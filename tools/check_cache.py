#!/usr/bin/env python3
"""Check that _cache/ matches the notebooks in src/ and the pinned Pluto version.

Run before the expensive export job so a stale cache fails in seconds instead
of after a six-hour build that gets killed by the GitHub job limit.

The cached notebook states are named by PlutoSliderServer as

    <pluto version><content hash>.plutostate      e.g. 0_20_13AbC...xyz.plutostate

    src/Export.jl:56    escapeuri(string(pluto_version, hash)), "." => "_"
    src/PlutoHash.jl    plutohash = base64urlencode . sha256   (of the file bytes)

Both halves matter:

  * content hash -- edit a notebook and only that notebook re-runs. Normal, cheap.
  * pluto version -- bump Pluto and EVERY entry is invalidated at once, so the
    next build cold-runs all notebooks. That does not fit in 6h. This is the
    case the script hard-fails on.

Deliberately plain Python: no Julia, no Pkg.instantiate, runs in ~2 seconds.
Keep it in sync with PlutoSliderServer if that naming ever changes.
"""

import base64
import fnmatch
import hashlib
import pathlib
import re
import sys

try:
    import tomllib  # stdlib since Python 3.11; ubuntu-latest runners ship 3.12+
except ModuleNotFoundError:
    sys.exit(f"FAIL: needs Python 3.11+ for tomllib, got {sys.version.split()[0]}")

ROOT = pathlib.Path(__file__).resolve().parent.parent
MANIFEST = ROOT / "pluto-deployment-environment" / "Manifest.toml"
DEPLOY_TOML = ROOT / "pluto-deployment-environment" / "PlutoDeployment.toml"
SRC = ROOT / "src"
CACHE = ROOT / "_cache"

PLUTO_HEADER = b"### A Pluto.jl notebook ###"


def pluto_version() -> str:
    """Read the exact Pluto version the Manifest pins."""
    text = MANIFEST.read_text()
    m = re.search(r'^\[\[deps\.Pluto\]\]$(?:(?!^\[\[).)*?^version = "([^"]+)"',
                  text, re.S | re.M)
    if not m:
        sys.exit(f"FAIL: no Pluto version found in {MANIFEST.relative_to(ROOT)}")
    return m.group(1)


def ignored_globs() -> list[str]:
    """Notebooks PlutoDeployment.toml tells PlutoSliderServer not to cache."""
    if not DEPLOY_TOML.is_file():
        return []
    with DEPLOY_TOML.open("rb") as f:
        return tomllib.load(f).get("Export", {}).get("ignore_cache", [])


def is_notebook(path: pathlib.Path) -> bool:
    """Same test PlutoPages uses: Pluto.is_pluto_notebook, i.e. the header."""
    with path.open("rb") as f:
        return f.read(len(PLUTO_HEADER)) == PLUTO_HEADER


def statefile(path: pathlib.Path, prefix: str) -> str:
    digest = hashlib.sha256(path.read_bytes()).digest()
    h = base64.urlsafe_b64encode(digest).decode().rstrip("=")
    return f"{prefix}{h}.plutostate"


def main() -> int:
    version = pluto_version()
    prefix = version.replace(".", "_")
    skip = ignored_globs()

    notebooks = [
        p for p in sorted(SRC.rglob("*.jl"))
        if is_notebook(p)
        and not any(fnmatch.fnmatch(str(p.relative_to(SRC)), g)
                    or fnmatch.fnmatch(p.name, g) for g in skip)
    ]
    cache = {p.name for p in CACHE.glob("*.plutostate")} if CACHE.is_dir() else set()
    want = {statefile(p, prefix): p for p in notebooks}

    missing = {fn: p for fn, p in want.items() if fn not in cache}
    orphans = sorted(cache - want.keys())
    wrong_version = [fn for fn in cache if not fn.startswith(prefix)]

    print(f"Pluto {version} | {len(notebooks)} notebooks | {len(cache)} cache entries")

    # Total invalidation: a Pluto bump, or an empty cache. Unrecoverable in CI.
    if cache and len(wrong_version) == len(cache):
        sys.exit(
            f"\nFAIL: every cache entry was built by a different Pluto version.\n"
            f"  Manifest pins Pluto {version}, so entries must start with '{prefix}'\n"
            f"  but _cache/ holds e.g. {sorted(wrong_version)[0]}\n\n"
            f"  Cold-running all {len(notebooks)} notebooks exceeds the 6h job limit.\n"
            f"  Either revert the Pluto bump, or regenerate the cache locally:\n"
            f"      rm -rf _cache && julia --project=pluto-deployment-environment generate.jl\n"
            f"  and commit the result. See website_maintenance.md."
        )
    if not cache:
        sys.exit(
            f"\nFAIL: _cache/ is empty or missing.\n"
            f"  Cold-running all {len(notebooks)} notebooks exceeds the 6h job limit.\n"
            f"  See website_maintenance.md for how to regenerate it."
        )

    for fn in orphans:
        print(f"  orphan  {fn}")
    for fn, p in sorted(missing.items(), key=lambda kv: str(kv[1])):
        print(f"  cold    {p.relative_to(ROOT)}")

    if orphans:
        print(f"\n{len(orphans)} orphaned entr{'y' if len(orphans) == 1 else 'ies'} "
              f"(notebook edited or deleted). Safe to delete; do it with the next commit "
              f"so _cache/ does not grow without bound.")
    if missing:
        print(f"\n{len(missing)} notebook(s) will run cold this build. That is expected "
              f"after editing a notebook. Commit the regenerated _cache/ entries to keep "
              f"later builds fast.")
    if not missing and not orphans:
        print("\nCache is complete and current.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
