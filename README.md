# Modelling and Simulation (ModSim)

Course material for *Modelling and Simulation of Biosystems* (I002445), Bachelor of Bioscience Engineering, Ghent University: course notes, exercise notebooks, examples and the course website.

## Launching Pluto

To launch Pluto, either:
- open a terminal in the top level folder, and run `julia launch_pluto.jl`
- open a Julia terminal in the top level folder, and run `include("launch_pluto.jl")`

## Repository layout

- `src/`: source of the course website (see [`website_maintenance.md`](website_maintenance.md)). Every file in here becomes a page.
- `exercises/`: the exercise notebooks.
  - `exercises/student_notebooks/` is the **single source of truth** for the exercises handed out to students.
  - `exercises/solved_notebooks/` contains the solved versions.
  - `src/exercises/*.jl` are **generated** copies for the website; never edit them by hand (see below).
- `notebook-checks/`: tooling around the notebooks.
  - `check_notebooks.jl` runs all solved notebooks; executed by CI (`.github/workflows/CheckNotebooks.yml`).
  - `sync_exercises.jl` regenerates `src/exercises/` from `exercises/student_notebooks/`; CI runs it with `--check` and fails when the website copies are out of sync.
- `project/`: the course project.
- `examples/`: worked examples used in the lectures.
- The course notes themselves (Typst sources, figures and the notebooks that generate them) live in the **private** repository `ModSim-course-notes`; they are not part of this repo.
- `pluto-deployment-environment/`: the Julia environment used to build the website (Julia 1.12). Keep its `Project.toml` and `Manifest.toml` up to date when you add packages to notebooks that are rendered on the website.

## Website

All website content lives in `src/`; the details are in [`website_maintenance.md`](website_maintenance.md).

- A push to `main` triggers `.github/workflows/ExportNotebooks.yml`, which builds the site and deploys it to the `gh-pages` branch (GitHub Pages). Pull requests get a preview under `previews/PR<n>`.
- To preview locally, run `julia develop.jl` in the top level folder (or use the VS Code task *PlutoPages: run development server*).
- Use Julia 1.12, the version the site environment was resolved with. With juliaup: `juliaup override set 1.12` inside this folder.

If you need a more in-depth example of the different pages, usage of tags, etc., check out release [v2425.1](https://github.com/Kermit-UGent/ModSim/releases/tag/v2425.1) and run the server on that code. Alternatively, have a look at the [original github source of the template](https://github.com/JuliaPluto/computational-thinking-template) or [the current website accompanying that repo](https://juliapluto.github.io/computational-thinking-template/).
