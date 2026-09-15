# Website maintenance

This document describes how the website works.

# Overview

This is the source code for the ModSim course website. It uses a site generation system inspired by [https://www.11ty.dev/](https://www.11ty.dev/), but there are only three template systems:
- **`.jlhtml` files** are rendered by [HypertextLiteral.jl](https://github.com/JuliaPluto/HypertextLiteral.jl)
- **`.jlmd` and `.md` files** are rendered by [MarkdownLiteral.jl](https://github.com/JuliaPluto/MarkdownLiteral.jl)
- **`.jl` files** (Pluto notebooks) are rendered by [PlutoSliderServer.jl](https://github.com/JuliaPluto/PlutoSliderServer.jl)

The `/src/` folder is scanned for files, and all files are turned into HTML pages. 

Paths correspond to URLs. For example, `src/cheat_sheets/intro_to_julia.jl` will become available at `https://kermit-ugent.github.io/ModSim/cheat_sheets/intro_to_julia/`. For files called *"index"*, the URL will point to its parent, e.g. `src/index.jlmd` becomes `https://kermit-ugent.github.io/ModSim/`. Remember that changing URLs is very bad! You can't share this site with your friends if the links break.

Because the site is served under the `/ModSim/` prefix on GitHub Pages, use **relative links** between pages (e.g. `../intro_to_julia/`) rather than root-relative ones (`/cheat_sheets/intro_to_julia/`), which break under the prefix.

> **To add something to our website, just create a new file!**

You can generate & preview the website locally (more on this later). The GitHub Action `.github/workflows/ExportNotebooks.yml` generates the website on every push to `main` and deploys the result to the `gh-pages` branch with GitHub Pages. Pull requests get a preview under `previews/PR<n>` on the same branch.

# Content

## Literal templates
We use *Julia* as our templating system! Because we use HypertextLiteral and MarkdownLiteral, you can write regular Markdown files and HTML files, but you can also include `$(interpolation)` to spice up your documents! For example:

```markdown
# Hey there!

This is some *text*. Here is a very big number: $(1 + 1).
```

Besides small inline values, you can also write big code blocks, with `$(begin ... end)`, and you can output HTML. Take a look at some of our files to learn more!

## Pluto notebooks

Pluto notebooks are included in the page, but they are **not executed** during the site build.

Each notebook is loaded, its plain markdown cells are prerendered (this happens inside the build process, no notebook is started), and the result is embedded as a static Pluto editor. So a page shows the prose and the source code of every cell, but no cell outputs: no plots, no numbers, no interactive sliders.

Visitors run the notebook themselves with the **"Edit or run this notebook"** button in the top right: on Binder, or on their own computer with the repository checked out.

Because nothing is executed, there is **no notebook output cache** anymore. Editing a notebook is picked up immediately — you never have to invalidate anything, and a full site build takes a few minutes instead of hours.

## Exercise notebooks (sync)

The exercise pages under `src/exercises/` are **generated**, do not edit them by hand.

1. Edit the notebooks in `exercises/student_notebooks/` (the single source of truth).
2. Run `julia notebook-checks/sync_exercises.jl` from the top level folder to regenerate `src/exercises/`.
3. Commit both `exercises/` and `src/exercises/`.

CI runs `julia notebook-checks/sync_exercises.jl --check` and fails when `src/exercises/` is out of sync with `exercises/student_notebooks/`.

## `.css`, `.html`, `.gif`, etc

Web assets go through the system unchanged.

# Front matter

Like many SSG systems, we use [*front matter*](https://www.11ty.dev/docs/data-frontmatter/) to add metadata to pages. In `.jlmd` files, this is done with a front matter block, e.g.:
```markdown
---
title: "🌼 How to install"
description: "Instructions to install Pluto.jl"
tags: ["docs", "introduction"]
layout: "md.jlmd"
---

# Let's install Pluto

here is how you do it
```

Every page **should probably** include:
- *`title`*: Will be used in the sidebar, on Google, in the window header, and on social media.
- *`description`*: Will be used on hover, on Google, and on social media.
- *`tags`*: List of *tags* that are used to create collections out of pages. Our sidebar uses collections to know which pages to list. (more details in `src/_data/sidebar.jl`)
- *`layout`*: The name of a layout file in `src/_includes`. For basic Markdown or HTML, you probably want `md.jlmd`. For Pluto notebooks, you should use `layout.jlhtml`. The homepage uses `welcome.jlmd`.

## How to write front matter
For `.jlmd` files, see the example above. 

For `.jl` notebooks, use the [Frontmatter GUI](https://github.com/fonsp/Pluto.jl/pull/2104) built into Pluto.

For `.jlhtml`, we still need to figure something out 😄.

# Running locally

## Developing *content, styles, etc.*

Open this repository in VS Code, and install the recommended extensions.

To start running the development server, open the VS Code *command palette* (press `Cmd+Shift+P`), and search for **`Tasks: Run Task`**, then **`PlutoPages: run development server`** (this runs `julia develop.jl`). The first run can take some time, as it precompiles the packages in `pluto-deployment-environment`. Leave it running.

Use the same Julia version as the CI workflow (`1.12`, the version that `pluto-deployment-environment/Manifest.toml` was resolved with). With juliaup: `juliaup override set 1.12` inside this folder.

This will start two things in parallel: the PlutoPages.jl notebook (which generates the website), and a static file server (with Deno_jll). It will open two tabs in your browser: one is the generation dashboard (PlutoPages), the other is the current site preview (Deno_jll).
 
Whenever you edit a file, PlutoPages will automatically regenerate! Refresh your browser tab. If it does not pick up the change, go to the generation dashboard and click the "Read input files again" button.

This workflow is recommended for writing static content, styles, and for site maintenance. Notebooks are never executed by the site build, so editing one is cheap — but it also means you have to check the notebook's output in Pluto itself.

## Developing PlutoPages itself

You need to manually run the notebook with Pluto:
1. Go to this folder, and run `julia --project=pluto-deployment-environment`. Then `import Pkg; Pkg.instantiate();`.
1. `import Pluto; Pluto.run()` and open the `PlutoPages.jl` notebook in this repository. The first run can take some time, as it precompiles packages. Leave it running.
2. In a second terminal, go to this folder, and run `julia --project=pluto-deployment-environment`, then:
    ```julia
	import Deno_jll
	run(`$(Deno_jll.deno()) run --allow-read --allow-net https://deno.land/std@0.102.0/http/file_server.ts _site`)
    ```
3. Go to the URL printed to your terminal. 
4. Whenever you edit a file, PlutoPages will automatically regenerate! Refresh your browser tab. If it does not pick up the change, go to the generation dashboard and click the "Read input files again" button.

# PlutoPages.jl

The site generator is a vendored copy of [PlutoPages.jl](https://github.com/JuliaPluto/PlutoPages.jl) (the `PlutoPages.jl` notebook in the top level folder); see that repository for its documentation and license.
