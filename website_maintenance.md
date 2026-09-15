# Website maintenance

This document describes how the website works.

# Overview

This is the source code for the computational thinking website! It uses a site generation system inspired by [https://www.11ty.dev/](https://www.11ty.dev/), but there are only three template systems:
- **`.jlhtml` files** are rendered by [HypertextLiteral.jl](https://github.com/JuliaPluto/HypertextLiteral.jl)
- **`.jlmd` files** are rendered by [MarkdownLiteral.jl](https://github.com/JuliaPluto/MarkdownLiteral.jl)
- **`.jl` files** are rendered by [PlutoSliderServer.jl](https://github.com/JuliaPluto/PlutoSliderServer.jl)

The `/src/` folder is scanned for files, and all files are turned into HTML pages. 

Paths correspond to URLs. For example, `src/data_science/pca.jl` will become available at `https://kermit-ugent.github.io/ModSim/data_science/pca/`. For files called *"index"*, the URL will point to its parent, e.g. `src/docs/index.jlmd` becomes `https://kermit-ugent.github.io/ModSim/docs/`. Remember that changing URLs is very bad! You can't share this site with your friends if the links break.

> **To add something to our website, just create a new file!** Fons will be happy to figure out the technical bits.

You can generate & preview the website locally (more on this later), and we have a github action generating the website when we push to the `Fall23` branch. The result (in the `Fall23-output` branch) is deployed with GitHub Pages.

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
- *`tags`*: List of *tags* that are used to create collections out of pages. Our sidebar uses collections to know which pages to list. (more details in `sidebar data.jl`)
- *`layout`*: The name of a layout file in `src/_includes`. For basic Markdown or HTML, you probably want `md.jlmd`. For Pluto, you should use `layout.jlhtml`.

## How to write front matter
For `.jlmd` files, see the example above. 

For `.jl` notebooks, use the [Frontmatter GUI](https://github.com/fonsp/Pluto.jl/pull/2104) built into Pluto.

For `.jlhtml`, we still need to figure something out 😄.

# Running locally

## Developing *content, styles, etc.*

Open this repository in VS Code, and install the recommended extensions.

To start running the development server, open the VS Code *command palette* (press `Cmd+Shift+P`), and search for **`Tasks: Run Task`**, then **`PlutoPages: run development server`**. The first run can take some time, as it precompiles the packages in `pluto-deployment-environment`. Leave it running.

Use the same Julia version as the CI workflow (`1.11.2`, the version that `pluto-deployment-environment/Manifest.toml` was resolved with). With juliaup: `juliaup override set 1.11.2` inside this folder.

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

# PlutoPages.jl?

The PlutoPages.jl is still in experimental stage, and I'm not sure if the Julia community is waiting for another SSG system. So right now it's sort of released in secret. If you use it, be sure to respect our LICENSE and be sure to share your feedback with fons@plutojl.org! Bug reports welcome.
