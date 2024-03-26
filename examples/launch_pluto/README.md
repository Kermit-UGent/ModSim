This folder explains the steps that are needed to run reproducible Pluto notebooks.
We define reproducible as:
- library versions are known on beforehand and fixed
- once installed and precompiled, the notebooks can be used as is and do not require redownloading.

To provide an example, we will use the [Measurements.jl](https://juliaphysics.github.io/Measurements.jl/stable/) package (used for quantification of measurement errors and error propagation).

Next up: launching the Pluto server.
Either via a terminal:

```shell
cd path/to/this/folder
julia launch_pluto.jl
```

or via a julia shell:

```julia
; cd /path/to/this/folder/
include("launch_pluto.jl")
```

In Pluto, open the `example_notebook.jl`.
The first cell in that notebook shows how to run that notebook in a specified environment instead of its own separate environment.
