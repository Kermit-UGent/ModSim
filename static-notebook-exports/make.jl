using Pkg; Pkg.activate("./exercises/solved_notebooks")
using PlutoStaticHTML
using Documenter

cd(@__DIR__)
bopts = BuildOptions("../exercises/solved_notebooks/P1_mtk", output_format = documenter_output)
files = ["ode_model_tank_h_mtk_sol.jl"]
build_notebooks(bopts, files)

