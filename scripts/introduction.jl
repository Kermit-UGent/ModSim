### A Pluto.jl notebook ###
# v0.20.10

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ f99b9e24-a7d8-48e5-b71e-a28e6abe2177
begin
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ 7b432fab-7e4a-49b7-9ccb-9487001f89a4
using PlutoUI, Plots

# ╔═╡ 5d71d7a9-4bef-4988-bcab-bce5ad08448f
using Catalyst, DifferentialEquations

# ╔═╡ a46b6d63-64de-4c8d-b447-d6023c77dfb4
using ModelingToolkit: t_nounits as t, D_nounits as D

# ╔═╡ 2e3f042e-42bb-11ef-326c-f913245163cf
md"""
# Introduction
"""

# ╔═╡ 2760c5e4-2d92-40bb-8a63-d2dff4498a94
plots = Dict()

# ╔═╡ 5c6c0291-ae53-4abb-8ed7-4a5abcd65787
md"## Example 1"

# ╔═╡ 8f52f63a-43b4-4633-bac4-2f2a5ceb4c05
oilfield = @reaction_network begin
	hill(R, I * v, 20, 1), R --> C  # turning resource into capital
	hill(R, i, 20, 2), C --> I   # invest capital into infrastructure
	d, C --> 0   # depreciation of capital
	r, I --> 0  # degradation of infrastructure
end

# ╔═╡ ea8b91c5-621a-416f-9423-53f3970b370c
@bind R0 Slider(100:100:1000, default=100, show_value=true)

# ╔═╡ 3bcd4ec4-a65d-4561-b57f-805623d8183c
oilfieldproblem = ODEProblem(oilfield, [:R => R0, :I=>0.0, :C=>10.], (0., 100.),
			[:v => 0.05, :i=>0.1, :d=>0.05, :r=>0.2])

# ╔═╡ 07e15719-a544-4d88-a42c-152a875c2858
plots["oilfield"] = oilfieldproblem |> solve |> plot

# ╔═╡ 0cc6be58-add9-4a58-b532-a43e33da414c
fishery = @reaction_network begin
	hill(R, I * v, 20, 1), R --> C  # turning resource into capital
	hill(R, i, 20, 2), C --> I   # invest capital into infrastructure
	g * (1-R/R0), R --> 2R  # growth of the resource
	d, C --> 0   # depreciation of capital
	r, I --> 0  # degradation of infrastructure
end

# ╔═╡ d5f1ed22-e9c7-4164-b266-18cea9d15d14
fisheryproblem = ODEProblem(fishery, [:R => R0, :I=>0.0, :C=>10.], (0., 100.),
			[:v => 0.01, :i=>0.03, :d=>0.05, :r=>0.2, :R0=>R0, :g=>0.1])

# ╔═╡ ed33e20f-cc2c-47e3-add2-298f4b46b8d2
plots["fishery"] = fisheryproblem |> solve |> plot

# ╔═╡ 8ef0561c-de57-4e47-a830-9f3b507110f0
md"## Example model building"

# ╔═╡ 2837e47a-d733-4eb5-be74-8486616862a6
using ModelingToolkit

# ╔═╡ 6a4f086a-0570-43d6-acc3-ba178177aadc
@variables y(t)  # position y as a function of time

# ╔═╡ 333dfcca-74ea-4853-bb2a-c7d311b8d330
@parameters g=9.81 m=2.0 k ν=0

# ╔═╡ cb200949-b7ae-45c5-82e2-d4c81df43afc
eq = m * D(D(y)) ~ - m * g - k * y - ν * D(y)

# ╔═╡ 4e4a90c3-1053-4473-b19d-0bf10ef64344
@mtkbuild model = ODESystem(eq, t)

# ╔═╡ d012505e-b8a2-4a70-a4f8-09b320073eb0
symbolic_linear_solve(substitute(eq, D(y)=>0), y)

# ╔═╡ b364c127-344d-4eea-97d1-2c964c956456
oprob = ODEProblem(model, [y => -1, D(y) => 0], (0.0, 25.0), [k=>1.2, ν=>0.2])

# ╔═╡ dc932c9f-6960-411a-a5b9-6cbb3c1d4cc1
sol = solve(oprob, Tsit5())

# ╔═╡ bf70bd70-34dc-49f1-9f64-a1f87b79aefe
plots["spring"] = plot(sol)

# ╔═╡ Cell order:
# ╟─2e3f042e-42bb-11ef-326c-f913245163cf
# ╠═f99b9e24-a7d8-48e5-b71e-a28e6abe2177
# ╠═7b432fab-7e4a-49b7-9ccb-9487001f89a4
# ╠═2760c5e4-2d92-40bb-8a63-d2dff4498a94
# ╟─5c6c0291-ae53-4abb-8ed7-4a5abcd65787
# ╠═5d71d7a9-4bef-4988-bcab-bce5ad08448f
# ╠═8f52f63a-43b4-4633-bac4-2f2a5ceb4c05
# ╠═ea8b91c5-621a-416f-9423-53f3970b370c
# ╠═3bcd4ec4-a65d-4561-b57f-805623d8183c
# ╠═07e15719-a544-4d88-a42c-152a875c2858
# ╠═0cc6be58-add9-4a58-b532-a43e33da414c
# ╠═d5f1ed22-e9c7-4164-b266-18cea9d15d14
# ╠═ed33e20f-cc2c-47e3-add2-298f4b46b8d2
# ╠═8ef0561c-de57-4e47-a830-9f3b507110f0
# ╠═2837e47a-d733-4eb5-be74-8486616862a6
# ╠═a46b6d63-64de-4c8d-b447-d6023c77dfb4
# ╠═6a4f086a-0570-43d6-acc3-ba178177aadc
# ╠═333dfcca-74ea-4853-bb2a-c7d311b8d330
# ╠═cb200949-b7ae-45c5-82e2-d4c81df43afc
# ╠═4e4a90c3-1053-4473-b19d-0bf10ef64344
# ╠═d012505e-b8a2-4a70-a4f8-09b320073eb0
# ╠═b364c127-344d-4eea-97d1-2c964c956456
# ╠═dc932c9f-6960-411a-a5b9-6cbb3c1d4cc1
# ╠═bf70bd70-34dc-49f1-9f64-a1f87b79aefe
