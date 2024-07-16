### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
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

# ╔═╡ 2e3f042e-42bb-11ef-326c-f913245163cf
md"""
# Introduction
"""

# ╔═╡ 8f52f63a-43b4-4633-bac4-2f2a5ceb4c05
oilfield = @reaction_network begin
	hill(R, I * v, 20, 1), R --> C  # turning resource into capital
	hill(R, i, 20, 2), C --> I   # invest capital into infrastructure
	d, C --> 0   # depreciation of capital
	r, I --> 0  # degradation of investment
end

# ╔═╡ ea8b91c5-621a-416f-9423-53f3970b370c
@bind R0 Slider(100:100:1000, default=100, show_value=true)

# ╔═╡ 3bcd4ec4-a65d-4561-b57f-805623d8183c
oilfieldproblem = ODEProblem(oilfield, [:R => R0, :I=>0.0, :C=>10.], (0., 100.),
			[:v => 0.05, :i=>0.1, :d=>0.05, :r=>0.2])

# ╔═╡ 07e15719-a544-4d88-a42c-152a875c2858
oilfieldproblem |> solve |> plot

# ╔═╡ 0cc6be58-add9-4a58-b532-a43e33da414c
fishery = @reaction_network begin
	hill(R, I * v, 20, 1), R --> C  # turning resource into capital
	hill(R, i, 20, 2), C --> I   # invest capital into infrastructure
	g * (1-R/R0), R --> 2R  # growth of the resource
	d, C --> 0   # depreciation of capital
	r, I --> 0  # degradation of investment
end

# ╔═╡ d5f1ed22-e9c7-4164-b266-18cea9d15d14
fisheryproblem = ODEProblem(fishery, [:R => R0, :I=>0.0, :C=>10.], (0., 100.),
			[:v => 0.01, :i=>0.03, :d=>0.05, :r=>0.2, :R0=>R0, :g=>0.1])

# ╔═╡ ed33e20f-cc2c-47e3-add2-298f4b46b8d2
fisheryproblem |> solve |> plot

# ╔═╡ Cell order:
# ╠═2e3f042e-42bb-11ef-326c-f913245163cf
# ╠═f99b9e24-a7d8-48e5-b71e-a28e6abe2177
# ╠═7b432fab-7e4a-49b7-9ccb-9487001f89a4
# ╠═5d71d7a9-4bef-4988-bcab-bce5ad08448f
# ╠═8f52f63a-43b4-4633-bac4-2f2a5ceb4c05
# ╠═ea8b91c5-621a-416f-9423-53f3970b370c
# ╠═3bcd4ec4-a65d-4561-b57f-805623d8183c
# ╠═07e15719-a544-4d88-a42c-152a875c2858
# ╠═0cc6be58-add9-4a58-b532-a43e33da414c
# ╠═d5f1ed22-e9c7-4164-b266-18cea9d15d14
# ╠═ed33e20f-cc2c-47e3-add2-298f4b46b8d2
