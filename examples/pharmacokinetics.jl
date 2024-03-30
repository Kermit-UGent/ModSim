### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 44489e9a-e767-11ee-3a73-39be5e6caedb
begin
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ 26a1ffe8-ac61-44ff-a77d-0fc37b72737d
using Catalyst, DifferentialEquations, Plots, PlutoUI

# ╔═╡ 9b0f4b7a-cd3f-445e-b852-1d8715a897c6
rs = @reaction_network begin
	# Note: adsorption not dependent on concentration!
    ka, Agut => Abody
	kfat, Abody → Afat
	kbody, Afat → Abody
	k, Abody → ∅
end

# ╔═╡ 103f67e3-223b-4cf2-9bca-6185698230f0
Graph(rs)

# ╔═╡ 76c2fd28-fbf1-4e11-b297-9ab1bf455f78
complexgraph(rs)

# ╔═╡ b87de75f-a680-426b-9739-2d6cb9e7f742
let
	p = (
		:ka => 0.1,
		:kfat => 0.2,
		:kbody => 0.5,
		:k => 0.05,
	)
	tspan = (0., 100.)
	u0 = [:Agut => 10., :Afat => 0., :Abody => 0.]
	# solve ODEs
	oprob = ODEProblem(rs, u0, tspan, p)
	# NOTE: callback is used to avoid negative concentrations
	osol  = solve(oprob, Tsit5(), callback = PositiveDomain([1.0, 1.0]; save = false, abstol = 1E-9))
	plot(osol)
end

# ╔═╡ Cell order:
# ╠═44489e9a-e767-11ee-3a73-39be5e6caedb
# ╠═26a1ffe8-ac61-44ff-a77d-0fc37b72737d
# ╠═9b0f4b7a-cd3f-445e-b852-1d8715a897c6
# ╠═103f67e3-223b-4cf2-9bca-6185698230f0
# ╠═76c2fd28-fbf1-4e11-b297-9ab1bf455f78
# ╠═b87de75f-a680-426b-9739-2d6cb9e7f742
