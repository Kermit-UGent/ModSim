### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 806ed5ad-788b-4f4b-b72f-3720766b6959
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ f3ae43dc-f7e1-11ee-3615-695b9e85b621
using Markdown

# ╔═╡ 28217175-2968-4219-81c9-eca0aa624196
using InteractiveUtils

# ╔═╡ 38b9cdf9-7d7a-44e5-aee4-e0b3a1c82464
using PlutoUI; TableOfContents()

# ╔═╡ 5d812731-8340-49c9-b181-aef5299388b5
using Catalyst

# ╔═╡ 63c243ff-649d-4494-98ef-e4dee5e88030
using DifferentialEquations, Plots

# ╔═╡ 2e8b2af9-2fed-49aa-a1cd-76c822664a5d
md"""
# Exercise: Irrigation experiment
"""

# ╔═╡ a7bb2db8-8e7e-45b6-843d-2db43f1a7ad1
md"""
An irrigation experiment is carried out on a soil column consisting of two layers of soil, each with specific soil characteristics. An adjustable volume of water per unit of time, $r$, is irrigated evenly over the soil column, starting with $5\;mm\,h^{-1}$ ($mm$ indicates a volume of water: $1\;mm = 10^{-3}\,m^3$). After $60\;h$ the added flow rate is increased to $10\;mm\,h^{-1}$. The water falls on the upper layer and percolates to the lower layer. The **relative moisture content** in both layers (i.e., relative to their residual moisture contents) is denoted by $S_1$ and $S_2$. Initially a moisture content of $30\;mm$ is present in the upper layer (cf. $S_1$) and of $25\;mm$ in the lower layer (cf. $S_2$). The residual moisture content in the upper layer is $S_{1,res}=10 \;mm$.

A model description of the relative moisture content in both soil layers is given by:

$$\begin{align}
\frac{dS_1}{dt} &= r\left(1-\cfrac{S_{1,res}}{S_{max}}\right) - \cfrac{r}{S_{max}}S_1 - \cfrac{k}{S_{max}}S_1 \\
\frac{dS_2}{dt} &= \cfrac{k}{S_{max}}S_1 - v \,S_2^2
\end{align}$$

Here $S_{max}$ ($150\;mm$) is the saturated water quantity for the top soil layer, $k$ is the percolation ratio ($3\;mm\,h^{-1}$) and $v$ is the flow factor into the groundwater ($10^{-3}\;h^{-1}\,mm^{-1}$).

Three measurements are made over the duration ($150\;h$) of the experiment:
- The excess running water (runoff): $r \cfrac{S_1 + S_{1,res}}{S_{max}}$,
- The underground outflow into groundwater: $v\,S_2^2$,
- The amount of percolation to deeper soil layers: $\cfrac{k}{S_{max}} S_1$.

The latter three are called *observables*.
"""

# ╔═╡ 072d96fa-4933-43ed-a4fd-162f64be5cdd
md"""
Create a *reaction network object* model for the aforementioned problem in order to simulate the evolution of the three afore mentioned measurements during $150\;h$. Name it `irrigation_mod`.

Tips:
- You can use any kind of expression for the reaction rates.
- The term $- v \,S_2^2$ is created by the reaction `v, 2S₂ --> 0`
"""

# ╔═╡ e557a067-52a1-42d5-98ae-40a0c63a2611
# Uncomment and complete the instruction
# irrigation_mod = @reaction_network begin
# 	missing
#   ...
# end
irrigation_mod = @reaction_network begin
    k/Smax, S₁ --> S₂
    v, 2S₂ --> 0
    r * (1 - S₁res / Smax), 0 --> S₁
    r/Smax, S₁ --> 0
end

# ╔═╡ 79277f44-ffca-44e3-867d-07a95dcb538c
md"""
Convert the system to a symbolic differential equation model and verify, by analyzing the differential equation, that your model is correctly implemented.
"""

# ╔═╡ 2f5b954b-1615-4b14-bc3b-3426c9296221
# osys = missing         # Uncomment and complete the instruction
osys = convert(ODESystem, irrigation_mod)

# ╔═╡ c163afdb-9411-4879-945a-8ee722ddb15a
md"""
Initialize a vector `u0` with the initial conditions:
"""

# ╔═╡ 32dcda4d-2e05-4451-86c2-22c2ea1d9b63
# u0 = missing           # Uncomment and complete the instruction
u0 = [:S₁ => 30, :S₂ => 25]

# ╔═╡ f974745b-5292-4bfa-aa8a-9836ae5870c3
md"""
Set the timespan for the simulation:
"""

# ╔═╡ 97c84665-c6ff-4c7a-8dfb-772244aa4986
# tspan = missing        # Uncomment and complete the instruction
tspan = (0.0, 150.0)

# ╔═╡ 0ceef9b8-4647-43a4-b9c3-11959834305b
md"""
Initialize a vector `params` with the parameter values:
"""

# ╔═╡ d172fc24-0fbe-46b4-a030-1ba44b4b57cd
# params = missing       # Uncomment and complete the instruction
params = [:k => 3.0, :Smax => 150.0, :v => 1.0e-3, :r => 5.0, :S₁res => 10.0]

# ╔═╡ 493d0d3d-cc43-4c2e-931f-f4acbc2f41ae
md"""
Unpack the variables and parameters so that we can use them in an intuitive way the calculate to observables.
"""

# ╔═╡ 67f62789-a5f0-4e05-a4fb-dbade2ad7891
# @unpack ..., ..., ..., ..., ..., ..., ... = ... # Uncomment and complete the instruction
@unpack S₁, S₂, k, Smax, v, S₁res, r = irrigation_mod

# ╔═╡ 7793eb61-3a2a-471a-b283-20db4a059b70
md"""
Create the *condition* that contains the timepoint for the sudden change in $r$. Store it in `condition`:
"""

# ╔═╡ c2bb92b9-ad9e-4bbc-8914-822848aea552
# condition = missing     # Uncomment and complete the instruction
condition = [60] => [irrigation_mod.r ~ irrigation_mod.r + 5.0]

# ╔═╡ 62ef97da-aedf-446d-ba2c-8c0ac72c9f05
md"""
Make a new *reaction system* where the discrete event is included. Name it `irrigation_mod_c`.
"""

# ╔═╡ 596b411b-2e76-48af-b3ce-8b73aa6f8700
# @named irrigation_mod_c = missing
@named irrigation_mod_c = ReactionSystem(equations(irrigation_mod), discrete_events=condition)

# ╔═╡ a38e9488-af43-49c6-9e7b-51f15e9a65ab
md"""
Complete the new *reaction system*. Name it `irrigation_mod_c_com`.
"""

# ╔═╡ a376c25d-059c-4a39-8c71-c8eb4310df72
# irrigation_mod_c_com = missing    # Uncomment and complete the instruction
irrigation_mod_c_com = complete(irrigation_mod_c)

# ╔═╡ 147b012a-298f-4885-9e71-453151c12bee
md"""
Create the ODE problem and store it in `oprob`:
"""

# ╔═╡ f222f889-d602-4654-9c90-a813738377c3
# oprob = missing        # Uncomment and complete the instruction
oprob = ODEProblem(deepcopy(irrigation_mod_c_com), u0, tspan, params)

# ╔═╡ 18717738-ae99-4b41-8414-4c1823307dde
md"""
Solve the ODE problem. Use `Tsit5()` and `saveat=0.5`. Store the solution in `osol`:
"""

# ╔═╡ dd5ffca7-4b39-4089-afee-dd73e4a9ac15
# osol = missing        # Uncomment and complete the instruction
osol = solve(deepcopy(oprob), Tsit5(), saveat=0.5)

# ╔═╡ ba5adc9e-5339-45e2-a126-fbdcdd0f283e
md"""
Calculate the observables.

- The excess running water (runoff): $r \cfrac{S_1 + S_{1,res}}{S_{max}}$,
- The underground outflow into groundwater: $v\,S_2^2$,
- The amount of percolation to deeper soil layers: $\cfrac{k}{S_{max}} S_1$.
"""

# ╔═╡ aeebc2d2-2052-4b80-a537-90d073ac3fe2
# Uncomment and complete the instruction
# begin
# 	runoff = ...;
# 	outflow = ...;
# 	percolation = ...;
# end
begin
	runoff = r*(S₁ + S₁res)/Smax;
	outflow = v*(S₂)^2;
	percolation = k/Smax*S₁;
end;

# ╔═╡ 3e1beea6-0fe4-4d77-9e5d-0ca78a662253
md"""
Plot the runoff, outflow and percolation.
"""

# ╔═╡ 15f14688-cf6e-4562-bcca-d3d0774c2d5c
# plot(osol; idxs=[..., ..., ...], labels=["..." "..." "..."])
plot(osol; idxs=[runoff, outflow, percolation], labels=["runoff" "outflow" "percolation"])

# ╔═╡ d1a9c363-3ef6-4bbb-a150-63b44f9e7cf1
md"""
Interpret the results. Ask yourself the following questions:

1. Can you clearly see the effect of the increase in $r$?
"""

# ╔═╡ 63ca9bf9-8b49-41ee-a062-eecc54f88e28
md"- Answer: missing"
#=
1. Yes, the effect of the increase in r can clearly been seen at t=60 in the three observables.
=#

# ╔═╡ 85e8869e-5859-4b8f-a72d-f0c979191ec7
md"""
2. Argue why the outflow and the percolation tend to the same value.
"""

# ╔═╡ 8b2459ca-9982-4c3a-b872-184c34289ec9
md"- Answer: missing"
#=
2. The water the percolates through the first soil layer will eventually as go through the second layer and flow out.
=#

# ╔═╡ Cell order:
# ╠═f3ae43dc-f7e1-11ee-3615-695b9e85b621
# ╠═28217175-2968-4219-81c9-eca0aa624196
# ╠═806ed5ad-788b-4f4b-b72f-3720766b6959
# ╠═38b9cdf9-7d7a-44e5-aee4-e0b3a1c82464
# ╟─2e8b2af9-2fed-49aa-a1cd-76c822664a5d
# ╟─a7bb2db8-8e7e-45b6-843d-2db43f1a7ad1
# ╠═5d812731-8340-49c9-b181-aef5299388b5
# ╟─072d96fa-4933-43ed-a4fd-162f64be5cdd
# ╠═e557a067-52a1-42d5-98ae-40a0c63a2611
# ╟─79277f44-ffca-44e3-867d-07a95dcb538c
# ╠═2f5b954b-1615-4b14-bc3b-3426c9296221
# ╠═63c243ff-649d-4494-98ef-e4dee5e88030
# ╟─c163afdb-9411-4879-945a-8ee722ddb15a
# ╠═32dcda4d-2e05-4451-86c2-22c2ea1d9b63
# ╟─f974745b-5292-4bfa-aa8a-9836ae5870c3
# ╠═97c84665-c6ff-4c7a-8dfb-772244aa4986
# ╟─0ceef9b8-4647-43a4-b9c3-11959834305b
# ╠═d172fc24-0fbe-46b4-a030-1ba44b4b57cd
# ╟─493d0d3d-cc43-4c2e-931f-f4acbc2f41ae
# ╠═67f62789-a5f0-4e05-a4fb-dbade2ad7891
# ╟─7793eb61-3a2a-471a-b283-20db4a059b70
# ╠═c2bb92b9-ad9e-4bbc-8914-822848aea552
# ╟─62ef97da-aedf-446d-ba2c-8c0ac72c9f05
# ╠═596b411b-2e76-48af-b3ce-8b73aa6f8700
# ╟─a38e9488-af43-49c6-9e7b-51f15e9a65ab
# ╠═a376c25d-059c-4a39-8c71-c8eb4310df72
# ╟─147b012a-298f-4885-9e71-453151c12bee
# ╠═f222f889-d602-4654-9c90-a813738377c3
# ╟─18717738-ae99-4b41-8414-4c1823307dde
# ╠═dd5ffca7-4b39-4089-afee-dd73e4a9ac15
# ╟─ba5adc9e-5339-45e2-a126-fbdcdd0f283e
# ╠═aeebc2d2-2052-4b80-a537-90d073ac3fe2
# ╟─3e1beea6-0fe4-4d77-9e5d-0ca78a662253
# ╠═15f14688-cf6e-4562-bcca-d3d0774c2d5c
# ╟─d1a9c363-3ef6-4bbb-a150-63b44f9e7cf1
# ╠═63ca9bf9-8b49-41ee-a062-eecc54f88e28
# ╟─85e8869e-5859-4b8f-a72d-f0c979191ec7
# ╠═8b2459ca-9982-4c3a-b872-184c34289ec9
