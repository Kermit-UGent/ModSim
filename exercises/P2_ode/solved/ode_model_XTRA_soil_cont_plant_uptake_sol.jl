### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 6c7911b4-fec2-4139-8b54-36a4fb5916a0
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 121df656-f57a-11ee-140e-dfb61e112370
using Markdown

# ╔═╡ 74a1a82f-8c30-45e5-a2b7-1c1ce4d5c523
using InteractiveUtils

# ╔═╡ 895c1016-a8ce-43ae-8094-5d6ea75a6053
using Catalyst

# ╔═╡ 28e9e96c-fce6-4507-97c9-37337a0731bc
using DifferentialEquations, Plots

# ╔═╡ 2f0c3dd4-9429-4f5b-9150-011970b003f4
md"""
# Exercise: Soil Contamination with Plant Uptake
"""

# ╔═╡ c0c35547-9eed-428b-b513-4b5166decb7e
md"""
The following system of differential equations models the decay of a pollutant in soil and its uptake by plants. The variable $C(t)$ (in $mg/kg$) is the concentration of the pollutant in the soil and $P(t)$ (in $mg/kg$) is the concentration of the pollutant in the plants at time $t$. 

$$\begin{align}
\cfrac{dC}{dt} &= r-k_1 C(t)-k_2 C(t) P(t)\\
\cfrac{dP}{dt} &= k_2 C(t) P(t)-k_3 P(t)
\end{align}$$

The interpratation of the parameters is the following:

-  $r$ represents the rate at which the pollutant enters the soil from external sources.
-  $k_1$ is the natural degradation rate of the pollutant in the soil.
-  $k_2$ is the uptake coefficient, representing the rate at which plants absorb pollutant from the soil.
-  $k_3$ is the natural degradation rate of the pollutant in the plant.

The natural degradation of the pollutant in the soil or plant could be accounted for by processes like radiation decay, microbial degradation, volatilization, or leaching.
"""

# ╔═╡ 34e7090b-66c6-4e07-a2ad-b4f82a3669a0
md"""
Model the aforementioned system of differential equations using a *reaction network object*. Name it `soil_cont_plant_uptake`.
"""

# ╔═╡ ba6c2c9b-6ba9-48b8-9137-6b43813815ec
# Uncomment and complete the instruction
# soil_cont_plant_uptake = @reaction_network begin
#     missing
# end
soil_cont_plant_uptake = @reaction_network begin
	# @parameters K=2
    r, ∅ --> C
    k₁, C --> ∅
    k₂, C + P --> 2P
	# hill(C,k₂,K,4), C + P --> 2P
    k₃, P --> ∅
end

# ╔═╡ 16402872-3590-44dd-922f-1846640c92fa
md"""
Convert the system to a symbolic differential equation model and verify that you get the same system of differential equations as given in the problem.
"""

# ╔═╡ 1a7c3080-5773-44b9-a5c5-bb16f25048a3
# osys = missing         # Uncomment and complete the instruction
osys = convert(ODESystem, soil_cont_plant_uptake)

# ╔═╡ 7c4767ac-ee72-4f46-8699-68f4bfb15d92
md"""
Suppose that we simulate the evoluation of the pollutant in the soil and plant during $400$ days. The inital pollutant concentrations in the soil and plant both have the value of $0.001\;mg/kg$. In the simulation, the soil is being contaminated at a rate $0.06\;mg/(kg \cdot day)$. The degradation rates and uptake coefficient have the following values: $k_1 = 4.1 \times 10^{-3}$, $k_2 = 1.9 \times 10^{-2}$, $k_3 = 2.2 \times 10^{-2}$. There units are consistent with the units of the aforementioned values.
"""

# ╔═╡ 340328cc-78ed-4c6c-a0bf-73fd1be21d21
md"""
Initialize a vector `u0` with the initial conditions:
"""

# ╔═╡ 9d4fc31d-32b4-49c7-9e4c-577530199513
# u0 = missing         # Uncomment and complete the instruction
u0 = [:C => 0.001, :P => 0.001]

# ╔═╡ 51ffec7c-6033-47a9-b65a-5f9a7ab96fb9
md"""
Set the timespan for the simulation:
"""

# ╔═╡ 20b43337-58fd-4b23-8e4d-c4fd8234bb5f
# tspan =  missing      # Uncomment and complete the instruction
tspan = (0.0, 400.0)

# ╔═╡ ad63f799-bc65-4f2c-ba6c-461fa10139d0
md"""
Initialize a vector `param` with the parameter values:
"""

# ╔═╡ 01e33f75-fdfe-4983-a3bd-4cf074152390
# params = missing       # Uncomment and complete the instruction
params = [:r => 0.06, :k₁ => 4.1e-3, :k₂ => 1.9e-2, :k₃ => 2.2e-2]

# ╔═╡ 0f326aa7-044c-4c6f-be71-acf5c032f796
md"""
Create the ODE problem and store it in `oprob`:
"""

# ╔═╡ 6f84b532-a718-4d42-9829-91366693b51c
# oprob = missing        # Uncomment and complete the instruction
oprob = ODEProblem(soil_cont_plant_uptake, u0, tspan, params)

# ╔═╡ ca65797f-a1dd-42dd-992c-ba067932a018
md"""
Solve the ODE problem. Use `Tsit5()` and `saveat=1.0`. Store the solution in `osol`:
"""

# ╔═╡ 3d66d40f-f627-4268-890f-ab662c0efdd6
# osol = missing         # Uncomment and complete the instruction
osol = solve(oprob, Tsit5(), saveat=1.0)

# ╔═╡ 22dc63ad-47d5-45c4-8902-e9d7abc0a4f6
md"""
Plot the solutions:
"""

# ╔═╡ 0f920caa-5993-448c-a449-6449feea121a
# missing              # Uncomment and complete the instruction
plot(osol)

# ╔═╡ 62c66f8a-6561-4d50-b6d9-4dfc43cef0a8
md"""
1. Interprate the simulation results (cf. peak in $C$ and increase of $P$) in terms of the used parameter values.
"""

# ╔═╡ 27c1e08f-5e6a-43e3-b8e3-bba78e093556
md"- Answer: missing"
#=
In the beginning the pollution in the soil C strong increases because of the high contamination rate. The pollution in the plant stays low in the beginning because the natural decay is a bit higher than the uptake rate from the soil. At a certain moment the pollution in the soil will be so high that pollution in the plant will significantly increase. Eventually C and P will evolve to non zeros steady states.
=#

# ╔═╡ 73b63510-e1a1-44a6-8902-67ee383e6582
md"""
2. How would you modify the basic model to make it a more realistic biological model (cf. hill, monod, ...).
"""

# ╔═╡ a7068013-ed98-486a-ba59-a57f26d12d1c
md"- Answer: missing"
#=
We could suppose that the pollutant uptake by the plant is low when soil contamination is low and high, but limited, when soil contamination is high. In this case it would be suitable to use the Hill equation hill(C,k₂,K,n) that moderates the uptake rate.
=#

# ╔═╡ ad11f590-aa18-45f6-99be-571039ccbbae
md"""
3. What are the units of the parameters k₁, k₂ and k₃?
"""

# ╔═╡ 1ebbac82-068a-4754-8eec-1afa662feb96
md"- Answer: missing"
#=
C and P are in mg/kg, Hence:
- k₁ has unit: 1/s
- k₂ has unit: kg/(mgs)
- k₃ has unit: 1/s
=#

# ╔═╡ 12af799d-54aa-4bde-8133-98c3b4dfe1e4
# ╠═╡ disabled = true
#=╠═╡
(osol[:C][end], osol[:P][end])
  ╠═╡ =#

# ╔═╡ 3647ca59-c09e-40f8-8a95-03741e198428
# ╠═╡ disabled = true
#=╠═╡
u_guess = [osol[:C][end], osol[:P][end]]
  ╠═╡ =#

# ╔═╡ bf1ee7fb-027f-4b8b-8ce4-57a46ac1b4f4
# ╠═╡ disabled = true
#=╠═╡
Ceq, Peq = solve(SteadyStateProblem(ODEProblem(soil_cont_plant_uptake, u_guess, tspan, params)))
  ╠═╡ =#

# ╔═╡ 542d505f-8d65-4f68-8ac3-ee897b957b5e
# ╠═╡ disabled = true
#=╠═╡
(Ceq, Peq)
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═121df656-f57a-11ee-140e-dfb61e112370
# ╠═74a1a82f-8c30-45e5-a2b7-1c1ce4d5c523
# ╠═6c7911b4-fec2-4139-8b54-36a4fb5916a0
# ╟─2f0c3dd4-9429-4f5b-9150-011970b003f4
# ╟─c0c35547-9eed-428b-b513-4b5166decb7e
# ╠═895c1016-a8ce-43ae-8094-5d6ea75a6053
# ╟─34e7090b-66c6-4e07-a2ad-b4f82a3669a0
# ╠═ba6c2c9b-6ba9-48b8-9137-6b43813815ec
# ╟─16402872-3590-44dd-922f-1846640c92fa
# ╠═1a7c3080-5773-44b9-a5c5-bb16f25048a3
# ╠═28e9e96c-fce6-4507-97c9-37337a0731bc
# ╟─7c4767ac-ee72-4f46-8699-68f4bfb15d92
# ╟─340328cc-78ed-4c6c-a0bf-73fd1be21d21
# ╠═9d4fc31d-32b4-49c7-9e4c-577530199513
# ╟─51ffec7c-6033-47a9-b65a-5f9a7ab96fb9
# ╠═20b43337-58fd-4b23-8e4d-c4fd8234bb5f
# ╟─ad63f799-bc65-4f2c-ba6c-461fa10139d0
# ╠═01e33f75-fdfe-4983-a3bd-4cf074152390
# ╟─0f326aa7-044c-4c6f-be71-acf5c032f796
# ╠═6f84b532-a718-4d42-9829-91366693b51c
# ╟─ca65797f-a1dd-42dd-992c-ba067932a018
# ╠═3d66d40f-f627-4268-890f-ab662c0efdd6
# ╟─22dc63ad-47d5-45c4-8902-e9d7abc0a4f6
# ╠═0f920caa-5993-448c-a449-6449feea121a
# ╟─62c66f8a-6561-4d50-b6d9-4dfc43cef0a8
# ╠═27c1e08f-5e6a-43e3-b8e3-bba78e093556
# ╟─73b63510-e1a1-44a6-8902-67ee383e6582
# ╠═a7068013-ed98-486a-ba59-a57f26d12d1c
# ╟─ad11f590-aa18-45f6-99be-571039ccbbae
# ╠═1ebbac82-068a-4754-8eec-1afa662feb96
# ╠═12af799d-54aa-4bde-8133-98c3b4dfe1e4
# ╠═3647ca59-c09e-40f8-8a95-03741e198428
# ╠═bf1ee7fb-027f-4b8b-8ce4-57a46ac1b4f4
# ╠═542d505f-8d65-4f68-8ac3-ee897b957b5e
