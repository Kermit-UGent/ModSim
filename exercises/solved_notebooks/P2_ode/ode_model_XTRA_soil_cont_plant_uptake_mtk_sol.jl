### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 6c7911b4-fec2-4139-8b54-36a4fb5916a0
using Pkg; Pkg.activate("..")

# ╔═╡ 121df656-f57a-11ee-140e-dfb61e112370
using Markdown; InteractiveUtils;

# ╔═╡ 2c6324f9-1b00-48b3-a1af-8ece8b302f7d
using StatsPlots, PlutoUI; TableOfContents()

# ╔═╡ c54fee70-5250-43ad-8582-3d4ac2215041
using OrdinaryDiffEq, ModelingToolkit

# ╔═╡ 2e963f28-6d41-48be-a13f-6cd4bc1f7143
using ModelingToolkit: t_nounits as t, D_nounits as D

# ╔═╡ 2f0c3dd4-9429-4f5b-9150-011970b003f4
md"""
# Exercise: Soil Contamination with Plant Uptake
"""

# ╔═╡ 78dc9e01-ec30-4720-8fad-13e69a3146c4
md"""
![](https://users.ugent.be/~gvhaelew/fig/soil_plant_cont_model.png)
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

# ╔═╡ 7755d7e9-f8fd-476b-9178-4719754e4ba5
md"""
Model the aforementioned system of differential equations using ModelingToolkit.
"""

# ╔═╡ eef7afc7-b4b2-4a45-8459-291ec5b75038
md"""
Define the variables with `@variables`.
"""

# ╔═╡ 20c55b50-df74-4410-bc73-aa527c4f7dcc
@variables C(t) P(t)

# ╔═╡ 76dd7081-6a14-4151-b22a-4f330686d593
md"""
Define the parameters with `@parameters`.
"""

# ╔═╡ f3412cb7-4a16-44ce-82e7-62a183fe95df
@parameters r k₁ k₂ k₃

# ╔═╡ 7e8856f9-b452-44a0-976e-f415037239ad
md"""
Set up the equations. Name the set of equations `eqs_scpu`.
"""

# ╔═╡ 646b7f8f-1ce9-4a41-b274-cde8d016b929
eqs_scpu = [
	D(C) ~ r - k₁*C - k₂*C*P,
	D(P) ~ k₂*C*P - k₃*P
]

# ╔═╡ d8b04e69-ebf5-4659-9f0b-79f85ecfe59c
md"""
Build a system of equations with `@mtkbuild`. Name it `sys_scpu`.
"""

# ╔═╡ 38c44a2a-f3f1-4c2e-921e-7dbaa8748323
@mtkbuild sys_scpu = ODESystem(eqs_scpu, t)

# ╔═╡ 7c4767ac-ee72-4f46-8699-68f4bfb15d92
md"""
Suppose that we simulate the evoluation of the pollutant in the soil and plant during $400$ days. The inital pollutant concentrations in the soil and plant both have the value of $0.001\;mg/kg$. In the simulation, the soil is being contaminated at a rate $0.06\;mg/(kg \cdot day)$. The degradation rates and uptake coefficient have the following values: $k_1 = 4.1 \times 10^{-3}$, $k_2 = 1.9 \times 10^{-2}$, $k_3 = 2.2 \times 10^{-2}$. There units are consistent with the units of the aforementioned values.
"""

# ╔═╡ 340328cc-78ed-4c6c-a0bf-73fd1be21d21
md"""
Initialize a vector `u0` with the initial conditions:
"""

# ╔═╡ 9d4fc31d-32b4-49c7-9e4c-577530199513
u0 = [C=>0.001, P=>0.001]

# ╔═╡ 51ffec7c-6033-47a9-b65a-5f9a7ab96fb9
md"""
Set the timespan for the simulation:
"""

# ╔═╡ 20b43337-58fd-4b23-8e4d-c4fd8234bb5f
tspan = (0.0, 400.0)

# ╔═╡ ad63f799-bc65-4f2c-ba6c-461fa10139d0
md"""
Initialize a vector `parms` with the parameter values:
"""

# ╔═╡ 01e33f75-fdfe-4983-a3bd-4cf074152390
parms = [r=>0.06, k₁=>4.1e-3, k₂=>1.9e-2, k₃=>2.2e-2]

# ╔═╡ 0f326aa7-044c-4c6f-be71-acf5c032f796
md"""
Create the ODE problem and store it in `oprob`:
"""

# ╔═╡ 6f84b532-a718-4d42-9829-91366693b51c
oprob_scpu = ODEProblem(sys_scpu, u0, tspan, parms);

# ╔═╡ ca65797f-a1dd-42dd-992c-ba067932a018
md"""
Solve the ODE problem. Use `Tsit5()` and `saveat=1.0`. Store the solution in `osol`:
"""

# ╔═╡ 3d66d40f-f627-4268-890f-ab662c0efdd6
osol_scpu = solve(oprob_scpu, Tsit5(), saveat=1.0)

# ╔═╡ 22dc63ad-47d5-45c4-8902-e9d7abc0a4f6
md"""
Plot the solutions:
"""

# ╔═╡ 0f920caa-5993-448c-a449-6449feea121a
plot(osol_scpu)

# ╔═╡ 3cef5b2e-8186-49ae-8baa-b29454fba215
md"""
!!! question
	Interprate the simulation results (cf. peak in $C$ and increase of $P$) in terms of the used parameter values.
"""

# ╔═╡ 27c1e08f-5e6a-43e3-b8e3-bba78e093556
md"- Answer: missing"
#=
In the beginning the pollution in the soil C strong increases because of the high contamination rate. The pollution in the plant stays low in the beginning because the natural decay is a bit higher than the uptake rate from the soil. At a certain moment the pollution in the soil will be so high that pollution in the plant will significantly increase. Eventually C and P will evolve to non zeros steady states.
=#

# ╔═╡ 73b63510-e1a1-44a6-8902-67ee383e6582
md"""
!!! question
	How would you modify the basic model to make it a more realistic biological model (cf. hill, monod, ...).
"""

# ╔═╡ a7068013-ed98-486a-ba59-a57f26d12d1c
md"- Answer: missing"
#=
We could suppose that the pollutant uptake by the plant is low when soil contamination is low and high, but limited, when soil contamination is high. In this case it would be suitable to use the Hill equation hill(C,k₂,K,n) that moderates the uptake rate.
=#

# ╔═╡ ad11f590-aa18-45f6-99be-571039ccbbae
md"""
!!! question
	What are the units of the parameters k₁, k₂ and k₃?
"""

# ╔═╡ 1ebbac82-068a-4754-8eec-1afa662feb96
md"- Answer: missing"
#=
C and P are in mg/kg, Hence:
- k₁ has unit: 1/s
- k₂ has unit: kg/(mgs)
- k₃ has unit: 1/s
=#

# ╔═╡ 3a9900fc-bb80-4153-ae64-898e5e9fba1e
md"""
Calculate the steady state values for $C$ and $P$.
"""

# ╔═╡ 27a12024-a07d-4171-addb-48297ef72437
md"""
Initialize a vector `u_guess` with the end values of $C$ and $P$ from the solution.
"""

# ╔═╡ 3647ca59-c09e-40f8-8a95-03741e198428
u_guess = [osol_scpu[C][end], osol_scpu[P][end]]

# ╔═╡ 2f3ec6fb-bbed-4d28-ac3f-9ebce0121761
md"""
Make a state state problem and solve immediatelly.
"""

# ╔═╡ bf1ee7fb-027f-4b8b-8ce4-57a46ac1b4f4
eq = solve(SteadyStateProblem(sys_scpu, u_guess, parms))

# ╔═╡ a53dc273-4170-4eb6-8a25-dccbd1ad2c10
md"""
Assign the steady state values of $C$ and $P$ to `Ceq` and `Peq` respectively.
"""

# ╔═╡ 9cd89fbe-1add-4ba1-beef-8936a83f6987
Ceq = eq[C]; Peq = eq[P];

# ╔═╡ 828a2244-7e7f-4765-a1ff-f3dce2966bb9
md"""
Show the steady state values and compare with the end values of the solution.
"""

# ╔═╡ 542d505f-8d65-4f68-8ac3-ee897b957b5e
(Ceq, Peq)

# ╔═╡ 554f63ba-e3eb-4394-a9a8-9b7c0aaab78b
md"""
!!! question
	If the simulation time were longer, would the end values of the solution be closer to the steady state values?
"""

# ╔═╡ 90d1e3bf-0f38-4266-b076-2f9f0c5fea61
md"- Answer: missing"
#=
Yes!
=#

# ╔═╡ Cell order:
# ╟─2f0c3dd4-9429-4f5b-9150-011970b003f4
# ╠═121df656-f57a-11ee-140e-dfb61e112370
# ╠═6c7911b4-fec2-4139-8b54-36a4fb5916a0
# ╠═2c6324f9-1b00-48b3-a1af-8ece8b302f7d
# ╠═c54fee70-5250-43ad-8582-3d4ac2215041
# ╠═2e963f28-6d41-48be-a13f-6cd4bc1f7143
# ╟─78dc9e01-ec30-4720-8fad-13e69a3146c4
# ╟─c0c35547-9eed-428b-b513-4b5166decb7e
# ╟─7755d7e9-f8fd-476b-9178-4719754e4ba5
# ╟─eef7afc7-b4b2-4a45-8459-291ec5b75038
# ╠═20c55b50-df74-4410-bc73-aa527c4f7dcc
# ╟─76dd7081-6a14-4151-b22a-4f330686d593
# ╠═f3412cb7-4a16-44ce-82e7-62a183fe95df
# ╟─7e8856f9-b452-44a0-976e-f415037239ad
# ╠═646b7f8f-1ce9-4a41-b274-cde8d016b929
# ╟─d8b04e69-ebf5-4659-9f0b-79f85ecfe59c
# ╠═38c44a2a-f3f1-4c2e-921e-7dbaa8748323
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
# ╟─3cef5b2e-8186-49ae-8baa-b29454fba215
# ╠═27c1e08f-5e6a-43e3-b8e3-bba78e093556
# ╟─73b63510-e1a1-44a6-8902-67ee383e6582
# ╠═a7068013-ed98-486a-ba59-a57f26d12d1c
# ╟─ad11f590-aa18-45f6-99be-571039ccbbae
# ╠═1ebbac82-068a-4754-8eec-1afa662feb96
# ╟─3a9900fc-bb80-4153-ae64-898e5e9fba1e
# ╟─27a12024-a07d-4171-addb-48297ef72437
# ╠═3647ca59-c09e-40f8-8a95-03741e198428
# ╟─2f3ec6fb-bbed-4d28-ac3f-9ebce0121761
# ╠═bf1ee7fb-027f-4b8b-8ce4-57a46ac1b4f4
# ╟─a53dc273-4170-4eb6-8a25-dccbd1ad2c10
# ╠═9cd89fbe-1add-4ba1-beef-8936a83f6987
# ╟─828a2244-7e7f-4765-a1ff-f3dce2966bb9
# ╠═542d505f-8d65-4f68-8ac3-ee897b957b5e
# ╟─554f63ba-e3eb-4394-a9a8-9b7c0aaab78b
# ╠═90d1e3bf-0f38-4266-b076-2f9f0c5fea61
