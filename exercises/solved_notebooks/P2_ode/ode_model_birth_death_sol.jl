### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 12fee37d-eef9-4e90-b0af-d892469fad08
using Pkg; Pkg.activate("..")

# ╔═╡ 4a3d066c-f5b1-11ee-0145-2da7c11147a5
using Markdown, InteractiveUtils

# ╔═╡ 34a2f237-848c-4d07-9cea-fa5505a9e215
using StatsPlots, PlutoUI; TableOfContents()

# ╔═╡ b721cc4d-43ba-4221-bbd1-15293aaf54b5
# ╠═╡ show_logs = false
using OrdinaryDiffEq, Catalyst

# ╔═╡ c701d64e-640c-473f-b0fa-688024962f28
md"""
# Exercise: Simple birth-death model for mice
"""

# ╔═╡ 1cadc3db-af98-4f0d-874b-32f5b7852ad2
md"""
![](https://users.ugent.be/~gvhaelew/fig/mice_model_part1.png)
"""

# ╔═╡ c025dfec-df32-466c-b652-5a95e5415155
md"""
In a simple birth-death model for mice, the birth rate of mice represents the
rate at which new individuals are added to the population through reproduction.
This rate is influenced by factors such as the number of reproductive females,
their fertility, and the frequency of reproduction cycles. Conversely,
the death rate reflects the rate at which individuals are removed from
the population due to mortality factors such as predation, disease,
and environmental stressors. Together, these rates interact dynamically
to shape the population dynamics of mice in their natural habitat.
Denote the number of mice by $X$, the average birth rate by $b$ ($mice / day$), and the average death rate by $d$ ($day^{-1}$). Hence, assume for this overly simplified model, that the birth of mice is a zeroth-order process and that the death of mice is a first-order process.
"""

# ╔═╡ b583efef-dc7d-4447-bf27-e6373f5872c4
md"""
Create a *reaction network object* model for the aforementioned problem in order to simulate the evolution of $X$ with time. Name it `birth_death`.
"""

# ╔═╡ 74c3f7ed-d705-4cb9-b52c-e06c9df5ca13
birth_death = @reaction_network begin
    b, 0 --> X
    d, X --> 0
end

# ╔═╡ c84f2a32-e66b-4a36-a257-7031ce799225
md"""
Convert the system to a symbolic differential equation model and verify, by analyzing the differential equation, that your model is correctly implemented.
"""

# ╔═╡ 124eb20e-6945-4eec-a05a-6b835efbdd2e
osys  = convert(ODESystem, birth_death)

# ╔═╡ 04c645ad-873d-44c1-9bc0-940a412630b9
md"""
## Part 1

Simulate the evolution of the number of mice **per day** during $10$ years starting off with $2$ mice. **Assume that per year 25 pups are born.** Suppose the death rate to be $0.0015\;day^{-1}$.
"""

# ╔═╡ 1cd5136a-f472-4329-980a-b8382d1c04ea
md"""
First, calculate the birth rate in $mice/day$.
"""

# ╔═╡ 6b64ec8b-d85e-4a4f-aa3f-f177f5e34180
round(25/365, digits=4)

# ╔═╡ 52ee74d2-4dca-47ef-8b97-f826ef31ddf9
md"""
Initialize a vector `u0` with the initial conditions:
"""

# ╔═╡ 23f282d8-2eeb-43dd-8ba9-84684daca7a4
u0 = [:X => 2.0]

# ╔═╡ 30fa6fc6-54e0-4233-9334-b3d15f630489
md"""
Set the timespan for the simulation:
"""

# ╔═╡ b461a225-bae6-45b8-bb42-54742a72b98f
tspan = (0.0, 365*10.0)

# ╔═╡ c753d65f-9ed3-4e60-8e36-6e4bca77c19b
md"""
Initialize a vector `parms` with the parameter values:
"""

# ╔═╡ e7942ba9-5434-4b22-b1d4-f18f6227320d
parms = [:b => 0.0685, :d => 0.0015]

# ╔═╡ 2af7f6b9-9d95-4bf6-af67-449803440639
md"""
Create the ODE problem and store it in `oprob`:
"""

# ╔═╡ 9432f90b-7141-4518-b6fd-d55a1389e14a
oprob = ODEProblem(birth_death, u0, tspan, parms);

# ╔═╡ ab197ac7-413b-4305-8352-2723bf9f2aff
md"""
Solve the ODE problem. Use `Tsit5()` and `saveat=1.0`. Store the solution in `osol`:
"""

# ╔═╡ 8e517114-7c74-4d7c-954c-2787d975a1df
osol = solve(oprob, Tsit5(), saveat=1.0)

# ╔═╡ 9eaf1ca5-3d65-4263-a307-991bc8bf62d2
md"""
Plot the results:
"""

# ╔═╡ 83e36e63-7e6b-4c17-8f7d-ed0e07f36fc5
plot(osol)

# ╔═╡ 5fc64787-795c-4e17-b751-cd3683f8016b
md"""
!!! question
	Interpret the results. Ask yourself the following questions:
	- What is the (approximate) steady state value for $X$?
"""

# ╔═╡ f5718c08-c8aa-4d74-b70f-8c08c6decab2
md"- Answer: missing"

# ╔═╡ e26db23b-82ac-4a67-9580-a3125cd10fc8
osol[:X][end]

# ╔═╡ 50ec2f87-c10f-42c9-99bd-ad6ca496a6f7
u_guess = [:X => osol[:X][end]]

# ╔═╡ be84871f-dfc5-4cf9-98b5-8eab1472cd52
sol_eq = solve(SteadyStateProblem(birth_death, u_guess, parms))

# ╔═╡ 8a768ce3-bc6f-4ce2-bd7f-34d6718bde1d
sol_eq[:X]

# ╔═╡ fe4526c8-4504-4d91-841a-a83aeef55fc7
md"""
## Part 2

Suppose that at $t = 3\;years$ the death rate of the mice population increases by $50\,\%$ due to a new predator species in the area. Use the same initial condition, timespan and parameter values. Simulate the evolution of the number of mice.
"""

# ╔═╡ 1e869025-c537-4bf1-9c38-93e24b598156
md"""
Create the *condition*. Store it in `condition2`:
"""

# ╔═╡ 1c5d8daa-5c7e-49ef-a6da-6459a1131320
condition2 = [365*3.0] => [birth_death.d ~ birth_death.d*(1 + 50/100)]

# ╔═╡ 4d9ba697-5db3-4879-be88-278765e766f1
md"""
Make a new *reaction system* where the discrete event is included. Name it `birth_death2`.
"""

# ╔═╡ d40c2598-5783-49c5-916d-296e946de3f6
@named birth_death2 = ReactionSystem(equations(birth_death), discrete_events=condition2)

# ╔═╡ afbc4f00-0b4a-427a-9630-3a279ccbc330
md"""
Complete the new *reaction system*. Name it `birth_death2_com`.
"""

# ╔═╡ 3dd3febc-fc3f-477f-a5be-335ae90223e3
birth_death2_com = complete(birth_death2)

# ╔═╡ cf05864f-37f0-4fbb-93a5-ef406a2ae014
md"""
Create the ODE problem and store it in `oprob2`:
"""

# ╔═╡ 65e24cee-06c3-4cb7-8f8c-61e02b7dcd29
oprob2 = ODEProblem(birth_death2_com, u0, tspan, parms);

# ╔═╡ 9cab94b5-6d2d-4b5a-96f3-3b8eee1c61a1
md"""
Solve the ODE problem. Make a deepcopy and use `Tsit5()` and `saveat=1.0`. Store the solution in `osol2`:
"""

# ╔═╡ 00750aa4-6a00-4d33-a0ab-2d7abe952078
osol2 = solve(deepcopy(oprob2), Tsit5(), saveat=1.0)

# ╔═╡ 694e5b50-f4c9-4ecc-8860-1b217696e1ea
md"""
Plot the results:
"""

# ╔═╡ c69d0e86-0d15-4915-a890-8abd378da59a
plot(osol2)

# ╔═╡ 86e7af42-acb4-45e4-9a05-1c3ee2733dae
osol2[:X][end]

# ╔═╡ 1f869579-7019-4120-8d14-b38f48d05bdd
md"""
!!! questions
	Interpret the results. Ask yourself the following questions:
	1. Can you clearly see the effect of the increase in the death rate?
	2. If the death rate increases at a different timepoint, would you reach the same steady state value for $X$? Explain.
"""

# ╔═╡ 775f2cf3-ffbb-491a-9002-3a39b9c70f3c
md"""
Answers:
1. missing
2. missing
"""

# ╔═╡ Cell order:
# ╟─c701d64e-640c-473f-b0fa-688024962f28
# ╠═4a3d066c-f5b1-11ee-0145-2da7c11147a5
# ╠═12fee37d-eef9-4e90-b0af-d892469fad08
# ╠═34a2f237-848c-4d07-9cea-fa5505a9e215
# ╠═b721cc4d-43ba-4221-bbd1-15293aaf54b5
# ╟─1cadc3db-af98-4f0d-874b-32f5b7852ad2
# ╟─c025dfec-df32-466c-b652-5a95e5415155
# ╟─b583efef-dc7d-4447-bf27-e6373f5872c4
# ╠═74c3f7ed-d705-4cb9-b52c-e06c9df5ca13
# ╟─c84f2a32-e66b-4a36-a257-7031ce799225
# ╠═124eb20e-6945-4eec-a05a-6b835efbdd2e
# ╟─04c645ad-873d-44c1-9bc0-940a412630b9
# ╟─1cd5136a-f472-4329-980a-b8382d1c04ea
# ╠═6b64ec8b-d85e-4a4f-aa3f-f177f5e34180
# ╟─52ee74d2-4dca-47ef-8b97-f826ef31ddf9
# ╠═23f282d8-2eeb-43dd-8ba9-84684daca7a4
# ╟─30fa6fc6-54e0-4233-9334-b3d15f630489
# ╠═b461a225-bae6-45b8-bb42-54742a72b98f
# ╟─c753d65f-9ed3-4e60-8e36-6e4bca77c19b
# ╠═e7942ba9-5434-4b22-b1d4-f18f6227320d
# ╟─2af7f6b9-9d95-4bf6-af67-449803440639
# ╠═9432f90b-7141-4518-b6fd-d55a1389e14a
# ╟─ab197ac7-413b-4305-8352-2723bf9f2aff
# ╠═8e517114-7c74-4d7c-954c-2787d975a1df
# ╟─9eaf1ca5-3d65-4263-a307-991bc8bf62d2
# ╠═83e36e63-7e6b-4c17-8f7d-ed0e07f36fc5
# ╟─5fc64787-795c-4e17-b751-cd3683f8016b
# ╠═f5718c08-c8aa-4d74-b70f-8c08c6decab2
# ╠═e26db23b-82ac-4a67-9580-a3125cd10fc8
# ╠═50ec2f87-c10f-42c9-99bd-ad6ca496a6f7
# ╠═be84871f-dfc5-4cf9-98b5-8eab1472cd52
# ╠═8a768ce3-bc6f-4ce2-bd7f-34d6718bde1d
# ╟─fe4526c8-4504-4d91-841a-a83aeef55fc7
# ╟─1e869025-c537-4bf1-9c38-93e24b598156
# ╠═1c5d8daa-5c7e-49ef-a6da-6459a1131320
# ╟─4d9ba697-5db3-4879-be88-278765e766f1
# ╠═d40c2598-5783-49c5-916d-296e946de3f6
# ╟─afbc4f00-0b4a-427a-9630-3a279ccbc330
# ╠═3dd3febc-fc3f-477f-a5be-335ae90223e3
# ╟─cf05864f-37f0-4fbb-93a5-ef406a2ae014
# ╠═65e24cee-06c3-4cb7-8f8c-61e02b7dcd29
# ╟─9cab94b5-6d2d-4b5a-96f3-3b8eee1c61a1
# ╠═00750aa4-6a00-4d33-a0ab-2d7abe952078
# ╟─694e5b50-f4c9-4ecc-8860-1b217696e1ea
# ╠═c69d0e86-0d15-4915-a890-8abd378da59a
# ╠═86e7af42-acb4-45e4-9a05-1c3ee2733dae
# ╟─1f869579-7019-4120-8d14-b38f48d05bdd
# ╠═775f2cf3-ffbb-491a-9002-3a39b9c70f3c
