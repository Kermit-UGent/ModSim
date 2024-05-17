### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 12fee37d-eef9-4e90-b0af-d892469fad08
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 4a3d066c-f5b1-11ee-0145-2da7c11147a5
using Markdown

# ╔═╡ 0ab4d97a-0c47-4f6e-a76b-9cec00ad410a
using InteractiveUtils

# ╔═╡ b721cc4d-43ba-4221-bbd1-15293aaf54b5
using Catalyst

# ╔═╡ 47e9c791-99cc-4a78-94b0-e0f5d4e0ecc5
using DifferentialEquations, Plots

# ╔═╡ c701d64e-640c-473f-b0fa-688024962f28
md"
### Exercise: Simple birth-death model for mice

In a simple birth-death model for mice, the birth rate of mice represents the
rate at which new individuals are added to the population through reproduction.
This rate is influenced by factors such as the number of reproductive females,
their fertility, and the frequency of reproduction cycles. Conversely,
the death rate reflects the rate at which individuals are removed from
the population due to mortality factors such as predation, disease,
and environmental stressors. Together, these rates interact dynamically
to shape the population dynamics of mice in their natural habitat.
Denote the number of mice by $X$, the average birth rate by $b$ ($mice / day$), and the average death rate by $d$ ($day^{-1}$).
"

# ╔═╡ b583efef-dc7d-4447-bf27-e6373f5872c4
md"
Create a *reaction network object* model for the aforementioned problem in order to simulate the evolution of $X$ with time. Name it `birth_death`.
"

# ╔═╡ 74c3f7ed-d705-4cb9-b52c-e06c9df5ca13
# Uncomment and complete the instruction
# birth_death = @reaction_network begin
# 	...
# end
birth_death = @reaction_network begin
    b, 0 --> X
    d, X --> 0
end

# ╔═╡ c84f2a32-e66b-4a36-a257-7031ce799225
md"
Convert the system to a symbolic differential equation model and verify, by analyzing the differential equation, that your model is correctly implemented.
"

# ╔═╡ 124eb20e-6945-4eec-a05a-6b835efbdd2e
# osys  = ...             # Uncomment and complete the instruction
osys  = convert(ODESystem, birth_death)

# ╔═╡ 04c645ad-873d-44c1-9bc0-940a412630b9
md"
#### Part 1

Simulate the evolution of the number of mice during 520 days starting off from 2 mice, with an average birth rate of $3\; mice/day$ and an average death rate of $0.015\;day^{-1}$.
"

# ╔═╡ 52ee74d2-4dca-47ef-8b97-f826ef31ddf9
md"
Initialize a vector `u₀` with the initial conditions:
"

# ╔═╡ 23f282d8-2eeb-43dd-8ba9-84684daca7a4
# u0 = ...         # Uncomment and complete the instruction
u0 = [:X => 2.0]

# ╔═╡ 30fa6fc6-54e0-4233-9334-b3d15f630489
md"
Set the timespan for the simulation:
"

# ╔═╡ b461a225-bae6-45b8-bb42-54742a72b98f
# tspan = ...      # Uncomment and complete the instruction
tspan = (0.0, 520.0)

# ╔═╡ c753d65f-9ed3-4e60-8e36-6e4bca77c19b
md"
Initialize a vector `param` with the parameter values:
"

# ╔═╡ e7942ba9-5434-4b22-b1d4-f18f6227320d
# params = ...     # Uncomment and complete the instruction
params = [:b => 3.0, :d => 0.015]

# ╔═╡ 2af7f6b9-9d95-4bf6-af67-449803440639
md"
Create the ODE problem and store it in `oprob`:
"

# ╔═╡ 9432f90b-7141-4518-b6fd-d55a1389e14a
# oprob = ...      # Uncomment and complete the instruction
oprob = ODEProblem(birth_death, u0, tspan, params)

# ╔═╡ ab197ac7-413b-4305-8352-2723bf9f2aff
md"
Solve the ODE problem. Use `Tsit5()` and `saveat=1.0`. Store the solution in `osol`:
"

# ╔═╡ 8e517114-7c74-4d7c-954c-2787d975a1df
# osol = ...         # Uncomment and complete the instruction
osol = solve(oprob, Tsit5(), saveat=1.0)

# ╔═╡ 9eaf1ca5-3d65-4263-a307-991bc8bf62d2
md"
Plot the results:
"

# ╔═╡ 83e36e63-7e6b-4c17-8f7d-ed0e07f36fc5
# ...           # Uncomment and complete the instruction
plot(osol)

# ╔═╡ 5fc64787-795c-4e17-b751-cd3683f8016b
md"
Interpret the results. Ask yourself the following questions:

1. What is the (approximate) steady state value for $X$?
- Answer: ...
"

# ╔═╡ fe4526c8-4504-4d91-841a-a83aeef55fc7
md"
#### Part 2

Suppose that at $t = 250\;days$ the death rate of the mice population increases by 30 % due to a new predator species in the area. Use the same initial condition, timespan and parameter values. Simulate the evolution of the number of mice.
"

# ╔═╡ 9ed6b4d1-16e1-483b-9984-c0c31effbdf2
md"
Create the *condition*. Store it in `condition2`:
"

# ╔═╡ 9252e205-c504-4c33-9cd6-f83d7f48f39e
# condition2 = ...      # Uncomment and complete the instruction
condition2 = [250.0]

# ╔═╡ 46c29f88-299f-4338-8d27-0e623fecdd8a
md"
Determine the order of the relevant parameter that needs to be modified at $t = 250 \;days$ in the model:
"

# ╔═╡ 170336ac-2db3-4afe-96ca-26ffe4873e73
# ...                # Uncomment and complete the instruction
parameters(birth_death)

# ╔═╡ 6a20f731-6006-49cd-8282-52c196e660fd
md"
Create the function called `affect2!`, that will be called by the solver at the timepoint(s) stored in `condition2` in order to alter the relevant parameter value:
"

# ╔═╡ 000722d1-e2dc-4ac0-9607-6d72055fd15e
# Uncomment and complete the instruction
# function affect2!(integrator)
#     ...
# end
function affect2!(integrator)
    integrator.p[2] *= (1 + 30/100)      # death rate d is the 2nd parameter
end

# ╔═╡ 92b47c79-9336-437e-82d3-d72f26078d64
md"
Create the callback function using `condition2` and `affect2!`. Store it in `cb2`:
"

# ╔═╡ 363a84ac-4cf3-411f-adb4-adaad8ce6cee
# cb2 = ...        # Uncomment and complete the instruction
cb2 = PresetTimeCallback(condition2, affect2!)

# ╔═╡ 83071c07-b753-439f-b47c-8467ab95745d
md"
Solve the ODE problem. Use `Tsit5()` and `saveat=1.0`. Store the solution in `osol2`:
"

# ╔═╡ 7dd44c1a-44cf-45b5-901a-84164ae733ab
# osol2 = ...         # Uncomment and complete the instruction
osol2 = solve(deepcopy(oprob), Tsit5(), saveat=1.0, callback=cb2)

# ╔═╡ 35a146ac-6c36-44f0-bef2-516e669e9e8b
md"
Plot the results:
"

# ╔═╡ 1bdef5d4-3bd8-44cf-89d1-c9f9c03f3382
# ...
plot(osol2)

# ╔═╡ 2663bcf6-8236-402b-af0a-b8639c356785
osol2[end]

# ╔═╡ c84cf2f3-f645-4dde-b5f2-d1180e781c24
md"
Interpret the results. Ask yourself the following questions:

1. Can you clearly see the effect of the increase in the death rate?
- Answer: ...
2. If the death rate increases at a different timepoint, would you reach the same steady state value for $X$? Explain.
- Answer: ...
"


# ╔═╡ Cell order:
# ╠═4a3d066c-f5b1-11ee-0145-2da7c11147a5
# ╠═0ab4d97a-0c47-4f6e-a76b-9cec00ad410a
# ╠═12fee37d-eef9-4e90-b0af-d892469fad08
# ╠═c701d64e-640c-473f-b0fa-688024962f28
# ╠═b721cc4d-43ba-4221-bbd1-15293aaf54b5
# ╠═b583efef-dc7d-4447-bf27-e6373f5872c4
# ╠═74c3f7ed-d705-4cb9-b52c-e06c9df5ca13
# ╠═c84f2a32-e66b-4a36-a257-7031ce799225
# ╠═124eb20e-6945-4eec-a05a-6b835efbdd2e
# ╠═47e9c791-99cc-4a78-94b0-e0f5d4e0ecc5
# ╠═04c645ad-873d-44c1-9bc0-940a412630b9
# ╠═52ee74d2-4dca-47ef-8b97-f826ef31ddf9
# ╠═23f282d8-2eeb-43dd-8ba9-84684daca7a4
# ╠═30fa6fc6-54e0-4233-9334-b3d15f630489
# ╠═b461a225-bae6-45b8-bb42-54742a72b98f
# ╠═c753d65f-9ed3-4e60-8e36-6e4bca77c19b
# ╠═e7942ba9-5434-4b22-b1d4-f18f6227320d
# ╠═2af7f6b9-9d95-4bf6-af67-449803440639
# ╠═9432f90b-7141-4518-b6fd-d55a1389e14a
# ╠═ab197ac7-413b-4305-8352-2723bf9f2aff
# ╠═8e517114-7c74-4d7c-954c-2787d975a1df
# ╠═9eaf1ca5-3d65-4263-a307-991bc8bf62d2
# ╠═83e36e63-7e6b-4c17-8f7d-ed0e07f36fc5
# ╠═5fc64787-795c-4e17-b751-cd3683f8016b
# ╠═fe4526c8-4504-4d91-841a-a83aeef55fc7
# ╠═9ed6b4d1-16e1-483b-9984-c0c31effbdf2
# ╠═9252e205-c504-4c33-9cd6-f83d7f48f39e
# ╠═46c29f88-299f-4338-8d27-0e623fecdd8a
# ╠═170336ac-2db3-4afe-96ca-26ffe4873e73
# ╠═6a20f731-6006-49cd-8282-52c196e660fd
# ╠═000722d1-e2dc-4ac0-9607-6d72055fd15e
# ╠═92b47c79-9336-437e-82d3-d72f26078d64
# ╠═363a84ac-4cf3-411f-adb4-adaad8ce6cee
# ╠═83071c07-b753-439f-b47c-8467ab95745d
# ╠═7dd44c1a-44cf-45b5-901a-84164ae733ab
# ╠═35a146ac-6c36-44f0-bef2-516e669e9e8b
# ╠═1bdef5d4-3bd8-44cf-89d1-c9f9c03f3382
# ╠═2663bcf6-8236-402b-af0a-b8639c356785
# ╠═c84cf2f3-f645-4dde-b5f2-d1180e781c24
