### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 2fb64158-6293-4a0b-b252-223307b472d9
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 71118b72-1db2-11ef-1f5b-a163b0b7c390
using Markdown

# ╔═╡ dbf68cf4-7a18-4bf2-90a3-e36216a41a70
using InteractiveUtils

# ╔═╡ b0b9313a-cc6b-44bd-90c9-5797f5d65d4c
using Catalyst

# ╔═╡ 99229de2-1e9b-470f-b532-ed1afb91c971
using DifferentialEquations, Plots

# ╔═╡ 5308b9d5-8068-48d7-b7a2-9e6ae043321b
md"
# Solving SDE problems with Catalyst
"

# ╔═╡ da37d5d2-4615-4be7-a34f-566f02ca2cc6
md"
**Stochastic Differential Equations** (SDEs) are mathematical equations used to model
systems influenced by random noise. They extend **Ordinary Differential Equations** (ODEs) by incorporating terms that represent **stochastic processes**, typically in the form of a Wiener process or Brownian motion. SDEs are widely used in various fields, such as physics, biology, finance, and engineering, to describe the evolution of systems under uncertainty or with inherent randomness.
"

# ╔═╡ 06e9a270-a116-4304-82a8-8733fd313b9c
md"
We will illustrate the concepts of SDE problems using the infection model that was used to introduce the Catalyst package and how to solve it as an ODE problem.\
In a previous notebook elaborated on the infection model in detail, hence, we will here limit ourselves to summarizing the variables, parameters and *reactions*.
"

# ╔═╡ 2b83ca01-0e45-4831-96ae-7a22f706a7c7
md"
Below we summarize the **variables** (**species**):
"

# ╔═╡ 7f3ffadd-5ae0-45d5-9894-49426f8c3af8
md"
| Variable | Unit | Meaning |
|:---------- |:---------- |:------------|
| ``S``    | *persons* | number of susceptible persons             |
| ``I``    | *persons* | number of infected persons             |
| ``D``    | *persons* | number of deceased persons             |
| ``R``    | *persons* | number of resistant persons             |
"

# ╔═╡ 7e72d0f8-50bb-460a-9820-91c38474e1a2
md"
Below we summarize the **parameters**:
"

# ╔═╡ 0ad755ec-9746-4716-9be3-aa92480c52c9
md"
| Variable | Unit | Meaning |
|:---------- |:---------- |:------------|
| ``\alpha`` | ``\frac{persons}{contact}`` | chances of getting infected after contact |
| ``\beta`` | ``\frac{contact}{persons^2\,day}`` | contact rate |
| ``r`` | ``\frac{1}{day}`` | rate of leaving infection period |
| ``m`` | ``\frac{person}{person}`` | fraction of persons deceasing |
| ``1-m`` | ``\frac{person}{person}`` | fraction of persons becoming resistant |
"

# ╔═╡ d46a18d4-f189-4be4-9bcb-edc6b05119cf
md"
Hence, the infection rate is ``\alpha \beta``. This means that a susceptible person meets an infected person: ``S+I``, this will result in ``2I`` at a rate ``\alpha \beta``. Futhermore, an infected person ``I`` will either become a deceased person ``D`` at a rate ``m r`` or become a resistant person ``R`` at rate ``(1-m) r``
"

# ╔═╡ 1e8b7f98-bede-4c43-9d4d-7fd94b1e73d3
md"
Our infection model has three reaction events:

- Infection, where a susceptible persons meets an infected persons and also becomes infected.
- Deceasing, where an infected person die.
- Recovery, where an infected person recovers and becomes resistant.
"

# ╔═╡ 274a11b8-09fc-46dd-9f27-ead38a41080d
md"
Each reaction is also associated with a specific rate:

- ``\alpha \beta``, the infection rate.
- ``m r``, the death rate.
- ``(1-m) r``, the recovery rate.
"

# ╔═╡ ef20da6f-c5af-4a94-bf6c-1ad68292992d
md"""
Hence, the following *infection reactions* are:

$$S + I \xrightarrow[]{\alpha \beta} 2I$$
$$I \xrightarrow[]{mr} D$$
$$I \xrightarrow[]{(1-m)r} R$$
"""

# ╔═╡ dcdadb6b-8577-4a59-80a5-ebad28b7d6b8
md"
We are going to implement this system of *reactions* using Catalyst.
"

# ╔═╡ 67ceffa7-2f2d-4719-b277-1eecf7c3ea57
md"
We first load the Catalyst package, which is required for the code in this introduction to run:
"

# ╔═╡ 5161d6a6-df87-4d50-9c74-5b7c780d4eb4
md"
#### Implementation of the system

The following code creates a so called *reaction network object*, that we have named `infection_sde_model`, that implements the aforementioned *reactions*.\
"

# ╔═╡ d1e14521-3365-43a5-b730-9b51ca359e12
infection_sde_model = @reaction_network begin
	@parameters η=50
	@default_noise_scaling η
	α * β, S + I --> 2I, [noise_scaling = 80.0]
	r * m, I --> D
	r * (1 - m), I --> R
end

# ╔═╡ f2b28b56-1369-4972-b771-7684397ad2ce
md"
Note that we have now introducted a **new parameter** (cf. `@parameters η` and `@default_noise_scaling η`). This parameter represent a default **noise scaling** parameter applying to all reactions. You can overwrite this default value for specific reactions by specifying `[noise_scaling = ...]` on the same line.\
They are in principle not necessary to solve the problem as a SDE problem but they can in many cases be very useful (see later).\
When solving SDE problems, some random noise will be introduced upon the reaction rates of all reactions.
"

# ╔═╡ 2a1afc01-5707-4368-ab69-cf2ffb8aa956
md"
Similarity as before, you can get a list of the *reaction* **species** with the function `species`, and a list of the **parameters** with the function `parameters`.
"

# ╔═╡ 2921c702-0ec3-4154-966c-7ba9c31095eb
parameters(infection_sde_model)

# ╔═╡ 619b6dc1-7486-4561-ab48-8daee7299734
md"
Note that the parameter $\eta$ also is present in the list.
"

# ╔═╡ 59c69709-f058-4312-9158-1de975ed20cd
md"
This *reaction model* can of course also be converted to a symbolic differential equation model:
"

# ╔═╡ 182667f0-e9cb-40b5-8b38-51a061c2c6de
osys = convert(ODESystem, infection_sde_model)

# ╔═╡ 56dce5a0-5917-4052-ae52-39683a5f8d31
md"
#### Simulating the system as an SDE problem

We first need to load the Differential and Plot package, which is required for simulating the system and plotting the results.
"

# ╔═╡ a4fd70ed-7767-471e-b1ff-58fd36a8f077
md"
Assume, as before, that there are $10\,000\,000$ people in the country, and that initially $1\,000$ person are infected. Hence, $I_0 = 1\,000$, $S_0 = 10\,000\,000-I_0 = 9\,999\,000$, $D_0 = 0$ and $R_0 = 0$.\
Furthermore, we take the following values for the parameters: $\alpha = 0.08\;person/contact$, $\beta = 10^{-6}\;contact/(person^2\,day)$, $r = 0.2\;day^{-1}$  and $m=0.4$.\
Finally, we want to run our simulation from day $0$ till day $90$.
"

# ╔═╡ b88e39ef-21a4-4a94-8145-3799348c80bd
md"
##### Setting initial conditions
"

# ╔═╡ e2ff6711-8e2b-43f7-9917-a567227a306c
u0 = [:S => 9999000, :I => 1000, :D => 0, :R => 0]

# ╔═╡ c0fad359-8a7f-4690-a405-e75b8625b868
md"
##### Setting the timespan
"

# ╔═╡ 1d0b68e1-b363-46e8-8445-9fc0765733eb
tspan = (0.0, 90.0)

# ╔═╡ b0a5e259-37d9-4480-971e-7365250706bb
md"
##### Setting parameter values

In the parameter list, you could also mention another default value for the default noise scaling parameter.
"

# ╔═╡ 38e28866-8c09-4675-a9c1-d9ba4012070b
params = [:α => 0.08, :β => 1.0e-6, :r => 0.2, :m => 0.4, :η => 50]

# ╔═╡ 5ef6414d-943c-4a1e-a528-af0fb61d4f9a
md"
##### Creating a SDEProblem

When creating the SDE problem you should literally state that the parameters `η1 η2 η3` are noise scaling parameters.
"

# ╔═╡ 1d5c6ff5-6a63-4491-981e-9c1ab5d37d60
# sprob = SDEProblem(infection_sde_model, u0, tspan, params, noise_scaling = @parameters η1 η2 η3)
sprob = SDEProblem(infection_sde_model, u0, tspan, params)

# ╔═╡ 0da84898-6170-4d71-acaf-23d97bfd7be9
md"
##### Solving the SDEProblem

There are many solving methods available for solving SDE problems. You can find a [list of methods here](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/#Full-List-of-Methods). We will simply use the first one in this list, namely `EM()`, with the time step option `dt=0.1` that will introduce some randomness at every time step.
"

# ╔═╡ f6ef847c-f9ea-464b-be95-e2bf5502a200
ssol = solve(sprob, EM(), dt=0.1)

# ╔═╡ 4fe4f1f8-b6a4-4e4a-b966-58bcdc2b2435
md"
Finally, we can plot the solution through the plot function.
"

# ╔═╡ b0253973-1a98-499d-a4b9-5db402a878bf
plot(ssol)

# ╔═╡ 96547ac1-d353-4a90-a6dd-38124d59294c
md"
You might notice that is you run the above instruction `ssol = solve(sprob, EM(), dt=0.1)` subsequent times, you will each time get different solutions (plots) due to the randomness introduced by treating the problem as SDE problem.
"

# ╔═╡ Cell order:
# ╠═71118b72-1db2-11ef-1f5b-a163b0b7c390
# ╠═dbf68cf4-7a18-4bf2-90a3-e36216a41a70
# ╠═2fb64158-6293-4a0b-b252-223307b472d9
# ╠═5308b9d5-8068-48d7-b7a2-9e6ae043321b
# ╠═da37d5d2-4615-4be7-a34f-566f02ca2cc6
# ╠═06e9a270-a116-4304-82a8-8733fd313b9c
# ╠═2b83ca01-0e45-4831-96ae-7a22f706a7c7
# ╠═7f3ffadd-5ae0-45d5-9894-49426f8c3af8
# ╠═7e72d0f8-50bb-460a-9820-91c38474e1a2
# ╠═0ad755ec-9746-4716-9be3-aa92480c52c9
# ╠═d46a18d4-f189-4be4-9bcb-edc6b05119cf
# ╠═1e8b7f98-bede-4c43-9d4d-7fd94b1e73d3
# ╠═274a11b8-09fc-46dd-9f27-ead38a41080d
# ╠═ef20da6f-c5af-4a94-bf6c-1ad68292992d
# ╠═dcdadb6b-8577-4a59-80a5-ebad28b7d6b8
# ╠═67ceffa7-2f2d-4719-b277-1eecf7c3ea57
# ╠═b0b9313a-cc6b-44bd-90c9-5797f5d65d4c
# ╠═5161d6a6-df87-4d50-9c74-5b7c780d4eb4
# ╠═d1e14521-3365-43a5-b730-9b51ca359e12
# ╠═f2b28b56-1369-4972-b771-7684397ad2ce
# ╠═2a1afc01-5707-4368-ab69-cf2ffb8aa956
# ╠═2921c702-0ec3-4154-966c-7ba9c31095eb
# ╠═619b6dc1-7486-4561-ab48-8daee7299734
# ╠═59c69709-f058-4312-9158-1de975ed20cd
# ╠═182667f0-e9cb-40b5-8b38-51a061c2c6de
# ╠═56dce5a0-5917-4052-ae52-39683a5f8d31
# ╠═99229de2-1e9b-470f-b532-ed1afb91c971
# ╠═a4fd70ed-7767-471e-b1ff-58fd36a8f077
# ╠═b88e39ef-21a4-4a94-8145-3799348c80bd
# ╠═e2ff6711-8e2b-43f7-9917-a567227a306c
# ╠═c0fad359-8a7f-4690-a405-e75b8625b868
# ╠═1d0b68e1-b363-46e8-8445-9fc0765733eb
# ╠═b0a5e259-37d9-4480-971e-7365250706bb
# ╠═38e28866-8c09-4675-a9c1-d9ba4012070b
# ╠═5ef6414d-943c-4a1e-a528-af0fb61d4f9a
# ╠═1d5c6ff5-6a63-4491-981e-9c1ab5d37d60
# ╠═0da84898-6170-4d71-acaf-23d97bfd7be9
# ╠═f6ef847c-f9ea-464b-be95-e2bf5502a200
# ╠═4fe4f1f8-b6a4-4e4a-b966-58bcdc2b2435
# ╠═b0253973-1a98-499d-a4b9-5db402a878bf
# ╠═96547ac1-d353-4a90-a6dd-38124d59294c
