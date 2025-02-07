### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ a2582acb-7d17-43ab-b883-d766b1a2c984
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ e5f8c320-eda0-11ee-37d0-458bdbd94f15
using Markdown

# ╔═╡ 2a1306f3-d811-459d-83ca-98cf62dc2db0
using InteractiveUtils

# ╔═╡ 4e03b93e-f63e-466c-9941-d66e62306010
using PlutoUI; TableOfContents()

# ╔═╡ 9e8fd818-a14f-41cf-b2fd-a7425141b283
using Catalyst

# ╔═╡ 0b3921c9-6d6e-4c52-8c21-d883ed493028
using DifferentialEquations, Plots

# ╔═╡ 62b185be-e327-4ef3-af39-819732d107bf
md"""
# Introduction to Catalyst (SSA)
"""

# ╔═╡ f5f32d5d-0c13-4865-8024-ca47208c9b8e
md"""
Catalyst.jl is a symbolic modeling package for analysis and high performance simulation of chemical reaction networks. Catalyst defines symbolic ReactionSystems, which can be created programmatically or easily specified using Catalyst's domain specific language (DSL).
"""

# ╔═╡ 8d9af65e-0499-4a6d-afbf-5afa9903a42e
md"""
This notebook describes the syntax for building chemical reaction network models using Catalyst's **D**omain-**S**pecific **L**anguage (DSL). We will illustrate this by implementing and solving an infection model by means of an SSA (**S**tochastic **S**imulation **A**lgorithm).
"""

# ╔═╡ 9cfa5b79-2128-4d45-aa22-51da0e74f320
md"""
## The infection model (revisited)
"""

# ╔═╡ e1583a47-9171-4db2-a6e9-d4889ee294c7
md"""
For the sake of clarity we restate the describtion of the previous infection model.

It is important to model the outbreak of infectious diseases in order to devise appropriate measures to avoid global epidemics. In this exercise we consider an isolated group of people in which a viral disease is spreading. An infection model (similar to the SIR-model but slightly extended) will be used for this purpose. We are interested in the evolution of the number of susceptible ($S$), infected ($I$), deceased ($D$) and resistant ($R$) persons.\
We make the following assumptions:
1. Transmission of the disease from an infected person to a susceptible person takes place through direct contact. The chance of any two inhabitants of the group coming into contact with each other is $\beta$, and the probability of infection after contact between an infected and a susceptible person is $\alpha$.
2. Note that the above assumption implicitly states that the probability of two neighbours coming into contact with each other is as high as the probability of two people living at two extremes of the territory coming into contact with each other.
3. A pereson leaves the infection period at a rate $r$ (hence, a person is contagious for an average of $1/r$ days. Without appropriate medication, a fraction $m$ of infected people die and a fraction $(1-m)$ of infected people acquire immunity after healing.
4. We assume that no one crosses the territory borders.
"""

# ╔═╡ 723c2c53-7f75-4f2d-8608-11ef0ef274d9
md"""
| Variable | Unit | Meaning |
|:---------- |:---------- |:------------|
| ``S``    | *persons* | number of susceptible persons             |
| ``I``    | *persons* | number of infected persons             |
| ``D``    | *persons* | number of deceased persons             |
| ``R``    | *persons* | number of resistant persons             |
"""

# ╔═╡ 0148d340-a072-49c0-9b3c-29249c21a334
md"""
| Variable | Unit | Meaning |
|:---------- |:---------- |:------------|
| ``\alpha`` | ``\frac{persons}{contact}`` | chances of getting infected after contact |
| ``\beta`` | ``\frac{contact}{persons^2\,day}`` | contact rate |
| ``r`` | ``\frac{1}{day}`` | rate of leaving infection period |
| ``m`` | ``\frac{person}{person}`` | fraction of persons deceasing |
| ``1-m`` | ``\frac{person}{person}`` | fraction of persons becoming resistant |
"""

# ╔═╡ 66d75717-3bcd-4e73-a8de-ed3879efb509
md"""
Hence, the infection rate is ``\alpha \beta``. This means that a susceptible person meets an infected person: ``S+I``, this will result in ``2I`` at a rate ``\alpha \beta``. Futhermore, an infected person ``I`` will either become a deceased person ``D`` at a rate ``m r`` or become a resistant person ``R`` at rate ``(1-m) r``
"""

# ╔═╡ e884ea42-87ac-4ccc-bf19-937837a2645d
md"""
Our infection model has three reaction events:

- Infection, where a susceptible persons meets an infected persons and also becomes infected.
- Deceasing, where an infected person die.
- Recovery, where an infected person recovers.
"""

# ╔═╡ 10c93224-ac8d-4512-be8f-6961b885712f
md"""
Each reaction is also associated with a specific rate:

- ``\alpha \beta``, the infection rate.
- ``m r``, the death rate.
- ``(1-m) r``, the recovery rate.
"""

# ╔═╡ 7064d3d4-7ccb-4af2-a1ed-64430c50651f
md"""
Hence, the following *infection reactions* are:

$$S + I \xrightarrow[]{\alpha \beta} 2I$$
$$I \xrightarrow[]{mr} D$$
$$I \xrightarrow[]{(1-m)r} R$$
"""

# ╔═╡ 4dd77999-2ed0-415e-b044-1af1de1b4ae3
md"""
We are going to implement this system of *reactions* using Catalyst.
"""

# ╔═╡ 30a0e2ec-ec30-401b-9270-939a767b54a8
md"""
We first load the Catalyst package, which is required for the code in this introduction to run:
"""

# ╔═╡ d888c734-044b-424d-b895-f8bc738346ed
md"""
### Implementation of the system

First we create a *reaction network object*, that we have named `infection_model`, that implements the aforementioned *reactions*.
"""

# ╔═╡ c792559f-db9b-4d9a-8e79-c7e8f82b4603
infection_model = @reaction_network begin
	α * β, S + I --> 2I
	r * m, I --> D
	r * (1 - m), I --> R
end

# ╔═╡ 2085fd23-28fc-4b68-a3e7-590c50f2c9c6
md"""
You can get a list of the different *reaction* **species** with the command `species`
"""

# ╔═╡ f3e1feb4-d36a-4095-8e25-df02c73078e9
species(infection_model)

# ╔═╡ 22145b7a-5a55-4db1-9061-ad51e907490f
md"""
The *reaction model* can be converted to a symbolic differential equation model via
"""

# ╔═╡ 7f15968c-dadc-4860-8071-9b072673e414
osys  = convert(ODESystem, infection_model)

# ╔═╡ f9feffe6-c7d6-4962-ada5-85a51664ff2a
md"""
You can get a list of the differential equations with the command `equations`:
"""

# ╔═╡ a0ee4837-074a-4587-a9d6-b8892a6c99f4
equations(osys)

# ╔═╡ 1d6f5b69-ab47-4344-985f-3e9c9430fddd
md"""
To get a list of the state variables, you can use the command `unknowns`:
"""

# ╔═╡ 6bc9574c-0f04-45c0-a457-6601f1333c0e
unknowns(osys)

# ╔═╡ b3584f8c-50f4-4800-96aa-776cfc2b8db3
md"""
To get a list of the parameters, you can use the command `parameters`:
"""

# ╔═╡ 98fa7987-b60b-4748-a3e8-38259bb0cd8c
parameters(osys)

# ╔═╡ d94a467a-3194-4573-a3e0-27a265147e66
md"""
### Simulating the system as a (Discrete) Jump problem
"""

# ╔═╡ d8b11407-e8c0-4ad0-a49d-8e422ccc3a9c
md"""
We first need to load the Differential and Plot package, which is required for simulating the system and plotting the results.
"""

# ╔═╡ 78e55389-3431-41de-9f4a-b6b45cb8988b
md"""
Instead of simulating our model with the species defined as decimal numbers, we will simulate the individual reaction events through the so-called **Gillespie algorithm**. This algorithm is a so-called **Stochastic Simulation Algorithm** (SSA).\
The Gillespie algorithm is a computational method used to simulate **discrete** and **stochastic** (random) processes. The algorithm models the changes in a system over time by considering individual events and their probabilities, this allows to understand how random fluctuations affect the system's behavior.
"""

# ╔═╡ 9bed9d84-600a-44f0-87a1-0997d547791f
md"""
To illustrate the simulation based on the Gillespie-algorithm, we will use the same infection model as before, but considering much less individuals. Hence, we will use different initial conditions, parameter values and timespan as with the ODE problem.
"""

# ╔═╡ a2349455-5850-4950-bb10-221ae813b26f
md"""
Assume in this example that there are $50$ people on the territory, and that initially 1 person is infected. Hence, $I_0 = 1$, $S_0 = 50-I_0 = 49$, $D_0 = 0$ and $R_0 = 0$.\
Furthermore, we take the following values for the parameters: $\alpha = 0.15\;person/contact$, $\beta = 0.1\;contact/(person^2\,day)$, $r = 0.2\;day^{-1}$ (i.e. a person is contagious for an average of $5\;days$) and $m=0.6$.\
Finally, we want to run our simulation from day $0$ till day $60$.
"""

# ╔═╡ 60b150fc-b9e5-426d-8d48-95efc907692a
md"""
### Setting initial conditions

The vector holding the initial conditions for $S$, $I$, $D$ and $R$ is:
"""

# ╔═╡ bca8459b-8f05-4e4e-92f1-43659ef652ba
u0 = [:S => 49, :I => 1, :D => 0, :R => 0]

# ╔═╡ ad8215f6-676e-4352-b0c9-8e9701da3bc6
md"""
### Setting parameter values

The vector holding the parameter values for $\alpha$, $\beta$, $r$ and $m$ is:
"""

# ╔═╡ da136f6e-7605-4f1f-81bc-f0e81ed7f528
params = [:α => 0.15, :β => 0.1, :r => 0.2, :m => 0.6]

# ╔═╡ cef064bf-2168-4904-b753-a012a9c9f070
md"""
### Setting the timespan
"""

# ╔═╡ 69b38ef3-d99e-4f4d-8c53-79e20d0094f0
tspan = (0.0, 60.0)

# ╔═╡ 57567d6d-6303-4ad8-b171-f904b594b3fe
md"""
### Creating an DiscreteProblem
"""

# ╔═╡ 92c77720-5120-4f55-8cf6-168ec8553638
md"""
Unlike the previous approach with ODEProblem (denoting a deterministic ordinary differential equation), we wish to simulate our model as a jump process (where each reaction event denotes a single jump in the state of the system). We do this by first creating a **DiscreteProblem**, and then using this as an input to a **JumpProblem**.
"""

# ╔═╡ 65e8b392-4a4c-4e09-a95e-a633a42e4bee
md"""
We create a DiscreteProblem by calling the `DiscreteProblem` function. Applying this function ensures that the problem is approached at a level of individual infections (reactions). Hence, the variable values will be integers. *Note that the order in which the input (the model name, the initial condition, the timespan, and the parameter values) is provided to* `DiscreteProblem` *matters!* Here, we save our DiscreteProblem in the `dprob` variable.
"""

# ╔═╡ e2f215e7-8d20-4e72-9e0b-cc0ce150e55b
dprob = DiscreteProblem(infection_model, u0, tspan, params)

# ╔═╡ d061c4e2-8fbb-4ac5-95cf-c37d2916e669
md"""
Next, we create a so-called JumpProblem by calling the `JumpProblem` function. Applying this function ensures that the infections (reactions) will happen stochastically. *Note again that the order in which the input (the model name, the DiscreteProblem variable, the simulation method) is provided to* `JumpProblem` *matters!* The simulation method is denoted by the option `Direct()`, which we recommend for now.
"""

# ╔═╡ b462661a-3908-48e7-b162-f3b78d37231f
jprob = JumpProblem(infection_model, dprob, Direct())

# ╔═╡ bb1538ec-99e4-4ad0-8aaa-d552c6291c6d
md"""
### Solving the DiscreteProblem
"""

# ╔═╡ 89c75cee-2488-4e9c-a470-b2111dcfc871
md"""
Finally, we can simulate our model using the solve function, and plot the solution using the `plot` function. Here, the `solve` function also has a second argument `SSAStepper()`, which we recommend for now. This is a time stepping algorithm that calls the `Direct` solver method to advance a simulation.
"""

# ╔═╡ 8fbb4942-888c-4798-b972-7860e695d5ba
dsol = solve(jprob, SSAStepper())

# ╔═╡ 5cfe9878-cbc2-4852-a757-ad8247d1f3d6
md"""
Note that at the different time points the variables values in the solution are  integer numbers and reflect the number of persons in either state ($S$, $I$, $D$ and $R$).\
Futhermore, note that executing the `solve` command at different occasions will result in other solutions because of the stochastic character of the applied method.
"""

# ╔═╡ ed8dce8b-d83b-40bc-b1eb-1ab369611cc1
md"""
Finally, we can plot the solution through the plot function.
"""

# ╔═╡ f855160a-2f04-4d11-97de-688113067c1c
plot(dsol)

# ╔═╡ 719af6d7-3c98-499e-b789-a6de5af2a027
md"""
Below is a piece of code that solves the problem a $1000$ times and stores the time values at which the number of infected persons becomes zero.
"""

# ╔═╡ bcd1e77b-0869-4b36-9309-598b8c771d80
begin
	times = []                      # make empty vector
	while length(times) < 1000      # while statement
		dsol2 = solve(jprob, SSAStepper())    # solve the problem
		j = findfirst(dsol2[:I] .== 0)        # find index of first 0
		if j != nothing                       # if index is a valid index
			append!(times, dsol2.t[j])        # append time to vector times
		end
	end
end

# ╔═╡ e97ca371-9ee5-4419-afcd-2a35ed7426d1
md"""
The vector `times` is now filled with time values at which the number of infected persons becomes zero.
"""

# ╔═╡ a0e70668-6c4e-4395-b1e4-925952a52dea
times

# ╔═╡ 764c7739-faad-4cf3-a64c-8ce3cae64074
md"""
With this vector we make a histogram so that you can have an idea of the distribution when the infected persons becomes zero.
"""

# ╔═╡ 1971359b-df64-49e5-9011-c992da211c74
histogram(times, bins=range(0, 60, length=61))
# histogram(times, bins=range(0, 60, length=61), normalize=:pdf)

# ╔═╡ Cell order:
# ╠═e5f8c320-eda0-11ee-37d0-458bdbd94f15
# ╠═2a1306f3-d811-459d-83ca-98cf62dc2db0
# ╠═a2582acb-7d17-43ab-b883-d766b1a2c984
# ╠═4e03b93e-f63e-466c-9941-d66e62306010
# ╟─62b185be-e327-4ef3-af39-819732d107bf
# ╟─f5f32d5d-0c13-4865-8024-ca47208c9b8e
# ╟─8d9af65e-0499-4a6d-afbf-5afa9903a42e
# ╟─9cfa5b79-2128-4d45-aa22-51da0e74f320
# ╟─e1583a47-9171-4db2-a6e9-d4889ee294c7
# ╟─723c2c53-7f75-4f2d-8608-11ef0ef274d9
# ╟─0148d340-a072-49c0-9b3c-29249c21a334
# ╟─66d75717-3bcd-4e73-a8de-ed3879efb509
# ╟─e884ea42-87ac-4ccc-bf19-937837a2645d
# ╟─10c93224-ac8d-4512-be8f-6961b885712f
# ╟─7064d3d4-7ccb-4af2-a1ed-64430c50651f
# ╟─4dd77999-2ed0-415e-b044-1af1de1b4ae3
# ╟─30a0e2ec-ec30-401b-9270-939a767b54a8
# ╠═9e8fd818-a14f-41cf-b2fd-a7425141b283
# ╟─d888c734-044b-424d-b895-f8bc738346ed
# ╠═c792559f-db9b-4d9a-8e79-c7e8f82b4603
# ╟─2085fd23-28fc-4b68-a3e7-590c50f2c9c6
# ╠═f3e1feb4-d36a-4095-8e25-df02c73078e9
# ╟─22145b7a-5a55-4db1-9061-ad51e907490f
# ╠═7f15968c-dadc-4860-8071-9b072673e414
# ╟─f9feffe6-c7d6-4962-ada5-85a51664ff2a
# ╠═a0ee4837-074a-4587-a9d6-b8892a6c99f4
# ╟─1d6f5b69-ab47-4344-985f-3e9c9430fddd
# ╠═6bc9574c-0f04-45c0-a457-6601f1333c0e
# ╟─b3584f8c-50f4-4800-96aa-776cfc2b8db3
# ╠═98fa7987-b60b-4748-a3e8-38259bb0cd8c
# ╟─d94a467a-3194-4573-a3e0-27a265147e66
# ╟─d8b11407-e8c0-4ad0-a49d-8e422ccc3a9c
# ╠═0b3921c9-6d6e-4c52-8c21-d883ed493028
# ╟─78e55389-3431-41de-9f4a-b6b45cb8988b
# ╟─9bed9d84-600a-44f0-87a1-0997d547791f
# ╟─a2349455-5850-4950-bb10-221ae813b26f
# ╟─60b150fc-b9e5-426d-8d48-95efc907692a
# ╠═bca8459b-8f05-4e4e-92f1-43659ef652ba
# ╟─ad8215f6-676e-4352-b0c9-8e9701da3bc6
# ╠═da136f6e-7605-4f1f-81bc-f0e81ed7f528
# ╟─cef064bf-2168-4904-b753-a012a9c9f070
# ╠═69b38ef3-d99e-4f4d-8c53-79e20d0094f0
# ╟─57567d6d-6303-4ad8-b171-f904b594b3fe
# ╟─92c77720-5120-4f55-8cf6-168ec8553638
# ╟─65e8b392-4a4c-4e09-a95e-a633a42e4bee
# ╠═e2f215e7-8d20-4e72-9e0b-cc0ce150e55b
# ╟─d061c4e2-8fbb-4ac5-95cf-c37d2916e669
# ╠═b462661a-3908-48e7-b162-f3b78d37231f
# ╟─bb1538ec-99e4-4ad0-8aaa-d552c6291c6d
# ╟─89c75cee-2488-4e9c-a470-b2111dcfc871
# ╠═8fbb4942-888c-4798-b972-7860e695d5ba
# ╟─5cfe9878-cbc2-4852-a757-ad8247d1f3d6
# ╟─ed8dce8b-d83b-40bc-b1eb-1ab369611cc1
# ╠═f855160a-2f04-4d11-97de-688113067c1c
# ╟─719af6d7-3c98-499e-b789-a6de5af2a027
# ╠═bcd1e77b-0869-4b36-9309-598b8c771d80
# ╟─e97ca371-9ee5-4419-afcd-2a35ed7426d1
# ╠═a0e70668-6c4e-4395-b1e4-925952a52dea
# ╟─764c7739-faad-4cf3-a64c-8ce3cae64074
# ╠═1971359b-df64-49e5-9011-c992da211c74
