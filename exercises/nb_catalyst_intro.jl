### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ e6e748a5-b692-4c7c-b39e-ef0fe9e46c35
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ be326550-25ea-4c5f-ac4a-dd74d12bc89a
using Markdown

# ╔═╡ cd166045-7788-4f5c-97d9-d1bedeceedae
using InteractiveUtils

# ╔═╡ b1bfe3b7-782c-46c9-86db-829942b6a9d0
using Catalyst

# ╔═╡ 19398992-dca7-441a-a5fd-7f9d0ad30be5
using DifferentialEquations, Plots

# ╔═╡ 9ea49d7a-e524-11ee-1b68-9d1d71aaba24
md"
# Introduction to Catalyst (ODE)
"

# ╔═╡ ba4dfc95-b74f-4d36-b34b-5eb2836a5cd6
md"
Catalyst.jl is a symbolic modeling package for analysis and high performance simulation of chemical reaction networks. Catalyst defines symbolic ReactionSystems, which can be created programmatically or easily specified using Catalyst's domain specific language (DSL).
"

# ╔═╡ 94440e99-7c0b-4d87-8038-141a2bc5fcb8
md"
This notebook describes the syntax for building chemical reaction network models using Catalyst's **D**omain-**S**pecific **L**anguage (DSL). We will illustrate this by implementing and solving an infection model by means of ODE (**O**rdinary **D**ifferential **E**quations).
"

# ╔═╡ 83fa478c-ac54-482b-92de-8eb6a01ddfaf
md"
## The infection model
"

# ╔═╡ d257102b-8480-4ce3-bdba-50995ccbdc26
md"
It is important to model the outbreak of infectious diseases in order to devise appropriate measures to avoid global epidemics. In this exercise we consider an isolated group of people in which a viral disease is spreading. An infection model (similar to the SIR-model but slightly extended) will be used for this purpose. We are interested in the evolution of the number of susceptible ($S$), infected ($I$), deceased ($D$) and resistant ($R$) persons.\
We make the following assumptions:
1. Transmission of the disease from an infected person to a susceptible person takes place through direct contact. The chance of any two inhabitants of the group coming into contact with each other is $\beta$, and the probability of infection after contact between an infected and a susceptible person is $\alpha$.
2. Note that the above assumption implicitly states that the probability of two neighbours coming into contact with each other is as high as the probability of two people living at two extremes of the territory coming into contact with each other.
3. A pereson leaves the infection period at a rate $r$ (hence, a person is contagious for an average of $1/r$ days. Without appropriate medication, a fraction $m$ of infected people die and a fraction $(1-m)$ of infected people acquire immunity after healing.
4. We assume that no one crosses the territory borders.
"

# ╔═╡ df9a46ab-56eb-4d3d-a3a8-e1f2da59b32e
md"
Below we summarize the **relevant variables** (**species**):
"

# ╔═╡ 7dd37224-071b-4ec7-88d5-39fd4482e615
md"
| Variable | Unit | Meaning |
|:---------- |:---------- |:------------|
| ``S``    | *persons* | number of susceptible persons             |
| ``I``    | *persons* | number of infected persons             |
| ``D``    | *persons* | number of deceased persons             |
| ``R``    | *persons* | number of recovered persons             |
"

# ╔═╡ 97d34ab6-a1c6-404d-af1f-7da85287757f
md"
Below we summarize the **parameters**:
"

# ╔═╡ b5e6cb25-da6f-475a-bc73-041c4c9256d5
md"
| Variable | Unit | Meaning |
|:---------- |:---------- |:------------|
| ``\alpha`` | ``\frac{persons}{contact}`` | chances of getting infected after contact |
| ``\beta`` | ``\frac{contact}{persons^2\,day}`` | contact rate |
| ``r`` | ``\frac{1}{day}`` | rate of leaving infection period |
| ``m`` | ``\frac{person}{person}`` | fraction of persons deceasing |
| ``1-m`` | ``\frac{person}{person}`` | fraction of persons becoming resistant |
"

# ╔═╡ 529ddac3-bb7a-4cc3-9bbd-ca31f5796e27
md"
Hence, the infection rate is ``\alpha \beta``. This means that a susceptible person meets an infected person: ``S+I``, this will result in ``2I`` at a rate ``\alpha \beta``. Futhermore, an infected person ``I`` will either become a deceased person ``D`` at a rate ``m r`` or become a resistant person ``R`` at rate ``(1-m) r``
"

# ╔═╡ f9ca0b86-1447-48c1-88e1-876180d7628a
md"
Our infection model has three reaction events:

- Infection, where a susceptible persons meets an infected persons and also becomes infected.
- Deceasing, where an infected person die.
- Recovery, where an infected person recovers.
"

# ╔═╡ bf2b06d8-884b-4a35-b308-83a2f3d6696f
md"
Each reaction is also associated with a specific rate:

- ``\alpha \beta``, the infection rate.
- ``m r``, the death rate.
- ``(1-m) r``, the recovery rate.
"

# ╔═╡ 1f405b09-53cd-43e0-b498-70eafef92acb
md"""
Hence, the following *infection reactions* are:

$$S + I \xrightarrow[]{\alpha \beta} 2I$$
$$I \xrightarrow[]{mr} D$$
$$I \xrightarrow[]{(1-m)r} R$$
"""

# ╔═╡ 39d80a3f-14d2-49d6-905b-171adf5930a0
md"
We are going to implement this system of *reactions* using Catalyst.
"

# ╔═╡ 44c6cdc0-8610-43bd-b3b3-a7441ee62615
md"
We first load the Catalyst package, which is required for the code in this introduction to run:
"

# ╔═╡ cd3da548-b24f-426f-8b46-2db70e0a979c
md"
#### Implementation of the system

The following code creates a so called *reaction network object*, that we have named `infection_model`, that implements the aforementioned *reactions*.
"

# ╔═╡ 7d2f8bf9-3a00-4631-84ff-1a36b6d7e19b
infection_model = @reaction_network begin
	α * β, S + I --> 2I
	r * m, I --> D
	r * (1 - m), I --> R
end

# ╔═╡ 4a504e36-8269-4822-8919-75b7f6dcaf6d
md"
Each line (between `begin` and `end`) corresponds to a *reaction*. Each *reaction* consists of:

- a **reaction rate** (the expression on the left-hand side of `,`),
- a set of **substrates** (the expression in-between `,` and `-->`),
- a set of **products** (the expression on the right-hand side of `-->`).

The substrates and the products may contain one or more reactants, separated by `+`.
"

# ╔═╡ af6b169d-55c3-4384-9f57-ad629e02bad6
md"
The *reaction model* is stored in the variable `infection_model` (the variable name can be chosen freely). It is a symbolic representation of the (chemical) network.
"

# ╔═╡ 52b2cc8f-1eac-47c1-8d42-09c1b91d1350
md"
You can get a list of the different *reaction* **species** with the command `species`
"

# ╔═╡ 5f1d6198-ccdc-4a63-bf98-2032a6d683ad
species(infection_model)

# ╔═╡ 59686253-3e0e-454e-8bc6-dde189ee2197
md"
To get a list of the *reaction* **parameters**, you can use the command `parameters`:
"

# ╔═╡ fdb79a26-fbc5-482a-8cdb-760d18f57a77
parameters(infection_model)

# ╔═╡ 64a48327-03ec-46a0-9a76-79c10acb6e20
md"
The *reaction model* can be converted to a symbolic differential equation model via
"

# ╔═╡ 78a22ba2-1e60-40b3-989c-922dcf9ca054
osys  = convert(ODESystem, infection_model)

# ╔═╡ eeb25e50-165b-4e93-8397-5b09fe8e7242
md"""
Note that the model equations are essencially:

$$\cfrac{dS(t)}{dt} = -\alpha \beta S(t) I(t)$$
$$\cfrac{dI(t)}{dt} = \alpha \beta S(t) I(t) - r I(t)$$
$$\cfrac{dD(t)}{dt} = m r I(t)$$
$$\cfrac{dR(t)}{dt} = (1-m) r I(t)$$
"""

# ╔═╡ 0d39fbbc-4fcb-45eb-b24b-14f2e093c98c
md"
You can get a list of the differential equations with the command `equations`:
"

# ╔═╡ ca18795a-d270-493f-807f-a3b8c78aa6d7
equations(osys)

# ╔═╡ 3f6f38a3-12ec-4b5b-8ebb-9e48190d3f7e
md"
To get a list of the state variables, you can use the command `states`:
"

# ╔═╡ 72f7f870-b2fc-4df1-9325-3c587be7e014
states(osys)

# ╔═╡ da5cfd6e-b57a-4235-b79d-a1c2f148328b
md"
To get a list of the parameters, you can use the command `parameters`:
"

# ╔═╡ a09b3a39-bedc-4180-ae5b-8cec708a64cb
parameters(osys)

# ╔═╡ 7fc8c671-b75c-4487-a013-779ee2422c8b
md"
#### Simulating the system as an ODE-problem

We first need to load the Differential and Plot package, which is required for simulating the system and plotting the results.
"

# ╔═╡ 3197244f-655b-4dca-80f3-794b30722551
md"
Now we wish to simulate our model. To do this, we need to provide some the following information:

- Initial conditions for the state variables $S$, $I$, $D$ and $R$.
- The parameter values for $\alpha$, $\beta$, $r$ and $m$.
- The timespan, which is the timeframe over which we wish to run the simulation.

Assume in this example that there are $10\,000\,000$ people in the country, and that initially $1\,000$ person are infected. Hence, $I_0 = 1\,000$, $S_0 = 10\,000\,000-I_0 = 9\,999\,000$, $D_0 = 0$ and $R_0 = 0$.\
Furthermore, we take the following values for the parameters: $\alpha = 0.08\;person/contact$, $\beta = 10^{-6}\;contact/(person^2\,day)$, $r = 0.2\;day^{-1}$ (i.e. a person is contagious for an average of $5\;days$) and $m=0.4$.\
Finally, we want to run our simulation from day $0$ till day $90$.
"

# ╔═╡ 979afe85-b910-44a0-8ac0-6e719cb9157e
md"
##### Setting initial conditions
"

# ╔═╡ 1ba859fa-46a5-434f-a99c-e710ba85caf8
md"
The initial conditions are given as a *Vector*. This is a type which collects several different values. To declare a vector, the values are specific within brackets, `[]`, and separated by `,`. Since we have four species, the vector holds four elements. E.g., we set the value of $I$ using the `:I => 1` syntax. Here, we first denote the name of the species (with a colon `:` pre-appended), next follows a `=>` and then the value of `I`.\
The vector holding the initial conditions for $S$, $I$, $D$ and $R$ can be created in the following way:
"

# ╔═╡ 35bd9a1a-bb4a-4285-98f0-853b03c95cb7
u0 = [:S => 9999000, :I => 1000, :D => 0, :R => 0]

# ╔═╡ 95c1c0ea-51d0-47ee-8948-a5b87c42a70d
md"
Note that the order of the vector elements doesn't matter, because the initial values of each of the species is indicated using its variable name.
"

# ╔═╡ 57036f49-2f1f-4327-89ed-2c96098a1c22
md"
##### Setting the timespan
"

# ╔═╡ 56ebb72d-2351-4e67-b268-1f48bbb77cb3
md"
The timespan sets the time point at which we start the simulation (typically `0.0` is used) and the final time point of the simulation. These are combined into a two-valued Tuple. Tuples are similar to vectors, but are enclosed by `()` and not `[]`. Again, we will let both time points be decimal valued.
"

# ╔═╡ a25d5925-a254-488b-b782-d29cff4470a2
tspan = (0.0, 90.0)

# ╔═╡ 0952b1d1-24b4-4540-91cd-94f7a4dcbd57
md"
##### Setting parameter values
"

# ╔═╡ a235c7ce-f14a-4c7c-86a1-08aa5f2d9c85
md"
Similarly, the parameters values are also given as a vector. We have four parameters, hence, the parameter vector will also contain four elements. We use a similar notation for setting the parameter values as the initial condition (first the colon, then the parameter name, then an arrow, then the value).
"

# ╔═╡ 9a9440fa-d8a3-44bc-8037-4bf1f8af40b0
params = [:α => 0.08, :β => 1.0e-6, :r => 0.2, :m => 0.4]

# ╔═╡ c4e83ef8-9490-4361-a2a9-5abc45e242be
md"
##### Creating an ODEProblem
"

# ╔═╡ aa1b904a-c8a9-41a4-9297-8d7c821d4b77
md"
Next, before we can simulate our model, we bundle all the required information together in a so-called **ODEProblem**. *Note that the order in which the input (the model name, the initial condition, the timespan, and the parameter values) is provided to the ODEProblem matters!* Here, we save our ODEProblem in the `oprob` variable.
"

# ╔═╡ c6d2dd69-8c61-4a40-894f-664b2d2d14be
oprob = ODEProblem(infection_model, u0, tspan, params)

# ╔═╡ 3c253bf3-886d-4e86-82ac-7751d23f342f
md"
##### Solving the ODEProblem
"

# ╔═╡ 14756171-9e8e-4cb0-b7af-74c2d649fe9f
md"
We can now simulate our model. We do this by providing the ODEProblem to the `solve` function. We save the output to the `sol` variable. Optionally, one can provide a [solver method](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/#Full-List-of-Methods) (e.g., `Tsit5`), and the time stepsize (with `saveat`).
"
# https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/#Full-List-of-Methods

# ╔═╡ 142e3e48-bf75-4498-ad0e-9f47cb921045
# osol = solve(oprob)
osol = solve(oprob, Tsit5(), saveat=0.5)

# ╔═╡ a3599781-a690-4fa3-b483-cd47727935cb
md"
Note that at the different time points the variables values in the solution are decimal numbers (and not integer numbers), despite the fact that we are applying the model to individuals. This is inherent to using an ODE approach. Later on, we will discretise the problem, and hence, work on the level of individual infections (reactions).\
Futhermore, note that executing the `solve` command at different occasions with an ODE problem will never modify the solution because ODE problems are **deterministic**. This will become different when simulating the individual infection (reaction) events by means of a stochastic (random) algorithm.
"

# ╔═╡ 0af3c166-46b6-455d-af6b-a72c4d2a5ce4
md"
Finally, we can plot the solution through the plot function.
"

# ╔═╡ 513c037b-c54c-47fa-b97a-06f69a983386
plot(osol)

# ╔═╡ b07f09c6-0515-4d27-9ca1-c45427a5988c
md"
If you want to see the final values of $S$, $I$, $D$ and $R$, type:
"

# ╔═╡ ab4940b8-b8ec-4835-8d6e-4e57a5e2e464
osol.u[end]

# ╔═╡ 764e2f1a-f974-4916-8573-cacba897cf07
md"
#### More advanced examples
"

# ╔═╡ 2ae76ddb-71f5-49d7-a250-429d6c0138f6
md"
In Example 1 we will show you one way of how you could analyze the simulation results for a limited range of parameter values.

In Examples 2 and 3 we will apply a new concept, namely **Callback**. The latter will basically affect, e.g., one or more parameter values or state variables during the solving process based in one or more *conditions* (also called *events*). These conditions can be either *time* or *state variable* related:

- A time related condition is a vector of one (or more) timepoint(s) for which the value of one (or more) parameter(s) or state variable(s) need to be altered. For this purpose we will use a so called `PresetTimeCallback`.
- A state variable related condition is usually a *condition function* that hits zero for a certain value of a state variable. For this purpose we will use a so called `ContinuousCallback`.
"

# ╔═╡ 590f1b49-7442-4a71-af8c-8acdea071448
md"
**Important remark:**\
You may have noticed that while using the Pluto notebooks, when you change the value of some variable (e.g., a parameter or an initial condition) that your results/plots will subsequently and automatically be altered based on the currect variable values in memory. In some cases this can be advantageous, in others not. For the latter reason, in this notebook, we will use slightly different variable names for some variables in order not to alter other results.
"

# ╔═╡ 45c1c238-a9f7-4f7b-a0ce-07b5bb4768d4
md"
##### Example 1 - Influence of $r$

Influence of the duration of infection $1/r$. We will check the effect of the average period in which infected persons are contagious if they are on average either $10$, $5$, $2$ days or $1$ day contagious. Hence, check the effect for $r=0.1, 0.2, 0.5$ and $1.0$.
"

# ╔═╡ ae38c663-0ee4-409e-bfca-5f13ed88b67d
md"
We will create een new parameter value vector, ODE problem and solution object by putting `1` at the end of the corresponding variable names. In that way, the previous simulation results will be unaffected! The model, the initial conditions and the timespan are identical as before.
"

# ╔═╡ d44da6c6-c93d-4c61-8125-9eee464c897e
params1 = [:α => 0.08, :β => 1.0e-6, :r => 0.1, :m => 0.4]

# ╔═╡ 00a72697-d36a-41cc-9eec-8e821829ce0e
# put semi-colon at end of instruction to avoid seeing its output.
oprob1 = ODEProblem(infection_model, u0, tspan, params1);

# ╔═╡ cf39b4cf-9cd0-4755-80db-ca4aea7c1084
osol1 = solve(oprob1, Tsit5(), saveat=0.5);

# ╔═╡ 52901bbf-e47b-4da1-95c4-f0869812398c
plot(osol1)

# ╔═╡ 102b4fbc-23b1-46ed-bb72-124eb88517ce
md"
Now, change the value of $r$ in the `param1` vector and analyze the effect in the plot.
"

# ╔═╡ f6cbafbf-98b4-4286-86d2-fe95821a5ff4
md"
Try to interpret the results yourself.
Ask yourself the following questions:

1. What are the trends in the results obtained?
- Answer: ...
2. How can this be explained from the model structure?
- Answer: ...
"

# ╔═╡ ade413c2-d7d9-4250-8490-75534900a389
md"
##### Example 2 - Using PresetTimeCallback

Suppose that regulations are such that on day 14, people need to reduce their contacts by 50%. Hence, this means that the parameter value $\beta$ needs to be divided by a factor of 2 at timepoint 14. In order to realize that we need to now the order of the parameters in the model because we will need to address the value of $\beta$ by means of an index.
"

# ╔═╡ e216061d-b1cc-43e9-b8b3-f6986041fcff
md"
First, check the order of the parameters in the model:
"

# ╔═╡ a173af17-6106-4065-9a8e-5d1954c77782
parameters(infection_model)

# ╔═╡ 50f51687-1ca2-4fc3-a0d2-35a67c09e597
md"
Hence, parameter $\beta$ has index `2` (2nd parameter).
"

# ╔═╡ dd62a738-64ef-4579-b20d-3443df87c985
md"
We will use the ODEProblem, stored in the variable `oprob`, that was previously created.
"

# ╔═╡ b5fe2703-6888-4519-8127-be872a8ffa76
md"
We now create the *condition*. We will store it in the variable `condition2` (from Example 2) in order not to interfere with the next example.
"

# ╔═╡ 5e1347e7-7268-43fd-98f0-7baae4199c60
condition2 = [14.0]

# ╔═╡ 5ad304c0-89ce-4cc6-a974-da938c30a396
md"
Next, we create a function that we will call `affect2!` (note the exclamation mark), that will be called by the solver at the timepoint(s) stored in `condition2` in order to alter the parameter value of $\beta$:
"

# ╔═╡ 182a46a9-5dce-4144-81cb-aa71a73c4ea0
function affect2!(integrator)
	integrator.p[2] = 0.5e-6     # β is the 2-nd parameter
end

# ╔═╡ 6dcd4258-d3b5-46a5-ad36-e4c7094e3fdb
md"
Next, we combine both `condition2` and `affect2!` with the function `PresetTimeCallback` in order to create the callback function, that we will name `cb2`:
"

# ╔═╡ 46a7e8aa-c07d-4e9c-9907-f72eb45aaecd
cb2 = PresetTimeCallback(condition2, affect2!)

# ╔═╡ 4dfc3e60-1475-47b9-a67a-e2d5aadf67a2
md"
Then, we solve the ODE problem, specifying the callback function `cb2`. Note, that we take a *deepcopy* of the `oprob` because otherwise the value of $\beta$ will remain altered in the original ODE problem!
"

# ╔═╡ b0864d84-f17f-4e58-ad63-a9610f8fc3dc
osol2 = solve(deepcopy(oprob), Tsit5(), saveat=0.5, callback=cb2)

# ╔═╡ 1e039147-7c48-4d5c-a432-0925eaa0872f
md"
Finally, we can plot the results:
"

# ╔═╡ bb2a80ee-81e7-441a-9fbc-fdb40cec6963
plot(osol2)

# ╔═╡ e4b03969-b316-40aa-a424-f7b95829ccaa
md"
If you want to see the final values of $S$, $I$, $D$ and $R$, type:
"

# ╔═╡ 5011ec50-f3d3-4b29-b1b2-d5a38a6af5c1
osol2.u[end]

# ╔═╡ 233a229c-56ed-4fab-ad3c-118b88972057
md"
Try to interpret the results yourself. Ask yourself the following questions:

1. What are the trends in the results obtained?
- Answer: ...
2. How much less casualties are there compared to not altering the contact rate?
- Answer: ...
"

# ╔═╡ cb8a6f77-f08c-4fc8-9445-bd1c17521fcc
md"
##### Example 3 - Using ContinuousCallback

Suppose that when the number of infected individuals reaches $1\,000\,000$, then $999\,000$ of them are promptly put into isolation (or removed from the population). Hence, a $1000$ individuals remain infected at some point. In order to realize that we need to now the order of the state variables in the model because we will need to address the value of $I$ by means of an index.\
Check the order of the state variables in the model with:
"

# ╔═╡ b000d4d1-c5f5-4471-81d4-b17d138cdabc
species(infection_model)

# ╔═╡ efce354e-e3d9-4bed-bb23-faadf68d80e7
md"
Hence, state variable $I$ has index `2` (2nd state variable).
"

# ╔═╡ 2206a68f-d3f0-484d-ae90-3c526cee1bc1
md"
The strategy will be to create a *condition function* with a condition that hits zero when $I$ equals $1\,000\,000$. Next when this is the case, we need to decrease the value of $I$ by $999\,000$. Hence, the condition could be coded as `u[2] - 1.0e6`, the latter number will become zero when $I$ (denoted as `u[2]`) reaches `1.0e6`.
"

# ╔═╡ 26983998-4dd2-4119-bd1d-531a759f9471
md"
First, we create a vector with the single number `1` in it. The purpose of this will become clear afterwards.
"

# ╔═╡ f3d4e76e-f74e-468d-8f41-15077c7031fc
proceed_with_condition = [1]

# ╔═╡ e84f85f2-25ac-4ba1-a246-495dae2d20a5
md"
Next, we create the *condition function*, that we will call `condition3` in the following way:
"

# ╔═╡ 48d2b79a-6cf8-4888-96f5-af7618747042
function condition3(u, t, integrator)
	u[2] - 1.0e6*proceed_with_condition[1]
end

# ╔═╡ e1dd4816-1cd2-4aac-ac1e-0dda5c2726aa
md"
Note that `proceed_with_condition[1]` refers to the first element (there is only one element) in the vector `proceed_with_condition`, which is `1` for now!
"

# ╔═╡ 35f83b24-3b0d-4e3a-95af-a22a5c31eab1
md"
Next, we create a function that we will call `affect3!` (note the exclamation mark), that will be called by the solver when the condition stored in `condition3` hits zero in order to alter the value of $I$:
"

# ╔═╡ 900f4d6a-9a32-406e-a63d-a0cd59bab2f3
function affect3!(integrator, param=proceed_with_condition)
	integrator.u[2] -= 0.999e6
	param[1] = 0
end

# ╔═╡ 196366f8-b3d6-4241-9fa6-dd253e929bd5
md"
Note that this function takes the vector `proceed_with_condition` as an argument (named `param`) and that the single element in `proceed_with_condition` is set to `0`. The purpose of that is to alter the condition in `condition3` in a way that it hits zero only once! If we don't take care of the condition in that way, then $I$ will be decreased by $999\,000$ **each time** when $I$ reaches $1\,000\,000$.
"

# ╔═╡ 1e325cf9-7495-4141-bd03-23349e96d666
md"
Next, we combine both `condition3` and `affect3!` with the function `ContinuousCallback` in order to create the callback function, that we will name `cb3`:
"

# ╔═╡ 8b6e3def-2b98-4b15-84fc-596bcf429d2f
cb3 = ContinuousCallback(condition3, affect3!)

# ╔═╡ d44ef0fe-585e-44cc-98be-cd5e8ba48c90
md"
Then, we solve the ODE problem, specifying the callback function `cb3`. Note, that again we play safe and take a *deepcopy* of the `oprob` to avoid alterations in the original ODE problem!
"

# ╔═╡ ec17e793-46c7-4337-9d13-056f62066cbe
osol3 = solve(deepcopy(oprob), Tsit5(), saveat=0.1, callback=cb3)

# ╔═╡ 2bc1768c-9525-4127-9ca3-334411d3abaf
md"
Finally, we can plot the results:
"

# ╔═╡ 08856000-3bc2-49de-9cf9-6f29f1c29097
plot(osol3)

# ╔═╡ 34fa1a17-4b44-4d41-9a05-5e72013538b7
md"
If you want to see the final values of $S$, $I$, $D$ and $R$, type:
"

# ╔═╡ b1e993a3-dade-43df-b9e9-fb132eb54d5c
osol3.u[end]

# ╔═╡ 47634eb2-4dc6-4fa7-9a17-df61eb4e22fa
md"
Try to interpret the results yourself. Ask yourself the following questions:

1. What are the trends in the results obtained?
- Answer: ...
2. How much less casualties are there compared to not putting $999\,000$ individuals into isolation? (Hint: you also need to take into account the casualties in the $999\,000$ individuals that had been put into isolation.)
- Answer: ...
"

# ╔═╡ Cell order:
# ╠═be326550-25ea-4c5f-ac4a-dd74d12bc89a
# ╠═cd166045-7788-4f5c-97d9-d1bedeceedae
# ╠═e6e748a5-b692-4c7c-b39e-ef0fe9e46c35
# ╠═9ea49d7a-e524-11ee-1b68-9d1d71aaba24
# ╠═ba4dfc95-b74f-4d36-b34b-5eb2836a5cd6
# ╠═94440e99-7c0b-4d87-8038-141a2bc5fcb8
# ╠═83fa478c-ac54-482b-92de-8eb6a01ddfaf
# ╠═d257102b-8480-4ce3-bdba-50995ccbdc26
# ╟─df9a46ab-56eb-4d3d-a3a8-e1f2da59b32e
# ╠═7dd37224-071b-4ec7-88d5-39fd4482e615
# ╟─97d34ab6-a1c6-404d-af1f-7da85287757f
# ╠═b5e6cb25-da6f-475a-bc73-041c4c9256d5
# ╠═529ddac3-bb7a-4cc3-9bbd-ca31f5796e27
# ╠═f9ca0b86-1447-48c1-88e1-876180d7628a
# ╠═bf2b06d8-884b-4a35-b308-83a2f3d6696f
# ╠═1f405b09-53cd-43e0-b498-70eafef92acb
# ╠═39d80a3f-14d2-49d6-905b-171adf5930a0
# ╠═44c6cdc0-8610-43bd-b3b3-a7441ee62615
# ╠═b1bfe3b7-782c-46c9-86db-829942b6a9d0
# ╠═cd3da548-b24f-426f-8b46-2db70e0a979c
# ╠═7d2f8bf9-3a00-4631-84ff-1a36b6d7e19b
# ╟─4a504e36-8269-4822-8919-75b7f6dcaf6d
# ╟─af6b169d-55c3-4384-9f57-ad629e02bad6
# ╠═52b2cc8f-1eac-47c1-8d42-09c1b91d1350
# ╠═5f1d6198-ccdc-4a63-bf98-2032a6d683ad
# ╠═59686253-3e0e-454e-8bc6-dde189ee2197
# ╠═fdb79a26-fbc5-482a-8cdb-760d18f57a77
# ╠═64a48327-03ec-46a0-9a76-79c10acb6e20
# ╠═78a22ba2-1e60-40b3-989c-922dcf9ca054
# ╟─eeb25e50-165b-4e93-8397-5b09fe8e7242
# ╠═0d39fbbc-4fcb-45eb-b24b-14f2e093c98c
# ╠═ca18795a-d270-493f-807f-a3b8c78aa6d7
# ╠═3f6f38a3-12ec-4b5b-8ebb-9e48190d3f7e
# ╠═72f7f870-b2fc-4df1-9325-3c587be7e014
# ╠═da5cfd6e-b57a-4235-b79d-a1c2f148328b
# ╠═a09b3a39-bedc-4180-ae5b-8cec708a64cb
# ╠═7fc8c671-b75c-4487-a013-779ee2422c8b
# ╠═19398992-dca7-441a-a5fd-7f9d0ad30be5
# ╠═3197244f-655b-4dca-80f3-794b30722551
# ╠═979afe85-b910-44a0-8ac0-6e719cb9157e
# ╠═1ba859fa-46a5-434f-a99c-e710ba85caf8
# ╠═35bd9a1a-bb4a-4285-98f0-853b03c95cb7
# ╟─95c1c0ea-51d0-47ee-8948-a5b87c42a70d
# ╠═57036f49-2f1f-4327-89ed-2c96098a1c22
# ╠═56ebb72d-2351-4e67-b268-1f48bbb77cb3
# ╠═a25d5925-a254-488b-b782-d29cff4470a2
# ╠═0952b1d1-24b4-4540-91cd-94f7a4dcbd57
# ╠═a235c7ce-f14a-4c7c-86a1-08aa5f2d9c85
# ╠═9a9440fa-d8a3-44bc-8037-4bf1f8af40b0
# ╠═c4e83ef8-9490-4361-a2a9-5abc45e242be
# ╟─aa1b904a-c8a9-41a4-9297-8d7c821d4b77
# ╠═c6d2dd69-8c61-4a40-894f-664b2d2d14be
# ╠═3c253bf3-886d-4e86-82ac-7751d23f342f
# ╠═14756171-9e8e-4cb0-b7af-74c2d649fe9f
# ╠═142e3e48-bf75-4498-ad0e-9f47cb921045
# ╠═a3599781-a690-4fa3-b483-cd47727935cb
# ╠═0af3c166-46b6-455d-af6b-a72c4d2a5ce4
# ╠═513c037b-c54c-47fa-b97a-06f69a983386
# ╠═b07f09c6-0515-4d27-9ca1-c45427a5988c
# ╠═ab4940b8-b8ec-4835-8d6e-4e57a5e2e464
# ╟─764e2f1a-f974-4916-8573-cacba897cf07
# ╠═2ae76ddb-71f5-49d7-a250-429d6c0138f6
# ╠═590f1b49-7442-4a71-af8c-8acdea071448
# ╠═45c1c238-a9f7-4f7b-a0ce-07b5bb4768d4
# ╠═ae38c663-0ee4-409e-bfca-5f13ed88b67d
# ╠═d44da6c6-c93d-4c61-8125-9eee464c897e
# ╠═00a72697-d36a-41cc-9eec-8e821829ce0e
# ╠═cf39b4cf-9cd0-4755-80db-ca4aea7c1084
# ╠═52901bbf-e47b-4da1-95c4-f0869812398c
# ╠═102b4fbc-23b1-46ed-bb72-124eb88517ce
# ╟─f6cbafbf-98b4-4286-86d2-fe95821a5ff4
# ╠═ade413c2-d7d9-4250-8490-75534900a389
# ╠═e216061d-b1cc-43e9-b8b3-f6986041fcff
# ╠═a173af17-6106-4065-9a8e-5d1954c77782
# ╠═50f51687-1ca2-4fc3-a0d2-35a67c09e597
# ╠═dd62a738-64ef-4579-b20d-3443df87c985
# ╠═b5fe2703-6888-4519-8127-be872a8ffa76
# ╠═5e1347e7-7268-43fd-98f0-7baae4199c60
# ╠═5ad304c0-89ce-4cc6-a974-da938c30a396
# ╠═182a46a9-5dce-4144-81cb-aa71a73c4ea0
# ╠═6dcd4258-d3b5-46a5-ad36-e4c7094e3fdb
# ╠═46a7e8aa-c07d-4e9c-9907-f72eb45aaecd
# ╠═4dfc3e60-1475-47b9-a67a-e2d5aadf67a2
# ╠═b0864d84-f17f-4e58-ad63-a9610f8fc3dc
# ╠═1e039147-7c48-4d5c-a432-0925eaa0872f
# ╠═bb2a80ee-81e7-441a-9fbc-fdb40cec6963
# ╠═e4b03969-b316-40aa-a424-f7b95829ccaa
# ╠═5011ec50-f3d3-4b29-b1b2-d5a38a6af5c1
# ╠═233a229c-56ed-4fab-ad3c-118b88972057
# ╠═cb8a6f77-f08c-4fc8-9445-bd1c17521fcc
# ╠═b000d4d1-c5f5-4471-81d4-b17d138cdabc
# ╠═efce354e-e3d9-4bed-bb23-faadf68d80e7
# ╠═2206a68f-d3f0-484d-ae90-3c526cee1bc1
# ╠═26983998-4dd2-4119-bd1d-531a759f9471
# ╠═f3d4e76e-f74e-468d-8f41-15077c7031fc
# ╠═e84f85f2-25ac-4ba1-a246-495dae2d20a5
# ╠═48d2b79a-6cf8-4888-96f5-af7618747042
# ╠═e1dd4816-1cd2-4aac-ac1e-0dda5c2726aa
# ╠═35f83b24-3b0d-4e3a-95af-a22a5c31eab1
# ╠═900f4d6a-9a32-406e-a63d-a0cd59bab2f3
# ╠═196366f8-b3d6-4241-9fa6-dd253e929bd5
# ╠═1e325cf9-7495-4141-bd03-23349e96d666
# ╠═8b6e3def-2b98-4b15-84fc-596bcf429d2f
# ╠═d44ef0fe-585e-44cc-98be-cd5e8ba48c90
# ╠═ec17e793-46c7-4337-9d13-056f62066cbe
# ╠═2bc1768c-9525-4127-9ca3-334411d3abaf
# ╠═08856000-3bc2-49de-9cf9-6f29f1c29097
# ╠═34fa1a17-4b44-4d41-9a05-5e72013538b7
# ╠═b1e993a3-dade-43df-b9e9-fb132eb54d5c
# ╠═47634eb2-4dc6-4fa7-9a17-df61eb4e22fa
