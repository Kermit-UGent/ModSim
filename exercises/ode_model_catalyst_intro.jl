### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

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

# ╔═╡ 3cfe381f-5871-4447-9149-fc023c787b42
using PlutoUI; TableOfContents()

# ╔═╡ b1bfe3b7-782c-46c9-86db-829942b6a9d0
using Catalyst

# ╔═╡ 19398992-dca7-441a-a5fd-7f9d0ad30be5
using DifferentialEquations, Plots

# ╔═╡ 9ea49d7a-e524-11ee-1b68-9d1d71aaba24
md"""
# Introduction to Catalyst (ODE)
"""

# ╔═╡ ba4dfc95-b74f-4d36-b34b-5eb2836a5cd6
md"""
Catalyst.jl is a symbolic modelling package for construction, analysis and high performance simulation of chemical reaction networks. Catalyst defines symbolic ReactionSystems, which can be created programmatically or easily specified using Catalyst's **D**omain **S**pecific **L**anguage (DSL).
"""

# ╔═╡ 94440e99-7c0b-4d87-8038-141a2bc5fcb8
md"""
This notebook describes the syntax for building chemical reaction network models using Catalyst's DSL. We will illustrate this by implementing and solving an infection model by means of ODE (**O**rdinary **D**ifferential **E**quations).
"""

# ╔═╡ 83fa478c-ac54-482b-92de-8eb6a01ddfaf
md"""
## The infection model
"""

# ╔═╡ d257102b-8480-4ce3-bdba-50995ccbdc26
md"""
It is important to model the outbreak of infectious diseases in order to devise appropriate measures to avoid global epidemics. In this exercise, we consider an isolated group of people in which a viral disease is spreading. We use an infection model (similar to the SIR-model but slightly extended) for this purpose. We are interested in the evolution of the number of susceptible ($S$), infected ($I$), deceased ($D$) and resistant ($R$) persons.\
We make the following assumptions:
1. Transmission of the disease from an infected person to a susceptible person takes place through direct contact. The chance of any two inhabitants of the group coming into contact with each other is $\beta$, and the probability of infection after contact between an infected and a susceptible person is $\alpha$.
2. Note that the above assumption implicitly states that the probability of two neighbours coming into contact with each other is as high as the probability of two people living at two extremes of the territory coming into contact with each other.
3. A person leaves the infection period at a rate $r$ (hence, a person is contagious for an average of $1/r$ days. Without appropriate medication, a fraction $m$ of infected people die and a fraction $(1-m)$ of infected people acquire immunity after healing.
4. We assume that there is no migration in or out of the population.
"""

# ╔═╡ df9a46ab-56eb-4d3d-a3a8-e1f2da59b32e
md"
Below we summarize the **relevant variables** (**species**):
"

# ╔═╡ 7dd37224-071b-4ec7-88d5-39fd4482e615
md"""
| Variable | Unit | Meaning |
|:---------- |:---------- |:------------|
| ``S``    | *persons* | number of susceptible persons             |
| ``I``    | *persons* | number of infected persons             |
| ``D``    | *persons* | number of deceased persons             |
| ``R``    | *persons* | number of resistant persons             |
"""

# ╔═╡ 97d34ab6-a1c6-404d-af1f-7da85287757f
md"
Below we summarize the **parameters**:
"

# ╔═╡ b5e6cb25-da6f-475a-bc73-041c4c9256d5
md"""
| Variable | Unit | Meaning |
|:---------- |:---------- |:------------|
| ``\alpha`` | ``\frac{persons}{contact}`` | chances of getting infected after contact |
| ``\beta`` | ``\frac{contact}{persons^2\,day}`` | contact rate |
| ``r`` | ``\frac{1}{day}`` | rate of leaving infection period |
| ``m`` | ``\frac{person}{person}`` | fraction of persons deceasing |
| ``1-m`` | ``\frac{person}{person}`` | fraction of persons becoming resistant |
"""

# ╔═╡ 529ddac3-bb7a-4cc3-9bbd-ca31f5796e27
md"""
Hence, the infection rate is ``\alpha \beta``. This means that a susceptible person meets an infected person: ``S+I``, this will result in ``2I`` at a rate ``\alpha \beta``. Futhermore, an infected person ``I`` will either become a deceased person ``D`` at a rate ``m r`` or become a resistant person ``R`` at rate ``(1-m) r``
"""

# ╔═╡ f9ca0b86-1447-48c1-88e1-876180d7628a
md"""
Our infection model has three reaction events:

1. Infection, where a susceptible persons meets an infected persons and also becomes infected.
2. Deceasing, where an infected person die.
3. Recovery, where an infected person recovers and becomes resitant.
"""

# ╔═╡ bf2b06d8-884b-4a35-b308-83a2f3d6696f
md"""
Each reaction is also associated with a specific rate:

1. ``\alpha \beta``, the infection rate.
2. ``m r``, the death rate.
3. ``(1-m) r``, the recovery rate.
"""

# ╔═╡ 1f405b09-53cd-43e0-b498-70eafef92acb
md"""
Hence, the following *infection reactions* are:

$$S + I \xrightarrow[]{\alpha \beta} 2I$$
$$I \xrightarrow[]{mr} D$$
$$I \xrightarrow[]{(1-m)r} R$$
"""

# ╔═╡ 39d80a3f-14d2-49d6-905b-171adf5930a0
md"""
We are going to implement this system of *reactions* using Catalyst.
"""

# ╔═╡ 44c6cdc0-8610-43bd-b3b3-a7441ee62615
md"""
We first load the Catalyst package, which is required for the code in this introduction to run:
"""

# ╔═╡ cd3da548-b24f-426f-8b46-2db70e0a979c
md"""
### Implementation of the system

The following code creates a so called *reaction network object*, that we have named `infection_model`, that implements the aforementioned *reactions*.
"""

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

# ╔═╡ bb5578b8-40d4-48c9-be87-cfaeadb8c5f9
md"""
!!! tip "Tip"
	The Greek letters can be visualized by typing a **backslash** followed by the **name of the Greek letter** and then the **TAB** key. For example `\alpha` followed by the TAB key results in in a list where you can choose `α`.
"""

# ╔═╡ af6b169d-55c3-4384-9f57-ad629e02bad6
md"
The *reaction model* is stored in the variable `infection_model` (the variable name can be chosen freely). It is a symbolic representation of the (chemical) network.
"

# ╔═╡ 52b2cc8f-1eac-47c1-8d42-09c1b91d1350
md"""
You can get a list of the different *reaction* **species** with the command `species`
"""

# ╔═╡ 5f1d6198-ccdc-4a63-bf98-2032a6d683ad
species(infection_model)

# ╔═╡ 59686253-3e0e-454e-8bc6-dde189ee2197
md"""
To get a list of the *reaction* **parameters**, you can use the command `parameters`:
"""

# ╔═╡ fdb79a26-fbc5-482a-8cdb-760d18f57a77
parameters(infection_model)

# ╔═╡ 2bd02577-7321-46b5-9f42-1484ed86d1d3
md"""
!!! important "Important"
	You can also get the different species and parameters using `@unpack` followed by comma separated species and/or parameter names followed by the equal sign and the name of the reaction network model. For example:
"""

# ╔═╡ 64455b26-4374-435f-b71b-8b9600cab263
@unpack S, I, D, R = infection_model;

# ╔═╡ 64a48327-03ec-46a0-9a76-79c10acb6e20
md"""
The *reaction model* can be converted to a symbolic differential equation model via
"""

# ╔═╡ 78a22ba2-1e60-40b3-989c-922dcf9ca054
osys = convert(ODESystem, infection_model)

# ╔═╡ eeb25e50-165b-4e93-8397-5b09fe8e7242
md"""
Note that the model equations are essentially:

$$\cfrac{dS(t)}{dt} = -\alpha \beta S(t) I(t)$$
$$\cfrac{dI(t)}{dt} = \alpha \beta S(t) I(t) - r I(t)$$
$$\cfrac{dD(t)}{dt} = m r I(t)$$
$$\cfrac{dR(t)}{dt} = (1-m) r I(t)$$
"""

# ╔═╡ 0d39fbbc-4fcb-45eb-b24b-14f2e093c98c
md"""
You can get a list of the differential equations with the command `equations`:
"""

# ╔═╡ ca18795a-d270-493f-807f-a3b8c78aa6d7
equations(osys)

# ╔═╡ 3f6f38a3-12ec-4b5b-8ebb-9e48190d3f7e
md"""
To get a list of the state variables, you can use the command `unknowns`:
"""

# ╔═╡ 72f7f870-b2fc-4df1-9325-3c587be7e014
unknowns(osys)

# ╔═╡ da5cfd6e-b57a-4235-b79d-a1c2f148328b
md"""
To get a list of the parameters, you can use the command `parameters`:
"""

# ╔═╡ a09b3a39-bedc-4180-ae5b-8cec708a64cb
parameters(osys)

# ╔═╡ 7fc8c671-b75c-4487-a013-779ee2422c8b
md"""
### Simulating the system as an ODE-problem

We first need to load the `DifferentialEquations` and `Plots` packages, which are required for simulating the system and plotting the results.
"""

# ╔═╡ 3197244f-655b-4dca-80f3-794b30722551
md"""
Now we wish to simulate our model. To do this, we need to provide the following information:

- Initial conditions for the state variables $S$, $I$, $D$ and $R$.
- The parameter values for $\alpha$, $\beta$, $r$ and $m$.
- The timespan, which is the timeframe over which we wish to run the simulation.

Assume in this example that there are $10\,000\,000$ people in the country, and that initially $1\,000$ persons are infected. Hence, $I_0 = 1\,000$, $S_0 = 10\,000\,000-I_0 = 9\,999\,000$, $D_0 = 0$ and $R_0 = 0$.\
Furthermore, we take the following values for the parameters: $\alpha = 0.08\;person/contact$, $\beta = 10^{-6}\;contact/(person^2\,day)$, $r = 0.2\;day^{-1}$ (i.e. a person is contagious for an average of $5\;days$) and $m=0.4$. The following table summarizes the above values:

|Initial conditions                   |Parameters          |
|:------------------------------------|:-------------------|
|$S_0 = 9\,999\,000$                  |$\alpha = 0.08$     |
|$I_0 = 1\,000$                       |$\beta = 10^{-6}$   |
|$D_0 = 0$                            |$r = 0.2$           |
|$R_0 = 0$                            |$m=0.4$             |

Finally, we want to run our simulation from day $0$ till day $90$.
"""

# ╔═╡ 979afe85-b910-44a0-8ac0-6e719cb9157e
md"""
### Setting initial conditions
"""

# ╔═╡ 1ba859fa-46a5-434f-a99c-e710ba85caf8
md"""
The initial conditions are given as a *Vector*. This is a type which collects several different values. To declare a vector, the values are specified within brackets, `[]`, and separated by `,`. Since we have four species, the vector holds four elements. E.g., we set the value of $I$ using the `:I => 1` syntax. Here, we first denote the name of the species (with a colon `:` pre-appended), next follows a `=>` and then the value of `I`.\
The vector holding the initial conditions for $S$, $I$, $D$ and $R$ can be created in the following way:
"""

# ╔═╡ 35bd9a1a-bb4a-4285-98f0-853b03c95cb7
u0 = [:S => 9_999_000.0, :I => 1_000.0, :D => 0.0, :R => 0.0]

# ╔═╡ 95c1c0ea-51d0-47ee-8948-a5b87c42a70d
md"
Note that the order of the vector elements doesn't matter here, since the initial values of each of the species is indicated using its variable name.
"

# ╔═╡ 57036f49-2f1f-4327-89ed-2c96098a1c22
md"""
### Setting the timespan
"""

# ╔═╡ 56ebb72d-2351-4e67-b268-1f48bbb77cb3
md"""
The timespan sets the time point at which we start the simulation (typically `0.0` is used) and the final time point of the simulation. These are combined into a two-valued *Tuple*. Tuples are similar to vectors, but are enclosed by `()` and not `[]`. Again, we will let both time points be decimal valued.
"""

# ╔═╡ a25d5925-a254-488b-b782-d29cff4470a2
tspan = (0.0, 90.0)

# ╔═╡ 0952b1d1-24b4-4540-91cd-94f7a4dcbd57
md"""
### Setting parameter values
"""

# ╔═╡ a235c7ce-f14a-4c7c-86a1-08aa5f2d9c85
md"""
Similarly, the parameter values are also given as a vector. We have four parameters, hence, the parameter vector will also contain four elements. We use a similar notation for setting the parameter values as the initial condition (first the colon, then the parameter name, then an arrow, then the value).
"""

# ╔═╡ 9a9440fa-d8a3-44bc-8037-4bf1f8af40b0
params = [:α => 0.08, :β => 1.0e-6, :r => 0.2, :m => 0.4]

# ╔═╡ c4e83ef8-9490-4361-a2a9-5abc45e242be
md"""
### Creating an ODEProblem
"""

# ╔═╡ aa1b904a-c8a9-41a4-9297-8d7c821d4b77
md"
Next, before we can simulate our model, we bundle all the required information together in a so-called **ODEProblem**. *Note that the order in which the input (the model name, the initial condition, the timespan, and the parameter values) is provided to the ODEProblem matters!* Here, we save our ODEProblem in the `oprob` variable.
"

# ╔═╡ c6d2dd69-8c61-4a40-894f-664b2d2d14be
oprob = ODEProblem(infection_model, u0, tspan, params)

# ╔═╡ 3c253bf3-886d-4e86-82ac-7751d23f342f
md"""
### Solving the ODEProblem
"""

# ╔═╡ 14756171-9e8e-4cb0-b7af-74c2d649fe9f
md"""
We can now simulate our model. We do this by providing the ODEProblem to the `solve` function. There are some [examples](https://docs.sciml.ai/DiffEqDocs/stable/getting_started/) online on how to solve ODE problems with the *DifferentialEquations.jl* package. We save the output to the `sol` variable. Optionally, one can provide a [solver method](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/#Full-List-of-Methods) (e.g., `Tsit5`), and the time stepsize (with `saveat`).
"""
# https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/#Full-List-of-Methods

# ╔═╡ 142e3e48-bf75-4498-ad0e-9f47cb921045
# osol = solve(oprob)
osol = solve(oprob, Tsit5(), saveat=0.5)

# ╔═╡ a3599781-a690-4fa3-b483-cd47727935cb
md"""
Note that at the different time points the variables values in the solution are decimal numbers (and not integer numbers), despite the fact that we are applying the model to individuals. This is inherent to using an ODE approach. Later on, we will see how we can discretise the problem, and hence, work on the level of individual infections (reactions).\
Futhermore, note that executing the `solve` command at different occasions with an ODE problem will never modify the solution because ODE problems are **deterministic**. This will become different when simulating the individual infection (reaction) events by means of a stochastic (random) algorithm.
"""

# ╔═╡ 8ac90e5c-ce85-4113-bc1b-9dc17cf7e928
md"""
### Plotting the results
"""

# ╔═╡ 0af3c166-46b6-455d-af6b-a72c4d2a5ce4
md"""
Finally, we can plot the solution through the `plot` function.
"""

# ╔═╡ 513c037b-c54c-47fa-b97a-06f69a983386
plot(osol)

# ╔═╡ 5672749b-283a-4e8e-8f90-76b5312b30ac
md"""
If you want to plot less species, like for example just $S$ and $I$, you can specify this with the option `idxs=[:S, :I]` (*notice the brackets*) in the plot function.
"""

# ╔═╡ 01188934-ce0f-4d65-b880-caac4837d796
plot(osol, idxs=[:S, :I])     # brackets [ ]

# ╔═╡ f3cf01a1-3c65-4a37-bf9a-cd2233cb470b
md"""
If you want a phase plot of, for example, just $I$ versus $S$, you can specify this with the option `idxs=(:S, :I)` (*notice the parentheses*) in the plot function. You can indicate the $S$ and $I$ axes with the additional options `xlab="S"` and `ylab="I"`.
"""

# ╔═╡ ff8f4c23-5695-48aa-9965-02677103f2c9
plot(osol, idxs=(:S, :I), xlab="S", ylab="I")     # parentheses ( )

# ╔═╡ b07f09c6-0515-4d27-9ca1-c45427a5988c
md"""
If you want to see the final values of $S$, $I$, $D$ and $R$, type:
"""

# ╔═╡ ab4940b8-b8ec-4835-8d6e-4e57a5e2e464
osol.u[end]

# ╔═╡ a8b68546-9401-4815-9f33-cc22b7a0579d
md"""
If you want the vector of, e.g., $S$ values separately, type:
"""

# ╔═╡ 89cc8970-693d-4885-b0b5-46110db4ca32
osol[:S]

# ╔═╡ dcdd693d-f2f3-406c-8647-31cdedcb823b
md"""
If you want the last value in the S vector, type:
"""

# ╔═╡ 0a2ed6a6-1535-4721-a6da-008f62b16942
osol[:S][end]

# ╔═╡ 757d9f1e-5cba-4ae8-9172-e0c2fceb8458
md"""
If you want the time vector separately, type:
"""

# ╔═╡ 54bedb4d-51c5-4f58-bc4b-50cc469b7b04
osol.t

# ╔═╡ 764e2f1a-f974-4916-8573-cacba897cf07
md"
## More advanced examples
"

# ╔═╡ 2ae76ddb-71f5-49d7-a250-429d6c0138f6
md"""
In Example 1 we will show you one way of how you could analyze the simulation results for a limited range of parameter values.

In Examples 2 and 3 we will apply some new concepts, namely **discrete** and **continuous** events. The latter will basically affect, e.g., one or more parameter values or state variables during the solving process based on one or more *conditions* (also called *events*). These conditions can be either *time* or *state variable* related:

- A time related condition is a vector of one (or more) timepoint(s) for which the value of one (or more) parameter(s) or state variable(s) need to be altered. We refer to them as **discrete events**.
- A state variable related condition is usually a condition for a certain value of a state variable. We will refer to them as **continuous events**.
"""

# ╔═╡ 590f1b49-7442-4a71-af8c-8acdea071448
md"""
**Important remark:**\
You may have noticed that while using the Pluto notebooks, when you change the value of some variable (e.g., a parameter or an initial condition) that your results/plots will subsequently and automatically be altered based on the currect variable values in memory. In some cases this can be advantageous, in others not. For the latter reason, in this notebook, we will use slightly different variable names for some variables in order not to alter other results.
"""

# ╔═╡ 45c1c238-a9f7-4f7b-a0ce-07b5bb4768d4
md"""
### Example 1 - Influence of $r$

Influence of the duration of infection $1/r$ for average infection periods of between $10$ days and $1$ day of being contagious ($r$ between $0.1$ and $1.0$, step $0.1$, default value $0.1$).
"""

# ╔═╡ b6baafc2-6d5e-43c3-8ef9-845961cdd20b
md"""
We will create a slider for the $r$-values between $0.1$ and $1.0$, stepsize $0.1$, default value $0.1$.
"""

# ╔═╡ cd32beba-67cf-4b12-a77b-99f96263f0a4
@bind r Slider(0.1:0.1:1, default=0.1, show_value=true)

# ╔═╡ ae38c663-0ee4-409e-bfca-5f13ed88b67d
md"""
We will create a new parameter value vector, ODE problem and solution object by putting `1` at the end of the corresponding variable names. In that way, the previous simulation results will be unaffected! The model, the initial conditions and the timespan are identical as before. In there we also use the variable `r` coupled to the slider.
"""

# ╔═╡ d44da6c6-c93d-4c61-8125-9eee464c897e
params1 = [:α => 0.08, :β => 1.0e-6, :r => r, :m => 0.4]

# ╔═╡ 00a72697-d36a-41cc-9eec-8e821829ce0e
# type a semi-colon at end of an instruction to avoid seeing its output
oprob1 = ODEProblem(infection_model, u0, tspan, params1);

# ╔═╡ cf39b4cf-9cd0-4755-80db-ca4aea7c1084
osol1 = solve(oprob1, Tsit5(), saveat=0.5);

# ╔═╡ 52901bbf-e47b-4da1-95c4-f0869812398c
plot(osol1, ylim=(0, 1e7))

# ╔═╡ 102b4fbc-23b1-46ed-bb72-124eb88517ce
md"""
Now, change the value of $r$ in the `param1` vector and analyze the effect in the plot.
"""

# ╔═╡ f6cbafbf-98b4-4286-86d2-fe95821a5ff4
md"""
Try to interpret the results yourself.
Ask yourself the following questions:

!!! question
	1. What are the trends in the results obtained?
"""

# ╔═╡ 4242ec68-b36d-4f7f-9001-2c90c66d5b9b
md"""
- Answer: 
"""

# ╔═╡ f5f98168-119c-4b49-a8e3-a3cc775faeb0
md"""
!!! question
	2. How can this be explained from the model structure?
"""

# ╔═╡ 30b0e2b3-2b8d-43a5-9587-9fcb0b89b258
md"""
- Answer: 
"""

# ╔═╡ ade413c2-d7d9-4250-8490-75534900a389
md"""
### Example 2 - Discrete Event

Suppose that regulations are such that on day 14, people need to reduce their contacts by 50%. Hence, this means that the parameter value $\beta$ needs to be divided by a factor of 2 at timepoint 14. 
"""

# ╔═╡ 58730ac6-d83b-420d-a000-2f50545f0d39
md"""
We need to state that the parameter $\beta$ needs to be reduced by $50\%$ at time $t=14$. We include this condition in a variable named `condition2` in the following way:
"""

# ╔═╡ 7411474c-fac7-4b7c-8ded-4c2df5956fb0
condition2 = [14.0] => [infection_model.β ~ infection_model.β/2]

# ╔═╡ e6949052-6d12-4414-ac67-cb9435b46290
md"""
The discrete time event needs to be now included in our model.
"""

# ╔═╡ 8752373f-a602-437b-9b3a-2602c8babf87
@named infection_model2 = ReactionSystem(equations(infection_model), discrete_events=condition2)

# ╔═╡ ee596a6c-29bf-4d18-b687-be943c128aa7
md"""
After that, we need to *complete* our *reaction network model*, so that the model can be simulated.
"""

# ╔═╡ 40f849cb-20f9-4bbe-806e-512abd6f3210
infection_model2_com = complete(infection_model2)

# ╔═╡ d4e537da-490b-4781-9a76-c89d762849ec
md"""
Then we need to create a new ODE problem.
"""

# ╔═╡ 4a6f357c-afd7-4b0d-b1d2-15f5c0e4298b
oprob2 = ODEProblem(infection_model2_com, u0, tspan, params)

# ╔═╡ e76e77be-9f15-4aa5-8506-cdb32a6ec9b1
md"""
Finally, the ODE problem can be solved. Notice that you need to make a deepcopy of the ODE problem, because otherwise changes to the parameter $\beta$ will remain after the first call to `solve`.
"""

# ╔═╡ 37b11fdb-8c16-4ec4-a7ad-e2fdbaf4ea5d
osol2 = solve(deepcopy(oprob2), Tsit5(), saveat=0.5)

# ╔═╡ 59f4f4dd-f0b5-4354-b53c-4554c878d5a8
md"""
Now we can plot the results.
"""

# ╔═╡ 0af84963-ca64-4958-a544-42d445da5a7c
plot(osol2)

# ╔═╡ 18151ab1-1e4e-48d8-be70-4fa4c8a43af1
md"""
If you want to see the final values of $S$, $I$, $D$ and $R$, type:
"""

# ╔═╡ 7d7ed974-7c9c-4c30-ab2a-e3028ce702dd
osol2.u[end]

# ╔═╡ 7375edbc-24a4-4300-bdef-2686cb377cfc
md"""
Try to interpret the results yourself. Ask yourself the following questions:

!!! question
	1. What are the trends in the results obtained?
"""

# ╔═╡ 6eb368f5-8d43-4fb8-b6eb-c783675d1dd9
md"""
- Answer: 
"""

# ╔═╡ 2ba54805-e901-4281-9b07-cc7137b8c809
md"""
!!! question
	2. How much less casualties are there compared to not altering the contact rate?
"""

# ╔═╡ fec6e165-bd1a-4e57-b7f3-ae712e03276c
md"""
- Answer: 
"""

# ╔═╡ cb8a6f77-f08c-4fc8-9445-bd1c17521fcc
md"""
### Example 3 - Continuous Event

Suppose that when the number of infected individuals reaches $1\,000\,000$, then $999\,000$ of them are promptly put into isolation (or removed from the population). Hence, a $1000$ individuals remain infected at some point.
"""

# ╔═╡ 5005a09d-a844-4f9d-a058-0d93583d5bab
md"""
Normally in a continuous event the value of one or more species can be changed when a certain condition is met. In our specific case we want the change in the species happening only once! So, if you want that the continuous event "when $I$ reaches $10^6$ then $999000$ is subtracted from $I$" happens only once, then we need to include a ficticious new *species* in our *reaction network model*. We will call this ficticious species `pwc` (a short for _**p**roceed **w**ith **c**ondition_) and we set its default value to `true`.
"""

# ╔═╡ f77d2d75-468f-4746-b81e-0bbbf33fc8d7
infection_model3 = @reaction_network begin
	@species pwc(t)=true
	α * β, S + I --> 2I
	r * m, I --> D
	r * (1 - m), I --> R
end

# ╔═╡ 51f0f4b0-0332-416f-8ade-582fa4ba0257
species(infection_model3)

# ╔═╡ 30bcea84-f7ca-4907-9a0f-54934e1731d0
md"""
We create the condition in the following way. When `pwc` is true then $I$ will be changed and also `pwc` will become `false`, so that the condition happens only once. We assume hereby that $I$ will never reach $0$!
"""

# ╔═╡ 87d9be87-3bc7-4442-a396-fb501355fe8c
condition3 = [infection_model3.I ~ 1e6*infection_model3.pwc] => [infection_model3.I ~ infection_model3.I - 0.999e6, infection_model3.pwc ~ false]

# ╔═╡ 467834b5-0063-4e91-b446-2828a1d44d78
md"""
The continuous event needs to be included in our model.
"""

# ╔═╡ 09ec4a21-c8fd-42fe-ae0c-f3a548c7032c
@named infection_model3_c = ReactionSystem(equations(infection_model3), continuous_events=condition3)

# ╔═╡ 7b93b776-caff-4608-9ec3-8519b9b95e51
md"""
After that, we need to *complete* our *reaction network model*.
"""

# ╔═╡ 96787705-f294-4f3c-9acf-8af7f6e0d579
infection_model3_c_com = complete(infection_model3_c)

# ╔═╡ 54d5f92c-de07-4c06-8bc5-e11e2b2c5764
md"""
Then we need to create a new ODE problem.
"""

# ╔═╡ e0562699-d065-4ae1-9c39-f9e0bf561fa3
oprob3 = ODEProblem(infection_model3_c_com, u0, tspan, params)

# ╔═╡ 6b3ea888-1238-46ae-9eaf-52e86c76bc1c
md"""
Finally, the ODE problem can be solved. Notice that you need to make a deepcopy of the ODE problem again.
"""

# ╔═╡ abe11178-7177-418c-98b4-da1afc56842e
osol3 = solve(deepcopy(oprob3), Tsit5(), saveat=0.1)

# ╔═╡ 23804c2a-d031-4025-92d3-5cdb74f87353
md"""
Now we can plot the results.
"""

# ╔═╡ 02173eaa-1178-4383-b57a-03e53ce38af8
plot(osol3, idxs=[:S, :I, :D, :R])

# ╔═╡ 6d1aa79b-8614-4a77-a327-f7e6962d1944
md"""
If you want to see the final values of $S$, $I$, $D$, $R$ and `pwc`, type:
"""

# ╔═╡ d04b8d97-e9d3-4428-a695-bfde4b44a291
osol3.u[end]

# ╔═╡ 70fd136e-5ded-4c30-818e-a5de61fbbf86
md"""
Try to interpret the results yourself. Ask yourself the following questions:

!!! question
	1. What are the trends in the results obtained?
"""

# ╔═╡ 33ff901b-6c7e-4b8c-a163-663d9443a213
md"""
- Answer: 
"""

# ╔═╡ 0d6c2233-35c5-4257-9dfd-398704080ba3
md"""
!!! question
	2. How much less casualties are there compared to not putting $999\,000$ individuals into isolation? (Hint: you also need to take into account the casualties in the $999\,000$ individuals that had been put into isolation.)
"""

# ╔═╡ f5d23dba-02a9-4d27-8eb1-ff89d2902cdf
md"""
- Answer: 
"""

# ╔═╡ Cell order:
# ╠═be326550-25ea-4c5f-ac4a-dd74d12bc89a
# ╠═cd166045-7788-4f5c-97d9-d1bedeceedae
# ╠═e6e748a5-b692-4c7c-b39e-ef0fe9e46c35
# ╠═3cfe381f-5871-4447-9149-fc023c787b42
# ╟─9ea49d7a-e524-11ee-1b68-9d1d71aaba24
# ╟─ba4dfc95-b74f-4d36-b34b-5eb2836a5cd6
# ╟─94440e99-7c0b-4d87-8038-141a2bc5fcb8
# ╟─83fa478c-ac54-482b-92de-8eb6a01ddfaf
# ╟─d257102b-8480-4ce3-bdba-50995ccbdc26
# ╟─df9a46ab-56eb-4d3d-a3a8-e1f2da59b32e
# ╟─7dd37224-071b-4ec7-88d5-39fd4482e615
# ╟─97d34ab6-a1c6-404d-af1f-7da85287757f
# ╟─b5e6cb25-da6f-475a-bc73-041c4c9256d5
# ╟─529ddac3-bb7a-4cc3-9bbd-ca31f5796e27
# ╟─f9ca0b86-1447-48c1-88e1-876180d7628a
# ╟─bf2b06d8-884b-4a35-b308-83a2f3d6696f
# ╟─1f405b09-53cd-43e0-b498-70eafef92acb
# ╟─39d80a3f-14d2-49d6-905b-171adf5930a0
# ╟─44c6cdc0-8610-43bd-b3b3-a7441ee62615
# ╠═b1bfe3b7-782c-46c9-86db-829942b6a9d0
# ╟─cd3da548-b24f-426f-8b46-2db70e0a979c
# ╠═7d2f8bf9-3a00-4631-84ff-1a36b6d7e19b
# ╟─4a504e36-8269-4822-8919-75b7f6dcaf6d
# ╟─bb5578b8-40d4-48c9-be87-cfaeadb8c5f9
# ╟─af6b169d-55c3-4384-9f57-ad629e02bad6
# ╟─52b2cc8f-1eac-47c1-8d42-09c1b91d1350
# ╠═5f1d6198-ccdc-4a63-bf98-2032a6d683ad
# ╟─59686253-3e0e-454e-8bc6-dde189ee2197
# ╠═fdb79a26-fbc5-482a-8cdb-760d18f57a77
# ╟─2bd02577-7321-46b5-9f42-1484ed86d1d3
# ╠═64455b26-4374-435f-b71b-8b9600cab263
# ╟─64a48327-03ec-46a0-9a76-79c10acb6e20
# ╠═78a22ba2-1e60-40b3-989c-922dcf9ca054
# ╟─eeb25e50-165b-4e93-8397-5b09fe8e7242
# ╟─0d39fbbc-4fcb-45eb-b24b-14f2e093c98c
# ╠═ca18795a-d270-493f-807f-a3b8c78aa6d7
# ╟─3f6f38a3-12ec-4b5b-8ebb-9e48190d3f7e
# ╠═72f7f870-b2fc-4df1-9325-3c587be7e014
# ╟─da5cfd6e-b57a-4235-b79d-a1c2f148328b
# ╠═a09b3a39-bedc-4180-ae5b-8cec708a64cb
# ╟─7fc8c671-b75c-4487-a013-779ee2422c8b
# ╠═19398992-dca7-441a-a5fd-7f9d0ad30be5
# ╟─3197244f-655b-4dca-80f3-794b30722551
# ╟─979afe85-b910-44a0-8ac0-6e719cb9157e
# ╟─1ba859fa-46a5-434f-a99c-e710ba85caf8
# ╠═35bd9a1a-bb4a-4285-98f0-853b03c95cb7
# ╟─95c1c0ea-51d0-47ee-8948-a5b87c42a70d
# ╟─57036f49-2f1f-4327-89ed-2c96098a1c22
# ╟─56ebb72d-2351-4e67-b268-1f48bbb77cb3
# ╠═a25d5925-a254-488b-b782-d29cff4470a2
# ╟─0952b1d1-24b4-4540-91cd-94f7a4dcbd57
# ╟─a235c7ce-f14a-4c7c-86a1-08aa5f2d9c85
# ╠═9a9440fa-d8a3-44bc-8037-4bf1f8af40b0
# ╟─c4e83ef8-9490-4361-a2a9-5abc45e242be
# ╟─aa1b904a-c8a9-41a4-9297-8d7c821d4b77
# ╠═c6d2dd69-8c61-4a40-894f-664b2d2d14be
# ╟─3c253bf3-886d-4e86-82ac-7751d23f342f
# ╟─14756171-9e8e-4cb0-b7af-74c2d649fe9f
# ╠═142e3e48-bf75-4498-ad0e-9f47cb921045
# ╟─a3599781-a690-4fa3-b483-cd47727935cb
# ╟─8ac90e5c-ce85-4113-bc1b-9dc17cf7e928
# ╟─0af3c166-46b6-455d-af6b-a72c4d2a5ce4
# ╠═513c037b-c54c-47fa-b97a-06f69a983386
# ╟─5672749b-283a-4e8e-8f90-76b5312b30ac
# ╠═01188934-ce0f-4d65-b880-caac4837d796
# ╟─f3cf01a1-3c65-4a37-bf9a-cd2233cb470b
# ╠═ff8f4c23-5695-48aa-9965-02677103f2c9
# ╟─b07f09c6-0515-4d27-9ca1-c45427a5988c
# ╠═ab4940b8-b8ec-4835-8d6e-4e57a5e2e464
# ╟─a8b68546-9401-4815-9f33-cc22b7a0579d
# ╠═89cc8970-693d-4885-b0b5-46110db4ca32
# ╟─dcdd693d-f2f3-406c-8647-31cdedcb823b
# ╠═0a2ed6a6-1535-4721-a6da-008f62b16942
# ╟─757d9f1e-5cba-4ae8-9172-e0c2fceb8458
# ╠═54bedb4d-51c5-4f58-bc4b-50cc469b7b04
# ╟─764e2f1a-f974-4916-8573-cacba897cf07
# ╟─2ae76ddb-71f5-49d7-a250-429d6c0138f6
# ╟─590f1b49-7442-4a71-af8c-8acdea071448
# ╟─45c1c238-a9f7-4f7b-a0ce-07b5bb4768d4
# ╟─b6baafc2-6d5e-43c3-8ef9-845961cdd20b
# ╠═cd32beba-67cf-4b12-a77b-99f96263f0a4
# ╟─ae38c663-0ee4-409e-bfca-5f13ed88b67d
# ╠═d44da6c6-c93d-4c61-8125-9eee464c897e
# ╠═00a72697-d36a-41cc-9eec-8e821829ce0e
# ╠═cf39b4cf-9cd0-4755-80db-ca4aea7c1084
# ╠═52901bbf-e47b-4da1-95c4-f0869812398c
# ╟─102b4fbc-23b1-46ed-bb72-124eb88517ce
# ╟─f6cbafbf-98b4-4286-86d2-fe95821a5ff4
# ╟─4242ec68-b36d-4f7f-9001-2c90c66d5b9b
# ╟─f5f98168-119c-4b49-a8e3-a3cc775faeb0
# ╟─30b0e2b3-2b8d-43a5-9587-9fcb0b89b258
# ╟─ade413c2-d7d9-4250-8490-75534900a389
# ╟─58730ac6-d83b-420d-a000-2f50545f0d39
# ╠═7411474c-fac7-4b7c-8ded-4c2df5956fb0
# ╟─e6949052-6d12-4414-ac67-cb9435b46290
# ╠═8752373f-a602-437b-9b3a-2602c8babf87
# ╟─ee596a6c-29bf-4d18-b687-be943c128aa7
# ╠═40f849cb-20f9-4bbe-806e-512abd6f3210
# ╟─d4e537da-490b-4781-9a76-c89d762849ec
# ╠═4a6f357c-afd7-4b0d-b1d2-15f5c0e4298b
# ╟─e76e77be-9f15-4aa5-8506-cdb32a6ec9b1
# ╠═37b11fdb-8c16-4ec4-a7ad-e2fdbaf4ea5d
# ╟─59f4f4dd-f0b5-4354-b53c-4554c878d5a8
# ╠═0af84963-ca64-4958-a544-42d445da5a7c
# ╟─18151ab1-1e4e-48d8-be70-4fa4c8a43af1
# ╠═7d7ed974-7c9c-4c30-ab2a-e3028ce702dd
# ╟─7375edbc-24a4-4300-bdef-2686cb377cfc
# ╟─6eb368f5-8d43-4fb8-b6eb-c783675d1dd9
# ╟─2ba54805-e901-4281-9b07-cc7137b8c809
# ╟─fec6e165-bd1a-4e57-b7f3-ae712e03276c
# ╟─cb8a6f77-f08c-4fc8-9445-bd1c17521fcc
# ╟─5005a09d-a844-4f9d-a058-0d93583d5bab
# ╠═f77d2d75-468f-4746-b81e-0bbbf33fc8d7
# ╠═51f0f4b0-0332-416f-8ade-582fa4ba0257
# ╟─30bcea84-f7ca-4907-9a0f-54934e1731d0
# ╠═87d9be87-3bc7-4442-a396-fb501355fe8c
# ╟─467834b5-0063-4e91-b446-2828a1d44d78
# ╠═09ec4a21-c8fd-42fe-ae0c-f3a548c7032c
# ╟─7b93b776-caff-4608-9ec3-8519b9b95e51
# ╠═96787705-f294-4f3c-9acf-8af7f6e0d579
# ╟─54d5f92c-de07-4c06-8bc5-e11e2b2c5764
# ╠═e0562699-d065-4ae1-9c39-f9e0bf561fa3
# ╟─6b3ea888-1238-46ae-9eaf-52e86c76bc1c
# ╠═abe11178-7177-418c-98b4-da1afc56842e
# ╟─23804c2a-d031-4025-92d3-5cdb74f87353
# ╠═02173eaa-1178-4383-b57a-03e53ce38af8
# ╟─6d1aa79b-8614-4a77-a327-f7e6962d1944
# ╠═d04b8d97-e9d3-4428-a695-bfde4b44a291
# ╟─70fd136e-5ded-4c30-818e-a5de61fbbf86
# ╟─33ff901b-6c7e-4b8c-a163-663d9443a213
# ╟─0d6c2233-35c5-4257-9dfd-398704080ba3
# ╟─f5d23dba-02a9-4d27-8eb1-ff89d2902cdf
