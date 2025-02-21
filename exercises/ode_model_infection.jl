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

# ╔═╡ 989fd8c8-25d9-47b9-ade6-6c7f21a7dceb
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ 18a4df05-0349-400d-a29e-b3fa71aa4d88
using Markdown

# ╔═╡ 65b3536b-c9b3-485e-87e1-4f6657381503
using InteractiveUtils

# ╔═╡ e968f6d5-962c-4de8-a9a3-ec116805d5e1
using PlutoUI; TableOfContents()

# ╔═╡ e97637b0-446c-44ac-bd08-c56632f9b57f
using Catalyst

# ╔═╡ c992a16d-80d3-48d6-bf43-2dbfa8e2f081
using DifferentialEquations, Plots

# ╔═╡ c54c7e3c-f3d1-40a5-88f7-bedcab3269ba
md"
# Exercises - infection model
"

# ╔═╡ fab49cb7-c41e-498a-98f6-33b4821ceb90
md"""
We will work here with the same infection model as in the **Introdution to Catalyst** (revisit the concerned notebook if necessary). We shortly summarize some important aspects of the model and give a condensed version of the solution method and the examples.
"""

# ╔═╡ 8ff0b57e-acc1-4ebd-a562-a82c9efa7ffc
md"""
The **state variables**:

| Variable | Unit | Meaning |
|:---------- |:---------- |:------------|
| ``S``    | *persons* | number of susceptible persons             |
| ``I``    | *persons* | number of infected persons             |
| ``D``    | *persons* | number of deceased persons             |
| ``R``    | *persons* | number of recovered persons             |
"""

# ╔═╡ e36b150a-b4fc-476f-bab7-e1162a4b263d
md"""
The **parameters**:

| Variable | Unit | Meaning |
|:---------- |:---------- |:------------|
| ``\alpha`` | ``\frac{persons}{contact}`` | chances of getting infected after contact |
| ``\beta`` | ``\frac{contact}{persons^2\,day}`` | contact rate |
| ``r`` | ``\frac{1}{day}`` | rate of leaving infection period |
| ``m`` | ``\frac{person}{person}`` | fraction of persons deceasing |
| ``1-m`` | ``\frac{person}{person}`` | fraction of persons becoming resistant |
"""

# ╔═╡ bdefa657-88bd-43d6-8e54-cc6321a9e720
md"""
The infection model has three reaction events:

- **Infection**, where a susceptible persons meets an infected person and also becomes infected. The infection rate is $\alpha \beta$.
- **Decease**, where an infected person dies. The death rate is $m r$.
- **Recovery**, where an infected person recovers. The recovery rate is $(1-m) r$.
"""

# ╔═╡ afb751b9-39b8-4430-af4b-02f010667518
md"""
The involved reactions are:

$$S + I \xrightarrow[]{\alpha \beta} 2I$$
$$I \xrightarrow[]{mr} D$$
$$I \xrightarrow[]{(1-m)r} R$$
"""

# ╔═╡ aab497eb-c0bd-49e3-a0c1-71e13f5b0ad5
md"""
Load the Catalyst package:
"""

# ╔═╡ a91f1024-4800-4ffe-8bb2-5f6cf3c2bde1
md"""
## Examples
"""

# ╔═╡ 47f4a7cc-1a74-4ff0-ae21-2a5e4c495935
md"
### Implementation of the system
"

# ╔═╡ 13d94ba6-425e-471b-aa99-2f321891d1f0
md"""
Implement the *reaction model* in Catalyst:
"""

# ╔═╡ 3f6080e1-ab43-47f7-a82b-95956bf4cafd
infection_model = @reaction_network begin
	α * β, S + I --> 2I
	r * m, I --> D
	r * (1 - m), I --> R
end

# ╔═╡ 5d5fae40-a4ed-4908-ba65-f3a16c38ab4f
md"""
The species:
"""

# ╔═╡ fd347a58-5cad-4073-aadb-176fa54dcdfa
species(infection_model)

# ╔═╡ bd5bd448-8741-4528-989f-9754b0e55b9b
md"""
Alternatively:
"""

# ╔═╡ a78b178e-6952-4057-bd8c-6d5b0c5e0518
@unpack S, I, D, R = infection_model;

# ╔═╡ cbfe8833-786c-48f7-b849-705006d4c41d
md"""
The parameters:
"""

# ╔═╡ 96f3dd68-da80-4f07-af2b-e3b9ec05cc0c
parameters(infection_model)

# ╔═╡ cc4d514a-dfed-4424-9982-8315954b38ff
md"
Convert the *reaction model* if you want to see the symbolic differential equation model:
"

# ╔═╡ 3aad1f32-7daf-49a5-9dd3-9ddbba1a50de
osys  = convert(ODESystem, infection_model)

# ╔═╡ a432a88a-8db3-4729-b9a3-8ca9322fc81e
md"""
We can get a list of the differential equations, the state variables and the parameters:
"""

# ╔═╡ 7e268605-8314-4e7b-8c38-a65212e148a9
equations(osys)

# ╔═╡ 0707b599-c1ca-459b-b87b-956f8f49564b
unknowns(osys)

# ╔═╡ 113ff67e-6dc8-45cb-9617-7768e132f6e9
parameters(osys)

# ╔═╡ 39e396e1-20fa-4e78-b0ab-99774ce55f0f
md"""
### Simulating the system as an ODE-problem

Load the packages Differential and Plot:
"""

# ╔═╡ 4a581cce-9a1b-4516-9cf1-6a41ea6566e7
md"""
### Setting initial conditions, timespan and parameter values
"""

# ╔═╡ c5ee2780-c652-4c9b-981b-53e5f91e1766
u0 = [:S => 9_999_000.0, :I => 1_000.0, :D => 0.0, :R => 0.0]

# ╔═╡ fc380000-33a0-419d-9226-7a0c136de8c6
tspan = (0.0, 90.0)

# ╔═╡ a8756838-93d7-4909-84c2-6d89be4fa711
params = [:α => 0.08, :β => 1.0e-6, :r => 0.2, :m => 0.4]

# ╔═╡ 6d2580f7-f882-449d-a9b8-629ab41e2671
md"""
### Creating and solving the ODE problem and plotting results
"""

# ╔═╡ d3f39b29-71f5-4674-a52f-3ba63279fca7
oprob = ODEProblem(infection_model, u0, tspan, params)

# ╔═╡ 6173506b-8f56-4100-a0a0-117c1ee8aa0f
osol = solve(oprob, Tsit5(), saveat=0.5)

# ╔═╡ e6b17756-a7e9-4f50-98dc-b69f9a54903f
plot(osol)

# ╔═╡ 006bfc98-3c95-4763-b0ed-f4714ad5bd02
plot(osol, idxs=[:S, :I])     # only S and I
# plot(osol, idxs=[S, I])     # also possible when S and I are unpacked

# ╔═╡ 4c5d2e5d-3237-4886-b852-3142fb996dd6
plot(osol, idxs=(:S, :I), xlab="S", ylab="I")    # phase plot I vs S

# ╔═╡ f9bd32f5-c473-44d4-af9d-4d99032a2d21
osol.u[end] 	# final values in the order of defined species

# ╔═╡ 000495a6-e91d-438e-ab6d-7a5e95a99f27
md"""
### Example 1 - Influence of $r$

Influence of the duration of infection, $1/r$, for average infection periods  between $1$ and $10$ days of being contagious ($r$ between $0.1$ and $1.0$, with a step of $0.1$, and default value $0.1$).
"""

# ╔═╡ d73f8ed8-ac25-4c1c-b834-989a2a929df3
@bind r Slider(0.1:0.1:1, default=0.1, show_value=true)

# ╔═╡ 4ca0e1df-7765-4f36-a442-2749ea493426
params1 = [:α => 0.08, :β => 1.0e-6, :r => r, :m => 0.4] 	# r specified by slider

# ╔═╡ d56b486e-44e5-4db7-9c4e-00db5c1b733f
oprob1 = ODEProblem(infection_model, u0, tspan, params1);

# ╔═╡ 298ada89-4581-41c5-a89f-f01ae431c5ec
osol1 = solve(oprob1, Tsit5(), saveat=0.5);

# ╔═╡ 5b4ca642-0f89-40ab-9a2f-c35fc603a91d
plot(osol1, ylim=(0, 1e7))

# ╔═╡ 171208da-69a3-4737-85e0-fae0805bbaa8
md"""
Change the value of $r$ in the `params1` vector to visualize the effect in the plot.
"""

# ╔═╡ 98aedcf5-e0de-442f-b19b-c5b3d9e3d475
md"""
### Example 2 - Discrete event

Suppose that regulations are such that on day 14, people need to reduce their contacts by 50%.
"""

# ╔═╡ 462b2786-98d6-40cb-887d-aff6b63881e8
condition2 = [14.0] => [infection_model.β ~ infection_model.β/2]

# ╔═╡ efd74885-c25c-42e5-bef5-c9654329c220
@named infection_model2 = ReactionSystem(equations(infection_model), discrete_events=condition2)

# ╔═╡ ecf20050-12d3-48e7-84f1-6919a9aab865
infection_model2_com = complete(infection_model2)

# ╔═╡ fdea5a92-e115-4dbe-a885-95772857c686
oprob2 = ODEProblem(infection_model2_com, u0, tspan, params)

# ╔═╡ 35d6bce7-9778-4f42-bdcd-271cedd8bad9
osol2 = solve(deepcopy(oprob2), Tsit5(), saveat=0.5)

# ╔═╡ 31ceb303-dbb1-4f4d-9082-14ea1d83ab66
# We can compare the result now to the solution without contact reduction
begin
	plot(osol2)
	plot!(osol; linestyle=:dash, label=:none, color=:grey, lw=0.5)
end

# ╔═╡ 14143212-53f1-4753-8c25-e22559cdbccd
osol2.u[end]

# ╔═╡ 9c668254-3c54-4d95-9e1e-e2bd2d3eb1df
md"""
!!! question
	How do we interpret this event? Is the number of deceased and infections reduced? How much?
"""

# ╔═╡ 130a6ce2-dd8c-4a05-8863-953a77a407aa
osol2[:D][end]/osol[:D][end] - 1 	# The number of deceased is reduced 13%

# ╔═╡ 9abcc2f9-5c24-4840-8a14-fa7a197ffef0
osol2[:S][end]/osol[:S][end] - 1 	# 6 times less infections with contact measures

# ╔═╡ 1f970a9e-abde-422b-a793-eb51ad2001bc
osol2[:R][end]/osol[:R][end] - 1 	# 13% less recoveries (due to fewer infections)

# ╔═╡ bdf9cfa0-9dc4-4ea5-9024-63813ec30d09
md"""
### Example 3 - Continuous event

Suppose that when the number of infected individuals reaches $1\,000\,000$, then $999\,000$ of them are promptly put into isolation (or removed from the population).
"""

# ╔═╡ cbf0a768-f6ac-43fb-8d05-e7dac51d4cc0
infection_model3 = @reaction_network begin
	@species pwc(t)=true
	α * β, S + I --> 2I
	r * m, I --> D
	r * (1 - m), I --> R
end

# ╔═╡ fee597e0-97cb-40bc-a5b9-166631e8b9f6
condition3 = [infection_model3.I ~ 1e6*infection_model3.pwc] => [infection_model3.I ~ infection_model3.I - 0.999e6, infection_model3.pwc ~ false]

# ╔═╡ 51363f3c-7aa9-48ed-808e-d1f7a4aadc0c
@named infection_model3_c = ReactionSystem(equations(infection_model3), continuous_events=condition3)

# ╔═╡ cd060e93-7a92-4824-b812-080391fbf554
infection_model3_c_com = complete(infection_model3_c)

# ╔═╡ a41878a2-2fd3-4a2a-ba85-d2a09c67da33
oprob3 = ODEProblem(infection_model3_c_com, u0, tspan, params)

# ╔═╡ fc2dc11e-7022-4c6f-b7ed-5b077fc21ac5
osol3 = solve(deepcopy(oprob3), Tsit5(), saveat=0.5)

# ╔═╡ 2b511237-7739-44e5-b6bd-26a2fdc76bb5
begin
	plot(osol3)
	plot!(osol; linestyle=:dash, label=:none, color=:grey, lw=0.5)
end

# ╔═╡ 3a131f92-aa09-4e8c-98a1-0c4c1c74878a
osol3.u[end]

# ╔═╡ 96792abd-2b86-4245-aaaf-15c43b4eabb0
md"""
!!! question
	How do we interpret this new event? Is this a better measure than contact reduction alone? Would you know how we call such an event?
"""

# ╔═╡ 1116b608-e6b9-4e6d-a0eb-1ac58ef9e2f7
0.4*0.999e6 	# Amount of people deceased during isolation (infected)

# ╔═╡ 4111ae36-ddef-42a8-b8a9-61f9f464a5c0
(osol3[:D][end] + 0.4*0.999e6)/osol[:D][end] - 1 	# Deceased are reduced only 1%!

# ╔═╡ 6c8a9e97-0a60-4e42-ad7f-78c0dece2b30
osol3[:S][end]/osol[:S][end] - 1 	# Infections have increased 56% (wrt. example 1)

# ╔═╡ 801a038b-7485-40b6-b7a5-92b0b3c11c0e
osol3[:S][end]/osol2[:S][end] - 1 	# Amount of population exposed decreased 79%

# ╔═╡ 62adbc28-cec2-4da8-9842-2153563497d0
md"""
- Answer: 
"""

# ╔═╡ 3bc2d36a-6efc-4631-baae-0d3e62bf7e3f
md"""
!!! hint
	Isolation alone is not effective to protect the part of the population that's already been infected. The peak of infections is not avoided, only delayed. However, the proportion of the population exposed is much lower thanks to isolation. *Would then a combined set of rules be best in that case?*
"""

# ╔═╡ ad83a851-e424-4e8f-bae9-38b9283db2fa
md"""
## Exercises
"""

# ╔═╡ b580ef9b-7d53-4bf0-8511-213920d59ee2
md"""
### Exercise 1 - Influence of $\alpha$

Evaluate the effect of a decreasing risk of infection after contact with an infected person, i.e. $r = 0.2$, $\beta = 0.1$ and $\alpha$ between $8\%$ and $20\%$.

Use the same initial values and timespan as before.
"""

# ╔═╡ 090a949f-c18f-4e07-8c94-21a71bb0ceed
md"""
Make a slider for $\alpha$ in the range of $0.08$ and $0.20$ with a step of $0.02$. Take a default value of $0.08$.
"""

# ╔═╡ 85e81822-2ea6-4d05-b879-1791bff255b9
# missing                  # Uncomment and complete the instruction

# ╔═╡ a2a109e3-ea74-4e74-91d0-6ea12ddf28d4
md"
Initialize vector `params_ex1` with parameter values:
"

# ╔═╡ 27a4c4f5-0b2b-4d23-9bb0-a1a8ff6cfc4d
# params_ex1 = missing      # Uncomment and complete the instruction

# ╔═╡ e0a8a397-4500-47de-8be3-49d3174648b1
md"""
Create the ODE problem and store it in `oprob_ex1`:
"""

# ╔═╡ 7955722d-fc35-4cce-8bf2-ee40ee7fbc82
# oprob_ex1 = missing;     # Uncomment and complete the instruction

# ╔═╡ 0590bc0b-45d6-4a8e-8078-af057984523f
md"""
Solve the ODE problem and store the solution in `osol_ex1`:
"""

# ╔═╡ f7b09e96-cf30-4bd0-841b-46621df693ab
# osol_ex1 = missing;       # Uncomment and complete the instruction

# ╔═╡ e7947165-b244-4ff6-bdf3-61d03aefe696
md"""
Plot the solutions:
"""

# ╔═╡ adc0aa15-e0a9-4549-9609-967a9dc78b78
# missing                   # Uncomment and complete the instruction

# ╔═╡ 7a1effe4-776c-4143-8085-f98a213a2cc3
md"""
Change the value of $\alpha$ in the `params_ex1` vector to visualize the effect in the plot.
"""

# ╔═╡ f85bd6be-8984-49f0-a36b-7aa7bb937c54
md"""
!!! warning "Tip"
	You can move the cell containing the slider next to the figure to visualize the changes in $\alpha$.
"""

# ╔═╡ c14d45b8-3af1-4e05-b10a-081ccd62a34d
md"""
Try to interpret the results yourself.
Ask yourself the following questions:
!!! question
	1. What are the trends in the obtained results?
"""

# ╔═╡ 04f89fb3-f856-4206-9871-4d92f2816332
md"""
- Answer: 
""" 

# ╔═╡ a39eeec9-d816-4d35-90ae-e602a24bb056
md"""
!!! question
	2. How can this be explained from the model structure?
"""

# ╔═╡ fcad8ede-8295-46d6-a637-728b3986a761
md"""
- Answer: 
"""

# ╔═╡ b3f08d61-187d-4f99-b95f-508ff1b3d52d
md"""
### Exercise 2 - Administration of medicinal products

Scientists have developed a medicine that heals sick people and makes them immune to
the disease. After administering medication, the infection duration is reduced to two
days. All treated patients heal and acquire immunity to the virus. The model will have
to be extended with two additional parameters.

- Parameter $b$: the fraction of infected persons undergoing treatment.
- Parameter $h$: the rate at which the infected persons treated are no longer contagious ($day^{-1}$).

Administering the drug to a fraction of the infected individuals affects two *reactions*: $I \rightarrow D$ and $I \rightarrow R$, with the following assumptions:

- The fraction of infected persons treated ($b$) has a reduced infection duration.
- The fraction of infected individuals not receiving treatment ($1 − b$) still has the same duration of infection.
- The mortality rate $m$ only affects the group of sick people who were not given any medication.
- All treated individuals recover.
- A fraction of the untreated individuals also heals.

Check the effect on the epidemic when $0\%$, $25\%$, $50\%$, $75\%$ and $100\%$ of infected individuals are treated with $h = 0.5$. Use the same initial conditions and timespan as before.
"""

# ╔═╡ 08ba4614-e4fa-4203-9e28-c5c4bd9fcdd5
md"""
Set-up the new *reaction network/model* and name it `infection_med`:
"""

# ╔═╡ e4087441-7d0d-4716-b9fa-eef45e9f3b87
# infection_med = @reaction_network begin  # Uncomment and complete the instruction
# 	α * β, S + I --> 2I
# 	..., I --> D
# 	(..., ...), I --> R
# end

# ╔═╡ 87b3796d-16c5-4fe9-a32e-d8bf7c754254
md"""
Convert to an ODE system. Check the differential equations and make sure you understand each term.
"""

# ╔═╡ 2a66417b-79a2-4909-bd89-5650c49b9411
# osys_ex2 = missing       # Uncomment and complete the instruction

# ╔═╡ 35956a99-5729-4bf7-bf09-b2232fc958da
md"""
Set-up parameter values:
"""

# ╔═╡ cc31567f-631b-41b9-aa2d-3408ca951bcf
# params_ex2 = missing       # Uncomment and complete the instruction

# ╔═╡ 82f8c50c-511d-4399-ac87-0890b678c6c1
md"""
Create the ODE problem and store it in `oprob_ex2`:
"""

# ╔═╡ ec08adb2-561c-4f17-a0b5-019d7b1f1098
# oprob_ex2 = missing;       # Uncomment and complete the instruction

# ╔═╡ c55ec157-6ec8-45db-9a54-0f39aceb9607
md"""
Solve the ODE problem and store the solution in `osol_ex2`:
"""

# ╔═╡ 67668918-624e-4adf-be5d-6fdbc62555f6
# osol_ex2 = missing;        # Uncomment and complete the instruction

# ╔═╡ 832f984a-176d-4187-8476-5026d29a63d8
md"""
Plot the solutions:
"""

# ╔═╡ 5f9a9838-099d-457f-9402-915f42d0cd33
# missing          # Uncomment and complete the instruction

# ╔═╡ ef3ca25d-c110-452d-877a-6304ae6cd5a7
md"""
Change the value of $b$ in the `params_ex2` vector to visualize the effect in the plot. Interpret the obtained plots.
"""

# ╔═╡ f8dbab24-55aa-466f-a978-5f8db93b4b93
md"""
Try to answer the following questions:

!!! question
	1. Why does the peak in the number of infected individuals shift to the right when the value of $b$ increases?
"""

# ╔═╡ 381d3a54-0e4f-4203-874f-4a4716d481d4
md"""
- Answer: 
"""

# ╔═╡ fd7d89f3-9ab8-4b17-b26f-571bad862718
md"""
!!! question
	2. Why does the number of recovered individuals first rise when the value of $b$ increases and then fall when the value of $b$ continues to increase?
"""

# ╔═╡ 96de4f0d-c2d1-4fad-bded-a314cbca848b
md"""
- Answer: 
"""

# ╔═╡ dedb1ac1-cb54-416f-a092-12a50080c69f
md"""
Check the number of fatalities:
"""

# ╔═╡ a912f449-6c9d-4595-863f-a2ba523bad95
# missing

# ╔═╡ 089707c3-1df6-4799-a87b-0d4872a0267d
md"""
### Exercise 3 - Adding vaccination to the model

Scientists have developed a vaccine that makes healthy people immediately immune to
the disease.

Vaccination affects several differential equations:
- Susceptible individuals are vaccinated at a rate of $v$ (with unit $day^{-1}$). These persons can therefore no longer be infected.
- The vaccinated persons become resistant.

We are going to use a vaccination rate $v$ so that the number of fatalities is about 10 times smaller after a period of 90 days compared to those in absence of vaccination (cf. Exercise 2).

The vaccination programme is launched $2$ days after the outbreak of the disease.

Assume that individuals are still being treated ($b = 0.2$ and $h = 0.5$). Extend the model obtained in the previous exercise for the launch of a vaccination campaign after the outbreak of the disease.

Find out via trial and error what the minimum vaccination rate need to be so that the number of fatalities is $10$ times smaller after a period of $90$ days compared to those in absence of vaccination (cf. Exercise 2). Consider an initial step size in $v$ of $0.01$ and then fine tune with a step size of $0.001$.

Use the same initial values and timespan as before.
"""

# ╔═╡ 381a714d-699e-45c6-909a-02689b2a7e6b
md"""
Set-up the new *reaction network/model* and name it `infection_med_vac`:
"""

# ╔═╡ a1795123-876b-4bd3-ac3a-711367b500d1
# Uncomment and complete the instruction
# infection_med_vac = @reaction_network begin
# 	α * β, S + I --> 2I
# 	..., I --> D
# 	(..., ...), I --> R
# 	..., ... --> ...
# end

# ╔═╡ 6b46d116-0da3-428a-b25f-221622dfa7e0
md"""
Convert to an ODE system. Check the differential equations and make sure you understand each term.
"""

# ╔═╡ 710fe13c-bd43-4d1e-867e-72330722ac1e
# osys_ex3 = missing       # Uncomment and complete the instruction

# ╔═╡ 2ee20493-2dc2-4bd4-8ec3-78032f1fcc2d
md"""
Make a slider and bind it to the variable `v`. Use a range $[0.0, 0.1]$, step size $0.001$ and default value of $0.0$.
"""

# ╔═╡ 2943b515-fb4b-45a2-88dd-d0267ec95b09
# missing                # Uncomment and complete the instruction

# ╔═╡ 0274fbd8-c025-4e36-bbb0-f65f21b962c7
md"""
Set-up parameter values:
"""

# ╔═╡ 34c4d746-45f7-43c4-a0f0-6e570373a35d
# params_ex3 = missing    # Uncomment and complete the instruction

# ╔═╡ 060a5c5b-9834-42e2-98e9-ac23e7403b60
md"""
Create the ODE problem and store it in `oprob_ex3`:
"""

# ╔═╡ c61ae8b6-2324-43fa-a5ae-274a919af559
# oprob_ex3 = missing;       # Uncomment and complete the instruction

# ╔═╡ 69983b90-b6ef-4732-a4af-a772bb1364e5
md"""
Solve the ODE problem for (step wise) increasing values of $v$ and store the solution in `osol_ex3_vac`. Consider an initial step size in $v$ of $0.01$ and then fine tune with a step size of $0.001$.
"""

# ╔═╡ 16211bc6-cd97-43c0-8faf-25bb490930b6
# osol_ex3_no_vac = missing;   # Uncomment and complete the instruction

# ╔═╡ 461eada9-5f9c-4439-8b65-226382b6d148
md"""
First, put the value of $b$ in Exercise 2 to $0.2$. Compare the latter with the number of fatalities when no vaccination is/was available (cf. Exercise 2) by setting up a condition (a boolean expression return either `true` or `false`) here below where the final number of fatalities (with vaccination) divided by 10 is compared with (use larger than or smaller than) the number of fatalities (without vaccination):
"""

# ╔═╡ 0994d790-0fe8-46ab-bba2-a0a4ea466f64
# missing                     # Uncomment and complete the instruction

# ╔═╡ 8bb5dafd-6c26-4634-8be4-c2963efba056
md"""
Once you have found the required value of $v$ launch the vaccination programme $2$ days after the outbreak. Set-up the $2$-day time condition and store it in `condition_ex3`:
"""

# ╔═╡ 26c77206-fe3c-45c5-a627-f2b346c66be7
# condition_ex3 = missing     # Uncomment and complete the instruction

# ╔═╡ 527e55a8-1744-4c40-a412-054bed6d77f0
md"""
Make a new *reaction system* where the discrete event is included. Name it `infection_med_vac_c`.
"""

# ╔═╡ ceae9947-ed0f-4975-a944-cc938cf22dde
# @named infection_med_vac_c = missing    # Uncomment and complete the instruction

# ╔═╡ 6e767e36-f54e-4472-881a-ab04ad9c9d09
md"""
Complete the new *reaction system*. Name it `infection_med_vac_c_com`.
"""

# ╔═╡ 870c9d65-2229-40b7-9959-59afecb4f3bf
# infection_med_vac_c_com = missing      # Uncomment and complete the instruction

# ╔═╡ 93ca3970-bebd-4192-bc62-6ca08d58a6c2
md"""
Create the ODE problem and store it in `oprob_ex3_c`:
"""

# ╔═╡ 7e4114d6-932a-4def-9bb5-b534072ea513
# oprob_ex3_c = missing                # Uncomment and complete the instruction

# ╔═╡ 1bf8bc58-6848-4783-a614-7ac11646de92
md"""
Solve the ODE problem. Make a deepcopy and use `Tsit5()` and `saveat=0.5`. Store the solution in `osol_ex3`. 
"""

# ╔═╡ 8814c42b-16c7-466a-8ad8-58ee1c0921b2
# osol_ex3 = missing;                   # Uncomment and complete the instruction

# ╔═╡ 25158ea5-1930-4fee-aab4-490b475d8635
md"""
Plot the solutions:
"""

# ╔═╡ ed9472a0-09c9-4d1c-8f1d-1b7d87b97640
# missing          # Uncomment and complete the instruction

# ╔═╡ 6e213afb-fd4b-4ab9-8986-a319bee6be0c
md"""
Check the number of fatalities now and compare to the case without vaccination.
"""

# ╔═╡ 3af06996-bd43-4783-a5e0-3c11c23ca280
# missing 		   # Uncomment and complete the instruction

# ╔═╡ 93400573-ad1f-408d-b514-854998e7836a
md"""
!!! question
	3. What can you say about the rate of infection in the case of no vaccination? How would you measure this in the plot?
"""

# ╔═╡ 008cc499-fce9-43e7-b969-2419b847d24c
md"""
- Answer: 
"""

# ╔═╡ 14c3b2f1-14b4-415c-85c6-922361423826
md"""
!!! hint
	The rate of infection has to do with how many people get infected during a certain period...
"""

# ╔═╡ Cell order:
# ╠═18a4df05-0349-400d-a29e-b3fa71aa4d88
# ╠═65b3536b-c9b3-485e-87e1-4f6657381503
# ╠═989fd8c8-25d9-47b9-ade6-6c7f21a7dceb
# ╠═e968f6d5-962c-4de8-a9a3-ec116805d5e1
# ╟─c54c7e3c-f3d1-40a5-88f7-bedcab3269ba
# ╟─fab49cb7-c41e-498a-98f6-33b4821ceb90
# ╟─8ff0b57e-acc1-4ebd-a562-a82c9efa7ffc
# ╟─e36b150a-b4fc-476f-bab7-e1162a4b263d
# ╟─bdefa657-88bd-43d6-8e54-cc6321a9e720
# ╟─afb751b9-39b8-4430-af4b-02f010667518
# ╟─aab497eb-c0bd-49e3-a0c1-71e13f5b0ad5
# ╠═e97637b0-446c-44ac-bd08-c56632f9b57f
# ╟─a91f1024-4800-4ffe-8bb2-5f6cf3c2bde1
# ╟─47f4a7cc-1a74-4ff0-ae21-2a5e4c495935
# ╟─13d94ba6-425e-471b-aa99-2f321891d1f0
# ╠═3f6080e1-ab43-47f7-a82b-95956bf4cafd
# ╟─5d5fae40-a4ed-4908-ba65-f3a16c38ab4f
# ╠═fd347a58-5cad-4073-aadb-176fa54dcdfa
# ╟─bd5bd448-8741-4528-989f-9754b0e55b9b
# ╠═a78b178e-6952-4057-bd8c-6d5b0c5e0518
# ╟─cbfe8833-786c-48f7-b849-705006d4c41d
# ╠═96f3dd68-da80-4f07-af2b-e3b9ec05cc0c
# ╟─cc4d514a-dfed-4424-9982-8315954b38ff
# ╠═3aad1f32-7daf-49a5-9dd3-9ddbba1a50de
# ╟─a432a88a-8db3-4729-b9a3-8ca9322fc81e
# ╠═7e268605-8314-4e7b-8c38-a65212e148a9
# ╠═0707b599-c1ca-459b-b87b-956f8f49564b
# ╠═113ff67e-6dc8-45cb-9617-7768e132f6e9
# ╟─39e396e1-20fa-4e78-b0ab-99774ce55f0f
# ╠═c992a16d-80d3-48d6-bf43-2dbfa8e2f081
# ╟─4a581cce-9a1b-4516-9cf1-6a41ea6566e7
# ╠═c5ee2780-c652-4c9b-981b-53e5f91e1766
# ╠═fc380000-33a0-419d-9226-7a0c136de8c6
# ╠═a8756838-93d7-4909-84c2-6d89be4fa711
# ╟─6d2580f7-f882-449d-a9b8-629ab41e2671
# ╠═d3f39b29-71f5-4674-a52f-3ba63279fca7
# ╠═6173506b-8f56-4100-a0a0-117c1ee8aa0f
# ╠═e6b17756-a7e9-4f50-98dc-b69f9a54903f
# ╠═006bfc98-3c95-4763-b0ed-f4714ad5bd02
# ╠═4c5d2e5d-3237-4886-b852-3142fb996dd6
# ╠═f9bd32f5-c473-44d4-af9d-4d99032a2d21
# ╟─000495a6-e91d-438e-ab6d-7a5e95a99f27
# ╠═d73f8ed8-ac25-4c1c-b834-989a2a929df3
# ╠═4ca0e1df-7765-4f36-a442-2749ea493426
# ╠═d56b486e-44e5-4db7-9c4e-00db5c1b733f
# ╠═298ada89-4581-41c5-a89f-f01ae431c5ec
# ╠═5b4ca642-0f89-40ab-9a2f-c35fc603a91d
# ╟─171208da-69a3-4737-85e0-fae0805bbaa8
# ╟─98aedcf5-e0de-442f-b19b-c5b3d9e3d475
# ╠═462b2786-98d6-40cb-887d-aff6b63881e8
# ╠═efd74885-c25c-42e5-bef5-c9654329c220
# ╠═ecf20050-12d3-48e7-84f1-6919a9aab865
# ╠═fdea5a92-e115-4dbe-a885-95772857c686
# ╠═35d6bce7-9778-4f42-bdcd-271cedd8bad9
# ╠═31ceb303-dbb1-4f4d-9082-14ea1d83ab66
# ╠═14143212-53f1-4753-8c25-e22559cdbccd
# ╟─9c668254-3c54-4d95-9e1e-e2bd2d3eb1df
# ╠═130a6ce2-dd8c-4a05-8863-953a77a407aa
# ╠═9abcc2f9-5c24-4840-8a14-fa7a197ffef0
# ╠═1f970a9e-abde-422b-a793-eb51ad2001bc
# ╟─bdf9cfa0-9dc4-4ea5-9024-63813ec30d09
# ╠═cbf0a768-f6ac-43fb-8d05-e7dac51d4cc0
# ╠═fee597e0-97cb-40bc-a5b9-166631e8b9f6
# ╠═51363f3c-7aa9-48ed-808e-d1f7a4aadc0c
# ╠═cd060e93-7a92-4824-b812-080391fbf554
# ╠═a41878a2-2fd3-4a2a-ba85-d2a09c67da33
# ╠═fc2dc11e-7022-4c6f-b7ed-5b077fc21ac5
# ╠═2b511237-7739-44e5-b6bd-26a2fdc76bb5
# ╠═3a131f92-aa09-4e8c-98a1-0c4c1c74878a
# ╟─96792abd-2b86-4245-aaaf-15c43b4eabb0
# ╠═1116b608-e6b9-4e6d-a0eb-1ac58ef9e2f7
# ╠═4111ae36-ddef-42a8-b8a9-61f9f464a5c0
# ╠═6c8a9e97-0a60-4e42-ad7f-78c0dece2b30
# ╠═801a038b-7485-40b6-b7a5-92b0b3c11c0e
# ╠═62adbc28-cec2-4da8-9842-2153563497d0
# ╟─3bc2d36a-6efc-4631-baae-0d3e62bf7e3f
# ╟─ad83a851-e424-4e8f-bae9-38b9283db2fa
# ╟─b580ef9b-7d53-4bf0-8511-213920d59ee2
# ╟─090a949f-c18f-4e07-8c94-21a71bb0ceed
# ╠═85e81822-2ea6-4d05-b879-1791bff255b9
# ╟─a2a109e3-ea74-4e74-91d0-6ea12ddf28d4
# ╠═27a4c4f5-0b2b-4d23-9bb0-a1a8ff6cfc4d
# ╟─e0a8a397-4500-47de-8be3-49d3174648b1
# ╠═7955722d-fc35-4cce-8bf2-ee40ee7fbc82
# ╟─0590bc0b-45d6-4a8e-8078-af057984523f
# ╠═f7b09e96-cf30-4bd0-841b-46621df693ab
# ╟─e7947165-b244-4ff6-bdf3-61d03aefe696
# ╠═adc0aa15-e0a9-4549-9609-967a9dc78b78
# ╟─7a1effe4-776c-4143-8085-f98a213a2cc3
# ╟─f85bd6be-8984-49f0-a36b-7aa7bb937c54
# ╟─c14d45b8-3af1-4e05-b10a-081ccd62a34d
# ╠═04f89fb3-f856-4206-9871-4d92f2816332
# ╟─a39eeec9-d816-4d35-90ae-e602a24bb056
# ╠═fcad8ede-8295-46d6-a637-728b3986a761
# ╟─b3f08d61-187d-4f99-b95f-508ff1b3d52d
# ╟─08ba4614-e4fa-4203-9e28-c5c4bd9fcdd5
# ╠═e4087441-7d0d-4716-b9fa-eef45e9f3b87
# ╟─87b3796d-16c5-4fe9-a32e-d8bf7c754254
# ╠═2a66417b-79a2-4909-bd89-5650c49b9411
# ╟─35956a99-5729-4bf7-bf09-b2232fc958da
# ╠═cc31567f-631b-41b9-aa2d-3408ca951bcf
# ╟─82f8c50c-511d-4399-ac87-0890b678c6c1
# ╠═ec08adb2-561c-4f17-a0b5-019d7b1f1098
# ╟─c55ec157-6ec8-45db-9a54-0f39aceb9607
# ╠═67668918-624e-4adf-be5d-6fdbc62555f6
# ╟─832f984a-176d-4187-8476-5026d29a63d8
# ╠═5f9a9838-099d-457f-9402-915f42d0cd33
# ╟─ef3ca25d-c110-452d-877a-6304ae6cd5a7
# ╟─f8dbab24-55aa-466f-a978-5f8db93b4b93
# ╠═381d3a54-0e4f-4203-874f-4a4716d481d4
# ╟─fd7d89f3-9ab8-4b17-b26f-571bad862718
# ╠═96de4f0d-c2d1-4fad-bded-a314cbca848b
# ╟─dedb1ac1-cb54-416f-a092-12a50080c69f
# ╠═a912f449-6c9d-4595-863f-a2ba523bad95
# ╟─089707c3-1df6-4799-a87b-0d4872a0267d
# ╟─381a714d-699e-45c6-909a-02689b2a7e6b
# ╠═a1795123-876b-4bd3-ac3a-711367b500d1
# ╟─6b46d116-0da3-428a-b25f-221622dfa7e0
# ╠═710fe13c-bd43-4d1e-867e-72330722ac1e
# ╟─2ee20493-2dc2-4bd4-8ec3-78032f1fcc2d
# ╠═2943b515-fb4b-45a2-88dd-d0267ec95b09
# ╟─0274fbd8-c025-4e36-bbb0-f65f21b962c7
# ╠═34c4d746-45f7-43c4-a0f0-6e570373a35d
# ╟─060a5c5b-9834-42e2-98e9-ac23e7403b60
# ╠═c61ae8b6-2324-43fa-a5ae-274a919af559
# ╟─69983b90-b6ef-4732-a4af-a772bb1364e5
# ╠═16211bc6-cd97-43c0-8faf-25bb490930b6
# ╟─461eada9-5f9c-4439-8b65-226382b6d148
# ╠═0994d790-0fe8-46ab-bba2-a0a4ea466f64
# ╟─8bb5dafd-6c26-4634-8be4-c2963efba056
# ╠═26c77206-fe3c-45c5-a627-f2b346c66be7
# ╟─527e55a8-1744-4c40-a412-054bed6d77f0
# ╠═ceae9947-ed0f-4975-a944-cc938cf22dde
# ╟─6e767e36-f54e-4472-881a-ab04ad9c9d09
# ╠═870c9d65-2229-40b7-9959-59afecb4f3bf
# ╟─93ca3970-bebd-4192-bc62-6ca08d58a6c2
# ╠═7e4114d6-932a-4def-9bb5-b534072ea513
# ╟─1bf8bc58-6848-4783-a614-7ac11646de92
# ╠═8814c42b-16c7-466a-8ad8-58ee1c0921b2
# ╟─25158ea5-1930-4fee-aab4-490b475d8635
# ╠═ed9472a0-09c9-4d1c-8f1d-1b7d87b97640
# ╟─6e213afb-fd4b-4ab9-8986-a319bee6be0c
# ╠═3af06996-bd43-4783-a5e0-3c11c23ca280
# ╟─93400573-ad1f-408d-b514-854998e7836a
# ╠═008cc499-fce9-43e7-b969-2419b847d24c
# ╟─14c3b2f1-14b4-415c-85c6-922361423826
