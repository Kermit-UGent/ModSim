### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 2552c020-1e29-451e-9f59-c4bde047faad
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ 84e21b44-a1b0-11ef-014d-c58a169e3de3
using Markdown

# ╔═╡ 3cd7556f-73b8-4176-a779-ad38909f464d
using InteractiveUtils

# ╔═╡ 896b4151-e26b-40ee-bfb9-c56dfc4e7048
using Catalyst

# ╔═╡ aaf6da21-60cc-478c-b447-3f33aa375240
using DifferentialEquations, Plots

# ╔═╡ e119ed70-fd81-4978-92d9-696cad71125f
using PlutoUI; TableOfContents()

# ╔═╡ bb10266a-1f6c-4fda-b276-6a9cf3a86e90
md"""
# Exercise - Anaerobic fermentation
"""

# ╔═╡ c3f85244-e0a8-4808-aa2d-15ab6bbb1b26
md"""
## Part 1
"""

# ╔═╡ 2e7c2eab-a3bc-4574-a0fd-0a45b59b803b
md"""
An operator at a beverage factory would like to model the anaerobic fermentation that occurs in one of his reactors. After a literature review, he finds that sucrose ($S$) is converted to ethanol ($E$) via glucose ($G$) under the action of an enzyme called invertase ($I$) from yeast using the following reaction stoichiometry (all components have the unit $mol\;L^{-1}$):

$$C_{12}H_{22}O_{11} + I \xrightarrow{r_1} 2C_{6}H_{12}O_{6} + I$$

$$C_{6}H_{12}O_{6} \xrightarrow{r_2} 2C_{2}H_{5}OH + 2CO_2$$

The operator knows that both reactions are carried out in isothermal conditions in a reactor with volume $V$ $[L]$.

The operator knows from literature that the reaction rate $r_1$ is first-order in both sucrose and invertase and with specific reaction rate $k_1$. The reaction rate $r_2$ is second-order with respect to glucose with specific reaction rate $k_2$. Additionally, the reaction is inhibited by ethanol itself according to $\cfrac{K}{E+K}$ where $K$ represents the ethanol concentration at which $r_2$ achieves half of its maximum reaction rate.

The initial concentrations and the parameter values are summarised in the following tables:

|   $S_0$   |   $I_0$   |   $G_0$   |   $E_0$   | $CO_{2,0}$ |
|:---------:|:---------:|:---------:|:---------:|:----------:|
|   $0.04$  |   $0.02$  |  $0.00$   |  $0.01$   |    $0.00$  |

|   $k_1$   |   $k_2$   |    $K$    |
|:---------:|:---------:|:---------:|
|  $0.40$   |  $0.65$   |   $0.50$  |

"""

# ╔═╡ 35f2015d-9c28-4bc8-83cd-c185757db0dc
md"
#### Implementation of the system
"

# ╔═╡ bfd7fc6f-57f6-40d8-8618-d514d1e5d9ad
md"""
Create a *reaction network object* model for the aforementioned problem in order to simulate the evolution of $S$, $I$, $G$, $E$ and $CO_2$ during $1440\;min$ ($=24\;h$). Name it `anaerobic_fermentation1`.

Tips:
- For the reaction with reaction rate $r_2$, in order to have a second-order reaction with respect to glucose, you need to double the stoichiometric coefficients, i.e., you reaction should be $2G \rightarrow 4E + 4CO_2$.
- For the inhibition factor $\cfrac{K}{E+K}$ you can use the function `mmr(..., ..., ...)`!
"""

# ╔═╡ 53c32175-4298-4450-ae85-132ea0cd6a9b
# Uncomment and complete the instruction
# anaerobic_fermentation1 = @reaction_network begin
# 	@species missing
# 	missing
# 	missing
# end
anaerobic_fermentation1 = @reaction_network begin
	@species S(t)=0.04 I(t)=0.02 G(t)=0.0 E(t)=0.01 CO2(t)=0.0
    k1, S + I --> 2G + I
	mmr(E, k2, K), 2G --> 4E + 4CO2
end

# ╔═╡ 485c2b25-f281-4ee9-bcb8-7ca010205930
md"""
Check out the species.
"""

# ╔═╡ 86200a5e-55ce-4d58-87ba-fe30f65c9120
# missing           # Uncomment and complete the instruction
species(anaerobic_fermentation1)

# ╔═╡ 579318a9-bad9-4392-9420-3f5474ccfff8
md"""
Convert the system to a symbolic differential equation model and inspect your differential equations.
"""

# ╔═╡ bc56d0fb-db19-4f63-9f34-617308da6c8a
# osys1 = missing         # Uncomment and complete the instruction
osys1 = convert(ODESystem, anaerobic_fermentation1)
# Same result with:
# osys1 = convert(ODESystem, anaerobic_fermentation1, combinatoric_ratelaws=true)

# ╔═╡ 5983b9d3-5119-47ac-9a0c-8a2755ad23e6
md"""
Initialize a vector `u01` with the initial conditions:
"""

# ╔═╡ c822fcbc-e8d3-4952-93bc-38a41f1786a7
# u01 = missing            # Uncomment and complete the instruction
u01 = [:S => 0.04, :I => 0.02, :G => 0.0, :E => 0.01, :CO2 => 0.0]

# ╔═╡ eb7fda8e-d1db-4bcc-91d6-964a6dcbbe90
md"""
Set the timespan for the simulation:
"""

# ╔═╡ 124571f0-fe33-44a2-9ab4-f6a418ac4f51
# tspan1 = missing          # Uncomment and complete the instruction
tspan1=(0.0, 1440.0)

# ╔═╡ 703d2eb0-8504-496f-9abe-8c6911a0fdbc
md"""
Initialize a vector `params1` with the parameter values:
"""

# ╔═╡ 9dd45f30-9076-45d5-827b-02ae2d73ef97
# params1 = missing          # Uncomment and complete the instruction
params1=[:k1 => 0.4, :k2 => 0.65, :K => 0.5]

# ╔═╡ 64f85604-9562-4a97-bedd-ea5d28ada096
md"""
Create the ODE problem and store it in `oprob1`:
"""

# ╔═╡ 04816732-3419-40e3-a927-d6e102949553
# oprob1 = missing             # Uncomment and complete the instruction
oprob1 = ODEProblem(anaerobic_fermentation1, u01, tspan1, params1);

# ╔═╡ 9a517512-a7c0-4a92-a104-b19a8bf786c7
md"""
Solve the ODE problem. Use `Tsit5()` and `saveat=0.5`. Store the solution in `osol1`:
"""

# ╔═╡ 69f9ce3a-0422-499a-8492-ea70de43d587
# osol1 = missing                    # Uncomment and complete the instruction
osol1 = solve(oprob1, Tsit5(), saveat=0.5)

# ╔═╡ f1c134c9-9127-4ed4-84a9-e1d218a82295
md"""
Plot the results. Use a line width of 2 (`linewidth=...`).
"""

# ╔═╡ e8eeefd3-f83a-4154-945d-568ac629cd94
# missing                      # Uncomment and complete the instruction
plot(osol1, linewidth=2)

# ╔═╡ 622f65d3-0592-419f-8101-bc6f8b8ea5ed
md"""
Interprete the results. Try to come up with an answer to the following questions:

1. Why is the concentration of invertase ($I$) constant zero, and the concentration of sucrose ($S$) becoming zero?
"""

# ╔═╡ 39e96a71-7483-4c85-b44d-76cafc3b22dc
md"- Answer: missing"
#=
1. Sucrose S is being converted to glucose G by means of the enzyme I. The enzyme I is functioning as a catalyst, so it is not being destroyed. On the otherhand sucrose S is converted to glucose G, so S is consumed to produce G until S is used up.
=#

# ╔═╡ 192ec8c6-ff60-4584-91f9-d257399bac76
md"""
2. Try to explain the peak in the glucose ($G$) concentration.
"""

# ╔═╡ e2681cdb-581d-41d7-93c6-b137ef54c083
md"- Answer: missing"
#=
2. G is being produced by the consumption of S and I. Initially there is no G in the reactor, so the production precedes. The convertion of G into E will start with some delay because it is proportional to G. Hence, when there is enough G in the reactor the convertion will take place and will compete with the production of G. Hence G will first increase and then decline.
=#

# ╔═╡ 9f5bf914-7ec6-48ba-8fed-2087f42b37d8
md"""
3. Why is the difference in ethanol ($E$) and $CO_2$ concentration constant?
"""

# ╔═╡ e50013b7-be37-4fa9-bdd6-2a902443f222
md"- Answer: missing"
#=
3. Both E and CO2 are being produced at exactly the same rate. Since the initial concentrations are 0.01 and 0.0 respectively, their difference will always be constant.
=#

# ╔═╡ 939716cd-08c9-44a0-9217-146f9e312d35
md"""
## Part 2
"""

# ╔═╡ 194a596b-22de-412e-8ead-94a8c02c98c6
md"""
Additionally, sucrose and glucose are added at a flow rate $Q_{in},[L\; min^{-1}]$ and respective concentrations $S_{in}$ and $G_{in}$. The same flow rate is removed from the reactor but the invertase $I$ stays in the reactor. Now the volume $V$ of the reactor will matter. Furthermore, the invertase enzyme degrades at a rate $d=0.003\;min^{-1}$.

The additional parameter values are summarised in the following table:

|  $Q_{in}$ |    $V$    | $S_{in}$  | $G_{in}$  | $d$ |
|:---------:|:---------:|:---------:|:---------:|:--------:|
|   $1.00$  |   $100$   |  $0.12$   |  $0.05$   |  $0.003$  |
"""

# ╔═╡ 76b1d8ba-cc55-48f7-bee4-09f899582e62
md"""
Make a copy of the content of the previous *reaction network object* and complement it with the new information. Name it `anaerobic_fermentation2`.
"""

# ╔═╡ 4165677e-a483-404c-a8c5-6cd85e5a8090
# Uncomment and complete the instruction
# anaerobic_fermentation2 = @reaction_network begin
#     @species missing
# 	  @parameters missing
#     missing
#     ...
#     missing
# end
anaerobic_fermentation2 = @reaction_network begin
    @species S(t)=0.04 I(t)=0.02 G(t)=0.0 E(t)=0.01 CO2(t)=0.0
	@parameters d=0.003
    k1, S + I --> 2G + I
    mm(E, k2, K), 2G --> 4E + 4CO2
    Q/V, (S, G, E, CO2) --> (0, 0, 0, 0)
	d, I --> 0
    Q/V*Sin, 0 --> S
    Q/V*Gin, 0 --> G
end

# ╔═╡ 3da7ee15-6d30-42bb-be03-6f9129febc60
md"""
Convert the system to a symbolic differential equation model and inspect your differential equations.
"""

# ╔═╡ 58b8545f-1699-45b4-91c0-86722668a428
# osys2 = missing               # Uncomment and complete the instruction
osys2 = convert(ODESystem, anaerobic_fermentation2)

# ╔═╡ f6d58bac-7256-4b57-836d-6d243bf292c0
md"""
Make an exact copy of `u01` and rename it to `u02` with the initial conditions:
"""

# ╔═╡ 75eb088a-33fd-47d7-892a-da934ec9896a
# u02 = missing                # Uncomment and complete the instruction
u02 = [:S => 0.04, :I => 0.02, :G => 0.0, :E => 0.01, :CO2 => 0.0]

# ╔═╡ b4c4fc78-b6e4-40db-9948-6be83ee57d87
md"""
Make an exact copy of `tspan1` and rename it to `tspan2`:
"""

# ╔═╡ 4a0e9abf-99e7-4fe6-82ae-aa958ab51dd0
# tspan2 = missing                # Uncomment and complete the instruction
tspan2=(0.0, 1440.0)

# ╔═╡ 0ec39752-6654-452b-98db-eadd8c02b27f
md"""
Make a copy of `params1`, rename it to `params2` and supplement it with the new parameter values:
"""

# ╔═╡ 3aea0dbe-1cb3-4b19-bb38-53e58db153dd
# param2 = missing                # Uncomment and complete the instruction
params2=[:Q => 1, :V => 100.0, :Sin => 0.12, :Gin => 0.05, :k1 => 0.4, :k2 => 0.65, :K => 0.5]

# ╔═╡ 49984e8e-a194-464a-994f-501342800026
md"""
Create the ODE problem and store it in `oprob2`:
"""

# ╔═╡ 18a3513f-8ecc-4c2f-bb39-5d0cb39a3d92
# oprob2 = missing                 # Uncomment and complete the instruction
oprob2 = ODEProblem(anaerobic_fermentation2, u02, tspan2, params2);

# ╔═╡ 5318708e-3410-4675-8162-827c0c6e039c
md"""
Solve the ODE problem. Use `Tsit5()` and `saveat=0.5`. Store the solution in `osol2`:
"""

# ╔═╡ e5476343-8376-4a9a-8286-f69e0d023939
# osol2 = missing                    # Uncomment and complete the instruction
osol2 = solve(oprob2, Tsit5(), saveat=0.5)

# ╔═╡ 29a3cd1b-fd8f-4e5c-9ec7-a69fe328d824
md"""
Plot the results. Use a line width of 2 (`linewidth=...`). If you only want to see the curves for, e.g., $E$, $S$ and $G$, you can use the option `idxs=[:E, :S, :G]` in the `plot` command.
"""

# ╔═╡ 0c179e6b-bcaf-4e89-8d04-59f24b456127
# missing                # Uncomment and complete the instruction
plot(osol2, linewidth=2)
# plot(osol2, linewidth=2, idxs=[:E, :S, :G])

# ╔═╡ e596f49e-1b3a-422e-be9a-f72133f042d5
md"""
Interprete the results.
"""

# ╔═╡ 3505c10d-270c-49fc-a83a-3753974d4d5a
md"- Answer: missing"

# ╔═╡ 5f865fec-09d5-4358-be57-0d568740968d
md"""
Check out the last concentrations (at the end time) for each of the species.

Tips:
- You can see the last values of all species with `osol2.u[end]`
- If later you need all last values separately, you can access the last value of $S$ with `osol2[:S][end]` and then you can put everything on one line separating the values with comma's.
"""

# ╔═╡ acc8acb3-f014-4108-80d7-bff893ce07d1
# missing                 # Uncomment and complete the instruction
osol2.u[end]

# ╔═╡ 053cca09-0e31-44c0-89db-26a581816744
# osol2[:S][end], ..., ..., ..., ...)  # Uncomment and complete the instruction
(osol2[:S][end], osol2[:I][end], osol2[:G][end], osol2[:E][end], osol2[:CO2][end])

# ╔═╡ 83cfd592-896d-42f5-ae2a-b1d0a9f14aac
md"""
Create a vector named `u_guess` in the same way as `u02`, but now with the end values of the species.
"""

# ╔═╡ f5a8f97b-7cbb-4eb4-b5ff-82c184b7de85
# u_guess2 = missing              # Uncomment and complete the instruction
u_guess2 = [:S=>osol2[:S][end], :I=>osol2[:I][end], :G=>osol2[:G][end], :E=>osol2[:E][end], :CO2=>osol2[:CO2][end]]

# ╔═╡ ce67b1f7-a82d-4fa3-a22b-f8ee3d2f4af0
md"""
Calculate the steady-state values of the species:
"""

# ╔═╡ 3a014879-883b-4cab-a7ac-0c3a708f5da6
# Sw2, Iw2, Gw2, Ew2, CO2w2 = missing # Uncomment and complete the instruction
Sw2, Iw2, Gw2, Ew2, CO2w2 = solve(SteadyStateProblem(ODEProblem(anaerobic_fermentation2, u_guess2, tspan2, params2)))

# ╔═╡ d7bde85c-05f1-4794-9e91-c4e5e120f876
md"""
Check ou the steady states:
"""

# ╔═╡ 876683d7-b3cd-4a06-94b2-d38a1d4a34b6
# missing
Sw2, Iw2, Gw2, Ew2, CO2w2

# ╔═╡ 55d303ec-29cb-4e8e-ab52-808d881ae0f2
md"""
## Part 3
"""

# ╔═╡ e0ef008c-9276-41e5-a85a-649a276c4666
md"""
We now want to keep a relatively high production of ethanol. Therefore, if the invertase decreases to $0.008$, then the invertase is instantaneously renewed to the initial concentration of $0.02$. Apply the change in the invertase concentration using a continuous event.
"""

# ╔═╡ 81db6771-ee45-42c8-ab23-fe1960370b86
md"""
Create the correct condition.
"""

# ╔═╡ 166a21f0-46a5-4d7a-9f54-788d7a51a487
# condition3 = missing    # Uncomment and complete the instruction
condition3 = [anaerobic_fermentation2.I ~ 0.008] => [anaerobic_fermentation2.I ~ 0.02]

# ╔═╡ 5e6c5431-6fad-4272-a539-5e3af83fa226
md"""
Include the condition into the *reaction network model*.
"""

# ╔═╡ 98be4b39-4097-4cfa-886e-171558da0df6
# Uncomment and complete the instruction
# @named anaerobic_fermentation3_c = missing
@named anaerobic_fermentation3_c = ReactionSystem(equations(anaerobic_fermentation2), continuous_events=condition3)

# ╔═╡ 66568f6a-b874-42cb-834d-ab9d4df815d1
md"""
Complete the *reaction network model*.
"""

# ╔═╡ 66441227-c202-48da-a6ed-58e25ffde3ce
# Uncomment and complete the instruction
# anaerobic_fermentation3_c_com = missing
anaerobic_fermentation3_c_com = complete(anaerobic_fermentation3_c)

# ╔═╡ 71b4bd6e-b906-4d88-a9a3-18700528ea12
md"""
Create a new ODE problem.
"""

# ╔═╡ 16d0d671-464c-4de3-bd73-593e1a73168b
# oprob3 = missing        # Uncomment and complete the instruction
oprob3 = ODEProblem(anaerobic_fermentation3_c_com, u02, tspan2, params2);

# ╔═╡ 66015536-5f23-45b5-8d57-f1d27d604f8c
md"""
Solve the new ODE problem. Make a `deepcopy`, use `Tsit5()` and `saveat=0.5`.
"""

# ╔═╡ 7f0777a2-6eb9-4c6d-9fb7-25a1f32f291d
# osol3 = missing
osol3 = solve(deepcopy(oprob3), Tsit5(), saveat=0.5)

# ╔═╡ 124d578f-42ea-488e-af65-312793f104ca
md"""
Plot the results.
"""

# ╔═╡ 9ad3aff4-76ea-452f-bd7b-fc379725df70
plot(osol3)

# ╔═╡ e07e1667-2f2c-4086-9818-90dff54296de
md"""
Interprete the results.
"""

# ╔═╡ 91af34c0-9b38-4a62-8387-5e3c60b4b145
md"- Answer: missing"

# ╔═╡ Cell order:
# ╠═84e21b44-a1b0-11ef-014d-c58a169e3de3
# ╠═3cd7556f-73b8-4176-a779-ad38909f464d
# ╠═2552c020-1e29-451e-9f59-c4bde047faad
# ╠═896b4151-e26b-40ee-bfb9-c56dfc4e7048
# ╠═aaf6da21-60cc-478c-b447-3f33aa375240
# ╠═e119ed70-fd81-4978-92d9-696cad71125f
# ╟─bb10266a-1f6c-4fda-b276-6a9cf3a86e90
# ╟─c3f85244-e0a8-4808-aa2d-15ab6bbb1b26
# ╟─2e7c2eab-a3bc-4574-a0fd-0a45b59b803b
# ╟─35f2015d-9c28-4bc8-83cd-c185757db0dc
# ╟─bfd7fc6f-57f6-40d8-8618-d514d1e5d9ad
# ╠═53c32175-4298-4450-ae85-132ea0cd6a9b
# ╟─485c2b25-f281-4ee9-bcb8-7ca010205930
# ╠═86200a5e-55ce-4d58-87ba-fe30f65c9120
# ╟─579318a9-bad9-4392-9420-3f5474ccfff8
# ╠═bc56d0fb-db19-4f63-9f34-617308da6c8a
# ╟─5983b9d3-5119-47ac-9a0c-8a2755ad23e6
# ╠═c822fcbc-e8d3-4952-93bc-38a41f1786a7
# ╟─eb7fda8e-d1db-4bcc-91d6-964a6dcbbe90
# ╠═124571f0-fe33-44a2-9ab4-f6a418ac4f51
# ╟─703d2eb0-8504-496f-9abe-8c6911a0fdbc
# ╠═9dd45f30-9076-45d5-827b-02ae2d73ef97
# ╟─64f85604-9562-4a97-bedd-ea5d28ada096
# ╠═04816732-3419-40e3-a927-d6e102949553
# ╟─9a517512-a7c0-4a92-a104-b19a8bf786c7
# ╠═69f9ce3a-0422-499a-8492-ea70de43d587
# ╟─f1c134c9-9127-4ed4-84a9-e1d218a82295
# ╠═e8eeefd3-f83a-4154-945d-568ac629cd94
# ╟─622f65d3-0592-419f-8101-bc6f8b8ea5ed
# ╠═39e96a71-7483-4c85-b44d-76cafc3b22dc
# ╟─192ec8c6-ff60-4584-91f9-d257399bac76
# ╠═e2681cdb-581d-41d7-93c6-b137ef54c083
# ╟─9f5bf914-7ec6-48ba-8fed-2087f42b37d8
# ╠═e50013b7-be37-4fa9-bdd6-2a902443f222
# ╟─939716cd-08c9-44a0-9217-146f9e312d35
# ╟─194a596b-22de-412e-8ead-94a8c02c98c6
# ╟─76b1d8ba-cc55-48f7-bee4-09f899582e62
# ╠═4165677e-a483-404c-a8c5-6cd85e5a8090
# ╟─3da7ee15-6d30-42bb-be03-6f9129febc60
# ╠═58b8545f-1699-45b4-91c0-86722668a428
# ╟─f6d58bac-7256-4b57-836d-6d243bf292c0
# ╠═75eb088a-33fd-47d7-892a-da934ec9896a
# ╟─b4c4fc78-b6e4-40db-9948-6be83ee57d87
# ╠═4a0e9abf-99e7-4fe6-82ae-aa958ab51dd0
# ╟─0ec39752-6654-452b-98db-eadd8c02b27f
# ╠═3aea0dbe-1cb3-4b19-bb38-53e58db153dd
# ╟─49984e8e-a194-464a-994f-501342800026
# ╠═18a3513f-8ecc-4c2f-bb39-5d0cb39a3d92
# ╟─5318708e-3410-4675-8162-827c0c6e039c
# ╠═e5476343-8376-4a9a-8286-f69e0d023939
# ╟─29a3cd1b-fd8f-4e5c-9ec7-a69fe328d824
# ╠═0c179e6b-bcaf-4e89-8d04-59f24b456127
# ╟─e596f49e-1b3a-422e-be9a-f72133f042d5
# ╟─3505c10d-270c-49fc-a83a-3753974d4d5a
# ╟─5f865fec-09d5-4358-be57-0d568740968d
# ╠═acc8acb3-f014-4108-80d7-bff893ce07d1
# ╠═053cca09-0e31-44c0-89db-26a581816744
# ╟─83cfd592-896d-42f5-ae2a-b1d0a9f14aac
# ╠═f5a8f97b-7cbb-4eb4-b5ff-82c184b7de85
# ╟─ce67b1f7-a82d-4fa3-a22b-f8ee3d2f4af0
# ╠═3a014879-883b-4cab-a7ac-0c3a708f5da6
# ╟─d7bde85c-05f1-4794-9e91-c4e5e120f876
# ╠═876683d7-b3cd-4a06-94b2-d38a1d4a34b6
# ╟─55d303ec-29cb-4e8e-ab52-808d881ae0f2
# ╟─e0ef008c-9276-41e5-a85a-649a276c4666
# ╟─81db6771-ee45-42c8-ab23-fe1960370b86
# ╠═166a21f0-46a5-4d7a-9f54-788d7a51a487
# ╟─5e6c5431-6fad-4272-a539-5e3af83fa226
# ╠═98be4b39-4097-4cfa-886e-171558da0df6
# ╟─66568f6a-b874-42cb-834d-ab9d4df815d1
# ╠═66441227-c202-48da-a6ed-58e25ffde3ce
# ╟─71b4bd6e-b906-4d88-a9a3-18700528ea12
# ╠═16d0d671-464c-4de3-bd73-593e1a73168b
# ╟─66015536-5f23-45b5-8d57-f1d27d604f8c
# ╠═7f0777a2-6eb9-4c6d-9fb7-25a1f32f291d
# ╟─124d578f-42ea-488e-af65-312793f104ca
# ╠═9ad3aff4-76ea-452f-bd7b-fc379725df70
# ╟─e07e1667-2f2c-4086-9818-90dff54296de
# ╟─91af34c0-9b38-4a62-8387-5e3c60b4b145
