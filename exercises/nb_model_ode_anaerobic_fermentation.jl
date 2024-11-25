### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 2552c020-1e29-451e-9f59-c4bde047faad
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 84e21b44-a1b0-11ef-014d-c58a169e3de3
using Markdown

# ╔═╡ 3cd7556f-73b8-4176-a779-ad38909f464d
using InteractiveUtils

# ╔═╡ 896b4151-e26b-40ee-bfb9-c56dfc4e7048
using Catalyst

# ╔═╡ aaf6da21-60cc-478c-b447-3f33aa375240
using DifferentialEquations, Plots

# ╔═╡ bb10266a-1f6c-4fda-b276-6a9cf3a86e90
md"""
### Exercise - Anaerobic fermentation
"""

# ╔═╡ c3f85244-e0a8-4808-aa2d-15ab6bbb1b26
md"""
#### Part 1
"""

# ╔═╡ 2e7c2eab-a3bc-4574-a0fd-0a45b59b803b
md"""
An operator at a beverage factory would like to model the anaerobic fermentation that occurs in one of his reactors. After a literature review, he finds that sucrose ($S$) is converted to ethanol ($E$) via glucose ($G$) under the action of an enzyme called invertase ($I$) from yeast using the following reaction stoichiometry (all components have the unit $mol\;L^{-1}$):

$$C_{12}H_{22}O_{11} + I \xrightarrow{r_1} 2C_{6}H_{12}O_{6}$$

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
##### Implementation of the system
"

# ╔═╡ bfd7fc6f-57f6-40d8-8618-d514d1e5d9ad
md"""
Create a *reaction network object* model for the aforementioned problem in order to simulate the evolution of $S$, $I$, $G$, $E$ and $CO_2$ during $720\;min$ ($=12\;h$). Name it `anaerobic_fermentation1`.

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
    k1, S + I --> 2G
	# mm(K, k2, E), 2G --> 4E + 4CO2
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
tspan1=(0.0, 720.0)

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
oprob1 = ODEProblem(anaerobic_fermentation1, u01, tspan1, params1)

# ╔═╡ 9a517512-a7c0-4a92-a104-b19a8bf786c7
md"""
Solve the ODE problem. Use `Tsit5()` and `saveat=0.5`. Store the solution in `osol1`:
"""

# ╔═╡ 69f9ce3a-0422-499a-8492-ea70de43d587
# osol1 = missing                    # Uncomment and complete the instruction
osol1 = solve(oprob1, Tsit5(), saveat=0.5)

# ╔═╡ f1c134c9-9127-4ed4-84a9-e1d218a82295
md"""
Plot the results.
"""

# ╔═╡ e8eeefd3-f83a-4154-945d-568ac629cd94
# missing                      # Uncomment and complete the instruction
plot(osol1, linewidth=2)

# ╔═╡ 622f65d3-0592-419f-8101-bc6f8b8ea5ed
md"""
Interprete the results. Try to come up with an answer to the following questions:

1. Why is the concentration of invertase ($I$) becoming zero, and the concentration of sucrose ($S$) becoming a fixed non-zero value?
    - Answer: missing
2. Try to explain the peak in the glucose ($G$) concentration.
    - Answer: missing
3. Why is the difference in ethanol ($E$) and $CO_2$ concentration constant?
    - Answer: missing
"""

# ╔═╡ 939716cd-08c9-44a0-9217-146f9e312d35
md"""
#### Part 2
"""

# ╔═╡ 194a596b-22de-412e-8ead-94a8c02c98c6
md"""
Additionally, sucrose, invertase and glucose are added at a flow rate $Q_{in},[L\; min^{-1}]$ and respective concentrations $S_{in}, I_{in}$ and $G_{in}$. The same flow rate is removed from the reactor. Now the volume $V$ of the reactor will matter.

The additional parameter values are summarised in the following table:

|  $Q_{in}$ |    $V$    | $S_{in}$  | $I_{in}$  | $G_{in}$ |
|:---------:|:---------:|:---------:|:---------:|:--------:|
|   $1.00$  |   $100$   |  $0.12$   |  $0.08$   |  $0.05$  |
"""

# ╔═╡ 76b1d8ba-cc55-48f7-bee4-09f899582e62
md"""
Make a copy of the content of the previous *reaction network object* and complement it with the new information. Name it `anaerobic_fermentation2`.
"""

# ╔═╡ 4165677e-a483-404c-a8c5-6cd85e5a8090
# Uncomment and complete the instruction
# anaerobic_fermentation2 = @reaction_network begin
# 	@species missing
#     missing
#     ...
#     missing
# end
anaerobic_fermentation2 = @reaction_network begin
	@species S(t)=0.04 I(t)=0.02 G(t)=0.0 E(t)=0.01 CO2(t)=0.0
    k1, S + I --> 2G
	# mm(K, k2, E), 2G --> 4E + 4CO2
	mm(E, k2, K), 2G --> 4E + 4CO2
	Q/V, (S, I, G, E, CO2) --> (0, 0, 0, 0, 0)
	Q/V*Sin, 0 --> S
	Q/V*Iin, 0 --> I
	Q/V*Gin, 0 --> G
end

# ╔═╡ 3da7ee15-6d30-42bb-be03-6f9129febc60
md"""
Convert the system to a symbolic differential equation model and inspect your differential equation.
"""

# ╔═╡ 58b8545f-1699-45b4-91c0-86722668a428
# osys2 = missing               # Uncomment and complete the instruction
osys2 = convert(ODESystem, anaerobic_fermentation2)

# ╔═╡ f6d58bac-7256-4b57-836d-6d243bf292c0
md"""
Make a copy of `u01` and rename it to `u02` with the initial conditions:
"""

# ╔═╡ 75eb088a-33fd-47d7-892a-da934ec9896a
# u02 = missing                # Uncomment and complete the instruction
u02 = [:S => 0.04, :I => 0.02, :G => 0.0, :E => 0.01, :CO2 => 0.0]

# ╔═╡ b4c4fc78-b6e4-40db-9948-6be83ee57d87
md"""
Make a copy of `tspan1` and rename it to `tspan2`:
"""

# ╔═╡ 4a0e9abf-99e7-4fe6-82ae-aa958ab51dd0
# tspan2 = missing                # Uncomment and complete the instruction
tspan2=(0.0, 720.0)

# ╔═╡ 0ec39752-6654-452b-98db-eadd8c02b27f
md"""
Make a copy of `params1`, rename it to `params2` and supplement it with the new parameter values:
"""

# ╔═╡ 3aea0dbe-1cb3-4b19-bb38-53e58db153dd
# param2 = missing                # Uncomment and complete the instruction
params2=[:Q => 1, :V => 100.0, :Sin => 0.12, :Iin => 0.08, :Gin => 0.05, :k1 => 0.4, :k2 => 0.65, :K => 0.5]

# ╔═╡ 49984e8e-a194-464a-994f-501342800026
md"""
Create the ODE problem and store it in `oprob2`:
"""

# ╔═╡ 18a3513f-8ecc-4c2f-bb39-5d0cb39a3d92
oprob2 = ODEProblem(anaerobic_fermentation2, u02, tspan2, params2)

# ╔═╡ 5318708e-3410-4675-8162-827c0c6e039c
md"""
Solve the ODE problem. Use `Tsit5()` and `saveat=0.5`. Store the solution in `osol2`:
"""

# ╔═╡ e5476343-8376-4a9a-8286-f69e0d023939
osol2 = solve(oprob2, Tsit5(), saveat=0.5)

# ╔═╡ 29a3cd1b-fd8f-4e5c-9ec7-a69fe328d824
md"""
Plot the results.
"""

# ╔═╡ 0c179e6b-bcaf-4e89-8d04-59f24b456127
# plot(osol, linewidth=2; idxs=:E)
plot(osol2, linewidth=2)

# ╔═╡ 5f865fec-09d5-4358-be57-0d568740968d
md"""
Check out the last concentrations (at the end time) for each of the species.

Tips:
- You can access the last value of $S$ with `osol2[:S][end]`.
- You can put everything on one line separating the values with comma's.
"""

# ╔═╡ 053cca09-0e31-44c0-89db-26a581816744
(osol2[:S][end], osol2[:I][end], osol2[:G][end], osol2[:E][end], osol2[:CO2][end])

# ╔═╡ 83cfd592-896d-42f5-ae2a-b1d0a9f14aac
md"""
Create a vector named `u_guess` in the same way as `u02`, but now with the end values of the species.
"""

# ╔═╡ f5a8f97b-7cbb-4eb4-b5ff-82c184b7de85
u_guess2 = [:S=>osol2[:S][end], :I=>osol2[:I][end], :G=>osol2[:G][end], :E=>osol2[:E][end], :CO2=>osol2[:CO2][end]]

# ╔═╡ ce67b1f7-a82d-4fa3-a22b-f8ee3d2f4af0
md"""
Calculate the steady-state values of the species:
"""

# ╔═╡ 3a014879-883b-4cab-a7ac-0c3a708f5da6
Seq2, Ieq2, Geq2, Eeq2, CO2eq2 = solve(SteadyStateProblem(ODEProblem(anaerobic_fermentation2, u_guess2, tspan2, params2)))

# ╔═╡ d7bde85c-05f1-4794-9e91-c4e5e120f876
md"""
Check ou the steady states:
"""

# ╔═╡ 876683d7-b3cd-4a06-94b2-d38a1d4a34b6
Seq2, Ieq2, Geq2, Eeq2, CO2eq2

# ╔═╡ Cell order:
# ╠═84e21b44-a1b0-11ef-014d-c58a169e3de3
# ╠═3cd7556f-73b8-4176-a779-ad38909f464d
# ╠═2552c020-1e29-451e-9f59-c4bde047faad
# ╠═896b4151-e26b-40ee-bfb9-c56dfc4e7048
# ╠═aaf6da21-60cc-478c-b447-3f33aa375240
# ╠═bb10266a-1f6c-4fda-b276-6a9cf3a86e90
# ╠═c3f85244-e0a8-4808-aa2d-15ab6bbb1b26
# ╠═2e7c2eab-a3bc-4574-a0fd-0a45b59b803b
# ╠═35f2015d-9c28-4bc8-83cd-c185757db0dc
# ╠═bfd7fc6f-57f6-40d8-8618-d514d1e5d9ad
# ╠═53c32175-4298-4450-ae85-132ea0cd6a9b
# ╠═485c2b25-f281-4ee9-bcb8-7ca010205930
# ╠═86200a5e-55ce-4d58-87ba-fe30f65c9120
# ╠═579318a9-bad9-4392-9420-3f5474ccfff8
# ╠═bc56d0fb-db19-4f63-9f34-617308da6c8a
# ╠═5983b9d3-5119-47ac-9a0c-8a2755ad23e6
# ╠═c822fcbc-e8d3-4952-93bc-38a41f1786a7
# ╠═eb7fda8e-d1db-4bcc-91d6-964a6dcbbe90
# ╠═124571f0-fe33-44a2-9ab4-f6a418ac4f51
# ╠═703d2eb0-8504-496f-9abe-8c6911a0fdbc
# ╠═9dd45f30-9076-45d5-827b-02ae2d73ef97
# ╠═64f85604-9562-4a97-bedd-ea5d28ada096
# ╠═04816732-3419-40e3-a927-d6e102949553
# ╠═9a517512-a7c0-4a92-a104-b19a8bf786c7
# ╠═69f9ce3a-0422-499a-8492-ea70de43d587
# ╠═f1c134c9-9127-4ed4-84a9-e1d218a82295
# ╠═e8eeefd3-f83a-4154-945d-568ac629cd94
# ╠═622f65d3-0592-419f-8101-bc6f8b8ea5ed
# ╠═939716cd-08c9-44a0-9217-146f9e312d35
# ╠═194a596b-22de-412e-8ead-94a8c02c98c6
# ╠═76b1d8ba-cc55-48f7-bee4-09f899582e62
# ╠═4165677e-a483-404c-a8c5-6cd85e5a8090
# ╠═3da7ee15-6d30-42bb-be03-6f9129febc60
# ╠═58b8545f-1699-45b4-91c0-86722668a428
# ╠═f6d58bac-7256-4b57-836d-6d243bf292c0
# ╠═75eb088a-33fd-47d7-892a-da934ec9896a
# ╠═b4c4fc78-b6e4-40db-9948-6be83ee57d87
# ╠═4a0e9abf-99e7-4fe6-82ae-aa958ab51dd0
# ╠═0ec39752-6654-452b-98db-eadd8c02b27f
# ╠═3aea0dbe-1cb3-4b19-bb38-53e58db153dd
# ╠═49984e8e-a194-464a-994f-501342800026
# ╠═18a3513f-8ecc-4c2f-bb39-5d0cb39a3d92
# ╠═5318708e-3410-4675-8162-827c0c6e039c
# ╠═e5476343-8376-4a9a-8286-f69e0d023939
# ╠═29a3cd1b-fd8f-4e5c-9ec7-a69fe328d824
# ╠═0c179e6b-bcaf-4e89-8d04-59f24b456127
# ╠═5f865fec-09d5-4358-be57-0d568740968d
# ╠═053cca09-0e31-44c0-89db-26a581816744
# ╠═83cfd592-896d-42f5-ae2a-b1d0a9f14aac
# ╠═f5a8f97b-7cbb-4eb4-b5ff-82c184b7de85
# ╠═ce67b1f7-a82d-4fa3-a22b-f8ee3d2f4af0
# ╠═3a014879-883b-4cab-a7ac-0c3a708f5da6
# ╠═d7bde85c-05f1-4794-9e91-c4e5e120f876
# ╠═876683d7-b3cd-4a06-94b2-d38a1d4a34b6
