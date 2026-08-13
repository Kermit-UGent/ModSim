### A Pluto.jl notebook ###
# v0.20.13

#> [frontmatter]
#> order = "7"
#> title = "1. ODE_model_temp_reactor"
#> tags = ["exercises"]
#> layout = "layout.jlhtml"
#> description = "modeling temperature in a CSTR"

using Markdown
using InteractiveUtils

# ╔═╡ 7d61614c-9e15-11f0-1ce8-af091a3beb9a
# Running this yourself? Point this at your own environment —
# we advise one shared project in the parent folder: Pkg.activate("..")
using Pkg; Pkg.activate("../../pluto-deployment-environment")

# ╔═╡ 2e418e46-84fd-4922-9985-a92a412c1fb1
using StatsPlots, PlutoUI; TableOfContents()

# ╔═╡ d0a12f36-2592-4c3f-9581-6e41fbe5117f
using OrdinaryDiffEq, ModelingToolkit

# ╔═╡ 225d31d2-a974-4494-8813-c93484ce721e
using ModelingToolkit: t_nounits as t, D_nounits as D

# ╔═╡ 4bf49e49-7675-46eb-acfd-603dd68ac027
md"""
# Exercise: Temperature evolution in a set of reactors
"""

# ╔═╡ 0e4065c2-281c-439a-a69b-e4eca3b5ae41
solution(text) = Markdown.MD(Markdown.Admonition("hint", "Solution", [text]));

# ╔═╡ 7e2a3f8c-5f1c-4357-ae92-f4c1c677fde3
md"""
## Deriving the model equations
"""

# ╔═╡ 9e8c3a6b-6d69-4f9b-a8ea-de732de886ca
md"""
Consider two reactors in series with constant volumes of liquid $V_{1}$ and $V_{2}$ in a room at a constant ambient temperature $T_a$. A fluid with temperature $T_{in}$, adjustable by the control engineer, flows through both reactors at a constant flow rate $Q$. The temperature of the fluid in the first and second reactors are $T_1$ and $T_2$ respectively. The heat loss through the walls of the reactors is proportional to the wall areas $A_1$ (first reactor) and $A_2$ (second reactor). The heat transfer coefficient of both reactor walls is $\lambda$. In the second reactor a heater with an adjustable power $P_h$ delivers heat to the fluid. The density of the fluid is $\rho$ and its thermal heat capacity is $c_p$.
"""

# ╔═╡ 6a5cc6f2-b6b8-40fe-9b52-1a3e2d5d4fb7
# PlutoUI.LocalResource("fig/temperature_reactors.png")
md"""
![Temperature reators](https://users.ugent.be/~gvhaelew/fig/temperature_reactors.png)
"""

# ╔═╡ d1480cba-5e0f-4767-8164-b8ed5291e888
md"""
The following quantities are important for setting-up the model equations:
-  $c_p \rho V_1 T_1$: the amount of heat-energy in reactor 1
-  $c_p \rho V_2 T_2$: the amount of heat-energy in reactor 2
-  $c_p \rho Q \left( T_{in} - T_1 \right)$: rate of heat flow entering reactor 1
-  $c_p \rho Q \left( T_{1} - T_2 \right)$: rate of heat flow entering reactor 2
-  $\lambda A_1 \left(T_1 - T_a\right)$: rate of heat flow dissipating from reactor 1
-  $\lambda A_2 \left(T_2 - T_a\right)$: rate of heat flow dissipating from reactor 2
-  $P_h$: rate of heat supplied to reactor 2
"""

# ╔═╡ 826b1cda-3fcc-46e6-b186-d031929a569b
md"""
!!! task
	Derive the model equations for the fluid temperatures $T_1$ and $T_2$.
"""

# ╔═╡ 409e4f75-0cc7-4f2b-96c6-ca94a10ee22f
md"""
!!! hint
	- Use the 'important' quantities mentioned above.
	- Set up:
	``` math
	\begin{align}
	c_p\,\rho\,V_1\,\cfrac{dT_1}{dt} &= \dots \\
	c_p\,\rho\,V_2\,\cfrac{dT_2}{dt} &= \dots
	\end{align}
	```
	- Think well before putting a $+$ or $-$ sign before the terms while setting-up the differential equations!
"""

# ╔═╡ 41ed50f3-9270-4563-b31f-30351f07468b
solution(md"""
``` math
\left\{
\begin{array}{l}
c_p\,\rho\,V_1\,\cfrac{dT_1}{dt} &= c_p\,\rho\,Q \left( T_{in} - T_1 \right)  -  \lambda\,A_1 \left(T_1 - T_a\right) \\
c_p\,\rho\,V_2\,\cfrac{dT_2}{dt} &= c_p\,\rho\,Q \left( T_1 - T_2 \right) - \lambda\,A_2 \left(T_2 - T_a\right) + P_h
\end{array}
\right.
```
""")

# ╔═╡ bdabfb0c-66ff-4cb7-a61f-f31fb49eeb34
md"""
Temperature sensors measures the fluid temperatures $T_1$ and $T_2$ in both reactors. We are interested in the evolution of $T_1$ and $T_2$ as a function of time.
"""

# ╔═╡ f4431fbe-2823-4ce5-ac35-16f51879d393
md"""
## Setting up the equations
"""

# ╔═╡ ceb468ea-c05d-4938-bc95-31256af53a24
md"""
The values of the parameters are listed here

| Parameter | Value | Unit | Meaning |
|:---------- |:---------- |:------------|:------------|
| `Q`    | 0.05e-3 | $m^3/s$  |    flow rate   |
| `V₁`    | 0.004 | $m^3$  |  volume of liquid in reactor 1  |
| `V₂`    | 0.006 | $m^3$  |  volume of liquid in reactor 2  |
| `Tin`    | 40.0 | $^{\circ}C$  |   temperature of incoming liquid |
| `λ`    | 240 | $W/(m^2\,^{\circ}C)$  |   heat transfer coefficient |
| `A₁`    | 0.1 | $m^2$  |  wall area of reactor 1  |
| `A₂`    | 0.2 | $m^2$  |  wall area of reactor 2  |
| `cₚ`    | 4200 | $J/(kg\,^{\circ}C$  |  thermal heat capacity of fluid  |
| `ρ`   | 1000 | $kg/m^3$  |  density of fluid  |
| `Tₐ`    | 20.0 | $^{\circ}C$  |   ambient temperature |
| `Pₕ`    | 5000 | $W$  |   power of heater (in reactor 2) |

The parameters $T_{in}$ and $P_{h}$ will be called *manipulable parameters* since they can be (suddenly) changed by a control engineer.

"""

# ╔═╡ c76afb2e-a12a-4540-a622-c00304cce2e4
md"""
Define the variables for the temperature in both reactors. Use `T₁` and `T₂`.
"""

# ╔═╡ b60730ff-a3de-4bd0-b074-6f4ee0b35957
# @variables 

# ╔═╡ 6b189311-3fa1-4a81-b65a-74e2d972b047
md"""
Define the parameters for this model and assign their corresponding values.
"""

# ╔═╡ 07e57de7-6ac1-4351-b2c2-fd60212c15b0
# @parameters 

# ╔═╡ 01a57415-a5f1-47be-b164-6b951d7ae6e7
md"""
Set up the equation that describes the change in $T_1$.
"""

# ╔═╡ ef9b10d8-c1df-41aa-a7c6-1b1fbc0b727c
eq_T1 = missing

# ╔═╡ bb57ddca-de84-46cf-9751-88b9501a30d5
md"""
Set up the equation that describes the change in $T_2$.
"""

# ╔═╡ 420b93da-59bc-4aa0-af4c-3e9df7d00107
eq_T2 = missing

# ╔═╡ 12f793e3-6c33-4394-8b24-2aa046320943
md"""
Bundle both equations.
"""

# ╔═╡ d2160d11-f0d3-40cb-a921-7668af9d5b31
eqns_reactors = missing

# ╔═╡ 18b41ac9-0795-46e5-a9bb-9abb7695c601
md"""
## Part 1: simple simulation
"""

# ╔═╡ 30bd37c4-16b1-48c8-ad51-f3fdcc34bda3
md"""
In this part we will assume that the manipulable parameters remain constant with values as listed before.
"""

# ╔═╡ 4995f295-ca60-4035-8cb3-f3607d4ad3c4
md"""
### Building the ODE system
"""

# ╔═╡ 8748f2eb-2a02-426e-b71b-166298a2667f
md"""
Build the model. Name it `sys1_reactors`.
"""

# ╔═╡ de8c6064-8809-4aee-8809-b7cc0b07122f
# @mtkbuild missing

# ╔═╡ 9bfc0267-074b-43cf-b199-52cc0e4a0f46
md"""
### Create and solve the ODE problem
"""

# ╔═╡ 8f9d348e-d1ea-484d-a5e8-6174a7ce858b
md"""
Create the ODE problem. Assume an initial temperature of 15.0 $^{\circ}C$ for $T_1$ and $T_2$ and simulate during 600.0 seconds. Use `[]` for the parameters argument since their values have been set before.
"""

# ╔═╡ 85c9b5e6-3ca5-4e93-bdce-4cdf690d44b8
oprob1_reactors = missing

# ╔═╡ b10799d3-3da8-4f8c-af02-9e57448cbcf4
md"""
Solve the ODE problem. Use `Tsit5()` and `saveat=1`.
"""

# ╔═╡ a5c54103-1fc9-4358-9a58-eb6bd88e7d18
sol1_reactors = missing

# ╔═╡ 09b21e52-1648-42e8-90a9-a0719a6440df
md"""
### Plotting results
"""

# ╔═╡ 01c1fced-1852-422b-ab14-32751f903495
md"""
Make a plot of $T_1$ and $T_2$ over time.
"""

# ╔═╡ 14701ae1-99c6-47f9-ab59-e880ae3ea0d4
missing

# ╔═╡ 9d1b4080-300f-496d-aa40-af8bc6c080b1
md"""
Retrieve the final values of $T_1$ and $T_2$ in the solution (i.e., at $t=600\;s$).
"""

# ╔═╡ b1077d9d-8523-4ad4-83bf-bcee0dc11b26
missing

# ╔═╡ 6d05657b-24f2-49a3-b568-bd25e52394e6
missing

# ╔═╡ 36d951d4-d354-450d-8dbe-1b3a088fcef0
md"""
### Calculating steady state values
"""

# ╔═╡ f0e79793-5862-4806-9865-f18323d4577a
md"""
Calculate the steady state values using `SteadyStateProblem` and `solve`. Use the above retrieved final values as a first guess.
"""

# ╔═╡ a2d2ae7a-9c55-4a3a-9428-8376fdf5ed8b
equil_T_vals = missing

# ╔═╡ 9c7a9403-c0f9-4500-b84c-1478d72269eb
md"""
Display the steady state values.
"""

# ╔═╡ 25876e54-adfe-4f2b-b844-7e19b3820868
missing

# ╔═╡ 572aa97b-ea95-453b-8c27-d7f579ee8823
missing

# ╔═╡ a2725289-11ff-408e-8375-17f2b54f3a4f
md"""
!!! questions
	1. Why is the steady state value of $T_1$ not equal to $T_{in}$?
	2. Why is the steady state value of $T_2$ much larger than $T_1$?
"""

# ╔═╡ 5af588fa-f2e5-4b25-b431-cb93e0683ba7
md"""
Answers:
1. missing
2. missing
"""

# ╔═╡ a76bb35a-f405-4fc8-b8fb-58d513f86b2e
md"""
## Part 2: sudden change in $T_{in}$ and $P_h$
"""

# ╔═╡ e71032db-1e3e-48a5-b660-2b15dcbf8480
md"""
In this part we will assume that both manipulable parameters $T_{in}$ and $P_h$ suddenly change at some point in time:
- `Tin` changes from 40.0 $^{\circ}C$ to 45.0 $^{\circ}C$, and
- `Pₕ` changes from 5000 $W$ to 2250 $W$,
at $t =$ 600 seconds.
"""

# ╔═╡ 15a37c5f-921c-4543-b605-cf89fdb898a2
md"""
### Building the ODE system
"""

# ╔═╡ 2dc76e23-723e-4c1d-a703-0cb7fd0edfa0
md"""
Build a new model where you include both discrete events. Name it `sys2_reactors`.
"""

# ╔═╡ cb30b1f0-6fbb-4c11-9cfa-ee8ed3ebc74a
# @mtkbuild missing

# ╔═╡ 0b13a875-6c68-4ca9-812e-fc84989d70d4
md"""
### Create and solve the ODE problem
"""

# ╔═╡ b9aa284a-8058-4b76-ac40-7aefbb43f377
md"""
Create the ODE problem. Take the same initial values for the temperatures as before but now simulate for 1200 seconds.
"""

# ╔═╡ 0ea3b443-96d3-4291-9cac-266735b23454
oprob2_reactors = missing

# ╔═╡ 5c59de0c-31d2-41b1-b89a-c77cdca2c322
md"""
Solve the ODE problem. Make a deepcopy of the ODE problem, use `Tsit5()` and `saveat=1`.
"""

# ╔═╡ 1eaa2ae5-4d1f-4c24-9358-3cd150ab7ccd
sol2_reactors = missing

# ╔═╡ 6820dc35-d1e4-45e2-b5b0-664a62e172b8
md"""
### Plotting results
"""

# ╔═╡ 19acc5b0-484a-4a36-8a3b-a59421848848
md"""
Make a plot of $T_1$ and $T_2$ over time.
"""

# ╔═╡ b9774204-0cb5-47c0-9a5c-f3508ff86435
missing

# ╔═╡ 21b22b69-1aeb-49b7-ab7f-075f978e1acb
md"""
!!! question
	1. Is the variable $T_1$ affected by (A) only $T_{in}$, (B) only $P_h$, (C) both $T_{in}$ and $P_h$, or (D) none of them? Try to reason using the model equation for the rate of change of $T_1$.
	2. Is the variable $T_2$ affected by (A) only $T_{in}$, (B) only $P_h$, (C) both $T_{in}$ and $P_h$, or (D) none of them? Try to reason using the model equation for the rate of change of $T_2$.
"""

# ╔═╡ 18f1bdd7-da96-49d7-96d7-fe456818cd4b
md"""
Answers:
1. missing
2. missing
"""

# ╔═╡ Cell order:
# ╟─4bf49e49-7675-46eb-acfd-603dd68ac027
# ╠═7d61614c-9e15-11f0-1ce8-af091a3beb9a
# ╠═2e418e46-84fd-4922-9985-a92a412c1fb1
# ╠═d0a12f36-2592-4c3f-9581-6e41fbe5117f
# ╠═225d31d2-a974-4494-8813-c93484ce721e
# ╟─0e4065c2-281c-439a-a69b-e4eca3b5ae41
# ╟─7e2a3f8c-5f1c-4357-ae92-f4c1c677fde3
# ╟─9e8c3a6b-6d69-4f9b-a8ea-de732de886ca
# ╟─6a5cc6f2-b6b8-40fe-9b52-1a3e2d5d4fb7
# ╟─d1480cba-5e0f-4767-8164-b8ed5291e888
# ╟─826b1cda-3fcc-46e6-b186-d031929a569b
# ╟─409e4f75-0cc7-4f2b-96c6-ca94a10ee22f
# ╟─41ed50f3-9270-4563-b31f-30351f07468b
# ╟─bdabfb0c-66ff-4cb7-a61f-f31fb49eeb34
# ╟─f4431fbe-2823-4ce5-ac35-16f51879d393
# ╟─ceb468ea-c05d-4938-bc95-31256af53a24
# ╟─c76afb2e-a12a-4540-a622-c00304cce2e4
# ╠═b60730ff-a3de-4bd0-b074-6f4ee0b35957
# ╟─6b189311-3fa1-4a81-b65a-74e2d972b047
# ╠═07e57de7-6ac1-4351-b2c2-fd60212c15b0
# ╟─01a57415-a5f1-47be-b164-6b951d7ae6e7
# ╠═ef9b10d8-c1df-41aa-a7c6-1b1fbc0b727c
# ╟─bb57ddca-de84-46cf-9751-88b9501a30d5
# ╠═420b93da-59bc-4aa0-af4c-3e9df7d00107
# ╟─12f793e3-6c33-4394-8b24-2aa046320943
# ╠═d2160d11-f0d3-40cb-a921-7668af9d5b31
# ╟─18b41ac9-0795-46e5-a9bb-9abb7695c601
# ╟─30bd37c4-16b1-48c8-ad51-f3fdcc34bda3
# ╟─4995f295-ca60-4035-8cb3-f3607d4ad3c4
# ╟─8748f2eb-2a02-426e-b71b-166298a2667f
# ╠═de8c6064-8809-4aee-8809-b7cc0b07122f
# ╟─9bfc0267-074b-43cf-b199-52cc0e4a0f46
# ╟─8f9d348e-d1ea-484d-a5e8-6174a7ce858b
# ╠═85c9b5e6-3ca5-4e93-bdce-4cdf690d44b8
# ╟─b10799d3-3da8-4f8c-af02-9e57448cbcf4
# ╠═a5c54103-1fc9-4358-9a58-eb6bd88e7d18
# ╟─09b21e52-1648-42e8-90a9-a0719a6440df
# ╟─01c1fced-1852-422b-ab14-32751f903495
# ╠═14701ae1-99c6-47f9-ab59-e880ae3ea0d4
# ╟─9d1b4080-300f-496d-aa40-af8bc6c080b1
# ╠═b1077d9d-8523-4ad4-83bf-bcee0dc11b26
# ╠═6d05657b-24f2-49a3-b568-bd25e52394e6
# ╟─36d951d4-d354-450d-8dbe-1b3a088fcef0
# ╟─f0e79793-5862-4806-9865-f18323d4577a
# ╠═a2d2ae7a-9c55-4a3a-9428-8376fdf5ed8b
# ╟─9c7a9403-c0f9-4500-b84c-1478d72269eb
# ╠═25876e54-adfe-4f2b-b844-7e19b3820868
# ╠═572aa97b-ea95-453b-8c27-d7f579ee8823
# ╟─a2725289-11ff-408e-8375-17f2b54f3a4f
# ╠═5af588fa-f2e5-4b25-b431-cb93e0683ba7
# ╟─a76bb35a-f405-4fc8-b8fb-58d513f86b2e
# ╟─e71032db-1e3e-48a5-b660-2b15dcbf8480
# ╟─15a37c5f-921c-4543-b605-cf89fdb898a2
# ╟─2dc76e23-723e-4c1d-a703-0cb7fd0edfa0
# ╠═cb30b1f0-6fbb-4c11-9cfa-ee8ed3ebc74a
# ╟─0b13a875-6c68-4ca9-812e-fc84989d70d4
# ╟─b9aa284a-8058-4b76-ac40-7aefbb43f377
# ╠═0ea3b443-96d3-4291-9cac-266735b23454
# ╟─5c59de0c-31d2-41b1-b89a-c77cdca2c322
# ╠═1eaa2ae5-4d1f-4c24-9358-3cd150ab7ccd
# ╟─6820dc35-d1e4-45e2-b5b0-664a62e172b8
# ╟─19acc5b0-484a-4a36-8a3b-a59421848848
# ╠═b9774204-0cb5-47c0-9a5c-f3508ff86435
# ╟─21b22b69-1aeb-49b7-ab7f-075f978e1acb
# ╠═18f1bdd7-da96-49d7-96d7-fe456818cd4b
