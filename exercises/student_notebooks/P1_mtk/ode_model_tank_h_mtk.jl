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

# ╔═╡ e63a0a5a-954d-11f0-1ed0-a3fcaefba117
using Pkg; Pkg.activate("..")

# ╔═╡ 1fe30acc-26dd-45e7-9264-1571c47d07c1
using StatsPlots, PlutoUI; TableOfContents()

# ╔═╡ 7d45ffff-e093-41c4-a0e3-50ebdba5f816
using OrdinaryDiffEq

# ╔═╡ d6f533fd-e593-4fec-8c87-53a157eeed48
using ModelingToolkit

# ╔═╡ 9c093210-cc9f-40a1-8e35-48f58a7fa491
using ModelingToolkit: t_nounits as t, D_nounits as D

# ╔═╡ 23d4f233-1305-47a3-8b8f-89873765d7bc
md"""
# Exercise: Cylindrical tank (water level)
"""

# ╔═╡ d99457cc-8223-4970-b52d-45d42bf7f030
solution(text) = Markdown.MD(Markdown.Admonition("hint", "Solution", [text]));

# ╔═╡ 15576faf-0986-4bcd-8e13-24afc9eb70d4
md"""
![](https://users.ugent.be/~gvhaelew/fig/tank_height.png)
"""

# ╔═╡ eee09dae-82eb-4e3b-98b3-daad57c42d49
md"""
## Deriving the model equations
"""

# ╔═╡ e638ef6e-00c3-4666-86d8-5803c21abe64
md"""
A cylindrical tank is filled with water (its density is denoted $\rho$) with a constant flow rate of $Q_{in}$. The outlet flow rate $Q_{out}$ depends on the square root of the water level (height $h$) in the tank in the following manner:
``` math
Q_{out} = \cfrac{\sqrt{h}}{R}
```
where $R$ is the resistence coefficient of the orifice. The cross section $A$ of the tank is constant, hence, the mass $M$ of water in the tank only varies with the height of the water:
``` math
M = \rho\,V = \rho\,A\,h
```
"""

# ╔═╡ cab795f1-f477-4d02-8428-3d26778b49ae
md"""
!!! task
	Derive as an exercise a single model equation for the height $h$ of the water inside the tank.
"""

# ╔═╡ 5ab3e018-e59c-4b80-ab19-3d2368a21bb5
md"""
!!! hint
	- Start by setting up the rate of change of the mass of water in the tank.
	``` math
	\cfrac{dM}{dt} = \cdots
	```
	work toward:
	``` math
	\cfrac{dh}{dt} = \cdots
	```

"""

# ╔═╡ d850e7ac-dfdf-4ccb-aeaf-01cd5ae0d8d9
solution(md"""
``` math
\cfrac{dh}{dt} = \cfrac{Q_{in}}{A} - \cfrac{\sqrt{h}}{AR}
```
""")

# ╔═╡ b37112ed-2436-42ba-a576-5bca422a0e27
md"""
Apart from the height of the water, we are also interested in the hydrostatic pressure $p$ (in $bar$, note that $1\;Pa = 10^{-5}\;bar$) at the bottom of the tank:
``` math
p = \rho\,g\,h
```
where $g$ is the gravitational constant.
"""

# ╔═╡ 3ae63528-8a78-4884-85f6-1c6dc3ae178e
md"""
Define variables for the height $h$ of the water in the tank and the hydrostatic pressure $p$. Use the following variable names: `h` and `p`. Mention the dependency on the time $t$.
"""

# ╔═╡ fac7ad9e-a4a2-4f7c-90f1-0111f9b3639c
# @variables missing

# ╔═╡ 1df30e1f-190e-4e26-9133-152afc60119d
md"""
Define the parameters for this model and assign their corresponding values.

| Parameter | Value | Unit | Meaning |
|:---------- |:---------- |:------------|:------------|
| $A$    | 0.152 | $m^2$  |    cross sectional area   |
| $\rho$    | 1000.0 | $kg/m^3$  |  water density  |
| $g$    | 9.81 | $m/s^2$  |   gravitational constant |
| $Q_{in}$    |  | $m^3/s$  |   inlet flow |
| $R$    |  | $s/m^{5/2}$  |   res. coeff. orifice |

The values of $Q_{in}$ and $R$ will be set later when creating the ODE problem using values from sliders.
"""

# ╔═╡ b12868c1-85ef-4a51-b146-ea7be5f119f6
# @parameters missing

# ╔═╡ b4de5802-7cd0-496f-be3d-816511ff0281
md"""
## Part 1: filling up
"""

# ╔═╡ 29b1088f-4675-4538-bb78-26cfadc6ce28
md"""
### Setting up the equations
"""

# ╔═╡ b22c45f3-bfa8-4715-97bc-4b34affe99be
md"""
In this part we will analyze the height of the water $h$ over time for different values of the initial height of the water $h_0$, the inlet flow $Q_{in}$ and the resistance coefficient $R$ of the orifice. The latter three will be set using values from sliders (see below, just above the plot of the variables).
"""

# ╔═╡ e5b1b418-8fba-4c58-aa8d-8a3b0e1b4e05
md"""
Set up the equation for the rate of change in height. Multiply the term $-\cfrac{\sqrt{h}}{R}$ with $(h > 0)$. *Think about why this must be done.*
"""

# ╔═╡ 4c14b63b-8ba3-4df3-b41b-8c5e02be8b45
change_height = missing

# ╔═╡ be764668-1a10-4e79-9cc0-585eeb25caf0
md"""
Set up the expression that will keep track of the hydrostatic pressure. Don't forget to include `*1e-5` in the term in order to have it in $bar$.
"""

# ╔═╡ 4d24b00a-ec24-4fe9-a31d-472904f2d64e
eq_pressure = missing

# ╔═╡ f3acc97b-e501-4797-9334-fd3fb7197cbe
md"""
Bundle the equation and the expression.
"""

# ╔═╡ 25973753-4374-41ae-a391-87c750289653
eqns_tank = missing

# ╔═╡ a54ebeb7-cd9d-4eac-ab2d-f516758fc045
md"""
### Building the ODE system
"""

# ╔═╡ a9ef1b26-e3be-45e9-8604-415dfb885a77
md"""
Build the model. Name it `sys1_tank`.
"""

# ╔═╡ 81069e6e-0785-41c2-bbcd-4b9bc25930e0
# @mtkbuild missing

# ╔═╡ 29e4b87b-00ee-4dfe-a0a8-2dbd18554f11
md"""
### Create and solve the ODE problem
"""

# ╔═╡ 98cf7e7c-c52a-49d8-8d00-5d3d85582c50
md"""
Create the ODE problem. Use the variables `h₀`, `QinLmin` and `R_val` to set the initial height, the inlet flow and the resistance coefficient, respectively.

**Important remarks:**
- The value of `QinLmin` is in liters per minute ($L/min$). Hence, you need to multiply it with `1e-3/60` to have it in $m^3/s$ before assigning it to `Qin`.
- You don't need to reassign the parameters `A`, `ρ` and `g` because they were set before when defining the actual parameters.
"""

# ╔═╡ 6782c5da-c3ce-4101-bf42-ed0849024162
oprob1_tank = missing

# ╔═╡ 401cb8fe-4b56-40ef-a763-741867c55685
md"""
Solve the ODE problem. Use `Tsit5()` and `saveat=1`.
"""

# ╔═╡ 44878548-9661-421e-a9ab-b82ca92a6310
sol1_tank = missing

# ╔═╡ e9ae80c7-7554-42dd-8cda-8a0bcb1f0a5a
md"""
### Plotting results
"""

# ╔═╡ 5ee777ac-71a2-449d-a678-ffe0da09287c
# @bind h₀ Slider(0:0.1:2.5, default=0, show_value=true)

# ╔═╡ a98352c7-9580-4c94-b250-1a1f24d78771
# @bind QinLmin Slider(0:10:120, default=100, show_value=true)

# ╔═╡ 533b1c79-3408-48d0-b3ee-04ec05e89b65
# @bind R_val Slider(400:10:800, default=600, show_value=true)

# ╔═╡ 4358d80a-8240-41b5-acfc-d4c0b1b22876
md"""
Plot the height and the pressure. Use the options `ylim=(0.0, 2.6)` and `idxs=[h, p]`. Play with the sliders above to see their effect.
"""

# ╔═╡ 8b15f793-3628-486a-8fd3-c590db3bda7e
missing

# ╔═╡ 8da6d0e4-01aa-428a-b1e8-e951087115f6
md"""
!!! questions
	1. On what parameters does the equilibrium value for the height depend?
	2. Can you calculate this equilibrium value analytically using the equation for the rate of change in height? How does the expression look like?
	3. To what values do you need to set `h₀` and `Qin` in order to completely drain the tank?
"""

# ╔═╡ 474e479c-757e-498f-b40c-e7e8a3839770
md"""
Answers:
1. missing
2. missing
3. missing
"""

# ╔═╡ c634e6c9-78c3-4501-91ea-b5905ef8f26b
md"""
### Calculating the steady state value
"""

# ╔═╡ 0e2f2d97-3091-41e2-ae02-2e8307fb6795
md"""
Calculate the steady state value using `SteadyStateProblem`. Use `1.0` as a first guess for `h`. Make this calculation only when `QinLmin` is greater than $30\;L/min$.
"""

# ╔═╡ 0308b949-5376-4d4e-9c2e-c4c0e5120991
equil_val = missing

# ╔═╡ 9a3c9576-a3cc-46a3-911f-906e7a019e9b
md"""
Display the steady state value.
"""

# ╔═╡ 91066873-7530-44d5-8fac-e508066e9b1a
missing

# ╔═╡ c3263780-eaf9-43fb-ba7d-1496cf59b0fe
md"""
## Part 2: filling up and draining
"""

# ╔═╡ bca47fff-a618-4fb2-a059-2b98fd1632a0
md"""
In this part the tank will be filled up to a certain height and then completely drained. In order to achieve this, we will use a continuous event.
"""

# ╔═╡ 5559a586-59f2-470a-9bf8-efba2ce42270
md"""
The tank is initially empty and filled up with an inlet flow `Qin` of `100` $L/min$. Take `600` for the value of `R`. The continuous event should state that when `h` becomes `0.99`, the value of `Qin` should be set to `0`.
"""

# ╔═╡ 5bd86c38-0338-4bd1-bdc7-a7c494b1f17a
md"""
### Building the ODE system
"""

# ╔═╡ 618590fc-7fb9-4020-858b-d9b2d43fbd70
md"""
Build the new model that includes the continuous event. Name it `sys2_tank`.
"""

# ╔═╡ 0a7dc8af-e7d2-40fd-8c60-d266e9d6c0b7
# @mtkbuild missing

# ╔═╡ 919b6c52-eba1-4a20-a1be-6cd1539a8164
md"""
### Create and solve the ODE problem
"""

# ╔═╡ 84f7c85f-08f6-4ba0-bf80-fddf86506143
md"""
Create the new ODE problem. Set the correct values for the initial height and the parameters. The simulation time is the same as in Part 1.
"""

# ╔═╡ 3c9acd41-8451-4a13-8dff-f74a41f0bd3c
oprob2_tank = missing

# ╔═╡ b9c8d9aa-2f9f-487c-87b4-048d4fbe8aa9
md"""
Solve the ODE problem. Use `Tsit5()` and `saveat=1`. Don't forget to take a `deepcopy` of the ODE problem.
"""

# ╔═╡ c53ecb42-cd0a-43a2-b244-22eb6833ea92
sol2_tank = missing

# ╔═╡ b52e47a1-d642-4529-97e6-3d3837be657f
md"""
### Plotting results
"""

# ╔═╡ 333cafcc-c6bb-433c-9f34-fdcf27d59508
md"""
Plot the height and the pressure. Use the option `idxs=[h, p]`.
"""

# ╔═╡ 0de3c62d-b02e-4f9d-8213-58c4199e2662
missing

# ╔═╡ b5196b36-9092-4c37-86a8-ebe5da4e14c1
md"""
!!! questions
	1. What is the time needed to drain the tank after the height became `0.99`? Derive this (numerically) using the results. 
	2. Can you calculate this time also analytically?
	3. Do both numerical and analytical draining times more or less correspond?

	Hints:
	- You can you `findfirst(==(x), V)` to find the first index of the value `x` in the vector `V`.
	- To slice a vector `V` from a certain index `i` to the `end`, do `V[i:end]`.
	- You can use `findfirst(<(0), V)` to find the first index of the value less than `0` in the vector `V`. Don't use `findfirst(==(0), V)` because the numerical values will never be exactly `0` in the tail of the height vector.

"""

# ╔═╡ 5d78c408-b10a-4685-8b16-1ee30315aae8
md"""
Numerical derivation:
"""

# ╔═╡ ea6725ef-f392-4161-96fa-0da3fcffbb94
i_max = missing

# ╔═╡ 1c21b251-7ef4-4f07-a7f1-cf434fdaaac3
Δi_zero = missing

# ╔═╡ 704c37e8-ee9e-4ae9-82a1-aa95fefb14d1
missing

# ╔═╡ 6a807ca4-6b20-49d4-b095-c249afca6772
md"""
Analytical calculation (optional):
"""

# ╔═╡ 66fe0cf1-fcdf-4d4d-b7f1-bbc1b7e558f0
missing

# ╔═╡ b4da3154-63ef-486b-8d2a-874f93624884
md"""
Answers:
1. missing
2. missing
3. missing
"""

# ╔═╡ d4890a27-dca4-4207-a827-ae1968fa9fc5
md"""
## Part 3: modifying Qin and R
"""

# ╔═╡ 39099e15-ed02-47bb-8801-8b59ced0b0db
md"""
In this part the values of $Q_{in}$ and $R$ will be modified at distinct moments in time. In order to achieve this, we will use discrete events. The tank is initially empty and filled up with an inlet flow `Qin` of `100` $L/min$. Take `600` for the value of `R`.

At the time instant `600.0` $s$ the value of `Qin` will increase by $10\;\%$ (i.e., multiplied by `1.1`). At the time instant `1200.0` $s$ the value of `R` will be decreased by $25\;\%$ (i.e., multiplied by `0.75`).
"""

# ╔═╡ 2b68cd82-3c81-4a4f-a27d-e5cb35a7b0f0
md"""
### Building the ODE system
"""

# ╔═╡ c1f4c7ab-88b2-4282-8a6f-47f9161951eb
md"""
Build the new model that includes both discrete events. Name it `sys3_tank`.
"""

# ╔═╡ 07cc0c6d-1514-4d49-bab3-12794e199729
# @mtkbuild missing

# ╔═╡ e395c4c4-3761-4b74-9b1d-4e3d32366333
md"""
### Create and solve the ODE problem
"""

# ╔═╡ b917ecfb-bda3-453b-a66e-25d83183d5ee
md"""
Create the new ODE problem. Set the correct values for the initial height and the parameters. The simulation time is the same as in Part 1 and 2.
"""

# ╔═╡ ffe3c2f4-dee4-46bb-9f94-d2a59b621ffc
oprob3_tank = missing

# ╔═╡ 0a4cebcd-88a9-4193-bce4-60c619e592d4
md"""
Solve the ODE problem. Use `Tsit5()` and `saveat=1`. Don't forget to take a `deepcopy` of the ODE problem.
"""

# ╔═╡ 27f8d05b-94d1-4855-b7ce-33d7995e10d3
sol3_tank = missing

# ╔═╡ 8abc6dd9-ad9e-4cc3-a4c8-cf6faa6b3b39
md"""
### Plotting results
"""

# ╔═╡ 45479aac-1f8d-4334-a6e5-38bad49131a1
missing

# ╔═╡ 9089bc00-7c8e-48a3-b5ed-ec1aa9ec11df
md"""
!!! question
	Is the evolution of the height accoding to your intuition? Think about what should happen to the height when you increase $Q_{in}$ and decrease $R$.
"""

# ╔═╡ 3b8c0bf9-9b2a-4fe5-8cf8-73653d193ca6
md"""
Answer: missing
"""

# ╔═╡ Cell order:
# ╟─23d4f233-1305-47a3-8b8f-89873765d7bc
# ╠═e63a0a5a-954d-11f0-1ed0-a3fcaefba117
# ╠═1fe30acc-26dd-45e7-9264-1571c47d07c1
# ╠═7d45ffff-e093-41c4-a0e3-50ebdba5f816
# ╠═d6f533fd-e593-4fec-8c87-53a157eeed48
# ╠═9c093210-cc9f-40a1-8e35-48f58a7fa491
# ╟─d99457cc-8223-4970-b52d-45d42bf7f030
# ╟─15576faf-0986-4bcd-8e13-24afc9eb70d4
# ╟─eee09dae-82eb-4e3b-98b3-daad57c42d49
# ╟─e638ef6e-00c3-4666-86d8-5803c21abe64
# ╟─cab795f1-f477-4d02-8428-3d26778b49ae
# ╟─5ab3e018-e59c-4b80-ab19-3d2368a21bb5
# ╟─d850e7ac-dfdf-4ccb-aeaf-01cd5ae0d8d9
# ╟─b37112ed-2436-42ba-a576-5bca422a0e27
# ╟─3ae63528-8a78-4884-85f6-1c6dc3ae178e
# ╠═fac7ad9e-a4a2-4f7c-90f1-0111f9b3639c
# ╟─1df30e1f-190e-4e26-9133-152afc60119d
# ╠═b12868c1-85ef-4a51-b146-ea7be5f119f6
# ╟─b4de5802-7cd0-496f-be3d-816511ff0281
# ╟─29b1088f-4675-4538-bb78-26cfadc6ce28
# ╟─b22c45f3-bfa8-4715-97bc-4b34affe99be
# ╟─e5b1b418-8fba-4c58-aa8d-8a3b0e1b4e05
# ╠═4c14b63b-8ba3-4df3-b41b-8c5e02be8b45
# ╟─be764668-1a10-4e79-9cc0-585eeb25caf0
# ╠═4d24b00a-ec24-4fe9-a31d-472904f2d64e
# ╟─f3acc97b-e501-4797-9334-fd3fb7197cbe
# ╠═25973753-4374-41ae-a391-87c750289653
# ╟─a54ebeb7-cd9d-4eac-ab2d-f516758fc045
# ╟─a9ef1b26-e3be-45e9-8604-415dfb885a77
# ╠═81069e6e-0785-41c2-bbcd-4b9bc25930e0
# ╟─29e4b87b-00ee-4dfe-a0a8-2dbd18554f11
# ╟─98cf7e7c-c52a-49d8-8d00-5d3d85582c50
# ╠═6782c5da-c3ce-4101-bf42-ed0849024162
# ╟─401cb8fe-4b56-40ef-a763-741867c55685
# ╠═44878548-9661-421e-a9ab-b82ca92a6310
# ╟─e9ae80c7-7554-42dd-8cda-8a0bcb1f0a5a
# ╠═5ee777ac-71a2-449d-a678-ffe0da09287c
# ╠═a98352c7-9580-4c94-b250-1a1f24d78771
# ╠═533b1c79-3408-48d0-b3ee-04ec05e89b65
# ╟─4358d80a-8240-41b5-acfc-d4c0b1b22876
# ╠═8b15f793-3628-486a-8fd3-c590db3bda7e
# ╟─8da6d0e4-01aa-428a-b1e8-e951087115f6
# ╠═474e479c-757e-498f-b40c-e7e8a3839770
# ╟─c634e6c9-78c3-4501-91ea-b5905ef8f26b
# ╟─0e2f2d97-3091-41e2-ae02-2e8307fb6795
# ╠═0308b949-5376-4d4e-9c2e-c4c0e5120991
# ╟─9a3c9576-a3cc-46a3-911f-906e7a019e9b
# ╠═91066873-7530-44d5-8fac-e508066e9b1a
# ╟─c3263780-eaf9-43fb-ba7d-1496cf59b0fe
# ╟─bca47fff-a618-4fb2-a059-2b98fd1632a0
# ╟─5559a586-59f2-470a-9bf8-efba2ce42270
# ╟─5bd86c38-0338-4bd1-bdc7-a7c494b1f17a
# ╟─618590fc-7fb9-4020-858b-d9b2d43fbd70
# ╠═0a7dc8af-e7d2-40fd-8c60-d266e9d6c0b7
# ╟─919b6c52-eba1-4a20-a1be-6cd1539a8164
# ╟─84f7c85f-08f6-4ba0-bf80-fddf86506143
# ╠═3c9acd41-8451-4a13-8dff-f74a41f0bd3c
# ╟─b9c8d9aa-2f9f-487c-87b4-048d4fbe8aa9
# ╠═c53ecb42-cd0a-43a2-b244-22eb6833ea92
# ╟─b52e47a1-d642-4529-97e6-3d3837be657f
# ╟─333cafcc-c6bb-433c-9f34-fdcf27d59508
# ╠═0de3c62d-b02e-4f9d-8213-58c4199e2662
# ╟─b5196b36-9092-4c37-86a8-ebe5da4e14c1
# ╟─5d78c408-b10a-4685-8b16-1ee30315aae8
# ╠═ea6725ef-f392-4161-96fa-0da3fcffbb94
# ╠═1c21b251-7ef4-4f07-a7f1-cf434fdaaac3
# ╠═704c37e8-ee9e-4ae9-82a1-aa95fefb14d1
# ╟─6a807ca4-6b20-49d4-b095-c249afca6772
# ╠═66fe0cf1-fcdf-4d4d-b7f1-bbc1b7e558f0
# ╠═b4da3154-63ef-486b-8d2a-874f93624884
# ╟─d4890a27-dca4-4207-a827-ae1968fa9fc5
# ╟─39099e15-ed02-47bb-8801-8b59ced0b0db
# ╟─2b68cd82-3c81-4a4f-a27d-e5cb35a7b0f0
# ╟─c1f4c7ab-88b2-4282-8a6f-47f9161951eb
# ╠═07cc0c6d-1514-4d49-bab3-12794e199729
# ╟─e395c4c4-3761-4b74-9b1d-4e3d32366333
# ╟─b917ecfb-bda3-453b-a66e-25d83183d5ee
# ╠═ffe3c2f4-dee4-46bb-9f94-d2a59b621ffc
# ╟─0a4cebcd-88a9-4193-bce4-60c619e592d4
# ╠═27f8d05b-94d1-4855-b7ce-33d7995e10d3
# ╟─8abc6dd9-ad9e-4cc3-a4c8-cf6faa6b3b39
# ╠═45479aac-1f8d-4334-a6e5-38bad49131a1
# ╟─9089bc00-7c8e-48a3-b5ed-ec1aa9ec11df
# ╠═3b8c0bf9-9b2a-4fe5-8cf8-73653d193ca6
