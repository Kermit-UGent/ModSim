### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ a6fbee1a-9863-11f0-1dfe-cf500206e163
using Pkg; Pkg.activate("..")

# ╔═╡ 21c2d9ad-9411-4e05-be3e-7f37deebc70d
using Plots, PlutoUI; TableOfContents()

# ╔═╡ eaa957df-d305-48fb-b2c2-4b4eb29a08be
using DifferentialEquations, ModelingToolkit

# ╔═╡ 5f3f3184-7bf8-4f40-80ca-bb48ad85e2f9
using ModelingToolkit: t_nounits as t, D_nounits as D

# ╔═╡ 0d12503e-5480-4ed8-9b4a-569b7801f855
md"""
# Exercise: Cylindrical tank with cooling jacket
"""

# ╔═╡ fdf5f147-5e6d-4730-8f16-ddc62ff9fcfe
solution(text) = Markdown.MD(Markdown.Admonition("hint", "Solution", [text]));

# ╔═╡ 82207f70-ca9a-409b-952d-520f3b8babc3
md"""
## Deriving the model equations
"""

# ╔═╡ defe58ce-d407-4ab4-a4e4-e37fc061149c
md"""
A cylindrical tank, surrounded by a cooling jacket, contains an initial amount of water (initial temperature $T_c$ and initial height $h_0$). The cooling jacket remains at a fixed temperature $T_c$. At a certain moment ($t=0$), a stream of warmer water with flow rate $Q_{in}$ and temperature $T_{in}$ enters the tank, and at the same moment a valve at the bottom of the tank is opened. The water in the tank is well mixed at all times. We are interested in the evolution of the temperature $T$ of the water in the tank and its height $h$.

We denote $\rho$ as the mass density of the water, $c_p$ as the specific heat capacity of the water, $U_c$ the heat transfer coefficient, $V (= h \cdot A)$ the volume of the water in the tank and $Q_{out} (= \sqrt{h}/R)$ the outgoing flow at the bottom of the tank. The total amount of heat-energy of the water in the tank is $H= c_p \rho V T$ [in $J$].

The following quantities are important for setting-up the model equations:
-  $c_p \rho Q_{in} T_{in}$ : rate of heat flow due to the incoming water [$Js^{-1}$]
-  $c_p \rho Q_{out} T$ : rate of heat flow due to the outgoing water [$Js^{-1}$]
-  $U_c (A + \pi d h) (T - T_c)$ : rate of heat transfer due to the cooling jacket [$Js^{-1}$]

For the last quantity, we have assumed that the heat transfer only occurs at the side-walls and bottom where the water is in direct contact with the tank. Hence, no heat transfer at the water surface inside the tank.
"""

# ╔═╡ 1992f00e-c2a6-479b-a81b-d272eeb2dbed
md"""
!!! task
	Derive the model equations for the total heat-energy $H$ and the volume $V$ of the water inside the tank.
"""

# ╔═╡ 10c49e36-f058-4803-8dd5-a3ea2c22602d
md"""
!!! hint
	- Use the 'important' quantities mentioned above.
	- Set up:
	``` math
	\begin{align}
	H &= c_p \rho V T  \\
	V &= h A \\
	\frac{d H}{dt} &= \dots \\
	\frac{d V}{dt} &= \dots
	\end{align}
	```
	- Think well before putting a $+$ or $-$ sign before the terms while setting-up the differential equations!
"""

# ╔═╡ 27539311-aaf7-4b94-87cf-dfc4f4d2bcc6
solution(md"""
``` math
\begin{align}
	H &= c_p \rho V T  \\
	V &= h A \\
	\frac{d H}{dt} &= c_p \rho Q_{in} T_{in} - c_p \rho \cfrac{\sqrt{h}}{R} T - U_c (A + \pi d h) (T - T_c)\\
	\cfrac{dV}{dt} &= Q_{in} - \cfrac{\sqrt{h}}{R}
\end{align}
```
""")

# ╔═╡ 50bbe739-85fb-4403-b49a-2c210ab95701
md"""
## Setting up the equations
"""

# ╔═╡ de867d6d-b327-4e6f-a03d-9fb8f45dab7e
md"""
Define variables for the temperature $T$, the height $h$, the total heat-energy $H$ and the volume $V$ of the water in the tank. Use the following variable names: `T`, `h`, `H` and `V`. Mention the dependency on the time $t$.
"""

# ╔═╡ 3efd0cbd-293e-421d-b469-802486709698
# @variables missing

# ╔═╡ bf1b6e89-5c90-4327-9a26-066feb2a51e4
md"""
Define the parameters for this model and assign their corresponding values.

| Parameter | Value | Unit | Meaning |
|:---------- |:---------- |:------------|:------------|
| $A$    | 0.283 | $m^2$  |    cross sectional area   |
| $d$    | 0.6 | $m$  |    diameter of the tank   |
| $\rho$    | 1000.0 | $kg/m^3$  |  water density  |
| $Q_{in}$    | 0.01 | $m^3/s$  |   inlet flow |
| $T_{in}$    | 30.0 | $^{\circ}C$  |   temperature of the incoming water stream |
| $U_c$    | 4000 | $J/(s m^2\,^{\circ}C)$  |   heat transfer coefficient |
| $T_c$    | 15.0 | $^{\circ}C$  |   temperature of the cooling water |
| $c_p$    | 4200 | $J/(kg ^{\circ}C)$  |   heat capacity |
| $R$    | 200 | $s/m^{5/2}$  |   res. coeff. orifice |
"""

# ╔═╡ cb5b651d-4e26-4136-a8e2-cadbd75234db
# @parameters missing

# ╔═╡ 3f94b45e-284f-4d09-9b63-8bd11215f148
md"""
Set up the equation for total heat-energy, the volume, the rate of change in the total heat energy and the rate of change in volume.
"""

# ╔═╡ f8cdf897-6fca-41e6-92c9-ae2a7cc36eb8
eq_heat = missing

# ╔═╡ 2fae0b7b-1440-4adb-9531-1652b89d6f64
eq_volume = missing

# ╔═╡ 1c5707fa-cb02-49d4-bc8a-d17bb7055ac7
change_heat = missing

# ╔═╡ 70e23f70-aac4-46e6-b20c-18947e0153d7
change_volume = missing

# ╔═╡ 2d03afaf-55b9-4eab-8a41-5581468d97e1
md"""
Bundle the equations.
"""

# ╔═╡ 6052c938-48d5-4bf4-82ac-43eb42aca0ff
eqns_tank = missing

# ╔═╡ bc764684-0ef8-4925-b874-8479affd6d12
md"""
## Part 1: simple simulation
"""

# ╔═╡ 370763e7-ff6c-4713-ac52-c67df0607da0
md"""
### Building the ODE system
"""

# ╔═╡ 016fcba0-d53a-4dd0-abfe-22460a70b8f9
md"""
Build the model. Name it `sys1_tank`.
"""

# ╔═╡ 38cd75e9-4820-4e5e-a692-5d684f883c24
# @mtkbuild missing

# ╔═╡ 6909277e-7c4e-4d79-aaa0-3eeff6841998
md"""
### Create and solve the ODE problem
"""

# ╔═╡ c1d1f8a4-8934-43ff-9b35-8c291dfe016d
md"""
Create the ODE problem. Take as initial temperature of the water `15.0` $^{\circ}C$ and initial height of the water `0.8` $m$. Simulate for `2400.0`$s$.
"""

# ╔═╡ 3ab78b19-7d18-42d7-9962-11da4d46f5c7
oprob1_tank = missing

# ╔═╡ fcae1777-2c2f-4ad0-ae9c-a60ff6504de7
md"""
Solve the ODE problem. Use `saveat=1` and `reltol=1e-9`.
"""

# ╔═╡ bcb70e32-49ac-4d04-91cf-92f8f5d6548e
sol1_tank = missing

# ╔═╡ 9b0988e7-4546-4a98-a0fe-325b1e874f93
md"""
### Plotting results
"""

# ╔═╡ 9117e5c7-0617-4cec-8c29-5d362fa2eb16
md"""
Plot the temperature (with the option `idxs=[T]`) and the height (with the option `idxs=[h]`). Use `twinx()` as first argument in the second plot instruction to plot the height on the right y-axis.
"""

# ╔═╡ ff0cd780-ea3c-407f-a18a-9404ba208b68
# begin
# 	plot(missing)
# 	plot!(missing)
# end

# ╔═╡ 3cee6515-f9b8-4435-aa1a-bf93f976a618
md"""
!!! question
	Why do you think the temperature first rises, goes to a maximum and then decreases to go to an equilibrium value?
"""

# ╔═╡ b27bbc52-9131-467e-85c5-03f9b0ee2e78
md"""
Answers: missing
"""

# ╔═╡ 75d483ad-d983-487e-901c-40d3e530780d
md"""
## Part 2: effect of the resistance coefficient of the orifice
"""

# ╔═╡ afb92363-2ae7-4106-a829-d20bbf8c49da
md"""
In this part the resistance coefficient of the orifice will be increased to `250` at the time instant of `1200.0` $s$. In order to achieve this a discrete event will be used.
"""

# ╔═╡ 55cd8777-593d-4733-b740-ab4e0ab19724
md"""
### Building the ODE system
"""

# ╔═╡ b968af64-6cbb-4e01-b290-1eb68b3564f9
md"""
Build the new model that includes the discrete event. Name it `sys2_tank`.
"""

# ╔═╡ 09cab085-e5b4-4030-a027-d1ffe9cdbdc0
# @mtkbuild missing

# ╔═╡ f1965bf2-6985-453e-bdba-617245c6d705
md"""
### Create and solve the ODE problem
"""

# ╔═╡ f3bc477c-3044-413b-b1d5-e84f485106a3
md"""
Create the new ODE problem. Use the same initial values and simulation time as in the previous part.
"""

# ╔═╡ 291c2204-5026-425c-a431-2b27946cf72a
oprob2_tank = missing

# ╔═╡ d762e7f2-c030-49a4-865a-139d648e42b1
md"""
Solve the ODE problem. Use `saveat=1` and `reltol=1e-9`. Don't forget to take a `deepcopy` of the ODFE problem.
"""

# ╔═╡ 50719433-3e82-4614-ad19-b1610fe23916
sol2_tank = missing

# ╔═╡ 671696da-42ed-4d08-8234-aeea74e29a11
md"""
### Plotting results
"""

# ╔═╡ e5b36ae7-0a44-4c90-bc67-8f35e51d6ea7
md"""
Plot the temperature and the height. Use `twinx()` as first argument in the second plot instruction to plot the height on the right y-axis.
"""

# ╔═╡ 41a61694-551d-4034-879e-67066385158c
# begin
# 	plot()
# 	plot!()
# end

# ╔═╡ 7f9b877b-eab3-4df5-a5a6-55fee459f272
md"""
!!! question
	Is the evolution of $T$ and $h$ after the increase of $R$ according to your exceptations? Explain it.
"""

# ╔═╡ cf9f3ddf-f1ce-4bf8-b42b-18c53ccdc6d6
md"""
Answer: missing
"""

# ╔═╡ Cell order:
# ╟─0d12503e-5480-4ed8-9b4a-569b7801f855
# ╠═a6fbee1a-9863-11f0-1dfe-cf500206e163
# ╠═21c2d9ad-9411-4e05-be3e-7f37deebc70d
# ╠═eaa957df-d305-48fb-b2c2-4b4eb29a08be
# ╠═5f3f3184-7bf8-4f40-80ca-bb48ad85e2f9
# ╟─fdf5f147-5e6d-4730-8f16-ddc62ff9fcfe
# ╟─82207f70-ca9a-409b-952d-520f3b8babc3
# ╟─defe58ce-d407-4ab4-a4e4-e37fc061149c
# ╟─1992f00e-c2a6-479b-a81b-d272eeb2dbed
# ╟─10c49e36-f058-4803-8dd5-a3ea2c22602d
# ╟─27539311-aaf7-4b94-87cf-dfc4f4d2bcc6
# ╟─50bbe739-85fb-4403-b49a-2c210ab95701
# ╟─de867d6d-b327-4e6f-a03d-9fb8f45dab7e
# ╠═3efd0cbd-293e-421d-b469-802486709698
# ╟─bf1b6e89-5c90-4327-9a26-066feb2a51e4
# ╠═cb5b651d-4e26-4136-a8e2-cadbd75234db
# ╟─3f94b45e-284f-4d09-9b63-8bd11215f148
# ╠═f8cdf897-6fca-41e6-92c9-ae2a7cc36eb8
# ╠═2fae0b7b-1440-4adb-9531-1652b89d6f64
# ╠═1c5707fa-cb02-49d4-bc8a-d17bb7055ac7
# ╠═70e23f70-aac4-46e6-b20c-18947e0153d7
# ╟─2d03afaf-55b9-4eab-8a41-5581468d97e1
# ╠═6052c938-48d5-4bf4-82ac-43eb42aca0ff
# ╟─bc764684-0ef8-4925-b874-8479affd6d12
# ╟─370763e7-ff6c-4713-ac52-c67df0607da0
# ╟─016fcba0-d53a-4dd0-abfe-22460a70b8f9
# ╠═38cd75e9-4820-4e5e-a692-5d684f883c24
# ╟─6909277e-7c4e-4d79-aaa0-3eeff6841998
# ╟─c1d1f8a4-8934-43ff-9b35-8c291dfe016d
# ╠═3ab78b19-7d18-42d7-9962-11da4d46f5c7
# ╟─fcae1777-2c2f-4ad0-ae9c-a60ff6504de7
# ╠═bcb70e32-49ac-4d04-91cf-92f8f5d6548e
# ╟─9b0988e7-4546-4a98-a0fe-325b1e874f93
# ╟─9117e5c7-0617-4cec-8c29-5d362fa2eb16
# ╠═ff0cd780-ea3c-407f-a18a-9404ba208b68
# ╟─3cee6515-f9b8-4435-aa1a-bf93f976a618
# ╠═b27bbc52-9131-467e-85c5-03f9b0ee2e78
# ╟─75d483ad-d983-487e-901c-40d3e530780d
# ╟─afb92363-2ae7-4106-a829-d20bbf8c49da
# ╟─55cd8777-593d-4733-b740-ab4e0ab19724
# ╟─b968af64-6cbb-4e01-b290-1eb68b3564f9
# ╠═09cab085-e5b4-4030-a027-d1ffe9cdbdc0
# ╟─f1965bf2-6985-453e-bdba-617245c6d705
# ╟─f3bc477c-3044-413b-b1d5-e84f485106a3
# ╠═291c2204-5026-425c-a431-2b27946cf72a
# ╟─d762e7f2-c030-49a4-865a-139d648e42b1
# ╠═50719433-3e82-4614-ad19-b1610fe23916
# ╟─671696da-42ed-4d08-8234-aeea74e29a11
# ╟─e5b36ae7-0a44-4c90-bc67-8f35e51d6ea7
# ╠═41a61694-551d-4034-879e-67066385158c
# ╟─7f9b877b-eab3-4df5-a5a6-55fee459f272
# ╠═cf9f3ddf-f1ce-4bf8-b42b-18c53ccdc6d6
