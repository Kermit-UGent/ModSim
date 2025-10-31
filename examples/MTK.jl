### A Pluto.jl notebook ###
# v0.20.10

using Markdown
using InteractiveUtils

# ╔═╡ 1c46508e-2354-11f0-08db-d373e43929bd
using Pkg; Pkg.activate("..")

# ╔═╡ 2d125aad-334c-4210-a921-1acacefc5cfa
using Plots, PlutoUI

# ╔═╡ d0f32979-01f6-4f07-acb3-627077c7c029
using DifferentialEquations, ModelingToolkit

# ╔═╡ 43c4f402-45a9-48d5-b1c6-8a6ea8107ce2
using ModelingToolkit: t_nounits as t, D_nounits as D

# ╔═╡ 3ce7be2a-51e8-419c-b8fb-764c15658903
md"This is a barebones tutorial to ModelingToolkit, helping you to build simple models based on ODEs. See the [documentations](https://docs.sciml.ai/ModelingToolkit/stable/) for more details."

# ╔═╡ f4db6481-80c9-4d56-b5a4-50e620d95143
md"""
## Case 1: conical tank

![](https://i.ibb.co/7xYxhD7f/conical-tank.png)
"""

# ╔═╡ 41af0644-4611-4880-950d-0b6ff2c58a15
@variables h(t) r(t) V(t)

# ╔═╡ 38602497-c591-4e22-af60-f9622096067f
@parameters k q H R

# ╔═╡ f98ecb14-9df2-4ba3-a1c7-3ab75e3e8055
eq_radius = r ~ h * R / H 


# ╔═╡ c18b1610-eadf-474e-8da9-9c1f21845bc1
eq_volume = V ~ π * r^2 * h / 3

# ╔═╡ e87ea53c-8629-4f49-af3b-1c753940e27d
change_vol = D(V) ~ q - h * k   # change in volume = inflow - outflow

# ╔═╡ 8e6efab4-ac91-40bd-9805-88e0dabb2cd1
eqns_tank = [eq_radius, eq_volume, change_vol]

# ╔═╡ c6bfc46c-7ca5-4d9e-88fb-3e7065ec6bb0
@mtkbuild tank = ODESystem(eqns_tank, t)

# ╔═╡ d6a153d4-d41f-47d4-88b0-b10ac2e22b04
tank_prob = ODEProblem(tank, [V=>0.1], (0, 100), [q=>0.1, H=>5, R=>2, k=>0.05], 				guesses=[r=>0.1, h=>0.1])

# ╔═╡ 9492a458-2e7c-4828-85a3-20ed266e794c
tank_sol = solve(tank_prob)

# ╔═╡ 66947904-6b80-4279-ae2b-a3c20886062d
length(tank_sol)

# ╔═╡ 99158cd8-753c-4d6b-9288-a78dfc1fa86a
plot(tank_sol, idxs=[V, r, h])

# ╔═╡ 8c9e806b-ee8b-4c84-b9bb-113e1d8ab4e4
plot(tank_sol, idxs=V/h)

# ╔═╡ b108c4a6-ee94-4773-8b88-fe02b70addfb
md"""
## Case 2: Mass-and-spring system

![](https://i.ibb.co/MDTGngfq/damped-spring.png)
"""

# ╔═╡ a0f5eeee-3b5f-4588-bd31-1cff6c993f6b
@variables y(t) E(t)  # position and energy

# ╔═╡ d3954593-63ae-4dfe-b752-01339b44adcd
@parameters μ [description="friction"]  m [description="mass"]

# ╔═╡ 268ab13e-2f15-4b1e-90b7-cbbcf5496163
D^2  # second-order derivative wrt time

# ╔═╡ ed48e70a-a608-46a1-9339-0bb3a93a4ea1
spring_eq = m * D(D(y)) ~ - k * y - μ * D(y)

# ╔═╡ dc7e26c6-2c88-475d-a4ab-64dbb4f9a956
energy = E ~ m*D(y)^2 /2 + y^2*k/2  # we can add quantities to keep track off. 

# ╔═╡ 75a3db3f-818a-4cb8-948a-dfbcd0582a45
@mtkbuild spring_ode = ODESystem([spring_eq, energy], t)

# ╔═╡ 30b3ed77-0709-4410-8d56-d576d69b5e0b
spring_prob = ODEProblem(spring_ode, [y=>2.0, D(y)=>-1.0], (0.0, 100.), [m=>3, k=>0.6, μ=>1e-1])

# ╔═╡ 15ba5472-72da-4ae4-9438-fe8a1d88d9c0
sol_spring = solve(spring_prob, Tsit5(), reltol=1e-9);

# ╔═╡ d29f022c-924b-4d26-ab8c-65f11eb85b6a
plot(sol_spring)

# ╔═╡ ff3645ab-12d1-46ac-b4d8-ce431ba27491
plot(sol_spring, idxs=E)

# ╔═╡ 27670f6a-9f2c-4768-a1ea-986f9ba5301c
md"""

## Case 3: Lotka-Volterra

Classical model for prey-predator relations
"""

# ╔═╡ 3100925e-d5dc-4d56-af84-ffab2059a4c5
@variables 🐰(t) 🦊(t)

# ╔═╡ 298c0f1b-4237-4c17-993b-a2969257790f
@parameters α β γ δ K

# ╔═╡ 130f2711-9ac8-491f-b580-28fe0889257a
LV_eqs = [
	D(🐰) ~ α * 🐰 * (1 - 🐰 / K) - β * 🐰 * 🦊,
	D(🦊) ~ γ *  🐰 * 🦊 - δ * 🦊 
]

# ╔═╡ 1a920b39-00f6-4455-b226-b6961af5c6d1
foxmanagment = [🐰 ~ 100.0] => [🦊 ~ 🦊+1] 

# ╔═╡ 1cdf417d-daea-4684-89bd-7deb951c6b3a
@mtkbuild lv_sys = ODESystem(LV_eqs, t, continuous_events=foxmanagment)

# ╔═╡ 08bfcddc-703f-4ef4-a44e-6386eb91de76
lv_prob = ODEProblem(lv_sys, [🐰=>1.0, 🦊=>1e-2], (0, 100), [α=>0.6, β=>0.6, γ=>0.04, δ=>0.5, K=>1000])

# ╔═╡ 5c927b82-1133-4508-996a-7ed7111cc4d5
sol_LV = solve(lv_prob, Tsit5());

# ╔═╡ 9db55fe4-9f67-4d23-85ec-fbafe21f09a7
plot(sol_LV)

# ╔═╡ cc65899f-f16c-4194-b2d8-bcfe3dcbd8c8
plot(sol_LV, idxs = 🐰 / 🦊, label="rabits/foxes")

# ╔═╡ 2ca7379b-6b0b-488b-b760-b969014ec110
TableOfContents()

# ╔═╡ Cell order:
# ╠═1c46508e-2354-11f0-08db-d373e43929bd
# ╠═2d125aad-334c-4210-a921-1acacefc5cfa
# ╠═d0f32979-01f6-4f07-acb3-627077c7c029
# ╠═43c4f402-45a9-48d5-b1c6-8a6ea8107ce2
# ╟─3ce7be2a-51e8-419c-b8fb-764c15658903
# ╟─f4db6481-80c9-4d56-b5a4-50e620d95143
# ╠═41af0644-4611-4880-950d-0b6ff2c58a15
# ╠═38602497-c591-4e22-af60-f9622096067f
# ╠═f98ecb14-9df2-4ba3-a1c7-3ab75e3e8055
# ╠═c18b1610-eadf-474e-8da9-9c1f21845bc1
# ╠═e87ea53c-8629-4f49-af3b-1c753940e27d
# ╠═8e6efab4-ac91-40bd-9805-88e0dabb2cd1
# ╠═c6bfc46c-7ca5-4d9e-88fb-3e7065ec6bb0
# ╠═d6a153d4-d41f-47d4-88b0-b10ac2e22b04
# ╠═9492a458-2e7c-4828-85a3-20ed266e794c
# ╠═66947904-6b80-4279-ae2b-a3c20886062d
# ╠═99158cd8-753c-4d6b-9288-a78dfc1fa86a
# ╠═8c9e806b-ee8b-4c84-b9bb-113e1d8ab4e4
# ╟─b108c4a6-ee94-4773-8b88-fe02b70addfb
# ╠═a0f5eeee-3b5f-4588-bd31-1cff6c993f6b
# ╠═d3954593-63ae-4dfe-b752-01339b44adcd
# ╠═268ab13e-2f15-4b1e-90b7-cbbcf5496163
# ╠═ed48e70a-a608-46a1-9339-0bb3a93a4ea1
# ╠═dc7e26c6-2c88-475d-a4ab-64dbb4f9a956
# ╠═75a3db3f-818a-4cb8-948a-dfbcd0582a45
# ╠═30b3ed77-0709-4410-8d56-d576d69b5e0b
# ╠═15ba5472-72da-4ae4-9438-fe8a1d88d9c0
# ╠═d29f022c-924b-4d26-ab8c-65f11eb85b6a
# ╠═ff3645ab-12d1-46ac-b4d8-ce431ba27491
# ╟─27670f6a-9f2c-4768-a1ea-986f9ba5301c
# ╠═3100925e-d5dc-4d56-af84-ffab2059a4c5
# ╠═298c0f1b-4237-4c17-993b-a2969257790f
# ╠═130f2711-9ac8-491f-b580-28fe0889257a
# ╠═1a920b39-00f6-4455-b226-b6961af5c6d1
# ╠═1cdf417d-daea-4684-89bd-7deb951c6b3a
# ╠═08bfcddc-703f-4ef4-a44e-6386eb91de76
# ╠═5c927b82-1133-4508-996a-7ed7111cc4d5
# ╠═9db55fe4-9f67-4d23-85ec-fbafe21f09a7
# ╠═cc65899f-f16c-4194-b2d8-bcfe3dcbd8c8
# ╠═2ca7379b-6b0b-488b-b760-b969014ec110
