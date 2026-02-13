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

# ╔═╡ 893b6948-b0be-11f0-13b3-f14c2a2208bb
using Pkg; Pkg.activate("..")

# ╔═╡ 47bc4a5d-7e5a-4278-b881-951737277dcf
using StatsPlots, PlutoUI; TableOfContents()

# ╔═╡ 9bb68aa2-14dd-4c28-b7c1-21f05f5671e2
using OrdinaryDiffEq, ModelingToolkit

# ╔═╡ fa5bb459-e754-401d-84a3-0ae75fb9b914
using ModelingToolkit: t_nounits as t, D_nounits as D

# ╔═╡ 00038789-08c2-464a-9941-1109ab47d811
md"""
# Exercise: Tractor seat
"""

# ╔═╡ f16b68e8-02ae-49b8-a077-c0100e2d3e89
solution(text) = Markdown.MD(Markdown.Admonition("hint", "Solution", [text]));

# ╔═╡ d6f7fe5d-155d-434c-a0c3-6e6e9aae55ba
md"""
![](https://users.ugent.be/~gvhaelew/fig/tractor_seat_all.png)
"""

# ╔═╡ 40b40bcf-498b-46f8-8f5d-9eb24375d229
md"""
## Deriving the model equation
"""

# ╔═╡ 815032ea-12ad-446c-836c-d75e2724bed4
md"""
Find out how the vertical position $y$ of a tractor seat with mass $m$ varies when, after 1 second, a person with mass $M$ suddenly sits in the seat. The tractor seat consists of a spring and a shock absorber. The spring constant is denoted $k$ and the natural length of the spring is denoted $L$. The shock absorber has a damping effect which is proportional (proportionality factor $b$), but opposite, to the velocity of the seat when moving. Pay attention to the sense of forces! Assume that the seat was initially in a rest position $y_0$ and not moving ($v_0 = 0$).
"""

# ╔═╡ 00c93c86-f862-4999-8b54-4edeb85628a0
md"""
!!! task
	Derive the model equation for the position $y$ of the seat. Use Newton's second law.
"""

# ╔═╡ 86940e79-05d5-43e5-afb1-590475774850
md"""
!!! hint
	- Set up:
	``` math
	\begin{align}
	m\cfrac{d^2 y}{dt^2} &= \cdots \\
	\end{align}
	```
"""

# ╔═╡ 5210f34f-928a-4a14-ab61-e2a9454508dd
solution(md"""
``` math
\begin{align}
	m\cfrac{d^2 y}{dt^2} &= -mg + k\left(L - y\right) - bv
\end{align}
```
""")

# ╔═╡ 5895728b-53cb-4717-8300-aedb1694c0cb
md"""
## Setting up the equations
"""

# ╔═╡ 20073595-1940-406f-b30b-57a1033e9c4f
md"""
Define the variable for the position $y$. Use the following variable name: `y`. Mention the dependency on the time $t$.
"""

# ╔═╡ eb2ac37e-91e7-4253-9274-ef8e3a2c666b
# @variables missing

# ╔═╡ 16b70e7d-e9a4-4406-a830-6e591fd868e7
md"""
Define the parameters for this model and assign their corresponding values.

| Parameter | Value | Unit | Meaning |
|:---------- |:---------- |:------------|:------------|
| $m$    | 23.0 | $kg$  |    mass of the seat   |
| $g$    | 9.81 | $m/s^2$  |   gravitational constant |
| $k$    | 22570 | $N/m$  |  spring constant  |
| $L$    | 0.40 | $m$  |   natural length of the spring |
| $b$    | 900 | $Ns/m$  |   damping constant |

Use the following names: `m`, `g`, `k`, `L` and `b`.
"""

# ╔═╡ 0472f2b5-394f-4670-81ec-3046c9826ee0
# @parameters missing

# ╔═╡ e68827a7-5f6e-45f7-b7fb-d0b41472b2c8
md"""
In order to know the initial position $y_0$ of the seat, calculate the steady state of $y$.

Hints:
- Set both $\cfrac{d^2 y}{dt^2}$ and $v = \cfrac{d y}{dt}$ to zero in the model equation and determine an expression for the steady state of $y$.
- Once you have the expression it is preferable to calculate the actual value it a `let` ... `end` block. Everything in this block will be local scope and won't interfere with variable names defined elsewhere.
- Round to two digits after the decimal point with the function `round(..., digits=2)`.
"""

# ╔═╡ d8ab108b-89af-4f3d-bb4d-dae83bd71aaf
# let
# 	m=23.0; g=9.81; k=22570; L=0.40
# 	missing
# end

# ╔═╡ 42a05adc-a5ac-43d5-9e0b-d8cb1f539c09
md"""
The person that will sit in the seat has a mass of `80.0` $kg$. Use a normal variable for the parameter `M`, the mass of the person.
"""

# ╔═╡ be8bff08-aa49-4ef4-9be6-d13a8aac01c5
missing

# ╔═╡ e918256e-3999-4ab2-8792-7da6ba32ab19
md"""
Set up the model equation.
"""

# ╔═╡ 115c48a4-4dc0-4cf0-ac34-e651e773e675
eq_position = missing

# ╔═╡ 6049886b-95a3-4d9e-a4e7-d05f1e9c504a
md"""
## Building the ODE system
"""

# ╔═╡ 5856753b-2753-4a83-b3ec-747a7c0b9c52
md"""
Build the model and include the (discrete) event that when the time is `1`, `m` needs to be incremented with `M`. Name your model `sys_tractor_seat`.
"""

# ╔═╡ 75e757ef-ed89-40f7-a8d6-12884f4ccb98
# @mtkbuild missing

# ╔═╡ d4bc6ff0-7383-4ea2-a7e0-5719a9dfe794
md"""
## Create and solve the ODE problem
"""

# ╔═╡ 9270475c-1f63-4d3e-98e0-ba7d2417ee81
md"""
Create the ODE problem. Use the steady state value that you calculated before for the initial position $y_0$. As mentioned before $v_0 = 0$. The symbol for $\frac{dy}{dt}$ in the vector of initial conditions is `D(y)`. Use a simulation time span of `5.0` seconds. Use `[b=>b_val]` for the parameters argument so that you can see the effect of the damping coefficient on the results. The variable `b_val` is defined later and bound to a slider.
"""

# ╔═╡ db295a4e-e877-46c6-9f7d-b3b9fa8d04d4
oprob_treactor_seat = missing

# ╔═╡ a6ca51e5-4019-4aa1-a667-66af70b77f8a
md"""
Solve the ODE problem. Make a deepcopy of the ODE problem, use `Tsit5()` and `saveat=0.01`.
"""

# ╔═╡ 7c1edab7-722d-4014-90c0-1b834e80aec7
sol_tractor_seat = missing

# ╔═╡ 195752fc-f480-4996-98ae-41c0daf74942
# @bind b_val Slider(200:50:4000, default=900, show_value=true)

# ╔═╡ dbe1a845-2cdc-45c0-9966-3187e67ef296
md"""
## Plotting results
"""

# ╔═╡ a4dbf8e2-044a-427a-b773-9b5b6281bffc
md"""
Plot the position $y$ of the seat over time. Use the option `idxs=[y]` and `ylim=(0.32, 0.392)`. Play with the sliders above to see their effect.
"""

# ╔═╡ 0861aabb-60b5-442b-bb0d-a26ed10e7c00
missing

# ╔═╡ 0d90986f-7c35-4145-aa50-9093725b1440
md"""
Check out the final value of $y$.
"""

# ╔═╡ d8db3d47-7a0d-4d59-8820-5dc81475b5e3
missing

# ╔═╡ 37c21c15-7e1e-4dc5-9cce-eaacba2be5d4
md"""
!!! questions
	1. Can you interprete the plot?
	2. What is the approximate value of $b$ in order to have a critically damped motion? (A critically damped motion is when a system returns to its equilibrium position as quickly as possible without oscillating or overshooting).
	3. Does the steady state value depend on the value of $b$?
"""

# ╔═╡ b7b5a8d9-5557-4fad-823c-1fb306894605
md"""
Answers:
1. missing
2. missing
3. missing
"""

# ╔═╡ 12e1b09d-7ec4-4efe-9a27-3373a8e22116
md"""
Plot the velocity $v$ of the seat over time. Use `idxs=[D(y)]`.
"""

# ╔═╡ c4e74bfe-c266-4040-81b4-f5d60737d406
missing

# ╔═╡ 54196cf6-8a5c-4db6-ab16-209c447316b0
md"""
!!! question
	If you use a considerable value for $b$, will the seat still be moving at the end of the time span?
"""

# ╔═╡ 24aeaf8b-fc8a-45d2-9972-f4866aa6d17c
md"""
Answer: missing
"""

# ╔═╡ b80c35e7-437e-4377-aee3-31294c680f30
md"""
## Calculating the steady state values
"""

# ╔═╡ 571afefb-1655-40e9-bfd2-d875db112fbe
md"""
Create a steady state problem and solve it. Hint: set the mass `m` to the correct value!
"""

# ╔═╡ 5a07d780-0526-482b-a445-fd21c86a3104
stst_val = missing

# ╔═╡ 8e9e6d22-641c-41fc-aef6-4aabbccbf573
md"""
Show the steady state value for `y` and `D(y)`.
"""

# ╔═╡ 8b0c836a-5bb3-4490-8a22-b7ed2fa97aaa
missing

# ╔═╡ 6d050df4-eae3-4aef-a730-ef343cee812d
missing

# ╔═╡ edd73175-b75f-4538-b189-8049abfe5abb
md"""
!!! question
	Do the steady state values correspond to the ones you can derive from the plots?
"""

# ╔═╡ b2e0a8ee-fddc-4c8f-af07-92fb5ec9aa34
md"""
Answer: missing
"""

# ╔═╡ Cell order:
# ╟─00038789-08c2-464a-9941-1109ab47d811
# ╠═893b6948-b0be-11f0-13b3-f14c2a2208bb
# ╠═47bc4a5d-7e5a-4278-b881-951737277dcf
# ╠═9bb68aa2-14dd-4c28-b7c1-21f05f5671e2
# ╠═fa5bb459-e754-401d-84a3-0ae75fb9b914
# ╠═f16b68e8-02ae-49b8-a077-c0100e2d3e89
# ╟─d6f7fe5d-155d-434c-a0c3-6e6e9aae55ba
# ╟─40b40bcf-498b-46f8-8f5d-9eb24375d229
# ╟─815032ea-12ad-446c-836c-d75e2724bed4
# ╟─00c93c86-f862-4999-8b54-4edeb85628a0
# ╟─86940e79-05d5-43e5-afb1-590475774850
# ╟─5210f34f-928a-4a14-ab61-e2a9454508dd
# ╟─5895728b-53cb-4717-8300-aedb1694c0cb
# ╟─20073595-1940-406f-b30b-57a1033e9c4f
# ╠═eb2ac37e-91e7-4253-9274-ef8e3a2c666b
# ╟─16b70e7d-e9a4-4406-a830-6e591fd868e7
# ╠═0472f2b5-394f-4670-81ec-3046c9826ee0
# ╟─e68827a7-5f6e-45f7-b7fb-d0b41472b2c8
# ╠═d8ab108b-89af-4f3d-bb4d-dae83bd71aaf
# ╟─42a05adc-a5ac-43d5-9e0b-d8cb1f539c09
# ╠═be8bff08-aa49-4ef4-9be6-d13a8aac01c5
# ╟─e918256e-3999-4ab2-8792-7da6ba32ab19
# ╠═115c48a4-4dc0-4cf0-ac34-e651e773e675
# ╟─6049886b-95a3-4d9e-a4e7-d05f1e9c504a
# ╟─5856753b-2753-4a83-b3ec-747a7c0b9c52
# ╠═75e757ef-ed89-40f7-a8d6-12884f4ccb98
# ╟─d4bc6ff0-7383-4ea2-a7e0-5719a9dfe794
# ╟─9270475c-1f63-4d3e-98e0-ba7d2417ee81
# ╠═db295a4e-e877-46c6-9f7d-b3b9fa8d04d4
# ╟─a6ca51e5-4019-4aa1-a667-66af70b77f8a
# ╠═7c1edab7-722d-4014-90c0-1b834e80aec7
# ╠═195752fc-f480-4996-98ae-41c0daf74942
# ╟─dbe1a845-2cdc-45c0-9966-3187e67ef296
# ╟─a4dbf8e2-044a-427a-b773-9b5b6281bffc
# ╠═0861aabb-60b5-442b-bb0d-a26ed10e7c00
# ╟─0d90986f-7c35-4145-aa50-9093725b1440
# ╠═d8db3d47-7a0d-4d59-8820-5dc81475b5e3
# ╟─37c21c15-7e1e-4dc5-9cce-eaacba2be5d4
# ╠═b7b5a8d9-5557-4fad-823c-1fb306894605
# ╟─12e1b09d-7ec4-4efe-9a27-3373a8e22116
# ╠═c4e74bfe-c266-4040-81b4-f5d60737d406
# ╟─54196cf6-8a5c-4db6-ab16-209c447316b0
# ╠═24aeaf8b-fc8a-45d2-9972-f4866aa6d17c
# ╟─b80c35e7-437e-4377-aee3-31294c680f30
# ╟─571afefb-1655-40e9-bfd2-d875db112fbe
# ╠═5a07d780-0526-482b-a445-fd21c86a3104
# ╟─8e9e6d22-641c-41fc-aef6-4aabbccbf573
# ╠═8b0c836a-5bb3-4490-8a22-b7ed2fa97aaa
# ╠═6d050df4-eae3-4aef-a730-ef343cee812d
# ╟─edd73175-b75f-4538-b189-8049abfe5abb
# ╠═b2e0a8ee-fddc-4c8f-af07-92fb5ec9aa34
