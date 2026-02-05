### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ e079823b-8b40-42a2-a63f-1645a97b33f0
using Pkg; 	Pkg.activate("..")

# ╔═╡ 62dc7706-f58a-11ee-2d3d-f78f7ceca914
using Markdown; InteractiveUtils;

# ╔═╡ 65571bb8-e260-4a82-b0d2-198e47c56271
using Plots, PlutoUI

# ╔═╡ 5e79d359-393b-405a-847e-dbdec15a1848
using DifferentialEquations, ModelingToolkit

# ╔═╡ cd21ad3c-94b8-4031-8fc0-17a5a39043f4
using ModelingToolkit: t_nounits as t, D_nounits as D

# ╔═╡ 66e8a12c-74b6-4077-b90e-3d85e5a61d6e
solution(text) = Markdown.MD(Markdown.Admonition("hint", "Solution", [text]));

# ╔═╡ f55e0ca3-07ff-43a1-907f-9ad640227822
md"""
# Exercise: Water evaporation and infiltration

Consider a water reservoir, such as a lake, where the water in the reservoir is in contact with the air as well as with the groundwater. We will denote the water level in the reservoir as $W$ and the groundwater level as $G$.

The water in the reservoir evaporates at a rate $k_1$ (i.e. the *evaporation coefficient*) and there can be infiltration into or from the groundwater at a rate $k_2$ (i.e., *infiltration coefficient*) depending on the difference in the water level in the reservoir and groundwater (cf. $(W-G)$)

There is a natural constant inflow of water into the reservoir at a rate $I$. At time $t=0\;s$ a pumping device is switched on such that the reservoir is rapidly being emptied at an outflow rate $O$ until the level of the water reservoir drops to zero. From then on, the pump is switched off.
"""

# ╔═╡ a551f3c5-2fc3-4236-ab4a-3d9c14afc62b
md"""
!!! task
	Set-up a system of differential equations modelling the above problem.
"""

# ╔═╡ 4eb95688-ba9d-4524-a522-3b8343f4e2be
solution(
md"""
The system of differential equations that models the water level in a reservoir ($W$) and the groundwater level ($G$) considering evaporation, infiltration, inlet flow and outlet flow can be written down as:

$$\begin{align}
\frac{dW}{dt} & = I - O - k_1 \cdot W - k_2 \cdot (W - G) \\
\frac{dG}{dt} & = k_2 \cdot (W - G)
\end{align}$$
"""
)

# ╔═╡ f7373868-80fc-4838-bbc2-bc1593a86d6e
md"""
Model the aforementioned system of differential equations using ModelingToolkit.
"""

# ╔═╡ 37faba08-5b07-4471-bd4c-296a9c371100
md"""
Define the variables with `@variables`.
"""

# ╔═╡ f86d8696-aa00-463d-a421-ea6a738e5bf8
missing

# ╔═╡ c08c6d18-400d-46db-bb5e-23b4cfe6f2db
md"""
Define the parameters with `@parameters`.
"""

# ╔═╡ db7a6a7f-3e46-433a-9eaf-b9c1061958e0
missing

# ╔═╡ e0145073-e173-4d7a-be0b-aa270c9b57dc
md"""
Set up the equations. Name the set of equations `eqs_wei`.
"""

# ╔═╡ 6900fb0a-430a-411e-a9b4-bf5ceaff36c6
eqs_wei = missing

# ╔═╡ 1e3c0058-a5da-47b8-b37c-f2820d20fe41
md"""
## Part 1
"""

# ╔═╡ 0c107a12-0ec3-4225-95ba-9afecdcf5b61
md"""
There is a natural constant inflow of water into the reservoir at a rate $I$. At time $t=0\;s$ a pumping device is switched on such that the reservoir is rapidly being emptied at an outflow rate $O$ until the level of the water reservoir drops to zero. From then on, the pump is switched off.

The latter can be handled using the continuous event. Use the following option when creating an `ODESystem`: `continuous_events = [[W ~ 0]=>[O ~ 0]]`. This will make sure that when $W$ hits zero, then $O$ will be put to zero.
"""

# ╔═╡ b49dc322-2171-4890-b070-ee83593b9336
md"""
Build a system of equations with `@mtkbuild`. Name it `sys1_wei`.
"""

# ╔═╡ b7e1dc94-c1b5-4dab-ba9b-f8b9a8fe1898
missing

# ╔═╡ acab2bf0-b792-4ccc-bee0-7611bedab23c
md"""
Both water levels are initially $6.75\;m$. The inflow rate is constant and is $2.7\; m/min$. The evaporation and infiltration coefficient are $0.4\;min^{-1}$ and $1.0\;min^{-1}$ respectively. The outflow rate due to the pump is $20\;m/min$ and the pump stops working when $W$ equals zero. We wish to simulate the evolution of $W$ and $G$ during $20\;min$.
"""

# ╔═╡ 44cf2ac6-75e2-4440-94f1-8fe6887ee1e0
md"""
Initialize a vector `u0` with the initial conditions:
"""

# ╔═╡ 733bdb56-fb4f-4bc6-b50c-e3245fd59730
u0 = missing

# ╔═╡ 521fcca5-3c88-463b-9395-e5871b9fc5a3
md"""
Set the timespan for the simulation:
"""

# ╔═╡ f2267cc9-4bfd-44af-984c-cf54aa855f91
tspan = missing

# ╔═╡ e12840ed-d448-4f7c-893d-cfe974f62f9a
md"""
Initialize a vector `parms` with the parameter values:
"""

# ╔═╡ 5a3d520d-fd4f-485c-a8b0-63858eca4bfc
parms = missing

# ╔═╡ aafb49a4-b468-49b2-838d-5ddfcc852d48
md"""
Create the ODE problem and store it in `oprob1_wei`:
"""

# ╔═╡ a25d3652-13f5-47ef-9f16-c6698547a734
oprob1_wei = missing

# ╔═╡ 451a3c66-5bcc-4161-af81-f89af33b5862
md"""
Solve the ODE problem. Make a deepcopy and use `Tsit5()`, `saveat=0.1` and `reltol=1e-9`. Store the solution in `osol1_wei`:
"""

# ╔═╡ b9d2c6cb-88f0-4a88-9e61-eecc905ff3e6
osol1_wei = missing

# ╔═╡ ab77c284-d379-47d1-bd86-88fa77749165
md"""
Plot the results:
"""

# ╔═╡ e29a3294-245b-445d-bb1e-12cafb2ec175
missing

# ╔═╡ 66588291-d399-4169-9383-ac6c05cdf906
md"""
Interpret the results. Ask yourself the following questions:

1. Can you clearly see the drop in $W$? To what value does $W$ drops?
"""

# ╔═╡ 89ece6c0-3690-4265-8b1d-c3a3cc8b095f
md"- Answer: missing"

# ╔═╡ 92fa2a5a-123b-434c-a97f-39d9727a5ab6
md"""
2. Why does $G$ also drop when $W$ drops? Explain.
"""

# ╔═╡ 71dde4bd-711c-40c3-8ee2-9f7a79d72aa2
md"- Answer: missing"

# ╔═╡ d04d906d-1432-4753-81b9-b03a06df9816
md"""
3. To what values are $W$ and $G$ tending to go? Was the system with the initial values for $W$ and $G$ and no outflow in equilibrium? Explain.
"""

# ╔═╡ 4ce3328e-3436-45ad-b899-9b901b53a8ea
md"- Answer: missing"

# ╔═╡ 97fc5a24-f5bf-4f62-b3c6-d3b3be1f7b71
md"""
## Part 2
"""

# ╔═╡ 21560a84-db0e-495c-bb87-ae8a81a8755b
md"""
Copy the above system of equations that you have built with `@mtkbuild`. Name it `sys2_wei`. Extend this system with two discrete events so that the pump is switched on ($20\;L/min$) at time 10 and time 15.

Use the following additional option when creating an `ODESystem`:\
`discrete_events = [[10]=>[O ~ 20], [15]=>[O ~ 20]]`.
"""

# ╔═╡ ef66afba-119f-4b07-ae18-875188d5caab
missing

# ╔═╡ 071425d1-a012-4bbc-99ce-e2a4137988e2
md"""
We will use the same values for the initial conditions, time span and parameter, so we don't need to redefine them.
"""

# ╔═╡ d08aea59-c2f3-4fcd-b201-bff8999cb7ef
md"""
Create the ODE problem and store it in `oprob2_wei`:
"""

# ╔═╡ 8e855e0d-58e7-4413-8124-e8ea1e092dd6
oprob2_wei = missing

# ╔═╡ 07f5a55b-f79b-41f9-8568-23bcfa325090
md"""
Solve the ODE problem. Make a deepcopy and use `Tsit5()`, `saveat=0.1` and `reltol=1e-9`. Store the solution in `osol2_wei`:
"""

# ╔═╡ 37cf1c26-41d6-4f89-922b-4d5125fd1001
osol2_wei = missing

# ╔═╡ 15482456-5b18-4453-87ab-57bf6a75d1de
md"""
Plot the results:
"""

# ╔═╡ fe8cc5dc-4e47-46af-b18c-a56812338988
missing

# ╔═╡ 47822a88-4a89-4466-870f-207f825fb82f
md"""
Interpret the results.
"""

# ╔═╡ e1355c7b-9d8a-4908-abc5-f0a5d0e8c1f4
md"""
missing
"""

# ╔═╡ 35bed614-5c58-4690-9341-276228561436
md"""
## Part 3
"""

# ╔═╡ a55e40f6-dda7-4705-a37b-191ab7b1f16c
md"""
Solve for the equilibrium values using a `SteadyStateProblem`. You can using either `sys1_wei` or `sys2_wei` to do that. Provide a vector with an initial guess for the equibibrium values and vector for the parameter values where $O$ is zero (cf. pump is switched off)!
"""

# ╔═╡ ad03bfce-ca51-4727-934b-dbfdc42e9443
eq_val = missing

# ╔═╡ a86346a2-0d31-4c25-88fb-b2ed3c3bd425
Weq = missing

# ╔═╡ de29ebf9-01fb-451e-8100-cc171c7b4843
Geq = missing

# ╔═╡ Cell order:
# ╠═62dc7706-f58a-11ee-2d3d-f78f7ceca914
# ╠═e079823b-8b40-42a2-a63f-1645a97b33f0
# ╠═65571bb8-e260-4a82-b0d2-198e47c56271
# ╠═5e79d359-393b-405a-847e-dbdec15a1848
# ╠═cd21ad3c-94b8-4031-8fc0-17a5a39043f4
# ╟─66e8a12c-74b6-4077-b90e-3d85e5a61d6e
# ╟─f55e0ca3-07ff-43a1-907f-9ad640227822
# ╟─a551f3c5-2fc3-4236-ab4a-3d9c14afc62b
# ╟─4eb95688-ba9d-4524-a522-3b8343f4e2be
# ╟─f7373868-80fc-4838-bbc2-bc1593a86d6e
# ╟─37faba08-5b07-4471-bd4c-296a9c371100
# ╠═f86d8696-aa00-463d-a421-ea6a738e5bf8
# ╟─c08c6d18-400d-46db-bb5e-23b4cfe6f2db
# ╠═db7a6a7f-3e46-433a-9eaf-b9c1061958e0
# ╟─e0145073-e173-4d7a-be0b-aa270c9b57dc
# ╠═6900fb0a-430a-411e-a9b4-bf5ceaff36c6
# ╟─1e3c0058-a5da-47b8-b37c-f2820d20fe41
# ╟─0c107a12-0ec3-4225-95ba-9afecdcf5b61
# ╟─b49dc322-2171-4890-b070-ee83593b9336
# ╠═b7e1dc94-c1b5-4dab-ba9b-f8b9a8fe1898
# ╟─acab2bf0-b792-4ccc-bee0-7611bedab23c
# ╟─44cf2ac6-75e2-4440-94f1-8fe6887ee1e0
# ╠═733bdb56-fb4f-4bc6-b50c-e3245fd59730
# ╟─521fcca5-3c88-463b-9395-e5871b9fc5a3
# ╠═f2267cc9-4bfd-44af-984c-cf54aa855f91
# ╟─e12840ed-d448-4f7c-893d-cfe974f62f9a
# ╠═5a3d520d-fd4f-485c-a8b0-63858eca4bfc
# ╟─aafb49a4-b468-49b2-838d-5ddfcc852d48
# ╠═a25d3652-13f5-47ef-9f16-c6698547a734
# ╟─451a3c66-5bcc-4161-af81-f89af33b5862
# ╠═b9d2c6cb-88f0-4a88-9e61-eecc905ff3e6
# ╟─ab77c284-d379-47d1-bd86-88fa77749165
# ╠═e29a3294-245b-445d-bb1e-12cafb2ec175
# ╟─66588291-d399-4169-9383-ac6c05cdf906
# ╠═89ece6c0-3690-4265-8b1d-c3a3cc8b095f
# ╟─92fa2a5a-123b-434c-a97f-39d9727a5ab6
# ╠═71dde4bd-711c-40c3-8ee2-9f7a79d72aa2
# ╟─d04d906d-1432-4753-81b9-b03a06df9816
# ╠═4ce3328e-3436-45ad-b899-9b901b53a8ea
# ╟─97fc5a24-f5bf-4f62-b3c6-d3b3be1f7b71
# ╟─21560a84-db0e-495c-bb87-ae8a81a8755b
# ╠═ef66afba-119f-4b07-ae18-875188d5caab
# ╟─071425d1-a012-4bbc-99ce-e2a4137988e2
# ╟─d08aea59-c2f3-4fcd-b201-bff8999cb7ef
# ╠═8e855e0d-58e7-4413-8124-e8ea1e092dd6
# ╟─07f5a55b-f79b-41f9-8568-23bcfa325090
# ╠═37cf1c26-41d6-4f89-922b-4d5125fd1001
# ╟─15482456-5b18-4453-87ab-57bf6a75d1de
# ╠═fe8cc5dc-4e47-46af-b18c-a56812338988
# ╟─47822a88-4a89-4466-870f-207f825fb82f
# ╠═e1355c7b-9d8a-4908-abc5-f0a5d0e8c1f4
# ╟─35bed614-5c58-4690-9341-276228561436
# ╟─a55e40f6-dda7-4705-a37b-191ab7b1f16c
# ╠═ad03bfce-ca51-4727-934b-dbfdc42e9443
# ╠═a86346a2-0d31-4c25-88fb-b2ed3c3bd425
# ╠═de29ebf9-01fb-451e-8100-cc171c7b4843
