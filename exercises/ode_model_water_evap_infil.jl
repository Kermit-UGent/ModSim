### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ e079823b-8b40-42a2-a63f-1645a97b33f0
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 62dc7706-f58a-11ee-2d3d-f78f7ceca914
using Markdown

# ╔═╡ 4071647d-3084-4c8b-9fb7-eca7255253a9
using InteractiveUtils

# ╔═╡ 65571bb8-e260-4a82-b0d2-198e47c56271
using PlutoUI

# ╔═╡ 2574879f-28d1-4d30-a1aa-a637dd1b216a
using Catalyst

# ╔═╡ 6c4b3d09-09c2-4439-9167-63b59b078104
using DifferentialEquations, Plots

# ╔═╡ 66e8a12c-74b6-4077-b90e-3d85e5a61d6e
hint(text) = Markdown.MD(Markdown.Admonition("hint", "Hint", [text]));

# ╔═╡ f55e0ca3-07ff-43a1-907f-9ad640227822
md"""
# Exercise: Water evaporation and infiltration

Consider a water reservoir, such as a lake, where the water in the reservoir is in contact with the air as well as with the groundwater. We will denote the water level in the reservoir as $W$ and the groundwater level as $G$.

The water in the reservoir evaporates at a rate $k_1$ (i.e. the *evaporation coefficient*) and there can be infiltration into or from the groundwater at a rate $k_2$ (i.e., *infiltration coefficient*) depending on the difference in the water level in the reservoir and groundwater (cf. $(W-G)$)

There is a natural constant inflow of water into the reservoir at a rate $I$. At time $t=0\;s$ a pumping device is switched on such that the reservoir is rapidly being emptied at an outflow rate $O$ until the level of the water reservoir drops to zero. From then on, the pump is switched off.
"""

# ╔═╡ a551f3c5-2fc3-4236-ab4a-3d9c14afc62b
md"""
!!! question
Set-up a system of differential equations modelling the above problem.
"""

# ╔═╡ 4eb95688-ba9d-4524-a522-3b8343f4e2be
hint(
md"""
The system of differential equations that models the water level in a reservoir ($W$) and the groundwater level ($G$) considering evaporation, infiltration, inlet flow and outlet flow can be written down as:

$$\begin{align}
\frac{dW}{dt} & = I - O - k_1 \cdot W - k_2 \cdot (W - G) \\
\frac{dG}{dt} & = k_2 \cdot (W - G)
\end{align}$$
"""
)

# ╔═╡ dc9d4427-0a9e-41ab-9074-3e4f10dbae7f
md"""
Model the aforementioned system of differential equations using a *reaction network object*. Name it `water_evap_infil`.
"""

# ╔═╡ 3c275181-503e-4806-b0e6-1731ae881a30
# Uncomment and complete the instruction
# water_evap_infil = @reaction_network begin
#     missing
# end

# ╔═╡ 608c61bc-c029-433e-ad8f-b13cfb40bc3d
md"""
Convert the system to a symbolic differential equation model and verify that you get the same system of differential equations as given in the problem.
"""

# ╔═╡ 5647d122-5e8d-4ff9-a798-4076ee93b771
# osys = missing         # Uncomment and complete the instruction

# ╔═╡ acab2bf0-b792-4ccc-bee0-7611bedab23c
md"""
Both water levels are initially $6.75\;m$. The inflow rate is constant and is $2.7\; m/min$. The evaporation and infiltration coefficient are $0.4\;min^{-1}$ and $1.0\;min^{-1}$ respectively. The outflow rate due to the pump is $20\;m/min$ and the pump stops working when $W$ equals zero. We wish to simulate the evolution of $W$ and $G$ during $20\;min$.
"""

# ╔═╡ 44cf2ac6-75e2-4440-94f1-8fe6887ee1e0
md"""
Initialize a vector `u0` with the initial conditions:
"""

# ╔═╡ 733bdb56-fb4f-4bc6-b50c-e3245fd59730
# u0 = missing          # Uncomment and complete the instruction

# ╔═╡ 521fcca5-3c88-463b-9395-e5871b9fc5a3
md"""
Set the timespan for the simulation:
"""

# ╔═╡ f2267cc9-4bfd-44af-984c-cf54aa855f91
# tspan =  missing      # Uncomment and complete the instruction

# ╔═╡ e12840ed-d448-4f7c-893d-cfe974f62f9a
md"""
Initialize a vector `param` with the parameter values:
"""

# ╔═╡ 5a3d520d-fd4f-485c-a8b0-63858eca4bfc
# params = missing      # Uncomment and complete the instruction

# ╔═╡ 666c3aa1-26ef-4d77-bb81-f8830e66eea2
md"""
Set-up a the *condition*, name it `condition`.
"""

# ╔═╡ ce612639-0791-4df9-bbd1-11da5ae8b247
# condition = missing      # Uncomment and complete the instruction

# ╔═╡ b4e50fbc-6780-4559-8cab-d7f2fd533eba
md"""
Make a new *reaction system* where the discrete event is included. Name it `water_evap_infil_c`.
"""

# ╔═╡ 82db28cb-d842-442f-a566-32c6fe3acc90
# @named water_evap_infil_c = missing    # Uncomment and complete the instruction

# ╔═╡ 7b649517-ccab-4f30-b75d-527295cd24a2
md"""
Complete the new *reaction system*. Name it `water_evap_infil_c_com`.
"""

# ╔═╡ 2d72257d-cf23-4a07-b659-4b886abe5abc
# water_evap_infil_c_com = missing     # Uncomment and complete the instruction

# ╔═╡ aafb49a4-b468-49b2-838d-5ddfcc852d48
md"""
Create the ODE problem and store it in `oprob`:
"""

# ╔═╡ a25d3652-13f5-47ef-9f16-c6698547a734
# oprob = missing           # Uncomment and complete the instruction

# ╔═╡ 451a3c66-5bcc-4161-af81-f89af33b5862
md"""
Solve the ODE problem. Make a deepcopy and use `Tsit5()` and `saveat=0.1`. Store the solution in `osol`:
"""

# ╔═╡ b9d2c6cb-88f0-4a88-9e61-eecc905ff3e6
# osol = missing                # Uncomment and complete the instruction

# ╔═╡ ab77c284-d379-47d1-bd86-88fa77749165
md"""
Plot the results:
"""

# ╔═╡ e29a3294-245b-445d-bb1e-12cafb2ec175
# missing                  # Uncomment and complete the instruction

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

# ╔═╡ Cell order:
# ╠═62dc7706-f58a-11ee-2d3d-f78f7ceca914
# ╠═4071647d-3084-4c8b-9fb7-eca7255253a9
# ╠═e079823b-8b40-42a2-a63f-1645a97b33f0
# ╠═65571bb8-e260-4a82-b0d2-198e47c56271
# ╟─66e8a12c-74b6-4077-b90e-3d85e5a61d6e
# ╟─f55e0ca3-07ff-43a1-907f-9ad640227822
# ╟─a551f3c5-2fc3-4236-ab4a-3d9c14afc62b
# ╟─4eb95688-ba9d-4524-a522-3b8343f4e2be
# ╠═2574879f-28d1-4d30-a1aa-a637dd1b216a
# ╟─dc9d4427-0a9e-41ab-9074-3e4f10dbae7f
# ╠═3c275181-503e-4806-b0e6-1731ae881a30
# ╟─608c61bc-c029-433e-ad8f-b13cfb40bc3d
# ╠═5647d122-5e8d-4ff9-a798-4076ee93b771
# ╠═6c4b3d09-09c2-4439-9167-63b59b078104
# ╟─acab2bf0-b792-4ccc-bee0-7611bedab23c
# ╟─44cf2ac6-75e2-4440-94f1-8fe6887ee1e0
# ╠═733bdb56-fb4f-4bc6-b50c-e3245fd59730
# ╟─521fcca5-3c88-463b-9395-e5871b9fc5a3
# ╠═f2267cc9-4bfd-44af-984c-cf54aa855f91
# ╟─e12840ed-d448-4f7c-893d-cfe974f62f9a
# ╠═5a3d520d-fd4f-485c-a8b0-63858eca4bfc
# ╟─666c3aa1-26ef-4d77-bb81-f8830e66eea2
# ╠═ce612639-0791-4df9-bbd1-11da5ae8b247
# ╟─b4e50fbc-6780-4559-8cab-d7f2fd533eba
# ╠═82db28cb-d842-442f-a566-32c6fe3acc90
# ╟─7b649517-ccab-4f30-b75d-527295cd24a2
# ╠═2d72257d-cf23-4a07-b659-4b886abe5abc
# ╟─aafb49a4-b468-49b2-838d-5ddfcc852d48
# ╠═a25d3652-13f5-47ef-9f16-c6698547a734
# ╟─451a3c66-5bcc-4161-af81-f89af33b5862
# ╠═b9d2c6cb-88f0-4a88-9e61-eecc905ff3e6
# ╟─ab77c284-d379-47d1-bd86-88fa77749165
# ╠═e29a3294-245b-445d-bb1e-12cafb2ec175
# ╟─66588291-d399-4169-9383-ac6c05cdf906
# ╟─89ece6c0-3690-4265-8b1d-c3a3cc8b095f
# ╟─92fa2a5a-123b-434c-a97f-39d9727a5ab6
# ╟─71dde4bd-711c-40c3-8ee2-9f7a79d72aa2
# ╟─d04d906d-1432-4753-81b9-b03a06df9816
# ╟─4ce3328e-3436-45ad-b899-9b901b53a8ea
