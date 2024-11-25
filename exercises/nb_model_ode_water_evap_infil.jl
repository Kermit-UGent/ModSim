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

# ╔═╡ 2574879f-28d1-4d30-a1aa-a637dd1b216a
using Catalyst

# ╔═╡ 6c4b3d09-09c2-4439-9167-63b59b078104
using DifferentialEquations, Plots

# ╔═╡ f55e0ca3-07ff-43a1-907f-9ad640227822
md"
### Exercise: Water evaporation and infiltration

Consider a water reservoir, such as a lake, where the water in the reservoir is in contact with the air as well as with the groundwater. We will denote the water level in the reservoir as $W$ and the groundwater level as $G$.

The water in the reservoir evaporates at a rate $k_1$ (i.e. the *evaporation coefficient*) and there can be infiltration into or from the groundwater at a rate $k_2$ (i.e., *infiltration coefficient*) depending on the magnitude of the water level in the reservoir and groundwater.

There is a natural constant inflow of water into the reservoir at a rate $I$. At time $t=0\;s$ a pumping device is switched on such that the reservoir is rapidly being emptied at an outflow rate $O$ until the level of the water reservoir drops to zero. From then on, the pump is switched off.

The system of differential equations that models the water level in a reservoir ($W$) and the groundwater level ($G$) considering evaporation, infiltration, inlet flow and outlet flow can be written down as:

$$\begin{align}
\frac{dW}{dt} & = I - O - k_1 \cdot W - k_2 \cdot (W - G) \\
\frac{dG}{dt} & = k_2 \cdot (W - G)
\end{align}$$

"

# ╔═╡ dc9d4427-0a9e-41ab-9074-3e4f10dbae7f
md"
Model the aforementioned system of differential equations using a *reaction network object*. Name it `water_evap_infil`.
"

# ╔═╡ 3c275181-503e-4806-b0e6-1731ae881a30
# Uncomment and complete the instruction
# water_evap_infil = @reaction_network begin
#     missing
# end
water_evap_infil = @reaction_network begin
    I, ∅ --> W
    -O, ∅ --> W
    k₁, W --> ∅
    k₂, W --> G
    k₂, G --> W
end

# ╔═╡ 608c61bc-c029-433e-ad8f-b13cfb40bc3d
md"
Convert the system to a symbolic differential equation model and verify that you get the same system of differential equations as given in the problem.
"

# ╔═╡ 5647d122-5e8d-4ff9-a798-4076ee93b771
# osys = missing         # Uncomment and complete the instruction
osys = convert(ODESystem, water_evap_infil)

# ╔═╡ acab2bf0-b792-4ccc-bee0-7611bedab23c
md"
Both water levels are initially $6.75\;m$. The inflow rate is constant and is $2.7\; m^3/min$. The evaporation and infiltration coefficient are $0.4\;min^{-1}$ and $1.0\;min^{-1}$ respectively. The outflow rate due to the pump is $20\;m^3/min$ and the pump stops working when $W$ equals zero. We wish to simulate the evolution of $W$ and $G$ during $20\;min$.
"

# ╔═╡ 44cf2ac6-75e2-4440-94f1-8fe6887ee1e0
md"
Initialize a vector `u₀` with the initial conditions:
"

# ╔═╡ 733bdb56-fb4f-4bc6-b50c-e3245fd59730
# u0 = missing
u0 = [:W => 6.75, :G => 6.75]

# ╔═╡ 521fcca5-3c88-463b-9395-e5871b9fc5a3
md"
Set the timespan for the simulation:
"

# ╔═╡ f2267cc9-4bfd-44af-984c-cf54aa855f91
# tspan =  missing      # Uncomment and complete the instruction
tspan = (0.0, 20)

# ╔═╡ e12840ed-d448-4f7c-893d-cfe974f62f9a
md"
Initialize a vector `param` with the parameter values:
"

# ╔═╡ 5a3d520d-fd4f-485c-a8b0-63858eca4bfc
# params = missing
params = [:I => 2.7, :O => 20.0, :k₁ => 0.4, :k₂ => 1.0]

# ╔═╡ 90c3b104-7a76-4986-9345-a544de0450d0
md"
Create the ODE problem and store it in `oprob`:
"

# ╔═╡ dec48829-0f6e-48c7-bb96-ee285065e6f9
# oprob = missing
oprob = ODEProblem(water_evap_infil, u0, tspan, params)

# ╔═╡ f31f4636-3077-4454-94c5-67b19850712e
md"
Check the order of the state variables and parameters:
"

# ╔═╡ 48301c71-e9da-427d-b0de-d184a7199193
# missing           # Uncomment and complete the instruction
species(water_evap_infil)

# ╔═╡ d895157c-5ee4-4dab-b621-a00d713d6c8e
# missing           # Uncomment and complete the instruction
parameters(water_evap_infil)

# ╔═╡ 666c3aa1-26ef-4d77-bb81-f8830e66eea2
md"
Set-up a the *condition function*, named `condition`, that contains the state variable that will hit zero.
"

# ╔═╡ 575bcb78-adf4-418d-b674-acd4a6eb7b5a
# function condition(u, t, integrator)
#   missing
# end
function condition(u, t, integrator)
  u[1]
end

# ╔═╡ d9d9ced3-27e2-4f93-a71e-4f6c8a5d1a1d
md"
Create the function called `affect!`, that will be called by the solver when the condition stored in `condition` hits zero in order to alter the value of $O$:
"

# ╔═╡ 5ab50b97-ed1d-4403-8665-8af53e33f916
# function affect!(integrator)
#   missing
# end
function affect!(integrator)
  integrator.ps[:O] = 0.0
end

# ╔═╡ daea3af6-1ce8-4550-9ee0-609ea990488e
md"
Combine both `condition` and `affect!` with the function `ContinuousCallback` in order to create the callback function, name it `cb`:
"

# ╔═╡ 40512cf6-f9a5-4d9a-b475-dd550a9679bf
# cb = missing
cb = ContinuousCallback(condition, affect!)

# ╔═╡ 9b6d4927-7502-4bdd-b711-d5145eb726d0
md"
Solve the ODE problem using the callback function `cb`. Name the solution `osol`:
"

# ╔═╡ c1c91323-a844-45be-b8b9-ec1e15c484dd
# osol = missing
osol = solve(deepcopy(oprob), Tsit5(), saveat=0.1, callback=cb)

# ╔═╡ d01e7bab-4aa9-485c-9518-8f4d814491b7
md"
Plot the results:
"

# ╔═╡ c6174712-65b4-4102-927f-9bf7a1dcd523
# missing
plot(osol)

# ╔═╡ 0600244d-139f-4f37-b952-267d068dceb6
md"
Interpret the results. Ask yourself the following questions:

1. Can you clearly see the drop in $W$? To what value does $W$ drops?
- Answer: missing
2. Why does $G$ also drop when $W$ drops? Explain.
- Answer: missing
3. To what values are $W$ and $G$ tending to go? Was the system with the initial values for $W$ and $G$ and no outflow in equilibrium? Explain.
- Answer: missing
"

# ╔═╡ Cell order:
# ╠═62dc7706-f58a-11ee-2d3d-f78f7ceca914
# ╠═4071647d-3084-4c8b-9fb7-eca7255253a9
# ╠═e079823b-8b40-42a2-a63f-1645a97b33f0
# ╠═f55e0ca3-07ff-43a1-907f-9ad640227822
# ╠═2574879f-28d1-4d30-a1aa-a637dd1b216a
# ╠═dc9d4427-0a9e-41ab-9074-3e4f10dbae7f
# ╠═3c275181-503e-4806-b0e6-1731ae881a30
# ╠═608c61bc-c029-433e-ad8f-b13cfb40bc3d
# ╠═5647d122-5e8d-4ff9-a798-4076ee93b771
# ╠═6c4b3d09-09c2-4439-9167-63b59b078104
# ╠═acab2bf0-b792-4ccc-bee0-7611bedab23c
# ╠═44cf2ac6-75e2-4440-94f1-8fe6887ee1e0
# ╠═733bdb56-fb4f-4bc6-b50c-e3245fd59730
# ╠═521fcca5-3c88-463b-9395-e5871b9fc5a3
# ╠═f2267cc9-4bfd-44af-984c-cf54aa855f91
# ╠═e12840ed-d448-4f7c-893d-cfe974f62f9a
# ╠═5a3d520d-fd4f-485c-a8b0-63858eca4bfc
# ╠═90c3b104-7a76-4986-9345-a544de0450d0
# ╠═dec48829-0f6e-48c7-bb96-ee285065e6f9
# ╠═f31f4636-3077-4454-94c5-67b19850712e
# ╠═48301c71-e9da-427d-b0de-d184a7199193
# ╠═d895157c-5ee4-4dab-b621-a00d713d6c8e
# ╠═666c3aa1-26ef-4d77-bb81-f8830e66eea2
# ╠═575bcb78-adf4-418d-b674-acd4a6eb7b5a
# ╠═d9d9ced3-27e2-4f93-a71e-4f6c8a5d1a1d
# ╠═5ab50b97-ed1d-4403-8665-8af53e33f916
# ╠═daea3af6-1ce8-4550-9ee0-609ea990488e
# ╠═40512cf6-f9a5-4d9a-b475-dd550a9679bf
# ╠═9b6d4927-7502-4bdd-b711-d5145eb726d0
# ╠═c1c91323-a844-45be-b8b9-ec1e15c484dd
# ╠═d01e7bab-4aa9-485c-9518-8f4d814491b7
# ╠═c6174712-65b4-4102-927f-9bf7a1dcd523
# ╠═0600244d-139f-4f37-b952-267d068dceb6
