### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 806ed5ad-788b-4f4b-b72f-3720766b6959
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate("..")
end
  ╠═╡ =#

# ╔═╡ f3ae43dc-f7e1-11ee-3615-695b9e85b621
using Markdown

# ╔═╡ 28217175-2968-4219-81c9-eca0aa624196
using InteractiveUtils

# ╔═╡ 5d812731-8340-49c9-b181-aef5299388b5
using Catalyst

# ╔═╡ 63c243ff-649d-4494-98ef-e4dee5e88030
using DifferentialEquations, Plots

# ╔═╡ 2e8b2af9-2fed-49aa-a1cd-76c822664a5d
md"
### Exercise: Irrigation experiment
"

# ╔═╡ a7bb2db8-8e7e-45b6-843d-2db43f1a7ad1
md"
An irrigation experiment is carried out on a soil column consisting of two layers of soil, each with specific soil characteristics. An adjustable volume of water per unit of time, $R$, is irrigated evenly over the soil column, starting with $5\;mm\,h^{-1}$ ($mm$ indicates a volume of water: $1\;mm = 10^{-3}\,m^3$). After $60\;h$ the added flow rate is increased to $10\;mm\,h^{-1}$. The water falls on the upper layer and percolates to the lower layer. The relative moisture content in both layers (i.e., relative to their residual moisture contents) is denoted by $S_1$ and $S_2$. Initially a moisture content of $30\;mm$ is present in the upper layer (cf. $S_1$) and of $25\;mm$ in the lower layer (cf. $S_2$). The residual moisture content in the upper layer is $S_{1,res}=10 \;mm$.

A model description of the relative moisture content in both soil layers is given by:

$$\begin{align}
\frac{dS_1}{dt} &= R\left(1-\cfrac{S_{1,res}}{S_{max}}\right) - \cfrac{R}{S_{max}}S_1 - \cfrac{k}{S_{max}}S_1 \\
\frac{dS_2}{dt} &= \cfrac{k}{S_{max}}S_1 - v \,S_2^2
\end{align}$$

Here $S_{max}$ ($150\;mm$) is the saturated water quantity for the top soil layer, $k$ is the percolation ratio ($3\;mm\,h^{-1}$) and $v$ is the flow factor into the groundwater ($10^{-3}\;h^{-1}\,mm^{-1}$).

Three measurements are made over the duration ($150\;h$) of the experiment:
- The excess running water (runoff): $R \cfrac{S_1 + S_{1,res}}{S_{max}}$,
- The underground outflow into groundwater: $v\,S_2^2$,
- The amount of percolation to deeper soil layers: $\cfrac{k}{S_{max}} S_1$.
"

# ╔═╡ 072d96fa-4933-43ed-a4fd-162f64be5cdd
md"
Create a *reaction network object* model for the aforementioned problem in order to simulate the evolution of the three afore mentioned measurements during $150\;h$. Name it `irrigation_mod`.

Tips:
- You can use any kind of expression for the reaction rates.
- The term $- v \,S_2^2$ is created by the reaction `v*S₂^2, S₂ => ∅`, meaning we overrule custom behavior and say that `S₂` is removed with a rate of $vS_2^2$.
"

# ╔═╡ 33c6a48a-4424-400a-b984-7b19cb19149e
R(t) = ifelse(t < 60, 5, 10)

# ╔═╡ e557a067-52a1-42d5-98ae-40a0c63a2611
# Uncomment and complete the instruction
# irrigation_mod = @reaction_network begin
# 	...
# end
irrigation_mod = @reaction_network begin
    k/Smax, S₁ --> S₂
    v*S₂^2, S₂ => ∅
    R(t) * (1 - S₁res / Smax), ∅ --> S₁
    R(t)/Smax, S₁ --> ∅
end

# ╔═╡ 79277f44-ffca-44e3-867d-07a95dcb538c
md"
Convert the system to a symbolic differential equation model and verify, by analyzing the differential equation, that your model is correctly implemented.
"

# ╔═╡ 2f5b954b-1615-4b14-bc3b-3426c9296221
# osys = ...         # Uncomment and complete the instruction
osys = convert(ODESystem, irrigation_mod)

# ╔═╡ c163afdb-9411-4879-945a-8ee722ddb15a
md"
Initialize a vector `u₀` with the initial conditions:
"

# ╔═╡ 32dcda4d-2e05-4451-86c2-22c2ea1d9b63
# u₀ = ...           # Uncomment and complete the instruction
u₀ = [:S₁ => 30, :S₂ => 25]

# ╔═╡ f974745b-5292-4bfa-aa8a-9836ae5870c3
md"
Set the timespan for the simulation:
"

# ╔═╡ 97c84665-c6ff-4c7a-8dfb-772244aa4986
# tspan = ...        # Uncomment and complete the instruction
tspan = (0.0, 150.0)

# ╔═╡ 0ceef9b8-4647-43a4-b9c3-11959834305b
md"
Initialize a vector `param` with the parameter values:
"

# ╔═╡ d172fc24-0fbe-46b4-a030-1ba44b4b57cd
# params = ...       # Uncomment and complete the instruction
params = [:k => 3.0, :Smax => 150.0, :v => 1.0e-3, :S₁res => 10.0]

# ╔═╡ 309fc7b7-ba9d-41ea-a834-e3b60eeb0226
md"
Create the ODE problem and store it in `oprob`:
"

# ╔═╡ d69a3772-a016-4555-88fe-cb5daf86e0ac
# oprob = ...        # Uncomment and complete the instruction
oprob = ODEProblem(irrigation_mod, u₀, tspan, params)

# ╔═╡ 7793eb61-3a2a-471a-b283-20db4a059b70
md"
Create the *condition* that contains the timepoint for the sudden change in $R$. Store it in `condition`:
"

# ╔═╡ c2bb92b9-ad9e-4bbc-8914-822848aea552
# condition = ...     # Uncomment and complete the instruction
condition = [60]

# ╔═╡ 0e469b9f-f825-4f64-8340-5047598197cd
md"
Determine the index number of the relevant parameter that needs to be modified in the model:
"

# ╔═╡ d6a89017-db51-4f19-a2fc-e919d4bc8815
#  ...                 # Uncomment and complete the instruction
parameters(irrigation_mod)

# ╔═╡ 62ef97da-aedf-446d-ba2c-8c0ac72c9f05
md"
Create a function called `affect!`, that will be called by the solver at the timepoint(s) stored in `condition` in order to alter the relevant parameter value:
"

# ╔═╡ 596b411b-2e76-48af-b3ce-8b73aa6f8700
# function affect!(integrator)
#     ...
# end
function affect!(integrator)
    integrator.p[4] += 5.0      # R is the 4th parameter !!!
end

# ╔═╡ a38e9488-af43-49c6-9e7b-51f15e9a65ab
md"
Create the callback function using `condition` and `affect!`. Store it in `cb`:
"

# ╔═╡ a376c25d-059c-4a39-8c71-c8eb4310df72
cb = PresetTimeCallback(condition, affect!)

# ╔═╡ 18717738-ae99-4b41-8414-4c1823307dde
md"
Solve the ODE problem. Use `Tsit5()` and `saveat=0.5`. Store the solution in `osol`:
"

# ╔═╡ dd5ffca7-4b39-4089-afee-dd73e4a9ac15
osol = solve(deepcopy(oprob), Tsit5(), saveat=0.5, callback=cb)

# ╔═╡ 64e02433-beea-4ca5-9b87-d118540d380e
md"
Plot the three measurements (runoff, outflow, percolation) in one figure.

Tips:
- You can access the (vector) results of $S_1$ and $S_2$, by `osol[:S₁]` and `osol[:S₂]` respectively.
- The time vector is accessible through `osol.t`
- To calculate the three measurements, use element-wise operations (place a dot in front of all operators, e.g., `.+`, `./`, `.*`, `.^`)

The `begin`-`end`-statement is used to execute multiple lines of code. With the exclamation mark `!` in `plot!` you can plot in the same figure.
"

# ╔═╡ 20442788-9c2f-4b0d-85bc-c380c6c839e4
# begin
# 	plot(osol.t, ..., xaxis="time [s]", label="runoff")
# 	plot!(osol.t, ..., label="outflow")
# 	plot!(osol.t, ..., label="percolation")
# end
begin
	plot(osol.t, (osol[:S₁] .+ 10) ./ 150 .* 5, xaxis="time [s]", label="runoff")
	plot!(osol.t, (osol[:S₂].^2) .* 1.0e-3, label="outflow")
	plot!(osol.t, 3.0 .* osol[:S₁] ./ 150, label="percolation")
end

# ╔═╡ d1a9c363-3ef6-4bbb-a150-63b44f9e7cf1
md"
Interpret the results. Ask yourself the following questions:

1. Can you clearly see the effect of the increase in $R$?
- Answer: ...
2. Argue why the outflow and the percolation tend to the same value.
- Answer: ...
"

# ╔═╡ Cell order:
# ╠═f3ae43dc-f7e1-11ee-3615-695b9e85b621
# ╠═28217175-2968-4219-81c9-eca0aa624196
# ╠═806ed5ad-788b-4f4b-b72f-3720766b6959
# ╠═2e8b2af9-2fed-49aa-a1cd-76c822664a5d
# ╠═a7bb2db8-8e7e-45b6-843d-2db43f1a7ad1
# ╠═5d812731-8340-49c9-b181-aef5299388b5
# ╠═072d96fa-4933-43ed-a4fd-162f64be5cdd
# ╠═33c6a48a-4424-400a-b984-7b19cb19149e
# ╠═e557a067-52a1-42d5-98ae-40a0c63a2611
# ╠═79277f44-ffca-44e3-867d-07a95dcb538c
# ╠═2f5b954b-1615-4b14-bc3b-3426c9296221
# ╠═63c243ff-649d-4494-98ef-e4dee5e88030
# ╠═c163afdb-9411-4879-945a-8ee722ddb15a
# ╠═32dcda4d-2e05-4451-86c2-22c2ea1d9b63
# ╠═f974745b-5292-4bfa-aa8a-9836ae5870c3
# ╠═97c84665-c6ff-4c7a-8dfb-772244aa4986
# ╠═0ceef9b8-4647-43a4-b9c3-11959834305b
# ╠═d172fc24-0fbe-46b4-a030-1ba44b4b57cd
# ╠═309fc7b7-ba9d-41ea-a834-e3b60eeb0226
# ╠═d69a3772-a016-4555-88fe-cb5daf86e0ac
# ╠═7793eb61-3a2a-471a-b283-20db4a059b70
# ╠═c2bb92b9-ad9e-4bbc-8914-822848aea552
# ╠═0e469b9f-f825-4f64-8340-5047598197cd
# ╠═d6a89017-db51-4f19-a2fc-e919d4bc8815
# ╠═62ef97da-aedf-446d-ba2c-8c0ac72c9f05
# ╠═596b411b-2e76-48af-b3ce-8b73aa6f8700
# ╠═a38e9488-af43-49c6-9e7b-51f15e9a65ab
# ╠═a376c25d-059c-4a39-8c71-c8eb4310df72
# ╠═18717738-ae99-4b41-8414-4c1823307dde
# ╠═dd5ffca7-4b39-4089-afee-dd73e4a9ac15
# ╠═64e02433-beea-4ca5-9b87-d118540d380e
# ╠═20442788-9c2f-4b0d-85bc-c380c6c839e4
# ╠═d1a9c363-3ef6-4bbb-a150-63b44f9e7cf1
