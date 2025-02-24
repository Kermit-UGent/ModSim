### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ e8dd4f61-4f51-4255-9193-62ea1368240b
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 9b9d35c4-9136-11ef-284a-d5964e231d9e
using Markdown

# ╔═╡ ae107329-e87d-4ff7-a068-362cb03fe815
using InteractiveUtils

# ╔═╡ c7c12802-f3e7-4e42-be3f-f91eabe432d8
using Catalyst

# ╔═╡ 92c50b0a-d9f1-45a9-9c8a-ced71a6ee967
using DifferentialEquations, Plots

# ╔═╡ 2a129f80-f4c2-4429-85e8-98d946fb3f36
using PlutoUI; TableOfContents()

# ╔═╡ 3a1a01f6-7f6e-411d-958b-0dce35b791ff
md"""
# Exercise - foxes and rabbits
## An ODE and a discrete (jump) problem
"""

# ╔═╡ 92a0efee-1e40-4080-aee2-831d98312176
md"""
Rabbits live on some secluded territory. Their maximum growth rate coefficient is $r$ $[year^{-1}]$, and their population capacity is $R_m$ $[\#rabbits]$. The rabbits die of old age or sickness with a dying rate coefficient $d$ $[year^{-1}]$.

At $t=0$, foxes intrude the territory and stay there. The foxes exclusively feed themselves with the rabbits. They hunt the rabbits at a rate proportional to the number of foxes (proportionality factor is $h$ $[year^{-1} \#foxes^{-1}]$). The population of foxes grows at a rate proportional to the number of rabbits (proportionality factor is $g$ $[year^{-1} \#rabbits^{-1}]$). The foxes die of old age or sickness with a dying rate coefficient $\delta$ $[year^{-1}]$.

The initial number of rabbits on the territory is 89, the initial number of
foxes intruding the territory is 2.
"""

# ╔═╡ 91e77fb0-cab3-4da7-8bc9-6f657df2497c
md"""
**Exercises:**

1. Make simulations of the evolution of rabbits and foxes as a ODE problem in the time interval $[0, 10]\;years$.

2. Make simulations of the evolution of rabbits and foxes as a discrete (jump) problem in the time interval $[0, 10]\;years$.

Assume the following parameter values: $r = 18.4\;year^{-1}$, $R_m = 120\;(\#rabbits)$, $d = 2.0\;year^{-1}$, $h = 1.4\;year^{-1}(\#foxes)^{-1}$, $g = 0.05\;year^{-1}(\#rabbits)^{-1}$ and $\delta = 1.0\;year^{-1}$.
"""

# ╔═╡ d11b4e97-df16-4641-9e32-c7f2a098ffb7
md"""
Create a *reaction network object* model for the aforementioned problem. Name it `foxes_rabbits_rn`.

Hitns:
- Use the variable names `R` and `F` for the rabbits and foxes respectively.
- Use the variable names `r`, `Rm`, `d`, `h`, `g` and `δ` for the parameters. 
"""

# ╔═╡ 6d49c334-76e0-4210-890b-e3e79097222d
# Uncomment and complete the instruction
# foxes_rabbits_rn = @reaction_network begin
#     @species missing
#     @parameters missing
# 	  missing              # natural population growth of the rabbits
#     missing              # deaths by age or sickness of the rabbits
#     missing              # hunting of rabbits by the foxes
#     missing              # gaining of foxes by hinting rabbits
#     missing              # deaths by age or sickness of the foxes
# end
foxes_rabbits_rn = @reaction_network begin
    @species R(t)=89 F(t)=2
    @parameters r=18.4 Rm=120 d=2.0 h=1.4 g=0.05 δ=1.0
	r*(1 - R/Rm), R --> 2R     # natural population growth of the rabbits
    d, R --> 0                 # deaths by age or sickness of the rabbits
	h*F, R --> 0               # hunting of rabbits by the foxes
	g*R, F --> 2F              # gaining of foxes by hinting rabbits
	δ, F --> 0                 # deaths by age or sickness of the foxes
end

# ╔═╡ 1c2c38f9-6c76-4274-8e87-178ca7813706
md"""
Convert the system to a symbolic differential equation model and verify, by analyzing the differential equation, that your model makes sense.
"""

# ╔═╡ 3ffbec86-a33a-4a71-9c8c-9f895a3a1f26
# osys = missing         # Uncomment and complete the instruction
osys = convert(ODESystem, foxes_rabbits_rn)

# ╔═╡ 577ae94c-543c-4717-981b-ce90c02d1828
md"""
Initialize a vector `u0` with the initial conditions:
"""

# ╔═╡ 55c208e5-f123-4c4e-8b03-c3e3d3eda605
# u0 = missing           # Uncomment and complete the instruction
u0 = [:R => 89, :F => 2]

# ╔═╡ 5270da4a-7d14-44cd-8672-9e481888704a
md"""
Set the timespan for the simulation:
"""

# ╔═╡ 10f878d8-d2da-4701-b190-7a2feaa0c009
# tspan = missing        # Uncomment and complete the instruction
tspan = (0.0, 10.0)

# ╔═╡ 1a835253-3677-4698-9347-d59937f9e7c0
md"""
Initialize a vector `params` with the parameter values:
"""

# ╔═╡ 2bd4dd93-1be0-427b-b372-88cb9849a124
# params = missing       # Uncomment and complete the instruction
params = [:r=>18.4, :Rm=>120, :d=>2.0, :h=>1.4, :g=>0.05, :δ=>1.0]

# ╔═╡ 4953f436-afed-4363-b840-956440fddc25
md"""
### Exercise 1 - Solve the problem as an ODE problem.
"""

# ╔═╡ a52b62fe-dd92-49ca-b09e-fadb29cc1f35
md"""
Create the ODE problem and store it in `oprob`:
"""

# ╔═╡ 629126b3-0b6b-43f0-9142-2d4d7d6270d9
# oprob = missing        # Uncomment and complete the instruction
oprob = ODEProblem(foxes_rabbits_rn, u0, tspan, params)

# ╔═╡ 3c071a33-5092-42f7-80b2-63066bdf3f22
md"""
Solve the ODE problem. Use `Tsit5()` and `saveat=0.05`. Store the solution in `osol`:
"""

# ╔═╡ 0f17c44c-a2e6-4b4a-be9c-e2e9763a3374
# osol = missing          # Uncomment and complete the instruction
osol = solve(oprob, Tsit5(), saveat=0.05)

# ╔═╡ e4bf9eca-ffee-4dca-8318-74875c636245
md"""
Plot the solution.
"""

# ╔═╡ a35eccae-7890-41a1-82a0-6cf98488d8a7
# missing               # Uncomment and complete the instruction
plot(osol)

# ╔═╡ c9ca9a36-784a-4560-93bf-4381812f41ac
md"""
### Exercise 2 - Solve the problem as a Discrete (jump) problem.
"""

# ╔═╡ 7d6c894e-90f2-44ad-9d2e-adf5efd4213e
md"""
Create a DiscreteProblem and store it in `dprob`:
"""

# ╔═╡ 581a5b6e-04a7-4fd1-9f69-a9cf8655a06c
# dprob = missing                 # Uncomment and complete the instruction
dprob = DiscreteProblem(foxes_rabbits_rn, u0, tspan, params)

# ╔═╡ 497fe4e8-b451-4dd7-b961-6e7968c6cb3d
md"""
Create a JumpProblem and store it in `jdprob`. Use the simulation method `Direct()`.
"""

# ╔═╡ 74e91495-998e-48a0-b14a-9d7be2c6410d
# jdprob = missing                # Uncomment and complete the instruction
jdprob = JumpProblem(foxes_rabbits_rn, dprob, Direct())

# ╔═╡ 9c9f097b-2a84-4383-843c-0cdc1b40e061
md"""
Solve the problem and store it in `jdsol`. Use the `SSAStepper()` stepping algorithm.
"""

# ╔═╡ c80337e8-6dbe-4c78-a793-d3a1e968c319
# jdsol = missing                # Uncomment and complete the instruction
jdsol = solve(jdprob, SSAStepper())

# ╔═╡ a743af00-cce6-4e6e-b426-107ad04960d8
md"""
Plot the solution.
"""

# ╔═╡ 46610efe-927a-4c5f-bd6a-fa8001949aa7
# missing               # Uncomment and complete the instruction
plot(jdsol)

# ╔═╡ 969bc08d-0709-4941-8fa8-3be707c5b0ed
md"""
Solve the problem several times by running the cell which solves the problem and see what happens in the plot.
"""

# ╔═╡ a039ee1b-3232-459c-885b-4166b973ab1b
md"""
Think of your class of probability theory, why does the SSA model alsways lead to extinction?
"""

# ╔═╡ 125809f6-46e4-49ec-8f96-f193d3036df6
md"- Answer: missing"
#=
Because R=0 is an absorbing state from which one cannot return.
=#

# ╔═╡ db421cb6-ba09-4476-9021-846a145bea5c
md"""
Write a piece of code that solves the problem a $1000$ times and stores the time values at which the rabbits die out.
"""

# ╔═╡ a6bc48a8-3bbc-4abc-918b-da70eab39017
# Uncomment and complete the instruction
# begin
# 	times = []                      # make empty vector
# 	# missing
# 	# ...
# end
begin
	times = []                      # make empty vector
	while length(times) < 1000      # while statement
		jdsol2 = solve(jdprob, SSAStepper())    # solve the problem
		j = findfirst(jdsol2[:R] .== 0)         # find index of first 0
		if j != nothing                         # if index is a valid index
			append!(times, jdsol2.t[j])         # append time to vector times
		end
	end
end

# ╔═╡ 2eb41790-23ef-43fc-8b69-7dd9dbb9f820
md"""
Make a histogram so that you can have an idea of the distribution when the rabbits die out. Use `bins=range(0, 10, length=121)`.
"""

# ╔═╡ 1c8188de-e8b8-44ae-b0b6-7d490f3a15e8
# missing              # Uncomment and complete the instruction
# histogram(times, bins=range(0, 10, length=121))
histogram(times, bins=range(0, 10, length=121), normalize=:pdf)

# ╔═╡ Cell order:
# ╠═9b9d35c4-9136-11ef-284a-d5964e231d9e
# ╠═ae107329-e87d-4ff7-a068-362cb03fe815
# ╠═e8dd4f61-4f51-4255-9193-62ea1368240b
# ╠═c7c12802-f3e7-4e42-be3f-f91eabe432d8
# ╠═92c50b0a-d9f1-45a9-9c8a-ced71a6ee967
# ╠═2a129f80-f4c2-4429-85e8-98d946fb3f36
# ╟─3a1a01f6-7f6e-411d-958b-0dce35b791ff
# ╟─92a0efee-1e40-4080-aee2-831d98312176
# ╟─91e77fb0-cab3-4da7-8bc9-6f657df2497c
# ╟─d11b4e97-df16-4641-9e32-c7f2a098ffb7
# ╠═6d49c334-76e0-4210-890b-e3e79097222d
# ╟─1c2c38f9-6c76-4274-8e87-178ca7813706
# ╠═3ffbec86-a33a-4a71-9c8c-9f895a3a1f26
# ╟─577ae94c-543c-4717-981b-ce90c02d1828
# ╠═55c208e5-f123-4c4e-8b03-c3e3d3eda605
# ╟─5270da4a-7d14-44cd-8672-9e481888704a
# ╠═10f878d8-d2da-4701-b190-7a2feaa0c009
# ╟─1a835253-3677-4698-9347-d59937f9e7c0
# ╠═2bd4dd93-1be0-427b-b372-88cb9849a124
# ╟─4953f436-afed-4363-b840-956440fddc25
# ╟─a52b62fe-dd92-49ca-b09e-fadb29cc1f35
# ╠═629126b3-0b6b-43f0-9142-2d4d7d6270d9
# ╟─3c071a33-5092-42f7-80b2-63066bdf3f22
# ╠═0f17c44c-a2e6-4b4a-be9c-e2e9763a3374
# ╟─e4bf9eca-ffee-4dca-8318-74875c636245
# ╠═a35eccae-7890-41a1-82a0-6cf98488d8a7
# ╟─c9ca9a36-784a-4560-93bf-4381812f41ac
# ╟─7d6c894e-90f2-44ad-9d2e-adf5efd4213e
# ╠═581a5b6e-04a7-4fd1-9f69-a9cf8655a06c
# ╟─497fe4e8-b451-4dd7-b961-6e7968c6cb3d
# ╠═74e91495-998e-48a0-b14a-9d7be2c6410d
# ╟─9c9f097b-2a84-4383-843c-0cdc1b40e061
# ╠═c80337e8-6dbe-4c78-a793-d3a1e968c319
# ╟─a743af00-cce6-4e6e-b426-107ad04960d8
# ╠═46610efe-927a-4c5f-bd6a-fa8001949aa7
# ╟─969bc08d-0709-4941-8fa8-3be707c5b0ed
# ╟─a039ee1b-3232-459c-885b-4166b973ab1b
# ╠═125809f6-46e4-49ec-8f96-f193d3036df6
# ╟─db421cb6-ba09-4476-9021-846a145bea5c
# ╠═a6bc48a8-3bbc-4abc-918b-da70eab39017
# ╟─2eb41790-23ef-43fc-8b69-7dd9dbb9f820
# ╠═1c8188de-e8b8-44ae-b0b6-7d490f3a15e8
