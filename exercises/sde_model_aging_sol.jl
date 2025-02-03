### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 85b6ee78-2f7b-4f83-908b-c48fade56fcb
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 0c457482-95fe-11ef-0873-65164704c7a6
using Markdown

# ╔═╡ 7a0d5c7f-fb37-4b38-9736-32e81c8014c4
using InteractiveUtils

# ╔═╡ 68588061-fedd-4f19-9e89-39b1bd791e33
using Catalyst

# ╔═╡ 1d0a8fb2-dab9-44ae-a230-f3e92cf3cd6a
using DifferentialEquations, Plots

# ╔═╡ a176adcc-bde0-4a25-aa80-33379e05e63f
md"""
# Exercise: Aging and saturated repair
"""

# ╔═╡ 818ba2af-83dd-4ddb-93e8-d0549ddd806b
md"""
Aging is ultimately correlated with damaged cells. These damaged cells are called **senescent cells**. Senescent cells are cells that eventually stop multiplying but don't die off when they should. They instead remain and secrete factors that cause **chronic inflammation** and **reduce regeneration, leading to disease and decline**. Let $X$ denote the number of senescent cells - or the **damage** - in a human body. Research shows that they are **produced** at a **rate proportional to age**. Fortunately, in living organisms, these senescent cells are removed by so-called **natural-killer cells**. However, like many biological processen, this biological process of removing senescent cells is **saturated**.

Hence, the model that we could adopt in order to predict the number of senescent cells (or damaged cells) $X$ in a human body, has two features:

1. **production of damage that rises linearly with age**, and
2. the **saturating removal of damage**.

A possible model is the following differential equation:

$$\cfrac{dX}{dt} = \mu t - \beta \cfrac{X}{X + \kappa}$$

Lets denote $X$ in trillions [$tn$]. The term $\mu t$ stands for the procution of senescent cells, and the term $- \beta \cfrac{X}{X + \kappa}$ for the removal of senescent cells. The time $t$ is in years [$y$]. The coefficient $\eta$ [$tn/y^2$] is a proportionality factor for the production, $\beta$ [$tn/y$] is the removal rate coefficient and $\cfrac{X}{X + \kappa}$ [$-$] is the corresponding saturation factor, with $\kappa$ [$tn$] the amount of $X$ at which they inhibit half of their own removal rate.

If this model was all there was, then all individuals would age at the same rate and die at the same age. The model does not explain why genetically identical organisms could differ in the number of senenscent cells. Therefore, we will introduce  noise in the model by treating it as a **Stochastic Differential Equation** (SDE) model, where noise will be added to both, production and removal processes.
"""

# ╔═╡ d917e75d-eefb-415f-ba45-2f1d7c62d63c
md"""
#### Implementation of the system
"""

# ╔═╡ 7bda135b-2e0a-425f-b30c-e5a09dbb8a4a
md"""
Implement the above ODE into a *reaction network object*, and name it `senescent_cells_rn`.

Take a default initial value $X(t=0)=0.0$ for the *species* $X$, and default values of $\mu=0.00558$, $\beta=0.4464$, $\kappa=1.116$ for the *parameters* in the model. In addition to the parameters, take $\eta=0.1$ as the *default noise scaling* factor, and, furthermore, set the noise scaling to $0.5$ for the process exhibiting the saturating removal of damage.
"""

# ╔═╡ 38edfb50-7d6f-4fcc-b328-95ecb26d6de1
# Uncomment and complete the instruction
# senescent_cells_rn = @reaction_network begin
#     @species missing
#     @parameters missing
#     @default_noise_scaling missing
#     missing
#     missing
# end
senescent_cells_rn = @reaction_network begin
    @species X(t)=0.0
    @parameters μ=0.00558 β=0.4464 κ=1.116 η=0.1
    @default_noise_scaling η
    μ*t, 0 --> X
    # β/(κ+X), X --> 0, [noise_scaling = 0.5]   # alternative
    mm(X, β, κ), X => 0, [noise_scaling = 0.5]
end

# ╔═╡ 51edb519-2086-45a5-a666-b5a58c51b5e8
md"""
Convert this *reaction model* into a symbolic differential equation model and verify that you get the correct differential equation as mentioned above.
"""

# ╔═╡ b3dd11e4-1209-46f3-a47a-07b54c5e9d77
# osys = missing                # Uncomment and complete the instruction
osys = convert(ODESystem, senescent_cells_rn)

# ╔═╡ 0524da26-35cb-4846-b7ef-f5842ae5db74
md"""
#### Setting initial condition, time span and parameters.
"""

# ╔═╡ 98c1f94e-aea5-40c5-bee7-8df8794fedb3
md"""
Initialize a vector `u0` with the default initial condition, set the timespan for the simulation (we will simulate from $0\;y$ to $120\;y$), and initialize a vector `param` with the default parameter values. In that way, later, you can change the initial condition and the parameter values if you want to try other values.
"""

# ╔═╡ 80db134e-bc57-41b0-8e1b-2ac6bff0c806
# u0 = missing                  # Uncomment and complete the instruction
u0 = [:X => 0.0]

# ╔═╡ ab8bef21-8c4f-4962-b9cb-e87ac1db049a
# tspan = missing               # Uncomment and complete the instruction
tspan = (0.0, 120.0)

# ╔═╡ ab302108-e9a7-4941-a762-c0dc109c9b1d
# params = missing              # Uncomment and complete the instruction
params = [:μ=>0.00558, :β=>0.4464, :κ=>1.116]

# ╔═╡ 3934ba8f-49c9-4f0c-8373-c8010f397d66
md"""
#### Simulating the system as an SDE problem
"""

# ╔═╡ 45fceec8-a641-43ed-9a17-72b795e09297
md"""
Create the SDE problem.
"""

# ╔═╡ 7d164cdf-fd63-4c85-b168-3279d6d658eb
# sprob = missing                # Uncomment and complete the instruction
sprob = SDEProblem(senescent_cells_rn, u0, tspan, params)

# ╔═╡ 9f9a616a-4353-4f9a-8da5-effe4d214ed0
md"""
Solve the SDE problem using `EM()`as solver and time step `dt=0.1`.
"""

# ╔═╡ e4b9393d-d33c-4e16-814e-347b35435a82
# ssol = missing                 # Uncomment and complete the instruction
ssol = solve(sprob, EM(), dt=0.1)

# ╔═╡ 91881bc3-a1f0-48f7-82e8-7d7acea2f4de
md"""
Plot the solutions. Use the option `ylim=(0, 6)` in order to limit the range of $X$.
"""

# ╔═╡ 7128ae66-c576-435b-8e33-4221f469ba1c
# missing                     # Uncomment and complete the instruction
plot(ssol, ylim=(0, 6))

# ╔═╡ 1bce56da-0c35-41df-adbd-b5dee116be11
md"""
Execute the cell, where the SDE problem is being solved, a few times and watch the (stochastic) changes in the solutions.
"""

# ╔═╡ 7ca7c177-40b3-42c0-a57d-2f7eb874abdf
md"""
#### Simulating the system as an EnsembleProblem.
"""

# ╔═╡ 33debdb3-2008-4055-a778-600bf7409361
md"""
In order to see to have an idea of the extend of the stochastic effect on the solutions, we can make a so-called *EnsembleProblem*. This allows us to plot many possible solutions in one plot.
"""

# ╔═╡ 85c7adcd-5a69-4900-849e-6222568bf68f
md"""
Create an `EnsembleProblem` based on `sprob`.
"""

# ╔═╡ 1f25730d-d93e-462a-ac02-f83381633849
# esprob = missing                # Uncomment and complete the instruction
esprob = EnsembleProblem(sprob)

# ╔═╡ 7d531f11-1a64-4b43-aa4f-04d282a615bd
md"""
Solve the ensemble problem. Use `EM()` as solver, take a time step `dt=0.1`, use the options `save_everystep=true`, and `trajectories=100`.
"""

# ╔═╡ 5aaf1395-9b7c-451f-b1a3-44e5ba3c4d6a
# essol = missing                  # Uncomment and complete the instruction
essol = solve(esprob, EM(), dt=0.1, save_everystep=true, trajectories=100)

# ╔═╡ c048528f-58a2-415d-a338-fc87302367b8
md"""
Plot the solutions. Use the option `ylim=(0, 6)` in order to limit the range of $X$.
"""

# ╔═╡ 147adea6-e7ed-4f7b-b753-127e41500109
# missing                     # Uncomment and complete the instruction
plot(essol, ylim=(0, 6))

# ╔═╡ b10a73a1-63ae-4a15-9a83-394f6ed4e36e
md"""
Interpret the results. Ask yourself the following question:

1. Suppose that $5$ trillion senescent cells is about the maximum a human body can bear. What is the (approximate) corresponding range of ages?
"""

# ╔═╡ a4f5a40e-c4ec-45ef-b44b-b64b71717c25
md"- Answer: missing"
#=
If we imagine a horizontal line through X = 5 and look at the crossings with the earliest and latest curves, then the range in age would be about [60, 120].
=#

# ╔═╡ 88c5c0f5-5bbf-4c54-9a9a-bf6369a2e895
md"""
2. What is the effect of halving the damage rate $\mu$?
"""

# ╔═╡ 369b91ee-68c0-406d-8832-363faec6feec
md"- Answer: missing"

# ╔═╡ 5ded917b-9524-4492-9b05-5a1ac1507dc2
md"""
3. What is the effect of doubling the damage removal rate $\beta$?
"""

# ╔═╡ 2560cf75-95ff-4760-a25c-baa812b40fc6
md"- Answer: missing"

# ╔═╡ 75b08cc8-20b2-4423-bea2-bf218c9f8202
md"""
4. What is the effect of halving the noise?
"""

# ╔═╡ 6123aaa2-2afb-46a3-bcab-1cf3d8335666
md"- Answer: missing"

# ╔═╡ Cell order:
# ╠═0c457482-95fe-11ef-0873-65164704c7a6
# ╠═7a0d5c7f-fb37-4b38-9736-32e81c8014c4
# ╠═85b6ee78-2f7b-4f83-908b-c48fade56fcb
# ╠═68588061-fedd-4f19-9e89-39b1bd791e33
# ╠═1d0a8fb2-dab9-44ae-a230-f3e92cf3cd6a
# ╟─a176adcc-bde0-4a25-aa80-33379e05e63f
# ╟─818ba2af-83dd-4ddb-93e8-d0549ddd806b
# ╟─d917e75d-eefb-415f-ba45-2f1d7c62d63c
# ╟─7bda135b-2e0a-425f-b30c-e5a09dbb8a4a
# ╠═38edfb50-7d6f-4fcc-b328-95ecb26d6de1
# ╟─51edb519-2086-45a5-a666-b5a58c51b5e8
# ╠═b3dd11e4-1209-46f3-a47a-07b54c5e9d77
# ╟─0524da26-35cb-4846-b7ef-f5842ae5db74
# ╟─98c1f94e-aea5-40c5-bee7-8df8794fedb3
# ╠═80db134e-bc57-41b0-8e1b-2ac6bff0c806
# ╠═ab8bef21-8c4f-4962-b9cb-e87ac1db049a
# ╠═ab302108-e9a7-4941-a762-c0dc109c9b1d
# ╟─3934ba8f-49c9-4f0c-8373-c8010f397d66
# ╟─45fceec8-a641-43ed-9a17-72b795e09297
# ╠═7d164cdf-fd63-4c85-b168-3279d6d658eb
# ╟─9f9a616a-4353-4f9a-8da5-effe4d214ed0
# ╠═e4b9393d-d33c-4e16-814e-347b35435a82
# ╟─91881bc3-a1f0-48f7-82e8-7d7acea2f4de
# ╠═7128ae66-c576-435b-8e33-4221f469ba1c
# ╟─1bce56da-0c35-41df-adbd-b5dee116be11
# ╟─7ca7c177-40b3-42c0-a57d-2f7eb874abdf
# ╟─33debdb3-2008-4055-a778-600bf7409361
# ╟─85c7adcd-5a69-4900-849e-6222568bf68f
# ╠═1f25730d-d93e-462a-ac02-f83381633849
# ╟─7d531f11-1a64-4b43-aa4f-04d282a615bd
# ╠═5aaf1395-9b7c-451f-b1a3-44e5ba3c4d6a
# ╟─c048528f-58a2-415d-a338-fc87302367b8
# ╠═147adea6-e7ed-4f7b-b753-127e41500109
# ╟─b10a73a1-63ae-4a15-9a83-394f6ed4e36e
# ╠═a4f5a40e-c4ec-45ef-b44b-b64b71717c25
# ╟─88c5c0f5-5bbf-4c54-9a9a-bf6369a2e895
# ╠═369b91ee-68c0-406d-8832-363faec6feec
# ╟─5ded917b-9524-4492-9b05-5a1ac1507dc2
# ╠═2560cf75-95ff-4760-a25c-baa812b40fc6
# ╟─75b08cc8-20b2-4423-bea2-bf218c9f8202
# ╠═6123aaa2-2afb-46a3-bcab-1cf3d8335666
