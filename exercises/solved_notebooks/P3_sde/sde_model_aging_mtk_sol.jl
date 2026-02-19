### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 85b6ee78-2f7b-4f83-908b-c48fade56fcb
using Pkg; Pkg.activate("..")

# ╔═╡ 0c457482-95fe-11ef-0873-65164704c7a6
using Markdown, InteractiveUtils

# ╔═╡ 1d0a8fb2-dab9-44ae-a230-f3e92cf3cd6a
using ModelingToolkit, StochasticDiffEq

# ╔═╡ 7a9fa5bc-7bb3-436e-b796-b946c855fe0a
using ModelingToolkit: t_nounits as t, D_nounits as D

# ╔═╡ 22e0eda2-a922-48bb-9b6f-ec0c6dfeb524
using StatsPlots, PlutoUI, StatsBase; TableOfContents()

# ╔═╡ a176adcc-bde0-4a25-aa80-33379e05e63f
md"""
# Exercise: Aging and saturated repair
"""

# ╔═╡ 3750e702-22b2-4e2f-aa76-d906f5fd315e
md"""
![Aging](https://www.biotechniques.com/wp-content/uploads/2024/11/aging-800x344.png)
"""

# ╔═╡ 818ba2af-83dd-4ddb-93e8-d0549ddd806b
md"""
Aging is ultimately correlated with damaged cells. These damaged cells are called **senescent cells**. Senescent cells are cells that eventually stop multiplying but don't die off when they should. They instead remain and secrete factors that cause **chronic inflammation** and **reduce regeneration, leading to disease and decline**. Let $X$ denote the number of senescent cells - or the **damage** - in a human body. Research shows that they are **produced** at a **rate proportional to age**. Fortunately, in living organisms, these senescent cells are removed by so-called **natural-killer cells**. However, like many biological processen, this biological process of removing senescent cells is **saturated**.

Hence, the model that we could adopt in order to predict the number of senescent cells (or damaged cells) $X$ in a human body, has two features:

1. **production of damage that rises linearly with age**, and
2. the **saturating removal of damage**.

A possible model is the following differential equation:

$$\cfrac{dX}{dt} = \mu t - \beta \cfrac{X}{X + \kappa} + \cfrac{dW}{dt}$$

Lets denote the amount of senescent cells as $X$ in trillions [$tn$]. The term $\mu t$ stands for the procution of senescent cells, and the term $- \beta \cfrac{X}{X + \kappa}$ for the removal of senescent cells. The time $t$ is in years [$y$]. The coefficient $\mu$ [$tn/y^2$] is a proportionality factor for the production, $\beta$ [$tn/y$] is the removal rate coefficient and $\cfrac{X}{X + \kappa}$ [$-$] is the corresponding saturation factor, with $\kappa$ [$tn$] the amount of $X$ at which they inhibit half of their own removal rate. $W$ is a Wiener process.

If this model was all there was, then all individuals would age at the same rate and die at the same age. The model does not explain why genetically identical organisms could differ in the number of senenscent cells. Therefore, we will introduce  noise in the model by treating it as a **Stochastic Differential Equation** (SDE) model, where noise will be added to both, production and removal processes.
"""

# ╔═╡ d917e75d-eefb-415f-ba45-2f1d7c62d63c
md"""
## Implementation of the system
"""

# ╔═╡ 7bda135b-2e0a-425f-b30c-e5a09dbb8a4a
md"""
Implement the above ODE using MTK.

Take a default initial value $X(t=0)=$`0.0` for the *species* `X`, and default values of `μ=0.00558`, `β=0.4364`, `κ=1.116` for the *parameters* in the model. Use a noise scaling parameter `n=0.2`.
"""

# ╔═╡ 36c678de-56ea-437f-b667-064149075fa1
md"""
Define the only variable $X$ with `@variables`.
"""

# ╔═╡ 72d60bcc-7d78-45e8-90ca-306bb429aede
# missing
@variables X(t)=0.0

# ╔═╡ e560a876-57af-42b7-a16a-8f18676fe8e6
md"""
Define the parameters with `@parameters`.
"""

# ╔═╡ 27783ea2-c055-4ee6-a5cc-dbc5b4f516f7
# missing
@parameters μ=0.00558 β=0.4464 κ=1.116 n=0.2

# ╔═╡ e881a30a-d622-468e-9be9-69523e5c5ae6
md"""
Instantiate a source `B` of Brownian noise (a Wiener process by default) with `@brownian`.
"""

# ╔═╡ 8ad575bd-abc1-41b5-8fe3-e28a5c227ef5
# missing
@brownian B

# ╔═╡ ef9af0de-b394-4a59-81e8-569716481862
md"""
Set up the model equation with the noise `n*B` added as a term in the equation and put it in a vector/array.
"""

# ╔═╡ 858cb14b-ae3d-4b71-94b1-ff7d66a532b9
# eq_scs = missing
eq_scs = [D(X) ~ μ*t - β*X/(X+κ) + n*B]

# ╔═╡ 9dc5cfed-fa3d-4637-b1da-37fae3429565
md"""
Define the expression for the diffusion/noise term and put it in a vector/array.
"""

# ╔═╡ 705d5c4f-6a7b-4067-80ec-8fb5ccaedb18
md"""
Build a system of equations with `@mtkbuild` and the function `System`. Provide the following arguments to `System`: the model equation and `t`. Name the system `sys_scs`.
"""

# ╔═╡ 4e5b490e-67a8-40de-a49f-9f88dcb4b558
# missing
@mtkbuild sys_scs = System(eq_scs, t)

# ╔═╡ 0524da26-35cb-4846-b7ef-f5842ae5db74
md"""
## Setting initial condition, time span and parameters.
"""

# ╔═╡ 98c1f94e-aea5-40c5-bee7-8df8794fedb3
md"""
Initialize a vector `u0` with the default initial condition, set the timespan for the simulation (we will simulate from `0.0` $y$ to `120.0` $y$), and initialize a vector `parms` with the default parameter values. In that way, later, you can change the initial condition and the parameter values if you want to try other values.
"""

# ╔═╡ 80db134e-bc57-41b0-8e1b-2ac6bff0c806
# u0 = missing
u0 = [X => 0.0]

# ╔═╡ ab8bef21-8c4f-4962-b9cb-e87ac1db049a
# tspan = missing
tspan = (0.0, 120.0)

# ╔═╡ ab302108-e9a7-4941-a762-c0dc109c9b1d
# parms = missing
parms = [μ=>0.00558, β=>0.4364, κ=>1.116]

# ╔═╡ 3934ba8f-49c9-4f0c-8373-c8010f397d66
md"""
## Simulating the system as an SDE problem
"""

# ╔═╡ 45fceec8-a641-43ed-9a17-72b795e09297
md"""
Create the SDE problem.
"""

# ╔═╡ 7d164cdf-fd63-4c85-b168-3279d6d658eb
# sprob_scs = missing
sprob_scs = SDEProblem(sys_scs, u0, tspan, parms)

# ╔═╡ 9f9a616a-4353-4f9a-8da5-effe4d214ed0
md"""
Solve the SDE problem using `EM()`as solver and time step `dt=0.1`.
"""

# ╔═╡ e4b9393d-d33c-4e16-814e-347b35435a82
# ssol_scs = missing
ssol_scs = solve(sprob_scs, EM(), dt=0.1)
# ssol_scs = solve(sprob_scs, SRIW1(), dt=0.1)

# ╔═╡ 91881bc3-a1f0-48f7-82e8-7d7acea2f4de
md"""
Plot the solutions. Use the option `ylim=(0, 6)` in order to limit the range of $X$.
"""

# ╔═╡ 7128ae66-c576-435b-8e33-4221f469ba1c
# missing
plot(ssol_scs, ylim=(0, 6))

# ╔═╡ 1bce56da-0c35-41df-adbd-b5dee116be11
md"""
Execute the cell, where the SDE problem is being solved, a few times and watch the (stochastic) changes in the solutions.
"""

# ╔═╡ 7ca7c177-40b3-42c0-a57d-2f7eb874abdf
md"""
## Simulating the system as an EnsembleProblem.
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
# esprob_scs = missing
esprob_scs = EnsembleProblem(sprob_scs)

# ╔═╡ 7d531f11-1a64-4b43-aa4f-04d282a615bd
md"""
Solve the ensemble problem. Use `EM()` as solver, take a time step `dt=0.1`, use the option `trajectories=100`.
"""
#=
Solve the ensemble problem. Use `EM()` as solver, take a time step `dt=0.1`, use the options `save_everystep=true`, and `trajectories=100`.
=#

# ╔═╡ 5aaf1395-9b7c-451f-b1a3-44e5ba3c4d6a
# essol_scs = missing
essol_scs = solve(esprob_scs, EM(), dt=0.1, trajectories=100)

# ╔═╡ c048528f-58a2-415d-a338-fc87302367b8
md"""
Plot the solutions. Use the option `ylim=(0, 6)` in order to limit the range of $X$.
"""

# ╔═╡ 147adea6-e7ed-4f7b-b753-127e41500109
# missing
plot(essol_scs, ylim=(0, 6))

# ╔═╡ 1db274cd-3b01-4f45-8aa0-ffe6878a1dd1
md"""
## Distribution of ages at 5 trillion senescent cells
"""

# ╔═╡ c024dd3a-8fa8-4dbe-9669-0f338416a6af
md"""
Set up a histogram that shows the distribution of ages once the 5 trillion senescent cells are present in the body.
"""

# ╔═╡ 3fe09628-cbfb-4e56-85e8-b97285b32e87
md"""
!!! hints
	- The number of senescent cells of the `i`-th trajoctory can be accessed with: `essol.u[i][:X]`.
	- The index of the first element in the `i`-th trajectory that is greater than 5 can be found with: `findfirst(>(5), essol.u[i][:X])`.
	- An index is a valid index when it is not `nothing`.
	- The time at index position `j` can be accessed with `essol.u[i].t[j]`
	- Appending an element, e.g., `x` to an array `times` can be done as follow: `append!(times, x)`
"""

# ╔═╡ 0f901112-e504-4b7d-b7ec-3a80153903ca
# begin
# 	times = missing                 # make empty vector
# 	for i = missing              # for loop from 1 to 100, default step is 1
# 		# find index of first element that is greater than 5
# 		j = missing
# 		if missing                       # if index is a valid index
# 			missing   # append time to vector times
# 		end
# 	end
# end
begin
	times = []                 # make empty vector
	for i = 1:100              # for loop from 1 to 100, default step is 1
		# find index of first element that is greater than 5
		j = findfirst(>(5), essol_scs.u[i][:X])
		if j != nothing                       # if index is a valid index
			append!(times, essol_scs.u[i].t[j])   # append time to vector times
		end
	end
end

# ╔═╡ 65edc40f-430f-4d1a-9062-17bf8e1d7d59
md"""
Make a histogram with the array `times`. Use `bins=range(0, 120, length=121)` or `bins=0:120`.
"""

# ╔═╡ 17a51f9e-0adf-4802-8d13-c43eb7801bc7
# missing
histogram(times, bins=range(0, 120, length=121), xlab="Age (years)", ylab="Frequency")
# histogram(times, bins=0:120, xlab="Age (years)", ylab="Frequency")

# ╔═╡ dad744a4-0fd2-406d-a600-429ad7619efd
md"""
Check the mean.
"""

# ╔═╡ a41d2ab4-ffad-4b27-9811-ab57c51526cb
# missing
mean(times)

# ╔═╡ 98ca03e1-9c4e-40db-9dd7-68f1898c1ef7
md"""
Check the standard deviation.
"""

# ╔═╡ 7ac54e2a-9c8e-4743-a1a0-38490b04f759
# missing
std(times)

# ╔═╡ 0f9118f5-93f4-42f8-a0cf-754ce9ce9dbe
md"""
Check the minimum value.
"""

# ╔═╡ a07565f7-31d9-4c50-8345-5777adb8a77a
# missing
minimum(times)

# ╔═╡ 6f4253ce-e686-4061-ae13-9b35065464d9
md"""
Check the maximum value.
"""

# ╔═╡ 70628044-aaa5-46d6-a2df-deabbd8f2df8
# missing
maximum(times)

# ╔═╡ 23178795-d7fb-4a23-9649-a1ee9d1706d6
md"""
!!! question
	Interpret the results. Ask yourself the following question:
	1. Suppose that $5$ trillion senescent cells is about the maximum a human body can bear. What is the (approximate) corresponding range of ages?
	2. What is the effect of halving the damage rate $\mu$?
	3. What is the effect of doubling the damage removal rate $\beta$?
	4. What is the effect of halving the noise?
"""

# ╔═╡ a4f5a40e-c4ec-45ef-b44b-b64b71717c25
md"""
Answers:
1. missing
2. missing
3. missing
4. missing
"""
#=
1. The range in age is about [70, 110].
=#

# ╔═╡ Cell order:
# ╠═0c457482-95fe-11ef-0873-65164704c7a6
# ╠═85b6ee78-2f7b-4f83-908b-c48fade56fcb
# ╠═1d0a8fb2-dab9-44ae-a230-f3e92cf3cd6a
# ╠═7a9fa5bc-7bb3-436e-b796-b946c855fe0a
# ╠═22e0eda2-a922-48bb-9b6f-ec0c6dfeb524
# ╟─a176adcc-bde0-4a25-aa80-33379e05e63f
# ╟─3750e702-22b2-4e2f-aa76-d906f5fd315e
# ╟─818ba2af-83dd-4ddb-93e8-d0549ddd806b
# ╟─d917e75d-eefb-415f-ba45-2f1d7c62d63c
# ╟─7bda135b-2e0a-425f-b30c-e5a09dbb8a4a
# ╟─36c678de-56ea-437f-b667-064149075fa1
# ╠═72d60bcc-7d78-45e8-90ca-306bb429aede
# ╟─e560a876-57af-42b7-a16a-8f18676fe8e6
# ╠═27783ea2-c055-4ee6-a5cc-dbc5b4f516f7
# ╟─e881a30a-d622-468e-9be9-69523e5c5ae6
# ╠═8ad575bd-abc1-41b5-8fe3-e28a5c227ef5
# ╟─ef9af0de-b394-4a59-81e8-569716481862
# ╠═858cb14b-ae3d-4b71-94b1-ff7d66a532b9
# ╟─9dc5cfed-fa3d-4637-b1da-37fae3429565
# ╟─705d5c4f-6a7b-4067-80ec-8fb5ccaedb18
# ╠═4e5b490e-67a8-40de-a49f-9f88dcb4b558
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
# ╟─1db274cd-3b01-4f45-8aa0-ffe6878a1dd1
# ╟─c024dd3a-8fa8-4dbe-9669-0f338416a6af
# ╟─3fe09628-cbfb-4e56-85e8-b97285b32e87
# ╠═0f901112-e504-4b7d-b7ec-3a80153903ca
# ╟─65edc40f-430f-4d1a-9062-17bf8e1d7d59
# ╠═17a51f9e-0adf-4802-8d13-c43eb7801bc7
# ╟─dad744a4-0fd2-406d-a600-429ad7619efd
# ╠═a41d2ab4-ffad-4b27-9811-ab57c51526cb
# ╟─98ca03e1-9c4e-40db-9dd7-68f1898c1ef7
# ╠═7ac54e2a-9c8e-4743-a1a0-38490b04f759
# ╟─0f9118f5-93f4-42f8-a0cf-754ce9ce9dbe
# ╠═a07565f7-31d9-4c50-8345-5777adb8a77a
# ╟─6f4253ce-e686-4061-ae13-9b35065464d9
# ╠═70628044-aaa5-46d6-a2df-deabbd8f2df8
# ╟─23178795-d7fb-4a23-9649-a1ee9d1706d6
# ╠═a4f5a40e-c4ec-45ef-b44b-b64b71717c25
