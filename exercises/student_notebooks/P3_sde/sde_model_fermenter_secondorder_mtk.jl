### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 8152f632-af15-4164-a8ff-07c33a9a49b3
using Pkg; Pkg.activate("..")

# ╔═╡ eb142900-1d94-11ef-12ed-6951b45f1817
using Markdown, InteractiveUtils

# ╔═╡ e317c2cc-c41a-4489-b06b-7f07a517e20b
using ModelingToolkit, OrdinaryDiffEq, StochasticDiffEq

# ╔═╡ 2e91e92c-5cd2-4018-8a26-9f0f9e01a175
using ModelingToolkit: t_nounits as t, D_nounits as D

# ╔═╡ e04f782d-67da-4e21-a3bf-d2ddff4bba0b
using StatsPlots

# ╔═╡ 3552fc05-a4dc-4242-91b3-3eb6fb2ea389
using PlutoUI; TableOfContents()

# ╔═╡ a28e7ddf-76e9-4628-888c-e1d838da75ce
md"""
# Exercise: Fermenter - 2nd order kinetics - SDE
"""

# ╔═╡ 959d6307-a30d-4ae7-970d-b2c7584c2c8f
md"""
In a fermenter reactor biomass grows on substrate. The reactor is fed with a inlet flow rate $Q_{in}$ [$L/h$], which consist of a (manipulable) input concentration of substrate $S_{in}$ [$g/L$]. Inside the reactor, biomass, with a concentration of $X$ [$g/L$], is produced through second-order kinetics resulting in the following system of ODEs:

$$\begin{eqnarray*}
\cfrac{dS}{dt} &= \cfrac{Q}{V}\,\left(S_{in} - S\right) - k\,X\,S \\
\cfrac{dX}{dt} &= -\cfrac{Q}{V}\,X + Y\,k\,X\,S
\end{eqnarray*}$$

with $k$ [$L\,gS^{-1}h^{-1}$] the reaction rate constant, and $Y$ [$gX/gS$] the yield coefficient which is defined here by the amount of produced biomass by consumption of one unit of substrate. Futhermore, the reactor is drained with an outlet flow $Q$ [$L/h$], which consist of the current concentrations of substrate $S$ [$g/L$] and biomass $X$ [$g/L$] inside the reactor. The volume $V$ [$L$] of the reactor content is kept constant by setting $Q_{in} = Q$.
"""

# ╔═╡ 6da59e21-6451-4f82-9387-63b90e88de00
md"""
## Implementation of the system
"""

# ╔═╡ 98858b9c-d4f9-451f-a7b9-fcaa012ee28e
md"""
Create a MTK model for the aforementioned problem in order to simulate the evolution of substrate $S$ and biomass $X$ with time as a Stochastic Differential Equation (SDE) problem with a noisy input concentration of the substrate $S_{in}$.
"""

# ╔═╡ f7941ddf-f1e7-41d1-9764-b644dfed68c0
md"""
The parameter values are `k=0.2`, `Y=0.76`, `Q=2.0`, `V=40.0` and `Sin=2.2`. Suppose that at $t=0\;h$ no substrate $S$ is present in the reactor but that there is initially some biomass with a concetration of `0.1`$\;g/L$. Simulate the evolution of $S$ and $X$ during `120.0` hours.
"""

# ╔═╡ 1dea7082-8562-4423-a03e-a8df45f8d57f
md"""
Define the variables `S` and `X` with `@variables`. You can also add the default initial conditions.
"""

# ╔═╡ c48a1118-817f-4c7f-ae7f-dd0460299036
missing

# ╔═╡ 5cd22cdf-708d-4597-8c20-5851af527e86
md"""
Define the parameters of the system with `@parameters`. Add a noise scaling parameter `n=0.5`.
"""

# ╔═╡ 5128e6f4-a9d4-4a03-a93d-949f1b21c150
missing

# ╔═╡ ea307a00-169c-417e-937a-aa00397f518a
md"""
Instantiate a source `B` of Brownian noise (a Wiener process by default) with `@brownian`.
"""

# ╔═╡ e8ea7103-684f-45e6-84e0-da763cab1386
missing

# ╔═╡ 1090d70a-a23b-40ce-9ab9-6a316c455286
md"""
Set up the model equations with the noise `n*B` added to `Sin` and put them in a vector/array. Name the array `eq_ferm`.
"""

# ╔═╡ 07e86c4d-9634-4304-bba3-3caed3fcda07
eqns_ferm = missing

# ╔═╡ 6b65a3a8-fde7-4586-98d6-3f262ec4e946
md"""
Build a system of equations with `@mtkbuild` and the function `System`. Provide the following arguments to `System`: the model equations and `t`. Name the system `sys_ferm`.
"""

# ╔═╡ 003fb585-74c2-4479-8f52-faa336d05b4c
missing

# ╔═╡ 8351a3a0-bbae-4f5c-9655-ded253b798ff
md"""
## Setting initial condition, time span and parameters.
"""

# ╔═╡ 7be20227-fb5e-4dc2-aa1e-1e94db474a69
md"""
Initialize a vector `u0` with the default initial condition, set the timespan for the simulation (we will simulate from `0.0` $y$ to `120.0` $y$), and initialize a vector `parms` with the default parameter values. In that way, later, you can change the initial condition and the parameter values if you want to try other values.
"""

# ╔═╡ ef3f7a41-61ab-4449-8ffc-784cf1e5cbe6
u0 = missing

# ╔═╡ 839a624b-ca62-4c57-9511-207b626ce864
tspan = missing

# ╔═╡ a4d28c40-e315-4bb9-87a5-2b45dd633e5f
parms = missing

# ╔═╡ 0545417b-e3f2-42a7-aa30-15c3bca9ed9c
md"""
## Simulating the system as an SDE problem
"""

# ╔═╡ b2973e86-d89f-476a-8c41-de25a9e5c69d
md"""
Create the SDE problem and store it in `sprob_ferm`.
"""

# ╔═╡ aeddc31e-9de2-4792-a2d8-59a14dfc8173
sprob_ferm = missing

# ╔═╡ 3a981326-2031-4c63-ad31-c44ddd7a88d5
md"""
Solve the SDE problem. Use `EM()` with `dt=0.1`. Store the solution in `ssol_ferm`:
"""

# ╔═╡ 6a69c369-4743-48c1-aed9-4f0ccb095707
ssol_ferm = missing

# ╔═╡ 6b6d2229-913b-41a2-8101-00e9fef0945a
md"""
Plot the results with the option `ylim=(0.0, 2.0)`:
"""

# ╔═╡ 593a0e0a-c4d8-4b12-b38b-15b47705f8a7
missing

# ╔═╡ 177acfd0-eed1-4a1d-925d-5b65e06c35bf
md"""
## Simulating the system as an EnsembleProblem.
"""

# ╔═╡ 08746e97-d794-4261-9ba6-9002cf17e4c1
md"""
Create an `EnsembleProblem` in order to visualize a multiple solutions. Store it in `esprob_ferm`.
"""

# ╔═╡ 3d07836a-60b1-4584-a89a-3d8bbc72b8cd
esprob_ferm = missing

# ╔═╡ 98c4ee2b-20ef-43a1-b420-7644d568810b
md"""
Solve the `EnsembleProblem` using the same solver (and time step) as before, for $100$ trajectories. Store the solution in `essol_ferm`.
"""

# ╔═╡ 683fe575-887b-4bd1-8960-10c04f68354d
essol_ferm = missing

# ╔═╡ 26746cab-d3e7-4a01-bbbd-9fcb49ef652f
md"""
Plot the results. Use as option again `ylim=(0.0,2.0)` and also `linealpha=0.5` (or `la=0.5`) to modify the line boldness.
"""

# ╔═╡ 4fcc1d0a-9d30-4056-b8a3-3d802edc42e5
missing

# ╔═╡ b7580824-8bea-4fd4-88c7-6d8c1eb3f00d
md"""
Check out the end value of $X$ in the first simulation with `essol_ferm.u[1][:X][end]`.
"""

# ╔═╡ 28da57db-ff57-4ade-887c-a132a616abd4
missing

# ╔═╡ 440c0d44-60c4-444f-b6b0-b53c6d9abb4f
md"""
Create an array/vector with all the end values of $X$. Use therefore a `for`-loop and the `append!` function. Put all values in `Xeq_values`.
"""

# ╔═╡ 2b61187f-12f4-4b80-b927-f6fb587b486a
# begin
# 	Xeq_values = missing
# 	for i=missing
# 		missing
# 	end
# end

# ╔═╡ 705e51c1-2cb9-4162-b356-9364715bbfa8
md"""
Plot a histogram in the range `[1.0, 1.8]` with a bin size of `0.01` (cf. `bins=1.0:0.01:1.8`).
"""

# ╔═╡ 89f2a414-ba53-43aa-bf61-1a7fbe252986
missing

# ╔═╡ 2d6202ca-fb4a-4b44-847f-07abe120e618
md"""
Calculate the mean.
"""

# ╔═╡ ac5a78ad-6575-4c0b-9f9e-a632b648552a
missing

# ╔═╡ 638304a0-1036-4c54-a40b-4cb8baf826c6
md"""
Calculate the standard deviation.
"""

# ╔═╡ fffe5c34-f0e0-400e-bb0f-44f9d9595919
missing

# ╔═╡ Cell order:
# ╠═eb142900-1d94-11ef-12ed-6951b45f1817
# ╠═8152f632-af15-4164-a8ff-07c33a9a49b3
# ╠═e317c2cc-c41a-4489-b06b-7f07a517e20b
# ╠═2e91e92c-5cd2-4018-8a26-9f0f9e01a175
# ╠═e04f782d-67da-4e21-a3bf-d2ddff4bba0b
# ╠═3552fc05-a4dc-4242-91b3-3eb6fb2ea389
# ╟─a28e7ddf-76e9-4628-888c-e1d838da75ce
# ╟─959d6307-a30d-4ae7-970d-b2c7584c2c8f
# ╟─6da59e21-6451-4f82-9387-63b90e88de00
# ╟─98858b9c-d4f9-451f-a7b9-fcaa012ee28e
# ╟─f7941ddf-f1e7-41d1-9764-b644dfed68c0
# ╟─1dea7082-8562-4423-a03e-a8df45f8d57f
# ╠═c48a1118-817f-4c7f-ae7f-dd0460299036
# ╟─5cd22cdf-708d-4597-8c20-5851af527e86
# ╠═5128e6f4-a9d4-4a03-a93d-949f1b21c150
# ╟─ea307a00-169c-417e-937a-aa00397f518a
# ╠═e8ea7103-684f-45e6-84e0-da763cab1386
# ╟─1090d70a-a23b-40ce-9ab9-6a316c455286
# ╠═07e86c4d-9634-4304-bba3-3caed3fcda07
# ╟─6b65a3a8-fde7-4586-98d6-3f262ec4e946
# ╠═003fb585-74c2-4479-8f52-faa336d05b4c
# ╟─8351a3a0-bbae-4f5c-9655-ded253b798ff
# ╟─7be20227-fb5e-4dc2-aa1e-1e94db474a69
# ╠═ef3f7a41-61ab-4449-8ffc-784cf1e5cbe6
# ╠═839a624b-ca62-4c57-9511-207b626ce864
# ╠═a4d28c40-e315-4bb9-87a5-2b45dd633e5f
# ╟─0545417b-e3f2-42a7-aa30-15c3bca9ed9c
# ╟─b2973e86-d89f-476a-8c41-de25a9e5c69d
# ╠═aeddc31e-9de2-4792-a2d8-59a14dfc8173
# ╟─3a981326-2031-4c63-ad31-c44ddd7a88d5
# ╠═6a69c369-4743-48c1-aed9-4f0ccb095707
# ╟─6b6d2229-913b-41a2-8101-00e9fef0945a
# ╠═593a0e0a-c4d8-4b12-b38b-15b47705f8a7
# ╟─177acfd0-eed1-4a1d-925d-5b65e06c35bf
# ╟─08746e97-d794-4261-9ba6-9002cf17e4c1
# ╠═3d07836a-60b1-4584-a89a-3d8bbc72b8cd
# ╟─98c4ee2b-20ef-43a1-b420-7644d568810b
# ╠═683fe575-887b-4bd1-8960-10c04f68354d
# ╟─26746cab-d3e7-4a01-bbbd-9fcb49ef652f
# ╠═4fcc1d0a-9d30-4056-b8a3-3d802edc42e5
# ╟─b7580824-8bea-4fd4-88c7-6d8c1eb3f00d
# ╠═28da57db-ff57-4ade-887c-a132a616abd4
# ╟─440c0d44-60c4-444f-b6b0-b53c6d9abb4f
# ╠═2b61187f-12f4-4b80-b927-f6fb587b486a
# ╟─705e51c1-2cb9-4162-b356-9364715bbfa8
# ╠═89f2a414-ba53-43aa-bf61-1a7fbe252986
# ╟─2d6202ca-fb4a-4b44-847f-07abe120e618
# ╠═ac5a78ad-6575-4c0b-9f9e-a632b648552a
# ╟─638304a0-1036-4c54-a40b-4cb8baf826c6
# ╠═fffe5c34-f0e0-400e-bb0f-44f9d9595919
