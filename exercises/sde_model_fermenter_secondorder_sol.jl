### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 8152f632-af15-4164-a8ff-07c33a9a49b3
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ eb142900-1d94-11ef-12ed-6951b45f1817
using Markdown

# ╔═╡ 7139a11b-64db-46b9-a41c-dca83a9eab26
using InteractiveUtils

# ╔═╡ 1ca1d1db-1fdd-4a76-b350-126add7e013c
using Catalyst

# ╔═╡ e04f782d-67da-4e21-a3bf-d2ddff4bba0b
using DifferentialEquations, Plots

# ╔═╡ b02df672-970c-4fa1-9107-cce8ad25d878
using Distributions

# ╔═╡ a28e7ddf-76e9-4628-888c-e1d838da75ce
md"""
# Exercise: Fermenter - 2nd order kinetics - SDE
"""

# ╔═╡ 959d6307-a30d-4ae7-970d-b2c7584c2c8f
md"""
In a fermenter reactor biomass grows on substrate. The reactor is fed with a inlet flow rate $Q_{in}$ [$L/h$], which consist of a (manipulable) input concentration of substrate $S_{in}$ [$g/L$]. Inside the reactor, biomass, with a concentration of $X$ [$g/L$], is produced through second-order kinetics:

$$\begin{eqnarray*}
%S \xrightarrow[\quad\quad]{\beta} Y \, X
S \xrightarrow[\quad\quad]{r} Y \, X \quad\quad\quad\quad r = k \, S\,X
\end{eqnarray*}$$

with $k$ [$L\,gS^{-1}h^{-1}$] the reaction rate constant, and $Y$ [$gX/gS$] the yield coefficient which is defined here by the amount of produced biomass by consumption of one unit of substrate. Futhermore, the reactor is drained with an outlet flow $Q$ [$L/h$], which consist of the current concentrations of substrate $S$ [$g/L$] and biomass $X$ [$g/L$] inside the reactor. The volume $V$ [$L$] of the reactor content is kept constant by setting $Q_{in} = Q$.
"""

# ╔═╡ 98858b9c-d4f9-451f-a7b9-fcaa012ee28e
md"""
Create a *reaction network object* model for the aforementioned problem in order to simulate the evolution of substrate $S$ and biomass $X$ with time as a Stochastic Differential Equation (SDE) problem with noise scaling. Name it `fermenter_sde_secondorder`.
"""

# ╔═╡ d2e2680c-01c6-449a-bb5f-7472bc1de243
md"""
Assign the following noise scaling values:
- `η = 0.10` for the main *reaction* (default value for `η`)
- noise scaling of `0.05` for the *reaction* describing the inlet $S_{in}$
- noise scaling of `0.0` for the remaining *reactions*
"""

# ╔═╡ f7941ddf-f1e7-41d1-9764-b644dfed68c0
md"""
The parameter values are $k =$ 0.2, $Y =$ 0.76, $Q =$ 2.0, $V =$ 40.0 and $S_{in} =$ 2.2$\;g/L$. Suppose that at $t=$0$\;h$ no substrate $S$ is present in the reactor but that there is initially some biomass with a concetration of 0.1$\;g/L$. Simulate the evolution of $S$ and $X$ during $120$ hours.
"""

# ╔═╡ 6b627d84-b6a5-444d-8163-40a4cab181bd
fermenter_sde_secondorder = @reaction_network begin
    @parameters η=0.10
    @default_noise_scaling η
	k, S + X --> (1 + Y)*X
	# Alternatively:
    # k*S*X, S => Y*X
    Q/V*Sin, 0 --> S, [noise_scaling = 0.05]
	Q/V, (S, X) --> 0, [noise_scaling = 0.00]
end

# ╔═╡ 47ca3573-691c-4127-85b9-d5b5a1a23fbb
md"""
Convert the system to a symbolic differential equation model and verify, by analyzing the differential equation, that your model is correctly implemented.

Hint:
- Use the option: `combinatoric_ratelaws=false`

For information about this option, click [here](https://docs.sciml.ai/Catalyst/stable/introduction_to_catalyst/introduction_to_catalyst/#introduction_to_catalyst_ratelaws).
"""

# ╔═╡ e1f11068-e489-4c54-a30d-981f5cb19b47
osys = convert(ODESystem, fermenter_sde_secondorder, combinatoric_ratelaws=false)

# ╔═╡ a65c84c8-d2ec-44e4-9a07-cbd59f190c57
md"""
Initialize a vector `u0` with the initial conditions:
"""

# ╔═╡ ef3f7a41-61ab-4449-8ffc-784cf1e5cbe6
u0 = [:S => 0.0, :X => 0.1]

# ╔═╡ 9d113e7f-8499-4b2a-a884-d34ce3da0b82
md"""
Set the timespan for the simulation:
"""

# ╔═╡ 839a624b-ca62-4c57-9511-207b626ce864
tspan = (0.0, 120.0)

# ╔═╡ a4d28c40-e315-4bb9-87a5-2b45dd633e5f
parms = [:k => 0.2, :Y => 0.76, :Q => 2, :V => 40, :Sin => 2.2, :η => 0.10]

# ╔═╡ b2973e86-d89f-476a-8c41-de25a9e5c69d
md"""
Create the SDE problem and store it in `sprob`.\
Hint:
- Use the option: `combinatoric_ratelaws=false`
"""

# ╔═╡ aeddc31e-9de2-4792-a2d8-59a14dfc8173
sprob = SDEProblem(fermenter_sde_secondorder, u0, tspan, parms, combinatoric_ratelaws=false);

# ╔═╡ 3a981326-2031-4c63-ad31-c44ddd7a88d5
md"""
Solve the SDE problem. Use `EM()` with `dt=0.1`. Store the solution in `ssol`:
"""

# ╔═╡ 6a69c369-4743-48c1-aed9-4f0ccb095707
ssol = solve(sprob, EM(), dt=0.1)

# ╔═╡ 6b6d2229-913b-41a2-8101-00e9fef0945a
md"""
Plot the results with the option `ylim=(0.0, 2.0)`:
"""

# ╔═╡ 593a0e0a-c4d8-4b12-b38b-15b47705f8a7
plot(ssol, ylim=(0.0, 2.0))

# ╔═╡ 08746e97-d794-4261-9ba6-9002cf17e4c1
md"""
Create an `EnsembleProblem` in order to visualize a multiple solutions. Store it in `esprob`.
"""

# ╔═╡ 3d07836a-60b1-4584-a89a-3d8bbc72b8cd
esprob = EnsembleProblem(sprob)

# ╔═╡ 98c4ee2b-20ef-43a1-b420-7644d568810b
md"""
Solve the `EnsembleProblem` using the same solver (and time step) as before, for $100$ trajectories. Store the solution in `essol`.
"""

# ╔═╡ 683fe575-887b-4bd1-8960-10c04f68354d
essol = solve(esprob, EM(), dt=0.1, trajectories=100)

# ╔═╡ 26746cab-d3e7-4a01-bbbd-9fcb49ef652f
md"""
Plot the results. Use as option again `ylim=(0.0,2.0)` and also `linealpha=0.5` (or `la=0.5`) to modify the line boldness.
"""

# ╔═╡ 4fcc1d0a-9d30-4056-b8a3-3d802edc42e5
plot(essol, ylim=(0.0,2.0), linealpha=0.5)
# plot(essol, ylim=(0.0,2.0), la=0.5)

# ╔═╡ 28da57db-ff57-4ade-887c-a132a616abd4
essol.u[1][:X][end]

# ╔═╡ 2b61187f-12f4-4b80-b927-f6fb587b486a
begin
	Xeq_values = []
	for i=1:100
		append!(Xeq_values, essol.u[i][:X][end])
	end
end

# ╔═╡ 89f2a414-ba53-43aa-bf61-1a7fbe252986
histogram(Xeq_values, bins=1.0:0.01:1.8)

# ╔═╡ ac5a78ad-6575-4c0b-9f9e-a632b648552a
mean(Xeq_values)

# ╔═╡ fffe5c34-f0e0-400e-bb0f-44f9d9595919
std(Xeq_values)

# ╔═╡ Cell order:
# ╠═eb142900-1d94-11ef-12ed-6951b45f1817
# ╠═7139a11b-64db-46b9-a41c-dca83a9eab26
# ╠═8152f632-af15-4164-a8ff-07c33a9a49b3
# ╠═1ca1d1db-1fdd-4a76-b350-126add7e013c
# ╠═e04f782d-67da-4e21-a3bf-d2ddff4bba0b
# ╠═b02df672-970c-4fa1-9107-cce8ad25d878
# ╟─a28e7ddf-76e9-4628-888c-e1d838da75ce
# ╟─959d6307-a30d-4ae7-970d-b2c7584c2c8f
# ╟─98858b9c-d4f9-451f-a7b9-fcaa012ee28e
# ╟─d2e2680c-01c6-449a-bb5f-7472bc1de243
# ╟─f7941ddf-f1e7-41d1-9764-b644dfed68c0
# ╠═6b627d84-b6a5-444d-8163-40a4cab181bd
# ╟─47ca3573-691c-4127-85b9-d5b5a1a23fbb
# ╠═e1f11068-e489-4c54-a30d-981f5cb19b47
# ╟─a65c84c8-d2ec-44e4-9a07-cbd59f190c57
# ╠═ef3f7a41-61ab-4449-8ffc-784cf1e5cbe6
# ╟─9d113e7f-8499-4b2a-a884-d34ce3da0b82
# ╠═839a624b-ca62-4c57-9511-207b626ce864
# ╠═a4d28c40-e315-4bb9-87a5-2b45dd633e5f
# ╟─b2973e86-d89f-476a-8c41-de25a9e5c69d
# ╠═aeddc31e-9de2-4792-a2d8-59a14dfc8173
# ╟─3a981326-2031-4c63-ad31-c44ddd7a88d5
# ╠═6a69c369-4743-48c1-aed9-4f0ccb095707
# ╟─6b6d2229-913b-41a2-8101-00e9fef0945a
# ╠═593a0e0a-c4d8-4b12-b38b-15b47705f8a7
# ╟─08746e97-d794-4261-9ba6-9002cf17e4c1
# ╠═3d07836a-60b1-4584-a89a-3d8bbc72b8cd
# ╟─98c4ee2b-20ef-43a1-b420-7644d568810b
# ╠═683fe575-887b-4bd1-8960-10c04f68354d
# ╟─26746cab-d3e7-4a01-bbbd-9fcb49ef652f
# ╠═4fcc1d0a-9d30-4056-b8a3-3d802edc42e5
# ╠═28da57db-ff57-4ade-887c-a132a616abd4
# ╠═2b61187f-12f4-4b80-b927-f6fb587b486a
# ╠═89f2a414-ba53-43aa-bf61-1a7fbe252986
# ╠═ac5a78ad-6575-4c0b-9f9e-a632b648552a
# ╠═fffe5c34-f0e0-400e-bb0f-44f9d9595919
