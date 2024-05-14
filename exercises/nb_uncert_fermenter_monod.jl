### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 9d8acf40-e635-4dc4-9938-ec63ab68e3bd
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 1ec7def2-0c43-11ef-0850-ddefa680d2a4
using Markdown

# ╔═╡ 2861bc17-7810-461e-8857-a55b6392977a
using InteractiveUtils

# ╔═╡ 093679cc-122b-4b0a-84c3-15cf3a6f5c1b
using Catalyst

# ╔═╡ 0720732e-105a-40c2-bf45-f8fbbc10e5bc
using DifferentialEquations, Plots

# ╔═╡ 28bb60c9-ab56-497f-b039-0d6577501be5
using Measurements

# ╔═╡ 437364b6-a953-428b-a891-e17ffd531656
md"
### Exercise: Fermenter - Monod kinetics - Uncertainty analysis
"

# ╔═╡ c44e10f3-8452-4a60-b840-f4a7e6ec3238
md"
In one of the previous practica we were introduced to a fermenter in which biomass $X$ [$g/L$] grows by breaking down substrate $S$ [$g/L$]. The reactor is fed with a inlet flow rate $Q_{in}$ [$L/h$], which consist of a (manipulable) input concentration of substrate $S_{in}$ [$g/L$]. This process was modelled using Monod kinetics, resulting in the model below:

$$\begin{eqnarray*}
%S \xrightarrow[\quad\quad]{\beta} Y \, X
S \xrightarrow[\quad\quad]{r} Y \, X \quad\quad\quad\quad r = \mu \, X \quad \textrm{with} \quad \mu = \mu_{max} \, \cfrac{S}{S + K_s}
\end{eqnarray*}$$
"

# ╔═╡ 356822df-bad3-41f9-8556-de04234b0a09
md"
The *reaction network object* for this model could be set-up as:
"

# ╔═╡ c551c33e-95b0-4ea8-85d4-6c68db72b8bc
fermenter_monod = @reaction_network begin
    X * mm(S, μmax, Ks), S --> Y*X
    Q/V, (S, X) --> ∅
    Q/V*Sin, ∅ --> S
end

# ╔═╡ 69f8c265-fe91-4e4b-971b-395325d4474e
md"
which resulted in the following differential equations:
"

# ╔═╡ 2806654f-7a24-4f8e-8b5a-3129b2431c95
md"
$$\begin{eqnarray*}
\cfrac{dS}{dt} &=& \cfrac{Q}{V} \left(S_{in} - S \right) - \mu_{max}\cfrac{S}{S + K_s} X\\
\cfrac{dX}{dt} &=& -\cfrac{Q}{V} X + Y \mu_{max}\cfrac{S}{S + K_s} X
\end{eqnarray*}$$
"

# ╔═╡ 5802078d-12bc-489e-8566-15f77a11e8ba
md"
Assume the uncertainties in the following parameter values:

-  $\mu_{max} = 0.30 \pm 0.06\;h^{-1}$
-  $K_s = 0.15 \pm 0.03 \;g/L$
-  $S_{in} = 2.2 \pm 0.4\;g/L$

and that the uncertainty in the other parameters values: $Y = 0.80$, $Q = 2.0\;L/h$, $V = 40.0\;L$ are negligible. Suppose that at $t$ no substrate $S$ is present in the reactor but that there is initially some biomass with a concetration of $0.01\;g/L$. Perform an uncertainty analysis by plotting the uncertainty bands on the simulation results of $S$ and $X$ in a timespan of $[0, 100]\,h$.

Interpret your results.
"

# ╔═╡ 562e0006-bb79-4895-98b4-c376ea6ce856
md"
Initialize a vector `u₀` with the initial conditions, and set the timespan:
"

# ╔═╡ 6afb23f2-3da7-43e0-8a3b-fc4bb5eb9bde
# u₀ = missing                 # Uncomment and complete the instruction
u₀ = [:S => 0.0, :X => 0.01]

# ╔═╡ 1b09c8cb-d3b7-4a15-b6ed-d03019c5f57b
# tspan = missing              # Uncomment and complete the instruction
tspan = (0.0, 100.0)

# ╔═╡ 40a8c697-56d1-4488-b089-00ddf41cdb5a
md"
We initialize a vector `params_uncert` with the parameter values and their corresponding uncertainty:
"

# ╔═╡ 62f4049b-617f-4268-beba-b2c634b46642
# params_uncert = missing       # Uncomment and complete the instruction
params_uncert = [:μmax => 0.30±0.06, :Ks => 0.15±0.03, :Y => 0.80, :Q => 2, :V => 40, :Sin => 2.2±0.4]
# params_uncert = [:μmax => 0.30±0.06, :Ks => 0.15, :Y => 0.80, :Q => 2, :V => 40, :Sin => 2.2]
# params_uncert = [:μmax => 0.30, :Ks => 0.15±0.03, :Y => 0.80, :Q => 2, :V => 40, :Sin => 2.2]
# params_uncert = [:μmax => 0.30, :Ks => 0.15, :Y => 0.80, :Q => 2, :V => 40, :Sin => 2.2±0.4]

# ╔═╡ 0e4c73d3-1630-463b-994d-073170c1a6c3
md"
We create the corresponding ODE problem and store it in `oprob_uncert`:
"

# ╔═╡ 3cad11d6-c2eb-4865-858a-d407aa25d780
# oprob_uncert = missing        # Uncomment and complete the instruction
oprob_uncert = ODEProblem(fermenter_monod, u₀, tspan, params_uncert, combinatoric_ratelaws=false)

# ╔═╡ a7b3b68f-d34e-4e9a-aa44-224919c44fec
md"
We solve the ODE problem. Use `Tsit5()` and `saveat=2.0`. Store the solution in `osol_uncert`:
"

# ╔═╡ 40168902-0780-4f92-8eca-23982c80511e
# osol_uncert = missing         # Uncomment and complete the instruction
osol_uncert = solve(oprob_uncert, Tsit5(), saveat=2.0)

# ╔═╡ 6b2abf2a-f1f3-4c2d-96e8-cc626848fa68
md"
Plot the results (simulation of the output variables $S$ and $X$ together with their uncertainty band):
"

# ╔═╡ 4d753fdb-089d-41dc-89ef-45a1733af786
plot(osol_uncert)

# ╔═╡ 4f3900d7-4ff8-4b4b-a59a-f7079ac6c457
md"
Try to relate the local sensitivity analysis to the uncertainty analysis. Hence, study the effect of the individual parameter uncertainties on the output variables $S$ and $X$ and compare with your local sensitivity results of the corresponding parameter.

In order to do that, analyse the effect on the uncertainty bands for $S$ and $X$ by taking one uncertainty on a parameter at a time. In other words, analyse the subsequent cases seperately:
- Assume uncertainty only in $\mu_{max}$
- Assume uncertainty only in $K_s$
- Assume uncertainty only in $S_{in}$

Draw your conclusions:
- Answer: missing
"

# ╔═╡ Cell order:
# ╠═1ec7def2-0c43-11ef-0850-ddefa680d2a4
# ╠═2861bc17-7810-461e-8857-a55b6392977a
# ╠═9d8acf40-e635-4dc4-9938-ec63ab68e3bd
# ╠═093679cc-122b-4b0a-84c3-15cf3a6f5c1b
# ╠═0720732e-105a-40c2-bf45-f8fbbc10e5bc
# ╠═28bb60c9-ab56-497f-b039-0d6577501be5
# ╠═437364b6-a953-428b-a891-e17ffd531656
# ╠═c44e10f3-8452-4a60-b840-f4a7e6ec3238
# ╠═356822df-bad3-41f9-8556-de04234b0a09
# ╠═c551c33e-95b0-4ea8-85d4-6c68db72b8bc
# ╠═69f8c265-fe91-4e4b-971b-395325d4474e
# ╠═2806654f-7a24-4f8e-8b5a-3129b2431c95
# ╠═5802078d-12bc-489e-8566-15f77a11e8ba
# ╠═562e0006-bb79-4895-98b4-c376ea6ce856
# ╠═6afb23f2-3da7-43e0-8a3b-fc4bb5eb9bde
# ╠═1b09c8cb-d3b7-4a15-b6ed-d03019c5f57b
# ╠═40a8c697-56d1-4488-b089-00ddf41cdb5a
# ╠═62f4049b-617f-4268-beba-b2c634b46642
# ╠═0e4c73d3-1630-463b-994d-073170c1a6c3
# ╠═3cad11d6-c2eb-4865-858a-d407aa25d780
# ╠═a7b3b68f-d34e-4e9a-aa44-224919c44fec
# ╠═40168902-0780-4f92-8eca-23982c80511e
# ╠═6b2abf2a-f1f3-4c2d-96e8-cc626848fa68
# ╠═4d753fdb-089d-41dc-89ef-45a1733af786
# ╠═4f3900d7-4ff8-4b4b-a59a-f7079ac6c457
