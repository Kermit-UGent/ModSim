### A Pluto.jl notebook ###
# v0.20.28

using Markdown
using InteractiveUtils

# ╔═╡ 9d8acf40-e635-4dc4-9938-ec63ab68e3bd
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ 1ec7def2-0c43-11ef-0850-ddefa680d2a4
using Markdown

# ╔═╡ 2861bc17-7810-461e-8857-a55b6392977a
using InteractiveUtils

# ╔═╡ 093679cc-122b-4b0a-84c3-15cf3a6f5c1b
using Catalyst

# ╔═╡ 0720732e-105a-40c2-bf45-f8fbbc10e5bc
using OrdinaryDiffEq

# ╔═╡ 28bb60c9-ab56-497f-b039-0d6577501be5
using StatsPlots, PlutoUI; TableOfContents()

# ╔═╡ 0a65afbc-643a-4fe4-b22d-c5fe54e3be1f
using Turing

# ╔═╡ d1b6739f-1e0c-452d-8763-717dd650b2ec
using ForwardDiff

# ╔═╡ 162b330b-02f7-4179-b950-d960eca68206
using Random

# ╔═╡ 437364b6-a953-428b-a891-e17ffd531656
md"""
# Exercise: Fermenter - Switching to a mutant strain?
"""

# ╔═╡ c44e10f3-8452-4a60-b840-f4a7e6ec3238
md"""
As a process engineer, you are responsible for a continuously stirred fermenter in which biomass $X$ [$g/L$] grows by breaking down substrate $S$ [$g/L$]. The reactor is fed with an inlet flow rate $Q$ [$L/h$] containing substrate at a concentration $S_{in}$ [$g/L$]. The yield $Y$ [$g$ biomass$/g$ substrate] tells how much biomass is formed per gram of substrate consumed. As in the previous practica, the growth is modelled with Monod kinetics:

$$\begin{eqnarray*}
S \xrightarrow[\quad\quad]{r} Y \, X \quad\quad\quad\quad r = \mu \, X \quad \textrm{with} \quad \mu = \mu_{max} \, \cfrac{S}{S + K_s}
\end{eqnarray*}$$
"""

# ╔═╡ f3065171-5526-476a-bfd7-a2b95bb24d61
md"""
## The question

A colleague has found a mutant strain with a **higher yield** $Y$, and management asks you a simple question: *should we switch to the mutant?*

The catch: the maximal growth rate $\mu_{max}$ of the mutant is not known. A higher yield often comes at the cost of slower growth, and a slower grower leaves more substrate unused in the outflow. So the mutant could be better or worse than the current strain.

You will answer this question with the tools of this practical:
1. **Uncertainty propagation**: what does our uncertainty about $\mu_{max}$ mean for the biomass? How sure are we that the mutant is better?
2. **Sensitivity analysis**: when does the biomass depend most strongly on $\mu_{max}$? That is the moment where a measurement tells us the most.
3. **Bayesian inference**: update our belief about $\mu_{max}$ with a few measurements, and propagate the reduced uncertainty again.
"""

# ╔═╡ 356822df-bad3-41f9-8556-de04234b0a09
md"""
The *reaction network object* for this model is:
"""

# ╔═╡ c551c33e-95b0-4ea8-85d4-6c68db72b8bc
fermenter_monod = @reaction_network begin
    X * mm(S, μmax, Ks), S => Y*X
    Q/V, (S, X) --> 0
    Q/V*Sin, ∅ --> S
end

# ╔═╡ 69f8c265-fe91-4e4b-971b-395325d4474e
md"""
which results in the following differential equations:
"""

# ╔═╡ 2806654f-7a24-4f8e-8b5a-3129b2431c95
md"""
$$\begin{eqnarray*}
\cfrac{dS}{dt} &=& \cfrac{Q}{V} \left(S_{in} - S \right) - \mu_{max}\cfrac{S}{S + K_s} X\\
\cfrac{dX}{dt} &=& -\cfrac{Q}{V} X + Y \mu_{max}\cfrac{S}{S + K_s} X
\end{eqnarray*}$$
"""

# ╔═╡ 5d6a881d-421e-4dea-9668-38b3f3378bf1
# osys = missing
osys = convert(ODESystem, fermenter_monod)

# ╔═╡ 5802078d-12bc-489e-8566-15f77a11e8ba
md"""
## The current strain

For the strain currently in use, the parameters are well known:
-  $\mu_{max} = 0.40\;h^{-1}$
-  $K_s = 0.015 \;g/L$
-  $Y = 0.67$
-  $S_{in} = 0.022 \;g/L$
-  $Q = 2.0\;L/h$
-  $V = 40.0\;L$

At $t=0$ no substrate is present in the reactor and the biomass concentration is $0.0005\;g/L$. We judge a strain by the biomass concentration in the outflow after $100\,h$.

Initialize a vector `u0` with the initial conditions, set the timespan `tspan` to $[0, 100]\,h$, the time step `dt` to `0.5`, and initialize a vector `parms` with the parameter values.
"""

# ╔═╡ 6afb23f2-3da7-43e0-8a3b-fc4bb5eb9bde
# u0 = missing
u0 = [:S => 0.0, :X => 0.0005]

# ╔═╡ 1b09c8cb-d3b7-4a15-b6ed-d03019c5f57b
# tspan = missing
tspan = (0.0, 100.0)

# ╔═╡ 9e3c427b-8a56-4787-9b28-43ec3932ddbe
# dt = missing
dt = 0.5

# ╔═╡ 62f4049b-617f-4268-beba-b2c634b46642
# parms = missing
parms = [:μmax => 0.40, :Ks => 0.015, :Y => 0.67, :Q => 2.0, :V => 40.0, :Sin => 0.022]

# ╔═╡ 0e4c73d3-1630-463b-994d-073170c1a6c3
md"""
Create the ODE problem and store it in `oprob`. Next, solve the ODE problem using `Tsit5()` and `saveat=dt`, and store the solution in `osol`. Finally, plot the results.
"""

# ╔═╡ 3cad11d6-c2eb-4865-858a-d407aa25d780
# oprob = missing
oprob = ODEProblem(fermenter_monod, u0, tspan, parms, combinatoric_ratelaws=false)

# ╔═╡ 40168902-0780-4f92-8eca-23982c80511e
# osol = missing
osol = solve(oprob, Tsit5(), saveat=dt)

# ╔═╡ 4d753fdb-089d-41dc-89ef-45a1733af786
# missing
plot(osol, title="Current strain", xlabel="Time (hours)", ylabel="concentration (g/l)", linewidth=2)

# ╔═╡ ca660e2c-9348-40b1-b92d-415068bea514
md"""
Store the biomass concentration at the end of the simulation in `X_current`. This is the value the mutant has to beat.
"""

# ╔═╡ ec042e61-04b8-4d51-8e7d-35e1646cc5ba
# X_current = missing
X_current = osol[:X][end]

# ╔═╡ 9ee22028-7a37-4510-9b61-198d3bd99a73
md"""
## The mutant strain
"""

# ╔═╡ cc9ac43b-448a-4a67-ba0b-8fdf21a4236e
md"""
Some batch tests show that the mutant has a yield of $Y = 0.70$ instead of $0.67$. All other parameters are the same as for the current strain, except $\mu_{max}$, which is unknown.

From the literature on mutants with an improved yield, you expect $\mu_{max}$ to be around $0.35\;h^{-1}$, give or take $15\%$. You express this **prior belief** with a log-normal distribution `LogNormal(log(0.35), 0.15)`.
"""

# ╔═╡ f5482531-b23f-472b-968f-20f7a6471778
begin
	plot(LogNormal(log(0.35), 0.15), xlabel="μmax (1/h)", ylabel="density", label="Prior μmax mutant", linewidth=2)
	vline!([0.40], label="μmax current strain", linestyle=:dash)
end

# ╔═╡ ff21b3cd-64b1-4d09-8eb1-afc617585fc8
md"""
Create an ODE problem `oprob_mut` for the mutant by remaking `oprob` with $Y = 0.70$ and the expected value $\mu_{max} = 0.35$.

!!! note
	Remember that you can remake an ODEProblem with new parameter values by using `remake(oprob, p = [...])`.
"""

# ╔═╡ 225eaaff-42a2-4aca-94e1-afcd988de4eb
# oprob_mut = missing
oprob_mut = remake(oprob, p = [:μmax => 0.35, :Y => 0.70])

# ╔═╡ 18ffe877-c579-4b4c-80bd-b5d50baf36bb
md"""
## Monte Carlo error propagation

We don't know $\mu_{max}$, but we do know which values are plausible. Using Monte Carlo simulation, we draw many values from the prior, simulate the mutant for each value, and look at the spread of the outputs. This shows how the uncertainty on the parameter **propagates** to the uncertainty on the output.

Make a Turing model `mutant_deviation` in which you implement the prior for `μmax_dev` and return the solution of the remade ODE problem.
"""

# ╔═╡ 7821dc21-3e6b-4ad1-829e-f29dcd020144
# @model function mutant_deviation()
# 	μmax_dev ~ missing
# 	new_prob = missing
# 	return missing
# end
@model function mutant_deviation()
	μmax_dev ~ LogNormal(log(0.35), 0.15)
	new_prob = remake(oprob_mut, p = [:μmax => μmax_dev])
	return solve(new_prob, Tsit5(), saveat=dt)
end

# ╔═╡ 3051d39a-fc2f-40e4-a688-4d6e75cfc2d1
md"""
Now sample 500 solutions from your `mutant_deviation` model.
"""

# ╔═╡ 33df05c9-08d3-4107-9240-931658cda2a4
# missing
begin
	dev_model = mutant_deviation()
	chain = sample(dev_model, Prior(), 500)
	solutions = generated_quantities(dev_model, chain)
end;

# ╔═╡ 02c5cf65-ebd3-477b-9a3e-8855bd46f55c
md"""
Plot the biomass of all Monte Carlo simulations together with the biomass of the current strain. Use `idxs=:X`, `label=false`, `color=:gray`, `alpha=0.4` and `linestyle=:dot`.
"""

# ╔═╡ aa833d92-567f-4971-a360-e60bce424c8e
# begin
# 	p1 = plot(title="Monte Carlo of the mutant", ylabel="X (g/l)", xlabel="Time (hours)")
# 	for missing in missing
# 		missing
# 	end
# 	plot!(osol, idxs=:X, label="Current strain", color=:red, linewidth=3)
# 	p1
# end
begin
	p1 = plot(title="Monte Carlo of the mutant", ylabel="X (g/l)", xlabel="Time (hours)")
	for solution in solutions
		plot!(solution, idxs=:X, label=false, color=:gray, alpha=0.4, linestyle=:dot)
	end
	plot!(osol, idxs=:X, label="Current strain", color=:red, linewidth=3)
	p1
end

# ╔═╡ ffc66a94-5435-4ca6-af87-76c85a599906
md"""
Create a vector `x1` that stores the biomass at the end of the simulations, make a histogram of it, and compute the probability `P_better_prior` that the mutant produces more biomass than the current strain.

!!! note
	The fraction of elements of a vector `x` that are larger than `a` can be computed with `mean(x .> a)`.
"""

# ╔═╡ d3afd5b4-034f-415c-a64a-776deed74d5c
# x1 = missing
x1 = [solution[:X][end] for solution in solutions];

# ╔═╡ 1d95728e-079e-4e9b-952a-c4f692661813
# missing
begin
	histogram(x1, title="Prior histogram of X(100 h)", xlabel="X (g/l)", ylabel="counts", label="Mutant", legend=:topleft)
	vline!([X_current], linewidth=3, color=:red, label="Current strain")
end

# ╔═╡ 02213469-a233-4759-a18d-29a76a4c13c6
# P_better_prior = missing
P_better_prior = mean(x1 .> X_current)

# ╔═╡ 97d2065c-6234-4b13-8e4c-6782b8301416
md"""
!!! questions
	- When is the spread on $X$ largest, and when is it smallest? Can you explain this?
	- How likely is it that the mutant beats the current strain? Would you advise management to switch based on this?
"""

# ╔═╡ d2acbfb9-a7ff-48b0-9d9a-9492ecb9c1bf
md"""
Answers:
- missing
- missing
"""
#=
- The spread is largest during the start-up phase (around 40-60 h). There, the growth rate decides how far the population has grown. Towards the end the spread shrinks again, because almost all substrate is consumed and X tends to Y*Sin, regardless of μmax.
- About 65-70%. That is better than a coin flip, but switching a production process is costly, and a 1 in 3 chance of ending up worse is too risky. We need more information about μmax. The question is: which information?
=#

# ╔═╡ eb613ead-6aea-41fe-acb7-f679b5a97c35
md"""
## Local Sensitivity Analysis (LSA)

The lab can run the mutant in a pilot reactor under the same conditions, but you can only afford **three measurements** of $X$, all at the same time point. Measurements are expensive, so we want them at the moment where $X$ depends most strongly on $\mu_{max}$. There, a small change in $\mu_{max}$ causes a large change in $X$. The other way around, a measurement of $X$ at that moment pins down $\mu_{max}$ best.

So, we will first perform a normalized local sensitivity analysis to pin down the moment, $X$ is most sensitive to a change in $μmax$.
"""

# ╔═╡ 67f572b0-3485-40bc-8f64-cc58471b38a9
md"""
Write a solution function `mutant_sim` with as argument a vector with the parameter `μmax` of the mutant, and that returns the solution. Based on it, write a function `mutant_sim_X` that only returns the output $X$.
"""

# ╔═╡ 565b6e18-012f-4b38-bf21-bdc683688107
# function mutant_sim(parms)
# 	missing
# 	...
# end
function mutant_sim(parms)
	μmax_val = parms[1]
	oprob = remake(oprob_mut, p = [:μmax => μmax_val])
	osol = solve(oprob, Tsit5(), saveat=dt)
	return osol
end

# ╔═╡ b7f14c48-ba19-4fbc-bded-829a34bd7ac8
# mutant_sim_X(parms) = missing
mutant_sim_X(parms) = mutant_sim(parms)[:X]

# ╔═╡ 73345fe4-6877-476d-8ce1-5a0056c73cc6
md"""
Make the time vector `t_vals` with timestep `dt` equal to $0.5$ and compute the output `X_sim` for the expected value $\mu_{max} = 0.35$, stored in `μmax_val`.
"""

# ╔═╡ 75ab68b2-1e88-421f-a089-4c0b3b72788d
# t_vals = missing
t_vals = tspan[1]:dt:tspan[2]

# ╔═╡ 82ae802b-2a7a-4525-9cc6-dfaa552ba990
# μmax_val = missing
μmax_val = 0.35

# ╔═╡ 8a655d11-0d45-4e7c-bcbf-8ea9659516aa
# X_sim = missing
X_sim = mutant_sim_X([μmax_val])

# ╔═╡ 9a97a7c7-bc3e-434e-9dce-7c859417aac2
md"""
Use `ForwardDiff.jacobian` to compute the sensitivity `sens_X_on_μmax`, and compute the normalized sensitivity `sens_X_on_μmax_rel`.
"""

# ╔═╡ d7d3028a-51c1-41e0-a865-b521ed235548
# sens_X_on_μmax = missing
sens_X_on_μmax = ForwardDiff.jacobian(mutant_sim_X, [μmax_val])[:, 1]

# ╔═╡ cb6047f2-3884-4fbd-96f2-f15cf684bea8
# sens_X_on_μmax_rel = missing
sens_X_on_μmax_rel = sens_X_on_μmax .* μmax_val ./ X_sim

# ╔═╡ f623b4b4-6010-45fa-b9a9-bc796f0fb106
md"""
Identify the time `t_star` at which the biomass is most sensitive to variations in $\mu_{max}$, and plot the normalized sensitivity function together with `t_star`.

!!! note
	You can find the index of the maximum value in a vector with the `argmax()` function. For example `argmax([1, 2, 5, 4])` will give you 3.
"""

# ╔═╡ 9c3fc6d8-8e5e-425d-9caa-56a9e0b03b99
# begin
# 	t_star_idx = missing
# 	t_star = t_vals[t_star_idx]
#
# 	plot(t_vals, abs.(sens_X_on_μmax_rel),
# 	    title="Sensitivity of X w.r.t. μmax", xlabel="Time (hours)",
# 	    ylabel="normalized sensitivity", label=false, linewidth=2)
# 	vline!([t_star], label="t* = $(t_star) h", linestyle=:dash, color=:red)
# end
begin
	t_star_idx = argmax(abs.(sens_X_on_μmax_rel))
	t_star = t_vals[t_star_idx]

	plot(t_vals, abs.(sens_X_on_μmax_rel),
	    title="Sensitivity of X w.r.t. μmax", xlabel="Time (hours)",
	    ylabel="normalized sensitivity", label=false, linewidth=2)
	vline!([t_star], label="t* = $(t_star) h", linestyle=:dash, color=:red)
end

# ╔═╡ b1825a21-e711-4f99-9a4b-51f3fb73ff82
md"""
!!! questions
	- Compare the sensitivity function with the Monte Carlo simulations above. What do you notice?
	- The sensitivity at $t = 100\,h$ is small, while that is exactly the time we care about. Why is it small? Is it a good idea to measure at $t = 100\,h$?
"""

# ╔═╡ 44234af3-6015-49b6-8709-1f3a3cc0f69d
md"""
Answers:
- missing
- missing
"""
#=
- The spread of the Monte Carlo simulations is largest where the sensitivity is largest. To first order, the uncertainty on the output is σ_X(t) ≈ |∂X/∂μmax(t)| σ_μmax: the sensitivity function tells us *where* parameter uncertainty propagates into output uncertainty.
- At the end, the substrate is almost completely consumed, so X ≈ Y(Sin - S) with S small, and μmax hardly matters. A measurement at 100 h therefore tells us very little about μmax. Yet this small dependency still decides whether the mutant beats the current strain, so we need to know μmax precisely, and we get that by measuring at t*.
=#

# ╔═╡ a00afabe-7bc1-4d8b-9b55-2953c7add9f5
md"""
## Measuring campaign

The lab ran the mutant in the pilot reactor and sampled it at your chosen time. We have defined a function `measure_X(t)` in this notebook that returns three hypothetical measurements of $X$ at the given timepoint `t`, with a relative measurement error of $5\%$. Use this function to get the measurements `X_meas` at the optimal time `t_star`.
"""

# ╔═╡ 635ff400-639c-4c80-93e6-e8721a836c16
begin
	μmax_real = 0.37
	real_solution = solve(remake(oprob_mut, p = [:μmax => μmax_real]), Tsit5(), saveat=dt)
	measure_X(t) = real_solution(t, idxs=:X) .* (1 .+ 0.05 .* randn(Xoshiro(round(Int, 10t)), 3))
end;

# ╔═╡ 93cc5a4f-f40a-4d28-abc3-8d6790fd60dc
# X_meas = missing
X_meas = measure_X(t_star)

# ╔═╡ 3e93a001-9401-4a1f-b283-c79fd0c43434
md"""
### Conditioning on the measurements
"""

# ╔═╡ 04ae6cd7-2d59-4bc0-a894-d9d5129e243d
md"""
Based on these measurements we can update what we believed before (the prior) with what the data tell us (the likelihood) into an updated belief (the posterior):

Create a Turing model `mutant_meas(t_meas)` with:
- the same prior for `μmax_dev` as before;
- a prior `Exponential()` for the measurement error `σ_X`;
- the simulated biomass `X_pred` at time `t_meas`, which you can get with `sol(t_meas, idxs=:X)`;
- the three measurements `X_meas`, each normally distributed around `X_pred` with standard deviation `σ_X`;
- the solution as return value.
"""

# ╔═╡ 02096fef-224e-4896-9f8d-3f7aba7334e1
# @model function mutant_meas(t_meas)
# 	μmax_dev ~ missing
# 	σ_X ~ missing
#
# 	sol = missing
# 	X_pred = missing
#
# 	X_meas = zeros(3)
# 	for i in 1:length(X_meas)
# 		X_meas[i] ~ missing
# 	end
# 	return sol
# end
@model function mutant_meas(t_meas)
	μmax_dev ~ LogNormal(log(0.35), 0.15)
	σ_X ~ Exponential()

	sol = solve(remake(oprob_mut, p = [:μmax => μmax_dev]), Tsit5(), saveat=dt)
	X_pred = sol(t_meas, idxs=:X)

	X_meas = zeros(3)
	for i in 1:length(X_meas)
		X_meas[i] ~ Normal(X_pred, σ_X)
	end
	return sol
end

# ╔═╡ 354f0768-efc0-4304-aa04-c6a68f961aa6
md"""
Condition the model on the measurements and store it in `mutant_cond_model`.
"""

# ╔═╡ 11e105de-270c-434a-bfb4-9118dfc9e461
# mutant_cond_model = missing
mutant_cond_model = mutant_meas(t_star) | (X_meas = X_meas,)

# ╔═╡ e46afb98-16de-4d19-86b4-b920b6b929c5
md"""
### MCMC
"""

# ╔═╡ f99ca771-8c14-4e19-a1ae-478a2103716a
# mutant_chain = missing
mutant_chain = sample(mutant_cond_model, NUTS(), 500)

# ╔═╡ e76a7bf7-4d8d-4f2f-8272-b63e49cd4ae9
md"""
Plot the chain to check if it properly converged to the posterior distribution.
"""

# ╔═╡ 50b30e65-72d3-4d8d-a739-ae969e0f50a8
# missing
plot(mutant_chain)

# ╔═╡ a65710f6-70ce-4348-923d-ed2ed2b4a469
md"""
Now, get the sampled posterior solutions `solutions_mutant` and plot them on top of the prior Monte Carlo simulations, together with the current strain and the measurements.
"""

# ╔═╡ 12c9f810-ef5b-49a4-8fdd-6efb1a040440
# solutions_mutant = missing
solutions_mutant = generated_quantities(mutant_cond_model, mutant_chain);

# ╔═╡ e2e11092-4498-437b-8c67-8d77302893fe
md"""
The posterior chance that the mutant is better is:
"""

# ╔═╡ cae4d9c4-3790-423e-8321-d44c0833aa40
begin
	sol_end = [solutions_mutant[i][:X][end] for i in eachindex(solutions_mutant)]
	P_better_posterior = mean(sol_end .> X_current)
	println("The chance the mutant is better is $(P_better_posterior*100)%")
end

# ╔═╡ 49307c36-3518-4118-b453-91e634f92daf
# begin
# 	p2 = plot(title="Prior vs posterior of the mutant", ylabel="X (g/l)", xlabel="Time (hours)")
# 	for missing in missing
# 		missing
# 	end
# 	for missing in missing
# 		missing
# 	end
# 	plot!(osol, idxs=:X, label="Current strain", color=:red, linewidth=3)
# 	scatter!(fill(t_star, 3), X_meas, label="Measurements", color=:black)
# 	p2
# end
begin
	p2 = plot(title="Prior vs posterior of the mutant", ylabel="X (g/l)", xlabel="Time (hours)")
	for solution in solutions
		plot!(solution, idxs=:X, label=false, color=:gray, alpha=0.2, linestyle=:dot)
	end
	for solution_mutant in solutions_mutant
		plot!(solution_mutant, idxs=:X, label=false, color=:orange, alpha=0.4)
	end
	plot!(osol, idxs=:X, label="Current strain", color=:red, linewidth=3)
	scatter!(fill(t_star, 3), X_meas, label="Measurements", color=:black)
	p2
end

# ╔═╡ aa479bcd-2314-47a5-83ba-b0d1a2654fa0
md"""
Create a vector `x2` with the posterior biomass at the end of the simulations, compare its histogram with the prior one, and compute the probability `P_better_post` that the mutant beats the current strain.
"""

# ╔═╡ 2f55a771-0d37-48e5-a9d9-24693bf1d177
# x2 = missing
x2 = [solution_mutant[:X][end] for solution_mutant in solutions_mutant];

# ╔═╡ 18059d0c-cf37-44a1-86a8-859ac66a66f0
# missing
begin
	histogram(x1, normalize=:pdf, alpha=0.5, label="Prior", title="Histogram of X(100 h)", xlabel="X (g/l)", ylabel="density", legend=:topleft)
	histogram!(x2, normalize=:pdf, alpha=0.5, label="Posterior")
	vline!([X_current], linewidth=3, color=:red, label="Current strain")
end

# ╔═╡ 40e73ddf-6c0d-4f83-860f-434d934193ed
# P_better_post = missing
P_better_post = mean(x2 .> X_current)

# ╔═╡ e1c87260-25ee-4d71-8744-ed0d7cb75f6a
md"""
### What if we had measured at the end?

Intuitively, you might have measured $X$ at $t = 100\,h$, since that is the quantity we care about. Repeat the inference with three measurements at $t = 100\,h$, and compute the probability `P_better_100` that the mutant is better.
"""

# ╔═╡ caa5c04d-49af-4576-83a2-605bae1c2494
# X_meas_100 = missing
X_meas_100 = measure_X(100.0)

# ╔═╡ 6a53aff6-29d4-4d41-8236-f621ce519361
# missing
begin
	mutant_cond_model_100 = mutant_meas(100.0) | (X_meas = X_meas_100,)
	mutant_chain_100 = sample(mutant_cond_model_100, NUTS(), 500)
	solutions_mutant_100 = generated_quantities(mutant_cond_model_100, mutant_chain_100)
	x3 = [solution_mutant[:X][end] for solution_mutant in solutions_mutant_100]
end;

# ╔═╡ 9292c593-1f0f-47a6-a2ec-91a8aadbe6e6
# P_better_100 = missing
P_better_100 = mean(x3 .> X_current)

# ╔═╡ b1814ea6-0065-4ac2-a012-505d2bbd8e40
md"""
Compare what we know about $\mu_{max}$ before and after both measuring campaigns:
"""

# ╔═╡ 44cb7615-86b2-45f0-b87a-69fcbfe863bb
begin
	density(vec(chain[:μmax_dev]), label="Prior", linewidth=2, title="What did we learn about μmax?", xlabel="μmax (1/h)", ylabel="density")
	density!(vec(mutant_chain_100[:μmax_dev]), label="Posterior, measured at t = 100 h", linewidth=2)
	density!(vec(mutant_chain[:μmax_dev]), label="Posterior, measured at t* = $(t_star) h", linewidth=2)
end

# ╔═╡ f364fa1f-5c47-4ca6-99df-17e2b735d98b
md"""
!!! questions
	- Compare `P_better_prior`, `P_better_100` and `P_better_post`.
"""

# ╔═╡ 5f9873a5-af70-4ac2-be08-4f249b71eb76
md"""
Answers:
- missing
- missing
"""
#=
- Measuring at t* shrinks the uncertainty on μmax from about 15% to about 3%, while measuring at 100 h only brings it down to about 10%. At t*, X is about ten times more sensitive to μmax than at 100 h, so the same 5% measurement error translates into a much smaller error on μmax. Sensitivity analysis tells you *where* to measure.
=#

# ╔═╡ 3980da6b-7036-44d2-9d8f-3ddc1062d0ad
md"""
## Conclusions

!!! questions
	- What do you advise management? Support your advice with numbers.
	- What did you learn/what did you use the information for from uncertainty propagation, sensitivity analysis and the Bayesian inference.
"""

# ╔═╡ 8021123d-1b19-47bf-9db0-0d7c587bc00e
md"""
Answers:
- missing
- missing
- missing
"""
#=
- Switch to the mutant. With the measurements, μmax ≈ 0.37 1/h and the probability that the mutant produces more biomass after 100 h is practically 1, with about 3-4% more biomass than the current strain.
- Uncertainty propagation turned "we don't know μmax" into a probability for the decision (about 65-70%), and showed that this was not good enough. Sensitivity analysis showed when a measurement is informative (t* ≈ 45 h, not at the end). Bayesian inference combined the prior knowledge with only three measurements into a posterior that is narrow enough to decide, and that posterior was propagated again in the same way.
- E.g. the dose of a drug for a patient with an unknown clearance rate, deciding whether a new fertilizer increases crop yield, choosing when to sample a river to estimate a degradation rate, ...
=#

# ╔═╡ Cell order:
# ╠═1ec7def2-0c43-11ef-0850-ddefa680d2a4
# ╠═2861bc17-7810-461e-8857-a55b6392977a
# ╠═9d8acf40-e635-4dc4-9938-ec63ab68e3bd
# ╠═093679cc-122b-4b0a-84c3-15cf3a6f5c1b
# ╠═0720732e-105a-40c2-bf45-f8fbbc10e5bc
# ╠═28bb60c9-ab56-497f-b039-0d6577501be5
# ╠═0a65afbc-643a-4fe4-b22d-c5fe54e3be1f
# ╠═d1b6739f-1e0c-452d-8763-717dd650b2ec
# ╠═162b330b-02f7-4179-b950-d960eca68206
# ╟─437364b6-a953-428b-a891-e17ffd531656
# ╟─c44e10f3-8452-4a60-b840-f4a7e6ec3238
# ╟─f3065171-5526-476a-bfd7-a2b95bb24d61
# ╟─356822df-bad3-41f9-8556-de04234b0a09
# ╠═c551c33e-95b0-4ea8-85d4-6c68db72b8bc
# ╟─69f8c265-fe91-4e4b-971b-395325d4474e
# ╟─2806654f-7a24-4f8e-8b5a-3129b2431c95
# ╠═5d6a881d-421e-4dea-9668-38b3f3378bf1
# ╟─5802078d-12bc-489e-8566-15f77a11e8ba
# ╠═6afb23f2-3da7-43e0-8a3b-fc4bb5eb9bde
# ╠═1b09c8cb-d3b7-4a15-b6ed-d03019c5f57b
# ╠═9e3c427b-8a56-4787-9b28-43ec3932ddbe
# ╠═62f4049b-617f-4268-beba-b2c634b46642
# ╟─0e4c73d3-1630-463b-994d-073170c1a6c3
# ╠═3cad11d6-c2eb-4865-858a-d407aa25d780
# ╠═40168902-0780-4f92-8eca-23982c80511e
# ╠═4d753fdb-089d-41dc-89ef-45a1733af786
# ╟─ca660e2c-9348-40b1-b92d-415068bea514
# ╠═ec042e61-04b8-4d51-8e7d-35e1646cc5ba
# ╟─9ee22028-7a37-4510-9b61-198d3bd99a73
# ╟─cc9ac43b-448a-4a67-ba0b-8fdf21a4236e
# ╟─f5482531-b23f-472b-968f-20f7a6471778
# ╟─ff21b3cd-64b1-4d09-8eb1-afc617585fc8
# ╠═225eaaff-42a2-4aca-94e1-afcd988de4eb
# ╟─18ffe877-c579-4b4c-80bd-b5d50baf36bb
# ╠═7821dc21-3e6b-4ad1-829e-f29dcd020144
# ╟─3051d39a-fc2f-40e4-a688-4d6e75cfc2d1
# ╠═33df05c9-08d3-4107-9240-931658cda2a4
# ╟─02c5cf65-ebd3-477b-9a3e-8855bd46f55c
# ╠═aa833d92-567f-4971-a360-e60bce424c8e
# ╟─ffc66a94-5435-4ca6-af87-76c85a599906
# ╠═d3afd5b4-034f-415c-a64a-776deed74d5c
# ╠═1d95728e-079e-4e9b-952a-c4f692661813
# ╠═02213469-a233-4759-a18d-29a76a4c13c6
# ╟─97d2065c-6234-4b13-8e4c-6782b8301416
# ╟─d2acbfb9-a7ff-48b0-9d9a-9492ecb9c1bf
# ╟─eb613ead-6aea-41fe-acb7-f679b5a97c35
# ╟─67f572b0-3485-40bc-8f64-cc58471b38a9
# ╠═565b6e18-012f-4b38-bf21-bdc683688107
# ╠═b7f14c48-ba19-4fbc-bded-829a34bd7ac8
# ╟─73345fe4-6877-476d-8ce1-5a0056c73cc6
# ╠═75ab68b2-1e88-421f-a089-4c0b3b72788d
# ╠═82ae802b-2a7a-4525-9cc6-dfaa552ba990
# ╠═8a655d11-0d45-4e7c-bcbf-8ea9659516aa
# ╟─9a97a7c7-bc3e-434e-9dce-7c859417aac2
# ╠═d7d3028a-51c1-41e0-a865-b521ed235548
# ╠═cb6047f2-3884-4fbd-96f2-f15cf684bea8
# ╟─f623b4b4-6010-45fa-b9a9-bc796f0fb106
# ╠═9c3fc6d8-8e5e-425d-9caa-56a9e0b03b99
# ╟─b1825a21-e711-4f99-9a4b-51f3fb73ff82
# ╟─44234af3-6015-49b6-8709-1f3a3cc0f69d
# ╟─a00afabe-7bc1-4d8b-9b55-2953c7add9f5
# ╟─635ff400-639c-4c80-93e6-e8721a836c16
# ╠═93cc5a4f-f40a-4d28-abc3-8d6790fd60dc
# ╟─3e93a001-9401-4a1f-b283-c79fd0c43434
# ╟─04ae6cd7-2d59-4bc0-a894-d9d5129e243d
# ╠═02096fef-224e-4896-9f8d-3f7aba7334e1
# ╟─354f0768-efc0-4304-aa04-c6a68f961aa6
# ╠═11e105de-270c-434a-bfb4-9118dfc9e461
# ╟─e46afb98-16de-4d19-86b4-b920b6b929c5
# ╠═f99ca771-8c14-4e19-a1ae-478a2103716a
# ╟─e76a7bf7-4d8d-4f2f-8272-b63e49cd4ae9
# ╠═50b30e65-72d3-4d8d-a739-ae969e0f50a8
# ╟─a65710f6-70ce-4348-923d-ed2ed2b4a469
# ╠═12c9f810-ef5b-49a4-8fdd-6efb1a040440
# ╟─e2e11092-4498-437b-8c67-8d77302893fe
# ╠═cae4d9c4-3790-423e-8321-d44c0833aa40
# ╠═49307c36-3518-4118-b453-91e634f92daf
# ╟─aa479bcd-2314-47a5-83ba-b0d1a2654fa0
# ╠═2f55a771-0d37-48e5-a9d9-24693bf1d177
# ╠═18059d0c-cf37-44a1-86a8-859ac66a66f0
# ╠═40e73ddf-6c0d-4f83-860f-434d934193ed
# ╟─e1c87260-25ee-4d71-8744-ed0d7cb75f6a
# ╠═caa5c04d-49af-4576-83a2-605bae1c2494
# ╠═6a53aff6-29d4-4d41-8236-f621ce519361
# ╠═9292c593-1f0f-47a6-a2ec-91a8aadbe6e6
# ╟─b1814ea6-0065-4ac2-a012-505d2bbd8e40
# ╟─44cb7615-86b2-45f0-b87a-69fcbfe863bb
# ╟─f364fa1f-5c47-4ca6-99df-17e2b735d98b
# ╟─5f9873a5-af70-4ac2-be08-4f249b71eb76
# ╟─3980da6b-7036-44d2-9d8f-3ddc1062d0ad
# ╟─8021123d-1b19-47bf-9db0-0d7c587bc00e
