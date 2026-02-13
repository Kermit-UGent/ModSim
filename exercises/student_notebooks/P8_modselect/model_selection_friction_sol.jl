### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ dd4af742-7f4a-4e38-9665-620a13b9d341
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ be891cc2-f306-4cc8-97d7-1efa7f97ecd6
using Markdown

# ╔═╡ ef0b5541-3f36-4e79-87eb-0d73ffbd80d4
using InteractiveUtils

# ╔═╡ bbd6918a-1d4b-442d-81a7-e25c8f48c904
using Catalyst, OrdinaryDiffEq

# ╔═╡ 09902aca-cb98-41f4-93f4-a7530c4373ba
using Turing, StatsPlots, StatsPlots.StatsBase

# ╔═╡ 579e2373-4836-4fd5-83b9-c6715984634c
using LinearAlgebra, Optim

# ╔═╡ f9c07cce-bada-11ef-36af-9f9bf85ccc38
using PlutoUI; TableOfContents()

# ╔═╡ eb43c9e6-d47e-414f-987e-3dcc7a7c3124
md"""
# Model selection for air friction terms
"""

# ╔═╡ 66243d50-9a2d-40db-b879-304f182eb46b
md"""
The table below contains data for how long it takes for a badminton projectile to fall a given distance when dropped. The goal is to determine which model for air resistance best describes it, namely:

$v'(t) = g - \frac{F_r(v)}{m}$

where $v$ is the object's speed, $m$ the mass and $F_r(v)$ some function of speed and extra parameters. The goal is also to determine, for such an object falling in air, the model parameters (calibration).
"""

# ╔═╡ 4a5e1692-171c-43fa-8309-b73c8ecd38c5
t = [0, 0.347, 0.47, 0.519, 0.582, 0.65, 0.674, 0.717, 0.766, 0.823, 0.87, 1.031, 1.193, 1.354, 1.501, 1.726, 1.873]

# ╔═╡ 8a71f11b-3dab-4ee3-bc61-dc5327f37eb9
d = [0, 0.61, 1.00, 1.22, 1.52, 1.83, 2.00, 2.13, 2.44, 2.74, 3.00, 4.00, 5.00, 6.00, 7.00, 8.50, 9.50]

# ╔═╡ 19bf959c-b385-44d3-8311-8f2d22e67584
scatter(t, d; xlabel="time (seconds)", ylabel="distance (meters)", legend=:none, title="measurements")

# ╔═╡ cac10b59-aba7-4c0c-a699-02251a349dd4
md"""
# 1. Mathematical model
"""

# ╔═╡ 16b137b5-466a-4206-8d33-de2bf4a8850a
md"""
The general model, obtained from Newton's second law and an opposing resistance force, follows:

$v'(t) = g - \frac{F_r(v)}{m}$

where $F_r$ will only have an influence when the object is falling, thus we assume $F_r=0$ for $v=0$.

"""

# ╔═╡ 1dca1d32-42c3-4636-aa86-66bcd3ea2dcd
md"""
From here, the different models can be developed. For example, for a linear model we have:

$v'(t) = g - \frac{k\,v(t)}{m}$

or

$v'(t) = g - \tilde k\,v(t)$

where the constant $m$ is absorbed in the definition of $k$ and calibrated jointly.
"""

# ╔═╡ 211988a4-e589-4b8f-98e0-c1db92fe5030
md"""
Parameters: $k$, $m$, or simply $k$ if we replace $k/m$ by $\tilde k$ \

Constants: $g$
"""

# ╔═╡ 7daa93a7-7968-4f7f-9e75-d27a412a7362
g = 9.80136 # m s²

# ╔═╡ 017e4b49-dda4-4be7-a6e9-e5bd09218c36
md"""
Integrating twice and solving for $v(0)=0$ and $x(0)=0$ we obtain:

$x(t) = \frac{g}{k^2} \Big(kt + e^{-kt}−1\Big)$

We can use this analytical solution to verify our implementation.
"""

# ╔═╡ 8145a2d1-c771-4530-8683-a057c4a4ad71
md"""
#### No air resistance
"""

# ╔═╡ cc0130e0-6fb1-482e-8b80-6726d4e7a08c
md"""
For $F(v) = 0$ and initial speed $v(0) = 0$, the position of the projectile would be $x(t) = gt^2/2$, according to classical laws and assuming $x(0)=0$. There are no parameters to estimate here, since $g$ is known. We can simply assess the fit to the data as a least squares problem.
"""

# ╔═╡ d8a642cc-9f5c-4696-8314-20a998498eee
md"""
Computing the sum of squares 

$SSR = \sum_{j=1}^N \big(x(t_j) - x_j\big)^2$

where $(t_j, x_j)$ denotes the $j$th time-distance pair from the table above, allows us to assess the fit.
"""

# ╔═╡ 35b9f77a-e939-42ef-9bfa-da2e30ffe659
md"""
## 1.1 Linear air resistance
"""

# ╔═╡ d2e7a649-2196-4f11-a790-96ce43396821
md"""
Let's implement the model above in Catalyst, where we 'create' the properties of the system:
"""

# ╔═╡ c27d4852-4852-41ba-a1b6-759b8efde91e
model_linear = @reaction_network begin
	@parameters k g
	g - k*v, 	∅ → v
	v, 			∅ → x
end

# ╔═╡ d4b00278-6933-45d8-bdfd-8d98682fd838
md"""
We verify the differential equations of the model:
"""

# ╔═╡ a324dda4-c785-4ad0-8f2f-0b1e5ae5445b
convert(ODESystem, model_linear)

# ╔═╡ f01a8cd4-87de-44a1-ba58-d1dd78e4f5f3
md"""
We will solve the model and plot the solution for different values of $k$ using a slider. See below.
"""

# ╔═╡ a75ce9c0-4e06-4950-a4f9-3baa5027a534
u0 = [:v => 0.0, :x => 0.0]

# ╔═╡ 4fe4f7bb-aee5-44d7-8503-916e3829ce65
tspan = (0.0, t[end])

# ╔═╡ 0fb380cf-ba93-4e01-918a-c052b725fa0e
@bind k_linear Slider(0:0.1:1.5, default=0, show_value=true)

# ╔═╡ bd6a3907-bc0f-4cc0-b35c-e33af1ce9a5a
prob_linear = ODEProblem(model_linear, u0, tspan, [:k => k_linear, :g => g])

# ╔═╡ 047bced8-8cdd-47c9-a683-4bb5adb208fd
sol_linear = solve(prob_linear);

# ╔═╡ 8c4f63ba-5424-44b0-aa85-3e32e7793fb0
begin
	td = 0:0.1:t[end] 	# vector to plot analytical solutions
	
	plot(sol_linear, label=["v, velocity (model)" "x, distance (model)"], lw=2)

	scatter!(t, d, label="d, measured distance", markersize=5)
	
	plot!(td, g*td.^2/2, ls=:dash, label="x, k=0 (analytical)", lw=2, c=:grey)
	plot!(td, @. g/k_linear^2*(k_linear*td + exp(-k_linear*td) - 1); ls=:dashdot, label="x, k>0 (analytical)", lw=1, c=:black)

	ylims!(0, 10)
	plot!(legend_position=:bottomright)
	title!("Predicted position vs data: linear model")
end

# ╔═╡ 32545342-1704-412f-a898-a3d1f5d3593d
md"""
!!! question
	How do we interpret this simulation? Can we decide on the approppriateness of such a model?
"""

# ╔═╡ a490335d-0142-48a8-b3be-656bef552139
md"""
Answer:
- For $k=0$, the model is equivalent to no air resistance. The velocity increases indefinitely, and thus the distance of the projectile increases quadratically. This doesn't correspond to the data.
- For increasing values of $k>0$, the modeled distance fits much better the measurements.
"""

# ╔═╡ 85471dd1-df76-4718-a7ee-d501bb8a2030
md"""
## 1.2 Calibration of the model
"""

# ╔═╡ 4bde3bb3-5dcd-42d4-ae48-b922e99a04bc
md"""
We create a Turing model with both parameters $k$ and $g$, where only the first one is calibrated:
"""

# ╔═╡ fa3383e4-f5fa-4d47-b21b-684c6aacbd46
@model function linear(t)
	σ_x ~ InverseGamma()
	k ~ LogNormal()
	params_linear = [:k => k, :g => g]
	oprob_linear = ODEProblem(model_linear, u0, tspan, params_linear)
	osol_linear = solve(oprob_linear, Tsit5(), saveat=t)
	x ~ MvNormal(osol_linear[:x], σ_x^2 * I)
end

# ╔═╡ 1f36c90e-ec20-4bdb-8bd0-5f44f03a1748
mod_linear_cond = linear(t) | (x = d,)

# ╔═╡ c4de5cb6-4b84-46b0-acc5-01aea794de3e
results_linear = optimize(mod_linear_cond, MLE(), NelderMead())

# ╔═╡ 3153b018-57f8-4d7a-bfb9-c4ca2886363a
coeftable(results_linear)

# ╔═╡ 7551a678-2873-4d10-bef7-81fd0e275b9f
k_opt_linear = coef(results_linear)[:k]

# ╔═╡ 621a2943-bb80-414e-babf-b5b2eac7269a
prob_opt_linear = remake(prob_linear, p=[:k => k_opt_linear])

# ╔═╡ cd8d78e3-18d3-4af1-ac69-71d007d41b7e
sol_opt_linear = solve(prob_opt_linear);

# ╔═╡ 0c9eb9d0-cda4-4f71-ba80-7f16680201d2
begin
	plot()
	plot!(sol_linear; idxs=:x, ls=:dash, label="candidate model, k=$k_linear", lw=2)
	plot!(sol_opt_linear; idxs=:x, label="calibrated model, k=$k_opt_linear", lw=3)
	scatter!(t, d, label="measured data", markersize=5)
	xlabel!("time (seconds)")
	ylabel!("distance (meters)")
	title!("Model prediction vs. measured position")
end

# ╔═╡ 1d5a95fb-2aeb-4a1e-808b-74b086970556
md"""
!!! note
	It seems we have found a good model for the air resistance. *But is this the best model of them all?*
"""

# ╔═╡ 3d2e0801-b075-4f91-a9db-55cdda05477d
md"""
# 2. Model selection criteria
"""

# ╔═╡ 8cc42a7b-ba84-43d0-aafa-344712e882c6
md"""
## 2.1 Akaike information criterion
"""

# ╔═╡ b78015f8-5edd-4a64-9e9b-b6856c525767
md"""
We implement here again a function to calculate the AIC criterion based on the calibration results.
"""

# ╔═╡ f5bdb7b0-5a32-42cf-9461-75a8b5b513c2
function AIC(results, measurements)
	L = results.lp
	k = length(coef(results))
	n = length(measurements)
	
	return round(2k - 2L; digits=3)		# L = log-likelihood
end

# ╔═╡ 1d8a59eb-67b2-4242-9ab3-96fec2856622
md"""
## 2.2 Bayesian information criterion
"""

# ╔═╡ d6b22a0e-f2d1-4869-8405-46d163686061
md"""
The same can be done for the BIC criterion.
"""

# ╔═╡ fb75ff46-58d4-4bb9-a7d9-36b66124ea90
function BIC(results, measurements)
	L = results.lp
	k = length(coef(results))
	n = length(measurements)
	
	return round(k*log(n) - 2L; digits=3) 	# L = log-likelihood
end

# ╔═╡ a6e960f6-fe3e-4f0e-9ee0-358010da6e39
md"""
## 2.3 Least squares problem
"""

# ╔═╡ abd8fec6-b1d9-467b-a543-e479506b134f
md"""
The least squares alternative form of both AIC and BIC can also be used.
"""

# ╔═╡ 9db0f7d9-1d0c-40f7-89ef-a4786f8978ea
function AIC_LS(SSR, n, k) 
	if n > 40
		return 2k + n*log(SSR/n) 
	else
		return 2k + n*log(SSR/n) + 2k*(k+1)/(n-k-1)
	end
end

# ╔═╡ a78f8bbb-52c5-45db-90c2-60837eefa6be
function BIC_LS(SSR, n, k) 
	return  k*log(n) + n*log(SSR/n)
end

# ╔═╡ 5e6bc47b-bdf3-4b49-babb-6f3bf795fb7a
md"""
We can thus obtain the squared sum of residuals from the calibrated model prediction and the data:
"""

# ╔═╡ 4d5cd521-bffd-4f30-8a23-a04b3c5b9624
function SSR(y_pred, y_data)
	return sum((y_pred - y_data).^2) 	# squared sum of residuals
end

# ╔═╡ 2393409e-b5ec-4ce9-a818-c2a89bf50885
md"""
## 2.4 Posterior probabilities
"""

# ╔═╡ 5c500c1a-a15a-49bf-bd74-53985a484f19
md"""
We can use the AIC to compute the posterior probabilities of the different candidate models:

$P(M_i|D) \propto \exp(-AIC(M_i)/2)$
"""

# ╔═╡ e56041b0-e722-40ae-8310-64d7069446a7
md"""
The following function will use the supplied AIC of several models to compute the normalized posterior probability that the model is the "true model", explaining the considered data set: 
"""

# ╔═╡ ca46ccb7-523c-4540-8daa-f6ee5f134f8f
function posterior(AICs) 	# AICs vector of AIC values
	AICmin = minimum(AICs)
	posterior = zeros(length(AICs))
	
	for i in eachindex(1:length(AICs))
		posterior[i] = exp((AICmin-AICs[i])/2)
	end
	
	return round.(posterior/sum(posterior); digits=3)	# normalized sum
end

# ╔═╡ c8f0cc9b-ff34-4baa-8234-b828df4fe894
md"""
# 3. Resistance models
"""

# ╔═╡ d0c9e1f1-ac1a-44c9-aaec-6c1cfd75f80d
md"""
We will examine several different models for the force of air resistance on the projectile as it falls. The models to be examined are:

(a) No air resistance. \
(b) Air resistance proportional to speed. \
(c) Air resistance proportional to the square of speed. \
(d) Air resistance proportional to a more general quadratic function of speed. \
(e) Air resistance proportional to an rth power of speed.

We’ll formulate and solve an ODE model for each possibility and compare how well each model ﬁts the data. To decide which model is best we’ll invoke several **model selection criteria** (AIC, BIC, etc.).
"""

# ╔═╡ 57f02857-1032-4ce2-a4fb-1cbeb6b6ac3d
md"""
$v'(t) = g - F\big(v(t)\big)$
"""

# ╔═╡ f2f496ab-4354-4d80-888e-89367efa00b0
md"""
The goal is to determine which type of function for $F$ best models air resistance in this situation by making use of the data in the table above. The choices for air resistance models listed above lead to the following possibilities:

| Model   | $F(v)$|
|:---|:---|
|(a) No air resistance 				|$F(v) = 0$   |
|(b) Linear air resistance 			|$F(v) = kv$  for some positive constant $k$ |
|(c) Pure quadratic resistance 		|$F(v) = kv^2$ for some positive constant $k$ |
|(d) General quadratic resistance 	|$F(v) = k_1v+k_2v^2$ for positive constants $k_1$ and $k_2$  |
|(e) General power law resistance 	|$F(v) = kv^r$ for positive constants with $k$ and $r$ |

"""

# ╔═╡ 06b3ede3-93b4-4784-8298-272f8646b1d5
md"""
We have already implemented the linear model and, by extension, the no resistance model ($k=0$).
"""

# ╔═╡ d6bba283-75ed-40a4-a7ef-48f74c305ff2
md"""
For the remaining models, we will need to estimate the different parameters if we want to assess the quality of the fit. This will already provide us some results for the choice of best resistance model.
"""

# ╔═╡ 5d829be3-9bc9-4889-82d5-dc3a755b6435
md"""
## 3.1 No air resistance
"""

# ╔═╡ 712f8974-ba9f-46b6-b569-fb9741567292
model_nr = @reaction_network begin
	@parameters g
	g, ∅ → v
	v, ∅ → x
end

# ╔═╡ 8ca0791a-3da5-4b08-bf80-8221e4646297
sys_nr = convert(ODESystem, model_nr)

# ╔═╡ 6652e22c-2600-4ab6-b84a-2205582de261
prob_nr = ODEProblem(model_nr, u0, tspan, [:g => g])

# ╔═╡ f72e9d76-76d4-4da6-8c9c-c4fe433d4611
sol_nr = solve(prob_nr);

# ╔═╡ 90b9a18d-eaa8-4674-a092-82e6cb3ae8b1
begin
	plot()
	plot!(sol_nr; idxs=[:x], label="no air resistance", lw=3)
	plot!(td, g*td.^2/2, linestyle=:dash, label="analytical solution")
	scatter!(t, d, label="measured distance", markersize=5)
	xlabel!("time (seconds)")
	ylabel!("distance (meters)")
	title!("No resistance model")
end

# ╔═╡ 381df36b-1766-4be3-aa32-01d995ab72c0
md"""
Let's now estimate the distance assuming a noisy error and no calibrated parameters.
"""

# ╔═╡ 965ebc6f-b440-42da-8aa7-00b22e6d0a3f
@model function no_resistance(t)
	σ_x ~ InverseGamma()
	params_nr = [:g => g]
	oprob_nr = ODEProblem(model_nr, u0, tspan, params_nr)
	osol_nr = solve(oprob_nr, Tsit5(), saveat=t)
	x ~ MvNormal(osol_nr[:x], σ_x^2 * I)
end

# ╔═╡ d7514616-736b-4cc3-a562-c23fb3d6924d
model_nr_cond = no_resistance(t) | (x = d,)

# ╔═╡ 4b953047-0f78-4c06-833f-3302a61882fc
results_nr = optimize(model_nr_cond, MLE(), NelderMead())

# ╔═╡ 14dca280-a67f-4760-ac3a-dd7d88155754
md"""
We can obtain the log-probability or likelihood `L` of this model explaining the data from the results:
"""

# ╔═╡ d451a113-249f-4a52-b753-45170c84a593
L_nr = results_nr.lp

# ╔═╡ 9fe87e80-709e-41dc-b2d2-0569c34aeb5d
md"""
If we obtain the maximum likelihood for the other candidate models, we can use it to compare them.
"""

# ╔═╡ a863f1ae-a06c-4829-9fb5-54c3a2f709aa
md"""
We can also calculate the AIC. For the no resistance model, we include the model error as parameter:
"""

# ╔═╡ eea8f1c3-a852-47d8-af77-588fa6499e84
k_nr = length(coef(results_nr)) 	# k here is no. of parameters

# ╔═╡ 4150f328-cf18-4e07-ba35-49f9db3d54fb
md"""
We calculate the AIC and BIC below:
"""

# ╔═╡ 79e107ca-cc86-494f-9f87-f594d9bb874c
AIC_nr = AIC(results_nr, d)

# ╔═╡ ee1909ff-7334-4874-baf6-1f5a4d5e897d
BIC_nr = BIC(results_nr, d)

# ╔═╡ 295db3bc-2007-4d18-928d-e7f67815c538
md"""
For the sum of square residuals, we need to make sure that model predictions are created at the times for which we have measurements. This is simply enforced by creating simulations at times `t`:
"""

# ╔═╡ 76159175-f38f-4979-9132-dfc4b12d3b30
sol_opt_nr = solve(prob_nr, Tsit5(), saveat=t);

# ╔═╡ 7cb8b9ac-bb09-424d-a661-12860446a03c
SSR_nr = SSR(sol_opt_nr[:x], d)

# ╔═╡ 8961def3-e608-4aa4-842d-4a4565aa57f2
n = length(d)

# ╔═╡ 9c083fda-ff0c-438e-abab-32e9845ee299
AIC_LS_nr = AIC_LS(SSR_nr, n, k_nr)

# ╔═╡ b2dc3596-b7b5-4791-8251-83ca8070afd8
md"""
## 3.2 Linear air resistance
"""

# ╔═╡ 2d4008e9-45fd-4bad-941d-51acb6c37a36
results_linear

# ╔═╡ 558e7a9b-85fd-44e0-8e36-5c32e3476964
coeftable(results_linear)

# ╔═╡ 3601dba1-43f9-4f50-b1fe-37a4f26985a8
k_lin = length(coef(results_linear))

# ╔═╡ 99a6cb5b-fe85-4476-abdc-471eadf6add6
L_linear = results_linear.lp

# ╔═╡ b55eaf74-9e69-4ab7-925c-db973cd0b667
AIC_linear = AIC(results_linear, d)

# ╔═╡ a85d7e97-9e0a-4e5e-af87-5816f2b0e44a
BIC_linear = BIC(results_linear, d)

# ╔═╡ 43bab8a0-fcb4-496a-8a0b-34beb3dd1766
sol_linear_t = solve(prob_opt_linear, Tsit5(), saveat=t);

# ╔═╡ 66d87f17-36c5-476f-8696-e69ae55091b4
SSR_linear = SSR(sol_linear_t[:x], d)

# ╔═╡ c59be115-75b8-4432-8cd0-61c7c4449ac0
AIC_LS_linear = AIC_LS(SSR_linear, n, k_linear)

# ╔═╡ a3c53584-7f2e-4324-ac45-ef3183a40c47
md"""
## 3.3 Pure quadratic resistance
"""

# ╔═╡ 17113056-e45f-4977-a022-22de535b1a49
model_quad = @reaction_network begin
	@parameters k g
	g - k*v^2, 	∅ → v
	v, 			∅ → x
end

# ╔═╡ a0227fba-9aca-46b1-aa78-c56742890e74
convert(ODESystem, model_quad)

# ╔═╡ 20d5a4d5-3d17-488a-bda6-62a180e459be
prob_quad = ODEProblem(model_quad, u0, tspan, [:k => 1, :g => g])

# ╔═╡ 250bd1af-97d6-4fd2-bd5f-2b10f3e3e212
sol_quad = solve(prob_quad);

# ╔═╡ b0908503-1523-4a35-bc2f-5087650e537e
@model function quadratic1(t)
	σ_x ~ InverseGamma()
	k ~ LogNormal()
	params_quad = [:k => k, :g => g]
	oprob_quad = ODEProblem(model_quad, u0, tspan, params_quad)
	osol_quad = solve(oprob_quad, Tsit5(), saveat=t)
	x ~ MvNormal(osol_quad[:x], σ_x^2 * I)
end

# ╔═╡ 8d95ddb9-bd58-4aea-9c8f-28e2f2119be2
quad_cond1 = quadratic1(t) | (x = d,)

# ╔═╡ 88c228fe-2aca-4058-803f-31d460ad1915
results_quad1 = optimize(quad_cond1, MLE(), NelderMead())

# ╔═╡ fec66ee6-88f5-4491-8b80-d5e12d4fec01
coeftable(results_quad1)

# ╔═╡ 6e804a51-3a8c-45c0-b078-33af8ca96040
k_opt_quad1 = coef(results_quad1)[:k]

# ╔═╡ 85c68c99-35a8-49e3-bc60-6138e45637fc
model_quad_opt1 = remake(prob_quad, p=[:k => k_opt_quad1])

# ╔═╡ 2a588744-d5c9-4e95-86e6-847fc9397663
sol_opt_quad1 = solve(model_quad_opt1);

# ╔═╡ 70369ddb-39b0-4ca0-9524-a12990e3478a
begin
	plot()
	plot!(sol_quad; idxs=:x, label="initial value, k=1", lw=2, ls=:dash)
	plot!(sol_opt_quad1; idxs=:x, label="calibrated model, k=$k_opt_quad1", lw=3)
	scatter!(t, d, label="measured distance", markersize=5)
	xlabel!("time (seconds)")
	ylabel!("distance (meters)")
	title!("Pure quadratic model")
end

# ╔═╡ d04ec167-203f-413b-aea5-a1a5cbdc70a2
k_quad1 = length(coef(results_quad1))

# ╔═╡ 6ce08886-9271-4628-9afb-5f7f6778fef8
L_quad1 = results_quad1.lp

# ╔═╡ 3cd13eef-7954-4245-ac76-74eb4cf57b2c
AIC_quad1 = AIC(results_quad1, d)

# ╔═╡ 263d4146-e574-4576-8db4-77aa799d7c1a
BIC_quad1 = BIC(results_quad1, d)

# ╔═╡ 32a9b18f-ed5b-498d-bc49-129f2427616b
sol_quad1 = solve(model_quad_opt1, Tsit5(), saveat=t);

# ╔═╡ 39a435b3-7e9c-4dd5-9027-a751d4b20167
SSR_quad1 = SSR(sol_quad1[:x], d)

# ╔═╡ 753a58ec-4a82-480b-bcd3-50448ae09503
AIC_LS_quad1 = AIC_LS(SSR_quad1, n, k_quad1)

# ╔═╡ 9e617891-4a1c-4e56-a7b2-177168d572e0
BIC_LS_quad1 = BIC_LS(SSR_quad1, n, k_quad1)

# ╔═╡ 42e9c861-a25a-4efc-bc6f-8d63c347c082
md"""
!!! question
	Is this the best resistance model that can fit the data? Can we find a better candidate model?
"""

# ╔═╡ 5d6f4d9a-fbe1-4b03-85de-d02d6f505c63
md"""
## 3.4 General quadratic resistance
"""

# ╔═╡ e387491c-a191-47cf-b5c4-8ae6ba4f561a
model_quad2 = @reaction_network begin
	@parameters k₁ k₂ g
	g - (k₁*v + k₂*v^2), 	∅ → v
	v, 						∅ → x
end

# ╔═╡ 990cd6d9-ef7e-4a62-83b8-7724515dfdcc
parameters(model_quad2)

# ╔═╡ 9c8f06fa-946a-4d2c-af3d-593a58bd0f92
convert(ODESystem, model_quad2)

# ╔═╡ 4e5f5f05-7f1a-41cb-a7e9-431f53e97bb5
prob_quad2 = ODEProblem(model_quad2, u0, tspan, [:k₁ => 1.0, :k₂ => 1.0, :g => g])

# ╔═╡ 41ea7cc3-fe7e-4114-8725-1180b8eea871
sol_quad2 = solve(prob_quad2);

# ╔═╡ 926154ef-34d7-452a-8627-bfb1c050a1b1
begin
	scatter(t, d, label="measured data", markersize=5)
	plot!(sol_quad2, idxs=:x, label="quadratic model", lw=3)
end

# ╔═╡ ba0a576e-4b9b-49e1-bfbc-33955dcec808
@model function quadratic2(t)
	σ_x ~ InverseGamma()
	k₁ ~ LogNormal()
	k₂ ~ LogNormal()
	params_quad2 = [:k₁ => k₁, :k₂ => k₂, :g => g]
	oprob_quad2 = ODEProblem(model_quad2, u0, tspan, params_quad2)
	osol_quad2 = solve(oprob_quad2, Tsit5(), saveat=t)
	x ~ MvNormal(osol_quad2[:x], σ_x^2 * I)
end

# ╔═╡ edd1c981-56f6-44d8-a805-0b85c963b39b
model_quad_cond2 = quadratic2(t) | (x = d,)

# ╔═╡ bc173213-4ece-4358-9d70-f1f9f7edea8e
results_quad2 = optimize(model_quad_cond2, MLE(), NelderMead())

# ╔═╡ a94febb7-163e-4f6e-91e5-1a1871e8ea39
coeftable(results_quad2)

# ╔═╡ ec661503-a2e6-4a3d-a244-a4cafebcff60
md"""
!!! question
	What conclusion can you draw from the calibration? Is the generic quadratic model any better?
"""

# ╔═╡ f8b9f59c-71ef-426c-bf25-880d6660e9aa
md"""
- Answer: 
"""

# ╔═╡ 2e8f3065-343f-482c-b3d3-0ace5337343d
k₁_opt = coef(results_quad2)[:k₁]

# ╔═╡ c73544a6-3b24-4f06-8b00-0e0e5059e5ee
k₂_opt = coef(results_quad2)[:k₂]

# ╔═╡ beb474e9-0694-457e-861a-e8ec9ddf865b
opt_quad2 = remake(prob_quad2, p=[:k₁ => k₁_opt, :k₂ => k₂_opt])

# ╔═╡ b032789f-9f7a-45a3-a3fd-9186598bab28
sol_opt_quad2 = solve(opt_quad2);

# ╔═╡ e0690c90-01f7-4e8d-9faf-6f9149a4031d
begin
	plot()
	plot!(sol_quad2, idxs=:x, label="initial values", ls=:dash, lw=2)
	plot!(sol_opt_quad2, idxs=[:x], label="optimal values", lw=2)
	scatter!(t, d, label="measured distance")
	xlabel!("time (seconds)")
	ylabel!("distance (meters)")
	title!("General quadratic model")
end

# ╔═╡ 8a6c6123-36a2-4769-8f22-876eeccf3963
k_quad2 = length(coef(results_quad2))

# ╔═╡ 65aef305-fdaa-4262-8bc5-a22f5a0b42ca
L_quad2 = results_quad2.lp

# ╔═╡ 4dcbb097-8eda-4d34-8d9f-88ff185ef108
AIC_quad2 = AIC(results_quad2, d)

# ╔═╡ a6e1fc57-2750-4765-925e-d321b23866e1
BIC_quad2 = BIC(results_quad2, d)

# ╔═╡ 50213066-8eb4-4d89-a660-da71911104cc
sol_quad2_t = solve(opt_quad2, Tsit5(), saveat=t);

# ╔═╡ a75c3dc3-03bf-4ed4-8fe8-59a63a022030
SSR_quad2 = SSR(sol_quad2_t[:x], d)

# ╔═╡ 6ec4751a-dd24-4da9-a864-1bde6aa5632d
AIC_LS_quad2 = AIC_LS(SSR_quad2, n, k_quad2)

# ╔═╡ c63dcc08-911f-49c5-9316-af07bc8c36b4
BIC_LS_quad2 = BIC_LS(SSR_quad2, n, k_quad2)

# ╔═╡ b88b0d16-a178-4395-a54f-836d4090ccec
md"""
## 3.5 General power law resistance
"""

# ╔═╡ 763b603a-f81b-4d9e-a6eb-7dc345e3e2df
model_power = @reaction_network begin
	@parameters r k g
	g - k*v^r, 	∅ → v
	v, 			∅ → x
end

# ╔═╡ 024b9f4e-600d-4d68-962d-a66c9b5a2f23
parameters(model_power)

# ╔═╡ 343945bd-203d-4a9b-b46d-2a12b46ae32c
prob_power = ODEProblem(model_power, u0, tspan, [:r => 1.0, :k => 1.0, :g => g])

# ╔═╡ 44a7a892-336f-4e52-a66d-230afbb75a7b
sol_power = solve(prob_power);

# ╔═╡ 4b53805c-c9ed-4708-b8a1-4922ce6d5a5d
@model function power(t)
	σ_x ~ InverseGamma()
	k ~ LogNormal()
	r ~ LogNormal()
	params_power = [:k => k, :r => r, :g => g]
	oprob_power = ODEProblem(model_power, u0, tspan, params_power)
	osol_power = solve(oprob_power, Tsit5(), saveat=t)
	x ~ MvNormal(osol_power[:x], σ_x^2 * I)
end

# ╔═╡ 4e579e09-2476-432d-8eb4-7a1080376998
model_power_cond = power(t) | (x = d,)

# ╔═╡ 02f2f82c-4c69-4fc3-90b1-ea0a8796fe12
results_power = optimize(model_power_cond, MLE(), NelderMead())

# ╔═╡ 1943119c-3007-4495-8148-1c346f052a72
coeftable(results_power)

# ╔═╡ d178630c-289d-4917-a182-8f95409d14c3
k_opt_power = coef(results_power)[:k]

# ╔═╡ b42769d2-16d2-449d-8925-fd5b71ddced3
r_opt_power = coef(results_power)[:r]

# ╔═╡ 06b60440-696e-4067-9772-2b91f4be9cbd
opt_power = remake(prob_power, p=[:k => k_opt_power, :r => r_opt_power])

# ╔═╡ 54e63222-6434-4f63-81fa-ab6d5d8fd7da
sol_opt_power = solve(opt_power);

# ╔═╡ 7ea9a771-d692-4ce1-88d8-0a31047dca3d
begin
	plot()
	plot!(sol_power, idxs=:x, label="initial values", ls=:dash, lw=2)
	plot!(sol_opt_power; idxs=:x, label="optimal values", lw=3)
	scatter!(t, d, label="measured distance", markersize=5)
	xlabel!("time (seconds)")
	ylabel!("distance (meters)")
	title!("Power law model")
end

# ╔═╡ f2f65b41-edb4-411f-aa4c-cef53e5aaefe
sol_power_t = solve(opt_power, Tsit5(), saveat=t);

# ╔═╡ 6ee18e30-996b-4e30-b0ec-39cf859d5d5a
k_power = length(coef(results_power))

# ╔═╡ 8740b8b8-06a9-420b-8526-2adf5f8219a3
L_power = results_power.lp

# ╔═╡ e49c777a-a2c3-4725-a9ba-abb7d8536d2e
AIC_power = AIC(results_power, d)

# ╔═╡ 8fe403f4-2f3b-482b-91b6-e49672f98d23
BIC_power = BIC(results_power, d)

# ╔═╡ 88fb0c51-9005-44de-bb15-660ddb9d3165
SSR_power = SSR(sol_power_t[:x], d)

# ╔═╡ bbd11d02-29ca-4b4d-b08f-2c4b460fd5c5
AIC_LS_power = AIC_LS(SSR_power, n, k_power)

# ╔═╡ e4c1292a-be82-466f-a586-430f7d34f9bb
BIC_LS_power = BIC_LS(SSR_power, n, k_power)

# ╔═╡ 3a1ce3ee-864b-4144-a9da-3d66bd1870e9
md"""
# 4. Model comparison
"""

# ╔═╡ da1dc114-8529-4ef8-b3af-2acb066fa02f
md"""
We can compare now all candidate models. Since we have found all calibrated parameters, we can plot all model predictions against the data.
"""

# ╔═╡ e84a0ef5-a6e9-4c66-a44c-81e7ef70351d
begin
	# k = k_lin_opt
	plot()
	plot!(sol_nr, idxs=:x, label="no air resistance", lw=2)
	plot!(sol_opt_linear, idxs=:x, label="linear", lw=2, ls=:dot)
	plot!(sol_opt_quad1, idxs=:x, label="quadratic 1", lw=2, ls=:dashdot)
	plot!(sol_opt_quad2, idxs=:x, label="quadratic 2", lw=2, ls=:dashdot)
	plot!(sol_opt_power, idxs=:x, label="power law", lw=2, ls=:dashdot)
	scatter!(t, d, label="measured distance", markersize=5)
	ylims!(0, 10)
	xlabel!("time (seconds)")
	ylabel!("distance (meters)")
	title!("Comparison all resistance models")
end

# ╔═╡ 9b8392a5-fc9f-433f-9b23-dab2dd472021
md"""
Both the no resistance and linear models fit worse. The more complex models are indistinguisable.
"""

# ╔═╡ 8d1c2906-f7bf-4557-81af-43b8d9b2a260
md"""
Since we cannot tell these very similar complex models apart, we will use the model selection criteria to rank them. The following graph should summarize everything once they have all been calculated.
"""

# ╔═╡ 520e22a5-e08f-4718-bf2f-f1e7e627a417
AIC_nr, AIC_linear, AIC_quad1, AIC_quad2, AIC_power

# ╔═╡ 8fd5f30a-ef31-442f-9124-8ca86a3b9b56
plot(
	bar(1:3, [AIC_quad1, AIC_quad2, AIC_power], title="AIC", ylims=(-78, -70)), 
	bar(1:3, [BIC_quad1, BIC_quad2, BIC_power], title="BIC", ylims=(-78, -70)), 
	bar(1:3, [L_quad1, L_quad2, L_power], title="Log-probability", ylims=(38, 40)), 
	bar(1:3, [k_quad1, k_quad2, k_power], title="k", ylims=(0, 4)), 
	xticks=(1:3, ["Quad (1)", "Quad (2)", "Power"]),
	legend=:none
)

# ╔═╡ 5e8a4a93-7316-4b61-b1e0-bce46597f614
md"""
!!! question
	Draw your conclusions for the best model based on the different criteria. What is your ranking? E.g.: for AIC, quad (1) = quad (2) > power; for Log-p, quad(1) = quad(2) ~ power, etc.
"""

# ╔═╡ dde10d51-3d75-4152-a7da-c05b1affaa6a
md"""
- Answer: 
"""

# ╔═╡ 8becb816-b96c-4df5-8af0-2464c7188d34
md"""
Since we have calculated the AIC for all models, we can estimate the model posterior probabilities:
"""

# ╔═╡ 9bb200cc-09db-4f19-bf29-5b4a1cc78574
posteriors = posterior([AIC_nr, AIC_linear, AIC_quad1, AIC_quad2, AIC_power])

# ╔═╡ 28932b7c-35b0-4e8e-91f2-4dfd5a41407e
md"""
The following table summarizes the model selection criteria for all models. Analyze the ranking.
"""

# ╔═╡ 9dd416bc-eead-4a95-820c-6f0cdf532243
md"""

| Model | Resistance term | k | SSR | log(L) | AIC | BIC | $P(M_i\|D)$ |
|:---|:---|:---|:---|:---|:---|:---|:---|
|(a) No resistance |$F(v) = 0$ | $k_nr | $(round(SSR_nr;digits=3)) | $(round(L_nr;digits=3)) | $AIC_nr | $BIC_nr | $(posteriors[1])
|(b) Linear |$F(v) = kv$ | $k_lin 	| $(round(SSR_linear;digits=3)) | $(round(L_linear;digits=3))	 | $AIC_linear | $BIC_linear | $(posteriors[2])
|(c) Pure quadratic |$F(v) = kv^2$ 	| $k_quad1 	| $(round(SSR_quad1;digits=3)) | $(round(L_quad1;digits=3)) | $AIC_quad1 | $BIC_quad1 | $(posteriors[3])
|(d) Gen. quadratic 	|$F(v) = k_1v+k_2v^2$ 	| $k_quad2  | $(round(SSR_quad2;digits=3)) | $(round(L_quad2;digits=3)) | $AIC_quad2 | $BIC_quad2 | $(posteriors[4])
|(e) Gen. power law 	|$F(v) = kv^r$ 	| $k_power 	| $(round(SSR_power;digits=3)) | $(round(L_power;digits=3))	 | $AIC_power | $BIC_power | $(posteriors[5])
"""

# ╔═╡ 077a5018-538c-4c33-b703-f7258849f2f8
md"""
!!! question
	Write your final ranking for the candidate models, e.g. linear > quadratic > ...
"""

# ╔═╡ d00abaf7-10de-4f7a-9bb9-a476bf4b8351
md"""
- Answer: 
"""

# ╔═╡ a4156697-08c0-418f-95cd-bffe4295eb82
md"""
!!! question
	What candidate model would you select and why? Why would you choose it over the second one?
"""

# ╔═╡ 07d655d8-91ba-49d4-8e1d-d738826b8099
md"""
- Answer: 
"""

# ╔═╡ b3fa82ec-6686-40bf-8709-4d281dcb518f
md"""
!!! question
	How would you combine the predictions of the different models using the posterior probability?
"""

# ╔═╡ 95aa37fe-3cb1-4a9a-952a-2bb2ea815a1b
md"""
- Answer: 
"""

# ╔═╡ Cell order:
# ╠═be891cc2-f306-4cc8-97d7-1efa7f97ecd6
# ╠═ef0b5541-3f36-4e79-87eb-0d73ffbd80d4
# ╠═dd4af742-7f4a-4e38-9665-620a13b9d341
# ╠═bbd6918a-1d4b-442d-81a7-e25c8f48c904
# ╠═09902aca-cb98-41f4-93f4-a7530c4373ba
# ╠═579e2373-4836-4fd5-83b9-c6715984634c
# ╠═f9c07cce-bada-11ef-36af-9f9bf85ccc38
# ╟─eb43c9e6-d47e-414f-987e-3dcc7a7c3124
# ╟─66243d50-9a2d-40db-b879-304f182eb46b
# ╠═4a5e1692-171c-43fa-8309-b73c8ecd38c5
# ╠═8a71f11b-3dab-4ee3-bc61-dc5327f37eb9
# ╠═19bf959c-b385-44d3-8311-8f2d22e67584
# ╟─cac10b59-aba7-4c0c-a699-02251a349dd4
# ╟─16b137b5-466a-4206-8d33-de2bf4a8850a
# ╟─1dca1d32-42c3-4636-aa86-66bcd3ea2dcd
# ╟─211988a4-e589-4b8f-98e0-c1db92fe5030
# ╠═7daa93a7-7968-4f7f-9e75-d27a412a7362
# ╟─017e4b49-dda4-4be7-a6e9-e5bd09218c36
# ╟─8145a2d1-c771-4530-8683-a057c4a4ad71
# ╟─cc0130e0-6fb1-482e-8b80-6726d4e7a08c
# ╟─d8a642cc-9f5c-4696-8314-20a998498eee
# ╟─35b9f77a-e939-42ef-9bfa-da2e30ffe659
# ╟─d2e7a649-2196-4f11-a790-96ce43396821
# ╠═c27d4852-4852-41ba-a1b6-759b8efde91e
# ╟─d4b00278-6933-45d8-bdfd-8d98682fd838
# ╠═a324dda4-c785-4ad0-8f2f-0b1e5ae5445b
# ╟─f01a8cd4-87de-44a1-ba58-d1dd78e4f5f3
# ╠═a75ce9c0-4e06-4950-a4f9-3baa5027a534
# ╠═4fe4f7bb-aee5-44d7-8503-916e3829ce65
# ╠═bd6a3907-bc0f-4cc0-b35c-e33af1ce9a5a
# ╠═047bced8-8cdd-47c9-a683-4bb5adb208fd
# ╠═0fb380cf-ba93-4e01-918a-c052b725fa0e
# ╟─8c4f63ba-5424-44b0-aa85-3e32e7793fb0
# ╟─32545342-1704-412f-a898-a3d1f5d3593d
# ╟─a490335d-0142-48a8-b3be-656bef552139
# ╟─85471dd1-df76-4718-a7ee-d501bb8a2030
# ╟─4bde3bb3-5dcd-42d4-ae48-b922e99a04bc
# ╠═fa3383e4-f5fa-4d47-b21b-684c6aacbd46
# ╠═1f36c90e-ec20-4bdb-8bd0-5f44f03a1748
# ╠═c4de5cb6-4b84-46b0-acc5-01aea794de3e
# ╠═3153b018-57f8-4d7a-bfb9-c4ca2886363a
# ╠═7551a678-2873-4d10-bef7-81fd0e275b9f
# ╠═621a2943-bb80-414e-babf-b5b2eac7269a
# ╠═cd8d78e3-18d3-4af1-ac69-71d007d41b7e
# ╟─0c9eb9d0-cda4-4f71-ba80-7f16680201d2
# ╟─1d5a95fb-2aeb-4a1e-808b-74b086970556
# ╟─3d2e0801-b075-4f91-a9db-55cdda05477d
# ╟─8cc42a7b-ba84-43d0-aafa-344712e882c6
# ╟─b78015f8-5edd-4a64-9e9b-b6856c525767
# ╠═f5bdb7b0-5a32-42cf-9461-75a8b5b513c2
# ╟─1d8a59eb-67b2-4242-9ab3-96fec2856622
# ╟─d6b22a0e-f2d1-4869-8405-46d163686061
# ╠═fb75ff46-58d4-4bb9-a7d9-36b66124ea90
# ╟─a6e960f6-fe3e-4f0e-9ee0-358010da6e39
# ╟─abd8fec6-b1d9-467b-a543-e479506b134f
# ╠═9db0f7d9-1d0c-40f7-89ef-a4786f8978ea
# ╠═a78f8bbb-52c5-45db-90c2-60837eefa6be
# ╟─5e6bc47b-bdf3-4b49-babb-6f3bf795fb7a
# ╠═4d5cd521-bffd-4f30-8a23-a04b3c5b9624
# ╟─2393409e-b5ec-4ce9-a818-c2a89bf50885
# ╟─5c500c1a-a15a-49bf-bd74-53985a484f19
# ╟─e56041b0-e722-40ae-8310-64d7069446a7
# ╠═ca46ccb7-523c-4540-8daa-f6ee5f134f8f
# ╟─c8f0cc9b-ff34-4baa-8234-b828df4fe894
# ╟─d0c9e1f1-ac1a-44c9-aaec-6c1cfd75f80d
# ╟─57f02857-1032-4ce2-a4fb-1cbeb6b6ac3d
# ╟─f2f496ab-4354-4d80-888e-89367efa00b0
# ╟─06b3ede3-93b4-4784-8298-272f8646b1d5
# ╟─d6bba283-75ed-40a4-a7ef-48f74c305ff2
# ╟─5d829be3-9bc9-4889-82d5-dc3a755b6435
# ╠═712f8974-ba9f-46b6-b569-fb9741567292
# ╠═8ca0791a-3da5-4b08-bf80-8221e4646297
# ╠═6652e22c-2600-4ab6-b84a-2205582de261
# ╠═f72e9d76-76d4-4da6-8c9c-c4fe433d4611
# ╟─90b9a18d-eaa8-4674-a092-82e6cb3ae8b1
# ╟─381df36b-1766-4be3-aa32-01d995ab72c0
# ╠═965ebc6f-b440-42da-8aa7-00b22e6d0a3f
# ╠═d7514616-736b-4cc3-a562-c23fb3d6924d
# ╠═4b953047-0f78-4c06-833f-3302a61882fc
# ╟─14dca280-a67f-4760-ac3a-dd7d88155754
# ╠═d451a113-249f-4a52-b753-45170c84a593
# ╟─9fe87e80-709e-41dc-b2d2-0569c34aeb5d
# ╟─a863f1ae-a06c-4829-9fb5-54c3a2f709aa
# ╠═eea8f1c3-a852-47d8-af77-588fa6499e84
# ╟─4150f328-cf18-4e07-ba35-49f9db3d54fb
# ╠═79e107ca-cc86-494f-9f87-f594d9bb874c
# ╠═ee1909ff-7334-4874-baf6-1f5a4d5e897d
# ╟─295db3bc-2007-4d18-928d-e7f67815c538
# ╠═76159175-f38f-4979-9132-dfc4b12d3b30
# ╠═7cb8b9ac-bb09-424d-a661-12860446a03c
# ╠═8961def3-e608-4aa4-842d-4a4565aa57f2
# ╠═9c083fda-ff0c-438e-abab-32e9845ee299
# ╟─b2dc3596-b7b5-4791-8251-83ca8070afd8
# ╠═2d4008e9-45fd-4bad-941d-51acb6c37a36
# ╠═558e7a9b-85fd-44e0-8e36-5c32e3476964
# ╠═3601dba1-43f9-4f50-b1fe-37a4f26985a8
# ╠═99a6cb5b-fe85-4476-abdc-471eadf6add6
# ╠═b55eaf74-9e69-4ab7-925c-db973cd0b667
# ╠═a85d7e97-9e0a-4e5e-af87-5816f2b0e44a
# ╠═43bab8a0-fcb4-496a-8a0b-34beb3dd1766
# ╠═66d87f17-36c5-476f-8696-e69ae55091b4
# ╠═c59be115-75b8-4432-8cd0-61c7c4449ac0
# ╟─a3c53584-7f2e-4324-ac45-ef3183a40c47
# ╠═17113056-e45f-4977-a022-22de535b1a49
# ╠═a0227fba-9aca-46b1-aa78-c56742890e74
# ╠═20d5a4d5-3d17-488a-bda6-62a180e459be
# ╠═250bd1af-97d6-4fd2-bd5f-2b10f3e3e212
# ╠═b0908503-1523-4a35-bc2f-5087650e537e
# ╠═8d95ddb9-bd58-4aea-9c8f-28e2f2119be2
# ╠═88c228fe-2aca-4058-803f-31d460ad1915
# ╠═fec66ee6-88f5-4491-8b80-d5e12d4fec01
# ╠═6e804a51-3a8c-45c0-b078-33af8ca96040
# ╠═85c68c99-35a8-49e3-bc60-6138e45637fc
# ╠═2a588744-d5c9-4e95-86e6-847fc9397663
# ╟─70369ddb-39b0-4ca0-9524-a12990e3478a
# ╠═d04ec167-203f-413b-aea5-a1a5cbdc70a2
# ╠═6ce08886-9271-4628-9afb-5f7f6778fef8
# ╠═3cd13eef-7954-4245-ac76-74eb4cf57b2c
# ╠═263d4146-e574-4576-8db4-77aa799d7c1a
# ╠═32a9b18f-ed5b-498d-bc49-129f2427616b
# ╠═39a435b3-7e9c-4dd5-9027-a751d4b20167
# ╠═753a58ec-4a82-480b-bcd3-50448ae09503
# ╠═9e617891-4a1c-4e56-a7b2-177168d572e0
# ╟─42e9c861-a25a-4efc-bc6f-8d63c347c082
# ╟─5d6f4d9a-fbe1-4b03-85de-d02d6f505c63
# ╠═e387491c-a191-47cf-b5c4-8ae6ba4f561a
# ╠═990cd6d9-ef7e-4a62-83b8-7724515dfdcc
# ╠═9c8f06fa-946a-4d2c-af3d-593a58bd0f92
# ╠═4e5f5f05-7f1a-41cb-a7e9-431f53e97bb5
# ╠═41ea7cc3-fe7e-4114-8725-1180b8eea871
# ╟─926154ef-34d7-452a-8627-bfb1c050a1b1
# ╠═ba0a576e-4b9b-49e1-bfbc-33955dcec808
# ╠═edd1c981-56f6-44d8-a805-0b85c963b39b
# ╠═bc173213-4ece-4358-9d70-f1f9f7edea8e
# ╠═a94febb7-163e-4f6e-91e5-1a1871e8ea39
# ╟─ec661503-a2e6-4a3d-a244-a4cafebcff60
# ╠═f8b9f59c-71ef-426c-bf25-880d6660e9aa
# ╠═2e8f3065-343f-482c-b3d3-0ace5337343d
# ╠═c73544a6-3b24-4f06-8b00-0e0e5059e5ee
# ╠═beb474e9-0694-457e-861a-e8ec9ddf865b
# ╠═b032789f-9f7a-45a3-a3fd-9186598bab28
# ╟─e0690c90-01f7-4e8d-9faf-6f9149a4031d
# ╠═8a6c6123-36a2-4769-8f22-876eeccf3963
# ╠═65aef305-fdaa-4262-8bc5-a22f5a0b42ca
# ╠═4dcbb097-8eda-4d34-8d9f-88ff185ef108
# ╠═a6e1fc57-2750-4765-925e-d321b23866e1
# ╠═50213066-8eb4-4d89-a660-da71911104cc
# ╠═a75c3dc3-03bf-4ed4-8fe8-59a63a022030
# ╠═6ec4751a-dd24-4da9-a864-1bde6aa5632d
# ╠═c63dcc08-911f-49c5-9316-af07bc8c36b4
# ╟─b88b0d16-a178-4395-a54f-836d4090ccec
# ╠═763b603a-f81b-4d9e-a6eb-7dc345e3e2df
# ╠═024b9f4e-600d-4d68-962d-a66c9b5a2f23
# ╠═343945bd-203d-4a9b-b46d-2a12b46ae32c
# ╠═44a7a892-336f-4e52-a66d-230afbb75a7b
# ╠═4b53805c-c9ed-4708-b8a1-4922ce6d5a5d
# ╠═4e579e09-2476-432d-8eb4-7a1080376998
# ╠═02f2f82c-4c69-4fc3-90b1-ea0a8796fe12
# ╠═1943119c-3007-4495-8148-1c346f052a72
# ╠═d178630c-289d-4917-a182-8f95409d14c3
# ╠═b42769d2-16d2-449d-8925-fd5b71ddced3
# ╠═06b60440-696e-4067-9772-2b91f4be9cbd
# ╠═54e63222-6434-4f63-81fa-ab6d5d8fd7da
# ╟─7ea9a771-d692-4ce1-88d8-0a31047dca3d
# ╠═f2f65b41-edb4-411f-aa4c-cef53e5aaefe
# ╠═6ee18e30-996b-4e30-b0ec-39cf859d5d5a
# ╠═8740b8b8-06a9-420b-8526-2adf5f8219a3
# ╠═e49c777a-a2c3-4725-a9ba-abb7d8536d2e
# ╠═8fe403f4-2f3b-482b-91b6-e49672f98d23
# ╠═88fb0c51-9005-44de-bb15-660ddb9d3165
# ╠═bbd11d02-29ca-4b4d-b08f-2c4b460fd5c5
# ╠═e4c1292a-be82-466f-a586-430f7d34f9bb
# ╟─3a1ce3ee-864b-4144-a9da-3d66bd1870e9
# ╟─da1dc114-8529-4ef8-b3af-2acb066fa02f
# ╟─e84a0ef5-a6e9-4c66-a44c-81e7ef70351d
# ╟─9b8392a5-fc9f-433f-9b23-dab2dd472021
# ╟─8d1c2906-f7bf-4557-81af-43b8d9b2a260
# ╠═520e22a5-e08f-4718-bf2f-f1e7e627a417
# ╠═8fd5f30a-ef31-442f-9124-8ca86a3b9b56
# ╟─5e8a4a93-7316-4b61-b1e0-bce46597f614
# ╠═dde10d51-3d75-4152-a7da-c05b1affaa6a
# ╟─8becb816-b96c-4df5-8af0-2464c7188d34
# ╠═9bb200cc-09db-4f19-bf29-5b4a1cc78574
# ╟─28932b7c-35b0-4e8e-91f2-4dfd5a41407e
# ╟─9dd416bc-eead-4a95-820c-6f0cdf532243
# ╟─077a5018-538c-4c33-b703-f7258849f2f8
# ╠═d00abaf7-10de-4f7a-9bb9-a476bf4b8351
# ╟─a4156697-08c0-418f-95cd-bffe4295eb82
# ╠═07d655d8-91ba-49d4-8e1d-d738826b8099
# ╟─b3fa82ec-6686-40bf-8709-4d281dcb518f
# ╠═95aa37fe-3cb1-4a9a-952a-2bb2ea815a1b
