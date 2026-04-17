### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 3ef93246-657d-4e77-9bf0-8380c64bfcfd
using Pkg; Pkg.activate("..")

# ╔═╡ 55cdebd2-0881-11ef-2722-91de1447877a
using Markdown, InteractiveUtils

# ╔═╡ a355b0ba-baaf-49f4-a5dc-965364a884f0
using Catalyst, OrdinaryDiffEq

# ╔═╡ 00fd6d49-f561-42e9-9413-d33af92f83dc
using StatsPlots, PlutoUI; TableOfContents()

# ╔═╡ 7ae714c4-d25d-4f9f-ab3d-cc067db9c156
using ForwardDiff

# ╔═╡ dfce0717-e33e-4b6b-bb5c-d8754a74aba4
using Turing

# ╔═╡ 31d294d1-3a1f-41db-abff-54f2a67c7ed9
md"""
# Exercise: Fermenter - Monod kinetics - Sensitivity analysis
"""

# ╔═╡ 5ffe7dcb-620d-4f22-95fe-2f77cda6fbe7
md"""
In one of the previous practica, we were introduced to a fermenter in which biomass $X$ [$g/L$] grows by breaking down substrate $S$ [$g/L$]. The reactor is fed with an inlet flow rate $Q_{in}$ [$L/h$], which consists of a (manipulable) input concentration of substrate $S_{in}$ [$g/L$]. This process was modelled using Monod kinetics, resulting in the model below:

$$\begin{eqnarray*}
S + X \xrightarrow[\quad\quad]{k} (1 + Y) \, X \quad\quad\quad\quad \textrm{with} \quad k = \cfrac{\mu_{max}}{S + K_s} \, .
\end{eqnarray*}$$
"""

# ╔═╡ d635d577-4d6d-40d7-a84c-3871981d59a4
md"""
Suppose that at $t=0$ no substrate $S$ is present in the reactor but that there is initially some biomass with a concentration of $0.0005\;g/L$. The (default) parameter values are $Q = 2.0$, $V = 40.0$ and $Y = 0.67$. The parameters $\mu_{max}$, $K_s$ and $S_{in}$ will be assigned later.
"""

# ╔═╡ 6ec6da23-853b-4129-94cf-67b5cadb1f95
md"""
The *reaction network object* model for this problem could be defined as:
"""

# ╔═╡ 935ca610-7a7a-4692-8908-fc26abb880b4
fermenter_monod = @reaction_network begin
	@species S(t)=0.0 X(t)=0.0005
	@parameters Q=2.0 V=40.0 Y=0.67 Sin μmax Ks
	μmax/(S + Ks), S + X --> (1 + Y)*X
	# Alternative:
	# mm(S, μmax, Ks)*X, S => Y*X
    Q/V, (S, X) --> ∅
    Q/V*Sin, ∅ --> S
end

# ╔═╡ 79e6056a-881c-442f-8989-5bc284d3d777
md"""
which resulted in the following differential equations:
"""

# ╔═╡ fa93e2c3-8b43-418e-ba24-406645b2e397
md"""
$$\begin{eqnarray*}
\cfrac{dS}{dt} &=& \cfrac{Q}{V} \left(S_{in} - S \right) - \mu_{max}\cfrac{S}{S + K_s} X\\
\cfrac{dX}{dt} &=& -\cfrac{Q}{V} X + Y \mu_{max}\cfrac{S}{S + K_s} X
\end{eqnarray*}$$
"""

# ╔═╡ 06730f54-7293-43f2-b772-84eec3e5528a
md"""
Convert the system to a symbolic differential equation model and verify, by analyzing the differential equation, that your model is correctly implemented.

Keep in mind that `mm(S, μmax, Ks)` stands for $\mu_{max} \, \cfrac{S}{S + K_s}$.
"""

# ╔═╡ 7f8b7a2e-bc65-4b51-ad59-bd7ac98604dd
# osys = missing

# ╔═╡ 1ee7359c-2757-44cd-8650-09b762cf30d0
md"""
## Goals of this exercise
"""

# ╔═╡ 55f1d688-0c53-481b-9965-5e92ca87ad83
md"""
Compute the following in a timespan of `(0.0, 100.0)`$\,h$:

- The sensitivities of $S$ and $X$ w.r.t. $\mu_{max}$, $K_s$ and $S_{in}$.

Plot the following:
- A figure with the sensitivity functions of $S$ and $X$ w.r.t. $S_{in}$.
- A figure with the sensitivity functions of $S$ w.r.t. $\mu_{max}$, $K_s$ and $S_{in}$.
- A figure with the sensitivity functions of $X$ w.r.t. $\mu_{max}$, $K_s$ and $S_{in}$.

Interpret your results.

Use the following parameter values $\mu_{max}$ = `0.40`, $K_s$ = `0.015` and $S_{in}$ = `0.022` for the calculations of the sensitivities.
"""

# ╔═╡ 6c4e3c09-4b84-4f5c-8739-2ac18e6f2af6
md"""
Initialize a vector `u0` with the initial conditions, set the timespan and initialize a vector `parms` with the parameter values:
"""

# ╔═╡ 2ee277e5-ce4a-4ade-be0e-9bba7a4dc08c
# u0 = missing   # in principle not necessary since we will use the default

# ╔═╡ 3fdc6b17-cdeb-4dc5-8886-9d3a62caac8d
# tspan = missing

# ╔═╡ 201dfb54-2056-4846-a64c-ff4951dc084d
md"""
For practical reasons, we will define the time step size as `dt = 0.5`. 
"""

# ╔═╡ 5627068f-6442-43bb-bacb-bcba917372f3
# dt = missing

# ╔═╡ 0139da85-02e3-4021-9b39-84af7e68d428
md"""
For practical reasons, we will use the variables `μmax_val`, `Ks_val`, and `Sin_val` to store the parameter values that are used for the calculation of the sensitivity functions.
"""

# ╔═╡ 0f995929-4d2b-4a7a-8da1-04e4d501385f
begin
	# μmax_val = missing
	# Ks_val = missing
	# Sin_val = missing
end;

# ╔═╡ 79b0eb65-5a0f-40b3-aa97-4088421c562e
# parms = missing   # no need to include Q, V, and Y since we will use the default

# ╔═╡ af882cf4-51fd-45ff-9b05-cfd52d6467b1
md"""
## Preliminary simulation
"""

# ╔═╡ f0f4fa14-6f99-4f21-a743-be61e08444a7
md"""
Create the ODE problem and store it in `oprob`. Next, solve the ODE problem using `Tsit5()` and `saveat=dt`, and store the solution in `osol`. Finally plot the results.
"""

# ╔═╡ 8b2f23f6-80b2-4e63-942e-e5cd17d8ba72
# oprob = missing

# ╔═╡ 89a31c32-88a4-479f-a688-ffcb75ee8e91
# osol = missing

# ╔═╡ 51a9b7e6-8ad9-477d-9596-ffd614df2c79
# missing

# ╔═╡ 570ebda9-b187-42ab-a761-bed0fa3ce097
md"""
## Local Sensitivity Analysis (LSA)
"""

# ╔═╡ 0343674d-32d1-49e1-a7c1-47196ae760ed
md"""
### Setting up function
"""

# ╔═╡ 693844d0-3858-4861-bae0-b47e78809f17
md"""
Write a solution function with as argument a vector of the parameters (that you want the sensitivity on), and that returns the outputs.
"""

# ╔═╡ 9622f7ca-f71a-4ad9-a309-d7d10a1c3e3b
# function fermenter_monod_sim(parms)
# 	missing
# 	...
# end

# ╔═╡ 4a5971b1-f4d0-43b6-805f-e17f5052ae92
md"""
Make two functions based on the solution function that each returns a single output, hence, one function that returns the output $S$, and another function that returns the output $X$.
"""

# ╔═╡ f40c6402-3c28-4a7d-b629-83507a9f29bd
# fermenter_monod_sim_S(parms) = missing

# ╔═╡ 3ae5bd00-2e06-4789-aab3-d897824d5e29
# fermenter_monod_sim_X(params) = missing

# ╔═╡ 6913fb1b-1986-4b37-819e-a6a7bd91f1a5
md"""
### Compute the outputs
"""

# ╔═╡ 4bd2bcca-9c42-4333-b062-2aaa9f7be3fe
md"""
Make the time vector.
"""

# ╔═╡ fbd98975-aa32-46ae-8db0-0e65cdf48309
# t_vals = missing

# ╔═╡ 93791eb3-1eaa-4146-90b5-c4811fb3485b
md"""
Compute the two outputs $S$ and $X$ for the given parameter values.
"""

# ╔═╡ dc0557d6-81b9-4759-8ed7-3129f60c6dc3
# S_sim = missing

# ╔═╡ 95bc683c-f6e6-4b42-b90b-b5a863edd4d5
# X_sim = missing

# ╔═╡ 35717c5a-9895-4e62-a890-03ef04aa07b9
md"""
### Compute the sensitivities
"""

# ╔═╡ fa970c0e-fb3b-486f-bbc1-345d44f8f0da
md"""
Using `ForwardDiff.jacobian` to compute the sensitivities for the single ouputs $S$ and $X$. Hence, you need to call `ForwardDiff.jacobian` twice.
"""

# ╔═╡ 07f2b8f8-5162-4077-a8a5-68fdec646644
md"""
#### Absolute sensitivities
"""

# ╔═╡ 64354302-f4cc-4592-9302-5db0f5bccb2e
# sens_S = missing

# ╔═╡ 49a94b9a-a543-495e-b4f1-c8579e59304d
# sens_X = missing

# ╔═╡ 9cace6c1-e678-4dd7-8705-92a55eb32fa9
md"""
Extract the (absolute) sensitivities of the outputs on the different parameters.
"""

# ╔═╡ a6dc2b60-6a0a-4140-892e-02cde8dc79d3
# begin
# 	sens_S_on_μmax = missing
# 	sens_S_on_Ks   = missing
# 	sens_S_on_Sin  = missing
# end

# ╔═╡ f806c243-9032-46b7-add3-4714344691c7
# begin
# 	sens_X_on_μmax = missing
# 	sens_X_on_Ks   = missing
# 	sens_X_on_Sin  = missing
# end

# ╔═╡ f4d0aaaa-bb3f-4335-901a-dbbcde643c78
md"""
#### Normalized sensitivities
"""

# ╔═╡ 5bf3a62d-d2aa-4653-8ee8-e90caa9504e8
md"""
Compute the normalized sensitivities.
"""

# ╔═╡ 76846731-929c-408f-a3de-970581c497e9
# begin
# 	sens_S_on_μmax_rel = missing
# 	sens_S_on_Ks_rel   = missing
# 	sens_S_on_Sin_rel  = missing
# end

# ╔═╡ b6c57444-547c-4e82-8526-6a30566e07c5
# begin
# 	sens_X_on_μmax_rel = missing
# 	sens_X_on_Ks_rel   = missing
# 	sens_X_on_Sin_rel  = missing
# end

# ╔═╡ 7f0b9e63-2856-4601-9f14-c20c6b1b3707
md"""
### Plotting + questions
"""

# ╔═╡ 5388c2a7-5a11-4da8-be09-46045cde8a4e
md"""
Plot the sensitivity functions of $S$ and $X$ on $S_{in}$. Provide a suitable title (`title="..."`), labels (`label=["..." "..."]`) and an x-label (`xlabel="..."`), and set the line width to 2 (`linewidth=...`).
"""

# ╔═╡ db840c76-a6c6-49fb-a0bb-d9149f947bc0
# missing

# ╔═╡ d41375ef-6958-4705-a417-4c6a491232ee
md"""
!!! questions
	Interpret your results. Try to answer the following question(s):
	- Which output variable $S$ or $X$ is most sensitive to $S_{in}$ in steady state?
	- Why is the sensitivity function of $S$ on $S_{in}$ at first positive but then becomes zero?
"""

# ╔═╡ 39c36eb4-5f15-4577-9b84-116e29d7891d
md"""
Answers:
- missing
- missing
"""

# ╔═╡ be89600a-4927-4afc-9813-d8a70adb2852
md"""
Plot the sensitivity functions of $S$ on $\mu_{max}$, $K_s$ and $S_{in}$. Provide a suitable title (`title="..."`), labels (`label=["..." "..." "..."]`) and an x-label (`xlabel="..."`), and set the line width to 2 (`linewidth=...`).
"""

# ╔═╡ c0223da4-9959-48d0-b607-633b2e82986c
# missing

# ╔═╡ ff86a29f-9308-473b-aa1c-dfd4af8179c7
md"""
!!! questions
	Interpret your results. Try to answer the following question(s):
	- Which parameter $\mu_{max}$, $K_s$ or $S_{in}$ affects the output $S$ the most in steady state?
	- Why is the sensitivity function of $S$ w.r.t. $K_s$ positive?
	- Why is the sensitivity function of $S$ w.r.t. $\mu_{max}$ negative?
"""


# ╔═╡ 430494da-402e-409d-85fc-c87cdfa3b427
md"""
Answers:
- missing
- missing
"""
#=
- It seems like μmax is affecting S the most in steady state, and its influence is negative; hence, the larger μmax, the smaller S.
- The sensitivity function of S w.r.t. Ks is positive, because the larger Ks, the smaller the reaction rate r (S => Y*X). Hence, less X will be produced; so less S will be consumed. Remember that r = μmax*S*X/(S + Ks) and Ks is in the denominator.
- The sensitivity function of S w.r.t. μmax is negative, because the larger μmax, the larger the reaction rate r (S => Y*X). Hence, more X will be produced; so more S will be consumed. Remember that r = μmax*S*X/(S + Ks) and μmax is in the numerator.
=#

# ╔═╡ 16a84fdb-8ce2-45b9-bfb7-7f4e1284a1d7
md"""
Plot the sensitivity functions of $X$ on $\mu_{max}$, $K_s$ and $S_{in}$. Provide a suitable title (`title="..."`), labels (`label=["..." "..." "..."]`) and an x-label (`xlabel="..."`), and set the line width to 2 (`linewidth=...`).
"""

# ╔═╡ 53134149-0bf7-41c1-9b35-e5037744211f
# missing

# ╔═╡ 355ca6a7-466b-4969-ab48-28e2257f9810
md"""
!!! questions
	Interpret your results. Try to answer the following question(s):
	- Which parameter $\mu_{max}$, $K_s$ or $S_{in}$ affects the output $X$ the most in steady state?
	- Why is the sensitivity function of $X$ w.r.t. $K_s$ negative?
	- Why is the sensitivity function of $X$ w.r.t. $\mu_{max}$ positive?
"""

# ╔═╡ 9a18a8bd-27d3-4771-a844-f6b65c7c8918
md"""
Answers:
- missing
- missing
- missing
"""
#=
- It seems like Sin is affecting X the most in steady state, and its influence is positive; hence, the larger Sin, the larger X.
- The sensitivity function of X w.r.t. Ks is negative, because the larger Ks, the smaller the reaction rate r (S => Y*X). Hence, less X will be produced. Remember that r = μmax*S*X/(S + Ks) and Ks is in the denominator.
- The sensitivity function of X w.r.t. μmax is positive, because the larger μmax, the larger the reaction rate r (S => Y*X). Hence, more X will be produced. Remember that r = μmax*S*X/(S + Ks) and μmax is in the numerator.
=#

# ╔═╡ 6536f76a-0f6d-4155-ab3b-38e2e3189886
md"""
## Monte Carlo error propagation
"""

# ╔═╡ 9fd59148-6993-4598-af95-ac85b6a4bb45
md"""
Let's take a look at how the uncertainty propagates in your model by using Monte Carlo simulations. Based on some literature search you can assume that the parameter `Ks` is normally distributed around the original value with a standard deviation of 20%.
"""

# ╔═╡ 5a559bc2-3d55-499c-93cf-4b0845ffc4b0
md"""
Start with making a Turing model in which you implement the prior and return the solution of the solved problem.
"""

# ╔═╡ ecc28b8e-0eaf-489f-8104-7c6383d4e01f
md"""
!!! note
	Instead of creating a new ODEProblem everytime, you can simply remake an old ODEProblem by using `new_prob = remake(oprob, p = [:Ks => Ks])`. 
"""

# ╔═╡ 26ae811d-11fd-44d8-b328-79381fa7472c
# @model function monod_deviation()
# 	Ks_dev ~ missing # 20% standard deviation
#   new_prob = missing
# 	return missing
# end

# ╔═╡ 62791c48-5e38-46fe-8336-3a21003f3fe3
md"""
Now get sample 100 solutions from your monod_deviation model.
"""

# ╔═╡ 9e049653-347a-4198-90fd-23f15c404d0c
# missing

# ╔═╡ f11827e5-de69-4a0c-b385-2c6ea031bb79
md"""
We can now visualise the results of our Monte Carlo simulation by looping through our solutions and plotting them. It is advised to put `label = false` (avoids 100 labels), `color=:gray` (neutral color) and `alpha=0.4` (makes it more transparent) and `linestyle=:dot`.
"""

# ╔═╡ 9d0f8147-3c82-4ace-a408-20e32b4d8ea0
# begin
#     p1 = plot(title="Monte Carlo of logistic growth") # make an empty plot
#     for missing in missing
#         missing
#     end
# 	  plot!(osol, label = "Original simulation", linewidth = 3) # Plots original result
#     p1
# end

# ╔═╡ 30c270f1-cbd5-401b-92db-7123efb01db6
md"""
As can be seen from the graph, the uncertainty is rather high, especially if we want to carefully monitor the concentration of biomass through time. Note that, since we have a sample of values for $S$ and $X$ at every time step, it is also perfectly possible here to quantify the uncertainty by calculating the standard deviation at every time step, but the graph already gives a clear indication.
"""

# ╔═╡ 16098df3-35c0-4d06-a482-35137fd53b97
md"""
Now create a vector x1 that stores the values of the biomass at the end of the simulations and make a histogram of it. Put `xlims = [0.010, 0.015]`.
"""

# ╔═╡ de5736d4-1c2f-4b0c-a84b-6d8b07d29a0a
# missing

# ╔═╡ 1b1d0997-b854-4445-8c19-e73db6784a9b
# missing

# ╔═╡ 43aaf563-f799-43e4-9f4c-8a19be1586dd
md"""
### Determininig optimal measurements
"""

# ╔═╡ f7f1f035-b48e-4108-bb84-0c00ab68a93a
md"""
As a process operator, you may want to reduce the uncertainty in biomass simulations. However, since taking measurements can be costly, it is important to perform them at the most informative time points.

Local sensitivity analysis is a useful tool in this context, as it indicates when the model output is most sensitive to changes in parameter values. At these time points, accurate measurements provide the greatest insight into the true parameter values.

Start off by identifying the time at which the biomass concentration is most sensitive to variations in `Ks`.
"""

# ╔═╡ 8dbadbc6-acd9-42f5-abf4-81fe6aef7fb7
md"""
!!! note
	You can find the index of the maximum value in a vector with the `argmax()` function. For example `argmax([1, 2, 5, 4])` will give you 3.
"""

# ╔═╡ ddfb530a-ba79-43c0-973a-4deb8cad1d33
# begin
# 	t_star_idx = missing
# 	t_star = t_vals[t_star_idx]
	
# 	plot(t_vals, abs.(sens_X_on_Ks_rel),
#	    title="Sensitivity of Ks wrt. X", xlabel="Time (h)",
#	    ylabel="normalized sensitivity")
# 	vline!([t_star], label="t* = $(round(t_star, digits=2)) h", linestyle=:dash, color=:red)
# end

# ╔═╡ 8359ea0e-a916-48f4-a07a-c905c3efaf15
begin 
	μmax_real = 0.42
	Ks_real = 0.0165
	Sin_real = 0.0207
	real_solution = solve(remake(oprob, p = [:μmax => μmax_real, :Ks => Ks_real,
											:Sin => Sin_real]));
end;

# ╔═╡ 9654dea2-e13e-4aa9-8a5e-c3c368c26e24
md"""
We have defined a function `real_solution(t)` in this notebook that returns a hypothetical measurement (`S_meas, X_meas`) at the given timepoint `t`. Use this function to get a measurement at the optimal time (`t_star`).
"""

# ╔═╡ 77eab4ec-5012-40dd-b702-9fb48982ff38
# S_meas, X_meas = real_solution(t_star)

# ╔═╡ 3e93a001-9401-4a1f-b283-c79fd0c43434
md"""
#### Conditioning on the measurements
"""

# ╔═╡ 684f9163-ac84-49be-beb1-66d7bdc21da8
md"""
Now we can use this measurement to calibrate our model. In this part we will use Monte Carlo Markov chain to show the decrease in uncertainty. 

We start by creating a Turing model, for which you can use the priors from the first part of this exercise. For the standard deviation `σ_X` you can assume a value of 0.0002 and a value of 0.00055 for `σ_S`.

Make sure to truncate the domain of your parameters to a reasonable domain or you might run into issues with the solver.
"""

# ╔═╡ 4fe68607-d67e-4e71-865a-b10a06e5908d
# @model function monod_meas()
#     σ_X = missing
#     σ_S = missing

#     Ks_dev ~ missing

#     sol = missing
    
#     S_pred, X_pred = sol(t_meas)

#     X_meas ~ missing
#     S_meas ~ missing
# end

# ╔═╡ 11e105de-270c-434a-bfb4-9118dfc9e461
# monod_cond_model = missing

# ╔═╡ e46afb98-16de-4d19-86b4-b920b6b929c5
md"""
### MCMC
"""

# ╔═╡ b8cd84de-57e4-4d8b-92da-e7c5b6d627bd
# monod_chain = missing

# ╔═╡ e76a7bf7-4d8d-4f2f-8272-b63e49cd4ae9
md"""
Plot the chain to check if it properly converged to the posterior distribution.
"""

# ╔═╡ 834c3500-7556-4999-b465-08a8a07c7f28
# missing

# ╔═╡ eeaf3255-95c8-4813-bdd7-38635e546606
md"""
Now, get the sampled posterior solutions and afterwards plot the resulting simulations.
"""

# ╔═╡ a66b6d91-e529-4c20-b5d6-91830eae1437
# missing

# ╔═╡ c4192ff2-fca0-4ab1-8ac1-ce9771bc69d3
md"""
Plot the resulting Monte Carlo simulation with the original simulations.
"""

# ╔═╡ 54f20cff-21fc-4169-9f47-5e670b7ec621
# begin
#     p = plot(ylabel = "concentration (g/l)", xlabel = "Time (hours)")
#     for missing in missing
#         missing
#     end
# 	plot!(osol, label = ["Original simulation S" "Original simulation X"], linewidth = 3)

# 	p
# end

# ╔═╡ 988aab6d-7193-4f3d-bd08-c3c36c920f81
md"""
Lastly, create a histogram of the value of X at the end of the simulations based on your posterior distribution of Ks. Put the `xlims = [0.010, 0.015]` to compare with the prior distribution. 
"""

# ╔═╡ 2f55a771-0d37-48e5-a9d9-24693bf1d177
# missing

# ╔═╡ 18059d0c-cf37-44a1-86a8-859ac66a66f0
# missing

# ╔═╡ Cell order:
# ╟─31d294d1-3a1f-41db-abff-54f2a67c7ed9
# ╠═55cdebd2-0881-11ef-2722-91de1447877a
# ╠═3ef93246-657d-4e77-9bf0-8380c64bfcfd
# ╠═a355b0ba-baaf-49f4-a5dc-965364a884f0
# ╠═00fd6d49-f561-42e9-9413-d33af92f83dc
# ╠═7ae714c4-d25d-4f9f-ab3d-cc067db9c156
# ╠═dfce0717-e33e-4b6b-bb5c-d8754a74aba4
# ╟─5ffe7dcb-620d-4f22-95fe-2f77cda6fbe7
# ╟─d635d577-4d6d-40d7-a84c-3871981d59a4
# ╟─6ec6da23-853b-4129-94cf-67b5cadb1f95
# ╠═935ca610-7a7a-4692-8908-fc26abb880b4
# ╟─79e6056a-881c-442f-8989-5bc284d3d777
# ╟─fa93e2c3-8b43-418e-ba24-406645b2e397
# ╟─06730f54-7293-43f2-b772-84eec3e5528a
# ╠═7f8b7a2e-bc65-4b51-ad59-bd7ac98604dd
# ╟─1ee7359c-2757-44cd-8650-09b762cf30d0
# ╟─55f1d688-0c53-481b-9965-5e92ca87ad83
# ╟─6c4e3c09-4b84-4f5c-8739-2ac18e6f2af6
# ╠═2ee277e5-ce4a-4ade-be0e-9bba7a4dc08c
# ╠═3fdc6b17-cdeb-4dc5-8886-9d3a62caac8d
# ╟─201dfb54-2056-4846-a64c-ff4951dc084d
# ╠═5627068f-6442-43bb-bacb-bcba917372f3
# ╟─0139da85-02e3-4021-9b39-84af7e68d428
# ╠═0f995929-4d2b-4a7a-8da1-04e4d501385f
# ╠═79b0eb65-5a0f-40b3-aa97-4088421c562e
# ╟─af882cf4-51fd-45ff-9b05-cfd52d6467b1
# ╟─f0f4fa14-6f99-4f21-a743-be61e08444a7
# ╠═8b2f23f6-80b2-4e63-942e-e5cd17d8ba72
# ╠═89a31c32-88a4-479f-a688-ffcb75ee8e91
# ╠═51a9b7e6-8ad9-477d-9596-ffd614df2c79
# ╟─570ebda9-b187-42ab-a761-bed0fa3ce097
# ╟─0343674d-32d1-49e1-a7c1-47196ae760ed
# ╟─693844d0-3858-4861-bae0-b47e78809f17
# ╠═9622f7ca-f71a-4ad9-a309-d7d10a1c3e3b
# ╟─4a5971b1-f4d0-43b6-805f-e17f5052ae92
# ╠═f40c6402-3c28-4a7d-b629-83507a9f29bd
# ╠═3ae5bd00-2e06-4789-aab3-d897824d5e29
# ╟─6913fb1b-1986-4b37-819e-a6a7bd91f1a5
# ╟─4bd2bcca-9c42-4333-b062-2aaa9f7be3fe
# ╠═fbd98975-aa32-46ae-8db0-0e65cdf48309
# ╟─93791eb3-1eaa-4146-90b5-c4811fb3485b
# ╠═dc0557d6-81b9-4759-8ed7-3129f60c6dc3
# ╠═95bc683c-f6e6-4b42-b90b-b5a863edd4d5
# ╟─35717c5a-9895-4e62-a890-03ef04aa07b9
# ╟─fa970c0e-fb3b-486f-bbc1-345d44f8f0da
# ╟─07f2b8f8-5162-4077-a8a5-68fdec646644
# ╠═64354302-f4cc-4592-9302-5db0f5bccb2e
# ╠═49a94b9a-a543-495e-b4f1-c8579e59304d
# ╟─9cace6c1-e678-4dd7-8705-92a55eb32fa9
# ╠═a6dc2b60-6a0a-4140-892e-02cde8dc79d3
# ╠═f806c243-9032-46b7-add3-4714344691c7
# ╟─f4d0aaaa-bb3f-4335-901a-dbbcde643c78
# ╟─5bf3a62d-d2aa-4653-8ee8-e90caa9504e8
# ╠═76846731-929c-408f-a3de-970581c497e9
# ╠═b6c57444-547c-4e82-8526-6a30566e07c5
# ╟─7f0b9e63-2856-4601-9f14-c20c6b1b3707
# ╟─5388c2a7-5a11-4da8-be09-46045cde8a4e
# ╠═db840c76-a6c6-49fb-a0bb-d9149f947bc0
# ╟─d41375ef-6958-4705-a417-4c6a491232ee
# ╠═39c36eb4-5f15-4577-9b84-116e29d7891d
# ╟─be89600a-4927-4afc-9813-d8a70adb2852
# ╠═c0223da4-9959-48d0-b607-633b2e82986c
# ╟─ff86a29f-9308-473b-aa1c-dfd4af8179c7
# ╟─430494da-402e-409d-85fc-c87cdfa3b427
# ╟─16a84fdb-8ce2-45b9-bfb7-7f4e1284a1d7
# ╠═53134149-0bf7-41c1-9b35-e5037744211f
# ╟─355ca6a7-466b-4969-ab48-28e2257f9810
# ╟─9a18a8bd-27d3-4771-a844-f6b65c7c8918
# ╟─6536f76a-0f6d-4155-ab3b-38e2e3189886
# ╟─9fd59148-6993-4598-af95-ac85b6a4bb45
# ╟─5a559bc2-3d55-499c-93cf-4b0845ffc4b0
# ╟─ecc28b8e-0eaf-489f-8104-7c6383d4e01f
# ╠═26ae811d-11fd-44d8-b328-79381fa7472c
# ╟─62791c48-5e38-46fe-8336-3a21003f3fe3
# ╠═9e049653-347a-4198-90fd-23f15c404d0c
# ╟─f11827e5-de69-4a0c-b385-2c6ea031bb79
# ╠═9d0f8147-3c82-4ace-a408-20e32b4d8ea0
# ╟─30c270f1-cbd5-401b-92db-7123efb01db6
# ╟─16098df3-35c0-4d06-a482-35137fd53b97
# ╠═de5736d4-1c2f-4b0c-a84b-6d8b07d29a0a
# ╠═1b1d0997-b854-4445-8c19-e73db6784a9b
# ╟─43aaf563-f799-43e4-9f4c-8a19be1586dd
# ╟─f7f1f035-b48e-4108-bb84-0c00ab68a93a
# ╟─8dbadbc6-acd9-42f5-abf4-81fe6aef7fb7
# ╠═ddfb530a-ba79-43c0-973a-4deb8cad1d33
# ╟─8359ea0e-a916-48f4-a07a-c905c3efaf15
# ╟─9654dea2-e13e-4aa9-8a5e-c3c368c26e24
# ╠═77eab4ec-5012-40dd-b702-9fb48982ff38
# ╟─3e93a001-9401-4a1f-b283-c79fd0c43434
# ╟─684f9163-ac84-49be-beb1-66d7bdc21da8
# ╠═4fe68607-d67e-4e71-865a-b10a06e5908d
# ╠═11e105de-270c-434a-bfb4-9118dfc9e461
# ╟─e46afb98-16de-4d19-86b4-b920b6b929c5
# ╠═b8cd84de-57e4-4d8b-92da-e7c5b6d627bd
# ╟─e76a7bf7-4d8d-4f2f-8272-b63e49cd4ae9
# ╠═834c3500-7556-4999-b465-08a8a07c7f28
# ╟─eeaf3255-95c8-4813-bdd7-38635e546606
# ╠═a66b6d91-e529-4c20-b5d6-91830eae1437
# ╟─c4192ff2-fca0-4ab1-8ac1-ce9771bc69d3
# ╠═54f20cff-21fc-4169-9f47-5e670b7ec621
# ╟─988aab6d-7193-4f3d-bd08-c3c36c920f81
# ╠═2f55a771-0d37-48e5-a9d9-24693bf1d177
# ╠═18059d0c-cf37-44a1-86a8-859ac66a66f0
