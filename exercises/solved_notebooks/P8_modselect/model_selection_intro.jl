### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ e1e7bc8e-7264-4cbc-98d2-aa73679fa2df
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ 5f4fea06-0632-11ef-102e-21f5606d2056
using Markdown

# ╔═╡ 489b5399-fe4c-481a-834f-0101bbe28cea
using InteractiveUtils

# ╔═╡ ab4a1911-461e-4878-9258-931fc2f1ea06
using PlutoUI; TableOfContents()

# ╔═╡ cdcebbb1-40e0-457f-a6ec-b769f6b1f2e9
using DifferentialEquations, Plots

# ╔═╡ 22489bd4-ab64-4bf9-ad03-5372ea273935
using Catalyst

# ╔═╡ 4ad72a6a-f541-4d0a-8392-56369e29ac96
using Turing, StatsPlots, StatsBase

# ╔═╡ 70590d28-058a-44df-8cf6-092d8c87438c
using LinearAlgebra, Optim

# ╔═╡ 427b509f-d08e-4d93-99ba-a79f9c244b28
md"""
# Introduction to model selection
"""

# ╔═╡ 49e7085c-3691-4182-9630-68dc9371ad18
md"""
## Goal of this practicum
"""

# ╔═╡ 84caa99f-32ab-433c-b620-f592d329d18d
md"""
In previous practicals, we have developed models to study phenomena and predict future behavior. We have also estimated the parameters associated with these models and we have also analyzed the sensitivity of the model predictions to changes in these parameters. We found that the mathematical structure of different models determines the sensitivity to errors in the parameters and errors in the model itself and we determined how these errors propagate through the model, allowing us to quantify the uncertainty in the predictions.
"""

# ╔═╡ 7cc1d1ce-a34a-4432-953d-4eb2df92696f
md"""
In this practical, we investigate how to make an objective choice between different candidate models by weighing the complexity of the models against the fit to the experimental data and the quality of the prediction. We will use two information criteria often used in practice to balance model quality and complexity: the Akaike information criterion (AIC) and the Bayesian information criterion (BIC).
"""

# ╔═╡ dcb51063-2b3f-4de9-90bc-0048e82cd4ab
md"""
In the Akaike information criterion, the fit or quality of the model (likelihood $L$) is compared against the number of model parameters ($k$), thus giving a measure of the balance between complexity and quality of the fit:

$AIC = 2k - 2\,\log(L)$
"""

# ╔═╡ 6690a844-8fd5-4ea1-9e62-a834bd454efb
md"""
The Bayesian information criterion gives similar information, but penalizes complexity more heavily:

$BIC = k\,\log(n) - 2\,\log(L)$

where $n$ is the number of data points considered.
"""

# ╔═╡ 61e5bead-6fc3-4c07-b2b4-832ce7a69198
md"""
In this notebook we will compare the different grass growth models and judge the quality of their fit to the calibration data set in order to select the simplest or least complex model that best represents the system.
"""

# ╔═╡ 04e95855-7c4a-4d2c-b836-c5dde291adad
md"""
## Grass growth models
"""

# ╔═╡ 4673bfdf-8f4a-42bd-a026-21a66d800f2b
md"""
In this notebook, three different models will be used, each modelling the yield of grass in a grassland:

- Logistic growth model: $\cfrac{dW}{dt} = \mu \left( 1 - \cfrac{W}{W_f} \right) W$
- Exponential growth model: $\cfrac{dW}{dt} = \mu \left( W_f - W \right)$
- Gompertz growth model: $\cfrac{dW}{dt} = \left( \mu - D \ln(W) \right) W$

with output $W$ the grass yield, and $W_f$, $\mu$ and $D$ parameters. The table below shows some typical values for the parameters:

|             | $\mu$      | $W_f$       | $D$          |
|:----------- |:----------:|:-----------:|:------------:|
| Logistic    |  0.07      | 10.0        |              |
| Exponential |  0.02      | 10.0        |              |
| Gompertz    |  0.09      |             | 0.04         |

We will use an initial condition of $W_0 = 2.0$ for each and a simulation time of $100$ days.
"""

# ╔═╡ 6e5f494e-ecfb-467b-a0f1-c85427da8215
md"""
In each of the three models we will use the following timespan:
"""

# ╔═╡ 0b3d35bb-c5b0-44b7-94b3-06fa571d339e
tspan = (0.0, 100.0)   # this will be the same for the three models

# ╔═╡ d1c5a659-70c4-4768-9542-a2988ea39c43
md"""
### The calibration data
"""

# ╔═╡ 583ed493-27e7-4fa7-b519-75cc71468fc6
md"""
Assume that the measured grass yields (of a certain plant type) over time are the following:
"""

# ╔═╡ d4ca6986-df24-48dd-8637-0974e0b16f2c
W_meas = [1.87, 2.45, 3.72, 4.32, 5.28, 7.01, 6.83, 8.62, 9.45, 10.31, 10.56, 11.72, 11.05, 11.53, 11.39, 11.7, 11.15, 11.49, 12.04, 11.95, 11.68]

# ╔═╡ 8215b6f0-017e-4a73-b4bf-aa61b8f3ca6f
md"""
They have been measured at the following corresponding time instances:
"""

# ╔═╡ 2ff5f464-38d7-4a90-ad31-d19d412781fc
t_meas = 0:5:100

# ╔═╡ b1d02f3c-f528-4297-a0c1-993cfe0d7c17
md"""
We can make a scatter plot of this data (including a title, a legend label, an X-axis label, X- and Y-axis limits) in the following way:
"""

# ╔═╡ 1499c727-ff19-4ef0-b467-bc71c1802979
scatter(t_meas, W_meas, title="Grass growth data",
                        label="Yield",
                        xlabel="t",
                        xlims=(0, 100),
                        ylims=(0, 14))

# ╔═╡ c692b627-8743-48f6-bb26-b7d9b3d3b9c8
md"""
### Logistic growth
"""

# ╔═╡ 676119ce-f4fe-41f9-8121-b2c21f0dd28c
md"""

$$\cfrac{dW}{dt} = \mu \left( 1 - \cfrac{W}{W_f} \right) W$$\
$W_0$ = 2.0, $\mu$ = 0.07 and $W_f$ = 10.0\
We will start by modelling our system and simulating using the aforementioned parameters values, initial condition and timespan in a way that we are familiar with.
"""

# ╔═╡ 5492210b-6ee6-4baf-b37e-d21358cdeb60
md"""
Implementation of the system:
"""

# ╔═╡ 18eece66-46fa-4458-aaed-f4c8fa002c20
growth_log = @reaction_network begin
	@species W(t)=2.0            # default initial condition
	@parameters μ=0.07 Wf=10.0   # default parameter values
	μ*(1-W/Wf), W --> 2W
end

# ╔═╡ 8bb0b24b-6cce-49cd-a625-5f375b92d9b7
md"""
Convert the *reaction model* to check that we work with the correct differential equation:
"""

# ╔═╡ ba181db8-d176-4b83-9168-d5939ffe9661
osys_log  = convert(ODESystem, growth_log)

# ╔═╡ cd8a7ba1-194a-4b42-8c87-2c9b1fe6b475
md"""
Setting initial conditions, timespan and parameter values:
"""

# ╔═╡ cbb2ca49-b019-495b-9310-83fcc00cad26
u0_log = [:W => 2.0]

# ╔═╡ be565a3c-31b6-4df1-b73b-08f308a8c09b
md"""
For the sake of clarity, we will use the variables `μ_log` and `Wf_log` to store the parameter values.
"""

# ╔═╡ f62806d1-77e1-470b-9711-33a924c788cc
μ_log = 0.07

# ╔═╡ 546ed163-26a7-4235-982f-7568ed609488
Wf_log = 10.0

# ╔═╡ 0da53fa2-5a42-46e6-8bd8-45d6aa903d46
params_log = [:μ => μ_log, :Wf => Wf_log]

# ╔═╡ b1298f40-4696-49d0-ac94-896e0cdbc996
md"""
Creating and solving the ODEProblem and plotting results:
"""

# ╔═╡ 270647d2-1371-4272-8bc1-3a6ad77bc716
oprob_log = ODEProblem(growth_log, u0_log, tspan, params_log);
# Also possible here if initial conditions and parameter values are defined in the catalyst model:
# oprob_log = ODEProblem(growth_mod_log, [], tspan, [])

# ╔═╡ ac235d86-1d93-4944-aa89-1b4fd38f0e6e
osol_log = solve(oprob_log, Tsit5(), saveat=0.5)

# ╔═╡ 3e13efa1-9bc6-456f-8e62-ecd3165e2a65
begin
	plot(osol_log, label="model", lw=2, ylabel="W")
	scatter!(t_meas, W_meas, title="Logistic growth model", label="data", xlabel="t",
	    xlims=(0, 100), ylims=(0, 14))
end

# ╔═╡ e2955af0-edf0-4702-9a8f-478141ffdc3b
md"""
We can see that the model does not predict well the data set for the considered parameter values. Thus we will use the data to both calibrate the model parameters and assess the quality of the fit.
"""

# ╔═╡ c73669c4-d7af-4877-b32d-6499f339e27a
md"""
### Parameter estimation
"""

# ╔═╡ 28e4a44f-e907-42ce-a748-6d21b60f0e33
md"""
We declare our Turing model function:
"""

# ╔═╡ 169a67ff-55fb-4d1d-b98b-126f4af47e77
@model function growth_log_fun(t_meas)
    σ_W ~ Exponential(10)
    W0 ~ LogNormal()
    μ ~ LogNormal()
    Wf ~ LogNormal()
	u0_log = [:W => W0]
	params_log = [:μ => μ, :Wf => Wf]
	oprob_log = ODEProblem(growth_log, u0_log, tspan, params_log)
    osol_log = solve(oprob_log, AutoTsit5(Rosenbrock23()), saveat=t_meas)
    W_s ~ MvNormal(osol_log[:W], σ_W^2 * I)
	return osol_log   # optionally, to be used with MCMC
end

# ╔═╡ bfb26b9d-6969-457c-8567-03422a7f2a93
md"""
!!! note
	As we will once again be calibrating our model, we use the more reliable Rosenbrock23 solver in our Turing model.
"""

# ╔═╡ aa3e553d-2731-42d6-b0e6-1821e4d7f4d4
md"""
We now provide the time measurements to the defined function (this results in the Turing model) and instantly condition the Turing model with the measurements of $W$:
"""

# ╔═╡ 51426716-03b8-4d54-9064-3943df282fa4
growth_log_cond_mod = growth_log_fun(t_meas) | (W_s = W_meas,)

# ╔═╡ f40873af-3ef9-420b-b0a8-ed9f56b17047
md"""
We are now ready to optimize the priors ($\sigma_W$, $W_0$, $\mu$ and $W_f$). This is done by calling the `optimize` function, providing the previously created object `growth_log_inf`, the method for estimating the parameters and (optionally) an algorithm (default: Nelder-Mead) to implement the method.
"""

# ╔═╡ 4536c6b7-f3bc-42e8-b044-0633f00b56bb
md"""
We will use the MLE (Maximum Likelihood Estimation) method here and store the optimization results in `results_log_mle`. If you get an error the first time, try running the optimization again.
"""

# ╔═╡ bc55fb07-1d84-4927-80d2-c1a012409400
results_log_mle = optimize(growth_log_cond_mod, MLE(), NelderMead())

# ╔═╡ c224f1b8-e5a1-4379-ae50-981ae89f9eaa
md"""
You can visualize a summary of the optimized parameters by piping them to `coeftable`:
"""

# ╔═╡ 4fe2b4cd-d0b0-489a-a867-14c58637327b
coeftable(results_log_mle)

# ╔═╡ 6cf5e2f2-344c-4a7a-8901-5a0ee9cb8689
md"""
You can obtain the actual optimized values using the function `coef` on the results object in conjunction by calling the parameters by name preceded by a colon. Here we assign the optimized parameter values to some suitable variable names:
"""

# ╔═╡ 551859cb-9cd1-4977-a528-b96423e83504
W0_opt_log = coef(results_log_mle)[:W0]

# ╔═╡ 955524a6-1a53-4d2d-ad16-a5a3a2a39b80
μ_opt_log = coef(results_log_mle)[:μ]

# ╔═╡ 118a79cc-afaf-4a99-8486-8e0e01fda40a
Wf_opt_log = coef(results_log_mle)[:Wf]

# ╔═╡ 8536cc9e-6cd4-4e38-9988-e45b06b75ea4
md"""
Now we can make a plot of $W$ simulated with the optimized initial condition and parameter values.
"""

# ╔═╡ 25db8d7b-e467-4a2f-b3f3-0ad7df64988a
md"""
Setting up initial condition with optimized initial condition:
"""

# ╔═╡ 01241971-cc9a-42d6-8e52-c840b91e6431
u0_opt_log = [:W => W0_opt_log]

# ╔═╡ e2db7492-4d68-4045-91b6-dca3aef8b514
md"""
Setting up parameter values with optimized parameter values:
"""

# ╔═╡ 7a890da6-7547-4854-b877-a5d935a9f2dd
params_opt_log = [:μ => μ_opt_log, :Wf => Wf_opt_log]

# ╔═╡ c95bcde6-ffb9-4ce6-9dfb-38c79cc2e4b8
md"""
Next, we create an ODEProblem and solve it:
"""

# ╔═╡ bf90a66a-8505-4b33-b821-d708a3b7f8b0
oprob_opt_log = ODEProblem(growth_log, u0_opt_log, tspan, params_opt_log)

# ╔═╡ 7af67f91-58e2-429d-a88f-b27a0258a805
osol_opt_log = solve(oprob_opt_log, Tsit5(), saveat=0.5);

# ╔═╡ 3ee0d4de-fdf5-40f0-a849-63b0759eb44b
md"""
Finally, we plot $W$ simulated with the optimized initial value and parameter values together with the measured data that was used to find the optimized values.
"""

# ╔═╡ 2c048fc0-5908-4f01-a9a5-480e12e098cc
begin
	plot(osol_opt_log, label="model", xlabel="t", ylabel="W",
		xlims=(0, 100), ylims=(0, 14), lw=2.0, title="Calibrated logistic growth model")
	scatter!(t_meas, W_meas, label="data")
end

# ╔═╡ 19ed4aa9-776e-4a80-b180-65d32f3f9f26
md"""
We can extract from the calibration results the log-likelihood or quality of the fit:
"""

# ╔═╡ 12a11f90-fa55-47f2-a75f-853cd54cbeab
L_log = results_log_mle.lp

# ╔═╡ f6b44166-3ba9-4770-80db-303bd28112cf
md"""
## Model selection criteria
"""

# ╔═╡ 02037c91-448e-4714-b932-86f4df631907
md"""
### Akaike information criterion
"""

# ╔═╡ cfea46a6-636c-4fcd-a04b-714e7a567cb3
md"""
To calculate the AIC, we can implement a function that uses the information from the calibration:
"""

# ╔═╡ cb9384a3-96b0-4694-a218-25ebf96187a9
function AIC(results, measurements)
	L = results.lp
	k = length(results.values)
	n = length(measurements)
	return 2k - 2L 			# L = log-likelihood
end

# ╔═╡ 73917400-43e0-4b9f-9448-9bbe02daebb9
md"""
This function uses the results from the calibration, from where we can extract as well the number of calibrated parameters, which includes the estimated prediction error:
"""

# ╔═╡ afba526d-0f71-4ba5-9a84-3c24804f959a
k_log = length(results_log_mle.values) 	# alternative: length(coef(results_log_mle))

# ╔═╡ 903380ce-1320-464e-a703-1222b4f381b1
md"""
The AIC will use this to balance the complexity with the quality of the fit. For the logistic model:
"""

# ╔═╡ 22bd4450-6c6d-4a4d-aa1d-a080d3b06eac
AIC_log = AIC(results_log_mle, W_meas)

# ╔═╡ 34f8cb7d-7482-4883-87e0-5c68be7ec0e9
md"""
### Bayesian information criterion
"""

# ╔═╡ 4d976b86-74d6-4bb5-af17-8a6bd658abe2
md"""
We can also calculate the BIC in a similar way to the AIC:
"""

# ╔═╡ ba7974ec-3965-4c92-ab37-08cdf5ab8125
function BIC(results, measurements)
	L = results.lp
	k = length(results.values)
	n = length(measurements)
	return k*log(n) - 2L 	# L = log-likelihood
end

# ╔═╡ 5e6837d5-35e9-434c-9337-7c90daec9c33
md"""
The BIC will additionally use the length of the data set for the complexity penalty term:
"""

# ╔═╡ b97ca16c-e6c3-4840-a855-1d29140eff59
n = length(W_meas)

# ╔═╡ 1af53492-f363-4d39-bfaf-2f668a441725
BIC_log = BIC(results_log_mle, W_meas)

# ╔═╡ e176e6ce-33da-4f3b-860e-4db99a931079
md"""
!!! question
	What conclusions can we extract from a comparison of the AIC or BIC for different models?
"""

# ╔═╡ 770e70d0-8a69-4585-a423-41335d453dcd
md"""
Conclusions:
- Lower values are better. For the same value of $L$, a simpler model would be preferred.
- BIC seems to penalize more complex models than AIC for the same values of $L$ and $k$.
"""

# ╔═╡ f4f3b1fb-165d-4f77-88f1-8c1a79f3fb42
md"""
### The posterior model probability
"""

# ╔═╡ 4a86d852-634c-43f5-90e5-3e675aeb2d69
md"""
We can use the AIC to compute the posterior probabilities of the different candidate models:

$P(M_i|D) \propto \exp(-AIC(M_i)/2)$
"""

# ╔═╡ 8a43d4b8-1902-4e02-9cd7-c36131ec9528
md"""
The following function will use the supplied AIC of several models to compute the normalized posterior probability that the model is the "true model", explaining the considered data set: 
"""

# ╔═╡ d149d829-7d3d-4f7a-bdad-b1debb2c5149
function posterior(AICs) 	# AICs vector of AIC values
	AICmin = minimum(AICs)
	posterior = zeros(length(AICs))
	
	for i in eachindex(1:length(AICs))
		posterior[i] = exp((AICmin-AICs[i])/2)
	end
	
	return round.(posterior/sum(posterior); digits=3)	# normalized sum
end

# ╔═╡ e2bec1ab-f161-45b6-a67e-1f4f197ae685
posterior([AIC_log])

# ╔═╡ 5db6b1a5-ed2e-4812-94db-37cfffb46227
md"""
!!! note
	This function will be used to compare the different candidate models (more than one).
"""

# ╔═╡ bc8d315a-ddfb-4c82-84d0-90a4fb05256c
md"""
### Least squares model fitting
"""

# ╔═╡ 51b348ef-0270-41d0-bbcb-a27eb10ff5d2
md"""
The Akaike information criterion can be reformulated in terms of **least squares** if we assume that the model residuals are normally (and independently) distributed with zero mean, giving rise to:

$AIC = 2k + n \log{\bigg(\frac{SSR}{n}}\bigg)$

where $SSR$ is the **squared sum of the model's residuals**. For small data sets, a correction is done:

$AIC_c = 2k + n \log{\bigg(\frac{SSR}{n}}\bigg) + \frac{2k(k+1)}{n-k-1}$

When the number of observations is large enough, the corrected $AIC_c$ and $AIC$ are identical.

The Bayesian information criterion can also be expressed in terms of the residuals:

$BIC = k\log(n) + n \log{\bigg(\frac{SSR}{n}}\bigg)$

Both criteria are implemented below and can be used to compare the fitness of different models. 
"""

# ╔═╡ 7bf8fa3a-ce98-4f40-8d6c-4e84b79ae1c0
function AIC_LS(SSR, n, k) 
	if n > 40
		return 2k + n*log(SSR/n) 
	else
		return 2k + n*log(SSR/n) + 2k*(k+1)/(n-k-1)
	end
end

# ╔═╡ 498e0b9b-17a9-4059-9196-c6c4d86f926e
function BIC_LS(SSR, n, k) 
	return  k*log(n) + n*log(SSR/n)
end

# ╔═╡ b57131ad-e9d6-4318-bca7-c59c5d8da4c3
md"""
We can thus obtain the squared sum of residuals from the calibrated model prediction and the data:
"""

# ╔═╡ f2fe206b-2832-42f6-84bd-20c3838fcb0b
function SSR(y_pred, y_data)
	return sum((y_pred - y_data).^2) 	# squared sum of residuals
end

# ╔═╡ d51fff9b-17d3-4a0d-aa5e-01d58821a246
md"""
We can now calculate the SSR and alternative AIC and BIC forms for the logistic model:
"""

# ╔═╡ 6d95d3b1-d964-4b61-adb6-f44875e23ecb
begin
	W_log = solve(oprob_opt_log, Tsit5(), saveat=t_meas)[:W] # model prediction
	SSR_log = SSR(W_log, W_meas)
	AIC_LS_log = AIC_LS(SSR_log, n, k_log)
	BIC_LS_log = BIC_LS(SSR_log, n, k_log)
end;

# ╔═╡ b1783355-0221-4f2d-a050-932b40920eab
AIC_log, BIC_log

# ╔═╡ 3c69edb3-ba7d-4ee0-996a-40fd985e57d3
AIC_LS_log, BIC_LS_log

# ╔═╡ 982e1de1-0820-42fb-8abb-31f1c99def46
md"""
!!! note
	See the exercises below to apply the different criteria for model selection to the other models.
"""

# ╔═╡ 79f0f1dd-850f-4dfd-b895-ff41a2d8adb8
md"""
## Exercises
"""

# ╔═╡ fa2270d0-2548-409c-a23f-4369d8bce8ec
md"""
### Exercise 1 - Compare the logistic and exponential models
"""

# ╔═╡ 545a317d-8bca-4b6c-9915-84b0b4ffbed5
md"""
Calibrate the initial condition and both parameters of the exponential growth model. Use the values mentioned in the Table as initials values for the optimization of the parameters. Then compare the fit to that of the logistic model by plotting both predictions in the same figure.
"""

# ╔═╡ 794b93c9-0a0a-4e77-b13c-0e2c06a9a0ec
md"""
$$\cfrac{dW}{dt} = \mu \left( W_f - W \right)$$\
$W_0$ = 2.0, $\mu$ = 0.02 and $W_f$ = 10.0
"""

# ╔═╡ cfca0095-81e7-4c88-84bd-773eb5af7aa4
growth_exp = @reaction_network begin
    μ*Wf, 0 --> W
    μ, W --> 0
end

# ╔═╡ ba9162dd-27f4-42b3-a6b5-ca9531614d5e
md"""
Use the same measurement data (`W_meas`, `t_meas`) as before.
"""

# ╔═╡ b1f03df5-4353-45ed-9d91-9e2c9960f706
md"""
Declare the Turing model function.
"""

# ╔═╡ 89bf91c6-117c-4647-bfb9-6fc9b8dcfb5f
@model function growth_exp_fun(t_meas)
    σ_W ~ Exponential(10)
    W0 ~ LogNormal()
    μ ~ LogNormal()
    Wf ~ LogNormal()
	u0_exp = [:W => W0]
	params_exp = [:μ => μ, :Wf => Wf]
	oprob_exp = ODEProblem(growth_exp, u0_exp, tspan, params_exp)
    osol_exp = solve(oprob_exp, AutoTsit5(Rosenbrock23()), saveat=t_meas)
    W_s ~ MvNormal(osol_exp[:W], σ_W^2 * I)
end

# ╔═╡ b9a8c0cb-fe38-448e-aa4c-e9422f24a4a4
md"""
Provide the time measurements to the defined function (this results in the Turing model) and instantly condition the Turing model with the measurements of $W$:
"""

# ╔═╡ 6b5fc852-1d40-45ce-a2e2-d22a6cc03156
growth_exp_cond_mod = growth_exp_fun(t_meas) | (W_s = W_meas,)

# ╔═╡ 29591722-bb44-47d0-a336-b2151502fd96
md"""
Optimize the priors ($\sigma_W$, $W_0$, $\mu$ and $W_f$). Do this with both the `MLE` and `MAP` methods and the Nelder-Mead algorithm. Store the optimization results in `results_exp_mle` and `results_exp_map`.
"""

# ╔═╡ ec032fb5-76a7-4a0c-a8be-b0f6bd166580
results_exp_mle = optimize(growth_exp_cond_mod, MLE(), NelderMead())

# ╔═╡ d44ed37a-d022-457c-befe-01bae5cd0ed9
md"""
Visualize a summary of the optimized parameters.
"""

# ╔═╡ 5d3f6390-1f71-4543-ae4e-d72f28b94f40
coeftable(results_exp_mle)

# ╔═╡ 9eb6e432-3f41-4b5f-9194-241e6ab3a039
md"""
Get the optimized values and assign them to `W0_opt_exp`, `μ_opt_exp` and `Wf_opt_exp`.
"""

# ╔═╡ cc0fc9ba-5564-4104-b0b3-a10ae39ef933
W0_opt_exp = coef(results_exp_mle)[:W0]

# ╔═╡ 6083305f-c838-4770-95c5-258fd7becbcb
μ_opt_exp = coef(results_exp_mle)[:μ]

# ╔═╡ 9777f3e8-785d-4ff5-a6c6-80a68a9c53ba
Wf_opt_exp = coef(results_exp_mle)[:Wf]

# ╔═╡ e80de22e-d5e6-4af7-af0e-feaf9ef660bf
md"""
Make a plot of $W$ simulated with the optimized initial condition and parameter values.
"""

# ╔═╡ 0020352d-52dc-4ca4-90dd-f53865306f9e
md"
Set up initial condition with optimized initial condition:
"

# ╔═╡ e58c49dd-5365-4dcb-9baa-c7a6048eb0ed
u0_opt_exp = [:W => W0_opt_exp]

# ╔═╡ 03a35cce-cacf-4b11-a0b8-3a26ad93c227
md"""
Set up parameter values with optimized parameter values:
"""

# ╔═╡ b3556b32-4665-439e-8751-a3ffa8ba1907
params_opt_exp = [:μ => μ_opt_exp, :Wf => Wf_opt_exp]

# ╔═╡ ae030bbd-6416-44b8-ba9e-3fe9be5c9c52
md"""
Create an ODEProblem and solve it. Solve it using `Tsit5()` and `saveat=0.5`.
"""

# ╔═╡ bfd9d570-7b12-4bbf-aaa9-449f2377a294
oprob_opt_exp = ODEProblem(growth_exp, u0_opt_exp, tspan, params_opt_exp)

# ╔═╡ f9e1ca7a-2a27-46aa-82a2-5e06882a9fad
osol_opt_exp = solve(oprob_opt_exp, Tsit5(), saveat=0.5);

# ╔═╡ 0e62cf17-f885-4c38-8df2-39d1719bfa04
md"""
Plot now $W$ simulated with the optimized initial value and parameter values of both logistic and exponential models together with the measured data that was used to find the optimized values.
"""

# ╔═╡ fe3c059f-cb4f-4ddf-ad54-c08aa70bc740
# Uncomment and complete the instruction
# begin
# plot()
# missing
# missing
# missing
# title!("Comparison logistic vs. exponential growth")
# end

# ╔═╡ fcfd18a2-154f-46fd-b4a5-01e406b500f1
md"""
!!! question
	By looking at the figure, how can you decide which candidate model is better?
"""

# ╔═╡ 78d5ad0e-318f-4af1-883b-b887f8d72591
md"""
- Answer: missing
"""

# ╔═╡ 58006718-5e76-4c6b-a968-aa5372bb4a15
md"""
Compare now the fit of both models by applying both the AIC and BIC criteria.
"""

# ╔═╡ 7606cd16-7463-462e-83e6-5adaea4632a2
md"""
Extract the log-probability and number of parameters from the calibration results of the exponential: 
"""

# ╔═╡ 2379169c-7f2c-415c-8955-e8b06fc9edf5
L_exp = results_exp_mle.lp

# ╔═╡ 577bb77e-d8c6-4121-84e6-1d3c000ce7b9
k_exp = length(coef(results_exp_mle))

# ╔═╡ ac2d3b13-54d2-4217-8b18-7f9e276f1b6f
md"""
Calculate the AIC and BIC for the exponential model:
"""

# ╔═╡ 2174d25b-b603-4d7b-b396-7ea00780143c
AIC_exp = AIC(results_exp_mle, W_meas)

# ╔═╡ 0e8a6614-3164-4eab-b28b-fde3a328cbc3
BIC_exp = BIC(results_exp_mle, W_meas)

# ╔═╡ f2805dab-cac4-4a31-8418-6763c14d6415
L_log, L_exp

# ╔═╡ 41d8cfe2-a997-4731-971e-3e0c6166177e
AIC_log, AIC_exp

# ╔═╡ 2b2110eb-b867-47a1-8238-3d43a64851bb
BIC_log, BIC_exp

# ╔═╡ 4fb53b80-4ea1-47a0-871a-adde769f6d18
md"""
!!! question
	Draw your conclusions.
"""

# ╔═╡ 22b58c97-1aa2-49c4-8a0a-488c63014a90
md"
- missing
"

# ╔═╡ 41abf8c4-e67e-4f66-a38a-7204a878d98d
md"""
### Exercise 2 - Comparison of the three models
"""

# ╔═╡ b02787ec-40eb-4000-9684-e69cbd912c6c
md"""
Perform the calibration of the Gompertz model and compare its fitness to the other two candidates.
"""

# ╔═╡ e34fed6f-6b58-4642-84eb-167467881bb2
md"""
$$\cfrac{dW}{dt} = \left( \mu - D \ln(W) \right) W$$\
$W_0$ = 2.0, $\mu$ = 0.09 and $D$ = 0.04.
"""

# ╔═╡ 68895842-0a1b-4236-b159-19a2805ea44d
growth_gom = @reaction_network begin
	μ-D*log(W), W --> 2W
end

# ╔═╡ ec12d3c4-6b48-4f3a-b463-0fc67d3adff3
md"""
Declare the Turing model. Take the same priors as before.
"""
# Take for $\sigma_W$ and $W_0$ the same priors (and distributions) as before, but take for $\mu$ a Uniform prior distribution in the range $[0, 2]$ and the same for $D$ but in the range $[0, 1]$.

# ╔═╡ 1f8b10ba-f380-492e-800c-2673774cb25c
@model function growth_gom_fun(t_meas)
    σ_W ~ Exponential(10)
    W0 ~ LogNormal()
    μ ~ LogNormal()
    D ~ LogNormal()
	u0_gom = [:W => W0]
	params_gom = [:μ => μ, :D => D]
	oprob_gom = ODEProblem(growth_gom, u0_gom, tspan, params_gom)
    osol_gom = solve(oprob_gom, AutoTsit5(Rosenbrock23()), saveat=t_meas)
    W_s ~ MvNormal(osol_gom[:W], σ_W^2 * I)
end

# ╔═╡ 82bee644-7b9a-42fa-a32f-76b3eb5038b5
md"""
Provide the time measurements to the defined function (this results in the Turing model) and instantly condition the Turing model with the measurements of $W$:
"""

# ╔═╡ e481b87c-c6f3-48c1-bfff-5d76de0dce36
growth_gom_cond_mod = growth_gom_fun(t_meas) | (W_s = W_meas,)

# ╔═╡ a438470b-e037-4b65-b6b0-fb9c202e74a4
md"""
Optimize the priors ($\sigma_W$, $W_0$, $\mu$ and $D$). Do this now with `MAP` method and Nelder-Mead. Store the optimization results in `results_gom_map`.
"""

# ╔═╡ f4c3dcd0-70ec-49ab-bbad-04b550da481c
results_gom_mle = optimize(growth_gom_cond_mod, MLE(), NelderMead())

# ╔═╡ df166b45-3c28-4bad-95b8-f43429a58fd0
md"""
Visualize a summary of the optimized parameters.
"""

# ╔═╡ a2738236-86a1-419a-9617-b1229f7c9240
coeftable(results_gom_mle)

# ╔═╡ 59f78a22-cac1-49a7-b3e6-9d74b694be64
md"""
Get the optimized values and assign them to `W0_opt_gom`, `μ_opt_gom` and `D_opt_gom`.
"""

# ╔═╡ 7c351d72-c5ce-4e5a-b2e5-87ff9dbc70e5
W0_opt_gom = coef(results_gom_mle)[:W0]

# ╔═╡ d7312308-7c95-4c14-9bd6-56a7ac99d49a
μ_opt_gom = coef(results_gom_mle)[:μ]

# ╔═╡ 8f8b23b3-2e33-46e4-917a-294934fa090e
D_opt_gom = coef(results_gom_mle)[:D]

# ╔═╡ 4ff1079f-8dfb-4069-9a98-0d263eba9920
md"""
Make a plot of $W$ simulated with the optimized initial condition and parameter values.
"""

# ╔═╡ 9edd77c3-eb1e-4be7-bb6c-83803e675102
md"""
Set up initial condition with optimized initial condition:
"""

# ╔═╡ c02a4f18-d51d-4954-b7a4-a4ae64b10fb7
u0_opt_gom = [:W => W0_opt_gom]

# ╔═╡ 060b5869-881e-4d57-8357-cff95912ba16
md"""
Set up parameter values with optimized parameter values:
"""

# ╔═╡ 0d225f18-bdf6-49d1-a827-2bc9551c158d
params_opt_gom = [:μ => μ_opt_gom, :D => D_opt_gom]

# ╔═╡ 5d27e199-ad09-454a-8966-1531371ce67b
md"""
Create an ODEProblem and solve it. Use the solver `Tsit5()` and `saveat=0.5`.
"""

# ╔═╡ 79c4439b-df8f-467a-b791-5581bd564996
oprob_opt_gom = ODEProblem(growth_gom, u0_opt_gom, tspan, params_opt_gom)

# ╔═╡ bcbe7409-e969-4c7a-9417-195d406617bf
osol_opt_gom = solve(oprob_opt_gom, Tsit5(), saveat=0.5);

# ╔═╡ 6dbc81c8-4bc7-4499-a9d1-9719df368176
md"""
Finally, we plot $W$ simulated with the optimized initial value and parameter values together with the measured data that was used to find the optimized values.
"""

# ╔═╡ ef5270de-da37-4a01-8734-dac1fd190a05
# Uncomment and complete the instruction
# begin
# plot()
# missing
# missing
# missing
# title!("Comparison logistic vs. exponential growth")
# end

# ╔═╡ b9030c3f-2cbb-44bf-90b7-23d116de4001
# L_gom = missing

# ╔═╡ 2e684fd0-8410-4037-a6a6-2c472d6592ab
# k_gom = missing

# ╔═╡ 79797d6d-9ca8-4f4f-9694-8045a8cf8cbb
# AIC_gom = missing

# ╔═╡ 2cd566bb-3088-43f8-8926-9952342b5e5c
# BIC_gom = missing

# ╔═╡ 750edde4-16b5-4025-82e0-078c0eacaca9
# AIC_log, AIC_exp, AIC_gom

# ╔═╡ 1bbbd1ae-3540-44da-bc53-c8182a85b611
# BIC_log, BIC_exp, BIC_gom

# ╔═╡ 3afb2a69-32d6-42c2-8c1e-51e64b83a6b4
md"""
!!! question
	Draw your conclusions.
"""

# ╔═╡ 01ddf3dd-5199-48cf-ab13-75686278f3dc
md"""
- Answer: missing
"""

# ╔═╡ 36f71777-ad33-4b27-9a0e-ff03b22ad79e
md"""
You can use the following graph with all the information calculated so far for your conclusions.
"""

# ╔═╡ d3b23128-2e7a-4246-a1f2-7325e0c2c432
# plot(
# 	bar(1:3, [AIC_log, AIC_exp, AIC_gom], title="AIC", ylims=(0, 50)), 
# 	bar(1:3, [BIC_log, BIC_exp, BIC_gom], title="BIC", ylims=(0, 50)), 
# 	bar(1:3, [L_log, L_exp, L_gom], title="Log-probability", ylims=(-20, 0)), 
# 	bar(1:3, [k_log, k_exp, k_gom], title="no. parameters", ylims=(0, 8)), 
# 	xticks=(1:3, ["Logistic", "Exponential", "Gompertz"]),
# 	legend=:none
# )

# ╔═╡ a9e5b236-67d5-400b-ba9a-8b565f6feeb1
md"""
### Exercise 3 - Calculation of the posterior probabilities
"""

# ╔═╡ 3bf0c957-f08b-4a45-90e0-40d56330c8e0
md"""
Use the above implemented function `posterior` to calculate the posterior model probabilities.
"""

# ╔═╡ 947e26fe-2010-48ad-b927-d4e5186f21e0
# posteriors = missing

# ╔═╡ 3506765b-a6cf-41ca-9201-fcc650c6f56d
# posteriors

# ╔═╡ 72dc27a0-0e61-4073-843a-67e1bfd1427a
# md"""
# We can summarize all calculated criteria so far in the following table:

# | Model | k | Log(L) | AIC | BIC | $P(M_i\|D)$ |
# |:---|:---|:---|:---|:---|:---|
# | Logistic | $k_log | $(round(L_log;digits=3)) | $(round(AIC_log;digits=3)) | $(round(BIC_log;digits=3)) | $(posteriors[1]) |
# | Exponential | $k_exp | $(round(L_exp;digits=3)) | $(round(AIC_exp;digits=3)) | $(round(BIC_exp;digits=3)) | $(posteriors[2]) |
# | Gompertz | $k_gom | $(round(L_gom;digits=3)) | $(round(AIC_gom;digits=3)) | $(round(BIC_gom;digits=3)) | $(posteriors[3]) |

# """

# ╔═╡ 6240553c-0c7c-45f4-b190-ebacf09ab632
md"""
!!! question
	Draw your conclusions. Does the posterior probability give the same raking as the other criteria?
"""

# ╔═╡ 12f64333-36d0-4e2a-9061-e3dcc7a4ae96
md"
- missing
"

# ╔═╡ 95ba48fe-f430-4f47-ab04-f5896d75343b
md"""
### Exercise 4 - Comparison with least squares
"""

# ╔═╡ ef293617-7a7a-4597-99d2-2c39c543c296
md"""
Repeat below the comparison to least squares for the exponential and Gompertz models.
"""

# ╔═╡ 8bbed35b-fec5-4682-8b45-99fd3bd55bf2
begin
	# Uncomment and complete the instruction
	# W_exp = missing
	# SSR_exp = missing
	# AIC_LS_exp = missing
	# BIC_LS_exp = missing
end;

# ╔═╡ 414ebb03-4fb7-4591-90e9-9ce3cf11e93d
begin
	# Uncomment and complete the instruction
	# W_gom = missing
	# SSR_gom = missing
	# AIC_LS_gom = missing
	# BIC_LS_gom = missing
end;

# ╔═╡ 7a633bd6-10f1-4755-94ac-40d1dc9580a7
# AIC_LS_log, AIC_LS_exp, AIC_LS_gom

# ╔═╡ aef172f1-8563-4e51-8a88-ec6b764a735e
# BIC_LS_log, BIC_LS_exp, BIC_LS_gom

# ╔═╡ 091e978a-fd89-455a-9b97-f3088d580683
# posterior([AIC_LS_log, AIC_LS_exp, AIC_LS_gom])

# ╔═╡ d3d6e099-66b0-444e-8142-ccd1718ac9ab
# plot(
# 	bar(1:3, [AIC_LS_log, AIC_LS_exp, AIC_LS_gom], title="AIC", ylims=(-40, 0)), 
# 	bar(1:3, [BIC_LS_log, BIC_LS_exp, BIC_LS_gom], title="BIC", ylims=(-40, 0)), 
# 	bar(1:3, [SSR_log, SSR_exp, SSR_gom], title="SSR", ylims=(0, 8)), 
# 	bar(1:3, [k_log, k_exp, k_gom], title="no. parameters", ylims=(0, 8)), 
# 	xticks=(1:3, ["Logistic", "Exponential", "Gompertz"]),
# 	legend=:none,
# 	suptitle="(Least squares)"
# )

# ╔═╡ 2de0f619-63e0-444b-916b-316e4ebfa13e
md"""
!!! question
	Draw your conclusions. Do the SSR and alternative AIC and BIC provide the same model ranking?
"""

# ╔═╡ ef5329d2-ac63-4b49-a4cb-e8d6281a3184
md"""
- Answer: missing
"""

# ╔═╡ 0b8e6775-7562-4c93-b64d-cb694531f250
md"""
## Additional exercises
"""

# ╔═╡ c99ac3a9-1d5e-43a3-a4ef-a22147187bb3
md"""
### 1. MAP estimation
"""

# ╔═╡ b620006e-d6ac-4038-960f-8d0d93c2de00
md"""
We can repeat the calibration and take into account the priors to obtain the MAP estimation.
"""

# ╔═╡ 00dbb9d0-2c9d-4559-91d8-1a54171cc257
results_log_map = optimize(growth_log_cond_mod, MAP(), NelderMead())

# ╔═╡ 3d5e9bb7-5e65-4de5-a521-6b1f1cb872b3
coeftable(results_log_map)

# ╔═╡ 3a3cd103-78b1-481a-8746-749d402346d1
md"""
!!! question
	How will this affect the different criteria for the model selection? How is log(L) compared to MLE?
"""

# ╔═╡ b9986617-923e-4744-9d44-1935770edbd7
md"""
### 2 - Watanabe-Akaike information criterion (WAIC)
"""

# ╔═╡ 6e0c869d-6f90-4104-8589-41c2fa1b7342
md"""
The AIC and BIC are easy to compute but do not take into account the uncertainty in the predictions for the assessment of the model. The more complex Widely Applicable Information Criterion (WAIC) or Watanabe-Akaike information criterion takes samples from the posterior distribution and provides a measure of uncertainty for each observation, which can be used for model assessment.
"""

# ╔═╡ 6ca83977-cb8c-4465-ab2e-c91971f16157
md"""
We can generate new samples from the posterior distribution with MCMC.
"""

# ╔═╡ fb76009b-3aaa-4eb5-8d69-45372a23af00
N = 200

# ╔═╡ b9c59973-512b-48eb-a2b1-831449e120a4
results_log_nuts = sample(growth_log_cond_mod, NUTS(), N)

# ╔═╡ cfeb28b8-2dd7-4c1d-8069-e8fe00d90567
plot(results_log_nuts)

# ╔═╡ 56784bdf-2a78-4986-9dd2-790a2a2ceea9
md"""
The log-pointwise-predictive-density (lppd) is the sum of the log-likelihood of all observations:
"""

# ╔═╡ bfb99943-1637-4dd6-adee-8a9e79274672
lppd_log = sum(results_log_nuts.value[:, :lp])

# ╔═╡ a08fb39b-3911-450e-8377-cc94398825a9
md"""
The second part of WAIC is the variance of the log-likelihood of each observation, also called the effective number of parameters, $p_{WAIC}$, considered here as a penalty term, similarly to AIC and BIC:
"""

# ╔═╡ 2575199f-a3ec-470f-87ea-d62be28c42dc
pWAIC_log = sum(results_log_nuts.value[:, :lp].^2)/N - (lppd_log/N)^2

# ╔═╡ 6391d316-8826-4f19-9a2e-7065db787bee
sum((results_log_nuts.value[:, :lp] .- lppd_log/N).^2)/N

# ╔═╡ 2879ff01-860f-4fb1-b52d-3a43e57ccb23
md"""
Finally, WAIC is defined as:

$WAIC = -2(\text{lppd} - p_{WAIC})$
"""

# ╔═╡ 3b78c1ff-a2f6-4c1c-8124-5fe58c5a312d
WAIC_log = -2*(lppd_log - pWAIC_log)

# ╔═╡ 75e221b4-fdfe-4ef0-a8ea-1d75700fc1e7
md"""
We repeat the calculation for the exponential and Gompertz models below.
"""

# ╔═╡ 0f228304-de82-48f2-8e43-6ef81f7a9830
results_exp_nuts = sample(growth_exp_cond_mod, NUTS(), N)

# ╔═╡ 9843e282-44b3-40d0-b2b5-506144c22b27
lppd_exp = sum(results_exp_nuts.value[:, :lp])

# ╔═╡ eb9fe668-0a0f-4a8f-a793-ba545cb02e85
pWAIC_exp = (sum(results_exp_nuts.value[:, :lp].^2)/N - (lppd_exp/N)^2)

# ╔═╡ ddf9ae6d-0abe-4c8c-9646-024dd9a41e02
WAIC_exp = -2*(lppd_exp - pWAIC_exp)

# ╔═╡ 8330fce5-13c9-447e-8368-01b7c7ebc68f
results_gom_nuts = sample(growth_gom_cond_mod, NUTS(), N)

# ╔═╡ 0d2cab0c-04a2-4c8b-bb4f-2115ae0cf871
lppd_gom = sum(results_gom_nuts.value[:, :lp])

# ╔═╡ 6900e703-de7d-4ea3-8ac5-0b38b21b083d
pWAIC_gom = sum(results_gom_nuts.value[:, :lp].^2)/N - (lppd_gom/N)^2

# ╔═╡ b353eab2-d78f-4c5e-8e38-bd57813b6816
WAIC_gom = -2*(lppd_gom - pWAIC_gom)

# ╔═╡ 48dd67e7-363a-4e7c-ad81-8f70fc6db847
[WAIC_log, WAIC_exp, WAIC_gom]

# ╔═╡ a7fa7870-f202-4694-a9c1-7021ed242966
posterior([WAIC_log, WAIC_exp, WAIC_gom])

# ╔═╡ 1699d3ae-d77a-4272-8c4b-ddd4eadd503f
md"""
!!! question
	Why is the WAIC criterion significantly better for model selection despite its complexity?
"""

# ╔═╡ 72878862-dcc0-4c79-9c91-7c0bc44359d3
md"""
#### References
1. [https://en.wikipedia.org/wiki/Watanabe%E2%80%93Akaike_information_criterion](https://en.wikipedia.org/wiki/Watanabe%E2%80%93Akaike_information_criterion)
2. [https://civil.colorado.edu/~balajir/CVEN6833/bayes-resources/RM-StatRethink-Bayes.pdf](https://civil.colorado.edu/~balajir/CVEN6833/bayes-resources/RM-StatRethink-Bayes.pdf)
"""

# ╔═╡ Cell order:
# ╠═5f4fea06-0632-11ef-102e-21f5606d2056
# ╠═489b5399-fe4c-481a-834f-0101bbe28cea
# ╠═e1e7bc8e-7264-4cbc-98d2-aa73679fa2df
# ╠═ab4a1911-461e-4878-9258-931fc2f1ea06
# ╠═cdcebbb1-40e0-457f-a6ec-b769f6b1f2e9
# ╠═22489bd4-ab64-4bf9-ad03-5372ea273935
# ╠═4ad72a6a-f541-4d0a-8392-56369e29ac96
# ╠═70590d28-058a-44df-8cf6-092d8c87438c
# ╟─427b509f-d08e-4d93-99ba-a79f9c244b28
# ╟─49e7085c-3691-4182-9630-68dc9371ad18
# ╟─84caa99f-32ab-433c-b620-f592d329d18d
# ╟─7cc1d1ce-a34a-4432-953d-4eb2df92696f
# ╟─dcb51063-2b3f-4de9-90bc-0048e82cd4ab
# ╟─6690a844-8fd5-4ea1-9e62-a834bd454efb
# ╟─61e5bead-6fc3-4c07-b2b4-832ce7a69198
# ╟─04e95855-7c4a-4d2c-b836-c5dde291adad
# ╟─4673bfdf-8f4a-42bd-a026-21a66d800f2b
# ╟─6e5f494e-ecfb-467b-a0f1-c85427da8215
# ╠═0b3d35bb-c5b0-44b7-94b3-06fa571d339e
# ╟─d1c5a659-70c4-4768-9542-a2988ea39c43
# ╟─583ed493-27e7-4fa7-b519-75cc71468fc6
# ╠═d4ca6986-df24-48dd-8637-0974e0b16f2c
# ╟─8215b6f0-017e-4a73-b4bf-aa61b8f3ca6f
# ╠═2ff5f464-38d7-4a90-ad31-d19d412781fc
# ╟─b1d02f3c-f528-4297-a0c1-993cfe0d7c17
# ╠═1499c727-ff19-4ef0-b467-bc71c1802979
# ╟─c692b627-8743-48f6-bb26-b7d9b3d3b9c8
# ╟─676119ce-f4fe-41f9-8121-b2c21f0dd28c
# ╟─5492210b-6ee6-4baf-b37e-d21358cdeb60
# ╠═18eece66-46fa-4458-aaed-f4c8fa002c20
# ╟─8bb0b24b-6cce-49cd-a625-5f375b92d9b7
# ╠═ba181db8-d176-4b83-9168-d5939ffe9661
# ╟─cd8a7ba1-194a-4b42-8c87-2c9b1fe6b475
# ╠═cbb2ca49-b019-495b-9310-83fcc00cad26
# ╟─be565a3c-31b6-4df1-b73b-08f308a8c09b
# ╠═f62806d1-77e1-470b-9711-33a924c788cc
# ╠═546ed163-26a7-4235-982f-7568ed609488
# ╠═0da53fa2-5a42-46e6-8bd8-45d6aa903d46
# ╟─b1298f40-4696-49d0-ac94-896e0cdbc996
# ╠═270647d2-1371-4272-8bc1-3a6ad77bc716
# ╠═ac235d86-1d93-4944-aa89-1b4fd38f0e6e
# ╠═3e13efa1-9bc6-456f-8e62-ecd3165e2a65
# ╟─e2955af0-edf0-4702-9a8f-478141ffdc3b
# ╟─c73669c4-d7af-4877-b32d-6499f339e27a
# ╟─28e4a44f-e907-42ce-a748-6d21b60f0e33
# ╠═169a67ff-55fb-4d1d-b98b-126f4af47e77
# ╟─bfb26b9d-6969-457c-8567-03422a7f2a93
# ╟─aa3e553d-2731-42d6-b0e6-1821e4d7f4d4
# ╠═51426716-03b8-4d54-9064-3943df282fa4
# ╟─f40873af-3ef9-420b-b0a8-ed9f56b17047
# ╟─4536c6b7-f3bc-42e8-b044-0633f00b56bb
# ╠═bc55fb07-1d84-4927-80d2-c1a012409400
# ╟─c224f1b8-e5a1-4379-ae50-981ae89f9eaa
# ╠═4fe2b4cd-d0b0-489a-a867-14c58637327b
# ╟─6cf5e2f2-344c-4a7a-8901-5a0ee9cb8689
# ╠═551859cb-9cd1-4977-a528-b96423e83504
# ╠═955524a6-1a53-4d2d-ad16-a5a3a2a39b80
# ╠═118a79cc-afaf-4a99-8486-8e0e01fda40a
# ╟─8536cc9e-6cd4-4e38-9988-e45b06b75ea4
# ╟─25db8d7b-e467-4a2f-b3f3-0ad7df64988a
# ╠═01241971-cc9a-42d6-8e52-c840b91e6431
# ╟─e2db7492-4d68-4045-91b6-dca3aef8b514
# ╠═7a890da6-7547-4854-b877-a5d935a9f2dd
# ╟─c95bcde6-ffb9-4ce6-9dfb-38c79cc2e4b8
# ╠═bf90a66a-8505-4b33-b821-d708a3b7f8b0
# ╠═7af67f91-58e2-429d-a88f-b27a0258a805
# ╟─3ee0d4de-fdf5-40f0-a849-63b0759eb44b
# ╠═2c048fc0-5908-4f01-a9a5-480e12e098cc
# ╟─19ed4aa9-776e-4a80-b180-65d32f3f9f26
# ╠═12a11f90-fa55-47f2-a75f-853cd54cbeab
# ╟─f6b44166-3ba9-4770-80db-303bd28112cf
# ╟─02037c91-448e-4714-b932-86f4df631907
# ╟─cfea46a6-636c-4fcd-a04b-714e7a567cb3
# ╠═cb9384a3-96b0-4694-a218-25ebf96187a9
# ╟─73917400-43e0-4b9f-9448-9bbe02daebb9
# ╠═afba526d-0f71-4ba5-9a84-3c24804f959a
# ╟─903380ce-1320-464e-a703-1222b4f381b1
# ╠═22bd4450-6c6d-4a4d-aa1d-a080d3b06eac
# ╟─34f8cb7d-7482-4883-87e0-5c68be7ec0e9
# ╟─4d976b86-74d6-4bb5-af17-8a6bd658abe2
# ╠═ba7974ec-3965-4c92-ab37-08cdf5ab8125
# ╟─5e6837d5-35e9-434c-9337-7c90daec9c33
# ╠═b97ca16c-e6c3-4840-a855-1d29140eff59
# ╠═1af53492-f363-4d39-bfaf-2f668a441725
# ╟─e176e6ce-33da-4f3b-860e-4db99a931079
# ╟─770e70d0-8a69-4585-a423-41335d453dcd
# ╟─f4f3b1fb-165d-4f77-88f1-8c1a79f3fb42
# ╟─4a86d852-634c-43f5-90e5-3e675aeb2d69
# ╟─8a43d4b8-1902-4e02-9cd7-c36131ec9528
# ╠═d149d829-7d3d-4f7a-bdad-b1debb2c5149
# ╠═e2bec1ab-f161-45b6-a67e-1f4f197ae685
# ╟─5db6b1a5-ed2e-4812-94db-37cfffb46227
# ╟─bc8d315a-ddfb-4c82-84d0-90a4fb05256c
# ╟─51b348ef-0270-41d0-bbcb-a27eb10ff5d2
# ╠═7bf8fa3a-ce98-4f40-8d6c-4e84b79ae1c0
# ╠═498e0b9b-17a9-4059-9196-c6c4d86f926e
# ╟─b57131ad-e9d6-4318-bca7-c59c5d8da4c3
# ╠═f2fe206b-2832-42f6-84bd-20c3838fcb0b
# ╟─d51fff9b-17d3-4a0d-aa5e-01d58821a246
# ╠═6d95d3b1-d964-4b61-adb6-f44875e23ecb
# ╠═b1783355-0221-4f2d-a050-932b40920eab
# ╠═3c69edb3-ba7d-4ee0-996a-40fd985e57d3
# ╟─982e1de1-0820-42fb-8abb-31f1c99def46
# ╟─79f0f1dd-850f-4dfd-b895-ff41a2d8adb8
# ╟─fa2270d0-2548-409c-a23f-4369d8bce8ec
# ╟─545a317d-8bca-4b6c-9915-84b0b4ffbed5
# ╟─794b93c9-0a0a-4e77-b13c-0e2c06a9a0ec
# ╠═cfca0095-81e7-4c88-84bd-773eb5af7aa4
# ╟─ba9162dd-27f4-42b3-a6b5-ca9531614d5e
# ╟─b1f03df5-4353-45ed-9d91-9e2c9960f706
# ╠═89bf91c6-117c-4647-bfb9-6fc9b8dcfb5f
# ╟─b9a8c0cb-fe38-448e-aa4c-e9422f24a4a4
# ╠═6b5fc852-1d40-45ce-a2e2-d22a6cc03156
# ╟─29591722-bb44-47d0-a336-b2151502fd96
# ╠═ec032fb5-76a7-4a0c-a8be-b0f6bd166580
# ╟─d44ed37a-d022-457c-befe-01bae5cd0ed9
# ╠═5d3f6390-1f71-4543-ae4e-d72f28b94f40
# ╟─9eb6e432-3f41-4b5f-9194-241e6ab3a039
# ╠═cc0fc9ba-5564-4104-b0b3-a10ae39ef933
# ╠═6083305f-c838-4770-95c5-258fd7becbcb
# ╠═9777f3e8-785d-4ff5-a6c6-80a68a9c53ba
# ╟─e80de22e-d5e6-4af7-af0e-feaf9ef660bf
# ╟─0020352d-52dc-4ca4-90dd-f53865306f9e
# ╠═e58c49dd-5365-4dcb-9baa-c7a6048eb0ed
# ╟─03a35cce-cacf-4b11-a0b8-3a26ad93c227
# ╠═b3556b32-4665-439e-8751-a3ffa8ba1907
# ╟─ae030bbd-6416-44b8-ba9e-3fe9be5c9c52
# ╠═bfd9d570-7b12-4bbf-aaa9-449f2377a294
# ╠═f9e1ca7a-2a27-46aa-82a2-5e06882a9fad
# ╟─0e62cf17-f885-4c38-8df2-39d1719bfa04
# ╠═fe3c059f-cb4f-4ddf-ad54-c08aa70bc740
# ╟─fcfd18a2-154f-46fd-b4a5-01e406b500f1
# ╠═78d5ad0e-318f-4af1-883b-b887f8d72591
# ╟─58006718-5e76-4c6b-a968-aa5372bb4a15
# ╟─7606cd16-7463-462e-83e6-5adaea4632a2
# ╠═2379169c-7f2c-415c-8955-e8b06fc9edf5
# ╠═577bb77e-d8c6-4121-84e6-1d3c000ce7b9
# ╟─ac2d3b13-54d2-4217-8b18-7f9e276f1b6f
# ╠═2174d25b-b603-4d7b-b396-7ea00780143c
# ╠═0e8a6614-3164-4eab-b28b-fde3a328cbc3
# ╠═f2805dab-cac4-4a31-8418-6763c14d6415
# ╠═41d8cfe2-a997-4731-971e-3e0c6166177e
# ╠═2b2110eb-b867-47a1-8238-3d43a64851bb
# ╟─4fb53b80-4ea1-47a0-871a-adde769f6d18
# ╠═22b58c97-1aa2-49c4-8a0a-488c63014a90
# ╟─41abf8c4-e67e-4f66-a38a-7204a878d98d
# ╟─b02787ec-40eb-4000-9684-e69cbd912c6c
# ╟─e34fed6f-6b58-4642-84eb-167467881bb2
# ╠═68895842-0a1b-4236-b159-19a2805ea44d
# ╟─ec12d3c4-6b48-4f3a-b463-0fc67d3adff3
# ╠═1f8b10ba-f380-492e-800c-2673774cb25c
# ╟─82bee644-7b9a-42fa-a32f-76b3eb5038b5
# ╠═e481b87c-c6f3-48c1-bfff-5d76de0dce36
# ╟─a438470b-e037-4b65-b6b0-fb9c202e74a4
# ╠═f4c3dcd0-70ec-49ab-bbad-04b550da481c
# ╟─df166b45-3c28-4bad-95b8-f43429a58fd0
# ╠═a2738236-86a1-419a-9617-b1229f7c9240
# ╟─59f78a22-cac1-49a7-b3e6-9d74b694be64
# ╠═7c351d72-c5ce-4e5a-b2e5-87ff9dbc70e5
# ╠═d7312308-7c95-4c14-9bd6-56a7ac99d49a
# ╠═8f8b23b3-2e33-46e4-917a-294934fa090e
# ╟─4ff1079f-8dfb-4069-9a98-0d263eba9920
# ╟─9edd77c3-eb1e-4be7-bb6c-83803e675102
# ╠═c02a4f18-d51d-4954-b7a4-a4ae64b10fb7
# ╟─060b5869-881e-4d57-8357-cff95912ba16
# ╠═0d225f18-bdf6-49d1-a827-2bc9551c158d
# ╟─5d27e199-ad09-454a-8966-1531371ce67b
# ╠═79c4439b-df8f-467a-b791-5581bd564996
# ╠═bcbe7409-e969-4c7a-9417-195d406617bf
# ╟─6dbc81c8-4bc7-4499-a9d1-9719df368176
# ╠═ef5270de-da37-4a01-8734-dac1fd190a05
# ╠═b9030c3f-2cbb-44bf-90b7-23d116de4001
# ╠═2e684fd0-8410-4037-a6a6-2c472d6592ab
# ╠═79797d6d-9ca8-4f4f-9694-8045a8cf8cbb
# ╠═2cd566bb-3088-43f8-8926-9952342b5e5c
# ╠═750edde4-16b5-4025-82e0-078c0eacaca9
# ╠═1bbbd1ae-3540-44da-bc53-c8182a85b611
# ╟─3afb2a69-32d6-42c2-8c1e-51e64b83a6b4
# ╠═01ddf3dd-5199-48cf-ab13-75686278f3dc
# ╟─36f71777-ad33-4b27-9a0e-ff03b22ad79e
# ╠═d3b23128-2e7a-4246-a1f2-7325e0c2c432
# ╟─a9e5b236-67d5-400b-ba9a-8b565f6feeb1
# ╟─3bf0c957-f08b-4a45-90e0-40d56330c8e0
# ╠═947e26fe-2010-48ad-b927-d4e5186f21e0
# ╠═3506765b-a6cf-41ca-9201-fcc650c6f56d
# ╠═72dc27a0-0e61-4073-843a-67e1bfd1427a
# ╟─6240553c-0c7c-45f4-b190-ebacf09ab632
# ╠═12f64333-36d0-4e2a-9061-e3dcc7a4ae96
# ╟─95ba48fe-f430-4f47-ab04-f5896d75343b
# ╟─ef293617-7a7a-4597-99d2-2c39c543c296
# ╠═8bbed35b-fec5-4682-8b45-99fd3bd55bf2
# ╠═414ebb03-4fb7-4591-90e9-9ce3cf11e93d
# ╠═7a633bd6-10f1-4755-94ac-40d1dc9580a7
# ╠═aef172f1-8563-4e51-8a88-ec6b764a735e
# ╠═091e978a-fd89-455a-9b97-f3088d580683
# ╠═d3d6e099-66b0-444e-8142-ccd1718ac9ab
# ╟─2de0f619-63e0-444b-916b-316e4ebfa13e
# ╠═ef5329d2-ac63-4b49-a4cb-e8d6281a3184
# ╟─0b8e6775-7562-4c93-b64d-cb694531f250
# ╟─c99ac3a9-1d5e-43a3-a4ef-a22147187bb3
# ╟─b620006e-d6ac-4038-960f-8d0d93c2de00
# ╠═00dbb9d0-2c9d-4559-91d8-1a54171cc257
# ╠═3d5e9bb7-5e65-4de5-a521-6b1f1cb872b3
# ╟─3a3cd103-78b1-481a-8746-749d402346d1
# ╟─b9986617-923e-4744-9d44-1935770edbd7
# ╟─6e0c869d-6f90-4104-8589-41c2fa1b7342
# ╟─6ca83977-cb8c-4465-ab2e-c91971f16157
# ╠═fb76009b-3aaa-4eb5-8d69-45372a23af00
# ╠═b9c59973-512b-48eb-a2b1-831449e120a4
# ╠═cfeb28b8-2dd7-4c1d-8069-e8fe00d90567
# ╟─56784bdf-2a78-4986-9dd2-790a2a2ceea9
# ╠═bfb99943-1637-4dd6-adee-8a9e79274672
# ╟─a08fb39b-3911-450e-8377-cc94398825a9
# ╠═2575199f-a3ec-470f-87ea-d62be28c42dc
# ╠═6391d316-8826-4f19-9a2e-7065db787bee
# ╟─2879ff01-860f-4fb1-b52d-3a43e57ccb23
# ╠═3b78c1ff-a2f6-4c1c-8124-5fe58c5a312d
# ╟─75e221b4-fdfe-4ef0-a8ea-1d75700fc1e7
# ╠═0f228304-de82-48f2-8e43-6ef81f7a9830
# ╠═9843e282-44b3-40d0-b2b5-506144c22b27
# ╠═eb9fe668-0a0f-4a8f-a793-ba545cb02e85
# ╠═ddf9ae6d-0abe-4c8c-9646-024dd9a41e02
# ╠═8330fce5-13c9-447e-8368-01b7c7ebc68f
# ╠═0d2cab0c-04a2-4c8b-bb4f-2115ae0cf871
# ╠═6900e703-de7d-4ea3-8ac5-0b38b21b083d
# ╠═b353eab2-d78f-4c5e-8e38-bd57813b6816
# ╠═48dd67e7-363a-4e7c-ad81-8f70fc6db847
# ╠═a7fa7870-f202-4694-a9c1-7021ed242966
# ╟─1699d3ae-d77a-4272-8c4b-ddd4eadd503f
# ╟─72878862-dcc0-4c79-9c91-7c0bc44359d3
