### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ f8a92690-990b-4341-89e1-322adbcb8d1b
using Pkg; Pkg.activate("..")

# ╔═╡ a09f814a-0c6a-11ef-0e79-a50b01287d63
using Markdown, InteractiveUtils

# ╔═╡ 015050b3-3339-4b1a-ad7d-c358cce73675
using Catalyst, ModelingToolkit, OrdinaryDiffEq

# ╔═╡ a83f424f-1bbb-4e71-af81-302c0ce68907
using ModelingToolkit: t_nounits as t, D_nounits as D

# ╔═╡ dbfe4800-0974-4ca1-bb0a-d8803409a98b
using Turing, StatsPlots, StatsBase

# ╔═╡ 6e227e07-166a-41ce-839a-4c4c72addb23
using LinearAlgebra, Optim

# ╔═╡ b992c080-a0ce-4188-b632-e734a141e67d
using PlutoUI; TableOfContents()

# ╔═╡ 37da8786-fea0-4c2f-a76f-6e6c68325a78
md"""
# Introduction to calibration
"""

# ╔═╡ 4623369d-8c5a-422d-9e40-0f1dd7586260
md"""
## Goal of this practicum
"""

# ╔═╡ 3cb0a166-ac53-4c3f-9832-e93742040cfb
md"""
In the models discussed in the previous sessions, we always knew the values of all parameters. In reality, the value of a parameter has to be calibrated, hence, estimated from experimental data. During this parameter estimation one attempts to find the set of parameter values for which the model predictions are as close as possible to the collected experimental data.
"""

# ╔═╡ 987f0a4d-e416-4ceb-adbe-3dcdca9d0996
md"""
The search of optimal parameter values usually involves a function, such as a loss function, a (log) likelihood function or a posterior distribution function. In this session we will be (mainly) using the MLE (Maximum Likelihood Estimation) and MAP (Maximum A Posteriori estimation) methods, which respectively try to maximize the likelihood function and the posterior probability of the data.

In order to understand the difference, recall Bayes' theorem applied to a set of parameters $θ$ and data $D$:
```math
P(θ \mid D) = \frac{P(D \mid θ) \, P(θ)}{P(D)}
```
In the MLE method, the likelihood function $P(D \mid θ)$ (= the probability of the data given the parameters) is maximized during the search of optimal parameter values of a model in order to fit experimental data. The parameter values are considered unknown but viewed as fixed points.

In the MAP method, the posterior probability $P(θ \mid D)$ (= the probability of the parameters given the data) is maximized instead. Instead of viewing the parameter values as fixed points, they are now treated as random variables in the model which follow a prior distribution. In other words, we have prior belief in which distribution these parameters come from (Normal, Beta, etc). Once new data comes in, we update our prior belief, leading to a posterior belief. Hence, we now have a better idea from which distribution these parameters come.

One caveat with the MLE and MAP methods are that they only return a point estimate of the optimal value, which gives no information on how certain we are about this value. A more informative, yet computationally exhaustive method is using the MCMC (Markov chain Monte Carlo) sampling methods, which approximate the entire posterior distribution of the unknown values.
"""

# ╔═╡ 75efff36-8da7-4d04-afa2-a2f8324bc103
md"""
In this notebook we will calibrate the different parameters involved in the grass growth models. To illustrate this concept, we first revisit the three simple models modelling the grass growth yield.
"""

# ╔═╡ 3dcb9c9d-370b-4031-b7c0-cee80742557a
md"""
## Grass growth models
"""

# ╔═╡ 7a14aa59-6e6f-4266-a0b3-84ab55f2efc5
md"""
In this notebook, three different models will be used, each modelling the yield of grass in a grassland:

- Logistic growth model: $\cfrac{dW}{dt} = \mu \left( 1 - \cfrac{W}{W_f} \right) W$
- Exponential growth model: $\cfrac{dW}{dt} = \mu \left( W_f - W \right)$
- Gompertz growth model: $\cfrac{dW}{dt} = \left( \mu - d \ln(W) \right) W$

with output $W$ the grass yield, and $W_f$, $\mu$ and $d$ parameters. The table below shows some typical parameter values and initial conditions for grasslands similar to the one we observed, which we can use as prior information.

|             | $\mu$      | $W_f$       | $d$          | $W_0$          |
|:----------- |:----------:|:-----------:|:------------:|:------------:|
| Logistic    |  0.07      | 10.0        |              | 2.0          |
| Exponential |  0.02      | 10.0        |              | 2.0          |
| Gompertz    |  0.09      |             | 0.040        | 2.0          |

Hence, for each grass growth model, we will optimize the parameter values together with the initial value.
"""

# ╔═╡ 85cd60a8-b448-4375-9b6d-399c4336c319
md"""
In each of the three models we will use the following timespan:
"""

# ╔═╡ 5b320989-3e0b-447b-bc9a-25fb221ce609
tspan = (0.0, 100.0)   # this will be the same for the three models

# ╔═╡ 2481cd4f-0efc-4450-ab3d-4a5492597f36
md"""
Variables containing the initial condition and parameters values will be defined later in the objective function.
"""

# ╔═╡ 1aa44f2b-6f33-437f-b9dd-89762d9f28ea
md"""
### The measurement data
"""

# ╔═╡ b2b433ed-0266-4bea-a7e8-32adba542d4c
md"""
Assume that the measured grass yields (of a certain plant type) are the following:
"""

# ╔═╡ 7c966a66-0091-4b81-9a7e-02ccd0d3db10
W_meas = [1.87, 2.45, 3.72, 4.32, 5.28, 7.01, 6.83, 8.62, 9.45, 10.31, 10.56, 11.72, 11.05, 11.53, 11.39, 11.7, 11.15, 11.49, 12.04, 11.95, 11.68]

# ╔═╡ 3edd2acc-a865-4675-afef-8868c68256f1
md"""
They have been measured at the following corresponding time instances:
"""

# ╔═╡ 877298e8-b61b-4c3a-ba2c-2827acdcfb50
t_meas = 0:5:100

# ╔═╡ ef06cc43-510b-4ff9-b0b7-1c7fc267e9b1
md"""
We can make a scatter plot of this data (including a title, a legend label, an X-axis label, X- and Y-axis limits) in the following way:
"""

# ╔═╡ cb2bc6ee-4211-47e1-9956-5cf1b0c0671d
scatter(t_meas, W_meas, title="Grass growth data",
                        label="Yield",
                        xlabel="t",
                        xlims=(0, 100),
                        ylims=(0, 14))

# ╔═╡ 50e6af5b-04f4-495c-853b-c746fd27254f
md"## Calibration with Turing"

# ╔═╡ a665740a-ea94-4452-8cbc-a1647319dfce
md"""
We will be using the familiar Turing framework to perform optimisation, which consists roughly of the following steps:
- Define your model as a Turing model, with the following components:
  - The observed input data as input of the model
  - Priors for all unknown values such as initial values, parameter values and the measurement error
  - An (ODE-based) model that predicts values of the output variable based on the observed input values and the unknown initial values, parameter values, etc.
  - The relationship between the observed outputs and the predicted outputs, which often comes down to specifying the measurement error
- Instantiate the Turing model and condition it on the observed values of the output variable
- Call an optimisation algorithm on the Turing model
- Extract the optimized parameters
- Visualise the results
"""

# ╔═╡ b2f20a9a-0613-4bb4-a1d3-df1620aef105
md"### The priors"

# ╔═╡ 6852f166-350e-4434-8304-382c01113ba1
md"""
During calibration, the prior distributions have two roles:
1. Define for every unknown quantity **the bounds**: the optimizer will only search within the domain of the prior distribution.
1. Define part of the **posterior probability** of observing the data given the parameter values. In other words, the theoretically optimal values depend on your prior distributions.

**Note that the second only applies to methods that take prior information into account, which Maximum Likelihood Estimation (MLE) does not.**
"""

# ╔═╡ 4f78b2e8-9132-4175-b9ae-9f5979649e90
md"#### Choice of priors for the grass growth models"

# ╔═╡ 557e69c5-d0c6-41ee-b3df-851bd1d16923
md"""
In general for our grass models, we will need to define priors for the following:
- the measurement error (standard deviation) $\sigma_W$.
- the initial condition $W_0$.
- the parameters $\mu$ and either $W_f$ or $d$ (depending on the model).
"""

# ╔═╡ 237934b7-c89a-4ef6-938f-4629c0fa4c95
md"""
For the measurement error $\sigma_W$, a distribution with most of its probability density around 0 and a long positive tail is often a good choice. We will use the Exponential distribution here for this purpose. As for its parameters, considering our yield is in the order of 1 to 10 [t/ha], an expected value (and therefore standard deviation) of around 1 seems appropriate.
"""

# ╔═╡ 513ee94d-5d49-43d2-a1db-638fe9e2de1a
@bind example_mean_error Slider(0.01:0.01:10.0, show_value = true, default = 1)

# ╔═╡ 8e163f77-b762-427a-b005-313284dc2e9e
plot(Exponential(example_mean_error), title = "Example prior for the measurement error", legend = false)

# ╔═╡ d9c38e9d-12b8-4e98-8546-bbc123960e6f
md"""
For the initial- and parameter values, we will default to log-normal distributions. These have the following properties that generally fit well with biological parameters, such as is the case in this notebook:
- A positive domain, which bounds possible values to the positive numbers
- Most of the probability density around the expected value
- A long postive tail, which allows for outliers

As the distribution of a variable whose logarithm is normally distributed, it does have some strange behaviour: the arguments of the distribution are the mean and standard deviation of the **logarithm**. To clarify, for $X \sim \mathrm{LogNormal(μ, σ)}$:
-  $μ = E[\mathrm{log}(X)]$
-  $σ = \sqrt{\mathrm{Var}(\mathrm{log}(X))}$

Practically, this means **you need to specify the logarithm of the desired expected value and play around with the standard deviation**. For this exercise, we will choose an expected value based on the table discussed earlier in this section and keep the standard deviation to the default value of 1.
"""

# ╔═╡ 39c70e22-86e0-467e-abab-24f99548b096
@bind example_log_μ Slider(0.01:0.01:10.0, default = 0.1, show_value = true)

# ╔═╡ ddef987d-59db-4adf-ae77-d4d32a41ddec
@bind example_log_σ Slider(0.01:0.01:3.0, default = 1.0, show_value = true)

# ╔═╡ 54fb706b-9036-42b5-833e-4081f4c2179c
begin
	example_param_prior = LogNormal(log(example_log_μ), example_log_σ)
		
	println("""Look at this strange distribution with 
			- mean = $(mean(example_param_prior))
			- median = $(median(example_param_prior))
			- σ = $(std(example_param_prior))""")
	plot(example_param_prior, xlims = (0, 5*example_log_μ), title = "Example prior for the model parameters", legend = false)
end

# ╔═╡ 9a5bc72b-346d-4e95-a873-783037ed98bc
md"""
## Example: the Logistic growth model
"""

# ╔═╡ 4b17f7e6-26b0-425a-9f53-ccf1689639fe
md"""
We will illustrate the calibration with the logistic growth model. The latter can be done via a Catalyst reaction network model or via a model built with ModelingToolkit. As an example, we will show both possibilities for the logistic growth model.
"""

# ╔═╡ 6736542c-5378-480e-a56b-65956b416225
md"""
#### 1) Catalyst based model
"""

# ╔═╡ ba56adb1-9405-40d5-be48-4273b42ab145
# growth_log = @reaction_network begin
# 	@species W(t)=2.0
# 	@parameters μ=0.07 Wf=10.0
# 	μ*(1-W/Wf), W --> 2W
# end

# ╔═╡ 1c311b6f-6b18-4170-b500-33a8e4d3cb2d
md"""
#### ̇2) ModelingToolkit based model
"""

# ╔═╡ 50220e1f-8e42-44da-be03-08c11df967f0
@variables W(t)

# ╔═╡ 449286b0-0210-4c6b-b28f-fca87b52d674
@parameters μ Wf d

# ╔═╡ f720eac9-cb29-4eef-ac94-471395941c0f
# IF YOU UNCOMMENT THE FOLLOWING LINE, PLEASE COMMENT THE CATALYST MODEL FIRST!!!
@mtkbuild growth_log = ODESystem([D(W) ~ μ*(1 - W/Wf)*W], t)

# ╔═╡ d75246d4-e03b-4684-be7d-4bcfb61ed7ef
md"""
### Declaration of the Turing model
"""

# ╔═╡ 8a9115eb-4044-4cab-a7db-39b5dd86c70d
@model function growth_log_inference(t_meas)
    σ_W ~ Exponential()
    W0 ~ LogNormal(log(1))
    μ ~ LogNormal(log(0.1))
    Wf ~ LogNormal(log(10))
	u0_log = [:W => W0]
	parms_log = [:μ => μ, :Wf => Wf]
	oprob_log = ODEProblem(growth_log, u0_log, tspan, parms_log)
    osol_log = solve(oprob_log, AutoTsit5(Rosenbrock23()), saveat=t_meas)
    W_s ~ MvNormal(osol_log[:W], σ_W)
end

# ╔═╡ 48c9f616-d298-40da-b917-225abd39b3d9
md"""
Some remarks:
- The time points are the ones from the measurements, therefore, we set: `saveat=t_meas`.
- Depending on the priors, the parameter values of our ODE may vary wildly during calibration. As this can influence the stiffness of the system, it can be beneficial to use an ODE solver that automatically detects the stiffness of the system and switches to a stiff solver if necessary, such as `AutoTsit5(Rosenbrock23())`. We don't expect you to know when this is necessary, but if your calibration regularly gives instability errors, this may help! 
"""

# ╔═╡ 35f158c1-858d-4e4d-ac3d-bf4807dad9a0
md"""
We will provide the measurements to the Turing model:
"""

# ╔═╡ 8b6534d6-776b-4285-8498-a9b34051facc
growth_log_inf = growth_log_inference(t_meas) | (W_s = W_meas, )

# ╔═╡ a6972aef-63ad-401c-acf5-6d59f9fc6698
md"""
We are now ready to optimize the priors ($\sigma_W$, $W_0$, $\mu$ and $W_f$). This is done by calling the `optimize` function, providing the previously created object `growth_log_inf`, the method for estimating the parameters and (optionally) an algorithm (default: Nelder-Mead) to implement the method.
"""

# ╔═╡ 73e35289-6dc0-4e2e-83eb-b56f83cdbbbf
md"""
### Method - Maximum Likelihood Estimation
"""

# ╔═╡ 47bd729c-4851-42f7-a03f-6ceacd3c717e
md"""
We will use the MLE (Maximum Likelihood Estimation) method here and store the optimization results in `results_log_mle`.
"""

# ╔═╡ f34bb7ac-1ed8-4dd9-b0b9-49bd6e0e1d71
results_log_mle = optimize(growth_log_inf, MLE(), NelderMead())

# ╔═╡ e55404ab-6762-4f39-bb42-9c195334a214
md"""
You can visualize a summary of the optimized parameters by piping them to `coeftable`. Beware that this can take a lot of time...
"""

# ╔═╡ 80e7f6b8-7592-48a3-8587-f1953d1bfcd8
# results_log_mle |> coeftable

# ╔═╡ a1ca7d0e-639c-42d4-be09-5c61a2008f29
md"""
You can obtain the actual optimized values using the function `coef` on the results object in conjunction by calling the parameters by name preceded by a colon. Here we assign the optimized parameter values to some suitable variable names:
"""

# ╔═╡ 30399b9a-1d77-4140-9ad3-5eed636a5b99
W0_opt1_log = coef(results_log_mle)[:W0]

# ╔═╡ e1b8e4ba-c1f1-48d5-87a2-edce19c9fe7a
μ_opt1_log = coef(results_log_mle)[:μ]

# ╔═╡ a943e0fe-1376-4ef1-9c45-25c6f95e3b96
Wf_opt1_log = coef(results_log_mle)[:Wf]

# ╔═╡ 72e065d4-7b1b-4f46-b373-935be8d801fc
md"""
Now we can make a plot of $W$ simulated with the optimized initial condition and parameter values.
"""

# ╔═╡ 0b2ffd6f-01cd-4f11-9062-d38b3c13a5b1
md"""
Setting up initial condition with optimized initial condition:
"""

# ╔═╡ 590b1006-0e37-4668-9b4a-3588fab45696
u0_opt1_log = [:W=>W0_opt1_log]

# ╔═╡ 8ff28a2e-185d-4dca-ad3b-a0b1507646d6
md"""
Setting up parameter values with optimized parameter values:
"""

# ╔═╡ 14d24cbd-3259-41bb-9013-b4fe25a3be4c
parms_opt1_log = [:μ=>μ_opt1_log, :Wf=>Wf_opt1_log]

# ╔═╡ 6c134677-3ec1-4e0a-88b5-01341a096675
md"""
Next, we create an ODEProblem and solve it:
"""

# ╔═╡ 7a81f4a0-8f7c-4e05-9c3f-2438eab9b691
oprob_opt1_log = ODEProblem(growth_log, u0_opt1_log, tspan, parms_opt1_log)

# ╔═╡ ac80099b-d8f0-4eba-809d-d482bd354d35
osol_opt1_log = solve(oprob_opt1_log, Tsit5(), saveat=0.5)

# ╔═╡ b58f2c24-e0ea-48a8-b0b7-d0faf9642340
md"""
Finally, we plot $W$ simulated with the optimized initial value and parameter values together with the measured data that was used to find the optimized values.
"""

# ╔═╡ 55eba435-6ca4-4f4f-b08a-be700d5bda91
begin
	plot(osol_opt1_log, label="Logistic growth", xlabel="t",
		xlims=(0, 100), ylims=(0, 14))
	scatter!(t_meas, W_meas, label="Yield")
end

# ╔═╡ 5d386b00-93b5-4a88-b4e3-e5c3eebd6dd5
md"""
### Method - Maximum A Posterior
"""

# ╔═╡ 6f0e91d0-6b99-4cf1-8145-523589a21e89
md"""
We will use the MAP (Maximum A Posterior) method here and store the optimization results in `results_log_map`.
"""

# ╔═╡ 209742ca-36bb-42a5-bf8a-291a40c47757
results_log_map = optimize(growth_log_inf, MAP(), NelderMead())

# ╔═╡ 07efb21a-0dac-46d5-b9e1-8f4757c6cedf
md"""
You can visualize a summary of the optimized parameters by piping them to `coeftable`. Beware that this can take a lot of time...
"""

# ╔═╡ a93a7eeb-aad1-49f3-a05f-1512adb08b6b
# results_log_map |> coeftable

# ╔═╡ ccc66805-efce-4f25-b76b-7cff23901ba4
md"""
You can obtain the actual optimized values using the function `coef` on the results object in conjunction by calling the parameters by name preceded by a colon. Here we assign the optimized parameter values to some suitable variable names:
"""

# ╔═╡ 57ef7d03-a6dc-40da-858d-5f9b9be613cd
W0_opt2_log = coef(results_log_map)[:W0]

# ╔═╡ 62bd2e58-343a-4155-a65e-cf326d0975f4
μ_opt2_log = coef(results_log_map)[:μ]

# ╔═╡ 4d6cba88-be68-412d-aeec-28300c26a9da
Wf_opt2_log = coef(results_log_map)[:Wf]

# ╔═╡ 4bc97af7-debc-4170-a2d5-4e6c928e5183
md"""
Now we can make a plot of $W$ simulated with the optimized initial condition and parameter values.
"""

# ╔═╡ a0033986-4946-45a2-a53e-29c38f35ed5c
md"""
Setting up initial condition with optimized initial condition:
"""

# ╔═╡ a570ebab-64de-4934-b887-77a8e2fb42e8
u0_opt2_log = [:W=>W0_opt2_log]

# ╔═╡ feee9cc1-ebe2-4c77-a11e-53ca00864ccf
md"""
Setting up parameter values with optimized parameter values:
"""

# ╔═╡ 73e86176-f314-463c-8821-51c0f3cc1a56
parms_opt2_log = [:μ=>μ_opt2_log, :Wf=>Wf_opt2_log]

# ╔═╡ c17c4e00-8688-4f9c-b0ae-741d70601aba
md"""
Next, we create an ODEProblem and solve it:
"""

# ╔═╡ f751b962-1810-429c-b34b-658b63fa9ba0
oprob_opt2_log = ODEProblem(growth_log, u0_opt2_log, tspan, parms_opt2_log)

# ╔═╡ 0fd058cc-9bd8-4e12-b5dd-79818519165b
osol_opt2_log = solve(oprob_opt2_log, Tsit5(), saveat=0.5)

# ╔═╡ 57c95382-69cc-47e7-aa94-4767136bb23b
md"""
Finally, we plot $W$ simulated with the optimized initial value and parameter values together with the measured data that was used to find the optimized values.
"""

# ╔═╡ 4e07c39b-502d-4b38-ae7f-f103cb4bae16
begin
	plot(osol_opt2_log, label="Logistic growth", xlabel="t",
		xlims=(0, 100), ylims=(0, 14))
	scatter!(t_meas, W_meas, label="Yield")
end

# ╔═╡ 29170e2a-9916-438e-92ca-9f4783397b5e
md"""
### Method - MCMC with NUTS
"""

# ╔═╡ 7ff9fe52-156b-4a92-9058-781670de3abb
md"""
We will use Markov chain Monte Carlo (MCMC) method in combination with the No U-Turn Sampler (NUTS) here and store the optimization results in `results_log_nuts`.
"""

# ╔═╡ 0c047043-3284-422a-9c88-2f4f4c170edf
results_log_nuts = sample(growth_log_inf, NUTS(), 500)

# ╔═╡ 78a608a4-0ba7-400e-abb9-a9285c60681c
plot(results_log_nuts) # check convergence

# ╔═╡ 19c362cb-2764-41c9-a571-2e8e2bfcde93
# summarize(results_log_nuts)

# ╔═╡ 93db47b2-34e8-43b4-beac-b5620fd444e7
W0_opt3_log = mean(results_log_nuts[:W0])

# ╔═╡ 68c71cd2-6cdc-4b80-b6e7-75bfea344295
μ_opt3_log = mean(results_log_nuts[:μ])

# ╔═╡ d6b9eaba-1d43-4e56-8fa1-bb2f87f6fe79
Wf_opt3_log = mean(results_log_nuts[:Wf])

# ╔═╡ b843ee5e-9618-4b50-93db-77e19b4be366
md"""
Now we can make a plot of $W$ simulated with the optimized initial condition and parameter values.
"""

# ╔═╡ 58014411-128c-41f4-b192-b815a3a8cc60
md"""
Setting up initial condition with optimized initial condition:
"""

# ╔═╡ ecf4b951-9b5f-441a-89f5-b0ccd040ba02
u0_opt3_log = [:W=>W0_opt3_log]

# ╔═╡ b92bf532-5610-4ce8-b469-154b20cc5956
md"""
Setting up parameter values with optimized parameter values:
"""

# ╔═╡ e629fbd1-a7c7-4dc3-8b7d-bdca380fe5ad
parms_opt3_log = [:μ=>μ_opt3_log, :Wf=>Wf_opt3_log]

# ╔═╡ a48df980-4847-418e-813d-1875d032ffd9
md"""
Next, we create an ODEProblem and solve it:
"""

# ╔═╡ 348e7454-c6dc-4d1e-b18b-6b34bc2fc0fc
oprob_opt3_log = ODEProblem(growth_log, u0_opt3_log, tspan, parms_opt3_log)

# ╔═╡ f71f0cf4-b63e-452c-b971-f30d7503f1e5
osol_opt3_log = solve(oprob_opt3_log, Tsit5(), saveat=0.5)

# ╔═╡ 7346a003-aa05-41e3-a420-c97c1994cf3c
md"""
Finally, we plot $W$ simulated with the optimized initial value and parameter values together with the measured data that was used to find the optimized values.
"""

# ╔═╡ 20a6f165-449a-42cc-9565-f342a7535422
begin
	plot(osol_opt3_log, label="Logistic growth", xlabel="t",
		xlims=(0, 100), ylims=(0, 14))
	scatter!(t_meas, W_meas, label="Yield")
end

# ╔═╡ 137bde23-76f2-4ebf-8bc2-ea8640001436
md"""
## Exercises
"""

# ╔═╡ 4aa71200-006b-4a15-ae75-67e36aa81522
md"""
### Exercise 1 - Calibration of the exponential growth model

Calibrate the initial condition and both parameters of the exponential growth model. Use the values mentioned in the *Table* as initials values for the optimization of the parameters.
"""

# ╔═╡ cdab3079-04b0-4a44-b770-468c20e321e4
md"""
Implement a *reaction network object* using Catalyst for the exponential growth model:
"""

# ╔═╡ cf1a144e-09e9-42a3-b2a3-b8676a200a39
# growth_exp = @reaction_network begin
#     missing
#     missing
# end
growth_exp = @reaction_network begin
    μ*Wf, 0 --> W
    μ, W --> 0
end

# ╔═╡ d69ed5dd-80cd-427e-9245-e42576a4688d
md"""
Convert the reaction network model into an ODE system to verify.
"""

# ╔═╡ 8b429656-1c57-4278-9633-71311f400ed8
# missing
convert(ODESystem, growth_exp)

# ╔═╡ b4300e8a-8052-419b-98c8-0508ebee2393
md"""
Declare the Turing model. Take the same priors (and distributions) as before.
"""

# ╔═╡ 2c6ae74c-2da4-4867-ad8a-f4e835101d63
# @model function growth_exp_inference(t_meas, W_meas)
#     σ_W ~ missing
#     W0 ~ missing
#     μ ~ missing
#     Wf ~ missing
#     u0_exp = missing
#     parms_exp = missing
#     oprob_exp = missing
#     osol_exp = missing
#     W_meas ~ missing
# end
@model function growth_exp_inference(t_meas)
    σ_W ~ Exponential()
    W0 ~ LogNormal(log(1))
    μ ~ LogNormal(log(0.1))
    Wf ~ LogNormal(log(10))
	u0_exp = [:W => W0]
	parms_exp = [:μ => μ, :Wf => Wf]
	oprob_exp = ODEProblem(growth_exp, u0_exp, tspan, parms_exp)
    osol_exp = solve(oprob_exp, AutoTsit5(Rosenbrock23()), saveat=t_meas)
    W_s ~ MvNormal(osol_exp[:W], σ_W)
end

# ╔═╡ c6a5d453-d610-4f65-847c-c878dd41726c
md"""
Provide the measurements to the Turing model.
"""

# ╔═╡ fff17cf7-173d-4f64-94a9-4bf46acc882d
# growth_exp_inf = missing
growth_exp_inf = growth_exp_inference(t_meas) | (W_s = W_meas, )

# ╔═╡ eee55784-a641-445e-be75-0b19e2a94754
md"""
Optimize the priors ($\sigma_W$, $W_0$, $\mu$ and $W_f$). Do this with `MLE` method and Nelder-Mead. Store the optimization results in `results_exp_mle`.
"""

# ╔═╡ 7844e4f5-3c7d-4b4b-beee-970c998c67a6
# results_exp_mle = missing
results_exp_mle = optimize(growth_exp_inf, MLE(), NelderMead())

# ╔═╡ c81d0140-3f4e-4eb4-8a77-1f48c5e0ecbf
md"""
Visualize a summary of the optimized parameters. Beware: this can take a lot of time...
"""

# ╔═╡ 7456455b-4f31-488f-990f-6ce534038e08
# missing
# results_exp_mle |> coeftable

# ╔═╡ b10c2ce4-d363-429c-a64c-ec29652137a5
md"""
Get the optimized values and assign them to `W0_opt_exp`, `μ_opt_exp` and `Wf_opt_exp`.
"""

# ╔═╡ 23b629a6-6866-416a-a768-9617ce6301db
# W₀_opt_exp = missing
W0_opt_exp = coef(results_exp_mle)[:W0]

# ╔═╡ 8788082d-f5d1-4385-8037-a0d360a841c7
# μ_opt_exp = missing
μ_opt_exp = coef(results_exp_mle)[:μ]

# ╔═╡ 03a4fa85-08db-46d4-bb53-c0ccea90a211
# Wf_opt_exp = missing
Wf_opt_exp = coef(results_exp_mle)[:Wf]

# ╔═╡ 8f5e1413-227c-44fa-bb2d-3653cbc27e38
md"""
Make a plot of $W$ simulated with the optimized initial condition and parameter values.
"""

# ╔═╡ 8a7f7aab-878e-41b5-b9da-d06747df042e
md"
Set up initial condition with optimized initial condition:
"

# ╔═╡ 30602ff1-041b-4fca-bf8e-55ff57df9e37
# u₀_opt_exp = missing
u0_opt_exp = [:W=>W0_opt_exp]

# ╔═╡ 881011be-6434-416f-915b-3333e8dea32f
md"""
Set up parameter values with optimized parameter values:
"""

# ╔═╡ 5602b88a-07f8-438b-994c-65f11e17a0ba
# parms_opt_exp = missing
parms_opt_exp = [:μ=>μ_opt_exp, :Wf=>Wf_opt_exp]

# ╔═╡ 25be4255-0888-4ecd-a2fd-d66402c5cb50
md"""
Create an ODEProblem and solve it. Use `Tsit5()` and `saveas=0.5`.
"""

# ╔═╡ 36e8a174-d526-45ee-b3c6-88d698ad5d5f
# oprob_opt_exp = missing
oprob_opt_exp = ODEProblem(growth_exp, u0_opt_exp, tspan, parms_opt_exp)

# ╔═╡ 7594147d-b3da-4e0d-896d-41baacb6d7be
# osol_opt_exp = missing
osol_opt_exp = solve(oprob_opt_exp, Tsit5(), saveat=0.5)

# ╔═╡ 8e047be9-f0f0-4a75-91b5-c523f55f8c67
md"""
Plot $W$ simulated with the optimized initial value and parameter values together with the measured data that was used to find the optimized values.
"""

# ╔═╡ 1a9587aa-2356-48df-abcc-2ce874fa5d24
# begin
# 	missing
# 	missing
# end
begin
	plot(osol_opt_exp, label="Exponential growth", xlabel="t",
		xlims=(0, 100),ylims=(0, 14))
	scatter!(t_meas, W_meas, label="Yield")
end

# ╔═╡ 785d500b-f8ea-446a-9952-2a5fd5d83d24
md"""
### Exercise 2 - Calibration of the Gompertz growth model

Calibrate the initial condition and both parameters of the Gompertz growth model. Use the values mentioned in the *Table* as initials values for the optimization of the parameters.
"""

# ╔═╡ e754826a-7411-4072-b0dc-a4bad7a15f98
md"""
Implement a system using MTK for the Gompertz growth model.\
*Hint*: no need to redefine the variable `W(t)` and the parameters `μ` and `d` because they are defined at the beginning of this notebook.
"""

# ╔═╡ 6e85c08e-9a97-4897-9c81-89517f55e254
@mtkbuild growth_gom = ODESystem([D(W) ~ (μ - d*log(W))*W], t)

# ╔═╡ e5081280-d226-4834-8932-c89becd8313c
md"""
Declare the Turing model. Take the same priors as before, and try to find a sensible prior distribution for $d$.
"""

# ╔═╡ c739a908-2353-4e7a-8fbd-f640dc8cabe0
# @model function growth_gom_inference(t_meas, W_meas)
#     σ_W ~ missing
#     W0 ~ missing
#     μ ~ missing
#     d ~ missing
#     u0_gom = missing
#     parms_gom = missing
#     oprob_gom = missing
#     osol_gom = missing
#     W_meas ~ missing
# end
@model function growth_gom_inference(t_meas)
    σ_W ~ Exponential()
    W0 ~ LogNormal(log(1.0))
    μ ~ LogNormal(log(0.1))
    d ~ LogNormal(log(0.1))
	u0_gom = [:W => W0]
	parms_gom = [:μ => μ, :d => d]
	oprob_gom = ODEProblem(growth_gom, u0_gom, tspan, parms_gom)
    osol_gom = solve(oprob_gom, AutoTsit5(Rosenbrock23()), saveat=t_meas)
    W_s ~ MvNormal(osol_gom[:W], σ_W)
end

# ╔═╡ 1d0383ad-54d6-4ff2-8555-def83bfff0e6
md"""
Provide the measurements to the Turing model.
"""

# ╔═╡ cd1cf2f8-9f7f-4ed4-9cb7-1a6efee68ab4
# growth_gom_inf = missing
growth_gom_inf = growth_gom_inference(t_meas) | (W_s = W_meas, )

# ╔═╡ aba74ee0-0163-4e15-8b49-d8dcad4839f7
md"""
Optimize the priors ($\sigma_W$, $W_0$, $\mu$ and $D$). Do this with `MLE` method and Nelder-Mead. Store the optimization results in `results_gom_mle`.
"""

# ╔═╡ 0eda4142-1aaf-4e17-bd78-857e13e94acd
# results_gom_mle = missing
results_gom_mle = optimize(growth_gom_inf, MLE(), NelderMead())

# ╔═╡ 50629194-98ed-4d45-86a2-95ac22daac29
md"""
Visualize a summary of the optimized parameters. Beware: this can take a lot of time...
"""

# ╔═╡ 9c239fc6-275c-4d64-9fa2-6fd57295b757
# missing
# results_gom_mle |> coeftable

# ╔═╡ dede17f2-655d-4871-b6de-5a32804947dd
md"""
Get the optimized values and assign them to `W0_opt_gom`, `μ_opt_gom` and `D_opt_gom`.
"""

# ╔═╡ e8c3f042-6058-4040-a67e-18563c04ee93
# W₀_opt_gom = missing
W0_opt_gom = coef(results_gom_mle)[:W0]

# ╔═╡ 665f4d03-3521-475a-a195-f861fd26bb69
# μ_opt_gom = missing
μ_opt_gom = coef(results_gom_mle)[:μ]

# ╔═╡ 13775d58-ac61-431c-a6c9-447c1eec7942
# D_opt_gom = missing
d_opt_gom = coef(results_gom_mle)[:d]

# ╔═╡ bc92f996-b626-4965-a286-d2c848eb1a21
md"""
Make a plot of $W$ simulated with the optimized initial condition and parameter values.
"""

# ╔═╡ e2fde8e9-1f87-4ffe-8dae-2794664bfaa4
md"""
Set up initial condition with optimized initial condition:
"""

# ╔═╡ 23fd9fd9-a13c-4c56-a0b4-daec6f1d2cd8
# u₀_opt_gom = missing
u0_opt_gom = [:W=>W0_opt_gom]

# ╔═╡ 6f414d15-af0a-452a-98a1-dc0b7e54d617
md"""
Set up parameter values with optimized parameter values:
"""

# ╔═╡ 57ee8a12-24df-4598-935c-f5e259b504cb
# parms_opt_gom = missing
parms_opt_gom = [:μ=>μ_opt_gom, :d=>d_opt_gom]

# ╔═╡ 48bc085c-9ce6-4752-a5a9-a814f803f571
md"""
Create an ODEProblem and solve it. Use `Tsit5()` and `saveas=0.5`.
"""

# ╔═╡ 5b9b2e5b-9deb-4ff6-a923-b15b8b08f0c9
# oprob_opt_gom = missing
oprob_opt_gom = ODEProblem(growth_gom, u0_opt_gom, tspan, parms_opt_gom)

# ╔═╡ 66ee9655-a006-4c0f-b1f1-5576efa8f896
# osol_opt_gom = missing
osol_opt_gom = solve(oprob_opt_gom, Tsit5(), saveat=0.5)

# ╔═╡ 52d7975a-5346-447f-9aad-4ecd11b6460a
md"""
Finally, we plot $W$ simulated with the optimized initial value and parameter values together with the measured data that was used to find the optimized values.
"""

# ╔═╡ e96efd6a-a666-4120-8480-9423e5d82ae1
# begin
# 	missing
# 	missing
# end
begin
	plot(osol_opt_gom, label="Gompertz growth", xlabel="t",
		xlims=(0, 100), ylims=(0, 14))
	scatter!(t_meas, W_meas, label="Yield")
end

# ╔═╡ f0b4772d-a72b-44e0-a3a1-ba9ad4c4dfeb
md"""
Which grass growth model fits best these data? How can you prove this numerically?
"""

# ╔═╡ 7e98a771-52bb-484e-82ab-2e42e7cb4053
md"- Answer: missing"

# ╔═╡ Cell order:
# ╟─37da8786-fea0-4c2f-a76f-6e6c68325a78
# ╠═a09f814a-0c6a-11ef-0e79-a50b01287d63
# ╠═f8a92690-990b-4341-89e1-322adbcb8d1b
# ╠═015050b3-3339-4b1a-ad7d-c358cce73675
# ╠═a83f424f-1bbb-4e71-af81-302c0ce68907
# ╠═dbfe4800-0974-4ca1-bb0a-d8803409a98b
# ╠═6e227e07-166a-41ce-839a-4c4c72addb23
# ╠═b992c080-a0ce-4188-b632-e734a141e67d
# ╟─4623369d-8c5a-422d-9e40-0f1dd7586260
# ╟─3cb0a166-ac53-4c3f-9832-e93742040cfb
# ╟─987f0a4d-e416-4ceb-adbe-3dcdca9d0996
# ╟─75efff36-8da7-4d04-afa2-a2f8324bc103
# ╟─3dcb9c9d-370b-4031-b7c0-cee80742557a
# ╟─7a14aa59-6e6f-4266-a0b3-84ab55f2efc5
# ╟─85cd60a8-b448-4375-9b6d-399c4336c319
# ╠═5b320989-3e0b-447b-bc9a-25fb221ce609
# ╟─2481cd4f-0efc-4450-ab3d-4a5492597f36
# ╟─1aa44f2b-6f33-437f-b9dd-89762d9f28ea
# ╟─b2b433ed-0266-4bea-a7e8-32adba542d4c
# ╠═7c966a66-0091-4b81-9a7e-02ccd0d3db10
# ╟─3edd2acc-a865-4675-afef-8868c68256f1
# ╠═877298e8-b61b-4c3a-ba2c-2827acdcfb50
# ╟─ef06cc43-510b-4ff9-b0b7-1c7fc267e9b1
# ╠═cb2bc6ee-4211-47e1-9956-5cf1b0c0671d
# ╟─50e6af5b-04f4-495c-853b-c746fd27254f
# ╟─a665740a-ea94-4452-8cbc-a1647319dfce
# ╟─b2f20a9a-0613-4bb4-a1d3-df1620aef105
# ╟─6852f166-350e-4434-8304-382c01113ba1
# ╟─4f78b2e8-9132-4175-b9ae-9f5979649e90
# ╟─557e69c5-d0c6-41ee-b3df-851bd1d16923
# ╟─237934b7-c89a-4ef6-938f-4629c0fa4c95
# ╠═513ee94d-5d49-43d2-a1db-638fe9e2de1a
# ╠═8e163f77-b762-427a-b005-313284dc2e9e
# ╟─d9c38e9d-12b8-4e98-8546-bbc123960e6f
# ╠═39c70e22-86e0-467e-abab-24f99548b096
# ╠═ddef987d-59db-4adf-ae77-d4d32a41ddec
# ╠═54fb706b-9036-42b5-833e-4081f4c2179c
# ╟─9a5bc72b-346d-4e95-a873-783037ed98bc
# ╟─4b17f7e6-26b0-425a-9f53-ccf1689639fe
# ╟─6736542c-5378-480e-a56b-65956b416225
# ╠═ba56adb1-9405-40d5-be48-4273b42ab145
# ╟─1c311b6f-6b18-4170-b500-33a8e4d3cb2d
# ╠═50220e1f-8e42-44da-be03-08c11df967f0
# ╠═449286b0-0210-4c6b-b28f-fca87b52d674
# ╠═f720eac9-cb29-4eef-ac94-471395941c0f
# ╟─d75246d4-e03b-4684-be7d-4bcfb61ed7ef
# ╠═8a9115eb-4044-4cab-a7db-39b5dd86c70d
# ╟─48c9f616-d298-40da-b917-225abd39b3d9
# ╟─35f158c1-858d-4e4d-ac3d-bf4807dad9a0
# ╠═8b6534d6-776b-4285-8498-a9b34051facc
# ╟─a6972aef-63ad-401c-acf5-6d59f9fc6698
# ╟─73e35289-6dc0-4e2e-83eb-b56f83cdbbbf
# ╟─47bd729c-4851-42f7-a03f-6ceacd3c717e
# ╠═f34bb7ac-1ed8-4dd9-b0b9-49bd6e0e1d71
# ╟─e55404ab-6762-4f39-bb42-9c195334a214
# ╠═80e7f6b8-7592-48a3-8587-f1953d1bfcd8
# ╟─a1ca7d0e-639c-42d4-be09-5c61a2008f29
# ╠═30399b9a-1d77-4140-9ad3-5eed636a5b99
# ╠═e1b8e4ba-c1f1-48d5-87a2-edce19c9fe7a
# ╠═a943e0fe-1376-4ef1-9c45-25c6f95e3b96
# ╟─72e065d4-7b1b-4f46-b373-935be8d801fc
# ╟─0b2ffd6f-01cd-4f11-9062-d38b3c13a5b1
# ╠═590b1006-0e37-4668-9b4a-3588fab45696
# ╟─8ff28a2e-185d-4dca-ad3b-a0b1507646d6
# ╠═14d24cbd-3259-41bb-9013-b4fe25a3be4c
# ╟─6c134677-3ec1-4e0a-88b5-01341a096675
# ╠═7a81f4a0-8f7c-4e05-9c3f-2438eab9b691
# ╠═ac80099b-d8f0-4eba-809d-d482bd354d35
# ╟─b58f2c24-e0ea-48a8-b0b7-d0faf9642340
# ╠═55eba435-6ca4-4f4f-b08a-be700d5bda91
# ╟─5d386b00-93b5-4a88-b4e3-e5c3eebd6dd5
# ╟─6f0e91d0-6b99-4cf1-8145-523589a21e89
# ╠═209742ca-36bb-42a5-bf8a-291a40c47757
# ╟─07efb21a-0dac-46d5-b9e1-8f4757c6cedf
# ╠═a93a7eeb-aad1-49f3-a05f-1512adb08b6b
# ╟─ccc66805-efce-4f25-b76b-7cff23901ba4
# ╠═57ef7d03-a6dc-40da-858d-5f9b9be613cd
# ╠═62bd2e58-343a-4155-a65e-cf326d0975f4
# ╠═4d6cba88-be68-412d-aeec-28300c26a9da
# ╟─4bc97af7-debc-4170-a2d5-4e6c928e5183
# ╟─a0033986-4946-45a2-a53e-29c38f35ed5c
# ╠═a570ebab-64de-4934-b887-77a8e2fb42e8
# ╟─feee9cc1-ebe2-4c77-a11e-53ca00864ccf
# ╠═73e86176-f314-463c-8821-51c0f3cc1a56
# ╟─c17c4e00-8688-4f9c-b0ae-741d70601aba
# ╠═f751b962-1810-429c-b34b-658b63fa9ba0
# ╠═0fd058cc-9bd8-4e12-b5dd-79818519165b
# ╟─57c95382-69cc-47e7-aa94-4767136bb23b
# ╠═4e07c39b-502d-4b38-ae7f-f103cb4bae16
# ╟─29170e2a-9916-438e-92ca-9f4783397b5e
# ╟─7ff9fe52-156b-4a92-9058-781670de3abb
# ╠═0c047043-3284-422a-9c88-2f4f4c170edf
# ╠═78a608a4-0ba7-400e-abb9-a9285c60681c
# ╠═19c362cb-2764-41c9-a571-2e8e2bfcde93
# ╠═93db47b2-34e8-43b4-beac-b5620fd444e7
# ╠═68c71cd2-6cdc-4b80-b6e7-75bfea344295
# ╠═d6b9eaba-1d43-4e56-8fa1-bb2f87f6fe79
# ╟─b843ee5e-9618-4b50-93db-77e19b4be366
# ╟─58014411-128c-41f4-b192-b815a3a8cc60
# ╠═ecf4b951-9b5f-441a-89f5-b0ccd040ba02
# ╟─b92bf532-5610-4ce8-b469-154b20cc5956
# ╠═e629fbd1-a7c7-4dc3-8b7d-bdca380fe5ad
# ╟─a48df980-4847-418e-813d-1875d032ffd9
# ╠═348e7454-c6dc-4d1e-b18b-6b34bc2fc0fc
# ╠═f71f0cf4-b63e-452c-b971-f30d7503f1e5
# ╟─7346a003-aa05-41e3-a420-c97c1994cf3c
# ╠═20a6f165-449a-42cc-9565-f342a7535422
# ╟─137bde23-76f2-4ebf-8bc2-ea8640001436
# ╟─4aa71200-006b-4a15-ae75-67e36aa81522
# ╟─cdab3079-04b0-4a44-b770-468c20e321e4
# ╠═cf1a144e-09e9-42a3-b2a3-b8676a200a39
# ╟─d69ed5dd-80cd-427e-9245-e42576a4688d
# ╠═8b429656-1c57-4278-9633-71311f400ed8
# ╟─b4300e8a-8052-419b-98c8-0508ebee2393
# ╠═2c6ae74c-2da4-4867-ad8a-f4e835101d63
# ╟─c6a5d453-d610-4f65-847c-c878dd41726c
# ╠═fff17cf7-173d-4f64-94a9-4bf46acc882d
# ╟─eee55784-a641-445e-be75-0b19e2a94754
# ╠═7844e4f5-3c7d-4b4b-beee-970c998c67a6
# ╟─c81d0140-3f4e-4eb4-8a77-1f48c5e0ecbf
# ╠═7456455b-4f31-488f-990f-6ce534038e08
# ╟─b10c2ce4-d363-429c-a64c-ec29652137a5
# ╠═23b629a6-6866-416a-a768-9617ce6301db
# ╠═8788082d-f5d1-4385-8037-a0d360a841c7
# ╠═03a4fa85-08db-46d4-bb53-c0ccea90a211
# ╟─8f5e1413-227c-44fa-bb2d-3653cbc27e38
# ╟─8a7f7aab-878e-41b5-b9da-d06747df042e
# ╠═30602ff1-041b-4fca-bf8e-55ff57df9e37
# ╟─881011be-6434-416f-915b-3333e8dea32f
# ╠═5602b88a-07f8-438b-994c-65f11e17a0ba
# ╟─25be4255-0888-4ecd-a2fd-d66402c5cb50
# ╠═36e8a174-d526-45ee-b3c6-88d698ad5d5f
# ╠═7594147d-b3da-4e0d-896d-41baacb6d7be
# ╟─8e047be9-f0f0-4a75-91b5-c523f55f8c67
# ╠═1a9587aa-2356-48df-abcc-2ce874fa5d24
# ╟─785d500b-f8ea-446a-9952-2a5fd5d83d24
# ╟─e754826a-7411-4072-b0dc-a4bad7a15f98
# ╠═6e85c08e-9a97-4897-9c81-89517f55e254
# ╟─e5081280-d226-4834-8932-c89becd8313c
# ╠═c739a908-2353-4e7a-8fbd-f640dc8cabe0
# ╟─1d0383ad-54d6-4ff2-8555-def83bfff0e6
# ╠═cd1cf2f8-9f7f-4ed4-9cb7-1a6efee68ab4
# ╟─aba74ee0-0163-4e15-8b49-d8dcad4839f7
# ╠═0eda4142-1aaf-4e17-bd78-857e13e94acd
# ╟─50629194-98ed-4d45-86a2-95ac22daac29
# ╠═9c239fc6-275c-4d64-9fa2-6fd57295b757
# ╟─dede17f2-655d-4871-b6de-5a32804947dd
# ╠═e8c3f042-6058-4040-a67e-18563c04ee93
# ╠═665f4d03-3521-475a-a195-f861fd26bb69
# ╠═13775d58-ac61-431c-a6c9-447c1eec7942
# ╟─bc92f996-b626-4965-a286-d2c848eb1a21
# ╟─e2fde8e9-1f87-4ffe-8dae-2794664bfaa4
# ╠═23fd9fd9-a13c-4c56-a0b4-daec6f1d2cd8
# ╟─6f414d15-af0a-452a-98a1-dc0b7e54d617
# ╠═57ee8a12-24df-4598-935c-f5e259b504cb
# ╟─48bc085c-9ce6-4752-a5a9-a814f803f571
# ╠═5b9b2e5b-9deb-4ff6-a923-b15b8b08f0c9
# ╠═66ee9655-a006-4c0f-b1f1-5576efa8f896
# ╟─52d7975a-5346-447f-9aad-4ecd11b6460a
# ╠═e96efd6a-a666-4120-8480-9423e5d82ae1
# ╟─f0b4772d-a72b-44e0-a3a1-ba9ad4c4dfeb
# ╠═7e98a771-52bb-484e-82ab-2e42e7cb4053
