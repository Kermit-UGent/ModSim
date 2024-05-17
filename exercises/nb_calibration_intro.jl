### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ f8a92690-990b-4341-89e1-322adbcb8d1b
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ a09f814a-0c6a-11ef-0e79-a50b01287d63
using Markdown

# ╔═╡ 7f521435-63ac-4178-a4aa-93d9c45fe820
using InteractiveUtils

# ╔═╡ 015050b3-3339-4b1a-ad7d-c358cce73675
using Catalyst, DifferentialEquations, Plots

# ╔═╡ dbfe4800-0974-4ca1-bb0a-d8803409a98b
using Turing

# ╔═╡ 7965b69b-fff6-4284-ad65-a92a58cda04a
using StatsPlots, StatsBase

# ╔═╡ 6e227e07-166a-41ce-839a-4c4c72addb23
using LinearAlgebra

# ╔═╡ 9f20022f-d74b-42ce-a751-4ea94eb896fe
using Optim

# ╔═╡ 37da8786-fea0-4c2f-a76f-6e6c68325a78
md"
# Introduction to calibration
"

# ╔═╡ 4623369d-8c5a-422d-9e40-0f1dd7586260
md"
## Goal of this practicum
"

# ╔═╡ 89b3701a-dc0c-4e8d-ba6b-b02e218ff79b
md"
In the models discussed in the previous sessions, we always knew the values of all parameters. In reality, the value of a parameter has to be calibrated, hence, estimated from experimental data. During this parameter estimation one attempts to find the set of parameter values for which the model predictions are as close as possible to the collected experimental data.

The search involves an objective function $J(\boldsymbol{\theta})$. This represents a measure for the distance (i.e. a single number, a scaler) between the model predictions and the data, as measured by the **S**um of **S**quared **R**esiduals (SSR):

$$\begin{equation}
J(\boldsymbol{\theta})=\frac{1}{\sigma^{2}_{1}}\sum_{j=1}^{N_1}\left(y_{1,j}-\hat{y}_{1,j}(\boldsymbol{\theta})\right)^{2}+\frac{1}{\sigma^{2}_{2}}\sum_{j=1}^{N_2}\left(y_{2,j}-\hat{y}_{2,j}(\boldsymbol{\theta})\right)^{2}+\ldots + \frac{1}{\sigma^{2}_{v}}\sum_{j=1}^{N_v}\left(y_{v, j}-\hat{y}_{v,j}(\boldsymbol{\theta})\right)^{2}.
\end{equation}$$

Here $y_{i,j}$ represents the $j$-th experimental data point for the $i$-th output variable; $\hat{y}_{i,j}(\boldsymbol{\theta})$ represents the model prediction associated with this data point, evaluated with the parameter values $\boldsymbol{\theta}$. The difference $(y_{i,j}-\hat{y}_{i,j}(\boldsymbol{\theta}))$ is called a *residual*. The number of data points collected for the $i$-th output variable is $N_i$, and the measurement error associated with the data collected for the $i$-th output variable is $\sigma_i$. The terms $1/\sigma_i^2$ thus ensure that the contribution of the different output variables to the objective function is weighted by the precision with which these output variables can be measured.

Thus, the estimation of the parameters is reduced to finding the set of parameter values $\hat{\boldsymbol{\theta}}$ for which the objective function is minimized, or mathematically:

$$\begin{equation*}
\hat{\boldsymbol{\theta}} = \textrm{argmin}_{\boldsymbol{\theta}} J(\boldsymbol{\theta}).
\end{equation*}$$

meaning that $\hat{\boldsymbol{\theta}}$ is the argument $\boldsymbol{\theta}$ that minimizes the objective function $J(\boldsymbol{\theta})$.
"

# ╔═╡ 75efff36-8da7-4d04-afa2-a2f8324bc103
md"
In this notebook we will calibrate the different parameters involved in the grass growth models. Therefore, we will implement an objective function for each model and then minimizing them using an optimization algorithm. To illustrate this concept, we first revisit the three simple models modelling the grass growth yield.
"

# ╔═╡ 3dcb9c9d-370b-4031-b7c0-cee80742557a
md"
## Grass growth models
"

# ╔═╡ 7a14aa59-6e6f-4266-a0b3-84ab55f2efc5
md"
In this notebook, three different models will be used, each modelling the yield of grass in a grassland:

- Logistic growth model: $\cfrac{dW}{dt} = \mu \left( 1 - \cfrac{W}{W_f} \right) W$
- Exponential growth model: $\cfrac{dW}{dt} = \mu \left( W_f - W \right)$
- Gompertz growth model: $\cfrac{dW}{dt} = \left( \mu - D \ln(W) \right) W$

with output $W$ the grass yield, and $W_f$, $\mu$ and $D$ parameters. The table below show the parameter values and the initial condition that will be used as initial values for the optimization algorithm:

|             | $\mu$      | $W_f$       | $D$          | $W_0$          |
|:----------- |:----------:|:-----------:|:------------:|:------------:|
| Logistic    |  0.07      | 10.0        |              | 2.0          |
| Exponential |  0.02      | 10.0        |              | 2.0          |
| Gompertz    |  0.09      |             | 0.040        | 2.0          |

Hence, for each grass growth model, we will optimize the parameter values together with the initial value.
"

# ╔═╡ 85cd60a8-b448-4375-9b6d-399c4336c319
md"
In each of the three models we will use the following timespan:
"

# ╔═╡ 2481cd4f-0efc-4450-ab3d-4a5492597f36
md"
Variables containing the initial condition and parameters values will be defined later in the objective function.
"

# ╔═╡ 9a5bc72b-346d-4e95-a873-783037ed98bc
md"
### Logistic growth model

We will illustrate the calibration with the logistic growth model:
"

# ╔═╡ ba56adb1-9405-40d5-be48-4273b42ab145
growth_log = @reaction_network begin
    μ*W, ∅ --> W
    μ/Wf*W, W --> ∅
end

# ╔═╡ f4748167-b635-47a1-9015-32e1258c0afa
md"
Check the order of the parameters:
"

# ╔═╡ a022b2ab-68a0-40ca-b914-7a2adcf4ae39
parameters(growth_log)

# ╔═╡ acccb2fa-12b2-4fc7-91e3-58a4b1a02892
md"
Next, we will need to create an `ODEProblem` in advance before we can optimize some of its parameters. We will provide the values in the aforementioned table as initial values for the problem.
"

# ╔═╡ e54ea8d1-0854-44fa-aed8-45d106e921e4
u₀_log = [:W => 2.0]

# ╔═╡ 8d96eb17-ce19-4523-916f-3cd0441a16ca
params_log = [:μ => 0.07, :Wf => 10.0]

# ╔═╡ 5c9db9df-0cbd-41ac-afe9-fb5616c967be
oprob_log = ODEProblem(growth_log, u₀_log, tspan, params_log)

# ╔═╡ 1aa44f2b-6f33-437f-b9dd-89762d9f28ea
md"
### The measurement data
"

# ╔═╡ b2b433ed-0266-4bea-a7e8-32adba542d4c
md"
Assume that the measured grass yields (of a certain plant type) are the following:
"

# ╔═╡ 7c966a66-0091-4b81-9a7e-02ccd0d3db10
W_meas = [1.87, 2.45, 3.72, 4.32, 5.28, 7.01, 6.83, 8.62, 9.45, 10.31, 10.56, 11.72, 11.05, 11.53, 11.39, 11.7, 11.15, 11.49, 12.04, 11.95, 11.68]

# ╔═╡ 3edd2acc-a865-4675-afef-8868c68256f1
md"
They have been measured at the following corresponding time instances:
"

# ╔═╡ 877298e8-b61b-4c3a-ba2c-2827acdcfb50
t_meas = 0:5:100

# ╔═╡ ef06cc43-510b-4ff9-b0b7-1c7fc267e9b1
md"
We can make a scatter plot of this data (including a title, a legend label, an X-axis label, X- and Y-axis limits) in the following way:
"

# ╔═╡ cb2bc6ee-4211-47e1-9956-5cf1b0c0671d
scatter(t_meas, W_meas, title="Grass growth data",
                        label="Yield",
                        xlabel="t",
                        xlims=(0, 100),
                        ylims=(0, 14))

# ╔═╡ d75246d4-e03b-4684-be7d-4bcfb61ed7ef
md"
### Declaration of the Turing model
"

# ╔═╡ dd3a32f1-bdb6-44a3-acbe-f4269725c9e4
md"
In the Turing model we will define our priors for the following magnitudes:
- the measurement error (standard deviation) $\sigma_W$, 
- the initial condition $W_0$, and
- the parameters $\mu$ and $W_f$.

We will thereby take an Inverse Gamma prior distribution for $\sigma_W$ and Uniform prior distributions for the initial condition and the parameters (assuming they lay in pre-defined range).
"

# ╔═╡ 8a9115eb-4044-4cab-a7db-39b5dd86c70d
@model function growth_log_inference(t_meas, W)
    σ_W ~ InverseGamma()
    W₀ ~ Uniform(0, 10)
    μ ~ Uniform(0, 1)
    Wf ~ Uniform(0, 100)
    osol_log = solve(remake(oprob_log; u0=[W₀]), Tsit5(), saveat=t_meas, p=[μ, Wf])
    W ~ MvNormal(osol_log[:W], σ_W^2 * I)
end

# ╔═╡ 48c9f616-d298-40da-b917-225abd39b3d9
md"
Some remarks:
- If one want to optimize the initial condition(s), the optimization algorithm will need to overwrite the initial condition, therefore, we use `remake(oprob_log; u0=[W₀])`. Note the **semi colon** between both arguments! If the initial condition is not being optimized you can just write `oprob_log`.
- The time points are the ones from the measurements, therefore, we set: `saveat=t_meas`.
- The same here, the optimization algorithm will need to overwrite the parameters, therefore we need to place all parameters in the option `p=[...]` in the `solve` function in the same order as in the reaction network model.
"

# ╔═╡ 35f158c1-858d-4e4d-ac3d-bf4807dad9a0
md"
We will provide the measurements to the Turing model:
"

# ╔═╡ 8b6534d6-776b-4285-8498-a9b34051facc
growth_log_inf = growth_log_inference(t_meas, W_meas)

# ╔═╡ a6972aef-63ad-401c-acf5-6d59f9fc6698
md"
We are now ready to optimize the priors ($\sigma_W$, $W_0$, $\mu$ and $W_f$). This is done by calling the `optimize` function, providing the previously created object `growth_log_inf`, the method for estimating the parameters and (optionally) an algorithm (default: Nelder-Mead) to implement the method.
"

# ╔═╡ 73e35289-6dc0-4e2e-83eb-b56f83cdbbbf
md"
#### Method - Maximum Likelihood Estimation
"

# ╔═╡ 47bd729c-4851-42f7-a03f-6ceacd3c717e
md"
We will use the MLE (Maximum Likelihood Estimation) method here and store the optimization results in `results_log_mle`.
"

# ╔═╡ f34bb7ac-1ed8-4dd9-b0b9-49bd6e0e1d71
results_log_mle = optimize(growth_log_inf, MLE(), NelderMead())

# ╔═╡ e55404ab-6762-4f39-bb42-9c195334a214
md"
You can visualize a summary of the optimized parameters by piping them to `coeftable`:
"

# ╔═╡ 80e7f6b8-7592-48a3-8587-f1953d1bfcd8
results_log_mle |> coeftable

# ╔═╡ a1ca7d0e-639c-42d4-be09-5c61a2008f29
md"
You can obtain the actual optimized values using the function `coef` on the results object in conjunction by calling the parameters by name preceded by a colon. Here we assign the optimized parameter values to some suitable variable names:
"

# ╔═╡ 30399b9a-1d77-4140-9ad3-5eed636a5b99
W₀_opt = coef(results_log_mle)[:W₀]

# ╔═╡ e1b8e4ba-c1f1-48d5-87a2-edce19c9fe7a
μ_opt = coef(results_log_mle)[:μ]

# ╔═╡ a943e0fe-1376-4ef1-9c45-25c6f95e3b96
Wf_opt = coef(results_log_mle)[:Wf]

# ╔═╡ 72e065d4-7b1b-4f46-b373-935be8d801fc
md"
Now we can make a plot of $W$ simulated with the optimized initial condition and parameter values.
"

# ╔═╡ 0b2ffd6f-01cd-4f11-9062-d38b3c13a5b1
md"
Setting up initial condition with optimized initial condition:
"

# ╔═╡ 590b1006-0e37-4668-9b4a-3588fab45696
u₀_opt_log = [:W => W₀_opt]

# ╔═╡ 8ff28a2e-185d-4dca-ad3b-a0b1507646d6
md"
Setting up parameter values with optimized parameter values:
"

# ╔═╡ 14d24cbd-3259-41bb-9013-b4fe25a3be4c
params_opt_log = [:μ => μ_opt, :Wf => Wf_opt]

# ╔═╡ 6c134677-3ec1-4e0a-88b5-01341a096675
md"
Next, we create an ODEProblem and solve it:
"

# ╔═╡ 7a81f4a0-8f7c-4e05-9c3f-2438eab9b691
oprob_opt_log = ODEProblem(growth_log, u₀_opt_log, tspan, params_opt_log)

# ╔═╡ ac80099b-d8f0-4eba-809d-d482bd354d35
osol_opt_log = solve(oprob_opt_log, Tsit5(), saveat=0.5)

# ╔═╡ b58f2c24-e0ea-48a8-b0b7-d0faf9642340
md"
Finally, we plot $W$ simulated with the optimized initial value and parameter values together with the measured data that was used to find the optimized values.
"

# ╔═╡ 55eba435-6ca4-4f4f-b08a-be700d5bda91
begin
plot(osol_opt_log, label="Logistic growth", xlabel="t", xlims=(0, 100), ylims=(0, 14))
scatter!(t_meas, W_meas, label="Yield")
end

# ╔═╡ 5d386b00-93b5-4a88-b4e3-e5c3eebd6dd5
md"
#### Method - Maximum A Posterior
"

# ╔═╡ 6f0e91d0-6b99-4cf1-8145-523589a21e89
md"
To be continued.
"

# ╔═╡ 29170e2a-9916-438e-92ca-9f4783397b5e
md"
#### Method - No-U-Turn Sampler
"

# ╔═╡ 7ff9fe52-156b-4a92-9058-781670de3abb
md"
To be continued.
"

# ╔═╡ 137bde23-76f2-4ebf-8bc2-ea8640001436
md"
## Exercises
"

# ╔═╡ 4aa71200-006b-4a15-ae75-67e36aa81522
md"
### Exercise 1 - Calibration of the exponential growth model

Calibrate the initial condition and both parameters of the exponential growth model. Use the values mentioned in the *Table* as initials values for the optimization of the parameters.
"

# ╔═╡ cdab3079-04b0-4a44-b770-468c20e321e4
md"
We have seen before that a possible *reaction network object* for the exponential growth model can be implemented as follows:
"

# ╔═╡ cf1a144e-09e9-42a3-b2a3-b8676a200a39
growth_exp = @reaction_network begin
    μ*Wf, ∅ --> W
    μ, W --> ∅
end

# ╔═╡ c6d373f4-f13c-4135-823d-ee8fbeb71b56
md"
Create an `ODEProblem`. Use the values in the aforementioned table as initial values for the problem. Use the same `tspan` as before.
"

# ╔═╡ a97abaa7-b642-4201-86f1-5c8995b07536
# u₀_exp = missing        # Uncomment and complete the instruction
u₀_exp = [:W => 2.0]

# ╔═╡ 387730b4-bd06-492f-94e6-231bd68b3436
# params_exp = missing    # Uncomment and complete the instruction
params_exp = [:μ => 0.02, :Wf => 10.0]

# ╔═╡ 290a7fe8-3b1e-423f-8b30-9bd8903d2e8f
# oprob_exp = missing     # Uncomment and complete the instruction
oprob_exp = ODEProblem(growth_exp, u₀_exp, tspan, params_exp)

# ╔═╡ febe2b67-2a8f-4575-946d-30877bd5f2d4
md"
Use the same measurement data (`W_meas`, `t_meas`) as before.
"

# ╔═╡ b4300e8a-8052-419b-98c8-0508ebee2393
md"
Declare the Turing model. Take the same priors (and distributions) as before.
"

# ╔═╡ 2c6ae74c-2da4-4867-ad8a-f4e835101d63
# Uncomment and complete the instruction
# @model function growth_exp_inference(t_meas, W)
#     σ_W ~ missing
#     W₀ ~ missing
#     μ ~ missing
#     Wf ~ missing
#     osol_exp = missing
#     W ~ missing
# end
@model function growth_exp_inference(t_meas, W)
    σ_W ~ InverseGamma()
    W₀ ~ Uniform(0, 10)
    μ ~ Uniform(0, 1)
    Wf ~ Uniform(0, 100)
    osol_exp = solve(remake(oprob_exp; u0=[W₀]), Tsit5(), saveat=t_meas, p=[μ, Wf])
    W ~ MvNormal(osol_exp[:W], σ_W^2 * I)
end

# ╔═╡ c6a5d453-d610-4f65-847c-c878dd41726c
md"
Provide the measurements to the Turing model.
"

# ╔═╡ fff17cf7-173d-4f64-94a9-4bf46acc882d
# growth_exp_inf = missing           # Uncomment and complete the instruction
growth_exp_inf = growth_exp_inference(t_meas, W_meas)

# ╔═╡ eee55784-a641-445e-be75-0b19e2a94754
md"
Optimize the priors ($\sigma_W$, $W_0$, $\mu$ and $W_f$). Do this with `MLE` method and Nelder-Mead. Store the optimization results in `results_exp_mle`.
"

# ╔═╡ 7844e4f5-3c7d-4b4b-beee-970c998c67a6
# results_exp_mle = missing          # Uncomment and complete the instruction
results_exp_mle = optimize(growth_exp_inf, MLE(), NelderMead())

# ╔═╡ c81d0140-3f4e-4eb4-8a77-1f48c5e0ecbf
md"
Visualize a summary of the optimized parameters.
"

# ╔═╡ 7456455b-4f31-488f-990f-6ce534038e08
# missing              # Uncomment and complete the instruction
results_exp_mle |> coeftable

# ╔═╡ b10c2ce4-d363-429c-a64c-ec29652137a5
md"
Get the optimized values and assign them to `W₀_opt_exp`, `μ_opt_exp` and `Wf_opt_exp`.
"

# ╔═╡ 23b629a6-6866-416a-a768-9617ce6301db
# W₀_opt_exp = missing            # Uncomment and complete the instruction
W₀_opt_exp = coef(results_exp_mle)[:W₀]

# ╔═╡ 8788082d-f5d1-4385-8037-a0d360a841c7
# μ_opt_exp = missing             # Uncomment and complete the instruction
μ_opt_exp = coef(results_exp_mle)[:μ]

# ╔═╡ 03a4fa85-08db-46d4-bb53-c0ccea90a211
# Wf_opt_exp = missing            # Uncomment and complete the instruction
Wf_opt_exp = coef(results_exp_mle)[:Wf]

# ╔═╡ 8f5e1413-227c-44fa-bb2d-3653cbc27e38
md"
Make a plot of $W$ simulated with the optimized initial condition and parameter values.
"

# ╔═╡ 8a7f7aab-878e-41b5-b9da-d06747df042e
md"
Set up initial condition with optimized initial condition:
"

# ╔═╡ 30602ff1-041b-4fca-bf8e-55ff57df9e37
# u₀_opt_exp = missing                 # Uncomment and complete the instruction
u₀_opt_exp = [:W => W₀_opt_exp]

# ╔═╡ 881011be-6434-416f-915b-3333e8dea32f
md"
Set up parameter values with optimized parameter values:
"

# ╔═╡ 5602b88a-07f8-438b-994c-65f11e17a0ba
# params_opt_exp = missing       # Uncomment and complete the instruction
params_opt_exp = [:μ => μ_opt_exp, :Wf => Wf_opt_exp]

# ╔═╡ 25be4255-0888-4ecd-a2fd-d66402c5cb50
md"
Create an ODEProblem and solve it:
"

# ╔═╡ 36e8a174-d526-45ee-b3c6-88d698ad5d5f
# oprob_opt_exp = missing        # Uncomment and complete the instruction
oprob_opt_exp = ODEProblem(growth_exp, u₀_opt_exp, tspan, params_opt_exp)

# ╔═╡ 7594147d-b3da-4e0d-896d-41baacb6d7be
# osol_opt_exp = missing       # Uncomment and complete the instruction
osol_opt_exp = solve(oprob_opt_exp, Tsit5(), saveat=0.5)

# ╔═╡ 8e047be9-f0f0-4a75-91b5-c523f55f8c67
md"
Plot $W$ simulated with the optimized initial value and parameter values together with the measured data that was used to find the optimized values.
"

# ╔═╡ 1a9587aa-2356-48df-abcc-2ce874fa5d24
# Uncomment and complete the instruction
# begin
# missing
# missing
# end
begin
plot(osol_opt_exp, label="Exponential growth", xlabel="t", xlims=(0, 100), ylims=(0, 14))
scatter!(t_meas, W_meas, label="Yield")
end

# ╔═╡ 785d500b-f8ea-446a-9952-2a5fd5d83d24
md"
### Exercise 2 - Calibration of the Gompertz growth model

Calibrate the initial condition and both parameters of the Gompertz growth model. Use the values mentioned in the *Table* as initials values for the optimization of the parameters.
"

# ╔═╡ e754826a-7411-4072-b0dc-a4bad7a15f98
md"
We have seen before that a possible *reaction network object* for the Gompertz growth model can be implemented as follows:
"

# ╔═╡ bc1edcbe-46eb-4531-9c5f-dee8d5dc2ff9
growth_gom = @reaction_network begin
    -μ, W --> ∅
    D*log(W), W --> ∅
end

# ╔═╡ 47fb9e4c-df6a-4811-9980-99d595a34908
md"
Create an `ODEProblem`. Use the values in the aforementioned table as initial values for the problem. Use the same `tspan` as before.
"

# ╔═╡ bbd150de-ff9a-4127-a0dd-2f9762f92b07
# u₀_gom = missing           # Uncomment and complete the instruction
u₀_gom = [:W => 2.0]

# ╔═╡ da5a0cbb-b033-46f1-a300-3954de138835
# params_gom = missing       # Uncomment and complete the instruction
params_gom = [:μ => 0.09, :D => 0.04]

# ╔═╡ 73da8f53-c3af-43b0-9b23-60471f1e3587
# oprob_gom = missing        # Uncomment and complete the instruction
oprob_gom = ODEProblem(growth_gom, u₀_gom, tspan, params_gom)

# ╔═╡ b0e67564-efe8-4fb2-bcf2-a711b770244e
md"
Use the same measurement data (`W_meas`, `t_meas`) as before.
"

# ╔═╡ e5081280-d226-4834-8932-c89becd8313c
md"
Declare the Turing model. Take for $\sigma_W$ and $W_0$ the same priors (and distributions) as before, but take for $\mu$ a Uniform prior distribution in the range $[0, 2]$ and the same for $D$ but in the range $[0, 1]$.
"

# ╔═╡ c739a908-2353-4e7a-8fbd-f640dc8cabe0
# Uncomment and complete the instruction
# @model function growth_gom_inference(t_meas, W)
#     σ_W ~ missing
#     W₀ ~ missing
#     μ ~ missing
#     D ~ missing
#     osol_gom = missing
#     W ~ missing
# end
@model function growth_gom_inference(t_meas, W)
    σ_W ~ InverseGamma()
    W₀ ~ Uniform(0, 10)
    μ ~ Uniform(0, 2)
    D ~ Uniform(0, 1)
    osol_gom = solve(remake(oprob_gom; u0=[W₀]), Tsit5(), saveat=t_meas, p=[μ, D])
    W ~ MvNormal(osol_gom[:W], σ_W^2 * I)
end

# ╔═╡ 1d0383ad-54d6-4ff2-8555-def83bfff0e6
md"
Provide the measurements to the Turing model.
"

# ╔═╡ cd1cf2f8-9f7f-4ed4-9cb7-1a6efee68ab4
# growth_gom_inf = missing         # Uncomment and complete the instruction
growth_gom_inf = growth_gom_inference(t_meas, W_meas)

# ╔═╡ aba74ee0-0163-4e15-8b49-d8dcad4839f7
md"
Optimize the priors ($\sigma_W$, $W_0$, $\mu$ and $D$). Do this with `MLE` method and Nelder-Mead. Store the optimization results in `results_gom_mle`.
"

# ╔═╡ 0eda4142-1aaf-4e17-bd78-857e13e94acd
# results_gom_mle = missing             # Uncomment and complete the instruction
results_gom_mle = optimize(growth_gom_inf, MLE(), NelderMead())

# ╔═╡ 50629194-98ed-4d45-86a2-95ac22daac29
md"
Visualize a summary of the optimized parameters.
"

# ╔═╡ 9c239fc6-275c-4d64-9fa2-6fd57295b757
# missing                       # Uncomment and complete the instruction
results_gom_mle |> coeftable

# ╔═╡ dede17f2-655d-4871-b6de-5a32804947dd
md"
Get the optimized values and assign them to `W₀_opt_gom`, `μ_opt_gom` and `D_opt_gom`.
"

# ╔═╡ e8c3f042-6058-4040-a67e-18563c04ee93
# W₀_opt_gom = missing          # Uncomment and complete the instruction
W₀_opt_gom = coef(results_gom_mle)[:W₀]

# ╔═╡ 665f4d03-3521-475a-a195-f861fd26bb69
# μ_opt_gom = missing           # Uncomment and complete the instruction
μ_opt_gom = coef(results_gom_mle)[:μ]

# ╔═╡ 13775d58-ac61-431c-a6c9-447c1eec7942
# D_opt_gom = missing           # Uncomment and complete the instruction
D_opt_gom = coef(results_gom_mle)[:D]

# ╔═╡ bc92f996-b626-4965-a286-d2c848eb1a21
md"
Make a plot of $W$ simulated with the optimized initial condition and parameter values.
"

# ╔═╡ e2fde8e9-1f87-4ffe-8dae-2794664bfaa4
md"
Set up initial condition with optimized initial condition:
"

# ╔═╡ 23fd9fd9-a13c-4c56-a0b4-daec6f1d2cd8
# u₀_opt_gom = missing          # Uncomment and complete the instruction
u₀_opt_gom = [:W => W₀_opt_gom]

# ╔═╡ 6f414d15-af0a-452a-98a1-dc0b7e54d617
md"
Set up parameter values with optimized parameter values:
"

# ╔═╡ 57ee8a12-24df-4598-935c-f5e259b504cb
# params_opt_gom = missing     # Uncomment and complete the instruction
params_opt_gom = [:μ => μ_opt_gom, :D => D_opt_gom]

# ╔═╡ 48bc085c-9ce6-4752-a5a9-a814f803f571
md"
Create an ODEProblem and solve it:
"

# ╔═╡ 5b9b2e5b-9deb-4ff6-a923-b15b8b08f0c9
# oprob_opt_gom = missing      # Uncomment and complete the instruction
oprob_opt_gom = ODEProblem(growth_gom, u₀_opt_gom, tspan, params_opt_gom)

# ╔═╡ 66ee9655-a006-4c0f-b1f1-5576efa8f896
# osol_opt_gom = missing        # Uncomment and complete the instruction
osol_opt_gom = solve(oprob_opt_gom, Tsit5(), saveat=0.5)

# ╔═╡ 52d7975a-5346-447f-9aad-4ecd11b6460a
md"
Finally, we plot $W$ simulated with the optimized initial value and parameter values together with the measured data that was used to find the optimized values.
"

# ╔═╡ e96efd6a-a666-4120-8480-9423e5d82ae1
# Uncomment and complete the instruction
# begin
# missing
# missing
# end
begin
plot(osol_opt_gom, label="Gompertz growth", xlabel="t", xlims=(0, 100), ylims=(0, 14))
scatter!(t_meas, W_meas, label="Yield")
end

# ╔═╡ f0b4772d-a72b-44e0-a3a1-ba9ad4c4dfeb
md"
Which grass growth model fits best these data? How can you prove this numerically?
- Answer: missing
"

# ╔═╡ 5b320989-3e0b-447b-bc9a-25fb221ce609
# ╠═╡ disabled = true
#=╠═╡
tspan = (0.0, 100.0)   # this will be the same for the three models
  ╠═╡ =#

# ╔═╡ 5f9005c7-5574-4717-81c8-d725a0fb2692
tspan = (0.0, 100.0)

# ╔═╡ Cell order:
# ╠═a09f814a-0c6a-11ef-0e79-a50b01287d63
# ╠═7f521435-63ac-4178-a4aa-93d9c45fe820
# ╠═f8a92690-990b-4341-89e1-322adbcb8d1b
# ╠═015050b3-3339-4b1a-ad7d-c358cce73675
# ╠═dbfe4800-0974-4ca1-bb0a-d8803409a98b
# ╠═7965b69b-fff6-4284-ad65-a92a58cda04a
# ╠═6e227e07-166a-41ce-839a-4c4c72addb23
# ╠═9f20022f-d74b-42ce-a751-4ea94eb896fe
# ╠═37da8786-fea0-4c2f-a76f-6e6c68325a78
# ╠═4623369d-8c5a-422d-9e40-0f1dd7586260
# ╠═89b3701a-dc0c-4e8d-ba6b-b02e218ff79b
# ╠═75efff36-8da7-4d04-afa2-a2f8324bc103
# ╠═3dcb9c9d-370b-4031-b7c0-cee80742557a
# ╠═7a14aa59-6e6f-4266-a0b3-84ab55f2efc5
# ╠═85cd60a8-b448-4375-9b6d-399c4336c319
# ╠═5b320989-3e0b-447b-bc9a-25fb221ce609
# ╠═2481cd4f-0efc-4450-ab3d-4a5492597f36
# ╠═9a5bc72b-346d-4e95-a873-783037ed98bc
# ╠═ba56adb1-9405-40d5-be48-4273b42ab145
# ╠═f4748167-b635-47a1-9015-32e1258c0afa
# ╠═a022b2ab-68a0-40ca-b914-7a2adcf4ae39
# ╠═acccb2fa-12b2-4fc7-91e3-58a4b1a02892
# ╠═e54ea8d1-0854-44fa-aed8-45d106e921e4
# ╠═5f9005c7-5574-4717-81c8-d725a0fb2692
# ╠═8d96eb17-ce19-4523-916f-3cd0441a16ca
# ╠═5c9db9df-0cbd-41ac-afe9-fb5616c967be
# ╠═1aa44f2b-6f33-437f-b9dd-89762d9f28ea
# ╠═b2b433ed-0266-4bea-a7e8-32adba542d4c
# ╠═7c966a66-0091-4b81-9a7e-02ccd0d3db10
# ╠═3edd2acc-a865-4675-afef-8868c68256f1
# ╠═877298e8-b61b-4c3a-ba2c-2827acdcfb50
# ╠═ef06cc43-510b-4ff9-b0b7-1c7fc267e9b1
# ╠═cb2bc6ee-4211-47e1-9956-5cf1b0c0671d
# ╠═d75246d4-e03b-4684-be7d-4bcfb61ed7ef
# ╠═dd3a32f1-bdb6-44a3-acbe-f4269725c9e4
# ╠═8a9115eb-4044-4cab-a7db-39b5dd86c70d
# ╠═48c9f616-d298-40da-b917-225abd39b3d9
# ╠═35f158c1-858d-4e4d-ac3d-bf4807dad9a0
# ╠═8b6534d6-776b-4285-8498-a9b34051facc
# ╠═a6972aef-63ad-401c-acf5-6d59f9fc6698
# ╠═73e35289-6dc0-4e2e-83eb-b56f83cdbbbf
# ╠═47bd729c-4851-42f7-a03f-6ceacd3c717e
# ╠═f34bb7ac-1ed8-4dd9-b0b9-49bd6e0e1d71
# ╠═e55404ab-6762-4f39-bb42-9c195334a214
# ╠═80e7f6b8-7592-48a3-8587-f1953d1bfcd8
# ╠═a1ca7d0e-639c-42d4-be09-5c61a2008f29
# ╠═30399b9a-1d77-4140-9ad3-5eed636a5b99
# ╠═e1b8e4ba-c1f1-48d5-87a2-edce19c9fe7a
# ╠═a943e0fe-1376-4ef1-9c45-25c6f95e3b96
# ╠═72e065d4-7b1b-4f46-b373-935be8d801fc
# ╠═0b2ffd6f-01cd-4f11-9062-d38b3c13a5b1
# ╠═590b1006-0e37-4668-9b4a-3588fab45696
# ╠═8ff28a2e-185d-4dca-ad3b-a0b1507646d6
# ╠═14d24cbd-3259-41bb-9013-b4fe25a3be4c
# ╠═6c134677-3ec1-4e0a-88b5-01341a096675
# ╠═7a81f4a0-8f7c-4e05-9c3f-2438eab9b691
# ╠═ac80099b-d8f0-4eba-809d-d482bd354d35
# ╠═b58f2c24-e0ea-48a8-b0b7-d0faf9642340
# ╠═55eba435-6ca4-4f4f-b08a-be700d5bda91
# ╠═5d386b00-93b5-4a88-b4e3-e5c3eebd6dd5
# ╠═6f0e91d0-6b99-4cf1-8145-523589a21e89
# ╠═29170e2a-9916-438e-92ca-9f4783397b5e
# ╠═7ff9fe52-156b-4a92-9058-781670de3abb
# ╠═137bde23-76f2-4ebf-8bc2-ea8640001436
# ╠═4aa71200-006b-4a15-ae75-67e36aa81522
# ╠═cdab3079-04b0-4a44-b770-468c20e321e4
# ╠═cf1a144e-09e9-42a3-b2a3-b8676a200a39
# ╠═c6d373f4-f13c-4135-823d-ee8fbeb71b56
# ╠═a97abaa7-b642-4201-86f1-5c8995b07536
# ╠═387730b4-bd06-492f-94e6-231bd68b3436
# ╠═290a7fe8-3b1e-423f-8b30-9bd8903d2e8f
# ╠═febe2b67-2a8f-4575-946d-30877bd5f2d4
# ╠═b4300e8a-8052-419b-98c8-0508ebee2393
# ╠═2c6ae74c-2da4-4867-ad8a-f4e835101d63
# ╠═c6a5d453-d610-4f65-847c-c878dd41726c
# ╠═fff17cf7-173d-4f64-94a9-4bf46acc882d
# ╠═eee55784-a641-445e-be75-0b19e2a94754
# ╠═7844e4f5-3c7d-4b4b-beee-970c998c67a6
# ╠═c81d0140-3f4e-4eb4-8a77-1f48c5e0ecbf
# ╠═7456455b-4f31-488f-990f-6ce534038e08
# ╠═b10c2ce4-d363-429c-a64c-ec29652137a5
# ╠═23b629a6-6866-416a-a768-9617ce6301db
# ╠═8788082d-f5d1-4385-8037-a0d360a841c7
# ╠═03a4fa85-08db-46d4-bb53-c0ccea90a211
# ╠═8f5e1413-227c-44fa-bb2d-3653cbc27e38
# ╠═8a7f7aab-878e-41b5-b9da-d06747df042e
# ╠═30602ff1-041b-4fca-bf8e-55ff57df9e37
# ╠═881011be-6434-416f-915b-3333e8dea32f
# ╠═5602b88a-07f8-438b-994c-65f11e17a0ba
# ╠═25be4255-0888-4ecd-a2fd-d66402c5cb50
# ╠═36e8a174-d526-45ee-b3c6-88d698ad5d5f
# ╠═7594147d-b3da-4e0d-896d-41baacb6d7be
# ╠═8e047be9-f0f0-4a75-91b5-c523f55f8c67
# ╠═1a9587aa-2356-48df-abcc-2ce874fa5d24
# ╠═785d500b-f8ea-446a-9952-2a5fd5d83d24
# ╠═e754826a-7411-4072-b0dc-a4bad7a15f98
# ╠═bc1edcbe-46eb-4531-9c5f-dee8d5dc2ff9
# ╠═47fb9e4c-df6a-4811-9980-99d595a34908
# ╠═bbd150de-ff9a-4127-a0dd-2f9762f92b07
# ╠═da5a0cbb-b033-46f1-a300-3954de138835
# ╠═73da8f53-c3af-43b0-9b23-60471f1e3587
# ╠═b0e67564-efe8-4fb2-bcf2-a711b770244e
# ╠═e5081280-d226-4834-8932-c89becd8313c
# ╠═c739a908-2353-4e7a-8fbd-f640dc8cabe0
# ╠═1d0383ad-54d6-4ff2-8555-def83bfff0e6
# ╠═cd1cf2f8-9f7f-4ed4-9cb7-1a6efee68ab4
# ╠═aba74ee0-0163-4e15-8b49-d8dcad4839f7
# ╠═0eda4142-1aaf-4e17-bd78-857e13e94acd
# ╠═50629194-98ed-4d45-86a2-95ac22daac29
# ╠═9c239fc6-275c-4d64-9fa2-6fd57295b757
# ╠═dede17f2-655d-4871-b6de-5a32804947dd
# ╠═e8c3f042-6058-4040-a67e-18563c04ee93
# ╠═665f4d03-3521-475a-a195-f861fd26bb69
# ╠═13775d58-ac61-431c-a6c9-447c1eec7942
# ╠═bc92f996-b626-4965-a286-d2c848eb1a21
# ╠═e2fde8e9-1f87-4ffe-8dae-2794664bfaa4
# ╠═23fd9fd9-a13c-4c56-a0b4-daec6f1d2cd8
# ╠═6f414d15-af0a-452a-98a1-dc0b7e54d617
# ╠═57ee8a12-24df-4598-935c-f5e259b504cb
# ╠═48bc085c-9ce6-4752-a5a9-a814f803f571
# ╠═5b9b2e5b-9deb-4ff6-a923-b15b8b08f0c9
# ╠═66ee9655-a006-4c0f-b1f1-5576efa8f896
# ╠═52d7975a-5346-447f-9aad-4ecd11b6460a
# ╠═e96efd6a-a666-4120-8480-9423e5d82ae1
# ╠═f0b4772d-a72b-44e0-a3a1-ba9ad4c4dfeb
