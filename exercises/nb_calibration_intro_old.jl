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

# ╔═╡ 8b6534d6-776b-4285-8498-a9b34051facc
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

# ╔═╡ 5b320989-3e0b-447b-bc9a-25fb221ce609
tspan = (0.0, 100.0)   # this will be the same for the three models

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
growth_mod_log = @reaction_network begin
    μ*W, ∅ --> W
    μ/Wf*W, W --> ∅
end

# ╔═╡ 1aa44f2b-6f33-437f-b9dd-89762d9f28ea
md"
### The measurement data
"

# ╔═╡ b2b433ed-0266-4bea-a7e8-32adba542d4c
md"
Assume that the measured grass yields (of a certain plant type) are the following:
"

# ╔═╡ 7c966a66-0091-4b81-9a7e-02ccd0d3db10
W_meas = [2.3, 4.5, 6.6, 7.6, 9.0, 9.1, 9.4]

# ╔═╡ 3edd2acc-a865-4675-afef-8868c68256f1
md"
They have been measured at the following corresponding time instances:
"

# ╔═╡ 877298e8-b61b-4c3a-ba2c-2827acdcfb50
t_meas = [0, 20, 29, 41, 50, 65, 72]

# ╔═╡ 68f608c1-607d-433d-8372-071a2bc541f8
md"
We will also assume that the error measurement is:
"

# ╔═╡ 34cc998b-2f96-4f30-b961-c3b584405982
W_sigma = 0.5

# ╔═╡ ef06cc43-510b-4ff9-b0b7-1c7fc267e9b1
md"
We can make a scatter plot of this data (including a title, a legend label, an X-axis label, X- and Y-axis limits) in the following way:
"

# ╔═╡ cb2bc6ee-4211-47e1-9956-5cf1b0c0671d
scatter(t_meas, W_meas, title="Grass growth data",
                        label="Yield",
                        xlabel="t",
                        xlims=(0, 80),
                        ylims=(0, 10))

# ╔═╡ d75246d4-e03b-4684-be7d-4bcfb61ed7ef
md"
### The objective function
"

# ╔═╡ dd3a32f1-bdb6-44a3-acbe-f4269725c9e4
md"
The objective function will now be implemented. We will implement very rudimentary objective functions and, hence, focus only on the essential concepts and not trying to seek for sophisticated and fancy implementations of this function.

For the logistic model, we will name our objective function as: `Jtheta_log`. The only input argument will be a vector of model parameters that need to be optimized.
"

# ╔═╡ 8a9115eb-4044-4cab-a7db-39b5dd86c70d
function Jtheta_log(thetas)           # function return to be minimized
    u₀ = [:W => thetas[1]]
    params = [:μ => thetas[2], :Wf => thetas[3]]
    oprob = ODEProblem(growth_mod_log, u₀, tspan, params, combinatoric_ratelaws=false)
    sol = solve(oprob, Tsit5(), saveat=t_meas)
    W_sol = sol[:W]
    J = (1/W_sigma^2)*sum(abs2, W_sol - W_meas)
    return J
end

# ╔═╡ 48c9f616-d298-40da-b917-225abd39b3d9
md"
Remarks on the function `Jtheta_log`:
- The input argument `thetas` will contain three values:
    - `thetas[1]` will be the initial condition for $W$
    - `thetas[2]` will be the parameter value for $\mu$
    - `thetas[3]` will be the parameter value for $W_f$
- Mind that vector indices always start with `1`.
- `u₀ = [:W => thetas[1]]` sets the initial condition.
- `params = [:μ => thetas[2], :Wf => thetas[3]]` sets the parameter values.
- You should already be familiar with the `ODEProblem` and `solve` commands.
- Be aware that the problem is being solved in the timepoints from the measurements, hence, in the `solve` command we use `saveat=t_meas`.
- In `W_sol = sol[:W]` the solution for the output variable $W$ is selected (in this example there is only one output variable) and assigned to `W_sol`.
- In `J = (1/sigma^2)*sum(abs2, W_sol - W_meas)` we compute the SSR weighted by $1/\sigma^{2}$ and assigned to `J` which will be returned by the function.
- In this example there is only one output variable, hence, only one term in the objective function $J(\theta)$. Thus, in this case it doesn't matter whether you weight the SSR by $1/\sigma^{2}$.
"

# ╔═╡ 35f158c1-858d-4e4d-ac3d-bf4807dad9a0
md"
In order to minimize the above objective function, we need to load the `Optim` package:
"

# ╔═╡ a6972aef-63ad-401c-acf5-6d59f9fc6698
md"
We are now ready to minimize the value of the objective function by optimizing the parameter values. This is done by calling the `optimize` function, providing the objective function `Jtheta_log`, initial values for the parameters to be optimized and (optionally) a minimization method (default: Nelder-Mead).

We will store the optimization results in `results_log`.
"

# ╔═╡ f34bb7ac-1ed8-4dd9-b0b9-49bd6e0e1d71
results_log = optimize(Jtheta_log, [2.0, 0.07, 10.0], NelderMead())

# ╔═╡ e55404ab-6762-4f39-bb42-9c195334a214
md"
You can obtain the optimized parameter values by calling `Optim.minimizer` on the results. We will store the optimized values in the variable `opt_log`.
"

# ╔═╡ 80e7f6b8-7592-48a3-8587-f1953d1bfcd8
opt_log = Optim.minimizer(results_log)

# ╔═╡ a1ca7d0e-639c-42d4-be09-5c61a2008f29
md"
You can obtain the minimum value of the objective function by calling `Optim.minimum` on the results.
"

# ╔═╡ 30399b9a-1d77-4140-9ad3-5eed636a5b99
Optim.minimum(results_log)

# ╔═╡ 72e065d4-7b1b-4f46-b373-935be8d801fc
md"
Now we can make a plot of $W$ simulated with the optimized initial value and parameter values.
"

# ╔═╡ 0b2ffd6f-01cd-4f11-9062-d38b3c13a5b1
md"
Setting up initial condition with optimized initial condition:
"

# ╔═╡ 590b1006-0e37-4668-9b4a-3588fab45696
u₀_opt_log = [:W => opt_log[1]]

# ╔═╡ 8ff28a2e-185d-4dca-ad3b-a0b1507646d6
md"
Setting up parameter values with optimized parameter values:
"

# ╔═╡ 14d24cbd-3259-41bb-9013-b4fe25a3be4c
params_opt_log = [:μ => opt_log[2], :Wf => opt_log[3]]

# ╔═╡ 6c134677-3ec1-4e0a-88b5-01341a096675
md"
Next, we create an ODEProblem and solve it:
"

# ╔═╡ 7a81f4a0-8f7c-4e05-9c3f-2438eab9b691
oprob_opt_log = ODEProblem(growth_mod_log, u₀_opt_log, tspan, params_opt_log)

# ╔═╡ ac80099b-d8f0-4eba-809d-d482bd354d35
osol_opt_log = solve(oprob_opt_log, Tsit5(), saveat=0.5)

# ╔═╡ b58f2c24-e0ea-48a8-b0b7-d0faf9642340
md"
Finally, we plot $W$ simulated with the optimized initial value and parameter values together with the measured data that was used to find the optimized values.
"

# ╔═╡ 55eba435-6ca4-4f4f-b08a-be700d5bda91
begin
	plot(osol_opt_log, label="Logistic growth", xlabel="t", xlims=(0, 80), ylims=(0, 10))
	scatter!(t_meas, W_meas, label="Yield")
end

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

# ╔═╡ febe2b67-2a8f-4575-946d-30877bd5f2d4
md"
Use the same measurement data (`W_meas`, `t_meas`, `W_sigma`) as before.
"

# ╔═╡ b4300e8a-8052-419b-98c8-0508ebee2393
md"
Implement the objective function for this case and name it `Jtheta_exp`.
"

# ╔═╡ 2c6ae74c-2da4-4867-ad8a-f4e835101d63
# Uncomment and complete the instruction
# function Jtheta_exp(thetas)
#     u₀ = missing
#     params = missing
#     oprob = missing
#     sol = missing
#     W_sol = missing
#     J = missing
#     return J
# end
function Jtheta_exp(thetas)
    u₀ = [:W => thetas[1]]
    params  = [:μ => thetas[2], :Wf => thetas[3]]
    oprob = ODEProblem(growth_exp, u₀, tspan, params, combinatoric_ratelaws=false)
    sol = solve(oprob, Tsit5(), saveat=t_meas)
    W_sol = sol[:W]
    J = (1/W_sigma^2)*sum(abs2, W_sol - W_meas)
    return J
end

# ╔═╡ eee55784-a641-445e-be75-0b19e2a94754
md"
Minimize the value of the objective function by optimizing the parameter values. Store the optimization results in `results_exp`.
"

# ╔═╡ 7844e4f5-3c7d-4b4b-beee-970c998c67a6
# results_exp = missing          # Uncomment and complete the instruction
results_exp = optimize(Jtheta_exp, [2.0, 0.02, 10.0], NelderMead())
# results_exp = optimize(Jtheta_exp, [2.0, 0.02, 1.0], NelderMead())

# ╔═╡ c81d0140-3f4e-4eb4-8a77-1f48c5e0ecbf
md"
Store the optimized values in the variable `opt_exp`:
"

# ╔═╡ 7456455b-4f31-488f-990f-6ce534038e08
# opt_exp = missing              # Uncomment and complete the instruction
opt_exp = Optim.minimizer(results_exp)

# ╔═╡ b10c2ce4-d363-429c-a64c-ec29652137a5
md"
Determine the minimum value of the objective function:
"

# ╔═╡ 23b629a6-6866-416a-a768-9617ce6301db
# missing                        # Uncomment and complete the instruction
Optim.minimum(results_exp)

# ╔═╡ 8a7f7aab-878e-41b5-b9da-d06747df042e
md"
Set up initial condition with optimized initial condition:
"

# ╔═╡ 30602ff1-041b-4fca-bf8e-55ff57df9e37
# missing                        # Uncomment and complete the instruction
u₀_opt_exp = [:W => opt_exp[1]]

# ╔═╡ 881011be-6434-416f-915b-3333e8dea32f
md"
Set up parameter values with optimized parameter values:
"

# ╔═╡ 5602b88a-07f8-438b-994c-65f11e17a0ba
# params_opt_exp = missing       # Uncomment and complete the instruction
params_opt_exp = [:μ => opt_exp[2], :Wf => opt_exp[3]]

# ╔═╡ 25be4255-0888-4ecd-a2fd-d66402c5cb50
md"
Create an ODEProblem and solve it:
"

# ╔═╡ 36e8a174-d526-45ee-b3c6-88d698ad5d5f
# oprob_opt_exp = missing        # Uncomment and complete the instruction
oprob_opt_exp = ODEProblem(growth_exp, u₀_opt_exp, tspan, params_opt_exp)

# ╔═╡ 7594147d-b3da-4e0d-896d-41baacb6d7be
# osol_opt_exp = missing
osol_opt_exp = solve(oprob_opt_exp, Tsit5(), saveat=0.5)

# ╔═╡ 8e047be9-f0f0-4a75-91b5-c523f55f8c67
md"
Finally, we plot $W$ simulated with the optimized initial value and parameter values together with the measured data that was used to find the optimized values.
"

# ╔═╡ 1a9587aa-2356-48df-abcc-2ce874fa5d24
begin
	plot(osol_opt_exp, label="Exponential growth",
		xlabel="t", xlims=(0, 80), ylims=(0, 10))
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

# ╔═╡ b0e67564-efe8-4fb2-bcf2-a711b770244e
md"
Use the same measurement data (`W_meas`, `t_meas`, `W_sigma`) as before.
"

# ╔═╡ e5081280-d226-4834-8932-c89becd8313c
md"
Implement the objective function for this case and name it `Jtheta_gom`.
"

# ╔═╡ c739a908-2353-4e7a-8fbd-f640dc8cabe0
# Uncomment and complete the instruction
# function Jtheta_gom(thetas)
#     u₀ = missing
#     params = missing
#     oprob = missing
#     sol = missing
#     W_sol = missing
#     J = missing
#     return J
# end
function Jtheta_gom(thetas)
    u₀ = [:W => thetas[1]]
    params  = [:μ => thetas[2], :D => thetas[3]]
    oprob = ODEProblem(growth_gom, u₀, tspan, params, combinatoric_ratelaws=false)
    sol = solve(oprob, Tsit5(), saveat=t_meas)
    W_sol = sol[:W]
    J = (1/W_sigma^2)*sum(abs2, W_sol - W_meas)
    return J
end

# ╔═╡ 1d0383ad-54d6-4ff2-8555-def83bfff0e6
md"
Minimize the value of the objective function by optimizing the parameter values. Store the optimization results in `results_gom`.
"

# ╔═╡ cd1cf2f8-9f7f-4ed4-9cb7-1a6efee68ab4
# results_gom = missing         # Uncomment and complete the instruction
results_gom = optimize(Jtheta_gom, [2.0, 0.09, 0.04], NelderMead())

# ╔═╡ aba74ee0-0163-4e15-8b49-d8dcad4839f7
md"
Store the optimized values in the variable `opt_gom`:
"

# ╔═╡ 0eda4142-1aaf-4e17-bd78-857e13e94acd
# opt_gom = missing             # Uncomment and complete the instruction
opt_gom = Optim.minimizer(results_gom)

# ╔═╡ 50629194-98ed-4d45-86a2-95ac22daac29
md"
Determine the minimum value of the objective function:
"

# ╔═╡ 9c239fc6-275c-4d64-9fa2-6fd57295b757
# missing                       # Uncomment and complete the instruction
Optim.minimum(results_gom)

# ╔═╡ 1904b8a6-5ff6-49d6-9f75-c2f524de181a
md"
Set up initial condition with optimized initial condition:
"

# ╔═╡ e8c3f042-6058-4040-a67e-18563c04ee93
# u₀_opt_gom = missing          # Uncomment and complete the instruction
u₀_opt_gom = [:W => opt_gom[1]]

# ╔═╡ 6f414d15-af0a-452a-98a1-dc0b7e54d617
md"
Set up parameter values with optimized parameter values:
"

# ╔═╡ 57ee8a12-24df-4598-935c-f5e259b504cb
# params_opt_gom = missing     # Uncomment and complete the instruction
params_opt_gom = [:μ => opt_gom[2], :D => opt_gom[3]]

# ╔═╡ 48bc085c-9ce6-4752-a5a9-a814f803f571
md"
Create an ODEProblem and solve it:
"

# ╔═╡ 5b9b2e5b-9deb-4ff6-a923-b15b8b08f0c9
# oprob_opt_gom = missing      # Uncomment and complete the instruction
oprob_opt_gom = ODEProblem(growth_gom, u₀_opt_gom, tspan, params_opt_gom)

# ╔═╡ 66ee9655-a006-4c0f-b1f1-5576efa8f896
# osol_opt_gom = missing
osol_opt_gom = solve(oprob_opt_gom, Tsit5(), saveat=0.5)

# ╔═╡ 52d7975a-5346-447f-9aad-4ecd11b6460a
md"
Finally, we plot $W$ simulated with the optimized initial value and parameter values together with the measured data that was used to find the optimized values.
"

# ╔═╡ e96efd6a-a666-4120-8480-9423e5d82ae1
begin
	plot(osol_opt_gom, label="Gompertz growth", xlabel="t",
		xlims=(0, 80), ylims=(0, 10))
	scatter!(t_meas, W_meas, label="Yield")
end

# ╔═╡ f0b4772d-a72b-44e0-a3a1-ba9ad4c4dfeb
md"
Which grass growth model fits best these data? How can you prove this numerically?
- Answer: missing

Conduct the calibration again for the exponential growth model, but now set
the initial value for $W_f$ equal to $1.0$. How can you explain this result?
- Answer: missing
"

# ╔═╡ Cell order:
# ╠═a09f814a-0c6a-11ef-0e79-a50b01287d63
# ╠═7f521435-63ac-4178-a4aa-93d9c45fe820
# ╠═f8a92690-990b-4341-89e1-322adbcb8d1b
# ╠═015050b3-3339-4b1a-ad7d-c358cce73675
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
# ╠═1aa44f2b-6f33-437f-b9dd-89762d9f28ea
# ╠═b2b433ed-0266-4bea-a7e8-32adba542d4c
# ╠═7c966a66-0091-4b81-9a7e-02ccd0d3db10
# ╠═3edd2acc-a865-4675-afef-8868c68256f1
# ╠═877298e8-b61b-4c3a-ba2c-2827acdcfb50
# ╠═68f608c1-607d-433d-8372-071a2bc541f8
# ╠═34cc998b-2f96-4f30-b961-c3b584405982
# ╠═ef06cc43-510b-4ff9-b0b7-1c7fc267e9b1
# ╠═cb2bc6ee-4211-47e1-9956-5cf1b0c0671d
# ╠═d75246d4-e03b-4684-be7d-4bcfb61ed7ef
# ╠═dd3a32f1-bdb6-44a3-acbe-f4269725c9e4
# ╠═8a9115eb-4044-4cab-a7db-39b5dd86c70d
# ╠═48c9f616-d298-40da-b917-225abd39b3d9
# ╠═35f158c1-858d-4e4d-ac3d-bf4807dad9a0
# ╠═8b6534d6-776b-4285-8498-a9b34051facc
# ╠═a6972aef-63ad-401c-acf5-6d59f9fc6698
# ╠═f34bb7ac-1ed8-4dd9-b0b9-49bd6e0e1d71
# ╠═e55404ab-6762-4f39-bb42-9c195334a214
# ╠═80e7f6b8-7592-48a3-8587-f1953d1bfcd8
# ╠═a1ca7d0e-639c-42d4-be09-5c61a2008f29
# ╠═30399b9a-1d77-4140-9ad3-5eed636a5b99
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
# ╠═137bde23-76f2-4ebf-8bc2-ea8640001436
# ╠═4aa71200-006b-4a15-ae75-67e36aa81522
# ╠═cdab3079-04b0-4a44-b770-468c20e321e4
# ╠═cf1a144e-09e9-42a3-b2a3-b8676a200a39
# ╠═febe2b67-2a8f-4575-946d-30877bd5f2d4
# ╠═b4300e8a-8052-419b-98c8-0508ebee2393
# ╠═2c6ae74c-2da4-4867-ad8a-f4e835101d63
# ╠═eee55784-a641-445e-be75-0b19e2a94754
# ╠═7844e4f5-3c7d-4b4b-beee-970c998c67a6
# ╠═c81d0140-3f4e-4eb4-8a77-1f48c5e0ecbf
# ╠═7456455b-4f31-488f-990f-6ce534038e08
# ╠═b10c2ce4-d363-429c-a64c-ec29652137a5
# ╠═23b629a6-6866-416a-a768-9617ce6301db
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
# ╠═b0e67564-efe8-4fb2-bcf2-a711b770244e
# ╠═e5081280-d226-4834-8932-c89becd8313c
# ╠═c739a908-2353-4e7a-8fbd-f640dc8cabe0
# ╠═1d0383ad-54d6-4ff2-8555-def83bfff0e6
# ╠═cd1cf2f8-9f7f-4ed4-9cb7-1a6efee68ab4
# ╠═aba74ee0-0163-4e15-8b49-d8dcad4839f7
# ╠═0eda4142-1aaf-4e17-bd78-857e13e94acd
# ╠═50629194-98ed-4d45-86a2-95ac22daac29
# ╠═9c239fc6-275c-4d64-9fa2-6fd57295b757
# ╠═1904b8a6-5ff6-49d6-9f75-c2f524de181a
# ╠═e8c3f042-6058-4040-a67e-18563c04ee93
# ╠═6f414d15-af0a-452a-98a1-dc0b7e54d617
# ╠═57ee8a12-24df-4598-935c-f5e259b504cb
# ╠═48bc085c-9ce6-4752-a5a9-a814f803f571
# ╠═5b9b2e5b-9deb-4ff6-a923-b15b8b08f0c9
# ╠═66ee9655-a006-4c0f-b1f1-5576efa8f896
# ╠═52d7975a-5346-447f-9aad-4ecd11b6460a
# ╠═e96efd6a-a666-4120-8480-9423e5d82ae1
# ╠═f0b4772d-a72b-44e0-a3a1-ba9ad4c4dfeb
