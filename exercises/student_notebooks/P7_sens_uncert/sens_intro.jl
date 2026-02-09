### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ e1e7bc8e-7264-4cbc-98d2-aa73679fa2df
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 5f4fea06-0632-11ef-102e-21f5606d2056
using Markdown

# ╔═╡ 489b5399-fe4c-481a-834f-0101bbe28cea
using InteractiveUtils

# ╔═╡ 22489bd4-ab64-4bf9-ad03-5372ea273935
using Catalyst

# ╔═╡ cdcebbb1-40e0-457f-a6ec-b769f6b1f2e9
using DifferentialEquations, Plots

# ╔═╡ 65885bbe-7b73-4efd-a8e7-984a34b80548
using ForwardDiff

# ╔═╡ 427b509f-d08e-4d93-99ba-a79f9c244b28
md"""
# Introduction to sensitivity analysis
"""

# ╔═╡ 49e7085c-3691-4182-9630-68dc9371ad18
md"""
## Goal of this practicum
"""

# ╔═╡ c41894dd-0f6c-483f-b19f-dbf8148f776f
md"""
Sensitivity functions indicate how sensitive the model output is to a change in parameter values. When a model output is very sensitive to a certain parameter, a small change in the value of this parameter will have a large influence on the value of the model output. Sensitivity functions thus provide important information about the model and are implicitely used to estimate parameters and explicitely in the context of optimal experimental design.
"""

# ╔═╡ 7b58f2f9-cb93-4a20-bf90-62a9006b57d6
md"""
The sensitivity function that measures how sensitive output $y_j$ is to changes in parameter $\theta_i$ is given by the partial derivative

$$\cfrac{\partial \hat{y}_i(\theta)}{\partial \theta_j}\tag{1}$$
"""

# ╔═╡ abd70f04-1bf8-454f-869a-0d8082d68453
md"""
Try to understand why the above expression $(1)$ does indeed give us the information we were promised in the first paragraph. How will we be able to see from the value determined by the above expression $(1)$ that whether or not output $y_j$ is sensitive to a change in $\theta_i$?
"""

# ╔═╡ e76aedec-69f2-4301-b015-6960e4503c42
md"""
Answer: missing
"""

# ╔═╡ a6a97279-0b64-45c3-8312-c22b1a8425d0
md"""
Sometimes expression $(1)$ can be evaluated analytically. Usually, however, we will have to approximate the partial derivative numerically. Expression $(1)$ can be made more specific:

$$\cfrac{\partial \hat{y}_i(\theta)}{\partial \theta_j} \approx \cfrac{\hat{y}_i(\theta_j+\Delta\theta_j)-\hat{y}_i(\theta_j)}{\Delta\theta_j}$$
"""

# ╔═╡ 25d3ef62-382c-455c-92bb-fadfc650c5a6
md"""
Thus, to calculate the sensitivity function numerically, the model is evaluated for the parameter values $\theta_i$ and $\theta_i+\Delta\theta_i$ and the difference between these evaluations is taken.
"""

# ╔═╡ 2a544bef-d7fe-4300-bc59-69e6f7429304
md"""
Since quantity $(1)$ is dependent on the units, a normalized variant is often used:

$$\cfrac{\partial \hat{y}_i(\theta)}{\partial \theta_j} \cdot \cfrac{\theta_j}{\hat{y}_i} \tag{2}$$

The interpretation of $(2)$ is how much the output changes per cent if the parameter is increased by one per cent. It assumes positive model outputs and parameters, which is often the case for biochemical models.

Using the normalized variant allows you to compare all possible sensitivity functions with each other.
"""

# ╔═╡ 7107a3a9-15ef-488e-8709-52a6444d0e1c
md"""
We now calculate and interpret sensitivity functions for some given models. To illustrate the concepts, we first consider three simple models modelling the growth of grass.
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

with output $W$ the grass yield, and $W_f$, $\mu$ and $D$ parameters. The table below show some typical values for the parameters:

|             | $\mu$      | $W_f$       | $D$          |
|:----------- |:----------:|:-----------:|:------------:|
| Logistic    |  0.07      | 10.0        |              |
| Exponential |  0.02      | 10.0        |              |
| Gompertz    |  0.09      |             | 0.04         |

We will use an initial condition of $W_0 = 2.0$ for each and a simulation time of $100$ days.
"""

# ╔═╡ 8d55bc42-1f8f-4bc1-b9fa-e8152fc125ce
md"""
We will illustrate how to compute the local sensitivity functions for the logistic model. The same will be left as exercises below for the exponential and gompertz models.

**Important:**
- We will use consequently `_log`, `_exp` and `_gom` appended to relevant variables names in order to indicate their model origin **and** to prevent cell-disabling that occurs when using the same variables names in these Notebooks.
"""

# ╔═╡ 676119ce-f4fe-41f9-8121-b2c21f0dd28c
md"""
### Modelling logistic growth

We will start by modelling our system and simulating using the aforementioned parameters values, initial condition and timespan.
"""

# ╔═╡ 5492210b-6ee6-4baf-b37e-d21358cdeb60
md"""
Implementation of the system:
"""

# ╔═╡ 18eece66-46fa-4458-aaed-f4c8fa002c20
growth_log = @reaction_network begin
	@species W(t)=2.0            # default initial condition
	@parameters μ=0.07 Wf=10.0   # default parameter values
    #μ*W, 0 --> W
    #μ/Wf*W, W --> 0
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

# ╔═╡ 0b3d35bb-c5b0-44b7-94b3-06fa571d339e
tspan = (0.0, 100.0)   # this will be the same for the three models

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
oprob_log = ODEProblem(growth_log, u0_log, tspan, params_log)
# Also possible:
# oprob_log = ODEProblem(growth_mod_log, [], tspan, [])

# ╔═╡ ac235d86-1d93-4944-aa89-1b4fd38f0e6e
osol_log = solve(oprob_log, Tsit5(), saveat=0.5)

# ╔═╡ 3e13efa1-9bc6-456f-8e62-ecd3165e2a65
plot(osol_log)

# ╔═╡ 79fe5411-fb78-490e-a9de-1868d1261aa6
md"""
### Local Sensitivity Analysis (LSA)
"""

# ╔═╡ c57938c8-566a-4502-a2e1-69dc77291500
md"""
In order to compute the local sensitivity functions, we will need to load the `ForwardDiff` package:
"""

# ╔═╡ 4aa73da9-a394-4ca3-a839-d076eb3c3d7f
md"""
We need to write a solution function with as argument a vector of the parameters (that you want the sensitivity on), and that returns the solution (time vector and outputs).
"""

# ╔═╡ bc201b61-0f58-49d2-a20b-4f18c42fcc96
function growth_sim_log(params)
	μ, Wf = params
    u0_log = [:W => 2.0]
    tspan = (0.0, 100.0)
    oprob_log = ODEProblem(growth_log, u0_log, tspan, [:μ=>μ, :Wf=>Wf])
	osol_log = solve(oprob_log, Tsit5(), saveat=0.5)
	return osol_log
end

# ╔═╡ 11f15a0e-2958-4075-928c-5a24bcc00c69
md"""
Next we will need to make a function based on the solution function that returns a single output.
"""

# ╔═╡ d0f3197f-3094-4f6e-b33c-e5b74e0947f6
growth_sim_W_log(params) = growth_sim_log(params)[:W]

# ╔═╡ ba4eaa90-b97a-4b71-b3e1-968340c6def5
md"""
Now make a time vector that is the same as the time vector from the solution.
"""

# ╔═╡ 71f897a6-1e91-412b-ad58-c5f1e7cd1adb
t_vals_log = 0:0.5:100.0
# Alternative:
# t_vals_log = tspan[1]:0.5:tspan[2]

# ╔═╡ b3b97b57-a4f7-4cb6-80c8-6ad893ade75d
md"""
Compute the single output with the given parameter values.
"""

# ╔═╡ 99864b21-3db0-4d60-a84e-8e96db4de4ae
W_log = growth_sim_W_log([μ_log, Wf_log])

# ╔═╡ 344eb3a7-8614-4349-a536-62d9acef6bda
md"""
Use the function `ForwardDiff.jacobian` to compute the sensitivities. This function takes two arguments: the solution function and a vector with the parameter values.
"""

# ╔═╡ 670ccfa4-24fa-4167-a46a-2a48dc19538b
sens_W_log = ForwardDiff.jacobian(growth_sim_W_log, [μ_log, Wf_log])

# ╔═╡ c262d3ed-9f62-45ad-a7ee-ec78ab16f35c
md"""
To get the sensitivities of $W$ on $\mu$, and of $W$ on $W_f$, you need to use indexing with `sens_W_log`:
- `sens_W_log[:,1]` gives the (absolute) sensitivity of $W$ on $\mu$.
- `sens_W_log[:,2]` gives the (absolute) sensitivity of $W$ on $W_f$.
"""

# ╔═╡ 6600d28e-2522-4069-a5e0-643be43f6117
sens_W_on_μ_log = sens_W_log[:,1]    # sensitivity of W on μ

# ╔═╡ d9e3a4ac-c138-4b93-a19c-737166a0f0ea
sens_W_on_Wf_log = sens_W_log[:,2]   # sensitivity of W on Wf

# ╔═╡ f07d8205-e172-45fa-b955-792bd95f3023
md"""
We now calculate the normalized sensitivities. For that we need to multiply by the parameter value and divide by the ouput. Beware that all element wise operations need a dot in front of the operator, e.g. as in `.*` and `./`.
"""

# ╔═╡ 1c142797-ada6-4d77-a22f-b981ffd38956
sens_W_on_μ_rel_log = sens_W_on_μ_log .* μ_log ./ W_log

# ╔═╡ cb397e0b-56f9-420f-a86e-bfab35286b44
sens_W_on_Wf_rel_log = sens_W_on_Wf_log .* Wf_log ./ W_log

# ╔═╡ 9e6946a4-f207-4bfb-9e19-aa7b77c2a05b
md"""
We are now ready to plot the two sensitivity functions. We provide the time vector (first argument) and a vector of the two sensitivity functions (second argument). Additionally, you can provide a title, (legend) labels and a x- and/or y-label.
"""

# ╔═╡ f58ce914-f366-4087-a7d1-8cfe69ac623b
plot(t_vals_log, [sens_W_on_μ_rel_log, sens_W_on_Wf_rel_log], title="Normalized sensitivities", label=["W on μ" "W on Wf"], xlabel="Time (day)")

# ╔═╡ 05972c7f-f64f-4865-b0a7-f33029d0a6fa
md"""
Notice that in the `label` option there is no comma separating the labels.
"""

# ╔═╡ 88aa13f6-e16c-40cd-ac10-3ede8b9fb429
# ╠═╡ disabled = true
#=╠═╡
md"""
#### Alternative

For each plot command we need to provide the time vector (first argument) and the column vector with the local sensitivity (second argument). Additionally, you can provide a title, a (legend) label and a x- and/or y-label.

Beware that if you want to execute multiple commands (here: `plot`) in a single cell, you need to put them in a `begin`-and-`end` block. Also if you want to visualize subsequent graphs in the same plot, the forthcoming plot command names should be followed by a `!`, as in `plot!(...)`.
"""
  ╠═╡ =#

# ╔═╡ dc0507ee-0f67-4556-bdc3-1294177c86b7
# ╠═╡ disabled = true
#=╠═╡
begin
plot(t_vals_log, sens_W_on_μ_rel_log, title="Normalized sensitivities", label="W on μ", xlabel="Time (day)")
plot!(t_vals_log, sens_W_on_Wf_rel_log, label="W on Wf")
end
  ╠═╡ =#

# ╔═╡ 04a1d1ad-53a9-4146-be09-65aa422e3730
md"""
Conclusions:
- From the sensitivity plot of $W$ on $\mu$ it can be seen that $W$ is most sensitive to $\mu$ in the time region $[0, 30]$. The latter corresponds to the region where the yield rate is largest (i.e., when the growth is largest). This makes sense because when looking at the differential equation, $\mu$ is approximately the growth rate for relatively small $W$ values.
- From the sensitivity plot of $W$ on $W_f$ it can be seen that $W$ is most sensitive to $W_f$ in the region where time values are large (cf. operating point). The latter corresponds to the region where the yield rate stagnates (i.e., when the yield reaches a steady value). This makes sense because when looking at the differential equation, $W_f$ is the steady state value.
"""

# ╔═╡ 79f0f1dd-850f-4dfd-b895-ff41a2d8adb8
md"""
## Exercises
"""

# ╔═╡ fa2270d0-2548-409c-a23f-4369d8bce8ec
md"""
### Exercise 1 - Sensitivity analysis of the exponential growth model
"""

# ╔═╡ dca50b37-b06c-4efe-881e-cf966ebc8fd7
md"""
Create a *reaction network object* for the exponential growth model. Name it `growth_exp`.
"""

# ╔═╡ 2836f231-3b90-4829-9012-3ed9a09239ff
# Uncomment and complete the instruction
# growth_exp = @reaction_network begin
#     @species missing
#     @parameters missing
#     missing
# end
growth_exp = @reaction_network begin
	@species W(t)=2.0
	@parameters μ=0.02 Wf=10.0
	(μ*Wf, μ), 0 <--> W
	# Alternative:
    # μ*Wf, 0 --> W
    # μ, W --> 0
end

# ╔═╡ e8db84be-31a0-415d-a194-064c4c87a293
md"""
Convert the system to a symbolic differential equation model (name it: `osys_exp`) and verify, by analyzing the differential equation, that your model has been correctly implemented.
"""

# ╔═╡ 78d175b4-a9b2-49b5-bbc8-eb348799985b
# osys_exp = missing           # Uncomment and complete the instruction
osys_exp = convert(ODESystem, growth_exp)

# ╔═╡ 63952c54-3304-4500-9536-b375c5c8f280
md"""
Initialize a vector `u0_exp` with the initial condition:
"""

# ╔═╡ 7dc013f7-7b35-42e4-aa77-0dca3755f389
# u0_exp = missing           # Uncomment and complete the instruction
u0_exp = [:W => 2.0]

# ╔═╡ 86785c55-b567-449c-ae76-bc15a16223bc
md"""
We will use the same timespan as before, so no need to redefine it.
"""

# ╔═╡ 05e96232-ce95-452b-98e4-e79817d45ae2
md"""
For the sake of clarity, we will use the variables `μ_exp` and `Wf_exp` to store the parameter values.
"""

# ╔═╡ 9a705b3f-3c45-4896-9b89-e9c37d47fd89
μ_exp = 0.02

# ╔═╡ 55a13ebc-09e7-433c-83a8-9b6805e2fbec
Wf_exp = 10.0

# ╔═╡ f1293675-b4d0-4827-8714-65a593ecc00e
md"""
Initialize a vector `params_exp` with the parameter values:
"""

# ╔═╡ 7005d88f-9e28-47fe-9a32-dec023903dc3
# params_exp = missing        # Uncomment and complete the instruction
params_exp = [:μ => μ_exp, :Wf => Wf_exp]

# ╔═╡ c8853554-26c5-491c-a536-4c07e8c6a986
md"""
Create the ODE problem and store it in `oprob_exp`:
"""

# ╔═╡ 025c1154-b2ae-4e1c-af1b-277b24d648a4
# oprob_exp = missing         # Uncomment and complete the instruction
oprob_exp = ODEProblem(growth_exp, u0_exp, tspan, params_exp)

# ╔═╡ f58825d2-55fc-43f1-b164-9555bf9f5b84
md"""
Solve the ODE problem. Use `Tsit5()` and `saveat=0.5`. Store the solution in `osol_exp`:
"""

# ╔═╡ 8044dd5f-9996-4062-a54a-8fb0eee30b46
# osol_exp = missing          # Uncomment and complete the instruction
osol_exp = solve(oprob_exp, Tsit5(), saveat=0.5)

# ╔═╡ a1752ad4-2c11-416c-8321-9da058f9aaea
md"""
Plot the result:
"""

# ╔═╡ 0f58fd3d-1495-46c7-9d92-5539395cf12e
# missing                     # Uncomment and complete the instruction
plot(osol_exp)

# ╔═╡ 6321774c-5bbe-42a7-bdad-0168c891b5ce
md"""
Write a solution function with as argument a vector of the parameters (that you want the sensitivity on), and that returns the outputs.
"""

# ╔═╡ eb1e5f97-3e09-43f0-b6ce-cc6cbc42ec2f
# Uncomment and complete the instruction
# function growth_sim_exp(params)
#     missing
#     ...
# end
function growth_sim_exp(params)
	μ, Wf = params
    u0_exp = [:W => 2.0]
    tspan = (0.0, 100.0)
    oprob_exp = ODEProblem(growth_exp, u0_exp, tspan, [:μ=>μ, :Wf=>Wf])
	osol_exp = solve(oprob_exp, Tsit5(), saveat=0.5)
	return osol_exp
end

# ╔═╡ b31b4d24-ae30-4cff-b6b7-348200feaa5d
md"""
Make a function based on the solution function that returns a single output.
"""

# ╔═╡ f2fa716f-d7b5-4c95-a8c6-6b16fbfe1499
# growth_sim_W_exp(params) = missing   # Uncomment and complete the instruction
growth_sim_W_exp(params) = growth_sim_exp(params)[:W]

# ╔═╡ 0b413b98-a800-44e9-a917-ba61467bd613
md"""
Make the time vector.
"""

# ╔═╡ 6be6fea7-561b-4eec-b249-fa4522a5b039
# t_vals = missing             # Uncomment and complete the instruction
t_vals_exp = 0:0.5:100.0

# ╔═╡ 98070d41-e8c5-49e6-9eea-04c718cbff65
md"""
Compute the output for the given parameter values.
"""

# ╔═╡ cef751aa-f9c8-46da-9436-9cc2e5d7515c
# W_exp = missing              # Uncomment and complete the instruction
W_exp = growth_sim_W_exp([μ_exp, Wf_exp])

# ╔═╡ 9e7bb7ec-8fe3-422e-920f-13c2ef055feb
md"""
Using `ForwardDiff.jacobian` to compute the sensitivities for the single ouput(s).
"""

# ╔═╡ b13b5e77-ca60-465b-9f63-be9e5da0482f
# sens_W_exp = missing           # Uncomment and complete the instruction
sens_W_exp = ForwardDiff.jacobian(growth_sim_W_exp, [μ_exp, Wf_exp])

# ╔═╡ 735f5191-259e-440f-a379-20a2a70ec72c
md"""
Extract the (absolute) sensitivities of the outputs on the different parameters.
"""

# ╔═╡ 689556c9-a690-4e9e-9052-7ac66e999d4d
# sens_W_on_μ_exp = missing      # Uncomment and complete the instruction
sens_W_on_μ_exp = sens_W_exp[:,1]

# ╔═╡ b85f0d09-a102-435b-a898-149ba4b29caa
# sens_W_on_Wf_exp = missing     # Uncomment and complete the instruction
sens_W_on_Wf_exp = sens_W_exp[:,2]

# ╔═╡ d067f814-1dd5-4e59-8cfd-3bc0b2d97612
md"""
Compute the normalized sensitivities.
"""

# ╔═╡ e01657f3-33ac-4eb8-b39e-e2441a165a8d
# sens_W_on_μ_rel_exp = missing   # Uncomment and complete the instruction
sens_W_on_μ_rel_exp = sens_W_on_μ_exp .* μ_exp ./ W_exp

# ╔═╡ 7231e54e-b293-48fc-b1f0-273f45f51539
# sens_W_on_Wf_rel_exp = missing   # Uncomment and complete the instruction
sens_W_on_Wf_rel_exp = sens_W_on_Wf_exp .* Wf_exp ./ W_exp

# ╔═╡ 0524a976-9d4d-4ef9-b088-bdd9102415e6
md"""
Plot both normalized sensitivity functions (with appropriate title and labels):
"""

# ╔═╡ 33031b53-6f2d-4256-a0aa-0ba97beaa7af
# missing     # Uncomment and complete the instruction
plot(t_vals_exp, [sens_W_on_μ_rel_exp, sens_W_on_Wf_rel_exp], title="Normalized sensitivities", label=["W on μ" "W on Wf"], xlabel="Time (day)")

# ╔═╡ 22b58c97-1aa2-49c4-8a0a-488c63014a90
md"
Draw your conclusions:

- missing
- missing
"

# ╔═╡ 41abf8c4-e67e-4f66-a38a-7204a878d98d
md"""
### Exercise 2 - Sensitivity analysis of the Gompertz growth model
"""

# ╔═╡ f311e943-c96b-4af3-8757-c662fd302a88
md"""
Create a *reaction network object* for the Gompertz growth model. Name it `growth_gom`.
"""

# ╔═╡ 1ae76036-42ef-46e9-88bf-d66d4267addc
# Uncomment and complete the instruction
# growth_gom = @reaction_network begin
#     @species missing
#     @parameters missing
#     missing
# end
growth_gom = @reaction_network begin
	@species W(t)=2.0
	@parameters μ=0.09 D=0.04
	μ-D*log(W), W --> 2*W
	# Alternative:
    # μ, W --> 2*W
    # D*log(W), W --> 0
end

# ╔═╡ e190ceca-30b1-49e2-baf5-ceb62929f4c0
md"""
Convert the system to a symbolic differential equation model (name it: `osys_gom`) and verify, by analyzing the differential equation, that your model has been correctly implemented.
"""

# ╔═╡ acaf5b91-accc-4dfa-9371-6495a47c8736
# osys_gom = missing                  # Uncomment and complete the instruction
osys_gom = convert(ODESystem, growth_gom)

# ╔═╡ ed1150f1-9c75-4867-b7f5-535d605810f4
md"""
Initialize a vector `u0_gom` with the initial condition:
"""

# ╔═╡ d35ef155-61aa-4b19-872f-d6e621f96572
# u0_gom = missing                    # Uncomment and complete the instruction
u0_gom = [:W => 2.0]

# ╔═╡ 3f546d58-15ff-442d-a4e8-7c0783ed2fe0
md"""
We will use the same timespan as before, so no need to redefine it.
"""

# ╔═╡ b778f786-87bc-4103-a792-2ae5cf004d50
md"""
For the sake of clarity, we will use the variables `μ_gom` and `D_gom` to store the parameter values.
"""

# ╔═╡ edd76a4e-e541-4a83-923f-6f593d8fadad
μ_gom = 0.09

# ╔═╡ 05f1ae6c-1229-4924-9c41-9554d6ca4261
D_gom = 0.04

# ╔═╡ 70df05e7-4527-4087-a762-942a33e89f74
md"""
Initialize a vector `params_gom` with the parameter values:
"""

# ╔═╡ bf9291d5-6fc5-4359-874e-167237a22147
# params_gom = missing                # Uncomment and complete the instruction
params_gom = [:μ => μ_gom, :D => D_gom]

# ╔═╡ 9e3158c5-d4b2-456c-b920-8431e2047770
md"""
Create the ODE problem and store it in `oprob_gom`:
"""

# ╔═╡ 3aab9073-cf4d-4ee5-9e8a-e3db586d9f68
# oprob_gom = missing                 # Uncomment and complete the instruction
oprob_gom = ODEProblem(growth_gom, u0_gom, tspan, params_gom)

# ╔═╡ 30c35371-d106-4673-8569-d07b39edbdfe
md"""
Solve the ODE problem. Use `Tsit5()` and `saveat=0.5`. Store the solution in `osol_gom`:
"""

# ╔═╡ c01fe052-3701-4e22-9a64-f681d52a445a
# osol_gom = missing                  # Uncomment and complete the instruction
osol_gom = solve(oprob_gom, Tsit5(), saveat=0.5)

# ╔═╡ 4f959d21-954b-4116-a498-bafc1aada47a
md"""
Plot the result:
"""

# ╔═╡ 5f427df2-1b7e-46c7-adb9-66b05deac6e8
# missing                             # Uncomment and complete the instruction
plot(osol_gom)

# ╔═╡ 958854fc-8e23-4364-88c5-1ec060b86a16
md"""
Write a solution function with as argument a vector of the parameters (that you want the sensitivity on), and that returns the outputs.
"""

# ╔═╡ be21e269-16c4-42b0-8930-2d009bd91161
# Uncomment and complete the instruction
# function growth_sim_gom(params)
#     missing
#     ...
# end
function growth_sim_gom(params)
	μ, D = params
    u0_gom = [:W => 2.0]
    tspan = (0.0, 100.0)
    oprob_gom = ODEProblem(growth_gom, u0_gom, tspan, [:μ=>μ, :D=>D])
	osol_gom = solve(oprob_gom, Tsit5(), saveat=0.5)
	return osol_gom
end

# ╔═╡ e57f5916-b753-429b-8434-8a85a67bd1fb
md"""
Make a function based on the solution function that returns a single output.
"""

# ╔═╡ 69f96b57-7c8d-4e64-9e6a-9804d001b0b1
# growth_sim_W_gom(params) = missing  # Uncomment and complete the instruction
growth_sim_W_gom(params) = growth_sim_gom(params)[:W]

# ╔═╡ f7f6c56b-4910-4641-a462-eacfa4b4d034
md"""
Make the time vector.
"""

# ╔═╡ 84c1bb38-23c3-4f70-be67-7da2300a737b
# t_vals_gom = missing      # Uncomment and complete the instruction
t_vals_gom = 0:0.5:100.0

# ╔═╡ 09eeab39-79a1-420a-8298-72b7002ec168
md"""
Compute the output for the given parameter values.
"""

# ╔═╡ 9e15601f-4ad2-41be-8168-db30548a1c3b
# W_gom = missing           # Uncomment and complete the instruction
W_gom = growth_sim_W_gom([μ_gom, D_gom])

# ╔═╡ 2b731e4a-d6eb-48f2-91f0-27e916d85683
md"""
Using `ForwardDiff.jacobian` to compute the sensitivities for the single ouput(s).
"""

# ╔═╡ 68ad7e05-7f79-4e8d-9031-fefb5f5c0897
# sens_W_gom = missing          # Uncomment and complete the instruction
sens_W_gom = ForwardDiff.jacobian(growth_sim_W_gom, [μ_gom, D_gom])

# ╔═╡ 362a8e24-8fac-4978-bc0e-13596a05d39e
md"""
Extract the (absolute) sensitivities of the outputs on the different parameters.
"""

# ╔═╡ e9114d64-fe3e-420d-9fde-4b0fcbbc8527
# sens_W_on_μ_gom = missing      # Uncomment and complete the instruction
sens_W_on_μ_gom = sens_W_gom[:,1]

# ╔═╡ 7c9e1112-0a63-417f-a65a-9c2bdfa8bc47
# sens_W_on_D_gom = missing      # Uncomment and complete the instruction
sens_W_on_D_gom = sens_W_gom[:,2]

# ╔═╡ 94e66042-b352-4bf2-a24b-b15789e10fe3
md"""
Compute the normalized sensitivities.
"""

# ╔═╡ 1711346a-c161-415a-96f9-1235d786a584
# sens_W_on_μ_rel_gom = missing   # Uncomment and complete the instruction
sens_W_on_μ_rel_gom = sens_W_on_μ_gom .* μ_gom ./ W_gom

# ╔═╡ 7a87bfbb-a8c0-4071-8f86-66513ae40968
# sens_W_on_D_rel_gom = missing    # Uncomment and complete the instruction
sens_W_on_D_rel_gom = sens_W_on_D_gom .* D_gom ./ W_gom

# ╔═╡ 7593a4bf-ad0f-4fd3-aa3c-8b6a5073cb02
md"
Plot both sensitivity functions (with appropriate title and labels):
"

# ╔═╡ ad70a1d5-9b68-4cd6-9259-a3c25cad706b
# missing        # Uncomment and complete the instruction
plot(t_vals_gom, [sens_W_on_μ_rel_gom, sens_W_on_D_rel_gom], title="Normalized sensitivities", label=["W on μ" "W on D"], xlabel="Time (day)")

# ╔═╡ 12f64333-36d0-4e2a-9061-e3dcc7a4ae96
md"
Draw your conclusions:

- missing
- missing
"

# ╔═╡ Cell order:
# ╠═5f4fea06-0632-11ef-102e-21f5606d2056
# ╠═489b5399-fe4c-481a-834f-0101bbe28cea
# ╠═e1e7bc8e-7264-4cbc-98d2-aa73679fa2df
# ╟─427b509f-d08e-4d93-99ba-a79f9c244b28
# ╟─49e7085c-3691-4182-9630-68dc9371ad18
# ╟─c41894dd-0f6c-483f-b19f-dbf8148f776f
# ╟─7b58f2f9-cb93-4a20-bf90-62a9006b57d6
# ╟─abd70f04-1bf8-454f-869a-0d8082d68453
# ╟─e76aedec-69f2-4301-b015-6960e4503c42
# ╟─a6a97279-0b64-45c3-8312-c22b1a8425d0
# ╟─25d3ef62-382c-455c-92bb-fadfc650c5a6
# ╟─2a544bef-d7fe-4300-bc59-69e6f7429304
# ╟─7107a3a9-15ef-488e-8709-52a6444d0e1c
# ╟─04e95855-7c4a-4d2c-b836-c5dde291adad
# ╟─4673bfdf-8f4a-42bd-a026-21a66d800f2b
# ╟─8d55bc42-1f8f-4bc1-b9fa-e8152fc125ce
# ╠═22489bd4-ab64-4bf9-ad03-5372ea273935
# ╠═cdcebbb1-40e0-457f-a6ec-b769f6b1f2e9
# ╟─676119ce-f4fe-41f9-8121-b2c21f0dd28c
# ╟─5492210b-6ee6-4baf-b37e-d21358cdeb60
# ╠═18eece66-46fa-4458-aaed-f4c8fa002c20
# ╟─8bb0b24b-6cce-49cd-a625-5f375b92d9b7
# ╠═ba181db8-d176-4b83-9168-d5939ffe9661
# ╟─cd8a7ba1-194a-4b42-8c87-2c9b1fe6b475
# ╠═cbb2ca49-b019-495b-9310-83fcc00cad26
# ╠═0b3d35bb-c5b0-44b7-94b3-06fa571d339e
# ╟─be565a3c-31b6-4df1-b73b-08f308a8c09b
# ╠═f62806d1-77e1-470b-9711-33a924c788cc
# ╠═546ed163-26a7-4235-982f-7568ed609488
# ╠═0da53fa2-5a42-46e6-8bd8-45d6aa903d46
# ╟─b1298f40-4696-49d0-ac94-896e0cdbc996
# ╠═270647d2-1371-4272-8bc1-3a6ad77bc716
# ╠═ac235d86-1d93-4944-aa89-1b4fd38f0e6e
# ╠═3e13efa1-9bc6-456f-8e62-ecd3165e2a65
# ╟─79fe5411-fb78-490e-a9de-1868d1261aa6
# ╟─c57938c8-566a-4502-a2e1-69dc77291500
# ╠═65885bbe-7b73-4efd-a8e7-984a34b80548
# ╟─4aa73da9-a394-4ca3-a839-d076eb3c3d7f
# ╠═bc201b61-0f58-49d2-a20b-4f18c42fcc96
# ╟─11f15a0e-2958-4075-928c-5a24bcc00c69
# ╠═d0f3197f-3094-4f6e-b33c-e5b74e0947f6
# ╟─ba4eaa90-b97a-4b71-b3e1-968340c6def5
# ╠═71f897a6-1e91-412b-ad58-c5f1e7cd1adb
# ╟─b3b97b57-a4f7-4cb6-80c8-6ad893ade75d
# ╠═99864b21-3db0-4d60-a84e-8e96db4de4ae
# ╟─344eb3a7-8614-4349-a536-62d9acef6bda
# ╠═670ccfa4-24fa-4167-a46a-2a48dc19538b
# ╟─c262d3ed-9f62-45ad-a7ee-ec78ab16f35c
# ╠═6600d28e-2522-4069-a5e0-643be43f6117
# ╠═d9e3a4ac-c138-4b93-a19c-737166a0f0ea
# ╟─f07d8205-e172-45fa-b955-792bd95f3023
# ╠═1c142797-ada6-4d77-a22f-b981ffd38956
# ╠═cb397e0b-56f9-420f-a86e-bfab35286b44
# ╟─9e6946a4-f207-4bfb-9e19-aa7b77c2a05b
# ╠═f58ce914-f366-4087-a7d1-8cfe69ac623b
# ╟─05972c7f-f64f-4865-b0a7-f33029d0a6fa
# ╠═88aa13f6-e16c-40cd-ac10-3ede8b9fb429
# ╠═dc0507ee-0f67-4556-bdc3-1294177c86b7
# ╟─04a1d1ad-53a9-4146-be09-65aa422e3730
# ╟─79f0f1dd-850f-4dfd-b895-ff41a2d8adb8
# ╟─fa2270d0-2548-409c-a23f-4369d8bce8ec
# ╟─dca50b37-b06c-4efe-881e-cf966ebc8fd7
# ╠═2836f231-3b90-4829-9012-3ed9a09239ff
# ╟─e8db84be-31a0-415d-a194-064c4c87a293
# ╠═78d175b4-a9b2-49b5-bbc8-eb348799985b
# ╟─63952c54-3304-4500-9536-b375c5c8f280
# ╠═7dc013f7-7b35-42e4-aa77-0dca3755f389
# ╟─86785c55-b567-449c-ae76-bc15a16223bc
# ╟─05e96232-ce95-452b-98e4-e79817d45ae2
# ╠═9a705b3f-3c45-4896-9b89-e9c37d47fd89
# ╠═55a13ebc-09e7-433c-83a8-9b6805e2fbec
# ╟─f1293675-b4d0-4827-8714-65a593ecc00e
# ╠═7005d88f-9e28-47fe-9a32-dec023903dc3
# ╟─c8853554-26c5-491c-a536-4c07e8c6a986
# ╠═025c1154-b2ae-4e1c-af1b-277b24d648a4
# ╟─f58825d2-55fc-43f1-b164-9555bf9f5b84
# ╠═8044dd5f-9996-4062-a54a-8fb0eee30b46
# ╟─a1752ad4-2c11-416c-8321-9da058f9aaea
# ╠═0f58fd3d-1495-46c7-9d92-5539395cf12e
# ╟─6321774c-5bbe-42a7-bdad-0168c891b5ce
# ╠═eb1e5f97-3e09-43f0-b6ce-cc6cbc42ec2f
# ╟─b31b4d24-ae30-4cff-b6b7-348200feaa5d
# ╠═f2fa716f-d7b5-4c95-a8c6-6b16fbfe1499
# ╟─0b413b98-a800-44e9-a917-ba61467bd613
# ╠═6be6fea7-561b-4eec-b249-fa4522a5b039
# ╟─98070d41-e8c5-49e6-9eea-04c718cbff65
# ╠═cef751aa-f9c8-46da-9436-9cc2e5d7515c
# ╟─9e7bb7ec-8fe3-422e-920f-13c2ef055feb
# ╠═b13b5e77-ca60-465b-9f63-be9e5da0482f
# ╟─735f5191-259e-440f-a379-20a2a70ec72c
# ╠═689556c9-a690-4e9e-9052-7ac66e999d4d
# ╠═b85f0d09-a102-435b-a898-149ba4b29caa
# ╟─d067f814-1dd5-4e59-8cfd-3bc0b2d97612
# ╠═e01657f3-33ac-4eb8-b39e-e2441a165a8d
# ╠═7231e54e-b293-48fc-b1f0-273f45f51539
# ╟─0524a976-9d4d-4ef9-b088-bdd9102415e6
# ╠═33031b53-6f2d-4256-a0aa-0ba97beaa7af
# ╟─22b58c97-1aa2-49c4-8a0a-488c63014a90
# ╟─41abf8c4-e67e-4f66-a38a-7204a878d98d
# ╟─f311e943-c96b-4af3-8757-c662fd302a88
# ╠═1ae76036-42ef-46e9-88bf-d66d4267addc
# ╟─e190ceca-30b1-49e2-baf5-ceb62929f4c0
# ╠═acaf5b91-accc-4dfa-9371-6495a47c8736
# ╟─ed1150f1-9c75-4867-b7f5-535d605810f4
# ╠═d35ef155-61aa-4b19-872f-d6e621f96572
# ╟─3f546d58-15ff-442d-a4e8-7c0783ed2fe0
# ╟─b778f786-87bc-4103-a792-2ae5cf004d50
# ╠═edd76a4e-e541-4a83-923f-6f593d8fadad
# ╠═05f1ae6c-1229-4924-9c41-9554d6ca4261
# ╟─70df05e7-4527-4087-a762-942a33e89f74
# ╠═bf9291d5-6fc5-4359-874e-167237a22147
# ╟─9e3158c5-d4b2-456c-b920-8431e2047770
# ╠═3aab9073-cf4d-4ee5-9e8a-e3db586d9f68
# ╟─30c35371-d106-4673-8569-d07b39edbdfe
# ╠═c01fe052-3701-4e22-9a64-f681d52a445a
# ╟─4f959d21-954b-4116-a498-bafc1aada47a
# ╠═5f427df2-1b7e-46c7-adb9-66b05deac6e8
# ╟─958854fc-8e23-4364-88c5-1ec060b86a16
# ╠═be21e269-16c4-42b0-8930-2d009bd91161
# ╟─e57f5916-b753-429b-8434-8a85a67bd1fb
# ╠═69f96b57-7c8d-4e64-9e6a-9804d001b0b1
# ╟─f7f6c56b-4910-4641-a462-eacfa4b4d034
# ╠═84c1bb38-23c3-4f70-be67-7da2300a737b
# ╟─09eeab39-79a1-420a-8298-72b7002ec168
# ╠═9e15601f-4ad2-41be-8168-db30548a1c3b
# ╟─2b731e4a-d6eb-48f2-91f0-27e916d85683
# ╠═68ad7e05-7f79-4e8d-9031-fefb5f5c0897
# ╟─362a8e24-8fac-4978-bc0e-13596a05d39e
# ╠═e9114d64-fe3e-420d-9fde-4b0fcbbc8527
# ╠═7c9e1112-0a63-417f-a65a-9c2bdfa8bc47
# ╟─94e66042-b352-4bf2-a24b-b15789e10fe3
# ╠═1711346a-c161-415a-96f9-1235d786a584
# ╠═7a87bfbb-a8c0-4071-8f86-66513ae40968
# ╟─7593a4bf-ad0f-4fd3-aa3c-8b6a5073cb02
# ╠═ad70a1d5-9b68-4cd6-9259-a3c25cad706b
# ╟─12f64333-36d0-4e2a-9061-e3dcc7a4ae96
