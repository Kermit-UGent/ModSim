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

# ╔═╡ d3a09030-90f4-4493-a15b-770abd4c55dd
using SciMLSensitivity

# ╔═╡ 427b509f-d08e-4d93-99ba-a79f9c244b28
md"
# Introduction to sensitivity analysis
"

# ╔═╡ 49e7085c-3691-4182-9630-68dc9371ad18
md"
## Goal of this practicum
"

# ╔═╡ c41894dd-0f6c-483f-b19f-dbf8148f776f
md"
Sensitivity functions indicate how sensitive the model output is to a change in parameter values. When a model output is very sensitive to a certain parameter, a small change in the value of this parameter will have a large influence on the value of the model output. Sensitivity functions thus provide important information about the model and are implicitely used to estimate parameters and explicitely in the context of optimal experimental design.
"

# ╔═╡ 7b58f2f9-cb93-4a20-bf90-62a9006b57d6
md"
The sensitivity function that measures how sensitive output $y_j$ is to changes in parameter $\theta_i$ is given by the partial derivative

$$\cfrac{\partial \hat{y}_i(\theta)}{\partial \theta_j}\tag{1}$$
"

# ╔═╡ abd70f04-1bf8-454f-869a-0d8082d68453
md"
Try to understand why the above expression $(1)$ does indeed give us the information we were promised in the first paragraph. How will we be able to see from the value determined by the above expression $(1)$ that whether or not output $y_j$ is sensitive to a change in $\theta_i$?
"

# ╔═╡ e76aedec-69f2-4301-b015-6960e4503c42
md"
Answer:
"

# ╔═╡ a6a97279-0b64-45c3-8312-c22b1a8425d0
md"
Sometimes expression $(1)$ can be evaluated analytically. Usually, however, we will have to approximate the partial derivative numerically. Expression $(1)$ can be made more specific:

$$\cfrac{\partial \hat{y}_i(\theta)}{\partial \theta_j} \approx \cfrac{\hat{y}_i(\theta_j+\Delta\theta_j)-\hat{y}_i(\theta_j)}{\Delta\theta_j}$$
"

# ╔═╡ 25d3ef62-382c-455c-92bb-fadfc650c5a6
md"
Thus, to calculate the sensitivity function numerically, the model is evaluated for the parameter values $\theta_i$ and $\theta_i+\Delta\theta_i$ and the difference between these evaluations is taken.
"

# ╔═╡ 2a544bef-d7fe-4300-bc59-69e6f7429304
md"
Since quantity $(1)$ is dependent on the units, a normalized variant is often used:

$$\cfrac{\partial \hat{y}_i(\theta)}{\partial \theta_j} \cdot \cfrac{\theta_j}{\hat{y}_i} \tag{2}$$

The interpretation of $(2)$ is how much the output changes per cent if the parameter is increased by one per cent. It assumes positive model outputs and parameters, which is often the case for biochemical models.

Using the normalized variant allows you to compare all possible sensitivity functions with each other.
"

# ╔═╡ 7107a3a9-15ef-488e-8709-52a6444d0e1c
md"
We now calculate and interpret sensitivity functions for some given models. To illustrate the concepts, we first consider three simple models modelling the growth of grass.
"

# ╔═╡ 04e95855-7c4a-4d2c-b836-c5dde291adad
md"
## Grass growth models
"

# ╔═╡ 4673bfdf-8f4a-42bd-a026-21a66d800f2b
md"
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
"

# ╔═╡ 8d55bc42-1f8f-4bc1-b9fa-e8152fc125ce
md"
We will illustrate how to compute the local sensitivity functions for the logistic model. The same will be left as exercises below for the exponential and gompertz models.

**Important:**
- We will use consequently `_log`, `_exp` and `_gom` appended to relevant variables names in order to indicate their model origin **and** to prevent cell-disabling that occurs when using the same variables names in these Notebooks.
"

# ╔═╡ 676119ce-f4fe-41f9-8121-b2c21f0dd28c
md"
### Modelling logistic growth

We will start by modelling our system and simulating using the aforementioned parameters values, initial condition and timespan.
"

# ╔═╡ 5492210b-6ee6-4baf-b37e-d21358cdeb60
md"
Implementation of the system:
"

# ╔═╡ 18eece66-46fa-4458-aaed-f4c8fa002c20
growth_mod_log = @reaction_network begin
	@species W(t)=2.0            # default initial condition
	@parameters μ=0.07 Wf=10.0   # default parameter values
    #μ*W, 0 --> W
    #μ/Wf*W, W --> 0
	μ*(1-W/Wf), W --> 2W
end

# ╔═╡ 8bb0b24b-6cce-49cd-a625-5f375b92d9b7
md"
Convert the *reaction model* to check that we work with the correct differential equation:
"

# ╔═╡ ba181db8-d176-4b83-9168-d5939ffe9661
osys_log  = convert(ODESystem, growth_mod_log)

# ╔═╡ cd8a7ba1-194a-4b42-8c87-2c9b1fe6b475
md"
Setting initial conditions, timespan and parameter values:
"

# ╔═╡ cbb2ca49-b019-495b-9310-83fcc00cad26
u0_log = [:W => 2.0]

# ╔═╡ 0b3d35bb-c5b0-44b7-94b3-06fa571d339e
tspan = (0.0, 100.0)   # this will be the same for the three models

# ╔═╡ 0da53fa2-5a42-46e6-8bd8-45d6aa903d46
params_log = [:μ => 0.07, :Wf => 10.0]

# ╔═╡ b1298f40-4696-49d0-ac94-896e0cdbc996
md"
Creating and solving the ODEProblem and plotting results:
"

# ╔═╡ 270647d2-1371-4272-8bc1-3a6ad77bc716
oprob_log = ODEProblem(growth_mod_log, u0_log, tspan, params_log)
# Also possible:
# oprob_log = ODEProblem(growth_mod_log, [], tspan, [])

# ╔═╡ ac235d86-1d93-4944-aa89-1b4fd38f0e6e
osol_log = solve(oprob_log, Tsit5(), saveat=0.5)

# ╔═╡ 3e13efa1-9bc6-456f-8e62-ecd3165e2a65
plot(osol_log)

# ╔═╡ 3d81d3b8-41f4-4d92-b7d5-3cf165c827d7
md"
### Local Sensitivity Analysis (LSA)

In order to compute the local sensitivity functions, we will need to load the `SciMLSensitivity` package:
"

# ╔═╡ a421375a-0204-406b-942f-a69f537d3473
md"
As we created an ODE problem using the command `ODEProblem` before solving it, we here need to create a similar *ODE forward sensitivity problem* with the command `ODEForwardSensitivityProblem` before we can retrieve the sensitivity functions. We need to provide (in this order) the following arguments:

- The function refering to the ODE model function. This is the function created by `ODEProblem` appended by `.f`.
- The initial condition(s) between square brackets `[` `]`.
- The timespan.
- A list of parameter values between square brackets `[` `]` with the order corresponding to the oder of the parameters in the model (cf. apply the function `parameters(`*model_name*`)`).

In this example this would be:
"

# ╔═╡ 30484d51-7e92-4a23-8a8c-75afa5fd1881
oprob_sens_log = ODEForwardSensitivityProblem(oprob_log.f, [2.0], tspan, [0.07, 10.0])

# ╔═╡ 19dfb850-b27d-4f81-b86d-47c8d653a645
md"
Next we need to solve the *ODE forward sensitivity problem* using the customary `solve` command with as first argument, the function created before:
"

# ╔═╡ cf5a134f-d4b3-4bcb-8cc6-c34228392ded
osol_sens_log = solve(oprob_sens_log, Tsit5(), saveat=0.5)

# ╔═╡ 5203ad2c-9d43-4f91-8840-3468ec8a8dc9
md"
From the `osol_sens_log` variable we can readily get the time vector with `osol_sens_log.t`, this will become useful when plotting the sensitivities versus time.
"

# ╔═╡ 137c0fc9-4aad-4a66-9214-c6edc4e04747
md"
Finally, the local sensitivities are preferably *extracted* from the solution by means of the function `extract_local_sensitivities` with the previous solution as sole argument:
"

# ╔═╡ 1bd6ed0b-b463-4f0f-8d32-c67479a70c30
u_log, dp_log = extract_local_sensitivities(osol_sens_log)

# ╔═╡ bbaa2ca8-f566-48b0-ab8e-36a2319395cb
u_log'

# ╔═╡ 2b13495a-993f-4f26-a1c7-05545d986052
dp_log[1]'

# ╔═╡ 8a8074c0-0f6d-4907-8907-cdb2b147d1a0
md"
The `extract_local_sensitivities` returns both the simulation results (here: `u_log`) as well as the local sensitivies (here: `dp_log`).

To get the sensitivities of $W$ to $\mu$, of $W$ and $W_f$, you need to use indexing with `dp_log`:
- `dp_log[1]` gives the (absolute) sensitivity of $W$ to $\mu$.
- `dp_log[2]` gives the (absolute) sensitivity of $W$ to $W_f$.

Beware that `u_log` and `dp_log[i]` are row vectors. If later you want to plot the sensitivities you will need to transpose these into a column vectors! Futhermore, we will compute the normalized sensitivities.

Beware that all element wise operations need a dot in front of the operator, e.g. as in `./` and `.*`.
"

# ╔═╡ f2eff54d-b19d-4e96-ba98-a54072f49118
sens_μ_log  = dp_log[1]'./u_log'.*0.07  # sensitivity for W on μ

# ╔═╡ 21a8558a-a048-40be-8a04-358ac6fa4f61
sens_Wf_log = dp_log[2]'./u_log'.*10.0  # sensitivity for W on Wf

# ╔═╡ d1d3645f-15dd-4b62-b149-79339cae8636
md"
We are now ready to plot the two sensitivity functions. For each plot command we need to provide the time vector (first argument) and the column vector with the local sensitivity (second argument). Additionally, you can provide a title, a (legend) label and a x- and/or y-label.

Beware that if you want to execute multiple commands (here: `plot`) in a single cell, you need to put them in a `begin`-and-`end` block. Also if you want to visualize subsequent graphs in the same plot, the forthcoming plot command names should be followed by a `!`, as in `plot!(...)`.
"

# ╔═╡ e2b40912-20f9-4ab9-8696-11a202c5aa74
begin
plot(osol_sens_log.t, sens_μ_log, title="Normalized sensitivities", label="W on μ", xlabel="Time (day)")
plot!(osol_sens_log.t, sens_Wf_log, label="W on Wf")
end

# ╔═╡ 04a1d1ad-53a9-4146-be09-65aa422e3730
md"
Conclusions:
- From the sensitivity plot of $W$ on $\mu$ it can be seen that $W$ is most sensitive to $\mu$ in the time region $[0, 30]$. The latter corresponds to the region where the yield rate is largest (i.e., when the growth is largest). This makes sense because when looking at the differential equation, $\mu$ is approximately the growth rate for relatively small $W$ values.
- From the sensitivity plot of $W$ on $W_f$ it can be seen that $W$ is most sensitive to $W_f$ in the region where time values are large (cf. operating point). The latter corresponds to the region where the yield rate stagnates (i.e., when the yield reaches a steady value). This makes sense because when looking at the differential equation, $W_f$ is the steady state value.
"

# ╔═╡ 79f0f1dd-850f-4dfd-b895-ff41a2d8adb8
md"
## Exercises
"

# ╔═╡ fa2270d0-2548-409c-a23f-4369d8bce8ec
md"
### Exercise 1 - Sensitivity analysis of the exponential growth model
"

# ╔═╡ dca50b37-b06c-4efe-881e-cf966ebc8fd7
md"
Create a *reaction network object* for the exponential growth model. Name it `growth_exp`.
"

# ╔═╡ 2836f231-3b90-4829-9012-3ed9a09239ff
# Uncomment and complete the instruction
# growth_exp = @reaction_network begin
#     missing
# end
growth_exp = @reaction_network begin
	@species W(t)=2.0
	@parameters μ=0.02 Wf=10.0
    μ*Wf, ∅ --> W
    μ, W --> ∅
end

# ╔═╡ e8db84be-31a0-415d-a194-064c4c87a293
md"
Convert the system to a symbolic differential equation model (name it: `osys_exp`) and verify, by analyzing the differential equation, that your model has been correctly implemented.
"

# ╔═╡ 78d175b4-a9b2-49b5-bbc8-eb348799985b
# osys_exp = missing           # Uncomment and complete the instruction
osys_exp = convert(ODESystem, growth_exp)

# ╔═╡ 63952c54-3304-4500-9536-b375c5c8f280
md"
Initialize a vector `u0_exp` with the initial condition:
"

# ╔═╡ 7dc013f7-7b35-42e4-aa77-0dca3755f389
# u₀_exp = missing           # Uncomment and complete the instruction
u0_exp = [:W => 2.0]

# ╔═╡ 86785c55-b567-449c-ae76-bc15a16223bc
md"
We will use the same timespan as before, so no need to redefine it.
"

# ╔═╡ f1293675-b4d0-4827-8714-65a593ecc00e
md"
Initialize a vector `params_exp` with the parameter values:
"

# ╔═╡ 7005d88f-9e28-47fe-9a32-dec023903dc3
# params_exp = missing        # Uncomment and complete the instruction
params_exp = [:μ => 0.02, :Wf => 10.0]

# ╔═╡ c8853554-26c5-491c-a536-4c07e8c6a986
md"
Create the ODE problem and store it in `oprob_exp`:
"

# ╔═╡ 025c1154-b2ae-4e1c-af1b-277b24d648a4
# oprob_exp = missing         # Uncomment and complete the instruction
oprob_exp = ODEProblem(growth_exp, u0_exp, tspan, params_exp)

# ╔═╡ f58825d2-55fc-43f1-b164-9555bf9f5b84
md"
Solve the ODE problem. Use `Tsit5()` and `saveat=0.5`. Store the solution in `osol_exp`:
"

# ╔═╡ 8044dd5f-9996-4062-a54a-8fb0eee30b46
# osol_exp = missing          # Uncomment and complete the instruction
osol_exp = solve(oprob_exp, Tsit5(), saveat=0.5)

# ╔═╡ a1752ad4-2c11-416c-8321-9da058f9aaea
md"
Plot the result:
"

# ╔═╡ 0f58fd3d-1495-46c7-9d92-5539395cf12e
# missing                     # Uncomment and complete the instruction
plot(osol_exp)

# ╔═╡ 70eb7e6a-73b8-4f12-b5c2-d8f3f8538d37
md"
Create the ODE forward sensitivity problem and store it in `oprob_sens_exp`:
"

# ╔═╡ a72b6b81-5c8b-48a8-b165-adf4f64d70c3
# oprob_sens_exp = missing    # Uncomment and complete the instruction
oprob_sens_exp = ODEForwardSensitivityProblem(oprob_exp.f, [2.0], tspan, [0.02, 10.0])

# ╔═╡ 7229b4aa-3f5e-4cae-9fa0-e8312006b51e
md"
Solve the ODE forward sensitivity problem. Use `Tsit5()` and `saveat=0.5`. Store the solution in `osol_sens_exp`:
"

# ╔═╡ ee11b8c2-5d9f-48a9-acc7-21b53e5e5821
# osol_sens_exp = missing      # Uncomment and complete the instruction
osol_sens_exp = solve(oprob_sens_exp, Tsit5(), saveat=0.5)

# ╔═╡ 7247431d-af78-4533-9d2b-0787076c73a0
md"
Extract the local sensitivities. Use `u_exp` and `dp_exp` for the simulation result and the local sensitivities, respectively:
"

# ╔═╡ 5b0f1fbc-7b62-4592-9b35-9413ac404769
# u_exp, dp_exp = missing      # Uncomment and complete the instruction
u_exp, dp_exp = extract_local_sensitivities(osol_sens_exp)

# ╔═╡ 6737df04-7ba9-4439-a5c6-39c89d111924
md"
Convert the output $W$ and the (absolute) sensitivities of both $W$ on $\mu$ and $W$ on $W_f$ into column vectors, compute the normalized sensitivities and name them `sens_μ_exp` and `sens_Wf_exp` respectively:
"

# ╔═╡ b5527e15-8b30-4253-82eb-8597086a591c
# sens_μ_exp  = missing         # Uncomment and complete the instruction
sens_μ_exp  = dp_exp[1]'./u_exp'.*0.02

# ╔═╡ b9235039-dba8-425c-afc5-016465f93f9f
# sens_Wf_exp = missing         # Uncomment and complete the instruction
sens_Wf_exp = dp_exp[2]'./u_exp'.*10.0

# ╔═╡ 0524a976-9d4d-4ef9-b088-bdd9102415e6
md"
Plot both normalized sensitivity functions (with appropriate title and label):
"

# ╔═╡ 2d51b318-8f5f-4328-973e-c6f4ae33a6ab
# Uncomment and complete the instruction
# begin
# missing
# missing
# end
begin
plot(osol_sens_exp.t, sens_μ_exp, title="Normalized sensitivities", label="W on μ")
plot!(osol_sens_exp.t, sens_Wf_exp, label="W on Wf", xlabel="Time (day)")
end

# ╔═╡ 22b58c97-1aa2-49c4-8a0a-488c63014a90
md"
Draw your conclusions:

- missing
- missing
"

# ╔═╡ 41abf8c4-e67e-4f66-a38a-7204a878d98d
md"
### Exercise 2 - Sensitivity analysis of the Gompertz growth model
"

# ╔═╡ f311e943-c96b-4af3-8757-c662fd302a88
md"
Create a *reaction network object* for the exponential growth model. Name it `growth_gom`.
"

# ╔═╡ 1ae76036-42ef-46e9-88bf-d66d4267addc
# Uncomment and complete the instruction
# growth_gom = @reaction_network begin
#     missing
# end
growth_gom = @reaction_network begin
	@species W(t)=2.0
	@parameters μ=0.09 D=0.04
    μ, W --> 2*W
    D*log(W), W --> 0
end

# ╔═╡ e190ceca-30b1-49e2-baf5-ceb62929f4c0
md"
Convert the system to a symbolic differential equation model (name it: `osys_gom`) and verify, by analyzing the differential equation, that your model has been correctly implemented.
"

# ╔═╡ acaf5b91-accc-4dfa-9371-6495a47c8736
# osys_gom = missing                  # Uncomment and complete the instruction
osys_gom = convert(ODESystem, growth_gom)

# ╔═╡ ed1150f1-9c75-4867-b7f5-535d605810f4
md"
Initialize a vector `u0_gom` with the initial condition:
"

# ╔═╡ d35ef155-61aa-4b19-872f-d6e621f96572
# u₀_gom = missing                    # Uncomment and complete the instruction
u0_gom = [:W => 2.0]

# ╔═╡ 3f546d58-15ff-442d-a4e8-7c0783ed2fe0
md"
We will use the same timespan as before, so no need to redefine it.
"

# ╔═╡ 70df05e7-4527-4087-a762-942a33e89f74
md"
Initialize a vector `params_gom` with the parameter values:
"

# ╔═╡ bf9291d5-6fc5-4359-874e-167237a22147
# params_gom = missing                # Uncomment and complete the instruction
params_gom = [:μ => 0.09, :D => 0.04]

# ╔═╡ 9e3158c5-d4b2-456c-b920-8431e2047770
md"
Create the ODE problem and store it in `oprob_gom`:
"

# ╔═╡ 3aab9073-cf4d-4ee5-9e8a-e3db586d9f68
# oprob_gom = missing                 # Uncomment and complete the instruction
oprob_gom = ODEProblem(growth_gom, u0_gom, tspan, params_gom)

# ╔═╡ 30c35371-d106-4673-8569-d07b39edbdfe
md"
Solve the ODE problem. Use `Tsit5()` and `saveat=0.5`. Store the solution in `osol_gom`:
"

# ╔═╡ c01fe052-3701-4e22-9a64-f681d52a445a
# osol_gom = missing                  # Uncomment and complete the instruction
osol_gom = solve(oprob_gom, Tsit5(), saveat=0.5)

# ╔═╡ 4f959d21-954b-4116-a498-bafc1aada47a
md"
Plot the result:
"

# ╔═╡ 5f427df2-1b7e-46c7-adb9-66b05deac6e8
# missing                             # Uncomment and complete the instruction
plot(osol_gom)

# ╔═╡ 0ee903e9-a6dd-44ee-af83-22917b3ec1e5
md"
Create the ODE forward sensitivity problem and store it in `oprob_sens_gom`:
"

# ╔═╡ 738c6544-a27d-4d1c-ab60-768933191546
# oprob_sens_gom = missing            # Uncomment and complete the instruction
oprob_sens_gom = ODEForwardSensitivityProblem(oprob_gom.f, [2.0], tspan, [0.09, 0.04])

# ╔═╡ 4210916b-ffa4-44da-9b3c-6a054ed50c4c
md"
Solve the ODE forward sensitivity problem. Use `Tsit5()` and `saveat=0.5`. Store the solution in `osol_sens_gom`:
"

# ╔═╡ 2afb7420-cc52-40fa-8397-1f89ed6a5103
# osol_sens_gom = missing              # Uncomment and complete the instruction
osol_sens_gom = solve(oprob_sens_gom, Tsit5(), saveat=0.5)

# ╔═╡ 9cad30aa-329f-4512-9c57-b53d78451882
md"
Extract the local sensitivities. Use `u_gom` and `dp_gom` for the simulation result and the local sensitivities, respectively:
"

# ╔═╡ 7916bd22-02d0-409d-ac4f-442b057502cf
# u_gom, dp_gom = missing              # Uncomment and complete the instruction
u_gom, dp_gom = extract_local_sensitivities(osol_sens_gom)

# ╔═╡ 48443235-b69f-47c0-8135-bec03a4b7013
md"
Convert the output $W$ and the sensitivities of both $W$ to $\mu$ and $W$ to $D$ into column vectors, compute the normalized sensitivities and name them `sens_μ_gom` and `sens_D_gom` respectively:
"

# ╔═╡ b14a4dfb-dba5-40c7-b7b4-e192b2aadb3a
# sens_μ_gom  = missing                # Uncomment and complete the instruction
sens_μ_gom = dp_gom[1]'./u_gom'.*0.09

# ╔═╡ 6ec3f325-9b75-49cb-bba1-bafd000d2fe5
# sens_D_gom = missing                 # Uncomment and complete the instruction
sens_D_gom = dp_gom[2]'./u_gom'.*0.04

# ╔═╡ 7593a4bf-ad0f-4fd3-aa3c-8b6a5073cb02
md"
Plot both sensitivity functions (with appropriate title and label):
"

# ╔═╡ 25ed859b-bb0e-4c72-8770-16a592bc63f1
# Uncomment and complete the instruction
# begin
# missing
# missing
# end
begin
plot(osol_sens_gom.t, sens_μ_gom, title="Normalized sensitivities", label="W on μ", xlabel="Time (day)")
plot!(osol_sens_gom.t, sens_D_gom, label="W on D")
end

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
# ╠═427b509f-d08e-4d93-99ba-a79f9c244b28
# ╠═49e7085c-3691-4182-9630-68dc9371ad18
# ╠═c41894dd-0f6c-483f-b19f-dbf8148f776f
# ╠═7b58f2f9-cb93-4a20-bf90-62a9006b57d6
# ╠═abd70f04-1bf8-454f-869a-0d8082d68453
# ╠═e76aedec-69f2-4301-b015-6960e4503c42
# ╠═a6a97279-0b64-45c3-8312-c22b1a8425d0
# ╠═25d3ef62-382c-455c-92bb-fadfc650c5a6
# ╠═2a544bef-d7fe-4300-bc59-69e6f7429304
# ╠═7107a3a9-15ef-488e-8709-52a6444d0e1c
# ╠═04e95855-7c4a-4d2c-b836-c5dde291adad
# ╠═4673bfdf-8f4a-42bd-a026-21a66d800f2b
# ╠═8d55bc42-1f8f-4bc1-b9fa-e8152fc125ce
# ╠═22489bd4-ab64-4bf9-ad03-5372ea273935
# ╠═cdcebbb1-40e0-457f-a6ec-b769f6b1f2e9
# ╠═676119ce-f4fe-41f9-8121-b2c21f0dd28c
# ╠═5492210b-6ee6-4baf-b37e-d21358cdeb60
# ╠═18eece66-46fa-4458-aaed-f4c8fa002c20
# ╠═8bb0b24b-6cce-49cd-a625-5f375b92d9b7
# ╠═ba181db8-d176-4b83-9168-d5939ffe9661
# ╠═cd8a7ba1-194a-4b42-8c87-2c9b1fe6b475
# ╠═cbb2ca49-b019-495b-9310-83fcc00cad26
# ╠═0b3d35bb-c5b0-44b7-94b3-06fa571d339e
# ╠═0da53fa2-5a42-46e6-8bd8-45d6aa903d46
# ╠═b1298f40-4696-49d0-ac94-896e0cdbc996
# ╠═270647d2-1371-4272-8bc1-3a6ad77bc716
# ╠═ac235d86-1d93-4944-aa89-1b4fd38f0e6e
# ╠═3e13efa1-9bc6-456f-8e62-ecd3165e2a65
# ╠═3d81d3b8-41f4-4d92-b7d5-3cf165c827d7
# ╠═d3a09030-90f4-4493-a15b-770abd4c55dd
# ╠═a421375a-0204-406b-942f-a69f537d3473
# ╠═30484d51-7e92-4a23-8a8c-75afa5fd1881
# ╠═19dfb850-b27d-4f81-b86d-47c8d653a645
# ╠═cf5a134f-d4b3-4bcb-8cc6-c34228392ded
# ╠═5203ad2c-9d43-4f91-8840-3468ec8a8dc9
# ╠═137c0fc9-4aad-4a66-9214-c6edc4e04747
# ╠═1bd6ed0b-b463-4f0f-8d32-c67479a70c30
# ╠═bbaa2ca8-f566-48b0-ab8e-36a2319395cb
# ╠═2b13495a-993f-4f26-a1c7-05545d986052
# ╠═8a8074c0-0f6d-4907-8907-cdb2b147d1a0
# ╠═f2eff54d-b19d-4e96-ba98-a54072f49118
# ╠═21a8558a-a048-40be-8a04-358ac6fa4f61
# ╠═d1d3645f-15dd-4b62-b149-79339cae8636
# ╠═e2b40912-20f9-4ab9-8696-11a202c5aa74
# ╠═04a1d1ad-53a9-4146-be09-65aa422e3730
# ╠═79f0f1dd-850f-4dfd-b895-ff41a2d8adb8
# ╠═fa2270d0-2548-409c-a23f-4369d8bce8ec
# ╠═dca50b37-b06c-4efe-881e-cf966ebc8fd7
# ╠═2836f231-3b90-4829-9012-3ed9a09239ff
# ╠═e8db84be-31a0-415d-a194-064c4c87a293
# ╠═78d175b4-a9b2-49b5-bbc8-eb348799985b
# ╠═63952c54-3304-4500-9536-b375c5c8f280
# ╠═7dc013f7-7b35-42e4-aa77-0dca3755f389
# ╠═86785c55-b567-449c-ae76-bc15a16223bc
# ╠═f1293675-b4d0-4827-8714-65a593ecc00e
# ╠═7005d88f-9e28-47fe-9a32-dec023903dc3
# ╠═c8853554-26c5-491c-a536-4c07e8c6a986
# ╠═025c1154-b2ae-4e1c-af1b-277b24d648a4
# ╠═f58825d2-55fc-43f1-b164-9555bf9f5b84
# ╠═8044dd5f-9996-4062-a54a-8fb0eee30b46
# ╠═a1752ad4-2c11-416c-8321-9da058f9aaea
# ╠═0f58fd3d-1495-46c7-9d92-5539395cf12e
# ╠═70eb7e6a-73b8-4f12-b5c2-d8f3f8538d37
# ╠═a72b6b81-5c8b-48a8-b165-adf4f64d70c3
# ╠═7229b4aa-3f5e-4cae-9fa0-e8312006b51e
# ╠═ee11b8c2-5d9f-48a9-acc7-21b53e5e5821
# ╠═7247431d-af78-4533-9d2b-0787076c73a0
# ╠═5b0f1fbc-7b62-4592-9b35-9413ac404769
# ╠═6737df04-7ba9-4439-a5c6-39c89d111924
# ╠═b5527e15-8b30-4253-82eb-8597086a591c
# ╠═b9235039-dba8-425c-afc5-016465f93f9f
# ╠═0524a976-9d4d-4ef9-b088-bdd9102415e6
# ╠═2d51b318-8f5f-4328-973e-c6f4ae33a6ab
# ╠═22b58c97-1aa2-49c4-8a0a-488c63014a90
# ╠═41abf8c4-e67e-4f66-a38a-7204a878d98d
# ╠═f311e943-c96b-4af3-8757-c662fd302a88
# ╠═1ae76036-42ef-46e9-88bf-d66d4267addc
# ╠═e190ceca-30b1-49e2-baf5-ceb62929f4c0
# ╠═acaf5b91-accc-4dfa-9371-6495a47c8736
# ╠═ed1150f1-9c75-4867-b7f5-535d605810f4
# ╠═d35ef155-61aa-4b19-872f-d6e621f96572
# ╠═3f546d58-15ff-442d-a4e8-7c0783ed2fe0
# ╠═70df05e7-4527-4087-a762-942a33e89f74
# ╠═bf9291d5-6fc5-4359-874e-167237a22147
# ╠═9e3158c5-d4b2-456c-b920-8431e2047770
# ╠═3aab9073-cf4d-4ee5-9e8a-e3db586d9f68
# ╠═30c35371-d106-4673-8569-d07b39edbdfe
# ╠═c01fe052-3701-4e22-9a64-f681d52a445a
# ╠═4f959d21-954b-4116-a498-bafc1aada47a
# ╠═5f427df2-1b7e-46c7-adb9-66b05deac6e8
# ╠═0ee903e9-a6dd-44ee-af83-22917b3ec1e5
# ╠═738c6544-a27d-4d1c-ab60-768933191546
# ╠═4210916b-ffa4-44da-9b3c-6a054ed50c4c
# ╠═2afb7420-cc52-40fa-8397-1f89ed6a5103
# ╠═9cad30aa-329f-4512-9c57-b53d78451882
# ╠═7916bd22-02d0-409d-ac4f-442b057502cf
# ╠═48443235-b69f-47c0-8135-bec03a4b7013
# ╠═b14a4dfb-dba5-40c7-b7b4-e192b2aadb3a
# ╠═6ec3f325-9b75-49cb-bba1-bafd000d2fe5
# ╠═7593a4bf-ad0f-4fd3-aa3c-8b6a5073cb02
# ╠═25ed859b-bb0e-4c72-8770-16a592bc63f1
# ╠═12f64333-36d0-4e2a-9061-e3dcc7a4ae96
