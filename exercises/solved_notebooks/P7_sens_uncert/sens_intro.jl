### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ e1e7bc8e-7264-4cbc-98d2-aa73679fa2df
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg; Pkg.activate("..")
	using Catalyst
	using OrdinaryDiffEq, StatsPlots
	using PlutoUI; TableOfContents()
end

# ╔═╡ 5f4fea06-0632-11ef-102e-21f5606d2056
using Markdown

# ╔═╡ 489b5399-fe4c-481a-834f-0101bbe28cea
using InteractiveUtils

# ╔═╡ 65885bbe-7b73-4efd-a8e7-984a34b80548
using ForwardDiff

# ╔═╡ a493c9e6-e250-4b7b-8cf2-03e3ac40cf19
using GlobalSensitivity

# ╔═╡ 427b509f-d08e-4d93-99ba-a79f9c244b28
md"""
# Introduction to sensitivity analysis
"""

# ╔═╡ 49e7085c-3691-4182-9630-68dc9371ad18
md"""
## Goal of this practical
"""

# ╔═╡ faf316bb-0ccc-4471-a53b-2476c3646705
md"### Sensitivity functions"

# ╔═╡ c41894dd-0f6c-483f-b19f-dbf8148f776f
md"""
Sensitivity functions indicate how sensitive the model output is to a change in parameter values. When a model output is very sensitive to a certain parameter, a small change in the value of this parameter will have a large influence on the value of the model output. Sensitivity functions thus provide important information about the model and are implicitely used to estimate parameters and explicitely in the context of optimal experimental design.
"""

# ╔═╡ 7b58f2f9-cb93-4a20-bf90-62a9006b57d6
md"""
The sensitivity function that measures how sensitive output $y_j$ is to changes in parameter $\theta_i$ is given by the partial derivative

$$\cfrac{\partial \hat{y}_i(\theta)}{\partial \theta_j}\tag{1}$$
"""

# ╔═╡ ec1ee750-1f35-4f51-a8c7-d42a18a7b233
md"""
Try to understand why the above expression $(1)$ does indeed give us the information we were promised in the first paragraph. How will we be able to see from the value determined by the above expression $(1)$ that whether or not output $y_j$ is sensitive to a change in $\theta_i$?
"""

# ╔═╡ f4ee83d1-259e-4dea-a86c-8c909ae4ec61
md"""
!!! example
	Consider the function $f(x) = a \, e^{b \, x}$, with parameters $\vec{θ} = (a, b)$.

	The sensitivity to $a$ is
	
	```math
	\dfrac{\partial f}{\partial a} = e^{b \, x} \, ,
	```
	and for $b$:
	
	```math
	\dfrac{\partial f}{\partial b} = a \, x \, e^{b \, x} \, .
	```

	If we now take a set of values for the input and parameters, such as $(x = 2, a = 1, b = 1)$, we can calculate the sensitivities at this point:
	- For $a$: $\frac{\partial f}{\partial a} = e^2$
	- For $b$: $\frac{\partial f}{\partial b} = 2 e^2$
	from which we learn the function is twice as sensitive to $b$ compared to $a$ here.
"""

# ╔═╡ a6a97279-0b64-45c3-8312-c22b1a8425d0
md"""
For simple functions, $(1)$ can be evaluated analytically. Usually, however, we will have to approximate the partial derivative numerically. Expression $(1)$ can be made more specific:

$$\cfrac{\partial \hat{y}_i(\theta)}{\partial \theta_j} \approx \cfrac{\hat{y}_i(\theta_j+\Delta\theta_j)-\hat{y}_i(\theta_j)}{\Delta\theta_j}$$
"""

# ╔═╡ 25d3ef62-382c-455c-92bb-fadfc650c5a6
md"""
Thus, to calculate the sensitivity function numerically, the model is evaluated for the parameter values $\theta_i$ and $\theta_i+\Delta\theta_i$ and the difference between these evaluations is taken.
"""

# ╔═╡ 33a71f15-dc39-45ed-b26c-855f42c55cae
md"### Normalized sensitivity"

# ╔═╡ 2a544bef-d7fe-4300-bc59-69e6f7429304
md"""
Since quantity $(1)$ is dependent on the units, a normalized variant is often used:

$$\cfrac{\partial \hat{y}_i(\theta)}{\partial \theta_j} \cdot \cfrac{\theta_j}{\hat{y}_i} \tag{2}$$

The interpretation of $(2)$ is how much the output changes per cent if the parameter is increased by one per cent. It assumes positive model outputs and parameters, which is often the case for biochemical models.

Using the normalized variant allows you to compare all possible sensitivity functions with each other.
"""

# ╔═╡ ec590df1-91be-40a2-844e-4ec3ea5b1077
md"### Global sensitivity analysis"

# ╔═╡ 3e2b5d93-4385-4261-9486-d93511ee7c26
md"""
The sensitivity function evaluates the partial derivative of your function **given some values for your parameters** (see also the given example). Therefore, the results of your sensitivity analysis will change based on the value of θ. As these methods are dependent on the current location of θ, they are called **local sensitivity analysis** methods.

**Global sensitivity methods** try to circumvent this problem by quantifying sensitivity using **multiple** values of the output $y$. We consider two methods in this course: the Sobol method and the Morris method. For more information regarding these methods, we refer to the theory.
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

For all simulations (regardless of model), we will use an initial condition of $W_0 = 2.0$, a simulation time of $100$ days and save our solution every $0.5$ days.
"""

# ╔═╡ 0b3d35bb-c5b0-44b7-94b3-06fa571d339e
begin
	const u0 = [:W => 2.0]
	const tspan = (0.0, 100.0)
	const Δt = 0.5
	const t_vals = tspan[1]:Δt:tspan[2] # time values of the simulation
end;

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
osys_log = convert(ODESystem, growth_log)

# ╔═╡ be565a3c-31b6-4df1-b73b-08f308a8c09b
md"""
For the sake of clarity, we will use the variables `μ_log` and `Wf_log` to store the parameter values.
"""

# ╔═╡ dca01fd2-28f9-4d36-9b8f-b209ab97d3db
μ_log = 0.07;

# ╔═╡ 292aff06-7dfa-424b-8a90-173c37b9e67b
Wf_log = 10.0;

# ╔═╡ 0da53fa2-5a42-46e6-8bd8-45d6aa903d46
params_log = [:μ => μ_log, :Wf => Wf_log];

# ╔═╡ b1298f40-4696-49d0-ac94-896e0cdbc996
md"""
Creating and solving the ODEProblem and plotting results:
"""

# ╔═╡ 270647d2-1371-4272-8bc1-3a6ad77bc716
oprob_log = ODEProblem(growth_log, u0, tspan, params_log);
# Also possible:
# oprob_log = ODEProblem(growth_mod_log, [], tspan, [])

# ╔═╡ ac235d86-1d93-4944-aa89-1b4fd38f0e6e
osol_log = solve(oprob_log, Tsit5(), saveat = Δt);

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
    oprob_log = ODEProblem(growth_log, u0, tspan, [:μ=>μ, :Wf=>Wf])
	osol_log = solve(oprob_log, Tsit5(), saveat = Δt)
	return osol_log
end

# ╔═╡ 11f15a0e-2958-4075-928c-5a24bcc00c69
md"""
Next we will need to make a function based on the solution function that returns the output of a single variable.
"""

# ╔═╡ d0f3197f-3094-4f6e-b33c-e5b74e0947f6
growth_sim_W_log(params) = growth_sim_log(params)[:W];

# ╔═╡ b3b97b57-a4f7-4cb6-80c8-6ad893ade75d
md"""
Compute the single output with the given parameter values.
"""

# ╔═╡ 99864b21-3db0-4d60-a84e-8e96db4de4ae
W_log = growth_sim_W_log([μ_log, Wf_log]);

# ╔═╡ 344eb3a7-8614-4349-a536-62d9acef6bda
md"""
Use the function `ForwardDiff.jacobian` to compute the sensitivities. This function takes two arguments: the solution function and a vector with the parameter values.
"""

# ╔═╡ 670ccfa4-24fa-4167-a46a-2a48dc19538b
sens_W_log = ForwardDiff.jacobian(growth_sim_W_log, [μ_log, Wf_log]);

# ╔═╡ c262d3ed-9f62-45ad-a7ee-ec78ab16f35c
md"""
To get the sensitivities of $W$ on $\mu$, and of $W$ on $W_f$, you need to use indexing with `sens_W_log`:
- `sens_W_log[:,1]` gives the (absolute) sensitivity of $W$ on $\mu$.
- `sens_W_log[:,2]` gives the (absolute) sensitivity of $W$ on $W_f$.
"""

# ╔═╡ 6600d28e-2522-4069-a5e0-643be43f6117
sens_W_on_μ_log = sens_W_log[:,1];    # sensitivity of W on μ

# ╔═╡ d9e3a4ac-c138-4b93-a19c-737166a0f0ea
sens_W_on_Wf_log = sens_W_log[:,2];   # sensitivity of W on Wf

# ╔═╡ f07d8205-e172-45fa-b955-792bd95f3023
md"""
We now calculate the normalized sensitivities. For that we need to multiply by the parameter value and divide by the ouput. Beware that all element wise operations need a dot in front of the operator, e.g. as in `.*` and `./`.
"""

# ╔═╡ 1c142797-ada6-4d77-a22f-b981ffd38956
sens_W_on_μ_rel_log = sens_W_on_μ_log .* μ_log ./ W_log;

# ╔═╡ cb397e0b-56f9-420f-a86e-bfab35286b44
sens_W_on_Wf_rel_log = sens_W_on_Wf_log .* Wf_log ./ W_log;

# ╔═╡ 9e6946a4-f207-4bfb-9e19-aa7b77c2a05b
md"""
We can now plot the two sensitivity functions as follows:
"""

# ╔═╡ dc0507ee-0f67-4556-bdc3-1294177c86b7
begin
	plot(
		t_vals, sens_W_on_μ_rel_log, title = "Normalized sensitivities",
		label = "W on μ", xlabel = "Time (day)"
	)
	plot!(t_vals, sens_W_on_Wf_rel_log, label="W on Wf")
end

# ╔═╡ 04a1d1ad-53a9-4146-be09-65aa422e3730
md"""
Conclusions:
- From the sensitivity plot of $W$ on $\mu$ it can be seen that $W$ is most sensitive to $\mu$ in the time region $[0, 30]$. The latter corresponds to the region where the yield rate is largest (i.e., when the growth is largest). This makes sense because when looking at the differential equation, $\mu$ is approximately the growth rate for relatively small $W$ values.
- From the sensitivity plot of $W$ on $W_f$ it can be seen that $W$ is most sensitive to $W_f$ in the region where time values are large (cf. operating point). The latter corresponds to the region where the yield rate stagnates (i.e., when the yield reaches a steady value). This makes sense because when looking at the differential equation, $W_f$ is the steady state value.
"""

# ╔═╡ 5152abfc-98fa-4e10-aa28-6455b996a256
md"### Global sensitivity analysis (GSA)"

# ╔═╡ 8af97e64-6a14-4c7d-b935-002af339c94a
md"""
To perform global sensitivity analysis, we use the `gsa` function from the `GlobalSensitivity` package. This function requires 4 inputs:
1. Similarly to local sensitivity analysis, we require a function that takes the model parameters as its input and returns a model output of interest. We can just reuse `growth_sim_log` here to get the model outputs at all timepoints.
1. A GSA method. The methods considered in this course are the Sobol method `Sobol()` and the Morris method `Morris()`. We will use Sobol here.
1. A range for every parameter from where values will be sampled. For our grass growth models, we will assume realistic values for $μ$ are $[0, 0.2]$ and for $W_f$ are $[1, 20]$.
1. The amount of samples. As the Sobol method is computationally expensive, you typically don't want to use too many for an exercise. For the Morris method on the other hand, you typically want a large amount of samples for stable results (> `1_000_000`).
"""

# ╔═╡ b86d3fe4-c3e5-4c39-8118-7970d609c17a
p_range = [[0.0, 0.2], [1.0, 20.0]];

# ╔═╡ 66083965-8402-47b6-9ccd-9ced9c7412e7
gsa_res = gsa(growth_sim_log, Sobol(), p_range, samples = 100)

# ╔═╡ 97eed321-4dfe-468b-8b4e-6748608a781d
md"""
!!! tip
	Does a function return some weird type you don't know what to do with? You can find out the fields of a variable `x` using `fieldnames(typeof(x))`.
"""

# ╔═╡ 8c923fea-2e75-42c1-ae73-0bdf8d1bee91
fieldnames(typeof(gsa_res))

# ╔═╡ 3016e5f5-2f13-4877-9ed8-e6a355e63d5e
md"""
The resulting variable contains the first, second and total order Sobol indices, which we can access using respectively `gsa_res.S1`, `gsa_res.S2` and `gsa_res.ST`. We will plot the total effects here:
"""

# ╔═╡ e2ec4060-342e-41c6-8bb2-1a20192ceeff
begin
	plot(t_vals, gsa_res.ST[1][1, :], label = "μ", xlabel = "Time (day)",
		title = "Global sensitivities", ylims = (0, 1)
	)
	plot!(t_vals, gsa_res.ST[2][1, :], label = "Wf")
end

# ╔═╡ 73990cc1-77d5-4f09-ab2e-fb67c7ff81a4
md"""
!!! question
	GSA returns an appreciably positive value for the effect of growth rate μ on final grass yield, while it was approximately $0$ in LSA. Why is this effect larger when we perform a global analysis where we consider multiple possible values for $μ$ and $W_f$?
"""

# ╔═╡ 51488d93-48ad-4a02-9b51-00ef26f64699
md"""
**Answer**: missing
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

# ╔═╡ 05e96232-ce95-452b-98e4-e79817d45ae2
md"""
For the sake of clarity, we will use the variables `μ_exp` and `Wf_exp` to store the parameter values.
"""

# ╔═╡ 9a705b3f-3c45-4896-9b89-e9c37d47fd89
μ_exp = 0.02;

# ╔═╡ 55a13ebc-09e7-433c-83a8-9b6805e2fbec
Wf_exp = 10.0;

# ╔═╡ f1293675-b4d0-4827-8714-65a593ecc00e
md"""
Initialize a vector `params_exp` with the parameter values:
"""

# ╔═╡ 7005d88f-9e28-47fe-9a32-dec023903dc3
# params_exp = missing        # Uncomment and complete the instruction
params_exp = [:μ => μ_exp, :Wf => Wf_exp];

# ╔═╡ c8853554-26c5-491c-a536-4c07e8c6a986
md"""
Create the ODE problem and store it in `oprob_exp`, using the same initial condition and timespan as before.
"""

# ╔═╡ 025c1154-b2ae-4e1c-af1b-277b24d648a4
# oprob_exp = missing         # Uncomment and complete the instruction
oprob_exp = ODEProblem(growth_exp, u0, tspan, params_exp);

# ╔═╡ f58825d2-55fc-43f1-b164-9555bf9f5b84
md"""
Solve the ODE problem. Use `Tsit5()` and `saveat=0.5`. Store the solution in `osol_exp`:
"""

# ╔═╡ 8044dd5f-9996-4062-a54a-8fb0eee30b46
# osol_exp = missing          # Uncomment and complete the instruction
osol_exp = solve(oprob_exp, Tsit5(), saveat=0.5);

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
    oprob_exp = ODEProblem(growth_exp, u0, tspan, [:μ=>μ, :Wf=>Wf])
	osol_exp = solve(oprob_exp, Tsit5(), saveat=0.5)
	return osol_exp
end

# ╔═╡ b31b4d24-ae30-4cff-b6b7-348200feaa5d
md"""
Make a function based on the solution function that returns the output of a single variable.
"""

# ╔═╡ f2fa716f-d7b5-4c95-a8c6-6b16fbfe1499
# growth_sim_W_exp(params) = missing   # Uncomment and complete the instruction
growth_sim_W_exp(params) = growth_sim_exp(params)[:W];

# ╔═╡ 98070d41-e8c5-49e6-9eea-04c718cbff65
md"""
Compute the output for the given parameter values.
"""

# ╔═╡ cef751aa-f9c8-46da-9436-9cc2e5d7515c
# W_exp = missing              # Uncomment and complete the instruction
W_exp = growth_sim_W_exp([μ_exp, Wf_exp]);

# ╔═╡ 9e7bb7ec-8fe3-422e-920f-13c2ef055feb
md"""
Using `ForwardDiff.jacobian` to compute the sensitivities for the single ouput(s).
"""

# ╔═╡ b13b5e77-ca60-465b-9f63-be9e5da0482f
# sens_W_exp = missing           # Uncomment and complete the instruction
sens_W_exp = ForwardDiff.jacobian(growth_sim_W_exp, [μ_exp, Wf_exp]);

# ╔═╡ 735f5191-259e-440f-a379-20a2a70ec72c
md"""
Extract the (absolute) sensitivities of the outputs on the different parameters.
"""

# ╔═╡ 689556c9-a690-4e9e-9052-7ac66e999d4d
# sens_W_on_μ_exp = missing      # Uncomment and complete the instruction
sens_W_on_μ_exp = sens_W_exp[:, 1];

# ╔═╡ b85f0d09-a102-435b-a898-149ba4b29caa
# sens_W_on_Wf_exp = missing     # Uncomment and complete the instruction
sens_W_on_Wf_exp = sens_W_exp[:, 2];

# ╔═╡ d067f814-1dd5-4e59-8cfd-3bc0b2d97612
md"""
Compute the normalized sensitivities.
"""

# ╔═╡ e01657f3-33ac-4eb8-b39e-e2441a165a8d
# sens_W_on_μ_rel_exp = missing   # Uncomment and complete the instruction
sens_W_on_μ_rel_exp = sens_W_on_μ_exp .* μ_exp ./ W_exp;

# ╔═╡ 7231e54e-b293-48fc-b1f0-273f45f51539
# sens_W_on_Wf_rel_exp = missing   # Uncomment and complete the instruction
sens_W_on_Wf_rel_exp = sens_W_on_Wf_exp .* Wf_exp ./ W_exp;

# ╔═╡ 0524a976-9d4d-4ef9-b088-bdd9102415e6
md"""
Plot both normalized sensitivity functions (with appropriate title and labels):
"""

# ╔═╡ 33031b53-6f2d-4256-a0aa-0ba97beaa7af
# missing     # Uncomment and complete the instruction
plot(t_vals, [sens_W_on_μ_rel_exp, sens_W_on_Wf_rel_exp], title="Normalized sensitivities", label=["W on μ" "W on Wf"], xlabel="Time (day)")

# ╔═╡ 22b58c97-1aa2-49c4-8a0a-488c63014a90
md"
Draw your conclusions:

- missing
- missing
"

# ╔═╡ ee17e112-a422-4b0e-be4c-c9e6abb040cb
md"""
Perform a global sensitivity analysis on this model, now using the **Morris method** and plot the results. You can use the same parameter range as before.
"""

# ╔═╡ a8c930a2-a119-4681-a063-7a6c99615c5b
gsa_res_exp = gsa(growth_sim_exp, Morris(num_trajectory = 50), p_range, samples = 1_000_000);

# ╔═╡ 1f9d7f84-fae7-4f0e-8783-7a8ce44051a3
begin
	plot(t_vals, gsa_res_exp.means[1][1, :], label = "μ", xlabel = "Time (day)",
		title = "Global sensitivities"
	)
	plot!(t_vals, gsa_res_exp.means[2][1, :], label = "Wf")
end

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

# ╔═╡ 70df05e7-4527-4087-a762-942a33e89f74
md"""
Initialize a vector `params_gom` with the parameter values:
"""

# ╔═╡ e9408c0f-cbcb-4703-8b85-e349e37b2025
μ_gom = 0.09;

# ╔═╡ 24c4e348-7322-4c70-b8da-ec70d25356a9
D_gom = 0.04;

# ╔═╡ bf9291d5-6fc5-4359-874e-167237a22147
# params_gom = missing                # Uncomment and complete the instruction
params_gom = [:μ => μ_gom, :D => D_gom];

# ╔═╡ 9e3158c5-d4b2-456c-b920-8431e2047770
md"""
Create the ODE problem and store it in `oprob_gom`, using the same initial conditions and timespan as before:
"""

# ╔═╡ 3aab9073-cf4d-4ee5-9e8a-e3db586d9f68
# oprob_gom = missing                 # Uncomment and complete the instruction
oprob_gom = ODEProblem(growth_gom, u0, tspan, params_gom);

# ╔═╡ 30c35371-d106-4673-8569-d07b39edbdfe
md"""
Solve the ODE problem. Use `Tsit5()` and `saveat=0.5`. Store the solution in `osol_gom`:
"""

# ╔═╡ c01fe052-3701-4e22-9a64-f681d52a445a
# osol_gom = missing                  # Uncomment and complete the instruction
osol_gom = solve(oprob_gom, Tsit5(), saveat=0.5);

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
    oprob_gom = ODEProblem(growth_gom, u0, tspan, [:μ=>μ, :D=>D])
	osol_gom = solve(oprob_gom, Tsit5(), saveat=0.5)
	return osol_gom
end

# ╔═╡ e57f5916-b753-429b-8434-8a85a67bd1fb
md"""
Make a function based on the solution function that returns a single output.
"""

# ╔═╡ 69f96b57-7c8d-4e64-9e6a-9804d001b0b1
# growth_sim_W_gom(params) = missing  # Uncomment and complete the instruction
growth_sim_W_gom(params) = growth_sim_gom(params)[:W];

# ╔═╡ 09eeab39-79a1-420a-8298-72b7002ec168
md"""
Compute the output for the given parameter values.
"""

# ╔═╡ 9e15601f-4ad2-41be-8168-db30548a1c3b
# W_gom = missing           # Uncomment and complete the instruction
W_gom = growth_sim_W_gom([μ_gom, D_gom]);

# ╔═╡ 2b731e4a-d6eb-48f2-91f0-27e916d85683
md"""
Using `ForwardDiff.jacobian` to compute the sensitivities for the single ouput(s).
"""

# ╔═╡ 68ad7e05-7f79-4e8d-9031-fefb5f5c0897
# sens_W_gom = missing          # Uncomment and complete the instruction
sens_W_gom = ForwardDiff.jacobian(growth_sim_W_gom, [μ_gom, D_gom]);

# ╔═╡ 362a8e24-8fac-4978-bc0e-13596a05d39e
md"""
Extract the (absolute) sensitivities of the outputs on the different parameters.
"""

# ╔═╡ e9114d64-fe3e-420d-9fde-4b0fcbbc8527
# sens_W_on_μ_gom = missing      # Uncomment and complete the instruction
sens_W_on_μ_gom = sens_W_gom[:,1];

# ╔═╡ 7c9e1112-0a63-417f-a65a-9c2bdfa8bc47
# sens_W_on_D_gom = missing      # Uncomment and complete the instruction
sens_W_on_D_gom = sens_W_gom[:,2];

# ╔═╡ 94e66042-b352-4bf2-a24b-b15789e10fe3
md"""
Compute the normalized sensitivities.
"""

# ╔═╡ 1711346a-c161-415a-96f9-1235d786a584
# sens_W_on_μ_rel_gom = missing   # Uncomment and complete the instruction
sens_W_on_μ_rel_gom = sens_W_on_μ_gom .* μ_gom ./ W_gom;

# ╔═╡ 7a87bfbb-a8c0-4071-8f86-66513ae40968
# sens_W_on_D_rel_gom = missing    # Uncomment and complete the instruction
sens_W_on_D_rel_gom = sens_W_on_D_gom .* D_gom ./ W_gom;

# ╔═╡ 7593a4bf-ad0f-4fd3-aa3c-8b6a5073cb02
md"
Plot both sensitivity functions (with appropriate title and labels):
"

# ╔═╡ ad70a1d5-9b68-4cd6-9259-a3c25cad706b
# missing        # Uncomment and complete the instruction
plot(t_vals, [sens_W_on_μ_rel_gom, sens_W_on_D_rel_gom], title="Normalized sensitivities", label=["W on μ" "W on D"], xlabel="Time (day)")

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
# ╟─faf316bb-0ccc-4471-a53b-2476c3646705
# ╟─c41894dd-0f6c-483f-b19f-dbf8148f776f
# ╟─7b58f2f9-cb93-4a20-bf90-62a9006b57d6
# ╟─ec1ee750-1f35-4f51-a8c7-d42a18a7b233
# ╟─f4ee83d1-259e-4dea-a86c-8c909ae4ec61
# ╟─a6a97279-0b64-45c3-8312-c22b1a8425d0
# ╟─25d3ef62-382c-455c-92bb-fadfc650c5a6
# ╟─33a71f15-dc39-45ed-b26c-855f42c55cae
# ╟─2a544bef-d7fe-4300-bc59-69e6f7429304
# ╟─ec590df1-91be-40a2-844e-4ec3ea5b1077
# ╟─3e2b5d93-4385-4261-9486-d93511ee7c26
# ╟─04e95855-7c4a-4d2c-b836-c5dde291adad
# ╟─4673bfdf-8f4a-42bd-a026-21a66d800f2b
# ╠═0b3d35bb-c5b0-44b7-94b3-06fa571d339e
# ╟─8d55bc42-1f8f-4bc1-b9fa-e8152fc125ce
# ╟─676119ce-f4fe-41f9-8121-b2c21f0dd28c
# ╟─5492210b-6ee6-4baf-b37e-d21358cdeb60
# ╠═18eece66-46fa-4458-aaed-f4c8fa002c20
# ╟─8bb0b24b-6cce-49cd-a625-5f375b92d9b7
# ╠═ba181db8-d176-4b83-9168-d5939ffe9661
# ╟─be565a3c-31b6-4df1-b73b-08f308a8c09b
# ╠═dca01fd2-28f9-4d36-9b8f-b209ab97d3db
# ╠═292aff06-7dfa-424b-8a90-173c37b9e67b
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
# ╠═dc0507ee-0f67-4556-bdc3-1294177c86b7
# ╟─04a1d1ad-53a9-4146-be09-65aa422e3730
# ╟─5152abfc-98fa-4e10-aa28-6455b996a256
# ╠═a493c9e6-e250-4b7b-8cf2-03e3ac40cf19
# ╟─8af97e64-6a14-4c7d-b935-002af339c94a
# ╠═b86d3fe4-c3e5-4c39-8118-7970d609c17a
# ╠═66083965-8402-47b6-9ccd-9ced9c7412e7
# ╟─97eed321-4dfe-468b-8b4e-6748608a781d
# ╠═8c923fea-2e75-42c1-ae73-0bdf8d1bee91
# ╟─3016e5f5-2f13-4877-9ed8-e6a355e63d5e
# ╠═e2ec4060-342e-41c6-8bb2-1a20192ceeff
# ╟─73990cc1-77d5-4f09-ab2e-fb67c7ff81a4
# ╠═51488d93-48ad-4a02-9b51-00ef26f64699
# ╟─79f0f1dd-850f-4dfd-b895-ff41a2d8adb8
# ╟─fa2270d0-2548-409c-a23f-4369d8bce8ec
# ╟─dca50b37-b06c-4efe-881e-cf966ebc8fd7
# ╠═2836f231-3b90-4829-9012-3ed9a09239ff
# ╟─e8db84be-31a0-415d-a194-064c4c87a293
# ╠═78d175b4-a9b2-49b5-bbc8-eb348799985b
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
# ╟─ee17e112-a422-4b0e-be4c-c9e6abb040cb
# ╠═a8c930a2-a119-4681-a063-7a6c99615c5b
# ╠═1f9d7f84-fae7-4f0e-8783-7a8ce44051a3
# ╟─41abf8c4-e67e-4f66-a38a-7204a878d98d
# ╟─f311e943-c96b-4af3-8757-c662fd302a88
# ╠═1ae76036-42ef-46e9-88bf-d66d4267addc
# ╟─e190ceca-30b1-49e2-baf5-ceb62929f4c0
# ╠═acaf5b91-accc-4dfa-9371-6495a47c8736
# ╟─70df05e7-4527-4087-a762-942a33e89f74
# ╠═e9408c0f-cbcb-4703-8b85-e349e37b2025
# ╠═24c4e348-7322-4c70-b8da-ec70d25356a9
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
