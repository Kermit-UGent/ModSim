### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 8349306c-e98c-4221-9a1b-322fde3e18cb
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ d2c4d230-0943-11ef-3aad-5719e74bb20e
using Markdown

# ╔═╡ 47790920-74d2-4f21-b0f2-45dc19d25156
using InteractiveUtils

# ╔═╡ 83f3d978-aefa-40c1-a647-b26e837aeed6
using Catalyst

# ╔═╡ 73d5f921-82e4-4160-a3c4-fc1935f4c58d
using DifferentialEquations, Plots

# ╔═╡ 5185d0eb-7392-4775-9336-3e0f9e1449ce
using Measurements

# ╔═╡ 241a8a65-c59f-44f1-be39-d5edd1321b49
md"
# Introduction to uncertainty analysis
"

# ╔═╡ 0f518e36-bc96-4599-807e-728504e5ca7b
md"
## Goal of this practicum
"

# ╔═╡ da3cdc89-b911-4af9-9d4a-526301cba581
md"
Parameter uncertainty plays a crucial role in shaping the behavior of output variables. Model equations describe how systems evolve over time, often incorporating parameters representing various aspects of the system's characteristics. However, these parameters are rarely known with absolute certainty and they often come with inherent uncertainty due to measurement errors, variability in real-world conditions, or incomplete knowledge about the system. This uncertainty can propagate through the model equations, leading to uncertainties in the predicted outcomes. Consequently, understanding the influence of parameter uncertainty becomes essential for assessing the reliability and robustness of the model predictions, as well as for making informed decisions based on these predictions. Techniques such as sensitivity analysis and uncertainty quantification are employed to explore and quantify the impact of parameter uncertainty on the output variables, providing insights into the system's behavior and guiding the refinement of models for improved accuracy and reliability.
"

# ╔═╡ 78afaded-5a19-4386-aa76-7974977ea354
md"
Uncertainty in model parameters manifests as variability in the predicted outcomes, resulting in error bars around the output variables. These error bars represent the range of potential values that the output variables could take due to the uncertainty in the parameters. As the uncertainty in parameters increases, the width of these error bars typically expands, reflecting the increased variability and unpredictability in the model's predictions.
"

# ╔═╡ eae11742-14c6-4b1b-a939-232710dfa10e
md"
We will now compute the variability in the output variables replected as error bars assuming some uncertainty in the model parameters. To illustrate this concept, we first revisit the three simple models modelling the growth of grass.
"

# ╔═╡ e5ab4490-6fd4-4b51-bfa7-1366438efefc
md"
## Grass growth models
"

# ╔═╡ 3e2750e7-220a-47b2-b445-eb3315734dec
md"
In this notebook, three different models will be used, each modelling the yield of grass in a grassland:

- Logistic growth model: $\cfrac{dW}{dt} = \mu \left( 1 - \cfrac{W}{W_f} \right) W$
- Exponential growth model: $\cfrac{dW}{dt} = \mu \left( W_f - W \right)$
- Gompertz growth model: $\cfrac{dW}{dt} = \left( \mu - D \ln(W) \right) W$

with output $W$ the grass yield, and $W_f$, $\mu$ and $D$ parameters. The table below show some typical values for the parameters together with their uncertainties:

|             | $\mu$      | $W_f$       | $D$          |
|:----------- |:----------:|:-----------:|:------------:|
| Logistic    |  0.07$\pm$0.02     | 10.0$\pm$0.15        |              |
| Exponential |  0.02$\pm$0.01      | 10.0$\pm$0.15        |              |
| Gompertz    |  0.09$\pm$0.01      |             | 0.040$\pm$0.002         |

We will use an initial condition of $W_0 = 2.0$ for each and a simulation time of $100$ days.
"

# ╔═╡ c874a08c-7d82-4631-b611-c598dcaded09
md"
We will illustrate how to compute the error bars, due to model parameter uncertainty, in conjunction with the output variable simulation for the logistic model. The same will be left as exercises below for the exponential and gompertz models.

**Important:**
- We will use consequently `_log`, `_exp` and `_gom` appended to relevant variables names in order to indicate their model origin **and** to prevent cell-disabling that occurs when using the same variables names in these Notebooks.
"

# ╔═╡ 9d3192bc-ee69-44a1-9273-06b8a55c60ef
md"
### Uncertainty analysis of the logistic growth model

We will start by modelling our system and simulating using the aforementioned parameters values and uncertainties, initial condition and timespan.
"

# ╔═╡ 8acb7ea2-54ef-428b-bb4f-364377d2c38d
md"
Implementation of the system:
"

# ╔═╡ 460d98ef-2177-49bf-87a4-35412f4183ed
growth_mod_log = @reaction_network begin
    #μ*W, 0 --> W
    #μ/Wf*W, W --> 0
	μ*(1-W/Wf), W --> 2W
end

# ╔═╡ 38e0f746-1d04-46ba-81cc-21522b9ce6a4
md"
Convert the *reaction model* to check that we work with the correct differential equation:
"

# ╔═╡ 2538bd6e-027f-47fe-813c-db0a02586d2c
osys_log  = convert(ODESystem, growth_mod_log)

# ╔═╡ 72b53ba6-dd05-4f2d-a65e-52783e76a1e9
md"
In order to use uncertainties in the parameter values, we will need to load the package `Measurements`:
"

# ╔═╡ 1a67db96-6973-446d-aa5b-253dd0008a73
md"
Setting initial conditions:
"

# ╔═╡ 048112b0-10f1-47c5-834f-c3851d3078e8
u0_log = [:W => 2.0]

# ╔═╡ 037cbf2b-dc0c-4b62-879b-aa702c964877
md"
Set the timespan for the simulation:
"

# ╔═╡ 40c1192c-bc76-4d3a-8144-dd76ba72fba2
tspan = (0.0, 100.0)

# ╔═╡ 0b5493eb-6ce6-43c1-9289-4c2e6f20d014
md"
In order to see the effect of the individual parameter uncertainties on the output variable, we will analyse this (only for didactical reasons) assuming different study cases:

1.  $\mu$ has a nominal value of $0.07$ and a standard deviation of $0.02$ while $W_f$ is exactly known and equal to $0.15$.
2.  $K_s$ has a nominal value of $10.0$ and a standard deviation of $0.15$ while $\mu$ is exactly known and equal to $0.07$.
3.  $\mu$ has a nominal value of $0.07$ and a standard deviation of $0.02$ and similarily $K_s$ has a nominal value of $10.0$ and a standard deviation of $0.15$.

The last study case is the true case because in reality all encertainties will contribute at once.
"

# ╔═╡ b4c76f21-9839-4db5-b952-43869d3a1efe
md"
#### Study case 1
"

# ╔═╡ 6ab73fc5-7315-4ec6-a3ee-9b5c689ce82e
md"
We initialize a vector `params1_uncert_log` with the parameter values and uncertainty (standard deviation) **only in the first parameter**:
"

# ╔═╡ b0f32f69-b6a2-4050-97e2-8018d07463d8
params1_uncert_log = [:μ => 0.07±0.02, :Wf => 10.0]

# ╔═╡ c1a0702c-6fa5-47d7-bd72-6bfcbdeba93f
md"
We create the corresponding ODE problem and store it in `oprob1_uncert_log`:
"

# ╔═╡ 2fc91add-8468-4402-bcd8-42da5544f612
oprob1_uncert_log = ODEProblem(growth_mod_log, u0_log, tspan, params1_uncert_log)

# ╔═╡ 5fa951d3-23ee-4510-873c-8158df4f2faf
md"
We solve the ODE problem. Use `Tsit5()` and `saveat=2.0`. Store the solution in `osol1_uncert_log`:
"

# ╔═╡ d54635e8-1693-42bc-91f4-7e504a8661c7
osol1_uncert_log = solve(oprob1_uncert_log, Tsit5(), saveat=2.0)

# ╔═╡ 707c6023-49ae-4a13-a220-d827a76a352c
md"
Plot the results (simulation of the output variable $W$ and uncertainty band):
"

# ╔═╡ 48e92497-f6b3-484c-bc5f-3ce69f4f1a91
plot(osol1_uncert_log)

# ╔═╡ 9b30b23b-6e92-49cf-bf93-e12a187fa2ba
md"
When thinking back of the sensitivity of $W$ to the parameter $\mu$, we saw that the corresponding sensitivity function had a maximum around $t=33\;s$. Looking at the above plot with error bars, we can see that de largest error bars (largest uncentrainty in the output variable) occurs at the timepoints where the sensitivity is strongest.
"

# ╔═╡ 86f8d937-c42a-47ff-a0f6-6cae9d637fcf
md"
#### Study case 2
"

# ╔═╡ 0755ea1c-50a0-4c8a-9adc-be2d63948a10
md"
We initialize a vector `params2_uncert_log` with the parameter values and uncertainty (standard deviation) **only in the second parameter**:
"

# ╔═╡ 7befd387-f3a6-4b34-90b0-5f04b5429252
params2_uncert_log = [:μ => 0.07, :Wf => 10.0±1.5]

# ╔═╡ 9e969176-f197-482f-9e2e-8b6d89de8f03
md"
We create the corresponding ODE problem and store it in `oprob2_uncert_log`:
"

# ╔═╡ 50cc4ee1-6601-442c-bac7-4d23b07db87e
oprob2_uncert_log = ODEProblem(growth_mod_log, u0_log, tspan, params2_uncert_log)

# ╔═╡ 6d69f64e-8e3a-46e5-bb2d-bfb6f55e129e
md"
We solve the ODE problem. Use `Tsit5()` and `saveat=2.0`. Store the solution in `osol2_uncert_log`:
"

# ╔═╡ ee548df5-1334-47dc-9258-a25062de3f67
osol2_uncert_log = solve(oprob2_uncert_log, Tsit5(), saveat=2.0)

# ╔═╡ afdb7450-e883-477b-b9a2-8e8585d16c40
md"
Plot the results (simulation of the output variable $W$ and uncertainty band):
"

# ╔═╡ 25bfbf46-69b1-47eb-a473-788aad066d4a
plot(osol2_uncert_log)

# ╔═╡ b8896997-ddec-4d4a-a7c4-0e1008089d1b
md"
When thinking back of the sensitivity of $W$ to the parameter $W_f$, we saw that the corresponding sensitivity function was strongest in the tail of the curve around the steady state value. Looking at the above plot with error bars, we can see that de largest error bars (largest uncentrainty in the output variable) occurs at the tail of the curve where the sensitivity is strongest.
"

# ╔═╡ 5522a5c2-ff7e-4093-86d0-ffb879fcc6e8
md"
#### Study case 3
"

# ╔═╡ 5e3da6ca-656e-41ad-ba04-bf3d6b51939d
md"
We initialize a vector `params_uncert_log` with the parameter values and uncertainty (standard deviation) **in all parameters**:
"

# ╔═╡ 7be29c69-7fa0-4f41-80c9-f69aa01e8ea4
params_uncert_log = [:μ => 0.07±0.02, :Wf => 10.0±1.5]

# ╔═╡ d48c378e-df1c-45f7-9c3a-c5e9640e25a5
md"
We create the corresponding ODE problem and store it in `oprob_uncert_log`:
"

# ╔═╡ f209f4fb-985c-4a9b-b4dd-d50cc3fae1cf
oprob_uncert_log = ODEProblem(growth_mod_log, u0_log, tspan, params_uncert_log)

# ╔═╡ de96646b-b752-4f72-ac78-2a503478bf50
md"
We solve the ODE problem. Use `Tsit5()` and `saveat=2.0`. Store the solution in `osol_uncert_log`:
"

# ╔═╡ 6b1a3c2a-0a88-41ad-8866-41a33da5cd90
osol_uncert_log = solve(oprob_uncert_log, Tsit5(), saveat=2.0)

# ╔═╡ 127ead82-2e2d-4c2e-9417-0490b6839514
md"
Plot the results (simulation of the output variable $W$ and uncertainty band):
"

# ╔═╡ a5c254ad-b629-4a9e-89e9-90460e84620c
plot(osol_uncert_log)

# ╔═╡ c181d973-4fcc-4fbe-b37c-5f71dd33b0d2
md"
Now you have the combined effect of the uncertainty in the parameter $\mu$ as well as in the parameter $W_f$.
"

# ╔═╡ f4406d4e-d9f9-45fe-865c-405dbbb3127f
md"
## Exercises
"

# ╔═╡ e9f49fbc-d877-4db7-b320-a2db964ecaaa
md"
### Exercise 1 - Uncertainty analysis of the exponential growth model

Perform an uncertainty analysis of the exponential growth model. Use the parameter uncertainties mentioned in the *Table* in the *Grass growth models* sections.
"

# ╔═╡ 5201995b-3298-4984-8839-dd25dc73a20f
md"
A possible *reaction network object* for the exponential growth model can be implemented as follows:
"

# ╔═╡ 97bdfce1-a670-469a-b8a9-34e1d1964409
growth_exp = @reaction_network begin
    μ*Wf, 0 --> W
    μ, W --> 0
end

# ╔═╡ 52991cca-0365-48fb-b157-62a71c5573ee
md"
The vector `u₀_exp` with the initial condition is:
"

# ╔═╡ f7241453-755e-466c-a963-a504bc9ea446
u0_exp = [:W => 2.0]

# ╔═╡ 11b81610-90b1-4d71-aeba-0bf33709dbda
md"
Initialize a vector `params_uncert_exp` with the parameter values and their uncertainties (standard deviation) **in all parameters**.\
**Remark**: you can use the same variable and leave a single uncertainty if you want to see the effect of the uncertainty in only one parameter later on.
"

# ╔═╡ 1af01af7-f5ee-4029-94ff-418c2f9748f2
# params_uncert_exp = missing     # Uncomment and complete the instruction
params_uncert_exp = [:μ => 0.02±0.01, :Wf => 10.0±1.5]

# ╔═╡ c25160d9-471a-414c-be99-10aa2efbb6d4
md"
Create the corresponding ODE problem and store it in `oprob_uncert_exp`:
"

# ╔═╡ 7fcec514-192b-4694-b149-14f5f87bae17
# oprob_uncert_exp = missing       # Uncomment and complete the instruction
oprob_uncert_exp = ODEProblem(growth_exp, u0_exp, tspan, params_uncert_exp)

# ╔═╡ bd120dbc-827a-499a-a11c-25423ec5f397
md"
Solve the ODE problem. Use `Tsit5()` and `saveat=2.0`. Store the solution in `osol_uncert_exp`:
"

# ╔═╡ 2dd58666-80ae-4346-ae15-34792eea2d25
# osol_uncert_exp = missing        # Uncomment and complete the instruction
osol_uncert_exp = solve(oprob_uncert_exp, Tsit5(), saveat=2.0)

# ╔═╡ 8059e00b-8e9c-4c11-a484-bea2a0bf155e
md"
Plot the results (simulation of the output variable $W$ and uncertainty band):
"

# ╔═╡ a89677f3-602b-4ef4-8f1a-ab203584f8df
# missing                          # Uncomment and complete the instruction
plot(osol_uncert_exp)

# ╔═╡ 2a470290-bca7-459f-b4f8-a5ac9fc031ed
md"
Draw your conclusions:

- missing
"

# ╔═╡ ded3aa76-ba09-499e-8c2d-52ff81ccfe02
md"
### Exercise 2 - Uncertainty analysis of the exponential Gompertz model

Perform an uncertainty analysis of the Gompertz growth model. Use the parameter uncertainties mentioned in the *Table* in the *Grass growth models* sections.
"

# ╔═╡ e245e5dd-bb11-455b-8de3-ae0ec1763399
md"
A possible *reaction network object* for the Gompertz growth model can be implemented as follows:
"

# ╔═╡ 3e16bd3a-898a-49f9-9ed1-8998be7dc045
growth_gom = @reaction_network begin
#    -μ, W --> ∅
#    D*log(W), W --> ∅
	μ-D*log(W), W --> 2W
end

# ╔═╡ db2f002f-20b7-4b5c-8d59-85b089fb3e59
md"
The vector `u₀_gom` with the initial condition is:
"

# ╔═╡ 22a0013b-04f9-4148-9344-10bee1def0cf
u0_gom = [:W => 2.0]

# ╔═╡ f4fff573-6533-4ac3-897d-840a9aaf55f1
md"
Initialize a vector `params_uncert_gom` with the parameter values and their uncertainties (standard deviation) **in all parameters**.\
**Remark**: you can use the same variable and leave a single uncertainty if you want to see the effect of the uncertainty in only one parameter later on.
"

# ╔═╡ 5ab3aca2-e570-4802-815b-5614680b427d
# params_uncert_gom = missing     # Uncomment and complete the instruction
params_uncert_gom = [:μ => 0.09±0.01, :D => 0.040±0.002]

# ╔═╡ 890fe099-28ae-4e4f-aa1e-8d9ab6240130
md"
Create the corresponding ODE problem and store it in `oprob_uncert_log`:
"

# ╔═╡ 1aed2273-7b3e-4a49-a8d0-91251c58b700
# oprob_uncert_gom = missing      # Uncomment and complete the instruction
oprob_uncert_gom = ODEProblem(growth_gom, u0_gom, tspan, params_uncert_gom)

# ╔═╡ 76da6edd-9e9f-4f5b-8673-b737c8f9c5f8
md"
Solve the ODE problem. Use `Tsit5()` and `saveat=2.0`. Store the solution in `osol_uncert_gom`:
"

# ╔═╡ 8b10af49-e3d0-426e-9147-7666262d1c54
# osol_uncert_gom = missing          # Uncomment and complete the instruction
osol_uncert_gom = solve(oprob_uncert_gom, Tsit5(), saveat=2.0)

# ╔═╡ 6e8318bb-c884-4726-bb3d-9e0cfb7f40c1
md"
Plot the results (simulation of the output variable $W$ and uncertainty band):
"

# ╔═╡ 92785455-14d2-4adb-9731-8575527c5217
# missing                   # Uncomment and complete the instruction
plot(osol_uncert_gom)

# ╔═╡ 9d11df31-57b2-4825-a9db-34b422d5f084
md"
Draw your conclusions:

- missing
"

# ╔═╡ Cell order:
# ╠═d2c4d230-0943-11ef-3aad-5719e74bb20e
# ╠═47790920-74d2-4f21-b0f2-45dc19d25156
# ╠═8349306c-e98c-4221-9a1b-322fde3e18cb
# ╠═83f3d978-aefa-40c1-a647-b26e837aeed6
# ╠═73d5f921-82e4-4160-a3c4-fc1935f4c58d
# ╠═241a8a65-c59f-44f1-be39-d5edd1321b49
# ╠═0f518e36-bc96-4599-807e-728504e5ca7b
# ╠═da3cdc89-b911-4af9-9d4a-526301cba581
# ╠═78afaded-5a19-4386-aa76-7974977ea354
# ╠═eae11742-14c6-4b1b-a939-232710dfa10e
# ╠═e5ab4490-6fd4-4b51-bfa7-1366438efefc
# ╠═3e2750e7-220a-47b2-b445-eb3315734dec
# ╠═c874a08c-7d82-4631-b611-c598dcaded09
# ╠═9d3192bc-ee69-44a1-9273-06b8a55c60ef
# ╠═8acb7ea2-54ef-428b-bb4f-364377d2c38d
# ╠═460d98ef-2177-49bf-87a4-35412f4183ed
# ╠═38e0f746-1d04-46ba-81cc-21522b9ce6a4
# ╠═2538bd6e-027f-47fe-813c-db0a02586d2c
# ╠═72b53ba6-dd05-4f2d-a65e-52783e76a1e9
# ╠═5185d0eb-7392-4775-9336-3e0f9e1449ce
# ╠═1a67db96-6973-446d-aa5b-253dd0008a73
# ╠═048112b0-10f1-47c5-834f-c3851d3078e8
# ╠═037cbf2b-dc0c-4b62-879b-aa702c964877
# ╠═40c1192c-bc76-4d3a-8144-dd76ba72fba2
# ╠═0b5493eb-6ce6-43c1-9289-4c2e6f20d014
# ╠═b4c76f21-9839-4db5-b952-43869d3a1efe
# ╠═6ab73fc5-7315-4ec6-a3ee-9b5c689ce82e
# ╠═b0f32f69-b6a2-4050-97e2-8018d07463d8
# ╠═c1a0702c-6fa5-47d7-bd72-6bfcbdeba93f
# ╠═2fc91add-8468-4402-bcd8-42da5544f612
# ╠═5fa951d3-23ee-4510-873c-8158df4f2faf
# ╠═d54635e8-1693-42bc-91f4-7e504a8661c7
# ╠═707c6023-49ae-4a13-a220-d827a76a352c
# ╠═48e92497-f6b3-484c-bc5f-3ce69f4f1a91
# ╠═9b30b23b-6e92-49cf-bf93-e12a187fa2ba
# ╠═86f8d937-c42a-47ff-a0f6-6cae9d637fcf
# ╠═0755ea1c-50a0-4c8a-9adc-be2d63948a10
# ╠═7befd387-f3a6-4b34-90b0-5f04b5429252
# ╠═9e969176-f197-482f-9e2e-8b6d89de8f03
# ╠═50cc4ee1-6601-442c-bac7-4d23b07db87e
# ╠═6d69f64e-8e3a-46e5-bb2d-bfb6f55e129e
# ╠═ee548df5-1334-47dc-9258-a25062de3f67
# ╠═afdb7450-e883-477b-b9a2-8e8585d16c40
# ╠═25bfbf46-69b1-47eb-a473-788aad066d4a
# ╠═b8896997-ddec-4d4a-a7c4-0e1008089d1b
# ╠═5522a5c2-ff7e-4093-86d0-ffb879fcc6e8
# ╠═5e3da6ca-656e-41ad-ba04-bf3d6b51939d
# ╠═7be29c69-7fa0-4f41-80c9-f69aa01e8ea4
# ╠═d48c378e-df1c-45f7-9c3a-c5e9640e25a5
# ╠═f209f4fb-985c-4a9b-b4dd-d50cc3fae1cf
# ╠═de96646b-b752-4f72-ac78-2a503478bf50
# ╠═6b1a3c2a-0a88-41ad-8866-41a33da5cd90
# ╠═127ead82-2e2d-4c2e-9417-0490b6839514
# ╠═a5c254ad-b629-4a9e-89e9-90460e84620c
# ╠═c181d973-4fcc-4fbe-b37c-5f71dd33b0d2
# ╠═f4406d4e-d9f9-45fe-865c-405dbbb3127f
# ╠═e9f49fbc-d877-4db7-b320-a2db964ecaaa
# ╠═5201995b-3298-4984-8839-dd25dc73a20f
# ╠═97bdfce1-a670-469a-b8a9-34e1d1964409
# ╠═52991cca-0365-48fb-b157-62a71c5573ee
# ╠═f7241453-755e-466c-a963-a504bc9ea446
# ╠═11b81610-90b1-4d71-aeba-0bf33709dbda
# ╠═1af01af7-f5ee-4029-94ff-418c2f9748f2
# ╠═c25160d9-471a-414c-be99-10aa2efbb6d4
# ╠═7fcec514-192b-4694-b149-14f5f87bae17
# ╠═bd120dbc-827a-499a-a11c-25423ec5f397
# ╠═2dd58666-80ae-4346-ae15-34792eea2d25
# ╠═8059e00b-8e9c-4c11-a484-bea2a0bf155e
# ╠═a89677f3-602b-4ef4-8f1a-ab203584f8df
# ╠═2a470290-bca7-459f-b4f8-a5ac9fc031ed
# ╠═ded3aa76-ba09-499e-8c2d-52ff81ccfe02
# ╠═e245e5dd-bb11-455b-8de3-ae0ec1763399
# ╠═3e16bd3a-898a-49f9-9ed1-8998be7dc045
# ╠═db2f002f-20b7-4b5c-8d59-85b089fb3e59
# ╠═22a0013b-04f9-4148-9344-10bee1def0cf
# ╠═f4fff573-6533-4ac3-897d-840a9aaf55f1
# ╠═5ab3aca2-e570-4802-815b-5614680b427d
# ╠═890fe099-28ae-4e4f-aa1e-8d9ab6240130
# ╠═1aed2273-7b3e-4a49-a8d0-91251c58b700
# ╠═76da6edd-9e9f-4f5b-8673-b737c8f9c5f8
# ╠═8b10af49-e3d0-426e-9147-7666262d1c54
# ╠═6e8318bb-c884-4726-bb3d-9e0cfb7f40c1
# ╠═92785455-14d2-4adb-9731-8575527c5217
# ╠═9d11df31-57b2-4825-a9db-34b422d5f084
