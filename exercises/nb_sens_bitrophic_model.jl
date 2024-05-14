### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ a3c1055e-b44f-4412-9ca0-f8ec5f972494
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ fd357fe0-0920-11ef-02e4-25a84575c6a2
using Markdown

# ╔═╡ 0ddc1082-bf8d-4357-a056-2fa861824b47
using InteractiveUtils

# ╔═╡ cce8fcde-069f-49b9-aadd-6abccbcfbace
using Catalyst

# ╔═╡ b48e4b46-3a07-4ff0-b617-6bc1125f6a1e
using DifferentialEquations, Plots

# ╔═╡ ddc035a2-3b03-4a54-807c-05a271384ca7
using SciMLSensitivity

# ╔═╡ ceec287b-5622-48d2-9fe1-26b3e84e3f45
md"
### Exercise: Bitrophic model - Sensitivity analysis
"

# ╔═╡ 2a732f39-0499-4c90-be62-a1e2278e0c8b
md"
The dynamic relationship between a field crop and a voracious insect population within an ecosystem can be represented by a bitrophic model. Such model typically consists of two variables: the abundance of the field crop, often representing a primary producer such as a plant species, and the population size of the voracious insect, which acts as a consumer feeding on the crop. The differential equations below describe how changes in the crop population affect the growth and behavior of the insect population, and vice versa, under the influence of an insecticide.

$$\begin{eqnarray*}
\frac{dC}{dt} &= \theta C \left(1-\frac{C}{k}\right)-fCA \\
\frac{dA}{dt} &= \phi f CA -(1 + p)\, \mu A
\end{eqnarray*}$$

Understanding this bitrophic interaction is crucial for predicting the impact of insect predation on crop yields and devising effective strategies for pest management in agriculture and ecological conservation efforts.

In these equations, $C$ and $A$ are both expressed in $kg/ha$, $\theta=0.2\;d^{-1}$, $k=4000\;kg/ha$, $f=0.001\;ha/(kg\,d)$, the efficiency ratio $\phi=0.2$, and the mortality ratio $\mu=0.1\;d^{-1}$. The crop can be treated with an insecticide which increases the insect's death coefficient by a factor of $p=3$. The factor $p$ depends on the applied insecticide concentration and can therefore be controlled externally. 

At the beginning of a season, $100\;kg$ of the crop and $0.5\;kg$ of insects per $ha$ are present. 
"

# ╔═╡ 87508f87-42d7-44e7-a0fb-415e6d7190c2
md"
Set up a *reaction network* model by analysing the terms in the above differential equations and simulate the evolution of $C$ and $A$ for $200$ days. Next, perform a sensitivity analysis of $C$ and $A$ to the parameters $\theta$, $\phi$ and $p$.
"

# ╔═╡ eeb52e5e-f57b-4746-bbfc-7f043fc50e96
md"
Set up a *reaction network* model and name it `bitrophic_model`.\
Hints:
-  $C$ is being created at a rate $\theta C \left(1-\frac{C}{k}\right)$
-  $C$ is being converted into a fraction $\phi$ of $A$ at a rate $f A$
-  $A$ is being degraded at a rate $(1 + p)\, \mu A$
"

# ╔═╡ 0b902d89-f0c2-4621-998e-7013931315ae
# bitrophic_model = @reaction_network begin
#     missing        # Uncomment and complete the instruction
#     missing        # Uncomment and complete the instruction
#     missing        # Uncomment and complete the instruction
# end
bitrophic_model = @reaction_network begin
    θ*(1-C/k)*C, ∅ --> C
    f*A, C --> ϕ*A
    (1+p)*μ, A --> ∅
end

# ╔═╡ 92027dcf-9765-4fad-a0a4-f19967ce6bce
md"
Check out the order of the parameters. This will be important later for selecting the right sensitivity functions.
"

# ╔═╡ 051f1b87-34ec-4e3f-8d4e-c4e42430be8f
# missing             # Uncomment and complete the instruction
parameters(bitrophic_model)

# ╔═╡ f11bae81-a855-48ec-90a7-16b4ea1fc30b
md"
Convert the system to a symbolic differential equations model, name it `osys` and verify, by analyzing the differential equations, that your model is correctly implemented.
"

# ╔═╡ 1a19c02c-5d0f-4271-b0e7-7da873030d7e
# osys = missing       # Uncomment and complete the instruction
osys = convert(ODESystem, bitrophic_model)

# ╔═╡ 55eba88f-d4cf-4ccf-b061-51feb3a4290c
md"
Initialize a vector `u₀` with the initial conditions, define the timespan in `tspan` and initialize a vector `param` with the parameter values:
"

# ╔═╡ 7769c7a5-5f19-4e1a-8721-41ee03c2024e
# u₀ = missing         # Uncomment and complete the instruction
u₀ = [:C => 100, :A => 0.5]

# ╔═╡ 00453a3c-88d1-4cc4-a8f2-d4fa2eed0b66
# tspan = missing      # Uncomment and complete the instruction
tspan = (0.0, 200.0)

# ╔═╡ 157eae84-9c8a-4d23-9ada-9d53ba520b0a
# params = missing     # Uncomment and complete the instruction
params = [:θ => 0.2, :k => 4000.0, :f => 0.001, :ϕ => 0.2, :p => 3, :μ => 0.1]

# ╔═╡ eddc01fa-3ab0-4fd0-963f-5704f48778a7
md"
Create the ODE problem and store it in `oprob`. Next, solve the ODE problem using `Tsit5()` and `saveat=0.5`, and store the solution in `osol`. Finally plot the results.
"

# ╔═╡ 1c185585-7c7a-4073-85f5-8300d3fe19e4
# oprob = missing       # Uncomment and complete the instruction
oprob = ODEProblem(bitrophic_model, u₀, tspan, params)

# ╔═╡ 2f1b4426-7235-4caf-84ed-ec6e04eaa34a
# osol = missing        # Uncomment and complete the instruction
osol = solve(oprob, Tsit5(), saveat=0.5)

# ╔═╡ 63b33638-b967-43a7-b6ff-f6912f9a12e3
# missing               # Uncomment and complete the instruction
plot(osol)

# ╔═╡ 81c6198f-5f69-4d78-b67c-b6e283a6c3b1
md"
Interpret the simulations results. Ask yourself the following questions:

1. What happens to $C$ and $A$ during the first 30 days?
    - Answer: missing
2. Why does $C$ starts to decline around day 40?
    - Answer: missing
3. What happen to $C$ and $A$ from day 50 on, do they finally reach steady state values?
    - Answer: missing
"

# ╔═╡ 18578789-f389-4d6e-a26b-9c8776f3d1a8
md"
Create the ODE forward sensitivity problem and store it in `oprob_sens`:
"

# ╔═╡ d17247da-0aeb-4f18-a9f1-0b8255e2be9c
# oprob_sens = missing         # Uncomment and complete the instruction
oprob_sens = ODEForwardSensitivityProblem(oprob.f, [100, 0.5], tspan, [0.2, 4000.0, 0.001, 0.2, 3, 0.1])

# ╔═╡ 3254b7a9-bf50-4d34-ae64-0ac79035db4b
md"
Solve the *ODE forward sensitivity problem* using `Tsit5()` and `saveat=0.5`, and store the solution in `osol_sens`:
"

# ╔═╡ c31fc27c-ca6f-45b7-9348-b16e05c44dfa
# osol_sens = missing           # Uncomment and complete the instruction
osol_sens = solve(oprob_sens, Tsit5(), saveat=0.5)

# ╔═╡ 45c00c7c-c57f-4727-8d5c-af599a1c5674
md"
Extract the sensitivity functions. Store the simulation results in the variable `u` and the sensitivities in the variable `dp`:
"

# ╔═╡ 2afbc8b2-e49e-40ca-ac14-d01d8d5b32bb
# u, dp = missing                 # Uncomment and complete the instruction
u, dp = extract_local_sensitivities(osol_sens)

# ╔═╡ 8b39af5a-2d86-48bd-abf0-8fd00a371c83
md"
Select by indexing and assigment the sensitivities of $C$ and $A$ to $\theta$, $\phi$ and $p$, to the variables `sens_theta`, `sens_phi` and `sens_p` respectively. Don't forget to transpose where necessary.

**Remark:**
- each of the elements `pd[i]` will contain two sensitivity functions (one of $C$ and one of $A$) to the `i`-th parameter. You can find the parameter indices of $\theta$, $\phi$ and $p$ by calling the function `parameters` in conjunction with the *reaction network* model name.
"

# ╔═╡ 6397e7a3-6265-4e84-9025-4617736c0ded
# Select sensitivity of C and A to θ
# sens_theta = missing             # Uncomment and complete the instruction
sens_theta = dp[1]'

# ╔═╡ ca36a17c-978a-43f1-ab19-06785ed6963f
# Select sensitivity of C and A to ϕ
# sens_phi   = missing             # Uncomment and complete the instruction
sens_phi   = dp[4]'

# ╔═╡ 17726224-7bf7-48f3-8763-e636de3b6775
# Select sensitivity of C and A to p
# sens_p     = missing             # Uncomment and complete the instruction
sens_p     = dp[5]'

# ╔═╡ 3814640e-ed0b-48a1-b4ac-6b7f04d899f3
md"
Plot both sensitivity functions in separate notebook cells (with appropriate title and label):
"

# ╔═╡ cda325f3-7086-4503-94d7-a954a0cbd90d
# Plot sensitivity of C and A to θ
# missing
plot(osol_sens.t, sens_theta, title="Absolute sensitivity of C and A to θ", label=["dC/dθ" "dA/dθ"])

# ╔═╡ 83903533-f8c1-4217-a930-4fddb885ff61
# Plot sensitivity of C and A to ϕ
# missing
plot(osol_sens.t, sens_phi, title="Absolute sensitivity of C and A to ϕ", label=["dC/dϕ" "dA/dϕ"])

# ╔═╡ 3eeb967d-ea9b-44da-af4c-d829a861c803
# Plot sensitivity of C and A to p
# missing
plot(osol_sens.t, sens_p, title="Absolute sensitivity of C and A to p", label=["dC/dp" "dA/dp"])

# ╔═╡ 3a4be902-38ec-412f-9d25-1472f1371b5f
md"
Analyse the sensitivity functions. Ask yourself the following questions:

1. In steady state, does $\phi$ have a positive or negative effect on $C$? Explain why this could be.
    - Answer: missing
2. In steady state, does $p$ have a positive or negative effect on $C$? Explain why this could be.
    - Answer: missing
"

# ╔═╡ Cell order:
# ╠═fd357fe0-0920-11ef-02e4-25a84575c6a2
# ╠═0ddc1082-bf8d-4357-a056-2fa861824b47
# ╠═a3c1055e-b44f-4412-9ca0-f8ec5f972494
# ╠═cce8fcde-069f-49b9-aadd-6abccbcfbace
# ╠═b48e4b46-3a07-4ff0-b617-6bc1125f6a1e
# ╠═ddc035a2-3b03-4a54-807c-05a271384ca7
# ╠═ceec287b-5622-48d2-9fe1-26b3e84e3f45
# ╠═2a732f39-0499-4c90-be62-a1e2278e0c8b
# ╠═87508f87-42d7-44e7-a0fb-415e6d7190c2
# ╠═eeb52e5e-f57b-4746-bbfc-7f043fc50e96
# ╠═0b902d89-f0c2-4621-998e-7013931315ae
# ╠═92027dcf-9765-4fad-a0a4-f19967ce6bce
# ╠═051f1b87-34ec-4e3f-8d4e-c4e42430be8f
# ╠═f11bae81-a855-48ec-90a7-16b4ea1fc30b
# ╠═1a19c02c-5d0f-4271-b0e7-7da873030d7e
# ╠═55eba88f-d4cf-4ccf-b061-51feb3a4290c
# ╠═7769c7a5-5f19-4e1a-8721-41ee03c2024e
# ╠═00453a3c-88d1-4cc4-a8f2-d4fa2eed0b66
# ╠═157eae84-9c8a-4d23-9ada-9d53ba520b0a
# ╠═eddc01fa-3ab0-4fd0-963f-5704f48778a7
# ╠═1c185585-7c7a-4073-85f5-8300d3fe19e4
# ╠═2f1b4426-7235-4caf-84ed-ec6e04eaa34a
# ╠═63b33638-b967-43a7-b6ff-f6912f9a12e3
# ╠═81c6198f-5f69-4d78-b67c-b6e283a6c3b1
# ╠═18578789-f389-4d6e-a26b-9c8776f3d1a8
# ╠═d17247da-0aeb-4f18-a9f1-0b8255e2be9c
# ╠═3254b7a9-bf50-4d34-ae64-0ac79035db4b
# ╠═c31fc27c-ca6f-45b7-9348-b16e05c44dfa
# ╠═45c00c7c-c57f-4727-8d5c-af599a1c5674
# ╠═2afbc8b2-e49e-40ca-ac14-d01d8d5b32bb
# ╠═8b39af5a-2d86-48bd-abf0-8fd00a371c83
# ╠═6397e7a3-6265-4e84-9025-4617736c0ded
# ╠═ca36a17c-978a-43f1-ab19-06785ed6963f
# ╠═17726224-7bf7-48f3-8763-e636de3b6775
# ╠═3814640e-ed0b-48a1-b4ac-6b7f04d899f3
# ╠═cda325f3-7086-4503-94d7-a954a0cbd90d
# ╠═83903533-f8c1-4217-a930-4fddb885ff61
# ╠═3eeb967d-ea9b-44da-af4c-d829a861c803
# ╠═3a4be902-38ec-412f-9d25-1472f1371b5f
