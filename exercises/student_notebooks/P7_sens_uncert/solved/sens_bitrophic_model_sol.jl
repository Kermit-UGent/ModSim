### A Pluto.jl notebook ###
# v0.19.46

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
using ForwardDiff

# ╔═╡ ceec287b-5622-48d2-9fe1-26b3e84e3f45
md"""
### Exercise: Bitrophic model - Sensitivity analysis
"""

# ╔═╡ 2a732f39-0499-4c90-be62-a1e2278e0c8b
md"""
The dynamic relationship between a field crop and a voracious insect population within an ecosystem can be represented by a bitrophic model. Such model typically consists of two variables: the abundance of the field crop, often representing a primary producer such as a plant species, and the population size of the voracious insect, which acts as a consumer feeding on the crop. The differential equations below describe how changes in the crop population affect the growth and behavior of the insect population, and vice versa, under the influence of an insecticide.

$$\begin{eqnarray*}
\frac{dC}{dt} &= \theta C \left(1-\frac{C}{k}\right)-fCA \\
\frac{dA}{dt} &= \phi f CA -(1 + p)\, \mu A
\end{eqnarray*}$$

Understanding this bitrophic interaction is crucial for predicting the impact of insect predation on crop yields and devising effective strategies for pest management in agriculture and ecological conservation efforts.

In these equations, $C$ and $A$ are both expressed in $kg/ha$, $\theta=0.2\;d^{-1}$, $k=4000\;kg/ha$, $f=0.001\;ha/(kg\,d)$, the efficiency ratio $\phi=0.2$, and the mortality ratio $\mu=0.1\;d^{-1}$. The crop can be treated with an insecticide which increases the insect's death coefficient by a factor of $p=3$. The factor $p$ depends on the applied insecticide concentration and can therefore be controlled externally. 

At the beginning of a season, $100\;kg$ of the crop and $0.5\;kg$ of insects per $ha$ are present. 
"""

# ╔═╡ 87508f87-42d7-44e7-a0fb-415e6d7190c2
md"""
Set up a *reaction network* model by analysing the terms in the above differential equations and simulate the evolution of $C$ and $A$ for $200$ days. Next, perform a sensitivity analysis of $C$ and $A$ to the parameters $\theta$, $\phi$ and $p$.
"""

# ╔═╡ eeb52e5e-f57b-4746-bbfc-7f043fc50e96
md"""
Set up a *reaction network* model and name it `bitrophic_model`.\
Hints:
-  $C$ is growing (i.e., $C \rightarrow 2C$) at a rate $\theta \left(1-\frac{C}{k}\right)$.
-  The insects $A$ eat crops $C$ (i.e., $C + A$) at a rate $f$ resulting in an increase of a factor of $\phi$ more insects (i.e., $(1+\phi)A$).
-  The insects $A$ are dying (i.e., $A \rightarrow 0$) at a rate $(1 + p)\, \mu$.
"""

# ╔═╡ 0b902d89-f0c2-4621-998e-7013931315ae
# bitrophic_model = @reaction_network begin
#     missing        # Uncomment and complete the instruction
#     missing        # Uncomment and complete the instruction
#     missing        # Uncomment and complete the instruction
# end
bitrophic_model = @reaction_network begin
	θ*(1-C/k), C --> 2C
	f, C+A --> (1+ϕ)*A
    (1+p)*μ, A --> 0
end

# ╔═╡ 92027dcf-9765-4fad-a0a4-f19967ce6bce
md"""
Check out the species and the parameters.
"""

# ╔═╡ 15dbf491-c5ae-4682-a21d-1a815a957993
# missing             # Uncomment and complete the instruction
species(bitrophic_model)

# ╔═╡ 051f1b87-34ec-4e3f-8d4e-c4e42430be8f
# missing             # Uncomment and complete the instruction
parameters(bitrophic_model)

# ╔═╡ f11bae81-a855-48ec-90a7-16b4ea1fc30b
md"""
Convert the system to a symbolic differential equations model, name it `osys` and verify, by analyzing the differential equations, that your model is correctly implemented.
"""

# ╔═╡ 1a19c02c-5d0f-4271-b0e7-7da873030d7e
# osys = missing       # Uncomment and complete the instruction
osys = convert(ODESystem, bitrophic_model)

# ╔═╡ 55eba88f-d4cf-4ccf-b061-51feb3a4290c
md"""
Initialize a vector `u0` with the initial conditions, define the timespan in `tspan` and initialize a vector `param` with the parameter values:
"""

# ╔═╡ 7769c7a5-5f19-4e1a-8721-41ee03c2024e
# u0 = missing         # Uncomment and complete the instruction
u0 = [:C => 100, :A => 0.5]

# ╔═╡ 00453a3c-88d1-4cc4-a8f2-d4fa2eed0b66
# tspan = missing      # Uncomment and complete the instruction
tspan = (0.0, 200.0)

# ╔═╡ 4e1acbdd-65ce-4f5e-b557-633af8a11aa9
md"""
For clarity, we will use the variables `θ`, `ϕ` and `p` to store the parameter values that are used for the calculation of the sensitivity functions.
"""

# ╔═╡ 9fda8191-e28d-4d10-b330-c386711a2a73
θ = 0.2

# ╔═╡ ae172c81-539c-4482-9cd2-ab74fd83c26a
ϕ = 0.2

# ╔═╡ 28417aea-9c05-4ba0-8b85-b8cccd808238
p = 3

# ╔═╡ 157eae84-9c8a-4d23-9ada-9d53ba520b0a
# params = missing     # Uncomment and complete the instruction
params = [:θ => θ, :k => 4000.0, :f => 0.001, :ϕ => ϕ, :p => p, :μ => 0.1]

# ╔═╡ eddc01fa-3ab0-4fd0-963f-5704f48778a7
md"""
Create the ODE problem and store it in `oprob`. Next, solve the ODE problem using `Tsit5()` and `saveat=0.5`, and store the solution in `osol`. Finally plot the results.
"""

# ╔═╡ 1c185585-7c7a-4073-85f5-8300d3fe19e4
# oprob = missing       # Uncomment and complete the instruction
oprob = ODEProblem(bitrophic_model, u0, tspan, params)

# ╔═╡ 2f1b4426-7235-4caf-84ed-ec6e04eaa34a
# osol = missing        # Uncomment and complete the instruction
osol = solve(oprob, Tsit5(), saveat=0.5)

# ╔═╡ 63b33638-b967-43a7-b6ff-f6912f9a12e3
# missing               # Uncomment and complete the instruction
plot(osol)

# ╔═╡ 81c6198f-5f69-4d78-b67c-b6e283a6c3b1
md"""
Interpret the simulations results. Ask yourself the following questions:

1. What happens to $C$ and $A$ during the first 30 days?
    - Answer: missing
2. Why does $C$ starts to decline around day 40?
    - Answer: missing
3. What happen to $C$ and $A$ from day 50 on, do they finally reach steady state values?
    - Answer: missing
"""
#=
1. C is increasing strongly while A remains relatively weak and constant.
2. Around t=40, A starts to increase (use the option ylim=(0, 300) in plot to see this). As A increase, then C will decline.
3. They will somehow start to oscillate around their steady state values. They finally reach steady state values.
=#

# ╔═╡ aa07ac11-1d6e-41ea-885d-1a6ae38bd58d
md"""
Write a solution function with as argument a vector of the parameters (that you want the sensitivity on), and that returns the outputs.
"""

# ╔═╡ 2a858811-9ffd-4dc7-84ba-f7a59421753f
# Uncomment and complete the instruction
# function bitrophic_model_sim(params)
# 	missing
# 	...
# end
function bitrophic_model_sim(params)
	θ, ϕ, p = params
    u0 = [:C => 100, :A => 0.5]
    tspan = (0.0, 200.0)
    params = [:θ => θ, :k => 4000.0, :f => 0.001, :ϕ => ϕ, :p => p, :μ => 0.1]
    oprob = ODEProblem(bitrophic_model, u0, tspan, params)
	osol = solve(oprob, Tsit5(), saveat=0.5)
	return osol
end

# ╔═╡ e911dc09-c56a-42ed-9a09-9c4bba1aa5f6
md"""
Make two functions based on the solution function that each returns a single output, hence, one function that returns the output $C$, and another function that returns the output $A$.
"""

# ╔═╡ 2e43110a-d948-4e6f-9ed9-a35101ba223a
# bitrophic_model_sim_C(params) = missing # Uncomment and complete the instruction
bitrophic_model_sim_C(params) = bitrophic_model_sim(params)[:C]

# ╔═╡ 298604da-3207-40db-8535-f628d17a0f39
# bitrophic_model_sim_A(params) = missing # Uncomment and complete the instruction
bitrophic_model_sim_A(params) = bitrophic_model_sim(params)[:A]

# ╔═╡ 76b256dc-bdec-4fb8-95c7-b32df03ba8b1
md"""
Make the time vector.
"""

# ╔═╡ 3d1f2c90-2e4e-4fa1-af41-9d44d00e5090
# t_vals = missing         # Uncomment and complete the instruction
t_vals = 0:0.5:200.0

# ╔═╡ c2ca4b68-39a6-45f7-a2a1-0f5c4903ed07
md"""
Compute the two outputs $C$ and $A$ for the given parameter values.
"""

# ╔═╡ 7bc637be-8238-43f4-b190-b46833fc9d34
# C_sim =  missing         # Uncomment and complete the instruction
C_sim = bitrophic_model_sim_C([θ, ϕ, p])

# ╔═╡ 127d8761-97d6-4fcb-b633-98ee30ebc617
# A_sim = missing           # Uncomment and complete the instruction
A_sim = bitrophic_model_sim_A([θ, ϕ, p])

# ╔═╡ cda34ccc-e018-4e2c-a4d3-8ed3d0c1624b
md"""
Using `ForwardDiff.jacobian` to compute the sensitivities for the single ouputs $C$ and $A$. Hence, you need to call `ForwardDiff.jacobian` twice.
"""

# ╔═╡ 55f00963-bfb0-420f-ae6f-d20f93f01b1e
# sens_C = missing           # Uncomment and complete the instruction
sens_C = ForwardDiff.jacobian(bitrophic_model_sim_C, [θ, ϕ, p])

# ╔═╡ 4d2c77f2-1a5f-437c-adb6-232155f2cafb
# sens_A = missing            # Uncomment and complete the instruction
sens_A = ForwardDiff.jacobian(bitrophic_model_sim_A, [θ, ϕ, p])

# ╔═╡ c261b580-2cc7-4aba-af15-1778f9293278
md"""
Extract the (absolute) sensitivities of the outputs on the different parameters.
"""

# ╔═╡ 730de2d2-104b-4f4a-96da-4f96b0445f95
# Uncomment and complete the instruction
# begin
# 	sens_C_on_θ = missing
# 	sens_C_on_ϕ = missing
# 	sens_C_on_p = missing
# end
begin
	sens_C_on_θ = sens_C[:, 1]
	sens_C_on_ϕ = sens_C[:, 2]
	sens_C_on_p = sens_C[:, 3]
end

# ╔═╡ 2d8df796-c362-49d9-b938-26f1917e800b
# Uncomment and complete the instruction
# begin
# 	sens_A_on_θ = missing
# 	sens_A_on_ϕ = missing
# 	sens_A_on_p = missing
# end
begin
	sens_A_on_θ = sens_A[:, 1]
	sens_A_on_ϕ = sens_A[:, 2]
	sens_A_on_p = sens_A[:, 3]
end

# ╔═╡ 6b3ee543-e3e2-4fe6-925d-1f201243d014
md"""
Compute the normalized sensitivities.
"""

# ╔═╡ 89830894-4481-480b-954d-b24971590eaf
# Uncomment and complete the instruction
# begin
# 	sens_C_on_θ_rel = missing
# 	sens_C_on_ϕ_rel = missing
# 	sens_C_on_p_rel = missing
# end
begin
	sens_C_on_θ_rel = sens_C_on_θ .* θ ./ C_sim
	sens_C_on_ϕ_rel = sens_C_on_ϕ .* ϕ ./ C_sim
	sens_C_on_p_rel = sens_C_on_p .* p ./ C_sim
end

# ╔═╡ c2b5c2c4-3aa6-466d-b106-d0d7aef43bdf
# Uncomment and complete the instruction
# begin
# 	sens_A_on_θ_rel = missing
# 	sens_A_on_ϕ_rel = missing
# 	sens_A_on_p_rel = missing
# end
begin
	sens_A_on_θ_rel = sens_A_on_θ .* θ ./ A_sim
	sens_A_on_ϕ_rel = sens_A_on_ϕ .* ϕ ./ A_sim
	sens_A_on_p_rel = sens_A_on_p .* p ./ A_sim
end

# ╔═╡ 3814640e-ed0b-48a1-b4ac-6b7f04d899f3
md"""
Plot the sensitivity functions of $C$ and $A$ on $\theta$. Provide a suitable title (`title="..."`), labels (`label=["..." "..."]`) and an x-label (`xlabel="..."`).
"""

# ╔═╡ bbe8980e-441d-4ac6-abc5-88ff77de9f24
# missing       # Uncomment and complete the instruction
plot(t_vals, [sens_C_on_θ_rel, sens_A_on_θ_rel], title="Normalized sensitivities", label=["C on θ" "A on θ"], xlabel="Time (days)")

# ╔═╡ a8257414-2d6a-4913-9f05-6369c62189e8
md"""
Interpret your results. Try to answer the following question(s):

- In steady state, does $\theta$ have any influence on $C$? Explain why this could be.
    - Answer: missing
- In steady state, why does $\theta$ have a positive effect on $A$? Explain why this could be.
    - Answer: missing
"""
#=
Answer for both questions: theta seems to have no influence on C in steady state. This could be because if theta increases, then the growth rate of C increases resulting in more crops, but then this will also result in more voracious insects (positive effect on A from θ) which will feed on the increased number of crops resulting in cancelling out the increase of the crops.
If you calculated analytically the steady state for C and A (call them Cw and Aw), then you will get:
Cw = (1+p)*μ/(ϕ*f) and Aw = (θ/f)*(1-Cw/k)
you can see that Cw is independend on θ, while Aw is proportional to θ
=#

# ╔═╡ 882f2f6b-eb58-4c8f-877e-76e460900c0c
md"""
Plot the sensitivity functions of $C$ and $A$ on $\phi$. Provide a suitable title (`title="..."`), labels (`label=["..." "..."]`) and an x-label (`xlabel="..."`).
"""

# ╔═╡ f76e5c58-12fc-4b2d-8718-e15144fbef98
# missing        # Uncomment and complete the instruction
plot(t_vals, [sens_C_on_ϕ_rel, sens_A_on_ϕ_rel], title="Normalized sensitivities", label=["C on ϕ" "A on ϕ"], xlabel="Time (days)")

# ╔═╡ 0bd84cf5-c0a2-4bb2-b5f4-f16da6eaab04
md"""
Plot the sensitivity functions of $C$ and $A$ on $p$. Provide a suitable title (`title="..."`), labels (`label=["..." "..."]`) and an x-label (`xlabel="..."`).
"""

# ╔═╡ cb02ebf8-aee5-47a2-9496-6f63fc123147
# missing            # Uncomment and complete the instruction
plot(t_vals, [sens_C_on_p_rel, sens_A_on_p_rel], title="Normalized sensitivities", label=["C on p" "A on p"], xlabel="Time (days)")

# ╔═╡ 869a9dc3-8227-4c6e-8901-198cafb489e2
md"""
Interpret your results. Try to answer the following question(s):

1. In steady state, does $\phi$ have a positive or negative effect on $C$? Explain why this could be.
    - Answer: missing
2. In steady state, does $p$ have a positive or negative effect on $C$? Explain why this could be.
    - Answer: missing
"""
#=
1. ϕ has a negative influence on C This makes sense because ϕ is the yield factor for the varocious insects, the larger ϕ, the more insects and hence less crops.
2. The factor p contributes to the amount of pesticides, hence if p is increased, A will decrease and hence C will increase.
=#

# ╔═╡ Cell order:
# ╠═fd357fe0-0920-11ef-02e4-25a84575c6a2
# ╠═0ddc1082-bf8d-4357-a056-2fa861824b47
# ╠═a3c1055e-b44f-4412-9ca0-f8ec5f972494
# ╠═cce8fcde-069f-49b9-aadd-6abccbcfbace
# ╠═b48e4b46-3a07-4ff0-b617-6bc1125f6a1e
# ╠═ddc035a2-3b03-4a54-807c-05a271384ca7
# ╟─ceec287b-5622-48d2-9fe1-26b3e84e3f45
# ╟─2a732f39-0499-4c90-be62-a1e2278e0c8b
# ╟─87508f87-42d7-44e7-a0fb-415e6d7190c2
# ╟─eeb52e5e-f57b-4746-bbfc-7f043fc50e96
# ╠═0b902d89-f0c2-4621-998e-7013931315ae
# ╟─92027dcf-9765-4fad-a0a4-f19967ce6bce
# ╠═15dbf491-c5ae-4682-a21d-1a815a957993
# ╠═051f1b87-34ec-4e3f-8d4e-c4e42430be8f
# ╟─f11bae81-a855-48ec-90a7-16b4ea1fc30b
# ╠═1a19c02c-5d0f-4271-b0e7-7da873030d7e
# ╟─55eba88f-d4cf-4ccf-b061-51feb3a4290c
# ╠═7769c7a5-5f19-4e1a-8721-41ee03c2024e
# ╠═00453a3c-88d1-4cc4-a8f2-d4fa2eed0b66
# ╟─4e1acbdd-65ce-4f5e-b557-633af8a11aa9
# ╠═9fda8191-e28d-4d10-b330-c386711a2a73
# ╠═ae172c81-539c-4482-9cd2-ab74fd83c26a
# ╠═28417aea-9c05-4ba0-8b85-b8cccd808238
# ╠═157eae84-9c8a-4d23-9ada-9d53ba520b0a
# ╟─eddc01fa-3ab0-4fd0-963f-5704f48778a7
# ╠═1c185585-7c7a-4073-85f5-8300d3fe19e4
# ╠═2f1b4426-7235-4caf-84ed-ec6e04eaa34a
# ╠═63b33638-b967-43a7-b6ff-f6912f9a12e3
# ╟─81c6198f-5f69-4d78-b67c-b6e283a6c3b1
# ╟─aa07ac11-1d6e-41ea-885d-1a6ae38bd58d
# ╠═2a858811-9ffd-4dc7-84ba-f7a59421753f
# ╟─e911dc09-c56a-42ed-9a09-9c4bba1aa5f6
# ╠═2e43110a-d948-4e6f-9ed9-a35101ba223a
# ╠═298604da-3207-40db-8535-f628d17a0f39
# ╟─76b256dc-bdec-4fb8-95c7-b32df03ba8b1
# ╠═3d1f2c90-2e4e-4fa1-af41-9d44d00e5090
# ╟─c2ca4b68-39a6-45f7-a2a1-0f5c4903ed07
# ╠═7bc637be-8238-43f4-b190-b46833fc9d34
# ╠═127d8761-97d6-4fcb-b633-98ee30ebc617
# ╟─cda34ccc-e018-4e2c-a4d3-8ed3d0c1624b
# ╠═55f00963-bfb0-420f-ae6f-d20f93f01b1e
# ╠═4d2c77f2-1a5f-437c-adb6-232155f2cafb
# ╟─c261b580-2cc7-4aba-af15-1778f9293278
# ╠═730de2d2-104b-4f4a-96da-4f96b0445f95
# ╠═2d8df796-c362-49d9-b938-26f1917e800b
# ╟─6b3ee543-e3e2-4fe6-925d-1f201243d014
# ╠═89830894-4481-480b-954d-b24971590eaf
# ╠═c2b5c2c4-3aa6-466d-b106-d0d7aef43bdf
# ╟─3814640e-ed0b-48a1-b4ac-6b7f04d899f3
# ╠═bbe8980e-441d-4ac6-abc5-88ff77de9f24
# ╟─a8257414-2d6a-4913-9f05-6369c62189e8
# ╟─882f2f6b-eb58-4c8f-877e-76e460900c0c
# ╠═f76e5c58-12fc-4b2d-8718-e15144fbef98
# ╟─0bd84cf5-c0a2-4bb2-b5f4-f16da6eaab04
# ╠═cb02ebf8-aee5-47a2-9496-6f63fc123147
# ╟─869a9dc3-8227-4c6e-8901-198cafb489e2
