### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 59243782-5b91-4afd-ab68-5283948781d3
using Pkg; Pkg.activate("..")

# ╔═╡ 27bd53ba-1ea7-11f0-3d13-65f334bf4250
using PlutoUI, Plots

# ╔═╡ b2f9ff92-f75f-479c-8d87-c34719661b3c
using Catalyst, DifferentialEquations

# ╔═╡ 3ebf762f-6e72-4811-8884-bf7ba98321e1
using ModelingToolkit

# ╔═╡ 92b24c3b-9a44-43aa-a446-998e2a1ed4c4
md"""
# Modelling mass and temperature changes

Given a degradation reaction 

`C --> ∅`

which has a kinetics of:

$$C' = - k C$$

that is temperature-dependent, with a reaction rate given by the *Arrhenius equation*:

$$k(T) = A e^{-\frac{E_a}{RT}}$$

Suppose our temperature is not constant but follows *Newton's cooling law*:

$$T' = - r (T-T_e)$$

How can we model this using our Julia tools?
"""

# ╔═╡ bf413087-6f78-47ba-bb8b-e33514ee1cfa
md"""
## 1. Abusing Catalyst.jl

Without using anything new, we can implement cooling as a first-order decay of the temperature difference. Note that we use `@observables` to keep track of the true temperature!
"""

# ╔═╡ 023cfc5d-6e2e-425e-a1fc-edb223e01589
sys1 = @reaction_network begin
	@observables T ~ Te + ΔT
	# ΔT = T - Te
	A * exp(-Ea/R/(Te+ΔT)), C --> 0
	r, ΔT --> 0
end

# ╔═╡ 134c7a48-01d1-4313-a9fe-157b7dab2a17
convert(ODESystem, sys1)

# ╔═╡ bc67c9d8-ae53-4195-9f35-5299f3dd7941
md"## 2. Decoupling

see [here](https://docs.sciml.ai/Catalyst/stable/model_creation/functional_parameters/)

Here, the temperature is independent of the concentration; we can solve for T (linear first-order ODE) and plug in a catalyst model. It requires a bit of computation and will fail when were have a more complex model.
"

# ╔═╡ 306f3a04-f53c-486f-8656-29ecfa56e40d
import ModelingToolkit: t  # import time

# ╔═╡ 3f2f4599-6b5e-46c8-84f9-dea9bb0183a5
md" ## 3. Using DifferentialEquations

Catalyst is not ideal for things beyond mass balances; you can just implement the ODE system directly.

This is precisely the same thing, though we build the ODE system here manually rather than automatically.
"

# ╔═╡ 509a36ed-fcd2-43ae-8768-f1b39f94a34d
function chemcooling!(du, u, (A, Ea, R, Te, r), t)
	C, T = u
	du[1] = dC = - A * exp(-Ea / R / T) * C
	du[2] = dT = - r * (T - Te)
end

# ╔═╡ b3c3bbe6-ee7f-40b8-be92-116eb1043368
md"""
## 4. Adding an extra equation using MTK

The most advanced and powerful way is adding additional equations to the system. This is the true scalable way. We can combine ODEs with reaction systems to create wildly complex models automatically.

see [here](https://docs.sciml.ai/Catalyst/stable/model_creation/constraint_equations/#constraint_equations_coupling_constraints)
"""

# ╔═╡ f9294388-e61d-4888-8548-dd6a1811944f
D = default_time_deriv()  # time derivative operator

# ╔═╡ 5fc2cf63-5bfa-49b9-acd9-b972bc5fbaee
sys4 = @reaction_network begin
	A * exp(-Ea / (R * T)), C --> 0
end

# ╔═╡ 63a4d6d6-d339-46a6-a728-931f015a997f
#sys4 = @reaction_network begin
#	$k, C --> 0   # alternative, splice in the code for k
#end

# ╔═╡ 01ea3432-a53f-4f16-9f13-cb5fa0245d08
md"## Parameters"

# ╔═╡ ee9550cf-4d00-4fd8-80cc-fe9cbf49e166
r = 0.018  # cooling rate

# ╔═╡ 25aaf503-95b3-44c4-a270-88b987ad1371
R = 8.314  # gas cst

# ╔═╡ ad3f22b6-a1a9-40c0-868a-3a6a63bea504
A = 0.17

# ╔═╡ 0e90bb09-0ae9-43d3-9d19-fe89f1d743ad
Ea = 89.3

# ╔═╡ 0eeafb45-5bbb-43ca-9cb3-a4ebe06dc64f
k_temp(T) = A * exp(-Ea / (R * T))

# ╔═╡ fac4b66e-533a-413f-8fbb-cf312244c10b
plot(k_temp, 20, 100, xlab="T [degree Celsius]", ylab="reaction rate", label="k")

# ╔═╡ 8952f810-79f8-42ec-948e-b78f847c154d
C0 = 50.0  # initial concentration

# ╔═╡ ec87305a-87ad-403b-a6a5-827f6c9f80e7
tspan = (0, 60)

# ╔═╡ a90c6b66-2cd7-4bd4-be1b-5a75a0088013
T0 = 100  # initial temperature

# ╔═╡ 61bce0a0-1a76-472c-886e-fc48c109df46
@variables T(t) = T0  # define temperature (function of time)

# ╔═╡ 1ecfbbda-f9a2-4382-b7f2-b8fbf1813419
k =  A * exp(-Ea / (R * T))  # this is an equation that uses the variable T

# ╔═╡ 90ad7025-500f-4526-8379-2aa0bf41b949
Te = 21  # env temp

# ╔═╡ 9b249ef1-eff5-41a3-aa93-3683560356a5
temperature(t) = (T0 - Te) * exp(-r*t) + Te

# ╔═╡ f57e5c92-aa19-4811-82ce-991ffa935374
plot(temperature, 0, 50, xlab="t", ylab="T")

# ╔═╡ 74876f57-8b41-4437-a9eb-4697e8cd6f53
# rather than model C & T, we eliminate T and make
# the reaction rate depend on the time
sys2 = @reaction_network begin
	A * exp(-Ea/R/(temperature(t))), C --> 0
end

# ╔═╡ 07e62e76-4271-415c-bfe4-cbf98f409208
odesys2 = ODEProblem(sys2, [:C => C0], tspan, [:A=>A, :R=>R, :Ea=>Ea]);

# ╔═╡ e5390072-ddfa-483d-9bbc-dbbf56a41480
sol2 = solve(odesys2);

# ╔═╡ 86f46366-9e6c-4315-898f-01480d937320
plot(sol2)

# ╔═╡ 7c0adb22-7daf-4110-9378-687843771cef
sys3 = ODEProblem(chemcooling!, [C0, T0], tspan, (A, Ea, R, Te, r));

# ╔═╡ 33228136-4e19-451d-96d8-c90eef259330
sol3 = solve(sys3);

# ╔═╡ 53cba252-77f9-463c-9f04-440a5bf7eba1
plot(sol3, label=["C" "T"])

# ╔═╡ 4c0c885e-6450-4329-9d4c-32d1b68dc687
newton_cooling_law = (D(T) ~ - r * (T - Te)) 

# ╔═╡ 8bd4b6e7-c7e6-4156-83d1-40aa0f6717d5
@named osys = ODESystem(newton_cooling_law, t)

# ╔═╡ 1295d9fc-f2ba-46a1-86e5-cf093fe31533
begin
	@named sys_with_cooling = extend(osys, sys4)
	sys_with_cooling = complete(sys_with_cooling)
end

# ╔═╡ 974311a5-515c-4dcf-bc05-63f4a2c85553
convert(ODESystem, sys_with_cooling)

# ╔═╡ 4b4e1bb0-d753-4c2e-aac4-ecf564be404f
odesys4 = ODEProblem(sys_with_cooling, [:C=>C0, :T=>T0], tspan, [:A=>A, :R=>R, :Ea=>Ea]);

# ╔═╡ eff6e2ff-99ad-4a82-8881-3f525d28c39e
sol4 = solve(odesys4);

# ╔═╡ 7f6a4cba-6656-4098-9852-3865cc68b70e
plot(sol4)

# ╔═╡ 5f0b6722-35bd-4818-be4e-9b79895ca509
ΔT0 = T0 - Te

# ╔═╡ 82cb2040-3b2f-4599-ab93-0b9fca9eee3b
odesys1 = ODEProblem(sys1, [:C => C0, :ΔT=>ΔT0], tspan, [:A=>A, :r=>r, :R=>R, :Ea=>Ea, :Te=>Te]);

# ╔═╡ 6e072039-2f16-478b-99c4-207505e85397
sol1 = solve(odesys1);

# ╔═╡ b3e54b72-4ab6-42c9-afda-ccb61844fa54
plot(sol1, idxs=[:C, :T])

# ╔═╡ f39e513d-f2ab-45fd-9621-fa893c793c47
TableOfContents()

# ╔═╡ Cell order:
# ╠═27bd53ba-1ea7-11f0-3d13-65f334bf4250
# ╠═59243782-5b91-4afd-ab68-5283948781d3
# ╟─92b24c3b-9a44-43aa-a446-998e2a1ed4c4
# ╠═0eeafb45-5bbb-43ca-9cb3-a4ebe06dc64f
# ╠═fac4b66e-533a-413f-8fbb-cf312244c10b
# ╟─bf413087-6f78-47ba-bb8b-e33514ee1cfa
# ╠═b2f9ff92-f75f-479c-8d87-c34719661b3c
# ╠═023cfc5d-6e2e-425e-a1fc-edb223e01589
# ╠═134c7a48-01d1-4313-a9fe-157b7dab2a17
# ╠═82cb2040-3b2f-4599-ab93-0b9fca9eee3b
# ╠═6e072039-2f16-478b-99c4-207505e85397
# ╠═b3e54b72-4ab6-42c9-afda-ccb61844fa54
# ╟─bc67c9d8-ae53-4195-9f35-5299f3dd7941
# ╠═306f3a04-f53c-486f-8656-29ecfa56e40d
# ╠═9b249ef1-eff5-41a3-aa93-3683560356a5
# ╠═f57e5c92-aa19-4811-82ce-991ffa935374
# ╠═74876f57-8b41-4437-a9eb-4697e8cd6f53
# ╠═07e62e76-4271-415c-bfe4-cbf98f409208
# ╠═e5390072-ddfa-483d-9bbc-dbbf56a41480
# ╠═86f46366-9e6c-4315-898f-01480d937320
# ╟─3f2f4599-6b5e-46c8-84f9-dea9bb0183a5
# ╠═509a36ed-fcd2-43ae-8768-f1b39f94a34d
# ╠═7c0adb22-7daf-4110-9378-687843771cef
# ╠═33228136-4e19-451d-96d8-c90eef259330
# ╠═53cba252-77f9-463c-9f04-440a5bf7eba1
# ╟─b3c3bbe6-ee7f-40b8-be92-116eb1043368
# ╠═3ebf762f-6e72-4811-8884-bf7ba98321e1
# ╠═f9294388-e61d-4888-8548-dd6a1811944f
# ╠═61bce0a0-1a76-472c-886e-fc48c109df46
# ╠═4c0c885e-6450-4329-9d4c-32d1b68dc687
# ╠═8bd4b6e7-c7e6-4156-83d1-40aa0f6717d5
# ╠═1ecfbbda-f9a2-4382-b7f2-b8fbf1813419
# ╠═5fc2cf63-5bfa-49b9-acd9-b972bc5fbaee
# ╠═63a4d6d6-d339-46a6-a728-931f015a997f
# ╠═1295d9fc-f2ba-46a1-86e5-cf093fe31533
# ╠═974311a5-515c-4dcf-bc05-63f4a2c85553
# ╠═4b4e1bb0-d753-4c2e-aac4-ecf564be404f
# ╠═eff6e2ff-99ad-4a82-8881-3f525d28c39e
# ╠═7f6a4cba-6656-4098-9852-3865cc68b70e
# ╟─01ea3432-a53f-4f16-9f13-cb5fa0245d08
# ╠═ee9550cf-4d00-4fd8-80cc-fe9cbf49e166
# ╠═25aaf503-95b3-44c4-a270-88b987ad1371
# ╠═ad3f22b6-a1a9-40c0-868a-3a6a63bea504
# ╠═0e90bb09-0ae9-43d3-9d19-fe89f1d743ad
# ╠═8952f810-79f8-42ec-948e-b78f847c154d
# ╠═ec87305a-87ad-403b-a6a5-827f6c9f80e7
# ╠═a90c6b66-2cd7-4bd4-be1b-5a75a0088013
# ╠═90ad7025-500f-4526-8379-2aa0bf41b949
# ╠═5f0b6722-35bd-4818-be4e-9b79895ca509
# ╠═f39e513d-f2ab-45fd-9621-fa893c793c47
