### A Pluto.jl notebook ###
# v0.20.28

using Markdown
using InteractiveUtils

# ╔═╡ a3c1055e-b44f-4412-9ca0-f8ec5f972494
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ fd357fe0-0920-11ef-02e4-25a84575c6a2
using Markdown

# ╔═╡ 0ddc1082-bf8d-4357-a056-2fa861824b47
using InteractiveUtils

# ╔═╡ cce8fcde-069f-49b9-aadd-6abccbcfbace
using ModelingToolkit

# ╔═╡ ab75267c-d28f-47ce-88bb-5aeb885398e4
using ModelingToolkit: t_nounits as t, D_nounits as D

# ╔═╡ 88aefc2c-5114-48ff-843b-e0d8f08ba3b5
using Turing

# ╔═╡ b48e4b46-3a07-4ff0-b617-6bc1125f6a1e
using OrdinaryDiffEq, StatsPlots

# ╔═╡ ddc035a2-3b03-4a54-807c-05a271384ca7
using ForwardDiff

# ╔═╡ aaa54f95-466f-463e-958a-fc02560fd7d1
using GlobalSensitivity

# ╔═╡ 2293569d-491b-4b7f-a065-0c17d43b56d5
md"""
# Hermione's PhD on plant insect interactions
"""

# ╔═╡ 9bc7b562-e1d2-4395-beb6-92f413589928
md"""
During her PhD Hermione has studied in her first year the plant Cucurbita pepo in a greenhouse. She determined the following parameters.
- the carrying capacity $k$ of her experimental setup is 4000 kg/ha
- the crop growth rate θ is 0.2 $d^{-1}$

In the second year she studied the leaf-cutter ant (Atta cephalotes) and the effect of a neurotoxin on the ant from which she determined the following parameters. 
- The mortality rate μ equals 0.1 $d^{-1}$
- the insecticide effect factor $p$ equals 3

In her last year, she is interested in the interactions between Cucurbita pepo and the leaf-cutter ant in her experimental setup. She knows that the relationship between a field crop (C) and a voracious insect (A) is described by the bitrophic model:

$$\begin{eqnarray*}
\frac{dC}{dt} &= \theta C \left(1-\frac{C}{k}\right)-fCA \\
\frac{dA}{dt} &= \phi f CA -(1 + p)\, \mu A
\end{eqnarray*}$$

| Parameter | Value | Description |
| :--- | :---: | :--- |
| $\theta$ | $0.2$ | Growth rate of $C$ |
| $k$ | $4000$ | Carrying capacity |
| $p$ | $3$ | Insecticide effect |
| $\mu$ | $0.1$ | Natural mortality rate |
| $f$ | **?** | feeding rate |
| $\phi$ | **?** | Conversion efficiency |

Based on a literature review she retrieves some initial ideas of the distributions of the feeding rate $f$ and the conversion $ϕ$.
- The feeding rate $f$ is LogNormal distributed around $f=0.001\;ha/(kg\,d)$ as `LogNormal(log(0.001), 0.5)`. In literature the values for $f$ are in the range of `[0.0005, 0.002]`.
- The conversion efficiency $ϕ$ is somewhere between 0.07 and 0.3 with a peak around 0.15.

**Her final goal is to be able to predict the pumpkin harvest in natural habitats with ants after 200 days.** She wants to know if in natural habitats the yield can reach over 2000 kg/ha, because that makes it profitable for farmers.
"""

# ╔═╡ 07495e3b-9b66-4f58-9b23-db461751e3ca
md"""
## Roadmap of Hermione's third year
"""

# ╔═╡ 39eeb5cf-2c10-40d3-a628-4778d593776a
md"""
1. Gain insights in the overall sensitivity of the model. What parameters have a major influence on the pumpkin harvest?
2. Determine the influence of the unknown parameters on the harvested pumpkin. What are the chances a farmer has less than 2000 kg/ha yield?
3. After the first year, the farmer reports a yield of 2400 kg/ha. Can you a priori know which parameter ($ϕ$ or $f$) will be better estimated using this extra measurement? Secondly use this measurement to update your beliefs about these parameters.
"""

# ╔═╡ 60a25144-e430-45cf-b250-3d196bdb4723
md"""
!!! tip "Step 0"
	Implement the model, use the default values `f=0.001`and `ϕ=0.15`. The initial crop density is `100 kg/ha` and the intial ant population is `0.5 kg/ha`. We want to simulate the harvest after `200 days`.
"""

# ╔═╡ fdd9a533-c441-4b77-bba8-687bddf3e4cb
@variables C(t) A(t)

# ╔═╡ 30ca6976-4336-4fe9-953c-2fc37a9f2288
@parameters θ k f ϕ p μ

# ╔═╡ 4f91272e-ddd9-4bbd-aab7-f6883566b037
eqns = [D(C) ~ θ*C*(1-C/k)-f*C*A,
	   D(A) ~ ϕ*f*C*A - (1+p)*μ*A]

# ╔═╡ 7769c7a5-5f19-4e1a-8721-41ee03c2024e
u0 = [:C => 100, :A => 0.5]

# ╔═╡ 00453a3c-88d1-4cc4-a8f2-d4fa2eed0b66
tspan = (0.0, 200.0)

# ╔═╡ 157eae84-9c8a-4d23-9ada-9d53ba520b0a
params = [:θ => 0.2, :k => 4000.0, :f => 0.001, :ϕ => 0.15, :p => 3, :μ => 0.1]

# ╔═╡ ab331db3-9eb6-4c6c-a0a1-63496e73a183
@mtkbuild ode_sys = ODESystem(eqns, t)

# ╔═╡ d6c1c97d-bd3b-4902-8224-31c119bc3026
ode_prob = ODEProblem(ode_sys, u0, tspan, params)

# ╔═╡ 1995aff0-c888-4a4c-8311-cc8792c41f98
ode_sol = solve(ode_prob, Tsit5(), saveat = 0.5)

# ╔═╡ 0f340b0c-468a-4570-a5f8-6f72bc384798
plot(ode_sol,lw = 2, idxs = [:A :C], title = "Bitrophic model")

# ╔═╡ f741a2eb-48dc-4ff4-817e-849a6dd839d3
md"""
## Part I: What determines the pumpkin harvest?
"""

# ╔═╡ ed3c3b7e-bb51-442b-a455-eebff351b227
md"""
!!! tip "Question 1"
    Gain insights in the overall sensitivity of the model. What is the potential influence of $f$ and $ϕ$ on the final harvest?
"""


# ╔═╡ 592a9822-f8e9-40ed-acdf-ee93975b6784
p_range = [[0.07, 0.3], [0.0005, 0.002]];

# ╔═╡ 2a858811-9ffd-4dc7-84ba-f7a59421753f
function bitrophic_model_sim(params)
	ϕ, f = params
	model_params = [:f => f, :ϕ => ϕ]
	oprob = remake(ode_prob, p = model_params)
	osol = solve(oprob, Tsit5(), saveat=0.5)
	return osol[:C]
end

# ╔═╡ fd71123c-dea5-49a9-91bd-32a59d97c0da
gsa_res = gsa(bitrophic_model_sim, Sobol(), p_range, samples = 200)

# ╔═╡ 48123288-495f-4a96-89ba-f35c6fc4dc44
bar(["ϕ", "f"], gsa_res.ST[end,:], colour = :orange,
   title = "global sensitivities on harvested mass of pumpkin",
   ylims = [0, 1])

# ╔═╡ 3d1f2c90-2e4e-4fa1-af41-9d44d00e5090
t_vals = 0:0.5:200.0

# ╔═╡ 06528d30-37cf-43b0-b0bd-82426df2f964
begin
    labels = ["ϕ", "f"]

    plt = plot(
        xlabel = "Time (day)",
        ylabel = "Total-order Sobol index",
        title = "Global sensitivity analysis",
        ylims = (0, 2),
        legend = :outerright
    )

    for i in eachindex(labels)
        plot!(plt, t_vals, gsa_res.ST[:, i], label = labels[i])
    end

    plt
end

# ╔═╡ ee3ff488-dfbc-410a-ab4c-23d854272dab
md"""
**Conclusion:** From the global sensitivity analysis we can conclude that both parameters explain approximately as much variance at harvest time, so we have to take both parameters into account. Looking at the global sensitivity through time, you can deduce that the feeding rate has a high influence on the growth phase of the crop (around t = 25).
"""

# ╔═╡ e296e434-017c-45ce-bfe0-76601e277612
md"""
## Part II: How reliable/variable are the outcomes?
"""

# ╔═╡ 8f8c1f5d-a5c5-4da0-8230-d81b9ce7dab5
md"""
!!! tip "Question 2"
	Determine the influence of the unknown parameters on the harvested pumpkin. What are the chances a farmer has less than 2000 kg/ha yield?
"""

# ╔═╡ e9f132c6-4828-43bc-9c7b-7b93515d4322
@model function uncertain_bitrophic()
	ϕ ~ TriangularDist(0.07, 0.3, 0.15)
	f ~ LogNormal(log(0.001), 0.5)
	oprob = remake(ode_prob, p = [:ϕ => ϕ, :f => f])
	return solve(oprob, Tsit5())
end

# ╔═╡ 3f419b58-52a9-4cca-ad3f-5a462e3174c8
bitrophic_chain = sample(uncertain_bitrophic(), Prior(), 200)

# ╔═╡ 69906d82-e61a-42eb-a370-b293c22b9b83
prior_sols = generated_quantities(uncertain_bitrophic(), bitrophic_chain);

# ╔═╡ 9fad293a-46e4-474a-9d71-e727772563a2
begin 
	p2 = plot(title = "uncertainty propagation")
	plot!(ode_sol,lw = 5, idxs = :C)
	for sol in prior_sols
		plot!(sol, colour = :gray, alpha = 0.3, label = false, ls = :dot, idxs = :C)
	end
	p2
end

# ╔═╡ 13331268-92f3-468d-9414-564b73b6c7c3
sol_end = collect([prior_sols[i][:C][end] for i in eachindex(prior_sols)])

# ╔═╡ ec8906b5-6388-46cf-844d-d49578e8f9fc
histogram(sol_end, title = "Distribution of yield at harvest time",
		 label = "counts")

# ╔═╡ 001649f9-f722-4e5a-a2c3-48b01b677c57
chance = mean(sol_end .< 2000.0)

# ╔═╡ ae341780-b8e9-4b73-a8fe-26efd8588e6a
println("The chance of having less than 2000 kg/ha is : $(chance*100)%")

# ╔═╡ 1b092702-1614-4fb6-90e1-df97aab1183c
md"""
**Conclusion:** Based on our Prior belief the farmer has 39% to have less than 2000 kg/ha yield.
"""

# ╔═╡ af93263f-fda7-4cad-8f3d-74ab4a5b345d
md"""
## Part III: Can we reduce uncertainty after year one?
"""

# ╔═╡ fe08b97c-08ad-4a19-8272-b047bc652ec6
md"""
!!! tip "Question 3"
	After the first year, the farmer reports a yield of 2400 kg/ha. Can you a priori know which parameter ($ϕ$ or $f$) will be better estimated using this extra measurement? Secondly use this measurement to update your beliefs about these parameters. What is the updated chance that the farmer will have less than 2000 kg/ha yield?

	Make use of the Metropolis-Hastings algorithm and use 5000 samples. You can assume a standard deviation of 50 kg/ha on the measured harvest yield.
"""

# ╔═╡ 81aa36db-b28f-4dff-9513-cfc3b80daaab
begin
	f_val = 0.001
	ϕ_val = 0.15
end

# ╔═╡ 4609ffec-3400-4f75-9ffa-a8c5006c3b48
begin
    p0 = [ϕ_val, f_val]
    C_nom = bitrophic_model_sim(p0)
    jac = ForwardDiff.jacobian(bitrophic_model_sim, p0)
    rel_sens_ϕ = jac[:, 1] ./ C_nom .* ϕ_val
    rel_sens_f = jac[:, 2] ./ C_nom .* f_val
end

# ╔═╡ 850feaf3-2b20-4ca9-988e-3f934a767a4a
begin
	plot(t_vals, rel_sens_f, label = "sens f", title = "normalized sensitivities",
		xlabel = "time")
	plot!(t_vals, rel_sens_ϕ, label = "sens ϕ")
end

# ╔═╡ 12803a13-ccbb-4a98-9a23-2c4e5e5559b7
md"""
**Conclusion**: A priori, it is not possible to define which parameter is most sensititive for the extra measurement. So based on this, we do not know which parameter will benefit the most. A second observation is that the parameters are less sensitive at the end of the simulation period. If you could choose to pick one measurement day, a measurement after 85 days would be more informative.
"""

# ╔═╡ 2c9fe511-e31c-4cd6-b60f-19798803b715
@model function uncertain_harvest()
	ϕ ~ TriangularDist(0.07, 0.3, 0.15)
	f ~ LogNormal(log(0.001), 0.5)
	osol = solve(remake(ode_prob, p = [:ϕ => ϕ, :f => f]), AutoTsit5(Rosenbrock23()))
	harvest ~ Normal(osol[:C][end], 50.0)
	return harvest
end

# ╔═╡ a8b2df6a-ec78-49cf-9226-e11a7cf6d522
harvest_cond = uncertain_harvest() | (harvest = 2400.0,)

# ╔═╡ 527ae9f4-8a0b-4e3c-9923-d240923f242c
cond_chain = sample(harvest_cond, MH(), 5000)

# ╔═╡ fef6de4a-e362-4422-8097-7298e8dbcfff
plot(cond_chain)

# ╔═╡ 1ccdb12d-a30a-49ae-b79a-2c12334ffe21
begin
	println("mean f is: $(mean(cond_chain[:f]))")
	histogram(cond_chain[:f], title = "Posterior distribution f")
end

# ╔═╡ e43afba4-0bc6-417b-b44b-244e1ab2436e
begin
	println("mean ϕ is: $(mean(cond_chain[:ϕ]))")
	histogram(cond_chain[:ϕ], title = "Posterior distribution ϕ")
end

# ╔═╡ 77df7791-2411-47b9-b891-af6c00fbde18
begin
	p3 = plot(title = "Updated beliefs model")
	harvest_end = []
	for (f, ϕ) in zip(cond_chain[:f][1:5:end], cond_chain[:ϕ][1:5:end])
		prob = remake(ode_prob, p = [:f => f, :ϕ => ϕ])
		osol_new = solve(prob, Tsit5())
		push!(harvest_end, osol_new[:C][end])
		plot!(osol_new, label = false, color = :gray, ls = :dot)
	end
	plot!(ode_sol, lw = 2, idxs = [:A :C])
	p3
end

# ╔═╡ 424871d9-dc81-49a7-804b-a25865962c09
println("chances that harvest is smaller than 2000 kg/ha: $(mean((harvest_end .< 2000)))")

# ╔═╡ 69b26cd0-7964-42f8-b471-e85118114867
md"""
**Conclusion:** Based on the simulation results above, the chance of having less than 2000 kg/ha yield has shrunk to approximately 0%.

**Side note:** Although we obtained a posterior distribution, a single harvest measurement cannot pin down $f$ and $ϕ$ individually. Both parameters have a similar influence on the yield at $t = 200$ d (see the sensitivity analysis), so many *combinations* of $f$ and $ϕ$ explain the observed 2400 kg/ha equally well. The data constrain the combination of parameters, not each parameter on its own. Yet the *prediction* of the harvest has become very certain.
"""

# ╔═╡ 9fafbd32-28de-4e8e-b53a-8b5fe5d48033
scatter(cond_chain[:f], cond_chain[:ϕ], xlabel = "f", ylabel = "ϕ", alpha = 0.3, label = false, title = "Posterior relation between ϕ and f")

# ╔═╡ Cell order:
# ╠═fd357fe0-0920-11ef-02e4-25a84575c6a2
# ╠═0ddc1082-bf8d-4357-a056-2fa861824b47
# ╠═a3c1055e-b44f-4412-9ca0-f8ec5f972494
# ╠═cce8fcde-069f-49b9-aadd-6abccbcfbace
# ╠═ab75267c-d28f-47ce-88bb-5aeb885398e4
# ╠═88aefc2c-5114-48ff-843b-e0d8f08ba3b5
# ╠═b48e4b46-3a07-4ff0-b617-6bc1125f6a1e
# ╠═ddc035a2-3b03-4a54-807c-05a271384ca7
# ╠═aaa54f95-466f-463e-958a-fc02560fd7d1
# ╟─2293569d-491b-4b7f-a065-0c17d43b56d5
# ╟─9bc7b562-e1d2-4395-beb6-92f413589928
# ╟─07495e3b-9b66-4f58-9b23-db461751e3ca
# ╟─39eeb5cf-2c10-40d3-a628-4778d593776a
# ╟─60a25144-e430-45cf-b250-3d196bdb4723
# ╠═fdd9a533-c441-4b77-bba8-687bddf3e4cb
# ╠═30ca6976-4336-4fe9-953c-2fc37a9f2288
# ╠═4f91272e-ddd9-4bbd-aab7-f6883566b037
# ╠═7769c7a5-5f19-4e1a-8721-41ee03c2024e
# ╠═00453a3c-88d1-4cc4-a8f2-d4fa2eed0b66
# ╠═157eae84-9c8a-4d23-9ada-9d53ba520b0a
# ╠═ab331db3-9eb6-4c6c-a0a1-63496e73a183
# ╠═d6c1c97d-bd3b-4902-8224-31c119bc3026
# ╠═1995aff0-c888-4a4c-8311-cc8792c41f98
# ╠═0f340b0c-468a-4570-a5f8-6f72bc384798
# ╟─f741a2eb-48dc-4ff4-817e-849a6dd839d3
# ╟─ed3c3b7e-bb51-442b-a455-eebff351b227
# ╠═592a9822-f8e9-40ed-acdf-ee93975b6784
# ╠═2a858811-9ffd-4dc7-84ba-f7a59421753f
# ╠═fd71123c-dea5-49a9-91bd-32a59d97c0da
# ╠═48123288-495f-4a96-89ba-f35c6fc4dc44
# ╠═3d1f2c90-2e4e-4fa1-af41-9d44d00e5090
# ╠═06528d30-37cf-43b0-b0bd-82426df2f964
# ╟─ee3ff488-dfbc-410a-ab4c-23d854272dab
# ╟─e296e434-017c-45ce-bfe0-76601e277612
# ╟─8f8c1f5d-a5c5-4da0-8230-d81b9ce7dab5
# ╠═e9f132c6-4828-43bc-9c7b-7b93515d4322
# ╠═3f419b58-52a9-4cca-ad3f-5a462e3174c8
# ╠═69906d82-e61a-42eb-a370-b293c22b9b83
# ╠═9fad293a-46e4-474a-9d71-e727772563a2
# ╠═13331268-92f3-468d-9414-564b73b6c7c3
# ╠═ec8906b5-6388-46cf-844d-d49578e8f9fc
# ╠═001649f9-f722-4e5a-a2c3-48b01b677c57
# ╠═ae341780-b8e9-4b73-a8fe-26efd8588e6a
# ╟─1b092702-1614-4fb6-90e1-df97aab1183c
# ╟─af93263f-fda7-4cad-8f3d-74ab4a5b345d
# ╟─fe08b97c-08ad-4a19-8272-b047bc652ec6
# ╠═81aa36db-b28f-4dff-9513-cfc3b80daaab
# ╠═4609ffec-3400-4f75-9ffa-a8c5006c3b48
# ╠═850feaf3-2b20-4ca9-988e-3f934a767a4a
# ╟─12803a13-ccbb-4a98-9a23-2c4e5e5559b7
# ╠═2c9fe511-e31c-4cd6-b60f-19798803b715
# ╠═a8b2df6a-ec78-49cf-9226-e11a7cf6d522
# ╠═527ae9f4-8a0b-4e3c-9923-d240923f242c
# ╠═fef6de4a-e362-4422-8097-7298e8dbcfff
# ╠═1ccdb12d-a30a-49ae-b79a-2c12334ffe21
# ╠═e43afba4-0bc6-417b-b44b-244e1ab2436e
# ╠═77df7791-2411-47b9-b891-af6c00fbde18
# ╠═424871d9-dc81-49a7-804b-a25865962c09
# ╟─69b26cd0-7964-42f8-b471-e85118114867
# ╠═9fafbd32-28de-4e8e-b53a-8b5fe5d48033
