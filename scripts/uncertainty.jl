### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ a981929c-2a53-11ef-16ce-f3918dd88313
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	using Pkg
	Pkg.activate("..")
end
  ╠═╡ =#

# ╔═╡ 1f5c60e2-e3ec-41ba-9c3e-44083987d710
using Catalyst, DifferentialEquations

# ╔═╡ 58eb5c0f-557d-42a1-ae5d-728de0c498bd
using Measurements, ForwardDiff, Turing

# ╔═╡ a5cbf15c-490e-4d63-9a99-1c579be36c7e
using Plots, LaTeXStrings, Latexify, PlutoUI, StatsPlots

# ╔═╡ d8818694-9270-4882-bc50-2eff71ea1c7a
using GlobalSensitivity

# ╔═╡ 9629b2f8-2c21-40f1-a48d-a206f8bad479
md"""
# Uncertainty & sensitivity analysis
"""

# ╔═╡ 288c20e5-f1b2-4fbc-bc52-c1ec0e9d8422
md"""
## Aleatoric vs epistemic uncertainty
"""

# ╔═╡ 2100c8b0-0f88-4228-a7f9-dc62f4a3cf9c
flogistic = (u, (r, K), t) -> r*u*(1-u/K)

# ╔═╡ 48856eee-8568-483a-bd0f-f3fb42c58906
logode = ODEProblem(flogistic, 2.0, (0.0, 50.0), (0.28, 100))

# ╔═╡ 6f6ee7e0-3413-43c5-a9bd-132267a96d82
logsde = SDEProblem(flogistic, (u,p,t)->0.25sqrt(u), 2.0, (0, 50), (0.28, 100))

# ╔═╡ e887a9ab-a4d8-4b05-a3f8-a8165d517a8a
β = 0.03

# ╔═╡ 733ae8b7-b3d6-443a-b4a1-f8fdbde0e0b6
γ = 0.3

# ╔═╡ 9813f4f9-1ae0-4480-bf9d-14ab547d01c2
md"## SIR example"

# ╔═╡ 9bed9078-0c22-4f29-96db-4274999adf9c
sir = @reaction_network begin
	@species S(t)=99 I(t)=1 R(t)=0
	@parameters β=$β γ=$γ
    β, S + I --> 2I  # susceptible persons get infected
    γ, I --> R       # infected persons become recovered
end

# ╔═╡ 34f47469-63b7-4e9a-98e8-e67827e4dbf4
convert(ODESystem, sir)

# ╔═╡ e34de9d2-6189-46fa-88a6-a37f94294384
sirprob = ODEProblem(sir, [], (0., 30.));

# ╔═╡ 80eeefec-1fd5-4880-9209-df27375f9e1e
sirsol = solve(sirprob)

# ╔═╡ bdf2dd3b-55c8-445f-b9ed-66805014aa47
md"""
## Uncertainty
"""

# ╔═╡ ea341f0b-b081-4eda-b7d8-5ca616261ee8
βu  = β ± 0.006

# ╔═╡ 6ed4121b-2b68-4446-bf24-e8ce61661f84
γu = γ ± 0.04

# ╔═╡ 2645bb50-bec6-4b5e-b423-9cbc2f4dbd24
sirprobu = remake(sirprob, p=[:β=>βu, :γ=>γu])

# ╔═╡ 02209c47-4ff3-4823-a998-fe1f173f749a
@model function sir_mc()
	β ~ Normal(0.03, 0.006)
	γ ~ TriangularDist(0.23, 0.37, 0.3)
	return solve(remake(sirprob, p=[:β=>β, :γ=>γ]))
end

# ╔═╡ 9b048268-30d7-4eaf-8a26-ab2e1c0404cd
throw = rand(sir_mc())

# ╔═╡ 4f1cc2c6-254d-4f9c-b97a-09d22cb3838a
generated_quantities(sir_mc(), throw)

# ╔═╡ 118f4f61-12fb-4025-bd01-91a7fa6f9bd7
md"[Using the Beer-Lambert Law to Calculate the Concentration of a Solution | Chemistry | Study.com](https://study.com/skill/learn/using-the-beer-lambert-law-to-calculate-the-concentration-of-a-solution-explanation.html)"

# ╔═╡ e6bfb7d1-e594-41c9-9ee2-9668cae70e71
Iobs = 0.2 ± 0.03;

# ╔═╡ ccf61cf5-33b2-4349-a422-16a7c883038d
I0 = 1.29 ± 0.06;

# ╔═╡ 5a9da0cd-9809-445e-aeb3-b91875064106
A = log10(I0 / Iobs)

# ╔═╡ ae807e33-3e05-4a3d-a591-2b0e04c9d6d8
ϵ = 8850 ± 50;

# ╔═╡ 260ba845-ae2f-4e9f-a2d0-f43c64648461
l = 3 ± 5e-4;

# ╔═╡ 8cf4e046-c979-4633-b922-43a46c914f3c
c = log10(I0 / Iobs) / (l * ϵ)

# ╔═╡ a3f82cf6-b25d-4de1-b264-93d6e047ec6d
md"Radius of a circle"

# ╔═╡ 090f89e0-d1e2-4363-830e-4d0635749890
r = 12.5 ± 0.2

# ╔═╡ 660c25d8-bcc5-4de3-9d04-013356acd22a
area = pi * r^2

# ╔═╡ 06929445-f0ac-430f-a09f-1743f2de6f91
md"## Local sensitivity"

# ╔═╡ 11c32674-9817-4b0c-bcb2-235c431a5519
function sir_sim(pars)
	β, γ = pars
	prob = remake(sirprob, p=[:β=>β, :γ=>γ])
	sol = solve(prob, Tsit5(), saveat=0.1)
	return sol
end

# ╔═╡ 1ed87d81-4b19-4de3-9735-fc36894b64c7
sir_sim([β, γ])

# ╔═╡ 633b06bf-ad63-418d-bf2d-899c51eee807
sir_sim_S(pars) = sir_sim(pars)[:S]

# ╔═╡ 42908d09-768d-4004-b304-cb136ed5b791
sir_sim_I(pars) = sir_sim(pars)[:I]

# ╔═╡ 449b2aa3-9efe-4399-8123-95fd952cb062
sir_sim_R(pars) = sir_sim(pars)[:R]

# ╔═╡ 0bb3f53f-50a6-40d3-95ac-be64a8dc79da
tvals = 0:.1:30

# ╔═╡ 10cafa44-8089-4aae-bd27-686be0496dc4
Sv, Iv, Rv = sir_sim_S([β, γ]), sir_sim_I([β, γ]), sir_sim_R([β, γ]);

# ╔═╡ 305b0b4a-4051-409e-9c64-f95ab1778a60
sens_S = ForwardDiff.jacobian(sir_sim_S, [β, γ])

# ╔═╡ 718a957c-74e6-455d-b847-e3d10e7f2016
sens_I = ForwardDiff.jacobian(sir_sim_I, [β, γ])

# ╔═╡ 20af2948-d688-45ee-804c-26e35a59162c
sens_R = ForwardDiff.jacobian(sir_sim_R, [β, γ])

# ╔═╡ 95d7287c-8008-4cad-9af6-ed3e77b9db41
sens_β = [sens_S[:,1] sens_I[:,1] sens_R[:,1]]

# ╔═╡ 4c23f7ac-bf2e-42bd-aa2d-75d98cb691f1
sens_γ = [sens_S[:,2] sens_I[:,2] sens_R[:,2]]

# ╔═╡ 611d9607-99e6-488a-95f2-18da8fb8759e
sens_β_rel = sens_β .* β ./ [Sv Iv Rv]

# ╔═╡ 34cfa9fc-d197-40aa-9987-4736493e533e
sens_γ_rel = sens_γ .* γ ./ [Sv Iv Rv]

# ╔═╡ e436ef36-1f14-48fd-93af-aec0eadc51ed
Ks = 25

# ╔═╡ 4c6849ca-723d-44da-b345-8373e6c67284
mymm = x-> mm(x, 1, Ks)

# ╔═╡ 93f93512-8337-40cd-8697-e46de3fcd945
sens_mm_Ks = x -> ForwardDiff.derivative(K->mm(x, 1, K), Ks)

# ╔═╡ f6aab879-37ad-4425-a5a5-82db2ca4d188
sens_mm_Ks_rel = x -> sens_mm_Ks(x) * Ks / mm(x, 1, Ks) 

# ╔═╡ bffa2e06-271f-4ddb-a6bf-cdf9903b5c31
sens_mm_mumax = x -> ForwardDiff.derivative(μ->mm(x, μ, Ks), 1)

# ╔═╡ 3cf65781-9369-4219-9520-88cfd38adf0e
sens_mm_mumax_rel = x -> sens_mm_mumax(x) / mm(x, 1, Ks)

# ╔═╡ 415cbdec-71f3-4ac7-8479-14d0c7f6481c
sir_summary = function(p)
	stepsize = 0.1
    # beta, gamma = pars
    prob_new_pars = remake(sirprob; p = p)
    sol = solve(prob_new_pars, Tsit5(); saveat = stepsize)
    # effect on total and maximum of infected
	tot_infected = stepsize * sum(sol[2, :])
	max_infected = maximum(sol[2, :])
    return [tot_infected, max_infected]
end

# ╔═╡ 95b84eb0-d3e2-4cc3-882b-d3a4e82c1c4f
β_vals = 1e-5:0.01:0.5

# ╔═╡ 0ac6b303-183b-4256-8436-052071d989cd
γ_vals = 1e-5:0.01:1

# ╔═╡ 3f6491cf-d9d2-4be0-a341-db1a1e891029
md"## Sobol sensitivity"

# ╔═╡ 6d1e0bc0-4b28-4353-8389-c5ff4d156b05
p_range = [[1e-5, 0.5], [1e-5, 2]] # parameter ranges beta and gamma

# ╔═╡ 2d58a23a-7354-4e07-8cbb-ec1351051845
sobol = gsa(sir_summary, Sobol(), p_range, samples = 1000)

# ╔═╡ 5d2ffa7b-1202-45b4-b3a7-536ea38e2bfc
# cols are parameters

# ╔═╡ 8c13eb13-4b54-40ab-8408-ede4f72d7af4
sobol.S1  # first order effects

# ╔═╡ 55b7b2ca-0e80-42a0-9aea-12b46947305c
sobol.ST  # total effects

# ╔═╡ f4869054-c5c6-4028-b9f3-0bb3e6b4305a
md"## Morris"

# ╔═╡ 9e4b6c92-03a9-4ddb-b4f1-97dddf1b7a30
morris = gsa(sir_summary, Morris(num_trajectory=100), p_range, samples = 1000_000)

# ╔═╡ 14676d4b-a43a-4e4f-8eb1-d93c098df156
morris.means

# ╔═╡ 1478f40c-f530-48de-9608-bdd24cab0f23
morris.means_star

# ╔═╡ b455c7c0-d9d5-402f-b49d-fcd1f4b3626a
morris.variances

# ╔═╡ 8928645a-c20e-407c-8ba0-92f3f27206f0
md"""

## Appendix 
"""

# ╔═╡ 3902e862-3bc9-4786-bdc5-0415bc17dbd7
TableOfContents()

# ╔═╡ 30d0f3e0-a6c0-44c8-bcfd-7b163871b4f0
plots = Dict()

# ╔═╡ 3ba4a491-c64d-48de-b263-83679b230b32
let
	p = plot(title="ODE throws of Verhulst model (epistemic)")
	for i in 1:100
		u0 = 10rand()
		r = max(0.28 + 0.1rand(), 0)
		K = 20rand() + 90
		prob = remake(logode, u0=u0, p=(r, K))
		sol = solve(prob, Tsit5())
		plot!(p, sol, lw=0.5, alpha=0.5, color=1, label="")
	end
	plots["epistemic"] = p
end

# ╔═╡ 36eb6a4d-a167-47eb-9208-52d34949f37e
let
	p = plot(title="SDE throws of Verhulst model (aleatoric)")
	for i in 1:100
		prob = remake(logsde)
		sol = solve(logsde, EM(), dt=0.1)
		plot!(p, sol, lw=0.5, alpha=0.5, color=1, label="")
	end
	plots["aleatoric"] = p
end

# ╔═╡ 7305bd7b-e823-451d-9eb1-15f5229fe431
plots["sir_sim"] = plot(sirsol, lw=2, title="SIR model with β=$β and γ=$γ", legend=:right)

# ╔═╡ ef6cb5d2-aa48-459d-931f-a5a79fcf2507
plots["sir_uncertainty"] =plot(solve(sirprobu, saveat=2), lw=2, title="SIR model with β=$βu and γ=$γu", legend=:right)

# ╔═╡ 5075c78a-56e3-4372-b605-de3dfec5cebc
let
	p = plot(sirsol, lw=2, title="Monte Carlo of SIR model", legend=:right)
	for i in 1:200
		throw = rand(sir_mc())
		sol = generated_quantities(sir_mc(), throw)
		plot!(sol, alpha=0.2, color=[1 2 3], label="", lw=0.5)
	end
	plots["sir_MC"] = p
end

# ╔═╡ 7c5a9bca-a568-43d5-88f1-313bcf19bc96
let
	plots["sir_sens_beta"] = pβ = plot(tvals, sens_β, label=["S" "I" "R"], lw=2, xlab="t",
			title="Absolute sensitivity SIR w.r.t. β")
end

# ╔═╡ bb1f08d7-d551-4226-9873-21fc903ed10c
let
	plots["sir_sens_gamma"] = pγ = plot(tvals, sens_γ, label=["S" "I" "R"], lw=2, xlab="t",
			title="Absolute sensitivity SIR w.r.t. γ")
end

# ╔═╡ 57107079-7baa-4bce-b801-b50d1dfc7a5f
let
	plots["sir_sens_beta_rel"] = pβ = plot(tvals, sens_β_rel, label=["S" "I" "R"], lw=2, xlab="t",
			title="Relative sensitivity SIR w.r.t. β")
end

# ╔═╡ 77b7b12f-1b91-4e78-af25-7de8996be2b7
let
	plots["sir_sens_gamma_rel"] = pγ = plot(tvals, sens_γ_rel, label=["S" "I" "R"], lw=2, xlab="t",
			title="Relative sensitivity SIR w.r.t. γ")
end

# ╔═╡ 136e14d5-2330-4e24-8085-4f13fa3b9e26
plots["sir_sens"] = plot(plots["sir_sens_beta"], plots["sir_sens_gamma"], plots["sir_sens_beta_rel"], plots["sir_sens_gamma_rel"], size=(1200, 800))

# ╔═╡ 9563199f-257b-4f89-a4c3-7cd1cc762bac
plots["mm"] = plot(mymm, 0, 100, lw=2, label="MM", xlab="X [mol/L]", ylab="Reaction rate", title="Michaelis-Menten")

# ╔═╡ 32cdba4c-84b9-4449-8a68-4a568592dce7
plots["mm_sens_Ks"] = plot(sens_mm_Ks, 0, 100, lw=2, title="Michaelis-Menten sensitivity", label="local sensititivity w.r.t. Ks", xlab="S [mol/L]", ylab="Sensitivity")

# ╔═╡ aef87a3e-fa60-4666-8517-8c53fb797402
plots["mm_sens_Ks_rel"] = plot(sens_mm_Ks_rel, 0, 100, lw=2, title="Michaelis-Menten sensitivity (relative)", label="local sensititivity w.r.t. Ks", xlab="S [mol/L]", ylab="Sensitivity")

# ╔═╡ a87311b2-bc02-48a8-9dea-11a0bcbb5377
plots["mm_sens_vmax"] = plot(sens_mm_mumax, 0, 100, lw=2, title="Michaelis-Menten sensitivity", label="local sensititivity w.r.t. vmax", xlab="S [mol/L]", ylab="Sensitivity")

# ╔═╡ af4aa315-1411-41f1-bd01-f342c9d51bc2
plots["mm_sens_vmax_rel"] = plot(sens_mm_mumax_rel, 0, 100, lw=2, title="Michaelis-Menten sensitivity (relative)", label="local sensititivity w.r.t. vmax", xlab="S [mol/L]", ylab="Sensitivity", ylims=(0, 2))

# ╔═╡ 7e7956b5-43fc-48fb-92ff-67d0fdb875b1
plots["mm_sens"] = plot(plots["mm_sens_Ks"], plots["mm_sens_vmax"], plots["mm_sens_Ks_rel"], plots["mm_sens_vmax_rel"], size=(1000, 800))

# ╔═╡ 71e3de0e-03d6-491b-9a1d-071e03f1f6ef
plots["cuminf"] = contourf(β_vals, γ_vals, (g, l)->sir_summary((g, l))[1], color=:speed,
	xlab="β", ylab="γ", title="Cumulative number of infected")

# ╔═╡ 164e1b01-18c8-4c23-b5a9-20c6c8d74ea8
plots["maxinf"] = contourf(β_vals, γ_vals, (g, l)->sir_summary((g, l))[2], color=:speed,
		xlab="β", ylab="γ", title="Maximum number of infected")

# ╔═╡ 2a10a3cc-cc24-4a62-9090-7ba424024706
plots["sir_heatmaps"] = plot(plots["cuminf"], plots["maxinf"], size=(800, 400))

# ╔═╡ b964e1e9-1dd2-4029-a4b9-66944d2c39ba
plots["sobol_sir_fo"] = groupedbar(["cumulative I", "maximum I"], sobol.S1, label=["β" "γ"], title="Sobol SIR first order", ylab="variance")

# ╔═╡ f45583a1-5872-4500-b6ca-1f8b3a423bdb
plots["sobol_sir_tot"] = groupedbar(["cumulative I", "maximum I"], sobol.ST, label=["β" "γ"], title="Sobol SIR total effects", ylab="variance")

# ╔═╡ 9c7ef318-38b7-49d5-a203-df0dccbe1e57
plots["morris_sir_mu"] = groupedbar(["cumulative I", "maximum I"], morris.means, label=["β" "γ"], title="Morris SIR μ")

# ╔═╡ 991e84f5-6711-45a7-b8ab-1e8199b314b9
plots["morris_sir_mu_star"] =groupedbar(["cumulative I", "maximum I"], morris.means_star, label=["β" "γ"], title="Morris SIR μ*")

# ╔═╡ 8866b3ee-de5d-448f-a995-9a4743cb511d
plots["morris_sir_var"] = groupedbar(["cumulative I", "maximum I"], morris.variances, label=["β" "γ"], title="Morris SIR Var", yscale=:log10)

# ╔═╡ e5e0a0db-2662-4886-9621-8e0032cb6e42
plots

# ╔═╡ Cell order:
# ╠═1f5c60e2-e3ec-41ba-9c3e-44083987d710
# ╠═58eb5c0f-557d-42a1-ae5d-728de0c498bd
# ╠═a5cbf15c-490e-4d63-9a99-1c579be36c7e
# ╠═a981929c-2a53-11ef-16ce-f3918dd88313
# ╠═9629b2f8-2c21-40f1-a48d-a206f8bad479
# ╠═288c20e5-f1b2-4fbc-bc52-c1ec0e9d8422
# ╠═2100c8b0-0f88-4228-a7f9-dc62f4a3cf9c
# ╠═48856eee-8568-483a-bd0f-f3fb42c58906
# ╠═3ba4a491-c64d-48de-b263-83679b230b32
# ╠═6f6ee7e0-3413-43c5-a9bd-132267a96d82
# ╠═36eb6a4d-a167-47eb-9208-52d34949f37e
# ╠═e887a9ab-a4d8-4b05-a3f8-a8165d517a8a
# ╠═733ae8b7-b3d6-443a-b4a1-f8fdbde0e0b6
# ╠═9813f4f9-1ae0-4480-bf9d-14ab547d01c2
# ╠═9bed9078-0c22-4f29-96db-4274999adf9c
# ╠═34f47469-63b7-4e9a-98e8-e67827e4dbf4
# ╠═e34de9d2-6189-46fa-88a6-a37f94294384
# ╠═80eeefec-1fd5-4880-9209-df27375f9e1e
# ╠═7305bd7b-e823-451d-9eb1-15f5229fe431
# ╠═bdf2dd3b-55c8-445f-b9ed-66805014aa47
# ╠═ea341f0b-b081-4eda-b7d8-5ca616261ee8
# ╠═6ed4121b-2b68-4446-bf24-e8ce61661f84
# ╠═2645bb50-bec6-4b5e-b423-9cbc2f4dbd24
# ╠═ef6cb5d2-aa48-459d-931f-a5a79fcf2507
# ╠═02209c47-4ff3-4823-a998-fe1f173f749a
# ╠═5075c78a-56e3-4372-b605-de3dfec5cebc
# ╠═9b048268-30d7-4eaf-8a26-ab2e1c0404cd
# ╠═4f1cc2c6-254d-4f9c-b97a-09d22cb3838a
# ╠═118f4f61-12fb-4025-bd01-91a7fa6f9bd7
# ╠═e6bfb7d1-e594-41c9-9ee2-9668cae70e71
# ╠═ccf61cf5-33b2-4349-a422-16a7c883038d
# ╠═5a9da0cd-9809-445e-aeb3-b91875064106
# ╠═ae807e33-3e05-4a3d-a591-2b0e04c9d6d8
# ╠═260ba845-ae2f-4e9f-a2d0-f43c64648461
# ╠═8cf4e046-c979-4633-b922-43a46c914f3c
# ╠═a3f82cf6-b25d-4de1-b264-93d6e047ec6d
# ╠═090f89e0-d1e2-4363-830e-4d0635749890
# ╠═660c25d8-bcc5-4de3-9d04-013356acd22a
# ╠═06929445-f0ac-430f-a09f-1743f2de6f91
# ╠═11c32674-9817-4b0c-bcb2-235c431a5519
# ╠═1ed87d81-4b19-4de3-9735-fc36894b64c7
# ╠═633b06bf-ad63-418d-bf2d-899c51eee807
# ╠═42908d09-768d-4004-b304-cb136ed5b791
# ╠═449b2aa3-9efe-4399-8123-95fd952cb062
# ╠═0bb3f53f-50a6-40d3-95ac-be64a8dc79da
# ╠═10cafa44-8089-4aae-bd27-686be0496dc4
# ╠═305b0b4a-4051-409e-9c64-f95ab1778a60
# ╠═718a957c-74e6-455d-b847-e3d10e7f2016
# ╠═20af2948-d688-45ee-804c-26e35a59162c
# ╠═95d7287c-8008-4cad-9af6-ed3e77b9db41
# ╠═4c23f7ac-bf2e-42bd-aa2d-75d98cb691f1
# ╠═611d9607-99e6-488a-95f2-18da8fb8759e
# ╠═34cfa9fc-d197-40aa-9987-4736493e533e
# ╠═7c5a9bca-a568-43d5-88f1-313bcf19bc96
# ╠═bb1f08d7-d551-4226-9873-21fc903ed10c
# ╠═57107079-7baa-4bce-b801-b50d1dfc7a5f
# ╠═77b7b12f-1b91-4e78-af25-7de8996be2b7
# ╠═136e14d5-2330-4e24-8085-4f13fa3b9e26
# ╠═e436ef36-1f14-48fd-93af-aec0eadc51ed
# ╠═4c6849ca-723d-44da-b345-8373e6c67284
# ╠═9563199f-257b-4f89-a4c3-7cd1cc762bac
# ╠═93f93512-8337-40cd-8697-e46de3fcd945
# ╠═f6aab879-37ad-4425-a5a5-82db2ca4d188
# ╠═bffa2e06-271f-4ddb-a6bf-cdf9903b5c31
# ╠═3cf65781-9369-4219-9520-88cfd38adf0e
# ╠═32cdba4c-84b9-4449-8a68-4a568592dce7
# ╠═aef87a3e-fa60-4666-8517-8c53fb797402
# ╠═a87311b2-bc02-48a8-9dea-11a0bcbb5377
# ╠═af4aa315-1411-41f1-bd01-f342c9d51bc2
# ╠═7e7956b5-43fc-48fb-92ff-67d0fdb875b1
# ╠═415cbdec-71f3-4ac7-8479-14d0c7f6481c
# ╠═95b84eb0-d3e2-4cc3-882b-d3a4e82c1c4f
# ╠═0ac6b303-183b-4256-8436-052071d989cd
# ╠═71e3de0e-03d6-491b-9a1d-071e03f1f6ef
# ╠═164e1b01-18c8-4c23-b5a9-20c6c8d74ea8
# ╠═2a10a3cc-cc24-4a62-9090-7ba424024706
# ╠═d8818694-9270-4882-bc50-2eff71ea1c7a
# ╠═3f6491cf-d9d2-4be0-a341-db1a1e891029
# ╠═6d1e0bc0-4b28-4353-8389-c5ff4d156b05
# ╠═2d58a23a-7354-4e07-8cbb-ec1351051845
# ╠═5d2ffa7b-1202-45b4-b3a7-536ea38e2bfc
# ╠═8c13eb13-4b54-40ab-8408-ede4f72d7af4
# ╠═b964e1e9-1dd2-4029-a4b9-66944d2c39ba
# ╠═55b7b2ca-0e80-42a0-9aea-12b46947305c
# ╠═f45583a1-5872-4500-b6ca-1f8b3a423bdb
# ╠═f4869054-c5c6-4028-b9f3-0bb3e6b4305a
# ╠═9e4b6c92-03a9-4ddb-b4f1-97dddf1b7a30
# ╠═14676d4b-a43a-4e4f-8eb1-d93c098df156
# ╠═1478f40c-f530-48de-9608-bdd24cab0f23
# ╠═b455c7c0-d9d5-402f-b49d-fcd1f4b3626a
# ╠═9c7ef318-38b7-49d5-a203-df0dccbe1e57
# ╠═991e84f5-6711-45a7-b8ab-1e8199b314b9
# ╠═8866b3ee-de5d-448f-a995-9a4743cb511d
# ╟─8928645a-c20e-407c-8ba0-92f3f27206f0
# ╠═3902e862-3bc9-4786-bdc5-0415bc17dbd7
# ╠═30d0f3e0-a6c0-44c8-bcfd-7b163871b4f0
# ╠═e5e0a0db-2662-4886-9621-8e0032cb6e42
