### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ eb128248-dfbc-11ee-2efb-cb60fad76024
begin
    using Pkg
	Pkg.activate("..")
    using Catalyst, DifferentialEquations, Plots, PlutoUI, LaTeXStrings
end

# ╔═╡ 0b52aede-1a71-4034-9921-2ffc9120c0ae
using Latexify

# ╔═╡ 6e94632a-cad9-49ea-8cdc-e4ec55871682
using Distributions

# ╔═╡ 4c420cd0-975a-4d3b-a6cc-9eca35490c62
md"## Stiff ODEs"

# ╔═╡ 081852cd-b589-4891-83c0-f19e4da2eb35
combustion_model = @reaction_network begin
	@species u(t)=1/100
	u^2, u => 2u  # growth proportional with surface
	u^3, u => 0   # decay proportional with volume (oxygen use)
end

# ╔═╡ 39ee1679-23ad-4448-ae20-b55345c950b4
convert(ODESystem, combustion_model)

# ╔═╡ 43183390-219c-4ade-a92a-80373a7c8585
combustion_problem = ODEProblem(combustion_model, [:u=>1/100], (0.0, 200.0))

# ╔═╡ b446a764-b80a-43c6-9fd1-72b099cfed0d
radicals = @reaction_network begin
	@species A(t)=10.0 B(t)=10.0 C(t)=10.0
	0.04, A --> B
	3e7, B + B --> C + B
	1e4, B + C --> A + C
end

# ╔═╡ c90f5869-89f8-486e-a4ad-9e4fbbd1eccd
latexify(radicals, form=:ode) |> println

# ╔═╡ 8372fad0-07ea-4905-8b8e-26316aa93f68


# ╔═╡ 312e9d28-ed88-47ae-92b1-b1334e2e9d7a
plot(solve(ODEProblem(radicals, [], (0.0, 40.0)), Rosenbrock23()))

# ╔═╡ 2aad9bca-736f-4aaa-ad50-74cc01661117
plot(solve(ODEProblem(radicals, [], (0.0, 1e12)), Rosenbrock23()))

# ╔═╡ 945febc8-2954-4d5c-a4b1-7f95c70ef4b6
md"## Stochastic differential equations"

# ╔═╡ 6501070e-9093-4928-89b9-b9dd34128810
brownian_motion = @reaction_network begin
	@species A(t)=1
	@parameters r₁=2 r₂=1
	r₁, ∅ --> A
	r₂, A --> ∅
end

# ╔═╡ 6c5f1ee1-ece0-4a44-8e92-c375f2bd40a7
sprob_bm = SDEProblem(brownian_motion, [], (0., 20.))

# ╔═╡ 7d0a5294-6fc2-45d5-b04f-2824649f5d81
println(latexify(brownian_motion))

# ╔═╡ cba57c64-a43f-4a4c-a793-86358131ef70
competition_model = @reaction_network begin
	@species A(t)=2.0 B(t)=2.5
	@parameters r=0.5 K=100 d=0.1
	r * (1 - (A+B) / K), A --> 2A
	r * (1 - (A+B) / K), B --> 2B
	d, (A, B) --> ∅
end

# ╔═╡ 44142e83-1f6c-4644-996e-246817748b7e
println(competition_model)

# ╔═╡ 1924514f-0112-4ec8-830a-4eed54d81392
latexify(competition_model, form=:ode) |> println

# ╔═╡ a807b275-99fa-47c8-abc0-e9afcde445b8


# ╔═╡ bc5adc9a-6f9c-4d55-ae8f-9bb8af5ff644
convert(ODESystem, competition_model)

# ╔═╡ ffa88278-4f16-4f4a-9b12-9a475254dc63
oprob_comp = ODEProblem(competition_model, [], (0., 50.))

# ╔═╡ bc48b429-6382-47a9-91d5-02a9941f5c79
sprob_comp = SDEProblem(competition_model, [], (0., 50.))

# ╔═╡ 23f1da8d-7a67-4dab-8928-f1f5480fadd0
sol_ode = solve(oprob_comp);

# ╔═╡ a72a568b-17c8-4350-b719-d72abbb10e84
sol_sde = solve(sprob_comp);

# ╔═╡ 4f671499-19b2-426c-a689-2ee4b084f3d0
md"## Discrete stochastic differential equations"

# ╔═╡ ade6e868-48ff-4b54-88ca-d180f115ca39
function gillespie_sir(S₀, I₀, R₀, α, β, t_end)
	S, I, R = S₀, I₀, R₀
	t = 0.0
	timesteps = [t]
	states = [(;S, I, R)]
	while t < t_end
		# compute propensities
		a_i = α * S * I  # infection
		a_r = β * I      # recovery
		# total propensity
		Rtot = a_i + a_r
		Rtot == 0 && break
		# sample time step
		τ = rand(Exponential(1 / Rtot))
		# update
		if rand() < (a_i / Rtot)  # infection
			S -= 1
			I += 1
		else                 # recovery
			I -= 1
			R += 1
		end
		t += τ
		push!(states, (;S, I, R))
		push!(timesteps, t)
	end
	return timesteps, states
end		

# ╔═╡ a9ba703a-c9f6-4c9a-930c-3d5c977bd559
S₀ = 50

# ╔═╡ a7049066-b990-470f-804d-00adb9201122
I₀ = 2

# ╔═╡ bfe45dd0-d429-40ef-859f-f5f6e2752904
R₀ = 0

# ╔═╡ 1e212d4d-5722-4c29-b27d-5d8ce128a1dd
α = 0.004

# ╔═╡ addfa406-0a45-4ba1-a6e1-345ef987ca3e
β = 0.03

# ╔═╡ 71d37d78-7a92-4971-aa6f-098030d618e4
timesteps, sir_states = gillespie_sir(S₀, I₀, R₀, α, β, 100)

# ╔═╡ aa9089b2-38b8-4ab1-9b34-87ee5b017c1d
sir_model  = @reaction_network begin
	α, S + I --> 2I
	β, I --> R
end

# ╔═╡ 3343b335-54c6-42b6-b2e1-f78bdf791b30
first.(sir_states)

# ╔═╡ 1e2caf48-30b3-4094-a159-a94585e248b3
pars_sir = [:α=>α, :β=>β]

# ╔═╡ 5ac23911-376f-4ee7-ba33-a120b80e1020
u0_sir = [:S=>S₀, :I=>I₀, :R=>R₀]

# ╔═╡ 06cadbd3-01ae-4f68-ab1a-7e8bc4a7b00a


# ╔═╡ f57de319-5025-4ec8-a437-163b8ab0b94e
sir_oprob = ODEProblem(sir_model, u0_sir, (0.0, 100.0), pars_sir)

# ╔═╡ d0aa9ea8-38cb-4f00-92d2-58740ab2d32a
sir_dprob = DiscreteProblem(sir_model, u0_sir, (0.0, 100.0), pars_sir)

# ╔═╡ 88a00ea8-9c7b-49fc-8771-2324576a0124
sir_jprob = JumpProblem(sir_model, sir_dprob, Direct())

# ╔═╡ 9265ba3b-5979-4755-a7b8-cd422040de8c


# ╔═╡ 5c3c6af6-add6-4f1a-91ca-0251b43d5e1a
plot(solve(sir_oprob))

# ╔═╡ c9fb3ebb-e37a-42dc-85f4-87b4fc5ca492
plot(solve(sir_jprob, SSAStepper()))

# ╔═╡ 3122f1ec-7c47-4873-bd6c-94f1a8616088
convert(ODESystem, sir_model)

# ╔═╡ 486c9c19-e2ed-4632-be5d-b76dadf433c4
hcat(sir_states...)

# ╔═╡ 26c3f4e0-dc66-43cc-9215-94e5c73531d8
timesteps

# ╔═╡ 04eeeae2-f525-4717-a397-39b895fa80eb
sir_states

# ╔═╡ 7055aa36-e4d8-42e6-a92c-f7d8ec5aa960
second(x) = x[2]

# ╔═╡ 627418f5-33fe-4515-8a0e-7f79774814d8
let
	S_gp = first.(sir_states)
	I_gp = second.(sir_states)
	R_gp = last.(sir_states)
	p_gp = plot(timesteps, S_gp, label="S", xlab="t")
	plot!(timesteps, I_gp, label="I")
	plot!(timesteps, R_gp, label="R")
end

# ╔═╡ 13036978-3008-4d7e-9661-a967381f4db6
plots = Dict()

# ╔═╡ e22ed4e6-6420-490c-9d55-fee19069f534
let
	p = plot(solve(combustion_problem, Tsit5()), lw=2, title="Combustion model (general solver)")
	plots["combustion_Tsi5"] = p
	p
end

# ╔═╡ f2201e7f-9aec-4962-a9b7-a7fd627d64e6
let
	p = plot(solve(combustion_problem, Rosenbrock23()), lw=2, title="Combustion model (stiff solver)")
	plots["combustion_Rosenbrock"] = p
	p	
end

# ╔═╡ 8fcdad45-8789-473e-9d1a-2692b62d526d
let
	# variant starting close to the equilibrium
	combustion_problem2 = ODEProblem(combustion_model, [:u=>0.999], (0.0, 200.0))
	sol_combustion_nonstiff = solve(combustion_problem2, Tsit5())
	sol_combustion_stiff = solve(combustion_problem2, Rosenbrock23())
	p = plot(sol_combustion_nonstiff, label="non-stiff solver (Tsit5)", lw=2)
	plot!(sol_combustion_stiff, label="stiff (Rosenbrock23)", lw=2)
	plots["combustion_solvers"] = p
	p
end

# ╔═╡ 5099fd5b-735e-4de3-918a-685dd4a82c22
let
	τ = 0.01  # stepsize
	tsteps = 0:τ:10
	X = randn(length(tsteps), 5)
	W = cumsum(X, dims=1)
	W .-= X[[1], :]  # start at 0
	W .*= √(τ)
	p = plot(tsteps, W, lw=2, xlab=L"t", ylab=L"W(t)",
		title="Five draws form a Wiener process", label="")
	plots["Wiener"] = p
	p
end

# ╔═╡ 18def41b-0827-4356-adfa-d3fc63e23cfa
let
	τ = 0.001  # stepsize
	tsteps = 0:τ:5
	mask = 2 .≤ tsteps .≤ 3

	x = randn(length(tsteps))
	w = cumsum(x)
	w .-= x[1]
	w .*= √(τ)
	plarge = plot(tsteps[1:10:end], w[1:10:end])
	vspan!(plarge, [2, 3], alpha=0.3)
	psmall = plot(tsteps[mask], w[mask])
	p = plot(plarge, psmall,
		layout=(2,1), xlab=L"t",ylab=L"W(t)", lw=2, label="")
	plots["Wiener_scalefree"] = p
	p
end

# ╔═╡ 3d549763-2ba2-48b6-b486-2ddd462977dd
let
	p = plot(solve(sprob_bm), lw=2)
	plots["brownian_motion"] = p
	p
end

# ╔═╡ 98b16b2b-e792-41ef-a601-08289236af98
let
	p_ode = plot(sol_ode, lw=2, title="simulation competition ODE")
	plots["competition_ode"] = p_ode
	p_ode
end

# ╔═╡ 23b18b26-18de-4e08-abcf-ae5b5ac5af3d
let
	p_sde = plot(sol_sde, lw=2, title="Monte Carlo simulation competition SDE")
	plots["competition_sde"] = p_sde
	p_sde
end

# ╔═╡ 53f1b6b4-a1e8-4e4d-b936-eebc5c37a712
let
	p_sde = plot(solve(sprob_comp), lw=2, title="Monte Carlo simulation competition SDE (second throw)")
	plots["competition_sde2"] = p_sde
	p_sde
end

# ╔═╡ b90e3157-3aed-4371-a396-283305a4d5b5
let

	p = plot(sol_ode, idxs=(:A, :B), lw=2, label="ODE", title="Phase plot competition model", color="black", ls=:dash)
	plot!(p, sol_sde, idxs=(1, 2), lw=2, label="SDE", xlab="A", ylab="B", color="grey")
	scatter!([2], [2.5], color="red", label="x₀")
	plots["competition_phaseplot"] = p
	p
end

# ╔═╡ bb562e40-daaf-47be-a1d0-1938bc227fc0
begin
	n₀ = 20
	r = 0.1
	ns = [n₀]
	ts = [0.0]
	while last(ns) > 0
		t, n = last(ts), last(ns)
		τ = rand(Exponential(1/(n*r)))
		push!(ns, n-1)
		push!(ts, t+τ)
	end

	p_discrete_decay = scatter(ts, ns, xlab=L"t", ylab=L"n", label="", title="Random decay of $(n₀) particles")
	p_decay_times = bar(0:n₀-1, diff(ts), xlab=L"n_0-n", ylabel=L"\tau",label="", title="Event times")
	plot!(p_decay_times, inv.((n₀:-1:1) .* r), label=L"1/nr", lw=2)
	plots["discrete_decay"] = p_discrete_decay
	plots["discrete_decay_eventtimes"] = p_decay_times
end;

# ╔═╡ e2dd1f82-031b-46eb-ac8e-9d03172db770
p_discrete_decay

# ╔═╡ 38802384-48b0-4d29-aab3-2578527f995b
p_decay_times

# ╔═╡ aa3056b8-7140-4271-885d-8392c3634af6
let
	p = plot(solve(sir_jprob, SSAStepper()), lw=2, title="Discrete SIR model")
	sol_sir_ode = solve(sir_oprob)
	plot!(sol_sir_ode, idxs=:S, label="", alpha=0.5, color=:blue, ls=:dash, lw=2)
	plot!(sol_sir_ode, idxs=:I, label="", alpha=0.5, color=:orange, ls=:dash, lw=2)
	plot!(sol_sir_ode, idxs=:R, label="", alpha=0.5, color=:green, ls=:dash, lw=2)
	plots["Gillespie_SIR"] = p
	p
end

# ╔═╡ bfbccda0-8648-435f-8b12-13957c574f54
plots

# ╔═╡ Cell order:
# ╠═eb128248-dfbc-11ee-2efb-cb60fad76024
# ╟─4c420cd0-975a-4d3b-a6cc-9eca35490c62
# ╠═081852cd-b589-4891-83c0-f19e4da2eb35
# ╠═39ee1679-23ad-4448-ae20-b55345c950b4
# ╠═43183390-219c-4ade-a92a-80373a7c8585
# ╟─e22ed4e6-6420-490c-9d55-fee19069f534
# ╟─f2201e7f-9aec-4962-a9b7-a7fd627d64e6
# ╟─8fcdad45-8789-473e-9d1a-2692b62d526d
# ╠═b446a764-b80a-43c6-9fd1-72b099cfed0d
# ╠═c90f5869-89f8-486e-a4ad-9e4fbbd1eccd
# ╠═8372fad0-07ea-4905-8b8e-26316aa93f68
# ╠═312e9d28-ed88-47ae-92b1-b1334e2e9d7a
# ╠═2aad9bca-736f-4aaa-ad50-74cc01661117
# ╠═945febc8-2954-4d5c-a4b1-7f95c70ef4b6
# ╟─5099fd5b-735e-4de3-918a-685dd4a82c22
# ╟─18def41b-0827-4356-adfa-d3fc63e23cfa
# ╠═6501070e-9093-4928-89b9-b9dd34128810
# ╠═6c5f1ee1-ece0-4a44-8e92-c375f2bd40a7
# ╠═7d0a5294-6fc2-45d5-b04f-2824649f5d81
# ╠═3d549763-2ba2-48b6-b486-2ddd462977dd
# ╠═cba57c64-a43f-4a4c-a793-86358131ef70
# ╠═44142e83-1f6c-4644-996e-246817748b7e
# ╠═0b52aede-1a71-4034-9921-2ffc9120c0ae
# ╠═1924514f-0112-4ec8-830a-4eed54d81392
# ╠═a807b275-99fa-47c8-abc0-e9afcde445b8
# ╠═bc5adc9a-6f9c-4d55-ae8f-9bb8af5ff644
# ╠═ffa88278-4f16-4f4a-9b12-9a475254dc63
# ╟─98b16b2b-e792-41ef-a601-08289236af98
# ╠═bc48b429-6382-47a9-91d5-02a9941f5c79
# ╟─23b18b26-18de-4e08-abcf-ae5b5ac5af3d
# ╟─53f1b6b4-a1e8-4e4d-b936-eebc5c37a712
# ╟─b90e3157-3aed-4371-a396-283305a4d5b5
# ╠═23f1da8d-7a67-4dab-8928-f1f5480fadd0
# ╠═a72a568b-17c8-4350-b719-d72abbb10e84
# ╠═4f671499-19b2-426c-a689-2ee4b084f3d0
# ╠═6e94632a-cad9-49ea-8cdc-e4ec55871682
# ╠═bb562e40-daaf-47be-a1d0-1938bc227fc0
# ╟─e2dd1f82-031b-46eb-ac8e-9d03172db770
# ╠═38802384-48b0-4d29-aab3-2578527f995b
# ╠═ade6e868-48ff-4b54-88ca-d180f115ca39
# ╠═a9ba703a-c9f6-4c9a-930c-3d5c977bd559
# ╠═a7049066-b990-470f-804d-00adb9201122
# ╠═bfe45dd0-d429-40ef-859f-f5f6e2752904
# ╠═1e212d4d-5722-4c29-b27d-5d8ce128a1dd
# ╠═addfa406-0a45-4ba1-a6e1-345ef987ca3e
# ╠═71d37d78-7a92-4971-aa6f-098030d618e4
# ╠═627418f5-33fe-4515-8a0e-7f79774814d8
# ╠═aa9089b2-38b8-4ab1-9b34-87ee5b017c1d
# ╠═3343b335-54c6-42b6-b2e1-f78bdf791b30
# ╠═1e2caf48-30b3-4094-a159-a94585e248b3
# ╠═5ac23911-376f-4ee7-ba33-a120b80e1020
# ╠═06cadbd3-01ae-4f68-ab1a-7e8bc4a7b00a
# ╠═f57de319-5025-4ec8-a437-163b8ab0b94e
# ╠═d0aa9ea8-38cb-4f00-92d2-58740ab2d32a
# ╠═88a00ea8-9c7b-49fc-8771-2324576a0124
# ╟─aa3056b8-7140-4271-885d-8392c3634af6
# ╠═9265ba3b-5979-4755-a7b8-cd422040de8c
# ╠═5c3c6af6-add6-4f1a-91ca-0251b43d5e1a
# ╠═c9fb3ebb-e37a-42dc-85f4-87b4fc5ca492
# ╠═3122f1ec-7c47-4873-bd6c-94f1a8616088
# ╠═486c9c19-e2ed-4632-be5d-b76dadf433c4
# ╠═26c3f4e0-dc66-43cc-9215-94e5c73531d8
# ╠═04eeeae2-f525-4717-a397-39b895fa80eb
# ╠═7055aa36-e4d8-42e6-a92c-f7d8ec5aa960
# ╠═13036978-3008-4d7e-9661-a967381f4db6
# ╠═bfbccda0-8648-435f-8b12-13957c574f54
