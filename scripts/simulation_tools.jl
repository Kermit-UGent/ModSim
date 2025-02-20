### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ eb128248-dfbc-11ee-2efb-cb60fad76024
# ╠═╡ skip_as_script = true
#=╠═╡
begin
    using Pkg
	Pkg.activate("..")
end
  ╠═╡ =#

# ╔═╡ 4f364a47-c45e-4264-9a00-fb705ff3f169
using Plots, PlutoUI, LaTeXStrings

# ╔═╡ fd2f6fb4-d6da-4b07-8044-d3e2a09a6b4d
using Catalyst, DifferentialEquations

# ╔═╡ 6e94632a-cad9-49ea-8cdc-e4ec55871682
using Distributions

# ╔═╡ eb30d4b4-ced7-473a-a5df-5ebed6a1c357
using Symbolics

# ╔═╡ 0b52aede-1a71-4034-9921-2ffc9120c0ae
using Latexify

# ╔═╡ ff86545c-64a2-4c6f-8c36-9066b270aa6b
md"""
## Events and callbacks
"""

# ╔═╡ 018d215a-41f4-4d21-b57f-fca3ac3a755d
md"### Bouncing ball"

# ╔═╡ 992d2b4f-6521-4368-910d-3e7ef07fb6df
function ball!(du, u, g, t)
	y, v = u
	du[1] = v
	du[2] = -g
	return du
end

# ╔═╡ 89d089b2-3220-4696-9ded-364c318b4a1e
md"### Dosed bioreactor"

# ╔═╡ 21e810e7-0879-4ffb-84f9-d67dafb900d6
bacterial_growth = @reaction_network begin
	@species X(t)=10 G(t)=8
	@parameters r=0.2 m=0.8
	r, X + G --> 2X
	m, X --> 0
end

# ╔═╡ 7b3ca635-5a3c-43f2-9690-5cdbd9074ed2
# ╠═╡ disabled = true
#=╠═╡
let
	dosetimes = 5:5:20
	affect!(integrator) = integrator.u[2] += 10
	cb = PresetTimeCallback(dosetimes, affect!)
	
	prob = ODEProblem(bacterial_growth, [], (0, 20))
	sol = solve(prob, Tsit5(), callback=cb)
	plots["undosed bioreactor"] = plot(solve(prob, Tsit5()), lw=2,
						title="Undosed bioreactor")
	plots["dosed_bioreactor"] = plot(sol, lw=2, title="Dosed bioreactor")
end
  ╠═╡ =#

# ╔═╡ 4c7feb9f-a039-48a8-a7d0-3633f1a635fe
dosetimes = 5:5:20

# ╔═╡ 5ae78d7c-fc07-45f3-a59b-eb032ff2974d
timed_feeding = [dosetimes] => [bacterial_growth.G ~ bacterial_growth.G + 10]

# ╔═╡ 90a998c4-5fdd-42c3-9086-267a46fe3999
#@named reactor_dosed = ReactionSystem(equations(bacterial_growth); discrete_events=timed_feeding)

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

# ╔═╡ 22ff7a8b-fbae-4c5a-86c7-92d9c4eb7935
plot(u->u^2 - u^3, 0, 1, xlab="u", ylab="u'", lw=2)

# ╔═╡ b446a764-b80a-43c6-9fd1-72b099cfed0d
radicals = @reaction_network begin
	@species A(t)=10.0 B(t)=0 C(t)=0
	0.04, A --> B
	3e7, B + B --> C + B
	1e4, B + C --> A + C
end

# ╔═╡ c90f5869-89f8-486e-a4ad-9e4fbbd1eccd
#latexify(radicals, form=:ode) |> println

# ╔═╡ 8372fad0-07ea-4905-8b8e-26316aa93f68


# ╔═╡ 312e9d28-ed88-47ae-92b1-b1334e2e9d7a
plot(solve(ODEProblem(radicals, [], (0.0, 10)), Tsit5()))

# ╔═╡ 2cba2928-c7a7-4304-97ab-2efb5d6a8297
@time solve(ODEProblem(radicals, [], (0.0, 10)), Tsit5()) |> length

# ╔═╡ 2aad9bca-736f-4aaa-ad50-74cc01661117
plot(solve(ODEProblem(radicals, [], (0.0, 10)), Rosenbrock23()))

# ╔═╡ 80404482-0534-43fd-8783-84ae6522706d
@time solve(ODEProblem(radicals, [], (0.0, 10)), Rosenbrock23())  |> length

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
#println(latexify(brownian_motion))

# ╔═╡ cba57c64-a43f-4a4c-a793-86358131ef70
competition_model = @reaction_network begin
	@species A(t)=2.0 B(t)=2.5
	@parameters r=0.5 K=100 d=0.1
	r * (1 - (A+B) / K), A --> 2A
	r * (1 - (A+B) / K), B --> 2B
	d, (A, B) --> ∅
end

# ╔═╡ bc5adc9a-6f9c-4d55-ae8f-9bb8af5ff644
convert(ODESystem, competition_model)

# ╔═╡ ffa88278-4f16-4f4a-9b12-9a475254dc63
oprob_comp = ODEProblem(competition_model, [], (0., 50.))

# ╔═╡ bc48b429-6382-47a9-91d5-02a9941f5c79
sprob_comp = SDEProblem(competition_model, [], (0., 50.))

# ╔═╡ 23f1da8d-7a67-4dab-8928-f1f5480fadd0
sol_ode = solve(oprob_comp)

# ╔═╡ a72a568b-17c8-4350-b719-d72abbb10e84
sol_sde = solve(sprob_comp);

# ╔═╡ 6826363e-d807-4c52-a26a-018e5c6868c9
eprob = EnsembleProblem(sprob_comp)

# ╔═╡ 26f8a52c-3f8d-4928-9f33-a15995323cb5
comp_ensemble = solve(eprob; trajectories = 20)

# ╔═╡ 3604baf4-7623-4842-aeb8-74785417b398
e_sumary = EnsembleAnalysis.EnsembleSummary(comp_ensemble)

# ╔═╡ 4f671499-19b2-426c-a689-2ee4b084f3d0
md"## Discrete stochastic differential equations"

# ╔═╡ 2ee4fe85-feed-471d-86e2-53212fc05650
0.1/log(2)

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

# ╔═╡ 5c3c6af6-add6-4f1a-91ca-0251b43d5e1a
plot(solve(sir_oprob), lw=2)

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

# ╔═╡ aec510c8-6f67-4bad-ac60-db8081dfc7f2
md"## Getting derivatives"

# ╔═╡ 1e7c7e9c-b286-48a2-920b-4ac2fe7e9533
# ╠═╡ skip_as_script = true
#=╠═╡
f(x) = log(x) + sin(x)^2 / x
  ╠═╡ =#

# ╔═╡ 459fe27d-bc53-469d-919d-724a5df98fe9
diff_complstep(f, x; h=1e-10) = imag(f(x+im*h)) / h

# ╔═╡ f61a7f2d-faa2-4a54-bd4a-c2fac6da888f
# ╠═╡ skip_as_script = true
#=╠═╡
g(x) = sin(cos(exp(x)) + x^2)
  ╠═╡ =#

# ╔═╡ 615370c3-4986-454e-85c5-11d447b51ff1
@variables x

# ╔═╡ 6ae135a1-72b5-434c-89e6-a0df76338d46
# ╠═╡ skip_as_script = true
#=╠═╡
g(x)
  ╠═╡ =#

# ╔═╡ 77e9ad42-87fb-4292-bb89-3b0cb271f1a2
# ╠═╡ skip_as_script = true
#=╠═╡
f(x)
  ╠═╡ =#

# ╔═╡ e0f9312f-dae9-49ea-9618-e786de59921f
# ╠═╡ skip_as_script = true
#=╠═╡
diff_fordiff(f, x; h=1e-10) = (f(x+h) - f(x)) / h
  ╠═╡ =#

# ╔═╡ 05bd415e-09db-4a19-8ccb-922096f218fc
diff_centrdiff(f, x; h=1e-10) = (f(x+h) - f(x-h)) / 2h

# ╔═╡ 52082e32-8c8f-4699-a44a-2fbe99b96319
a = 2

# ╔═╡ 5bbf9db9-4537-416d-9edc-6487961b1690
Dx = Differential(x) 

# ╔═╡ a99a342d-5c69-411d-8bd6-55ed6f14f438
#=╠═╡
Dx(f(x))
  ╠═╡ =#

# ╔═╡ 8b36f54a-420a-40ae-b655-d25c99fe7182
# ╠═╡ skip_as_script = true
#=╠═╡
df_sym = expand_derivatives(Dx(f(x)))  # this expands the derviatve operator
  ╠═╡ =#

# ╔═╡ 856890fb-a420-4ae0-8f6d-c40763d11cc2
# ╠═╡ skip_as_script = true
#=╠═╡
df = build_function(df_sym, x) |> eval  #builds an expression and turns it into a function
  ╠═╡ =#

# ╔═╡ 9d84b996-a3e2-48c1-ae51-2bfc343bff0f
# ╠═╡ skip_as_script = true
#=╠═╡
df(a)
  ╠═╡ =#

# ╔═╡ 9be6aa61-46fe-4d53-941b-6cb0b21133a0
#=╠═╡
d2gdt2 = Dx(Dx(g(x))) |> expand_derivatives |> simplify
  ╠═╡ =#

# ╔═╡ 4aa42d80-81b5-4144-8b4a-6c3ebf29d84b
# ╠═╡ skip_as_script = true
#=╠═╡
dgdt = Dx(g(x)) |> expand_derivatives |> simplify
  ╠═╡ =#

# ╔═╡ 6f6ba9c5-bc87-4b71-82a0-a4c64750a2ab
md"## Appendix"

# ╔═╡ 2bf198cd-4463-4f97-975e-f805a37fa780
TableOfContents()

# ╔═╡ 13036978-3008-4d7e-9661-a967381f4db6
plots = Dict()

# ╔═╡ 0de4a715-89d6-403b-ac65-adaed6062371
let
	
	function condition(u, t, integrator)
		u[1]  # check when u[1] (i.e. x) == 0
	end
	
	function affect!(integrator)
		# nearly elastic collision
		integrator.u[2] = -0.9integrator.u[2]
	end
	
	cb = ContinuousCallback(condition, affect!)
	
	u0 = [50.0, 0.0]
	tspan = (0.0, 15.0)
	g = 9.81
	prob = ODEProblem(ball!, u0, tspan, g)
	sol = solve(prob, Tsit5(), callback = cb)
	plots["bouncing_ball"] = plot(sol, lw=2, label=[L"y(t)" L"v(t)"], title="Bouncing ball with callbacks")
end

# ╔═╡ 5e7453ce-721b-4e4e-8f75-bb69aac8cf42
let
	@unpack G, X = bacterial_growth
	dosing = [5, 10, 15, 20] => [G ~ G + 10]

	@named dosed_reactor = ReactionSystem(equations(bacterial_growth),
				discrete_events=dosing)

	dosed_reactor = complete(dosed_reactor)

	oprob = ODEProblem(dosed_reactor, [], (0, 20))
	
	sol = solve(oprob, Tsit5())

	prob = ODEProblem(bacterial_growth, [], (0, 20))

	plots["undosed bioreactor"] = plot(solve(prob, Tsit5()), lw=2,
						title="Undosed bioreactor")

	plots["dosed_bioreactor"] = plot(sol, lw=2, title="Dosed bioreactor")

	#R = G / X
	#plot(sol, idxs=[R])
	
end

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
		title="Five draws from a Wiener process", label="")
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

# ╔═╡ dc583986-3b76-44a7-ad5a-909aa392c38b
let
	τ = 0.01  # stepsize
	tsteps = 0:τ:15
	X = randn(length(tsteps), 10)
	W = cumsum(X, dims=1)
	W .-= X[[1], :]  # start at 0
	W .*= √(τ)
	p = plot(W[:,1:2:end], W[:,2:2:end], lw=2, xlab=L"x", ylab=L"y",
		title="Five draws from a Wiener process (2D)", label="", alpha=0.8, aspect_ratio=:equal)
	plots["Wiener_2D"] = p
	p
end

# ╔═╡ f3ec3c15-aa56-4b20-ad83-2acec0aad6b6
let
	τ = 0.02  # stepsize
	tsteps = 0:τ:20
	n=5
	X = randn(length(tsteps), 15)
	W = cumsum(X, dims=1)
	W .-= X[[1], :]  # start at 0
	W .*= √(τ)
	
	p = plot(lw=2, xlab=L"x", ylab=L"y", zlab=L"z",
		title="Five draws from a Wiener process (3D)", label="", alpha=0.8, aspect_ratio=:equal)
	for i in 1:3:15
		plot3d!(W[:,i], W[:,i+1], W[:,i+2], alpha=0.8, label="")
	end
	plots["Wiener_3D"] = p
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

# ╔═╡ ce7b1335-e43d-40dc-8ee8-d75b6df80dce
plots["competition_ensemble"] = plot(plot(comp_ensemble, idxs=:A, title="Species A"), plot(comp_ensemble, idxs=:B , title="Species B"), layout=(2,1))

# ╔═╡ 6286a56e-90c4-41b0-8b62-f7189fd5505d
plots["competition_ensemble_summary"] = plot(e_sumary, xlab=:t, title="Competition ensemble summary")

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

# ╔═╡ 6e433b23-3915-466a-85ca-31948896e3cb
#=╠═╡
plots["example_diff"] = plot(f, 1, 5, label=L"f(x)", lw=2, xlab=L"x")
  ╠═╡ =#

# ╔═╡ 525ae0f8-026a-43a3-a8b4-3a101c8ad6a7
# ╠═╡ skip_as_script = true
#=╠═╡
if !ismissing(f(a))
plot(f, 1, 10, label="\$f(x)\$", xlabel="\$x\$", lw=2)
plots["fdiff_example"] = plot!(df, 1, 10, label="\$f'(x)\$", lw=2)
end
  ╠═╡ =#

# ╔═╡ ae5be024-863d-46f1-b5d9-fc691b10e63b
#=╠═╡
let
	fexamp(x) = 64x*(1-x)*(1-2x)^2*(1-8x+8x^2)^2
	#dfexamp = diff(fexamp(x), x)
	dfexamp = build_function(expand_derivatives(Dx(fexamp(x))), x) |> eval
	error(diff, h; x=1.0) = max(abs(Float64(dfexamp(x)) - diff(fexamp, x, h=h)), 1e-50)
	stepsizes = map(t->10.0^t, -20:0.1:-1);
	p = plot(stepsizes, error.(diff_fordiff, stepsizes), label="forward difference",
    xscale=:log10, yscale=:log10, lw=2, legend=:bottomright)
	plot!(stepsizes, error.(diff_centrdiff, stepsizes), label="central difference", 		lw=2)
	plot!(stepsizes, error.(diff_complstep, stepsizes), label="complex step", lw=2)
	xlabel!("\$h\$")
	ylabel!("absolute error")
	plots["numdiff_error"] = p
end
  ╠═╡ =#

# ╔═╡ 23c38172-8916-42e3-a186-5cdd16939b5b
plots

# ╔═╡ Cell order:
# ╠═eb128248-dfbc-11ee-2efb-cb60fad76024
# ╠═4f364a47-c45e-4264-9a00-fb705ff3f169
# ╠═fd2f6fb4-d6da-4b07-8044-d3e2a09a6b4d
# ╠═ff86545c-64a2-4c6f-8c36-9066b270aa6b
# ╠═018d215a-41f4-4d21-b57f-fca3ac3a755d
# ╠═992d2b4f-6521-4368-910d-3e7ef07fb6df
# ╠═0de4a715-89d6-403b-ac65-adaed6062371
# ╟─89d089b2-3220-4696-9ded-364c318b4a1e
# ╠═21e810e7-0879-4ffb-84f9-d67dafb900d6
# ╠═5e7453ce-721b-4e4e-8f75-bb69aac8cf42
# ╠═7b3ca635-5a3c-43f2-9690-5cdbd9074ed2
# ╠═4c7feb9f-a039-48a8-a7d0-3633f1a635fe
# ╠═5ae78d7c-fc07-45f3-a59b-eb032ff2974d
# ╠═90a998c4-5fdd-42c3-9086-267a46fe3999
# ╟─4c420cd0-975a-4d3b-a6cc-9eca35490c62
# ╠═081852cd-b589-4891-83c0-f19e4da2eb35
# ╠═39ee1679-23ad-4448-ae20-b55345c950b4
# ╠═43183390-219c-4ade-a92a-80373a7c8585
# ╟─e22ed4e6-6420-490c-9d55-fee19069f534
# ╟─f2201e7f-9aec-4962-a9b7-a7fd627d64e6
# ╟─8fcdad45-8789-473e-9d1a-2692b62d526d
# ╠═22ff7a8b-fbae-4c5a-86c7-92d9c4eb7935
# ╠═b446a764-b80a-43c6-9fd1-72b099cfed0d
# ╠═c90f5869-89f8-486e-a4ad-9e4fbbd1eccd
# ╠═8372fad0-07ea-4905-8b8e-26316aa93f68
# ╠═312e9d28-ed88-47ae-92b1-b1334e2e9d7a
# ╠═2cba2928-c7a7-4304-97ab-2efb5d6a8297
# ╠═2aad9bca-736f-4aaa-ad50-74cc01661117
# ╠═80404482-0534-43fd-8783-84ae6522706d
# ╠═945febc8-2954-4d5c-a4b1-7f95c70ef4b6
# ╟─5099fd5b-735e-4de3-918a-685dd4a82c22
# ╟─18def41b-0827-4356-adfa-d3fc63e23cfa
# ╟─dc583986-3b76-44a7-ad5a-909aa392c38b
# ╠═f3ec3c15-aa56-4b20-ad83-2acec0aad6b6
# ╠═6501070e-9093-4928-89b9-b9dd34128810
# ╠═6c5f1ee1-ece0-4a44-8e92-c375f2bd40a7
# ╠═7d0a5294-6fc2-45d5-b04f-2824649f5d81
# ╠═3d549763-2ba2-48b6-b486-2ddd462977dd
# ╠═cba57c64-a43f-4a4c-a793-86358131ef70
# ╠═bc5adc9a-6f9c-4d55-ae8f-9bb8af5ff644
# ╠═ffa88278-4f16-4f4a-9b12-9a475254dc63
# ╟─98b16b2b-e792-41ef-a601-08289236af98
# ╠═bc48b429-6382-47a9-91d5-02a9941f5c79
# ╟─23b18b26-18de-4e08-abcf-ae5b5ac5af3d
# ╟─53f1b6b4-a1e8-4e4d-b936-eebc5c37a712
# ╟─b90e3157-3aed-4371-a396-283305a4d5b5
# ╠═23f1da8d-7a67-4dab-8928-f1f5480fadd0
# ╠═a72a568b-17c8-4350-b719-d72abbb10e84
# ╠═6826363e-d807-4c52-a26a-018e5c6868c9
# ╠═26f8a52c-3f8d-4928-9f33-a15995323cb5
# ╠═ce7b1335-e43d-40dc-8ee8-d75b6df80dce
# ╠═3604baf4-7623-4842-aeb8-74785417b398
# ╠═6286a56e-90c4-41b0-8b62-f7189fd5505d
# ╟─4f671499-19b2-426c-a689-2ee4b084f3d0
# ╠═6e94632a-cad9-49ea-8cdc-e4ec55871682
# ╠═bb562e40-daaf-47be-a1d0-1938bc227fc0
# ╠═2ee4fe85-feed-471d-86e2-53212fc05650
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
# ╠═5c3c6af6-add6-4f1a-91ca-0251b43d5e1a
# ╠═c9fb3ebb-e37a-42dc-85f4-87b4fc5ca492
# ╠═3122f1ec-7c47-4873-bd6c-94f1a8616088
# ╠═486c9c19-e2ed-4632-be5d-b76dadf433c4
# ╠═26c3f4e0-dc66-43cc-9215-94e5c73531d8
# ╠═04eeeae2-f525-4717-a397-39b895fa80eb
# ╠═7055aa36-e4d8-42e6-a92c-f7d8ec5aa960
# ╠═bfbccda0-8648-435f-8b12-13957c574f54
# ╠═6e433b23-3915-466a-85ca-31948896e3cb
# ╠═a99a342d-5c69-411d-8bd6-55ed6f14f438
# ╠═aec510c8-6f67-4bad-ac60-db8081dfc7f2
# ╠═1e7c7e9c-b286-48a2-920b-4ac2fe7e9533
# ╠═6ae135a1-72b5-434c-89e6-a0df76338d46
# ╠═459fe27d-bc53-469d-919d-724a5df98fe9
# ╠═77e9ad42-87fb-4292-bb89-3b0cb271f1a2
# ╠═eb30d4b4-ced7-473a-a5df-5ebed6a1c357
# ╠═525ae0f8-026a-43a3-a8b4-3a101c8ad6a7
# ╠═f61a7f2d-faa2-4a54-bd4a-c2fac6da888f
# ╠═615370c3-4986-454e-85c5-11d447b51ff1
# ╠═e0f9312f-dae9-49ea-9618-e786de59921f
# ╠═8b36f54a-420a-40ae-b655-d25c99fe7182
# ╠═9be6aa61-46fe-4d53-941b-6cb0b21133a0
# ╠═856890fb-a420-4ae0-8f6d-c40763d11cc2
# ╠═05bd415e-09db-4a19-8ccb-922096f218fc
# ╠═0b52aede-1a71-4034-9921-2ffc9120c0ae
# ╠═9d84b996-a3e2-48c1-ae51-2bfc343bff0f
# ╠═4aa42d80-81b5-4144-8b4a-6c3ebf29d84b
# ╟─ae5be024-863d-46f1-b5d9-fc691b10e63b
# ╠═52082e32-8c8f-4699-a44a-2fbe99b96319
# ╠═5bbf9db9-4537-416d-9edc-6487961b1690
# ╟─6f6ba9c5-bc87-4b71-82a0-a4c64750a2ab
# ╠═2bf198cd-4463-4f97-975e-f805a37fa780
# ╠═13036978-3008-4d7e-9661-a967381f4db6
# ╠═23c38172-8916-42e3-a186-5cdd16939b5b
