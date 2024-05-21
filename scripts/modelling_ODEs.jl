### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 093b722d-28af-4219-8546-39a3262146b2
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	using Pkg
	Pkg.activate("..")
end
  ╠═╡ =#

# ╔═╡ a52da2c2-f7df-11ee-033f-8500edb3c03f
using Plots, PlutoUI, LaTeXStrings, Latexify

# ╔═╡ 8cff27a7-fde1-4b49-8ad6-513302997a4e
using Catalyst, DifferentialEquations

# ╔═╡ 0686fc66-5428-451f-aa72-c0250ad4bf67
using Symbolics

# ╔═╡ dfca2f9f-0134-461c-a18b-f66f2bf02943
md"# Modelling with ordinary differential equations"

# ╔═╡ 34bec0a1-40e8-48a2-9109-94872aaff1b9
md"## Balance equations"

# ╔═╡ abebedae-b977-43ae-aaa0-6b00990a5de4
md"""
### Water tank example

Water flows with a constant flow $q$ (in m$^3$/h) into a cylindrical with a height $H$ (in m) and a floor area $A$ (in m$^2$). Water leaves the tank with a rate that is proportional to the height of the water: 
$$q_\text{out}=rh(t)\,,$$
Describe the volume $V(t)$ (in m$^3$) using an ODE and solve it when you know that at $t=0$, the tank is empty.

In this system, we have a "conservation of water": the volume of water in the system is determined by the in-and outgoing flows. Let us consider what happens in a small time step $\Delta t$ and how this impacts a change in volume $\Delta V$:

$$\Delta V = q \Delta t - rh(t)\Delta t$$

We see that the change in water volume in a small time interval is determined by:
- the ingoing flow $q$, which is constant here;
- the outgoing flow, $rh(t)$, which depends on the height, which in turn depends on how full the tank is. Note that, for this geometry, $h(t)=V(t)/A$.

Putting things together and rearranging, we have:

$$\Delta V / \Delta t= q  - rV(t)/A\,.$$

If we take the limit of $\Delta t\rightarrow 0$, we obtain:

$$\frac{\mathrm{d}V(t)}{\mathrm{d}t}= q  - rV(t)/A\,.$$

This is a linear first-order differential equation. Its general solution is 

$$V(t)=Ce^{-r/At} + q/r\,,$$

and filling in the initial condition $V(0)=0$ allows us to obtain a specific solution to the initial value problem:

$$V(t)=qA/re^{-rt} + qA/r\,,$$
"""

# ╔═╡ 0d09ba2e-3cac-4051-b98f-26b79736b225
q = 10

# ╔═╡ fdc534e1-e334-48c1-aac4-5a89c47484e0
A = 1^2 * π

# ╔═╡ 104ad25e-743b-4f47-8b10-0d3d6715f95c
r = 0.2

# ╔═╡ 73eb7d0a-5433-4e3d-a008-748db66b8ef9
md"### Coffee example"

# ╔═╡ d3b84441-ed9f-436d-a690-660c5f4b8fbd
function coffee!(du, u, (q, Tmilk, Tenv, k), t)
	V, T = u
	du[1] = dV = q(t)
	du[2] = dT = q(t) * Tmilk - dV * T - (T - Tenv) * k 
	return du
end

# ╔═╡ ddd43577-eb2e-4c72-b829-d7195c165ddf
# nog niet correct, want Cp niet in rekening voor wet Newton

# ╔═╡ 10a26b97-8b0a-454e-afbf-141aef4aa04f
qin = t -> 2 ≤ t < 3 ? 5e-2 : 0.0

# ╔═╡ 206dbd4b-2ec0-4591-b7f9-d8e78d568c2f
coffee_prob = ODEProblem(coffee!, [1.5e-1, 80], (0.0, 8.0), (qin, 5, 20, 0.1
))

# ╔═╡ bf062835-538d-437a-bae6-6309c66ebd19
coffee_sol = solve(coffee_prob, saveat=0.1, tstops=1.9:0.001:3.3)

# ╔═╡ f3737612-5458-4c6e-a634-246b2cb8cb05
begin
	plot(qin, 0:0.01:8, ls=:dash, label="q", color=:red, lw=2, xlab="t", title="Debit of milk added")
	vspan!([2, 3], alpha=0.4, color=:pink, label="")
end

# ╔═╡ 071b8f85-75cc-487b-a13f-64428bab7147
md"## Law of mass action"

# ╔═╡ 8ae850c8-ec2c-4a6c-9c37-1c1f93bb56e9
function lma_illustration(n, m; ϵ=0.05, kwargs...)
	x1, y1 = rand(n), rand(n)
	x2, y2 = rand(m), rand(m)
	p = plot(aspect_ratio=:equal, xlim=(0,1), ylim=(0, 1),
			xticks=[], yticks=[];kwargs...)
	title!(p, "[A] = $n, [B] = $m")
	scatter!(x1, y1, label="A")
	scatter!(x2, y2, label="B", m=:^)
	D = (x1 .- x2').^2 .+ (y1 .- y2').^2 .|> sqrt
	n_events = count(<(ϵ), D) ÷ 2
	title!(p, "[A] = $n, [B] = $m\n $(n_events) reaction events")
	d = minimum(D, dims=2)
	θs = 0:0.1:2π
	for i in 1:n
		if d[i] < ϵ
			x, y = x1[i], y1[i]
			plot!(p, x.+ϵ.*cos.(θs), y.+ϵ.*sin.(θs),
					ls=:dash, color=:red, alpha=0.5, label="")
		end
	end
	return p
end

# ╔═╡ c01904be-7ee0-4b43-bd29-6fac1e38b3c6
rn = @reaction_network begin
	k1, A + B --> C
	k2, C --> B
	k3, A + A --> B + C
end

# ╔═╡ 689f9c80-5680-4cc5-a60d-bacd846f0925
species(rn)

# ╔═╡ ad511e07-230f-4405-9613-40f2e55676e6
parameters(rn)

# ╔═╡ 730255e3-ba36-482c-bff3-8435d348fe40
reactions(rn)

# ╔═╡ 959b8b64-959a-4869-b499-159eccfa0768
latexify(rn, type=:ode)

# ╔═╡ 8a752cbe-8220-41f8-997f-f2eb321b2a2b
convert(ODESystem, rn) |> latexify |> clipboard

# ╔═╡ 063f7b0d-c8e6-4128-87f9-d3b5b5139670
reactionsys = @reaction_network begin
	(r1, r1), 2NO <--> N2O2
	r2, N2O2 + H2 --> N2O + H2O
	r3, N2O + H2 --> N2 + H2O
end

# ╔═╡ 691b4a76-fc29-4cb9-a3ca-26d03bd84cec
latexify(reactionsys, kind=:ode) |> clipboard

# ╔═╡ 7c4e736b-fac7-47c1-b73c-ef7d46151e15
convert(ODESystem, reactionsys) |> latexify |> clipboard

# ╔═╡ 2b80012d-78c8-4cd3-9b8b-6d0177afa963
reactsyst_prob = ODEProblem(reactionsys,
		[:NO=>5.2, :H2O=>0, :H2=>5.1, :N2O2=>0, :N2O=>0, :N2=>0],  # initial values
		(0.0, 60.0),  # time interval
		[:r1=>1e2, :r2=>0.1, :r3=>50])  # parameter values

# ╔═╡ 7ad626ab-12a5-4a53-9991-7512de5092fb
#g = Graph(reactionsys)

# ╔═╡ ca2304c8-9a09-4063-9478-2207925d6444
reactionsys2 = @reaction_network begin
	r * NO^2 * H2, 2H2 + 2NO => 0
end

# ╔═╡ bbfd0201-c25e-446e-828e-6632a7e0c116
convert(ODESystem, reactionsys2) |> latexify |> clipboard

# ╔═╡ c77660bc-a478-45d4-8d8b-c5147d70b123
tank = @reaction_network begin
	@species V(t)=0
	@parameters q=1 A=0.5^2*pi
	q, 0 --> V  # incoming water
	r / A, V --> 0  # emptying
end

# ╔═╡ a5177857-d6e4-4004-92e1-34bbfb42de53
growth1 = @reaction_network begin
	@species B(t)=1
	@parameters r=log(2)
	r, B --> 2B
end

# ╔═╡ e6b6eefb-6b29-4ac1-b3b7-769040252f3e
growth2 = @reaction_network begin
	@species B(t)=1
	@parameters r=log(2) K=1e3
	r * (1-B/K), B --> 2B
end

# ╔═╡ b16be4b8-8ae7-4edf-9912-3c828d85d0c7
growth3 = @reaction_network begin
	@species B(t)=1
	@parameters r=log(2) K=1e3
	r, B --> 2B
	r/K, B + B --> 0  
end

# ╔═╡ fa828b69-5b17-47cf-8461-7cc025d64501
plot(solve(ODEProblem(growth3, [], (0., 20.))), lw=2, title="Logistic growth")

# ╔═╡ d71fccaa-e8ed-4943-9659-d51d34a10bba
lotka_volterra = @reaction_network begin
	α, x --> 2x  # reproduction prey
	β, y --> 0  # mortality pred
	γ, x + y --> 1.1y  # predation
end

# ╔═╡ dfb27358-b5ab-443f-8f61-ff4bafc25b02
lv_sys = convert(ODESystem, lotka_volterra, combinatoric_ratelaws=false)

# ╔═╡ cd632750-68e4-48de-91be-1fc6c06ffe02
md"## Process dynamics"

# ╔═╡ 82d65028-e90a-4e91-bf67-d2f8da74134b


# ╔═╡ 8001cf7b-ffce-48ba-923e-1e630106bf4b
md"### Saturated processes"

# ╔═╡ 7d119d0f-7e10-474e-b802-e7f91eb01ec4
michaelis_menten_kinetics = @reaction_network begin
	(k₁, k₋₁), E + S <--> ES
	k₂, ES --> E + P
end

# ╔═╡ e1e6362e-bf87-46a8-9191-97e5a121a454
@bind Ks Slider(0:0.1:5, default=1, show_value=true)

# ╔═╡ 8578c4f3-d72f-402e-b886-d34b14b95705
@bind vmax Slider(0:0.1:5, default=1, show_value=true)

# ╔═╡ b8c26bea-b775-4ae2-a65d-edefa2f11c4f
convert(ODESystem, michaelis_menten_kinetics)

# ╔═╡ 52e846e3-b507-4c84-bedb-ce55207c6e2a
latexify(michaelis_menten_kinetics, form=:ode) |> clipboard

# ╔═╡ 5086bbf0-005c-47f8-81bb-029f5cf40614
@variables k₁ k₋₁ ES Eₜ S k₂ v_m K_s

# ╔═╡ 4fe3971f-0786-48f7-adbb-147c37b7dc24
eqmm  = (k₁ + k₋₁) * ES ~ k₂ * (Eₜ - ES) * S

# ╔═╡ b38ed246-abf9-4808-b9eb-0f003224d98d
Symbolics.solve_for([eqmm], [ES])

# ╔═╡ 03a01f2e-cf36-426f-9df2-b8d848be7526
DS = Differential(S)

# ╔═╡ 35e2a2fe-c74a-4f61-bdac-46865f595a56
v = v_m * S / (K_s+S)

# ╔═╡ b37068e8-1426-40d1-9028-347839605450
substitute(expand_derivatives(DS(v)), S=>0)

# ╔═╡ f15f67bf-5c82-4503-b292-151549271917
michealis_menten_direct = @reaction_network begin
	mm(S, vmax, Ks), S --> P
end

# ╔═╡ 77830fad-7719-43c0-aeb2-0037b9017c4e
convert(ODESystem, michealis_menten_direct)

# ╔═╡ ba662f0d-7a45-4e90-b067-a57da5069b2b
md"### Hill function"

# ╔═╡ fe75d100-73b1-4f80-a2cb-bea922521365
@bind n Slider(1:10, default=5)

# ╔═╡ b3217136-69b5-42f2-a4b5-f16371fcf827
# time in minutes
repressor = @reaction_network begin
	@species R(t)=0 mRNA(t)=0
	hillr(R, vtranscr, Ki, 4), 0 --> mRNA  # transcription
	d, mRNA --> 0  # degradation or mRNA
	vtransl, mRNA --> mRNA + R  # translation
	r, R --> 0  # degradation of R
end

# ╔═╡ ef74adcd-c1c2-49ba-b424-8d35881db4f7
convert(ODESystem, repressor)

# ╔═╡ 19b752fa-aad6-4ee5-a28c-29f36e8496c9
let
	prob = ODEProblem(repressor, [], [0, 200], [:vtranscr=>1.2, :Ki=>5, :d=>1/100, :r=>1/10, :vtransl=>1.2])
	plot(solve(prob), lw=2)
end

# ╔═╡ 46b39ff2-7087-44ad-86b2-f630ec2195cc
md"### Logistic growth"

# ╔═╡ 4c73138a-fab1-405b-a442-007881dfd094
logistic = @reaction_network begin
	@species P(t)=2
	@parameters r=1 K=100
	r, P --> 2P
	r/K, 2P --> 0
end	

# ╔═╡ 1955523c-7d60-4e7a-845b-8977dc9b01fb
convert(ODESystem, logistic)

# ╔═╡ 5140228d-a479-4232-8303-6f5b0da339d8
md"## Compartmental models"

# ╔═╡ 350eb38a-4345-460f-89ca-63982ceadcb1
md"### Pharmacokinetic model"

# ╔═╡ 43140694-0a7d-4fe6-af17-54662c856073
pharkin = @reaction_network begin
	@species B(t)=0 T(t)=0
	k₁, G --> B
	k₂, B --> T
end

# ╔═╡ ff620b5c-8aa8-4b48-83bd-de533c53f52e
md"### Reactor with dead zone"

# ╔═╡ 17689a8f-6da2-470a-b015-c55098c8a723
reactor2 = @reaction_network begin
	@species Ad(t)=0 B(t)=0
	@parameters V=100 k=0.3 r=0.05 q=1 c=.1
	q*c, 0 --> Am  # amount of A entering the bulk
	q/(V*(1-f)), Am --> 0  # amount of A leaving the bulk
	r, Am --> B  # reaction
	(k/(V*(1-f)), k/(V*f)), Am <--> Ad  # exchange bulk-dead zone
end

# ╔═╡ 2f2398ed-355f-49ff-b5e0-1fa8215b51f6
convert(ODESystem, reactor2)

# ╔═╡ 40aedd69-a6b0-4f5e-a220-197d2e72d21a
md"f: $(@bind f Slider(0.1:0.1:0.9, default=0.2, show_value=true))"

# ╔═╡ af9e366a-53e7-49d4-9354-44a05fb67888
let
	pars = [:f=>f]
	u0 = [:Am=>0]
	prob = ODEProblem(reactor2, u0, (0.0, 100.), pars)
	sol = solve(prob)
	plot(sol, lw=2, ls=:auto, ylab="amount [mol]", title="Reactor with dead zone\nf=$f")
end

# ╔═╡ 2f829ebd-9a8c-41d0-854c-a8407d2b161f
md"### SIR model"

# ╔═╡ 67422772-6ffe-499a-ba86-4cf70fb2fd53
sir = @reaction_network begin
	@species S(t)=50 I(t)=5 R(t)=0
	β, S + I --> 2I
	γ, I --> R
end

# ╔═╡ e7c84c65-7f49-410e-99d8-b9345b4559d7
md"### Leslie matrix model"

# ╔═╡ f4acd112-3341-4df2-b7f3-4738bf8b3bb3
butterfly = @reaction_network begin
	@species E(t)=10 C(t)=0 P(t)=0 B(t)=0
	@parameters f=0.5 s=0.05 m=0.12
	s, E --> C
	s, C --> P
	s, P --> B
	m, (E, C, P, B) --> 0
	f, B --> B + E
end

# ╔═╡ 79b13806-4d55-4524-9a43-09e3dc6603a5
# oefening: draagkracht voor 100 rupsen

# ╔═╡ e0203b11-b75c-4e79-87d0-c74105894aa8
convert(ODESystem, butterfly)

# ╔═╡ 92d35635-cc1c-4485-8a40-9ca16a3dbee3
let
	pars = [:f=>5, :s=>0.05, :m=>0.12]
	prob_sir = ODEProblem(butterfly, [], (0.0, 50.0), pars)
	plot(solve(prob_sir), lw=2)
end

# ╔═╡ 380b6a1d-f5da-4ba9-b7b0-3dbca0e267b8
md"## Appendix"

# ╔═╡ 07afed5f-6306-4620-98dc-ec729750850b
TableOfContents()

# ╔═╡ c7ee808d-3ec0-4130-b6d1-1fb993178f41
plots = Dict()

# ╔═╡ e2b7dd60-327c-45fd-a8d7-683b5d4f1274
plots["tank"] = plot(t->q*A/r - q*A/r*exp(-t*r/A), 0, 50, lw=2, label=L"V(t)", title="Tank filling problem")

# ╔═╡ 34906c24-e13c-422f-b330-bb905c356276
let
	p = plot(coffee_sol, idxs=1, label="V", title="Coffee volume", lw=2)
	vspan!([2, 3], alpha=0.4, color=:pink, label="")
	plots["coffeevol"] = p
end

# ╔═╡ 0cb4e475-676b-4d9c-9a75-ea0a0659c5b0
let
	p = plot(coffee_sol, idxs=2, label="T", title="Coffee temperature", color=:orange, lw=2)
	vspan!([2, 3], alpha=0.4, color=:pink, label="")
	plots["coffeetemp"] = p
end

# ╔═╡ 92e9cb14-2db1-423d-ab97-e99c726844f7
plots["LMA_50_50"] = lma_illustration(50, 50)

# ╔═╡ a4fb682e-0ae4-4d78-9f27-a3d0a4b6031e
plots["LMA_10_10"] = lma_illustration(10, 10)

# ╔═╡ e57db960-151e-4a49-9e0b-ab60b08a09f4
plots["LMA_100_100"] = lma_illustration(100, 100)

# ╔═╡ 7ac71bef-c7f8-4a07-bc1f-7416757e70ac
plots["LMA_50_150"] = lma_illustration(50, 150)

# ╔═╡ c8da6437-9981-4198-bcb3-9cd081181aa9
plots["NO_H2"] = plot(solve(reactsyst_prob, Rosenbrock23()), lw=2, ylabel="concentration [mol/L]", ls=:auto, idxs=[1,2,3,4])

# ╔═╡ ff4a8011-1d9d-4a5a-a89f-7f032f1619e1
plots["exp_growth"] = plot(solve(ODEProblem(growth1, [], (0., 20.))), lw=2, title="Exponential growth")

# ╔═╡ 8687a35e-5efd-45b7-b36e-acbe6d3932a3
plots["log_growth"] = plot(solve(ODEProblem(growth2, [], (0., 20.))), lw=2, title="Logistic growth")

# ╔═╡ d5e019dc-c8b4-4296-a3c0-921f7dc4a388
let
	lvp = ODEProblem(lotka_volterra, [:x=>1.0, :y=>1.0], (0., 50.), 
			[:α=>1.2, :β=>0.5, :γ=>0.5], combinatoric_ratelaws=false)
	plots["LK"] = plot(solve(lvp), lw=2, label=["prey" "predator"], title="Lotka-Volterra model")
end

# ╔═╡ 90e1b464-b025-4b11-8f97-1b6ef20714ed
let
	p = plot(x->mm(x, vmax, Ks), 0, 5, lw=2, xlab=L"[S]", ylab=L"v", title="Michaelis-Menten kinetics", label="reaction rate")
	hline!([vmax], ls=:auto, lw=2, label=L"v_{max}", alpha=0.6)
	plot!([0, Ks, Ks], [.5vmax, .5vmax, 0], ls=:auto, lw=2, label=L"K_s", alpha=0.6)
	plot!([0, 1/3], [0, (1/3)*vmax/Ks], label="slope at S=0", ls=:auto, lw=2, alpha=0.6)
	plots["MM_default"] = p
end

# ╔═╡ 251cc4df-e800-4945-8923-ffade8a109d0
let
	p = plot(x->mm(x, 1, 1), 0, 5, lw=2, xlab=L"[S]", ylab=L"v", title="Michaelis-Menten kinetics\n parameters", label=L"v_{max}=1, K_s=1")
	plot!(x->mm(x, 2, 1), 0, 5, lw=2, ls=:auto, label=L"v_{max}=2, K_s=1")
	plot!(x->mm(x, 1, 2), 0, 5, lw=2, ls=:auto, label=L"v_{max}=1, K_s=2")
	plot!(x->mm(x, 1, 1/2), 0, 5, lw=2, ls=:auto, label=L"v_{max}=1, K_s=1/2")
	plots["MM_pars"] = p
end

# ╔═╡ f1a46e45-bc97-443e-b7a6-9706a24d6ab0
let
	p = plot(x->mmr(x, vmax, Ks), 0, 5, lw=2, xlab=L"[X]", ylab=L"v", title="repressing Michaelis-Menten\nkinetics", label="reaction rate")
	hline!([vmax], ls=:dot, lw=2, label=L"v_{max}", alpha=0.6)
	plot!([0, Ks, Ks], [.5vmax, .5vmax, 0], ls=:dash, lw=2, label=L"K_s", alpha=0.6)
	plots["MMr_default"] = p
end

# ╔═╡ 213e3087-19c5-4226-828b-773fa2dcdbde
let
	p = plot(x->hill(x, vmax, Ks, n), 0, 5, lw=2, xlab=L"[S]", ylab=L"v", title="Hill equation", label=("v_max=$vmax, Ks=$Ks, n=$(n)"))
	#plot!(x->mm(x, 2, 1), 0, 5, lw=2, ls=:auto, label=L"v_{max}=2, K_s=1")
	#plot!(x->mm(x, 1, 2), 0, 5, lw=2, ls=:auto, label=L"v_{max}=1, K_s=2")
	#plot!(x->mm(x, 1, 1/2), 0, 5, lw=2, ls=:auto, label=L"v_{max}=1, K_s=1/2")
	plots["hill"] = p
end

# ╔═╡ 3402d9fc-b636-4edd-aa42-68f87ca9c012
let
	p = plot(lw=2, xlab=L"[S]", ylab=L"v", title="Hill equation\nvmax=$vmax and Ks=$Ks")
	for n in [1, 2, 5, 10]
		plot!(x->hill(x, vmax, Ks, n), 0, 5, ls=:auto, label="n=$n", lw=2)
	end
	plot!(S->vmax*>(S, Ks), 0, 5, label="n=∞", lw=2)
	plot!([0, Ks, Ks], [.5vmax, .5vmax, 0], ls=:dash, lw=2, label=L"K_s", alpha=0.6)

	plots["hill_n"] = p
end

# ╔═╡ 6732f82b-79e0-4e31-9b8d-3188246b692d
let
	p = plot(x->hillr(x, vmax, Ks, n), 0, 5, lw=2, xlab=L"[S]", ylab=L"v", title="Reverse Hill equation", label=("v_max=$vmax, Ks=$Ks, n=$(n)"))
	#plot!(x->mm(x, 2, 1), 0, 5, lw=2, ls=:auto, label=L"v_{max}=2, K_s=1")
	#plot!(x->mm(x, 1, 2), 0, 5, lw=2, ls=:auto, label=L"v_{max}=1, K_s=2")
	#plot!(x->mm(x, 1, 1/2), 0, 5, lw=2, ls=:auto, label=L"v_{max}=1, K_s=1/2")
	plots["hillr"] = p
end

# ╔═╡ 6a5ba438-a5aa-4481-b5cb-472b94dac066
let
	p = plot(lw=2, xlab=L"[S]", ylab=L"v", title="Reverse Hill equation\nvmax=$vmax and Ks=$Ks")
	for n in [1, 2, 5, 10]
		plot!(x->hillr(x, vmax, Ks, n), 0, 5, ls=:auto, label="n=$n", lw=2)
	end
	plot!(S->vmax*<(S, Ks), 0, 5, label="n=∞", lw=2)
	plot!([0, Ks, Ks], [.5vmax, .5vmax, 0], ls=:dash, lw=2, label=L"K_s", alpha=0.6)

	plots["hillr_n"] = p
end

# ╔═╡ 1d0e45e5-4ef8-46a1-9a7b-f76302889716
let
	p = plot(y->(1-y/100) * y, 0, 120, lw=2, xlab=L"P", ylab="growth", title="Logistic growth", label=L"r=1, K=100")
	
	vline!([50], ls=:auto, lw=2, label=L"K/2", alpha=0.6)
	vline!([100], ls=:auto, lw=2, label=L"K", alpha=0.6)

	plots["logistic"] = p
end

# ╔═╡ 30bbcf76-f92e-4291-aec2-1bfdd89bdc38
let
	prob = ODEProblem(logistic, [], (0.0, 10.0))
	plots["logistic_sol"] = plot(solve(prob), lw=2, title="The logistic equation")
end

# ╔═╡ fdcfd565-349e-4047-842e-32d2dd8329e6
let
	for f in [0.1, 0.33, 0.66]
		pars = [:f=>f]
		u0 = [:Am=>0]
		prob = ODEProblem(reactor2, u0, (0.0, 100.), pars)
		sol = solve(prob)
		plots["reactor_dead_zone_$(f)"] = plot(sol, lw=2, ls=:auto, ylab="amount [mol]", title="Reactor with dead zone\nf=$f")
	end
end

# ╔═╡ 0e2b9af8-9701-4d1c-838d-c289669522d1
let
	β, γ = 0.03, 0.3
	prob_sir = ODEProblem(sir, [], (0.0, 30.0), (;β, γ))
	plots["sir_beta=$(β)_gamma=$v"] = plot(solve(prob_sir), lw=2, title="SIR model\nβ=$β γ=$γ", ls=:auto)
end

# ╔═╡ ac6c7513-3476-4968-9d6e-eec6d09e019d
let
	β, γ = 0.3, 0.3
	prob_sir = ODEProblem(sir, [], (0.0, 30.0), (;β, γ))
	plots["sir_beta=$(β)_gamma=$v"] = plot(solve(prob_sir), lw=2, title="SIR model\nβ=$β γ=$γ", ls=:auto)
end

# ╔═╡ d408c7dc-4a57-400d-bad2-7ac02cc4fa14
let
	β, γ = 0.01, 0.3
	prob_sir = ODEProblem(sir, [], (0.0, 30.0), (;β, γ))
	plots["sir_beta=$(β)_gamma=$v"] = plot(solve(prob_sir), lw=2, title="SIR model\nβ=$β γ=$γ", ls=:auto)
end

# ╔═╡ 2f4153cb-5aab-495b-98e8-cbdd8ae99816
plots

# ╔═╡ 4d4d8bea-faf7-4045-a8a7-fc0e6b92d7ea
length(plots)

# ╔═╡ Cell order:
# ╠═dfca2f9f-0134-461c-a18b-f66f2bf02943
# ╠═093b722d-28af-4219-8546-39a3262146b2
# ╠═a52da2c2-f7df-11ee-033f-8500edb3c03f
# ╠═8cff27a7-fde1-4b49-8ad6-513302997a4e
# ╠═34bec0a1-40e8-48a2-9109-94872aaff1b9
# ╠═abebedae-b977-43ae-aaa0-6b00990a5de4
# ╠═0d09ba2e-3cac-4051-b98f-26b79736b225
# ╠═fdc534e1-e334-48c1-aac4-5a89c47484e0
# ╠═104ad25e-743b-4f47-8b10-0d3d6715f95c
# ╠═e2b7dd60-327c-45fd-a8d7-683b5d4f1274
# ╠═73eb7d0a-5433-4e3d-a008-748db66b8ef9
# ╠═d3b84441-ed9f-436d-a690-660c5f4b8fbd
# ╠═ddd43577-eb2e-4c72-b829-d7195c165ddf
# ╠═10a26b97-8b0a-454e-afbf-141aef4aa04f
# ╠═206dbd4b-2ec0-4591-b7f9-d8e78d568c2f
# ╠═bf062835-538d-437a-bae6-6309c66ebd19
# ╠═34906c24-e13c-422f-b330-bb905c356276
# ╠═f3737612-5458-4c6e-a634-246b2cb8cb05
# ╠═0cb4e475-676b-4d9c-9a75-ea0a0659c5b0
# ╠═071b8f85-75cc-487b-a13f-64428bab7147
# ╠═8ae850c8-ec2c-4a6c-9c37-1c1f93bb56e9
# ╠═92e9cb14-2db1-423d-ab97-e99c726844f7
# ╠═a4fb682e-0ae4-4d78-9f27-a3d0a4b6031e
# ╠═e57db960-151e-4a49-9e0b-ab60b08a09f4
# ╠═7ac71bef-c7f8-4a07-bc1f-7416757e70ac
# ╠═c01904be-7ee0-4b43-bd29-6fac1e38b3c6
# ╠═689f9c80-5680-4cc5-a60d-bacd846f0925
# ╠═ad511e07-230f-4405-9613-40f2e55676e6
# ╠═730255e3-ba36-482c-bff3-8435d348fe40
# ╠═959b8b64-959a-4869-b499-159eccfa0768
# ╠═8a752cbe-8220-41f8-997f-f2eb321b2a2b
# ╠═063f7b0d-c8e6-4128-87f9-d3b5b5139670
# ╠═691b4a76-fc29-4cb9-a3ca-26d03bd84cec
# ╠═7c4e736b-fac7-47c1-b73c-ef7d46151e15
# ╠═2b80012d-78c8-4cd3-9b8b-6d0177afa963
# ╠═c8da6437-9981-4198-bcb3-9cd081181aa9
# ╠═7ad626ab-12a5-4a53-9991-7512de5092fb
# ╠═ca2304c8-9a09-4063-9478-2207925d6444
# ╠═bbfd0201-c25e-446e-828e-6632a7e0c116
# ╠═c77660bc-a478-45d4-8d8b-c5147d70b123
# ╠═a5177857-d6e4-4004-92e1-34bbfb42de53
# ╠═ff4a8011-1d9d-4a5a-a89f-7f032f1619e1
# ╠═e6b6eefb-6b29-4ac1-b3b7-769040252f3e
# ╠═8687a35e-5efd-45b7-b36e-acbe6d3932a3
# ╠═b16be4b8-8ae7-4edf-9912-3c828d85d0c7
# ╠═fa828b69-5b17-47cf-8461-7cc025d64501
# ╠═d71fccaa-e8ed-4943-9659-d51d34a10bba
# ╠═dfb27358-b5ab-443f-8f61-ff4bafc25b02
# ╟─d5e019dc-c8b4-4296-a3c0-921f7dc4a388
# ╠═cd632750-68e4-48de-91be-1fc6c06ffe02
# ╠═82d65028-e90a-4e91-bf67-d2f8da74134b
# ╠═8001cf7b-ffce-48ba-923e-1e630106bf4b
# ╠═7d119d0f-7e10-474e-b802-e7f91eb01ec4
# ╠═e1e6362e-bf87-46a8-9191-97e5a121a454
# ╠═8578c4f3-d72f-402e-b886-d34b14b95705
# ╟─90e1b464-b025-4b11-8f97-1b6ef20714ed
# ╟─251cc4df-e800-4945-8923-ffade8a109d0
# ╠═b8c26bea-b775-4ae2-a65d-edefa2f11c4f
# ╠═52e846e3-b507-4c84-bedb-ce55207c6e2a
# ╠═0686fc66-5428-451f-aa72-c0250ad4bf67
# ╠═5086bbf0-005c-47f8-81bb-029f5cf40614
# ╠═4fe3971f-0786-48f7-adbb-147c37b7dc24
# ╠═b38ed246-abf9-4808-b9eb-0f003224d98d
# ╠═03a01f2e-cf36-426f-9df2-b8d848be7526
# ╠═35e2a2fe-c74a-4f61-bdac-46865f595a56
# ╠═b37068e8-1426-40d1-9028-347839605450
# ╟─f1a46e45-bc97-443e-b7a6-9706a24d6ab0
# ╠═f15f67bf-5c82-4503-b292-151549271917
# ╠═77830fad-7719-43c0-aeb2-0037b9017c4e
# ╠═ba662f0d-7a45-4e90-b067-a57da5069b2b
# ╠═fe75d100-73b1-4f80-a2cb-bea922521365
# ╠═213e3087-19c5-4226-828b-773fa2dcdbde
# ╠═3402d9fc-b636-4edd-aa42-68f87ca9c012
# ╠═6732f82b-79e0-4e31-9b8d-3188246b692d
# ╠═6a5ba438-a5aa-4481-b5cb-472b94dac066
# ╠═b3217136-69b5-42f2-a4b5-f16371fcf827
# ╠═ef74adcd-c1c2-49ba-b424-8d35881db4f7
# ╠═19b752fa-aad6-4ee5-a28c-29f36e8496c9
# ╠═46b39ff2-7087-44ad-86b2-f630ec2195cc
# ╠═1d0e45e5-4ef8-46a1-9a7b-f76302889716
# ╠═4c73138a-fab1-405b-a442-007881dfd094
# ╠═1955523c-7d60-4e7a-845b-8977dc9b01fb
# ╠═30bbcf76-f92e-4291-aec2-1bfdd89bdc38
# ╠═5140228d-a479-4232-8303-6f5b0da339d8
# ╠═350eb38a-4345-460f-89ca-63982ceadcb1
# ╠═43140694-0a7d-4fe6-af17-54662c856073
# ╠═ff620b5c-8aa8-4b48-83bd-de533c53f52e
# ╠═17689a8f-6da2-470a-b015-c55098c8a723
# ╠═2f2398ed-355f-49ff-b5e0-1fa8215b51f6
# ╟─40aedd69-a6b0-4f5e-a220-197d2e72d21a
# ╠═af9e366a-53e7-49d4-9354-44a05fb67888
# ╠═fdcfd565-349e-4047-842e-32d2dd8329e6
# ╠═2f829ebd-9a8c-41d0-854c-a8407d2b161f
# ╠═67422772-6ffe-499a-ba86-4cf70fb2fd53
# ╠═0e2b9af8-9701-4d1c-838d-c289669522d1
# ╠═ac6c7513-3476-4968-9d6e-eec6d09e019d
# ╠═d408c7dc-4a57-400d-bad2-7ac02cc4fa14
# ╠═e7c84c65-7f49-410e-99d8-b9345b4559d7
# ╠═f4acd112-3341-4df2-b7f3-4738bf8b3bb3
# ╠═79b13806-4d55-4524-9a43-09e3dc6603a5
# ╠═e0203b11-b75c-4e79-87d0-c74105894aa8
# ╠═92d35635-cc1c-4485-8a40-9ca16a3dbee3
# ╟─380b6a1d-f5da-4ba9-b7b0-3dbca0e267b8
# ╠═07afed5f-6306-4620-98dc-ec729750850b
# ╠═c7ee808d-3ec0-4130-b6d1-1fb993178f41
# ╠═2f4153cb-5aab-495b-98e8-cbdd8ae99816
# ╠═4d4d8bea-faf7-4045-a8a7-fc0e6b92d7ea
