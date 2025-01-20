### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 7bc363b0-9415-4954-807f-81a308bde531
begin
	import Pkg
	Pkg.activate("..")
end

# ╔═╡ 52b28a4c-b0bb-11ef-2841-17ecfb596676
using Plots, PlutoUI, DifferentialEquations, ForwardDiff, Catalyst

# ╔═╡ 77860f25-1b30-4cf1-812a-03ff2100bacd
using Turing

# ╔═╡ c897d6e6-0b33-4dfb-89f3-f7e05bb9f73b
using StatsPlots, Optim

# ╔═╡ e96779da-57c0-4c11-b3ab-40dc5ae81bee
using StatsBase

# ╔═╡ bfcc4b4e-073e-401e-851c-d01ef828028a
md"""
# The minimal glucose model and dynamic compensation

The Minimal Model of Glucose Regulation is a mathematical model used to describe how the body regulates glucose (sugar) levels in the blood. It was developed by Richard Bergman and Claudio Cobelli in the late 1970s and has become a cornerstone in diabetes research.

We will use this exercise to study insulin sensitivity.

The basic model considers only the concentration of glucose $G(t)$ in nM and the concentration of insulin $I(t)$ in nM:

- glucose is added to a system with a zeroth-order rate of $m$ (later $m(t)$ if we model a non-fixed input);
- glucose is removed from the blood with a rate of $sGI$, where $s$ is the insulin sensitivity;
- insulin decays according to first-order kinetics with a rate parameter $\gamma$
- $\beta$-cells produce insulin as a response to higher glucose concentrations. This is according to a saturated process, so it is well approximated using a Hill function ($n=2$). The  rate of insulin production is given by $qBf(G)$, with $B$ the amount of $\beta$-cells and $q$ the maximal rate of insulin production/unit of cells.

The following plot gives a fairly realistic response of insulin production as a function of glucose concentration in the blood. 
"""

# ╔═╡ e93bc1a5-e0b1-425c-b069-d5ce80756e60
f_insulin(G) = hill(G, 1, 5, 2)

# ╔═╡ 2ac43184-a393-4d39-9c3b-c299c3371b09
plot(f_insulin, 0, 30, xlab="G nm", ylab="f(G)", title="Insulin production rate")

# ╔═╡ 75a5b5d9-f8d7-4728-83dc-03dde502dcfb
md"Below is a reaction network implementing this model. All parameters are set to 1 for didactic purposes." 

# ╔═╡ f62c2b4a-a9e8-472b-b07b-44072da08ee8
glucose_insuline_circuit = @reaction_network begin
	@parameters q=1 s=1 γ=1 m=1 B=1 Ks=1.0
	m, 0 --> G
	s * I, G --> 0
	B * hill(G, q, Ks, 2), 0 --> I
	γ, I --> 0
end

# ╔═╡ 99ccbd1b-b9e4-4afc-8d9a-8c1055bcce0f
convert(ODESystem, glucose_insuline_circuit)

# ╔═╡ c319ae6e-2f02-4a7e-a0b0-4a3f57fdcb66
md"You can set initial values of glucose and insulin here."

# ╔═╡ 7623a596-d6d4-4e0a-b29c-f2f7de0d2940
G0 = 5.0

# ╔═╡ 0621968c-b39a-4346-a86a-42fa2cde724c
I0 = 0.0

# ╔═╡ 3c1465d2-0886-4b74-baee-55c7dad41f9a
md"Simulate the system over a time interval of 0 to 10 hours with $m=1$. Plot the results. What are the steady state concentrations for the two species? Does this depend on initial glucose levels (given enough time)?"

# ╔═╡ bbbc3335-4190-4bc7-91ba-e96476a75751
@unpack t, G, I = glucose_insuline_circuit

# ╔═╡ afc40964-0a97-43a6-ab27-5fe96ecd319b
oprob1 = ODEProblem(glucose_insuline_circuit, [G=>G0, I=>I0], (0., 10.), [:m=>1.0])

# ╔═╡ f050b32e-1a14-4d5e-9fb3-5dbc4dfbf79c
sol1 = solve(oprob1, Tsit5());

# ╔═╡ 7e35b56b-44a0-46df-a812-f5fb437b672e
plot(sol1)

# ╔═╡ 00db9d76-c4dd-450d-b7e3-cf87f8bfeffb
plot(sol1, idxs=(G, I))

# ╔═╡ aab02fef-6d84-4f53-aff5-6f94bf2bd460
md"Now simulate the system but rather than with $m$ being a constant glucose input, we give in a pulse of glucose (i.e., drinking a soda) with a peak at $t=5$h. Note that our parameter now depends on the time!"

# ╔═╡ 8a5ce760-55bb-4610-9226-32c3d8376748
glucose_pulse(t) = .5 + exp(-(t-5)^2)

# ╔═╡ a2fb8c21-c594-44f3-92df-99e6d3afeba1
plot(glucose_pulse, 0, 10, label="G [nm]", xlab="t")

# ╔═╡ ba68f740-ba8b-4efa-be2d-ac50c77dfa1a
glucose_pulse(t)

# ╔═╡ 31aaced7-1fe4-4117-80f3-0fa70035822a
md"Up to now, we set $s$, the insulin sensitivity to 1. This parameter represents how sensitive the body is to insuline in taking up glucose. Aging and obesity increase glucose resistance ($1/s$), resulting in diabetes! Explore the effect of this parameter on your plots below."

# ╔═╡ 69e21e9b-7d0e-44e3-8b33-ee7102dadadf
@bind s Slider(0.1:0.1:5, default=1, show_value=true)

# ╔═╡ 6890e3ec-0875-423e-89c1-988411c869e3
oprob2 = ODEProblem(glucose_insuline_circuit, [G=>0., I=>0.0], (0., 10.), [:m=>glucose_pulse(t), :s=>s])

# ╔═╡ fb439f9d-ee88-43cd-b91b-44ed8013b27d
sol2 = solve(oprob2, Tsit5());

# ╔═╡ 095e359a-2ec4-4637-8477-93510bef23ed
plot(sol2)

# ╔═╡ b4073e81-6db8-4be6-a2d1-afa01640dbcf
s

# ╔═╡ 4d58a832-9157-42fe-b578-c3fc29f01ea8
md"We have made a function here that computes the steady state glucose concentration after 100 hours. Plot this as a function of s."

# ╔═╡ 7cf11e6f-7ef4-478f-8cc7-eae93374cc58
function glucose_steady_state(s)
	oprob = ODEProblem(glucose_insuline_circuit, [G=>0., I=>0.0], (0., 100.), [:m=>1/2, :s=>s])
	sol = solve(oprob, Tsit5())
	return sol[G][end]
end

# ╔═╡ 5eda8c4a-2baa-4987-b408-dbe7383abe95
glucose_steady_state(s)

# ╔═╡ 34788b1d-1530-4b88-8690-03fe9af52a8b
plot(glucose_steady_state, 0.1, 5, xlab="s", ylab="Gss")

# ╔═╡ f616da2f-b51e-4328-bd3c-47e6d7377cbb
md"Using automatic differentiation (ForwardDiff), compute the absolute and relative sensitity indices. Is this system sensitive to the insuline sensitivity?"

# ╔═╡ 97636a73-4732-4377-b760-c36cef13904b
fsens(s) = ForwardDiff.derivative(glucose_steady_state, s)

# ╔═╡ a67d70ee-c8b5-4f22-926c-f48c5dcf9815
fsens(s)

# ╔═╡ 1bb5834a-6a4e-40ec-82c1-f2c505ecfeb8
frelsens(s) = ForwardDiff.derivative(glucose_steady_state, s) * glucose_steady_state(s) / s

# ╔═╡ a3156073-07bd-4f13-be5e-d19483fbad28
frelsens(s)

# ╔═╡ 0ae57436-0b8c-41ce-9d9d-8c48ed59b820
md"""
We see that the final glucose concentration is highly dependent on $s$! This seems to be a flaw in the model, as we can imagine that the physiological parameters can greatly differ from person to person (for example, a person can have a large pancreas). 

A mechanism that stabilizes this is called *dynamic compensation*. Simply put, we have assumed here that the amount of beta cells ($B$) is fixed. However, in practice, these cells are capable of dividing, growing, and thus producing more insulin. Their growth rate depends on the concentration of glucose, creating an additional feedback loop that stabilizes the physiological circuit. The growth rate of the $\beta$-cells follows a sigmoid shape, being negative when $G$ is smaller than a threshold and positive if $G$ exceeds this threshold.

```julia
μ(G), B --> 2B
```

This curve is plotted below.
"""

# ╔═╡ 87d9c826-56af-409a-a884-d784cfa16a64
μ(G) = 0.3atan(0.5(G-1))

# ╔═╡ 4d95e77b-fab0-454b-8319-1fe83e12f9e9
plot(μ, 0, 30, xlab="G", label="μ(G)", title="Glucose-dependend growth rate")

# ╔═╡ ea46f5f9-47a8-4641-a7e3-08a1a29cd04f
md"Add dynamic compensation to the model and show that this greatly reduces the sentitivty w.r.t. $s$."

# ╔═╡ 870e9229-374a-40b8-97f7-1a84d1a857e9
md"""
Robust version with dynamic compensation
```julia
glucose_insuline_circuit = @reaction_network begin
	@parameters q=1 s=1 γ=1 m=1 Ks=1.0
	@species B(t)=1
	m, 0 --> G
	s * I, G --> 0
	q * B * hill(G, 1, Ks, 2), 0 --> I
	γ, I --> 0
	μ(G), B --> 2B
end
```
"""

# ╔═╡ e37ce327-504b-4fe4-baea-493cda834b4a
md" this does not work yet. Split and use different model: https://allendowney.github.io/ModSimPy/chap18.html"

# ╔═╡ e66378f7-8ec8-4f7f-ab40-5c64479f865e
tmeas = 0:2:8

# ╔═╡ bcafb5c2-4f3c-47f4-b786-8f656969abe6
length(tmeas)

# ╔═╡ 08d3cf6e-dc0b-4a56-a241-0d0033e79fad
Gmeas = [92, 350, 287, 251, 240]

# ╔═╡ c0d39fc8-d779-465e-91c1-3a6f9b694ee3
Imeas = [11, 26, 130, 85, 51]

# ╔═╡ c2c18077-add0-4ff9-a122-f248dbd1da79
begin
	scatter(tmeas, Gmeas, xlab="t [h]", label="G [mg/dL]")
	scatter!(tmeas, Imeas, label="I [μU/mL]")
end

# ╔═╡ e93a2495-aa5b-46f6-ae3c-b94559583c60
@model function glucose_inference(tmeas, Gmeas, Imeas)
	@assert length(tmeas) == length(Gmeas) == length(Imeas)
	N = length(tmeas)
	σGsq ~ Gamma(5)
	σIsq ~ Gamma(3)
	q~10LogNormal()
	s~10LogNormal()
	γ~10LogNormal()
	B=1
 	m~10LogNormal()
	Ks~10LogNormal()
	G0~100LogNormal()
	I0~1LogNormal()
	prob = remake(oprob1; u0=[G=>G0, I=>I0], tspan=(0., 8.), 
									p=[:m=>m, :q=>q, :s=>s,
										:γ=>γ, :Ks=>Ks, :B=>B])
	sol = solve(prob, saveat=tmeas)
	for i in 1:N
		Gmeas[i] ~ Normal(sol[G][i], sqrt(σGsq))
		Imeas[i] ~ Normal(sol[I][i], sqrt(σIsq))
	end
	return sol
end

# ╔═╡ 55e4db93-f376-4cbb-829d-fea9ce3f3998
glucose_inference(tmeas, Gmeas, Imeas)()

# ╔═╡ 8989e653-e3e8-4e4f-8fc2-596e931ed21c
mod_data = glucose_inference(tmeas, Gmeas, Imeas)

# ╔═╡ b763bb0c-581e-46ec-8104-53df25369642
mod_data() |> plot

# ╔═╡ 5e187f36-3883-4f08-b08c-a55d24e409cd
chain1 = sample(mod_data, NUTS(), 100)

# ╔═╡ 2bf57ee8-ac1d-412f-a4c8-6e4ce87a74c5
summarize(chain1)

# ╔═╡ aa32db53-c944-442a-88df-5937f380b54c
#=╠═╡
plot(chain1)
  ╠═╡ =#

# ╔═╡ 7b247ff1-f873-4b00-90c3-4c7a87485fb2
#=╠═╡
scatter(chain1[:q], chain1[:B])
  ╠═╡ =#

# ╔═╡ e591eff0-fd70-43d3-bf2a-34309a8ecd62
#=╠═╡
begin
	p = scatter(tmeas, Gmeas, xlab="t [h]", label="G [mg/dL]")
	scatter!(tmeas, Imeas, label="I [μU/mL]")
	N = length(sample_sols)
	for i in 1:10
		n = rand(1:N)
		plot!(sample_sols[n], lw=0.5, alpha=0.8, label="", color=[:blue :orange])
	end
	p
end
  ╠═╡ =#

# ╔═╡ 04d6078e-f46b-43f2-988c-d91252b26747
#=╠═╡
sample_sols = generated_quantities(mod_data, chain1);
  ╠═╡ =#

# ╔═╡ 851e6a2f-3569-489d-b375-149329c68a7c
mapGI = optimize(mod_data, MAP(), NelderMead())

# ╔═╡ 246a06ba-ef48-47ba-aef0-51f07ceb96e4


# ╔═╡ 403f6e68-e827-4885-830b-cb0d7f872ee0
coeftable(mapGI)

# ╔═╡ 2f992357-4608-40e6-9954-8306f31e9768


# ╔═╡ 6ecda2c4-2a1d-4d24-b7b8-ec956a3e8d70
#coeftable(mapG)

# ╔═╡ Cell order:
# ╠═52b28a4c-b0bb-11ef-2841-17ecfb596676
# ╠═7bc363b0-9415-4954-807f-81a308bde531
# ╠═bfcc4b4e-073e-401e-851c-d01ef828028a
# ╠═e93bc1a5-e0b1-425c-b069-d5ce80756e60
# ╠═2ac43184-a393-4d39-9c3b-c299c3371b09
# ╟─75a5b5d9-f8d7-4728-83dc-03dde502dcfb
# ╠═f62c2b4a-a9e8-472b-b07b-44072da08ee8
# ╠═99ccbd1b-b9e4-4afc-8d9a-8c1055bcce0f
# ╠═c319ae6e-2f02-4a7e-a0b0-4a3f57fdcb66
# ╠═7623a596-d6d4-4e0a-b29c-f2f7de0d2940
# ╠═0621968c-b39a-4346-a86a-42fa2cde724c
# ╠═3c1465d2-0886-4b74-baee-55c7dad41f9a
# ╠═afc40964-0a97-43a6-ab27-5fe96ecd319b
# ╠═bbbc3335-4190-4bc7-91ba-e96476a75751
# ╠═f050b32e-1a14-4d5e-9fb3-5dbc4dfbf79c
# ╠═7e35b56b-44a0-46df-a812-f5fb437b672e
# ╠═00db9d76-c4dd-450d-b7e3-cf87f8bfeffb
# ╠═aab02fef-6d84-4f53-aff5-6f94bf2bd460
# ╠═8a5ce760-55bb-4610-9226-32c3d8376748
# ╟─a2fb8c21-c594-44f3-92df-99e6d3afeba1
# ╠═ba68f740-ba8b-4efa-be2d-ac50c77dfa1a
# ╠═6890e3ec-0875-423e-89c1-988411c869e3
# ╠═fb439f9d-ee88-43cd-b91b-44ed8013b27d
# ╠═095e359a-2ec4-4637-8477-93510bef23ed
# ╠═31aaced7-1fe4-4117-80f3-0fa70035822a
# ╠═69e21e9b-7d0e-44e3-8b33-ee7102dadadf
# ╠═b4073e81-6db8-4be6-a2d1-afa01640dbcf
# ╠═4d58a832-9157-42fe-b578-c3fc29f01ea8
# ╠═7cf11e6f-7ef4-478f-8cc7-eae93374cc58
# ╠═5eda8c4a-2baa-4987-b408-dbe7383abe95
# ╠═34788b1d-1530-4b88-8690-03fe9af52a8b
# ╠═f616da2f-b51e-4328-bd3c-47e6d7377cbb
# ╠═97636a73-4732-4377-b760-c36cef13904b
# ╠═a67d70ee-c8b5-4f22-926c-f48c5dcf9815
# ╠═1bb5834a-6a4e-40ec-82c1-f2c505ecfeb8
# ╠═a3156073-07bd-4f13-be5e-d19483fbad28
# ╠═0ae57436-0b8c-41ce-9d9d-8c48ed59b820
# ╠═87d9c826-56af-409a-a884-d784cfa16a64
# ╠═4d95e77b-fab0-454b-8319-1fe83e12f9e9
# ╠═ea46f5f9-47a8-4641-a7e3-08a1a29cd04f
# ╟─870e9229-374a-40b8-97f7-1a84d1a857e9
# ╠═e37ce327-504b-4fe4-baea-493cda834b4a
# ╠═e66378f7-8ec8-4f7f-ab40-5c64479f865e
# ╠═bcafb5c2-4f3c-47f4-b786-8f656969abe6
# ╠═08d3cf6e-dc0b-4a56-a241-0d0033e79fad
# ╠═c0d39fc8-d779-465e-91c1-3a6f9b694ee3
# ╠═c2c18077-add0-4ff9-a122-f248dbd1da79
# ╠═77860f25-1b30-4cf1-812a-03ff2100bacd
# ╠═e93a2495-aa5b-46f6-ae3c-b94559583c60
# ╠═55e4db93-f376-4cbb-829d-fea9ce3f3998
# ╠═8989e653-e3e8-4e4f-8fc2-596e931ed21c
# ╠═b763bb0c-581e-46ec-8104-53df25369642
# ╠═5e187f36-3883-4f08-b08c-a55d24e409cd
# ╠═2bf57ee8-ac1d-412f-a4c8-6e4ce87a74c5
# ╠═c897d6e6-0b33-4dfb-89f3-f7e05bb9f73b
# ╠═aa32db53-c944-442a-88df-5937f380b54c
# ╠═7b247ff1-f873-4b00-90c3-4c7a87485fb2
# ╠═e591eff0-fd70-43d3-bf2a-34309a8ecd62
# ╠═04d6078e-f46b-43f2-988c-d91252b26747
# ╠═851e6a2f-3569-489d-b375-149329c68a7c
# ╠═246a06ba-ef48-47ba-aef0-51f07ceb96e4
# ╠═e96779da-57c0-4c11-b3ab-40dc5ae81bee
# ╠═403f6e68-e827-4885-830b-cb0d7f872ee0
# ╠═2f992357-4608-40e6-9954-8306f31e9768
# ╠═6ecda2c4-2a1d-4d24-b7b8-ec956a3e8d70
