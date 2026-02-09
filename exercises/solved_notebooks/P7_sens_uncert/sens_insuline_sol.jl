### A Pluto.jl notebook ###
# v0.19.46

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

# ╔═╡ 7bc363b0-9415-4954-807f-81a308bde531
begin
	import Pkg
	Pkg.activate(".")
end

# ╔═╡ 52b28a4c-b0bb-11ef-2841-17ecfb596676
using Plots, PlutoUI, DifferentialEquations, ForwardDiff, Catalyst

# ╔═╡ bfcc4b4e-073e-401e-851c-d01ef828028a
md"""
### Exercise: The minimal glucose model and dynamic compensation - sensitivity

The Minimal Model of Glucose Regulation is a mathematical model used to describe how the body regulates glucose (sugar) levels in the blood. It was developed by Richard Bergman and Claudio Cobelli in the late 1970s and has become a cornerstone in diabetes research.

We will use this exercise to study insulin sensitivity.

The basic model considers only the concentration of glucose $G(t)$ in $mmol/L$ and the concentration of insulin $I(t)$ in $mmol/L$:

- Glucose is added to a system with a zeroth-order rate of $m$ (later $m(t)$ if we model a non-fixed input).
- Glucose is removed from the blood with a rate of $sGI$, where $s$ is the insulin sensitivity.
- Insulin decays according to first-order kinetics with a rate parameter $\gamma$
-  $\beta$-cells produce insulin as a response to higher glucose concentrations. This is according to a saturated process, so it is well approximated using a Hill function ($n=2$). The  rate of insulin production is given by $qBf(G)$, with $B$ the amount of $\beta$-cells, $q$ the maximal rate of insulin production/unit of cells and $f(G)$ the Hill function.
"""

# ╔═╡ a5789d9d-a690-4c6b-8b8d-fbfe63f659d0
md"""
The Hill function is defined as:

$$hill(X, v, K, n) = \cfrac{v\,X^n}{X^n + K^n}$$

In order to have an idea of how it looks like, lets define it as `f_insuline(G)` for $v=1$, $K=5$ and $n=2$:
"""

# ╔═╡ e93bc1a5-e0b1-425c-b069-d5ce80756e60
f_insulin(G) = hill(G, 1, 5, 2)

# ╔═╡ be33a8a3-a3f5-448d-9390-16f7c4ca099d
md"""
The following plot gives a fairly realistic response of insulin production as a function of glucose concentration in the blood.
"""

# ╔═╡ 2ac43184-a393-4d39-9c3b-c299c3371b09
plot(f_insulin, 0, 30, xlab="G (mmol/L)", ylab="f(G)", title="Insulin production rate")

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

# ╔═╡ 882e106f-d3a5-4ccb-ae65-cfb5e4e5de57
md"""
Convert the system to a symbolic differential equation model and inspect the differential equations. You do not need to make a new variable, just call `convert` with the right arguments.
"""

# ╔═╡ 99ccbd1b-b9e4-4afc-8d9a-8c1055bcce0f
# missing            # Uncomment and complete the instruction
convert(ODESystem, glucose_insuline_circuit)

# ╔═╡ 3c1465d2-0886-4b74-baee-55c7dad41f9a
md"""
Simulate the system over a time interval of $0$ to $10$ hours with $m=1$ for various initial glucose concentrations (e.g., between $0.1$ and $5$) by means of the variable `G0` bound to the slider just here below. Use an initial insuline concentration of $0.0$.
"""

# ╔═╡ 5bac69d1-0bd4-41bd-b15d-5ca6292a31e4
@bind G0 Slider(0.1:0.1:5, default=1, show_value=true)

# ╔═╡ 6ccd5ecc-5579-4ce5-9425-d1bbae0bc5ba
G0

# ╔═╡ afc40964-0a97-43a6-ab27-5fe96ecd319b
# Putting a semi-colon (;) after an instruction will hide its return value.
# oprob1 = missing;     # Uncomment and complete the instruction
oprob1 = ODEProblem(glucose_insuline_circuit, [:G=>G0, :I=>0.0], (0., 10.), []);

# ╔═╡ f050b32e-1a14-4d5e-9fb3-5dbc4dfbf79c
# sol1 = missing;       # Uncomment and complete the instruction
sol1 = solve(oprob1, Tsit5(), saveat=0.01);

# ╔═╡ cfbbfe48-a2c3-4691-9fac-b0a18f99222c
md"""
Plot the results. Use thereby `ylim=(0.0, 5.5)`.
"""

# ╔═╡ 7e35b56b-44a0-46df-a812-f5fb437b672e
# missing              # Uncomment and complete the instruction
plot(sol1, ylim=(0.0, 5.5))

# ╔═╡ 0133ab41-2c8c-4147-9594-1ad238833a9e
md"""
What are the steady state concentrations for the two species? Does this depend on initial glucose levels (given enough time)?

- Answer: missing
"""

# ╔═╡ 479e0f11-e18e-4523-8d24-b2026acdeaef
md"""
Check out the final glucose and insuline concentrations (at the end time).
"""

# ╔═╡ 824677b4-6533-453f-909f-bf04fe4bf5d9
# missing               # Uncomment and complete the instruction
(sol1[:G][end], sol1[:I][end])

# ╔═╡ cef7f30d-43ef-4f2c-8bf3-35c140d237ce
md"""
Create a vector named `u1_guess` with the previous final values.
"""

# ╔═╡ 9d3ba8bd-e9c3-4bb6-957d-c7f8211a2cda
# u1_guess = missing    # Uncomment and complete the instruction
u1_guess = [:G=>sol1[:G][end], :I=>sol1[:I][end]]

# ╔═╡ 3f5e2195-da93-40fb-88ec-9e85a7bb6c24
md"""
Calculate the steady-state values of glucose and insuline.
"""

# ╔═╡ 3dae761f-5fb6-4f1c-a440-d2c21a6edb17
# Gw, Iw = missing      # Uncomment and complete the instruction
Gw, Iw = solve(SteadyStateProblem(ODEProblem(glucose_insuline_circuit, u1_guess, (0., 10.), [])))

# ╔═╡ 314aa77c-5906-4772-94c5-f0b7cf648de0
md"""
Check ou the steady state values for glucose and insulin.
"""

# ╔═╡ 3b56f88c-0ca8-4f33-8516-11e639cb6a1a
# missing               # Uncomment and complete the instruction
(Gw, Iw)

# ╔═╡ 00db9d76-c4dd-450d-b7e3-cf87f8bfeffb
# ╠═╡ disabled = true
#=╠═╡
plot(sol1, idxs=[:G, :I])
  ╠═╡ =#

# ╔═╡ 5d1134d5-734e-4851-a91c-615dc37e3c68
# ╠═╡ disabled = true
#=╠═╡
plot(sol1, idxs=(:G, :I), xlabel="G", ylabel="I", title="I versus G")
  ╠═╡ =#

# ╔═╡ aab02fef-6d84-4f53-aff5-6f94bf2bd460
md"""
Now simulate the system but rather than with $m$ being a constant glucose input, we give in a pulse of glucose (i.e., drinking a soda) with a peak at $t=5$h. **Note that our parameter now depends on the time!**
"""

# ╔═╡ 8a5ce760-55bb-4610-9226-32c3d8376748
glucose_pulse(t) = .5 + exp(-(t-5)^2)

# ╔═╡ a2fb8c21-c594-44f3-92df-99e6d3afeba1
plot(glucose_pulse, 0, 10, label="G [mmol/L]", xlabel="t")

# ╔═╡ d95198d9-fb96-4389-bf66-f039a673d1a2
md"""
For purpose of solving the ODE problem we will need to define $t$ as a (default) time variable with the command below.
"""

# ╔═╡ 8bb8feaa-63e7-4777-ab79-92d98f84142b
t = default_t();

# ╔═╡ d7d91377-2b6b-4daa-aeb6-d0b20c20cd34
md"""
Redo the ODE problem but now with `:m=>glucose_pulse(t)` as parameter. Use an initial value of $0.0$ for both glucose and insuline concentrations.
"""

# ╔═╡ 6890e3ec-0875-423e-89c1-988411c869e3
# oprob2 = missing;             # Uncomment and complete the instruction
oprob2 = ODEProblem(glucose_insuline_circuit, [:G=>0.0, :I=>0.0], (0., 10.), [:m=>glucose_pulse(t)]);

# ╔═╡ cfcd040b-f762-4cf4-8b9c-1d77ac2837cd
md"""
Solve the new ODE problem using `Tsit5()` and `saveat=0.01`.
"""

# ╔═╡ fb439f9d-ee88-43cd-b91b-44ed8013b27d
# sol2 = missing;                   # Uncomment and complete the instruction
sol2 = solve(oprob2, Tsit5(), saveat=0.01);

# ╔═╡ 83595165-3d71-43b3-93c0-4ba880e0b86e
md"""
Plot the results.
"""

# ╔═╡ 095e359a-2ec4-4637-8477-93510bef23ed
# missing                          # Uncomment and complete the instruction
plot(sol2)

# ╔═╡ 31aaced7-1fe4-4117-80f3-0fa70035822a
md"""
Up to now, we set $s$, the insulin sensitivity to $1$. This parameter represents how sensitive the body is to insuline in taking up glucose. Aging and obesity increase glucose resistance ($1/s$), resulting in diabetes! Explore the effect of this parameter on your plots below.
"""

# ╔═╡ 53addc28-7c87-4c79-ad8f-4a41a3ca5f19
md"""
We make a slider so that s can get values between 0.1 and 5 in step of 0.1.
"""

# ╔═╡ 69e21e9b-7d0e-44e3-8b33-ee7102dadadf
@bind s Slider(0.1:0.1:5.0, default=1, show_value=true)

# ╔═╡ b4073e81-6db8-4be6-a2d1-afa01640dbcf
s

# ╔═╡ 4d58a832-9157-42fe-b578-c3fc29f01ea8
md"""
Below we have made a function that computes/shows the steady state glucose concentration after $100$ hours. Initial values for $G$ and $I$ are $0.0$ and $m$ was set to $0.5$.
"""

# ╔═╡ 7cf11e6f-7ef4-478f-8cc7-eae93374cc58
function glucose_steady_state(s)
	oprob = ODEProblem(glucose_insuline_circuit, [:G=>0.0, :I=>0.0], (0., 100.), [:m=>0.5, :s=>s])
	sol = solve(oprob, Tsit5(), saveat=0.01)
	return sol[:G][end]   # final (steady state) glucose concentration is returned
end

# ╔═╡ 581b82c2-0957-4d20-a664-6ff80b5c40d9
md"""
Calling the function results in the final (steady state) glucose concentration for a specific valu of $s$ set by the slider above.
"""

# ╔═╡ 5eda8c4a-2baa-4987-b408-dbe7383abe95
glucose_steady_state(s)

# ╔═╡ 1c22a730-776e-4508-bdee-d8861c8949e4
md"""
Below is a plot of the steady state value of $G$ as a function of $s$.
"""

# ╔═╡ 34788b1d-1530-4b88-8690-03fe9af52a8b
plot(glucose_steady_state, 0.1, 5, xlabel="s", ylabel="Gss", label="Gss", ylim=(0, 5))

# ╔═╡ f616da2f-b51e-4328-bd3c-47e6d7377cbb
md"""
Use automatic differentiation with `ForwardDiff.derivative(..., ...)` to compute the absolute and relative sensitivity index. Is this system sensitive to the insuline sensitivity $s$?
"""

# ╔═╡ 9be6f5b3-d280-4222-ad64-4b4a7b682310
md"""
Calculate the absolute sensitivity.
"""

# ╔═╡ 97636a73-4732-4377-b760-c36cef13904b
# Uncomment and complete the instruction
# sens_G(s) = ForwardDiff.derivative(..., ...) 
sens_G(s) = ForwardDiff.derivative(glucose_steady_state, s)

# ╔═╡ e0318ae7-d27d-4d73-bf44-632d813e8aa5
md"""
Display the absolute sensitivity for the current $s$ value (cf. slider).
"""

# ╔═╡ a67d70ee-c8b5-4f22-926c-f48c5dcf9815
# missing                     # Uncomment and complete the instruction
sens_G(s)

# ╔═╡ 407437a0-36b1-452f-910c-f244aeeedc4c
md"""
Calculate the normalized (total relative) sensitivity.
"""

# ╔═╡ 1bb5834a-6a4e-40ec-82c1-f2c505ecfeb8
# sens_G_rel(s) = missing     # Uncomment and complete the instruction
sens_G_rel(s) = sens_G(s) * glucose_steady_state(s) / s

# ╔═╡ 66105045-a8f6-4513-a255-a9e1b0c78d37
md"""
Display the normalized (total relative) sensitivity for the current $s$ value (cf. slider).
"""

# ╔═╡ a3156073-07bd-4f13-be5e-d19483fbad28
# missing                     # Uncomment and complete the instruction
sens_G_rel(s)

# ╔═╡ 0ae57436-0b8c-41ce-9d9d-8c48ed59b820
md"""
We see that the final glucose concentration is highly dependent on $s$! This seems to be a flaw in the model, as we can imagine that the physiological parameters can greatly differ from person to person (for example, a person can have a large pancreas). The final glucose concentration should not depend on the insuline sensitivity $s$.

A mechanism that stabilizes this is called *dynamic compensation*. Simply put, we have assumed here that the amount of beta cells ($B$) is fixed. However, in practice, these cells are capable of dividing, growing, and thus producing more insulin. Their growth rate depends on the concentration of glucose, creating an additional feedback loop that stabilizes the physiological circuit.

```julia
μ(G), B --> 2B        # dynamic compensation
```

"""

# ╔═╡ 92d689a4-9f2b-4a0f-81a5-6194af06940c
md"""
Growth rate function depending on the glucose concentration.
"""

# ╔═╡ 87d9c826-56af-409a-a884-d784cfa16a64
μ(G) = 0.3atan(0.5(G-1))

# ╔═╡ 2f784fff-90fc-4273-8b9a-db4593849652
md"""
The growth rate of the $\beta$-cells follows a sigmoid shape, being negative when $G$ is smaller than a threshold and positive if $G$ exceeds this threshold. This curve is plotted below.
"""

# ╔═╡ 4d95e77b-fab0-454b-8319-1fe83e12f9e9
plot(μ, 0, 30, xlab="G", label="μ(G)", title="Glucose-dependend growth rate")

# ╔═╡ ea46f5f9-47a8-4641-a7e3-08a1a29cd04f
md"""
Add dynamic compensation to the model and show that this greatly reduces the sentitivty w.r.t. $s$.
"""

# ╔═╡ 870e9229-374a-40b8-97f7-1a84d1a857e9
md"""
Robust version of a *reaction network object* with dynamic compensation.
```julia
glucose_insuline_circuit_robust = @reaction_network begin
	@parameters q=1 s=1 γ=1 m=1 Ks=1.0
	@species B(t)=1
	m, 0 --> G
	s * I, G --> 0
	B * hill(G, q, Ks, 2), 0 --> I
	γ, I --> 0
	μ(G), B --> 2B
end
```
"""

# ╔═╡ a8814cfd-8fd0-4628-b1ae-80211386a402
md"""
Create the aforementioned robust version of a *reaction network object*.
"""

# ╔═╡ 246a06ba-ef48-47ba-aef0-51f07ceb96e4
# Uncomment and complete the instruction
# glucose_insuline_circuit_robust = @reaction_network begin
# 	@parameters missing
# 	@species missing
# 	missing
# 	...
# end
glucose_insuline_circuit_robust = @reaction_network begin
	@parameters q=1 s=1 γ=1 m=1 Ks=1.0
	@species B(t)=1
	m, 0 --> G
	s * I, G --> 0
	B * hill(G, q, Ks, 2), 0 --> I
	γ, I --> 0
	μ(G), B --> 2B
end

# ╔═╡ 69965b61-043b-44a0-998b-f5f34deefdd2
md"""
Convert the system to a symbolic differential equation model and inspect the differential equations. You do not need to make a new variable, just call `convert` with the right arguments.
"""

# ╔═╡ 2f992357-4608-40e6-9954-8306f31e9768
# missing                         # Uncomment and complete the instruction
convert(ODESystem, glucose_insuline_circuit_robust)

# ╔═╡ f87451c3-1707-4f13-a373-d8949ce4ca08
md"""
Simulate the system over a time interval of $0$ to $10$ hours with default parameter values, and initial glucose and insuline concentrations of $5.0$ and $0.0$, repectively. Use `Tsit5()` and `saveat=0.01` to solve.
"""

# ╔═╡ 93457418-a032-48b6-9c63-94ce49645c68
# oprob_robust = missing;         # Uncomment and complete the instruction
oprob_robust = ODEProblem(glucose_insuline_circuit_robust, [:G=>5.0, :I=>0.0], (0., 10.), []);

# ╔═╡ cac1644a-7604-4716-85c8-63b657ff95b2
# sol_robust = missing;           # Uncomment and complete the instruction
sol_robust = solve(oprob_robust, Tsit5(), saveat=0.01);

# ╔═╡ 06645013-b63b-49b4-822b-f08059f3a394
md"""
Plot the results. Use thereby `ylim=(0.0, 5.5)`.
"""

# ╔═╡ de7919d3-04a8-4b4f-ad76-2f835275fc6b
# missing                            # Uncomment and complete the instruction
plot(sol_robust, ylim=(0.0, 5.5))

# ╔═╡ 318e83fc-340b-4fd1-9e72-dbfdb891ff5f
md"""
Implement a function that computes/shows the steady state glucose concentration after $100$ hours. Set initial values for $G$ and $I$ to $0.0$ and set $m$ to $0.5$.

- Tip: copy the *body* of the former function `glucose_steady_state` and adapt.
"""

# ╔═╡ bcf597ef-c57f-402d-8209-511a7fc16437
# Uncomment and complete the instruction
# function glucose_steady_state_robust(s)
# 	oprob = missing
# 	sol = missing
# 	return missing
# end
function glucose_steady_state_robust(s)
	oprob = ODEProblem(glucose_insuline_circuit_robust, [:G=>0.0, :I=>0.0], (0., 100.), [:m=>1/2, :s=>s])
	sol = solve(oprob, Tsit5(), saveat=0.01)
	return sol[:G][end]
end

# ╔═╡ 5c2d0d0c-7bce-4696-9807-e793b6f8033e
md"""
Create a new slider object with a range between $0.1$ and $5.0$ and step size $0.1$, and bind it to the new variable `s_robust`.

- Tip: copy the previous slider and adapt.
"""

# ╔═╡ 8ef6d8a4-be72-49dc-be44-0e4eb06336be
# missing                   # Uncomment and complete the instruction
@bind s_robust Slider(0.1:0.1:5, default=1, show_value=true)

# ╔═╡ fd56b95e-59a3-46e7-8216-e50181c22a58
md"""
Call the function `glucose_steady_state_robust` with `s_robust` as argument and observe the new steady state glucose concentration for different insulin sensitivity values.
"""

# ╔═╡ daf82017-c427-480c-94d2-938e675004f0
# missing                   # Uncomment and complete the instruction
glucose_steady_state_robust(s_robust)

# ╔═╡ b95899e9-3af9-4fc5-820b-f72164beb5db
md"""
Plot of the new steady state value of $G$ as a function of $s$ (in the range $[0.1, 5.0]$). You might need to use `ylim=(0.98, 1.02)`.
"""

# ╔═╡ 4adeac81-8724-4a2c-a004-17be8f92f95e
# missing                     # Uncomment and complete the instruction
plot(glucose_steady_state_robust, 0.1, 5, xlab="s", ylab="Gss", ylim=(0.98, 1.02))

# ╔═╡ 041f6717-2709-47cb-ad54-8ef6b1a61fcc
md"""
Use automatic differentiation with `ForwardDiff.derivative(..., ...)` to compute the normalized (total relative) sensitivity index. Is this new system sensitive to the insuline sensitivity $s$?

- Answer: missing
"""

# ╔═╡ 2180719e-ad6b-4b63-a835-3c30b811d35a
# sens_G_rel_robust(s) = ...      # Uncomment and complete the instruction
sens_G_rel_robust(s) = ForwardDiff.derivative(glucose_steady_state_robust, s) * glucose_steady_state_robust(s) / s

# ╔═╡ d4af98b9-40dc-486b-bbc0-9538a193cf19
md"""
Display the new normalized (total relative) sensitivity for the current $s$ value (cf. slider).
"""

# ╔═╡ 6ca12c8b-a535-43ca-9d6a-86e324eac51e
sens_G_rel_robust(s_robust)

# ╔═╡ Cell order:
# ╠═52b28a4c-b0bb-11ef-2841-17ecfb596676
# ╠═7bc363b0-9415-4954-807f-81a308bde531
# ╟─bfcc4b4e-073e-401e-851c-d01ef828028a
# ╟─a5789d9d-a690-4c6b-8b8d-fbfe63f659d0
# ╠═e93bc1a5-e0b1-425c-b069-d5ce80756e60
# ╟─be33a8a3-a3f5-448d-9390-16f7c4ca099d
# ╠═2ac43184-a393-4d39-9c3b-c299c3371b09
# ╟─75a5b5d9-f8d7-4728-83dc-03dde502dcfb
# ╠═f62c2b4a-a9e8-472b-b07b-44072da08ee8
# ╟─882e106f-d3a5-4ccb-ae65-cfb5e4e5de57
# ╠═99ccbd1b-b9e4-4afc-8d9a-8c1055bcce0f
# ╟─3c1465d2-0886-4b74-baee-55c7dad41f9a
# ╠═5bac69d1-0bd4-41bd-b15d-5ca6292a31e4
# ╠═6ccd5ecc-5579-4ce5-9425-d1bbae0bc5ba
# ╠═afc40964-0a97-43a6-ab27-5fe96ecd319b
# ╠═f050b32e-1a14-4d5e-9fb3-5dbc4dfbf79c
# ╟─cfbbfe48-a2c3-4691-9fac-b0a18f99222c
# ╠═7e35b56b-44a0-46df-a812-f5fb437b672e
# ╟─0133ab41-2c8c-4147-9594-1ad238833a9e
# ╟─479e0f11-e18e-4523-8d24-b2026acdeaef
# ╠═824677b4-6533-453f-909f-bf04fe4bf5d9
# ╟─cef7f30d-43ef-4f2c-8bf3-35c140d237ce
# ╠═9d3ba8bd-e9c3-4bb6-957d-c7f8211a2cda
# ╟─3f5e2195-da93-40fb-88ec-9e85a7bb6c24
# ╠═3dae761f-5fb6-4f1c-a440-d2c21a6edb17
# ╟─314aa77c-5906-4772-94c5-f0b7cf648de0
# ╠═3b56f88c-0ca8-4f33-8516-11e639cb6a1a
# ╠═00db9d76-c4dd-450d-b7e3-cf87f8bfeffb
# ╠═5d1134d5-734e-4851-a91c-615dc37e3c68
# ╟─aab02fef-6d84-4f53-aff5-6f94bf2bd460
# ╠═8a5ce760-55bb-4610-9226-32c3d8376748
# ╠═a2fb8c21-c594-44f3-92df-99e6d3afeba1
# ╟─d95198d9-fb96-4389-bf66-f039a673d1a2
# ╠═8bb8feaa-63e7-4777-ab79-92d98f84142b
# ╟─d7d91377-2b6b-4daa-aeb6-d0b20c20cd34
# ╠═6890e3ec-0875-423e-89c1-988411c869e3
# ╟─cfcd040b-f762-4cf4-8b9c-1d77ac2837cd
# ╠═fb439f9d-ee88-43cd-b91b-44ed8013b27d
# ╟─83595165-3d71-43b3-93c0-4ba880e0b86e
# ╠═095e359a-2ec4-4637-8477-93510bef23ed
# ╟─31aaced7-1fe4-4117-80f3-0fa70035822a
# ╟─53addc28-7c87-4c79-ad8f-4a41a3ca5f19
# ╠═69e21e9b-7d0e-44e3-8b33-ee7102dadadf
# ╠═b4073e81-6db8-4be6-a2d1-afa01640dbcf
# ╟─4d58a832-9157-42fe-b578-c3fc29f01ea8
# ╠═7cf11e6f-7ef4-478f-8cc7-eae93374cc58
# ╟─581b82c2-0957-4d20-a664-6ff80b5c40d9
# ╠═5eda8c4a-2baa-4987-b408-dbe7383abe95
# ╟─1c22a730-776e-4508-bdee-d8861c8949e4
# ╠═34788b1d-1530-4b88-8690-03fe9af52a8b
# ╟─f616da2f-b51e-4328-bd3c-47e6d7377cbb
# ╟─9be6f5b3-d280-4222-ad64-4b4a7b682310
# ╠═97636a73-4732-4377-b760-c36cef13904b
# ╟─e0318ae7-d27d-4d73-bf44-632d813e8aa5
# ╠═a67d70ee-c8b5-4f22-926c-f48c5dcf9815
# ╟─407437a0-36b1-452f-910c-f244aeeedc4c
# ╠═1bb5834a-6a4e-40ec-82c1-f2c505ecfeb8
# ╟─66105045-a8f6-4513-a255-a9e1b0c78d37
# ╠═a3156073-07bd-4f13-be5e-d19483fbad28
# ╟─0ae57436-0b8c-41ce-9d9d-8c48ed59b820
# ╟─92d689a4-9f2b-4a0f-81a5-6194af06940c
# ╠═87d9c826-56af-409a-a884-d784cfa16a64
# ╟─2f784fff-90fc-4273-8b9a-db4593849652
# ╠═4d95e77b-fab0-454b-8319-1fe83e12f9e9
# ╟─ea46f5f9-47a8-4641-a7e3-08a1a29cd04f
# ╟─870e9229-374a-40b8-97f7-1a84d1a857e9
# ╟─a8814cfd-8fd0-4628-b1ae-80211386a402
# ╠═246a06ba-ef48-47ba-aef0-51f07ceb96e4
# ╟─69965b61-043b-44a0-998b-f5f34deefdd2
# ╠═2f992357-4608-40e6-9954-8306f31e9768
# ╟─f87451c3-1707-4f13-a373-d8949ce4ca08
# ╠═93457418-a032-48b6-9c63-94ce49645c68
# ╠═cac1644a-7604-4716-85c8-63b657ff95b2
# ╟─06645013-b63b-49b4-822b-f08059f3a394
# ╠═de7919d3-04a8-4b4f-ad76-2f835275fc6b
# ╟─318e83fc-340b-4fd1-9e72-dbfdb891ff5f
# ╠═bcf597ef-c57f-402d-8209-511a7fc16437
# ╟─5c2d0d0c-7bce-4696-9807-e793b6f8033e
# ╠═8ef6d8a4-be72-49dc-be44-0e4eb06336be
# ╟─fd56b95e-59a3-46e7-8216-e50181c22a58
# ╠═daf82017-c427-480c-94d2-938e675004f0
# ╟─b95899e9-3af9-4fc5-820b-f72164beb5db
# ╠═4adeac81-8724-4a2c-a004-17be8f92f95e
# ╟─041f6717-2709-47cb-ad54-8ef6b1a61fcc
# ╠═2180719e-ad6b-4b63-a835-3c30b811d35a
# ╟─d4af98b9-40dc-486b-bbc0-9538a193cf19
# ╠═6ca12c8b-a535-43ca-9d6a-86e324eac51e
