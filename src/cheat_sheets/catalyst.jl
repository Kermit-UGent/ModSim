### A Pluto.jl notebook ###
# v0.20.4

#> [frontmatter]
#> order = "4"
#> title = "Catalyst Cheat Sheet"
#> date = "2025-01-26"
#> tags = ["welcome"]
#> description = "Catalyst Cheat Sheet"
#> layout = "layout.jlhtml"
#> 
#>     [[frontmatter.author]]
#>     name = "Daan Van Hauwermeiren"
#>     [[frontmatter.author]]
#>     name = "Michiel Stock"

using Markdown
using InteractiveUtils

# ╔═╡ c3cccc9c-dbfc-11ef-31c5-07bc34711b2d
using Pkg; Pkg.activate("../../pluto-deployment-environment")

# ╔═╡ fe5e5f00-638f-4ca3-a331-c746af53b6b4
using Catalyst, Plots, DifferentialEquations

# ╔═╡ 22cc8f5e-da25-48bb-8cbc-952e205d5640
md"""## Notes before the cheat sheet

In Pluto, there can only be one statement per cell (because of the syntax tree that is generated to determine the order of execution). But sometimes we want to group multiple expressions for clarity. We can do that in two ways, using the same syntax, but with a different keyword: both call for wrapping the code in a block.

The first one is for simply wrapping multiple statements (the variable names can be accessed in the rest of the notebook):
```julia
begin
	...
end
```
In the second one, the variables only live within the scope of the block. We will use this to illustrate behaviour, but we do not need the rest in any other part of the notebook, and we want to avoid having to add numerous postfixes to the variable names.
```julia
let
    ...
end
```

Note that the indentation is optional, and only used for readability.
"""


# ╔═╡ df1d7710-f666-46c4-9bca-2707a997eb83
md"# `Catalyst` cheat sheet"

# ╔═╡ 19b4c4a9-8e6c-4a4e-ae8f-3506ef7d15ec
md"""
`Catalyst.jl` is a Julia package that provides a clean interface to
building reaction networks, which can be turned into ODESystems
that `DifferentialEquations.jl` can simulate. Since almost all mass
transfer problems can be written using reactions (or, more generally
speaking, processes), `Catalyst.jl` will be our main tool for building
mechanistic models.
"""

# ╔═╡ 6c756d57-075b-446d-b88b-9934a57179bd
md"## Defining a reaction system"

# ╔═╡ 20f23ad3-ce33-4e12-b3fb-b9376b1a4962
md"""
Use the `@reaction_network` macro to define a reaction network.
This macro allows you to specify reactions using a simple syntax.
"""

# ╔═╡ dc557b21-e296-446b-b7ab-0462707569c6
mm = @reaction_network begin
    (kB, kD), S + E <--> ES  # reversible binding
    kP, ES --> P + E         # conversion of substrate by enzyme
end

# ╔═╡ f62bff24-00a4-4b7f-97ef-c2dce24b70fe
species(mm) # check the species (states)

# ╔═╡ 063a77c2-5404-4729-b1c1-679c3326d636
parameters(mm) # check the parameters

# ╔═╡ 96703b25-43f5-4fa0-8283-c43fc5a670f2
equations(mm) # check the equations

# ╔═╡ ea2a4e2b-a546-442c-abbe-43dc52f68c85
let
	@unpack S = mm  # extract parameter / variable
end

# ╔═╡ 10e37455-df3d-443c-ac22-529c63624b59
osys = convert(ODESystem, mm) # convert ReactionSystem in an ODE system

# ╔═╡ b4ead481-5cdd-4f5d-a995-d034cc3a51d5
md"""
It is also possible to add reactions using functions, such as building big reaction networks using code. For this, we refer to the documentations.
"""

# ╔═╡ 7b2dbf80-68ed-47ad-b324-a4fdb003a765
md"## Simulation"

# ╔═╡ e88f863f-1d13-4ebd-8169-7d88a3fde1c2
md"""
Simulation of a reaction network builds upon DifferentialEquations.jl. Reaction networks can directly be transformed in ODE systems (if you want to see the differential equations) or a problem (needed to solve numerically).
"""

# ╔═╡ 83a462e1-15be-4395-b6a5-55214fcfa6b7
u0map = [:S => 10.0, :E => 0.1, :P => 0.0, :ES => 0] # define intial values

# ╔═╡ e388f03c-37a3-4841-bc36-1d1c506ee941
pmap = [:kB => 0.5, :kD => 0.1, :kP => 2.2] # define parameters

# ╔═╡ cbcd8520-20ce-4585-a03c-83572d996d10
let
	# alternative, just as good!
	u0map = [mm.S => 10.0, mm.E => 0.1, mm.P => 0.0, mm.ES => 0]
	pmap = [mm.kB => 0.5, mm.kD => 0.1, mm.kP => 2.2]
end

# ╔═╡ 0a6d9520-acb7-4bc1-9f90-65738a4d0b1c
tspan = (0.0, 100.)

# ╔═╡ 4f4c1e47-c65c-48b9-b57f-3ba0a031b95d
oprob = ODEProblem(mm, u0map, tspan, pmap)

# ╔═╡ c4420960-a3fb-4d0d-a095-3a0d656236f7
sol = solve(oprob, Tsit5())

# ╔═╡ 0a2632e7-8942-4007-a8a4-38cedd53eb1d
plot(sol)  # plot all variables

# ╔═╡ a8abeb40-cd52-4193-b484-4a56b6cda714
plot(sol, idxs=[:S, :P])  # plot only S and P

# ╔═╡ a6397a2a-5a81-43ac-9dad-8c8766872eec
begin
	# alternative, unpack variables
	@unpack P, S = mm  # extract product and substrate
	plot(sol, idxs=[P, S])
end

# ╔═╡ a6eec47f-bf65-4fa6-ac1e-8c4f3a5d7ea4
plot(sol, idxs=P/(S + P), title="Fraction of substrate converted")

# ╔═╡ 574a345b-350c-41dc-8bd7-680b2055e3c5
md"""
You see that reaction systems have access to the variables and parameters, which are also used for plotting.
"""

# ╔═╡ bb7827c6-2f44-49ba-bbb2-043f70a115ae
md"## Using default options"

# ╔═╡ 3d4a712c-0729-4112-897c-bbf2be20bb76
md"You can already specify initial values and parameters directly in the reaction network."

# ╔═╡ 78c45bff-a134-4bc6-9003-31fe016e01cc
mm_defaults = @reaction_network begin
	@species S(t)=10. E(t)=0.1 P(t)=0 ES(t)=0
	@parameters kB=0.5 kD=0.1 kP=2.2
	(kB, kD), S + E <--> ES
	kP, ES --> P + E
end

# ╔═╡ cd4c918d-bca1-4fda-89a1-b980311aabc9
oprob_defaults = ODEProblem(mm_defaults, [], tspan); 
	# no need to specify initial states and parameters

# ╔═╡ ccd1f502-fdd4-44dd-985a-dfde703e1154
plot(solve(oprob_defaults, Tsit5()))

# ╔═╡ 8ba6c0b8-1dcc-4a38-b040-a155e3c63414
ODEProblem(mm_defaults, [:E=>0.2], tspan, [:kP=>1.5]);
	# overwriting defaults

# ╔═╡ b64be898-7bd2-4266-a28b-425c45749df2
mm_ann = @reaction_network begin
	@species S(t)=10. [description="substrate"] P(t)=0 [description="product"]
	@parameters kB=0.5 [description="binding rate"] kD=0.1 [description="dissociated rate rate"] kP=2.2 [description="conversion rate"]
	(kB, kD), S + E <--> ES
	kP, ES --> P + E
end
	# adding annotation to the parameters

# ╔═╡ 14ae384c-4ebd-4264-8e27-5ee7c00fd8fb
md"## Observables"

# ╔═╡ e8e66e87-f678-4ac1-9683-75a7330eb0ad
md"""
Often, we are interested in the states or species in the system. However, sometimes we want to track something that is not a state but computed based on the states. This is an observable. Even though you can always compute these afterward, adding them to the system is likely useful so you have access to them while plotting. 

Consider the total amount of enzyme in the system. It is clear that
this is conserved.
"""

# ╔═╡ 773e72a8-1041-428b-9ae4-96da03aac835
mm_obs = @reaction_network begin
    @observables begin
        Etot ~ E + ES
    end
    (kB, kD), S + E <--> ES
    kP, ES --> P + E
end

# ╔═╡ 23cd6fbb-4edb-4498-b97b-6e6cf2db8ad9
observed(mm_obs)

# ╔═╡ fada3709-04d3-4794-b3ee-44e4e9583860
oprob_obs = ODEProblem(mm_obs, u0map, tspan, pmap);

# ╔═╡ 54355a46-6d7c-4a18-9aa6-01428d2ed69a
sol_obs = solve(oprob_obs, Tsit5());

# ╔═╡ ed65589e-4618-4496-882d-411bd272d509
plot(sol_obs, idxs=:Etot, ylims=(0,1))
	# flat line, enzyme is conserved

# ╔═╡ 8ac9f9a3-ecbb-4f25-a7b6-88f96e256f8e
md"## Removing conserved quantities"

# ╔═╡ c8a14bda-547f-40b6-be81-a6d5689a044f
equations(osys)

# ╔═╡ 778232f3-3a5d-42f9-b698-941328866833
mm_conserved = convert(ODESystem, mm, remove_conserved = true)

# ╔═╡ 183c9ba5-9dbb-42b9-9800-3b3e99ef6890
equations(mm_conserved) # enzyme (E + ES) is conserved

# ╔═╡ db0390cc-8a4c-4bd5-a911-93e009c21fa8
observed(mm_conserved) # still in observables

# ╔═╡ 702d7e17-4d7f-4401-873c-c7d7ea04b450
md"## Events"

# ╔═╡ c8050c1d-4684-440d-9b97-0f371b8d548b
md"""
Events are external perturbations of the system, either by changing the states or the parameters. There are two types of events:
- discrete events: which happen at fixed time points;
- continuous events: which happen when a variable reaches a certain value.
"""

# ╔═╡ b5501abe-088b-44c5-96c0-9b732f691e4a
md"### Discrete Events (time-based)"

# ╔═╡ 584e3f6a-ad10-4292-bb05-dd66534d4621
dilluting = [20.0] => [
	mm.E ~ mm.E/2, mm.S ~ mm.S/2, mm.ES ~ mm.ES / 2, mm.P ~ mm.P/2
]
	# at t=20, we double the volume of the reactor, hence halving the concentration

# ╔═╡ 7ef543ea-2974-4576-a2f4-ca876590b2e2
@named mm_dil = ReactionSystem(equations(mm), discrete_events=dilluting);

# ╔═╡ 0373d49d-c2fa-4a75-8826-62b389c74e02
mm_dillute = complete(mm_dil);

# ╔═╡ 97dc4594-b3e3-4392-b107-5d71612378c4
oprob_dil = ODEProblem(mm_dillute, u0map, tspan, pmap);

# ╔═╡ b876f4b4-e25d-4fdf-b87c-27058facc68d
sol_dil = solve(oprob_dil, Tsit5(), tstops=20);

# ╔═╡ c9b59f48-8585-45d9-91e9-d9b3ef024736
plot(sol_dil)

# ╔═╡ b62081ba-aa94-4403-a300-1b08234d36c3
md"### Continuous events"

# ╔═╡ 496f16b2-f776-456b-9920-3206afe18b6c
substrate_feeding = [mm.S ~ 5] => [mm.S ~ mm.S + 5]
	# when the substrate reaches a concentration of 5, we add new substrate

# ╔═╡ 0a3f23a6-6a1b-482f-ab45-dbc91a99eafb
@named mm_f = ReactionSystem(equations(mm), continuous_events=substrate_feeding);

# ╔═╡ 2de7223f-d2dc-4d64-bff4-fac0677c083e
mm_fed = complete(mm_f);

# ╔═╡ df317585-fcab-4898-a2d3-3f0d73a663db
oprob_fed = ODEProblem(mm_fed, u0map, tspan, pmap);

# ╔═╡ 2a1b3ea1-bc2b-437a-822a-7199d532007c
sol_fed = solve(oprob_fed, Tsit5(), tstops=20);

# ╔═╡ 6641a5eb-5bac-4442-8610-710c1ea826bc
plot(sol_fed)

# ╔═╡ 44fa725e-f4e5-43db-8647-4fa9f4e5e24b
md"## Computing steady state"

# ╔═╡ 70d30f17-e919-4651-926c-09a3f9f38567
md"It is possible to numerically compute when the system reaches a steady state."

# ╔═╡ 529b385b-738c-4601-9fc7-ff3c678656ad
mm_continuous = @reaction_network begin
	    (kB, kD), S + E <--> ES  # reversible binding
	    kP, ES --> P + E         # conversion of substrate by enzyme
	    1, 0 --> S  # constant inflow
	    1, (S, P) --> 0  # constant outflow
end;

# ╔═╡ c85dd29e-2eac-4b1a-b23d-a84ff36b7853
ssprob = SteadyStateProblem(mm_continuous, u0map, pmap);

# ╔═╡ bb9f02a8-c4db-45d2-9dd7-163037d8c91d
print(solve(ssprob))

# ╔═╡ 88c50285-6073-4532-a28d-d3373a2b632e
md"## Discrete jump equations"

# ╔═╡ ba0eadd7-41ed-4558-8323-d4f386aabe3a
md"""
When modelling discrete systems, we can turn the problem into a discrete jump equation. Here, the species are described by an integer; there are a discrete number of molecules. This simulation is now stochastic.
"""

# ╔═╡ fcd24970-b0ca-4d4e-809f-ea9078cc81b7
u0map_int = [:S => 500, :E => 5, :P => 0, :ES => 0];
	# starting with 500 substrate molecules and 5 enzyme molecules

# ╔═╡ 09932682-b17b-4fea-9e0d-cc29cd6f63d6
jsys = JumpInputs(mm, u0map_int, tspan, pmap);

# ╔═╡ e0de6197-a122-40f3-9bc4-c398eaf5e713
jprob = JumpProblem(jsys);

# ╔═╡ 74bd42b4-c866-4b28-ac87-493db87a7ed6
jsol = solve(jprob);

# ╔═╡ b762e40b-b208-43dd-8395-f241ad3a875b
plot(jsol)

# ╔═╡ 952ee553-567c-4a70-bb24-5db4c7043395
md"## Stochastic Differential Equations"

# ╔═╡ 7c98aa51-0853-40f4-8cb1-fe763cda8698
md"""
We can also simulate a reaction system as an SDE. Be wary that no mechanism prevents concentrations from being zero!
"""

# ╔═╡ aae67739-e1ac-4949-8924-087a63881118
let
	sys = @reaction_network begin
	    (k1, k2), A <--> B
	end
	
	u0map = [:A => 10., :B=> 200.]
	pmap = [:k1 => 2, :k2=> 5]
	
	
	sprob = SDEProblem(sys, u0map, (0., 20.), pmap)
	
	sol = solve(sprob, STrapezoid(), dt=0.02)
	plot(sol)
end

# ╔═╡ 73632fbb-88f7-4ff4-a7e6-b5e32780764c
let
	sys = @reaction_network begin
	    @parameters η
	    @default_noise_scaling η
	    (k1, k2), A <--> B
	end
	
	pmap = [:k1 => 2, :k2=> 5, :η=>0.1]
	u0map = [:A => 10., :B=> 200.]
	
	sprob = SDEProblem(sys, u0map, (0., 20.), pmap)
	sol = solve(sprob, STrapezoid(), dt=0.02)
	plot(sol)
end

# ╔═╡ 5e0acdd7-318b-4f7e-b4d3-051acc423f69
let
	sys = @reaction_network begin
	    @parameters η
	    @default_noise_scaling η
	    k1, A --> B, [noise_scaling = 0.0]
	    k2, B --> A, [noise_scaling = η]
	end
	
	pmap = [:k1 => 2, :k2=> 5, :η=>1]
	u0map = [:A => 10., :B=> 200.]
	sprob = SDEProblem(sys, u0map, (0., 20.), pmap)
	sol = solve(sprob, STrapezoid(), dt=0.02)
	plot(sol)
end

# ╔═╡ Cell order:
# ╟─22cc8f5e-da25-48bb-8cbc-952e205d5640
# ╟─df1d7710-f666-46c4-9bca-2707a997eb83
# ╠═c3cccc9c-dbfc-11ef-31c5-07bc34711b2d
# ╠═fe5e5f00-638f-4ca3-a331-c746af53b6b4
# ╟─19b4c4a9-8e6c-4a4e-ae8f-3506ef7d15ec
# ╟─6c756d57-075b-446d-b88b-9934a57179bd
# ╟─20f23ad3-ce33-4e12-b3fb-b9376b1a4962
# ╠═dc557b21-e296-446b-b7ab-0462707569c6
# ╠═f62bff24-00a4-4b7f-97ef-c2dce24b70fe
# ╠═063a77c2-5404-4729-b1c1-679c3326d636
# ╠═96703b25-43f5-4fa0-8283-c43fc5a670f2
# ╠═ea2a4e2b-a546-442c-abbe-43dc52f68c85
# ╠═10e37455-df3d-443c-ac22-529c63624b59
# ╟─b4ead481-5cdd-4f5d-a995-d034cc3a51d5
# ╟─7b2dbf80-68ed-47ad-b324-a4fdb003a765
# ╟─e88f863f-1d13-4ebd-8169-7d88a3fde1c2
# ╠═83a462e1-15be-4395-b6a5-55214fcfa6b7
# ╠═e388f03c-37a3-4841-bc36-1d1c506ee941
# ╠═cbcd8520-20ce-4585-a03c-83572d996d10
# ╠═0a6d9520-acb7-4bc1-9f90-65738a4d0b1c
# ╠═4f4c1e47-c65c-48b9-b57f-3ba0a031b95d
# ╠═c4420960-a3fb-4d0d-a095-3a0d656236f7
# ╠═0a2632e7-8942-4007-a8a4-38cedd53eb1d
# ╠═a8abeb40-cd52-4193-b484-4a56b6cda714
# ╠═a6397a2a-5a81-43ac-9dad-8c8766872eec
# ╠═a6eec47f-bf65-4fa6-ac1e-8c4f3a5d7ea4
# ╟─574a345b-350c-41dc-8bd7-680b2055e3c5
# ╟─bb7827c6-2f44-49ba-bbb2-043f70a115ae
# ╟─3d4a712c-0729-4112-897c-bbf2be20bb76
# ╠═78c45bff-a134-4bc6-9003-31fe016e01cc
# ╠═cd4c918d-bca1-4fda-89a1-b980311aabc9
# ╠═ccd1f502-fdd4-44dd-985a-dfde703e1154
# ╠═8ba6c0b8-1dcc-4a38-b040-a155e3c63414
# ╠═b64be898-7bd2-4266-a28b-425c45749df2
# ╟─14ae384c-4ebd-4264-8e27-5ee7c00fd8fb
# ╟─e8e66e87-f678-4ac1-9683-75a7330eb0ad
# ╠═773e72a8-1041-428b-9ae4-96da03aac835
# ╠═23cd6fbb-4edb-4498-b97b-6e6cf2db8ad9
# ╠═fada3709-04d3-4794-b3ee-44e4e9583860
# ╠═54355a46-6d7c-4a18-9aa6-01428d2ed69a
# ╠═ed65589e-4618-4496-882d-411bd272d509
# ╟─8ac9f9a3-ecbb-4f25-a7b6-88f96e256f8e
# ╠═c8a14bda-547f-40b6-be81-a6d5689a044f
# ╠═778232f3-3a5d-42f9-b698-941328866833
# ╠═183c9ba5-9dbb-42b9-9800-3b3e99ef6890
# ╠═db0390cc-8a4c-4bd5-a911-93e009c21fa8
# ╟─702d7e17-4d7f-4401-873c-c7d7ea04b450
# ╟─c8050c1d-4684-440d-9b97-0f371b8d548b
# ╟─b5501abe-088b-44c5-96c0-9b732f691e4a
# ╠═584e3f6a-ad10-4292-bb05-dd66534d4621
# ╠═7ef543ea-2974-4576-a2f4-ca876590b2e2
# ╠═0373d49d-c2fa-4a75-8826-62b389c74e02
# ╠═97dc4594-b3e3-4392-b107-5d71612378c4
# ╠═b876f4b4-e25d-4fdf-b87c-27058facc68d
# ╠═c9b59f48-8585-45d9-91e9-d9b3ef024736
# ╟─b62081ba-aa94-4403-a300-1b08234d36c3
# ╠═496f16b2-f776-456b-9920-3206afe18b6c
# ╠═0a3f23a6-6a1b-482f-ab45-dbc91a99eafb
# ╠═2de7223f-d2dc-4d64-bff4-fac0677c083e
# ╠═df317585-fcab-4898-a2d3-3f0d73a663db
# ╠═2a1b3ea1-bc2b-437a-822a-7199d532007c
# ╠═6641a5eb-5bac-4442-8610-710c1ea826bc
# ╟─44fa725e-f4e5-43db-8647-4fa9f4e5e24b
# ╟─70d30f17-e919-4651-926c-09a3f9f38567
# ╠═529b385b-738c-4601-9fc7-ff3c678656ad
# ╠═c85dd29e-2eac-4b1a-b23d-a84ff36b7853
# ╠═bb9f02a8-c4db-45d2-9dd7-163037d8c91d
# ╟─88c50285-6073-4532-a28d-d3373a2b632e
# ╟─ba0eadd7-41ed-4558-8323-d4f386aabe3a
# ╠═fcd24970-b0ca-4d4e-809f-ea9078cc81b7
# ╠═09932682-b17b-4fea-9e0d-cc29cd6f63d6
# ╠═e0de6197-a122-40f3-9bc4-c398eaf5e713
# ╠═74bd42b4-c866-4b28-ac87-493db87a7ed6
# ╠═b762e40b-b208-43dd-8395-f241ad3a875b
# ╟─952ee553-567c-4a70-bb24-5db4c7043395
# ╟─7c98aa51-0853-40f4-8cb1-fe763cda8698
# ╠═aae67739-e1ac-4949-8924-087a63881118
# ╠═73632fbb-88f7-4ff4-a7e6-b5e32780764c
# ╠═5e0acdd7-318b-4f7e-b4d3-051acc423f69
