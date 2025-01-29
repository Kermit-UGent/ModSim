### A Pluto.jl notebook ###
# v0.20.4

#> [frontmatter]
#> order = "3"
#> title = "Catalyst"
#> date = "2025-01-26"
#> tags = ["cheat sheets"]
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
begin
	using Pkg
	Pkg.activate("../../pluto-deployment-environment")
	using Catalyst, Plots, DifferentialEquations
end

# ╔═╡ df1d7710-f666-46c4-9bca-2707a997eb83
md"# The `Catalyst` cheat sheet"

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


# ╔═╡ 6c756d57-075b-446d-b88b-9934a57179bd
md"## defining a reaction system"

# ╔═╡ dc557b21-e296-446b-b7ab-0462707569c6
mm = @reaction_network begin
    (kB, kD), S + E <--> ES  # reversible binding
    kP, ES --> P + E         # conversion of substrate by enzyme
end

# ╔═╡ f62bff24-00a4-4b7f-97ef-c2dce24b70fe
species(mm)  # check the species

# ╔═╡ 063a77c2-5404-4729-b1c1-679c3326d636
parameters(mm)  # check the parameters

# ╔═╡ 96703b25-43f5-4fa0-8283-c43fc5a670f2
equations(mm)  # check the equationd

# ╔═╡ ea2a4e2b-a546-442c-abbe-43dc52f68c85
let
	@unpack S = mm  # extract parameter / variable
end

# ╔═╡ 10e37455-df3d-443c-ac22-529c63624b59
osys = convert(ODESystem, mm) # convert ReactionSystem in an ODE system

# ╔═╡ 7b2dbf80-68ed-47ad-b324-a4fdb003a765
md"## Simulation"

# ╔═╡ de0ac799-a461-4009-a923-df79340baf67
begin
	# define paramters and intial values
	u0map = [:S => 10.0, :E => 0.1, :P => 0.0, :ES => 0]
	pmap = [:kB => 0.5, :kD => 0.1, :kP => 2.2]
end

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

# ╔═╡ bb7827c6-2f44-49ba-bbb2-043f70a115ae
md"## Default options"

# ╔═╡ 2fa02aa4-d035-49e0-b19f-85fe65587914
begin
	mm_defaults = @reaction_network begin
	    @species S(t)=10. E(t)=0.1 P(t)=0 ES(t)=0
	    @parameters kB=0.5 kD=0.1 kP=2.2
	    (kB, kD), S + E <--> ES
	    kP, ES --> P + E
	end
	
	# no need to specify initial states and parameters
	oprob_defaults = ODEProblem(mm_defaults, [], tspan)
	
	plot(solve(oprob_defaults, Tsit5()))
end

# ╔═╡ 8ba6c0b8-1dcc-4a38-b040-a155e3c63414
# overwriting defaults
ODEProblem(mm_defaults, [:E=>0.2], tspan, [:kP=>1.5])

# ╔═╡ b64be898-7bd2-4266-a28b-425c45749df2
let
	# adding annotation to the parameters
	mm = @reaction_network begin
	    @species S(t)=10. [description="substrate"] P(t)=0 [description="product"]
	    @parameters kB=0.5 [description="binding rate"] kD=0.1 [description="dissociated rate rate"] kP=2.2 [description="conversion rate"]
	    (kB, kD), S + E <--> ES
	    kP, ES --> P + E
	end
end

# ╔═╡ 14ae384c-4ebd-4264-8e27-5ee7c00fd8fb
md"## Observables"

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

# ╔═╡ 6223afa4-c38f-46bc-b6e4-f5eff8a9efdc
let
	oprob = ODEProblem(mm_obs, u0map, tspan, pmap)

	sol = solve(oprob, Tsit5())
	
	plot(sol, idxs=:Etot, ylims=(0,1))
end

# ╔═╡ 8ac9f9a3-ecbb-4f25-a7b6-88f96e256f8e
md"## removing conserved quantities"

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

# ╔═╡ b5501abe-088b-44c5-96c0-9b732f691e4a
md"### Discrete Events"

# ╔═╡ fe8e9068-3b21-4450-a3c2-95833093a1f1
let 
	dilluting = [20.0] => [mm.E ~ mm.E/2, mm.S ~ mm.S/2,
                            mm.ES ~ mm.ES / 2, mm.P ~ mm.P/2]

	@named mm_dillute = ReactionSystem(equations(mm), discrete_events=dilluting)
	
	mm_dillute = complete(mm_dillute)
	
	oprob = ODEProblem(mm_dillute, u0map, tspan, pmap)
	
	sol = solve(oprob, Tsit5(), tstops=20)
	
	plot(sol)
end

# ╔═╡ b62081ba-aa94-4403-a300-1b08234d36c3
md"### Continuous events"

# ╔═╡ 4f69d58a-8adc-4c0a-8143-97a41a0ff203
let
	substrate_feeding = [mm.S ~ 5] => [mm.S ~ mm.S + 5]
	
	@named mm_fed = ReactionSystem(equations(mm), continuous_events=substrate_feeding)
	
	mm_fed = complete(mm_fed)
	
	oprob = ODEProblem(mm_fed, u0map, tspan, pmap)
	
	sol = solve(oprob, Tsit5(), tstops=20)
	
	plot(sol)
end

# ╔═╡ 44fa725e-f4e5-43db-8647-4fa9f4e5e24b
md"## Steady State"

# ╔═╡ 937dbffa-c710-4f3a-80d3-c47eafabb9af
let 
	mm_continuous = @reaction_network begin
	    (kB, kD), S + E <--> ES  # reversible binding
	    kP, ES --> P + E         # conversion of substrate by enzyme
	    1, 0 --> S  # constant inflow
	    1, (S, P) --> 0  # constant outflow
	end
	
	ssprob = SteadyStateProblem(mm_continuous, u0map, pmap)
	
	sol = solve(ssprob)

	sol
end

# ╔═╡ 88c50285-6073-4532-a28d-d3373a2b632e
md"## JUMP Simulation"

# ╔═╡ 7664ad01-9426-4b81-a3ef-d9041a8826e3
let
	# starting with 500 substrate molecules and 5 enzyme molecules
	u0map_int = [:S => 500, :E => 5, :P => 0, :ES => 0]
	
	jsys = JumpInputs(mm, u0map_int, tspan, pmap)
	jprob = JumpProblem(jsys)
	
	jsol = solve(jprob)
	plot(jsol)
end

# ╔═╡ 952ee553-567c-4a70-bb24-5db4c7043395
md"## Stochastic Differential Equations"

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
# ╟─c3cccc9c-dbfc-11ef-31c5-07bc34711b2d
# ╟─df1d7710-f666-46c4-9bca-2707a997eb83
# ╟─22cc8f5e-da25-48bb-8cbc-952e205d5640
# ╟─6c756d57-075b-446d-b88b-9934a57179bd
# ╠═dc557b21-e296-446b-b7ab-0462707569c6
# ╠═f62bff24-00a4-4b7f-97ef-c2dce24b70fe
# ╠═063a77c2-5404-4729-b1c1-679c3326d636
# ╠═96703b25-43f5-4fa0-8283-c43fc5a670f2
# ╠═ea2a4e2b-a546-442c-abbe-43dc52f68c85
# ╠═10e37455-df3d-443c-ac22-529c63624b59
# ╟─7b2dbf80-68ed-47ad-b324-a4fdb003a765
# ╠═de0ac799-a461-4009-a923-df79340baf67
# ╠═cbcd8520-20ce-4585-a03c-83572d996d10
# ╠═0a6d9520-acb7-4bc1-9f90-65738a4d0b1c
# ╠═4f4c1e47-c65c-48b9-b57f-3ba0a031b95d
# ╠═c4420960-a3fb-4d0d-a095-3a0d656236f7
# ╠═0a2632e7-8942-4007-a8a4-38cedd53eb1d
# ╠═a8abeb40-cd52-4193-b484-4a56b6cda714
# ╠═a6397a2a-5a81-43ac-9dad-8c8766872eec
# ╠═a6eec47f-bf65-4fa6-ac1e-8c4f3a5d7ea4
# ╟─bb7827c6-2f44-49ba-bbb2-043f70a115ae
# ╠═2fa02aa4-d035-49e0-b19f-85fe65587914
# ╠═8ba6c0b8-1dcc-4a38-b040-a155e3c63414
# ╠═b64be898-7bd2-4266-a28b-425c45749df2
# ╠═14ae384c-4ebd-4264-8e27-5ee7c00fd8fb
# ╠═773e72a8-1041-428b-9ae4-96da03aac835
# ╠═23cd6fbb-4edb-4498-b97b-6e6cf2db8ad9
# ╠═6223afa4-c38f-46bc-b6e4-f5eff8a9efdc
# ╠═8ac9f9a3-ecbb-4f25-a7b6-88f96e256f8e
# ╠═c8a14bda-547f-40b6-be81-a6d5689a044f
# ╠═778232f3-3a5d-42f9-b698-941328866833
# ╠═183c9ba5-9dbb-42b9-9800-3b3e99ef6890
# ╠═db0390cc-8a4c-4bd5-a911-93e009c21fa8
# ╟─702d7e17-4d7f-4401-873c-c7d7ea04b450
# ╟─b5501abe-088b-44c5-96c0-9b732f691e4a
# ╠═fe8e9068-3b21-4450-a3c2-95833093a1f1
# ╟─b62081ba-aa94-4403-a300-1b08234d36c3
# ╠═4f69d58a-8adc-4c0a-8143-97a41a0ff203
# ╟─44fa725e-f4e5-43db-8647-4fa9f4e5e24b
# ╠═937dbffa-c710-4f3a-80d3-c47eafabb9af
# ╠═88c50285-6073-4532-a28d-d3373a2b632e
# ╠═7664ad01-9426-4b81-a3ef-d9041a8826e3
# ╟─952ee553-567c-4a70-bb24-5db4c7043395
# ╠═aae67739-e1ac-4949-8924-087a63881118
# ╠═73632fbb-88f7-4ff4-a7e6-b5e32780764c
# ╠═5e0acdd7-318b-4f7e-b4d3-051acc423f69
