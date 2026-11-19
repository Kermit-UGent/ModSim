### A Pluto.jl notebook ###
# v0.20.4

#> [frontmatter]
#> order = "6"
#> title = "example ladybugs"
#> date = "2025-02-07"
#> tags = ["project"]
#> description = "Project example ladybugs"
#> layout = "layout.jlhtml"
#> 
#>     [[frontmatter.author]]
#>     name = "Michiel Stock"

using Markdown
using InteractiveUtils

# ╔═╡ 21357c48-f35d-11ee-23f8-2534bb1d82f4
begin

    using Pkg
    Pkg.activate("../../pluto-deployment-environment")
    
	# make this cell invisible when you are finished
	title = "APHID ANNIHILATION 🤘🔥🔥"
	names = ["Vo Orbeeld", "Pro Ject"]

	x = 4
	academic_year = "202$(x)_202$(x+1)"

	email_main_person = "mail@domain.be"

	using PlutoUI  # interactivity
	using Random # set seed
    using Catalyst, JumpProcesses, OrdinaryDiffEq # modeling
	using Optim
	using StatsPlots # plots
	TableOfContents()
end;

# ╔═╡ 14c7e803-c0ff-4211-8b60-2c9c246934dd
md"""
# $title

**$(join(names, ", ", " and "))**
"""

# ╔═╡ 6948149f-854e-4b3c-b51c-099dd221ab83
md"""
## Abstract

Belgium's native ladybugs, or lady beetles, have had it rough the last few decades. In 1995, the Asian lady beetle _Harmonia axyridis_ was introduced in Belgium as a biological pest-control agent, and quickly established itself as the new kingpin (R.L. Koch, 2003). While useful for controlling agricultural pests such as aphids, the beetle has been found to have a negative effect on biodiversity by outcompeting native ladybugs, essentially bullying the poor things (Brown et al., 2008).

The goal of this project is to create a simple model for the dynamics of a ladybug population living off aphids and subsequently to use it to investigate what feeding strategy is optimal for the ladybugs. We then intend to help the native ladybugs by passing this information on to them.

Since we are dealing with relatively small population sizes that can go to 0, we have chosen to model this system as a discrete stochastic jump process. We optimized the ladybug's aphid predation rate for an ideal deterministic case and then investigated whether how this optimal feeding rate performed in a stochastic environment.
"""

# ╔═╡ 89551690-500d-4e37-ae20-5beb71cc87ac
md"""
## Model
"""

# ╔═╡ aed58771-86c1-4928-a1ba-b7d2f503b188
md"""
First we introduce our notations for the different components of the model.

Variables:
- 🦗: The aphid population size
- 🐞: The ladybug population size

Parameters:
- 🚩🦗: The environment's carrying capacity for aphids.
- 🍼🦗: The aphid's birth rate.
- 👌🦗: The aphid's target population level as designated by the ladybugs.
- 🍴🐞: The ladybug's aphid predation rate. This is the rate at which native ladybugs eat aphids.
- 💀🐞: The ladybug's mortality rate.


We'll now introduce the model step by step in the following section. Let's start by defining the time span over which we'll simulate our model: six months, equating about one growing season
"""

# ╔═╡ eaac6193-f4b6-4901-bd7b-1ff6025c739b
tspan = (0.0, 6*30.0)

# ╔═╡ 98451eb8-8755-47fe-b8c6-ce2959256dd8
md"""
### Aphids
"""

# ╔═╡ ccfd3dbc-a9f0-4dfb-9a04-bdc69c95796b
md"""
We start with the aphids, which will serve as the ladybugs' food source. We will use a simple Gompertz model for their growth, as aphids are known for their explosive (exponential) population growth, yet their maximum population size is limited by how many aphids their host plant can feed. In reality aphids are of course able to switch host plants, but for simplicity's sake we will consider an aphid population with one single host plant (which does not deteriorate).
"""

# ╔═╡ d5b50429-b4d7-4d1a-aa24-d7aa6dabba4e
aphid_rn = @reaction_network begin
	@parameters 🚩🦗 = 10_000 🍼🦗 = 0.1
	@species 🦗(t) = 1000
	
    🍼🦗*(🚩🦗-🦗)/🚩🦗, 🦗 --> 2🦗
end

# ╔═╡ ec92e20d-bd34-4d2e-a828-281cdfff5184
sol_aph = JumpInputs(aphid_rn, [], tspan, []) |> JumpProblem |> solve;

# ╔═╡ 4273bd12-3df5-452d-a439-b2a612f8a138
plot(sol_aph, label = "Aphids", palette = :okabe_ito, title = "Evolution of an unbothered aphid population", size = (700, 400), margin = 3Plots.mm)

# ╔═╡ aa32c5e2-0f67-4910-a612-d958fa13169c
md"""
### Get eaten by ladybugs
"""

# ╔═╡ 9dcdbf69-0750-4cab-b174-0399095dc819
md"""
Next we add the ladybugs. They eat aphids to make more ladybugs, and then die of old age after living a fulfilling ladybug life.

We assume two ladybugs need to eat 30 aphids to produce a child. Additionally, the rate at which aphids get eaten scales linearly with the amount of ladybugs, but not the amount of aphids. This is because we assume the number of ladybugs will be the limiting factor: 20 ladybugs will eat twice as many aphids as 10 ladybugs, but it does not matter whether 1000 or 2000 aphids are present (we thus also assume the ladybugs have no difficulty finding the aphids). Finally, we also assume ladybugs will always leave a certain number of aphids alive lest their food source goes extinct.

To model ladybug mortality, we also assume a first-order process. The more ladybugs are present, the more ladybugs will die on any given day by reaching old age, or perhaps getting eaten by a bird.
"""

# ╔═╡ 7806ad6b-3a18-417c-98e0-a3657b9e8dd1
ladybug_rn = @reaction_network begin
	@parameters 🚩🦗 = 10_000 🍼🦗 = 0.1 👌🦗 = 300 🍴🐞 = 0.3 💀🐞 = 0.1
	@species 🦗(t) = 1000 🐞(t) = 10
	
    🍼🦗*(🚩🦗-🦗)/🚩🦗, 🦗 --> 2🦗
    🍴🐞*🐞*(🦗 - 👌🦗)/🦗, 2🐞 + 30🦗 => 3🐞
	💀🐞, 🐞 --> 0
end

# ╔═╡ d14a7e77-eb62-44dd-bc70-0853ccd41f08
md"""
Let's double check everything is in order:
"""

# ╔═╡ 8959ff62-c0fd-42b8-b03a-f58439470fb1
convert(ODESystem, ladybug_rn)

# ╔═╡ 0bd9835b-cfcc-440b-829e-666fd4899411
sol_bugged = JumpInputs(ladybug_rn, [], tspan, []) |> JumpProblem |> solve;

# ╔═╡ 5f78a353-90d2-451e-ad9d-c55b9124f3ae
# ╠═╡ show_logs = false
plot(sol_bugged, label = ["Aphids" "Ladybugs"], palette = :okabe_ito, yaxis = :log, ylims = (1, 1e4), title = "Evolution of a very much bothered aphid population (& ladybugs)", size = (800, 500), margin = 3Plots.mm)

# ╔═╡ 61667718-20fc-40b3-87b8-a7546e2713fb
md"""
And that's it for the structure of our ladybug model!
"""

# ╔═╡ f284fda5-aaf7-4121-a8f0-b996d501bec5
md"## Simulation and analysis"

# ╔═╡ 2a0b0c1f-9510-4f71-9d65-b4fb7c854a98
md"""
### Optimisation of ladybug behaviour
"""

# ╔═╡ eacb84c4-94dc-48d0-95c4-05636650057a
md"""
Given our simple ladybug model, we'd like to find the optimal ladybug behaviour so that they have the largest possible population size at the end of the growing season.

To this end we'll optimize a loss function:
- that returns the negative ladybug population size at the end of the growing season: maximize population size => minimize the inverse
- in function of the aphid predation rate $🍴🐞$ and the desired minimal aphid population level $👌🦗$: the other 3 parameters are outside of the ladybugs' control
"""

# ╔═╡ 797a59c5-d120-4008-901d-14db799ab558
md"""
For the callibration process, we'll assume a perfectly deterministic system by modeling it as an ODE problem rather than a jump problem. This is because trying to optimize a strongly stochastic process is a pain and we'd rather not deal with that.
"""

# ╔═╡ eaa2c110-9ba9-4b28-bf41-13943f54a8e6
ladyprob = ODEProblem(ladybug_rn, [], tspan, []);

# ╔═╡ ace97b1f-04ef-401a-93d8-4e6bffdb8985
function ladybug_loss(params)
	sol = remake(ladyprob, p = [:🍴🐞 => params[1], :👌🦗 => params[2]]) |> 
		x -> solve(x, reltol = 1e-9, abstol = 1e-9) # set tolerance of the solver very low or the optimizer finds the right parameters to break it and make the aphid population go below 0 for even more ladybugs

	return -sol[:🐞][end]
end

# ╔═╡ 698811dd-81b0-4a17-98f4-29109912d85f
md"""
We use constrained optimisation since both parameters must be positive. In addition, we set a minimum value of 100 for the aphids' target aphid population level. This is an estimation for how many aphids are needed to start a new population the next growing season.
"""

# ╔═╡ 88e28b17-7ed5-4da5-8310-860f1d1e8cff
res = optimize(ladybug_loss, [0, 100], [1, Inf], [0.3, 300], Fminbox(NelderMead()));

# ╔═╡ 00d095f7-e90f-4924-b6f1-9703ff26987f
params_opt = [:🍴🐞 => res.minimizer[1], :👌🦗 => res.minimizer[2]]

# ╔═╡ 3574468f-6996-4e6a-b86b-21f8ef1457a6
# ╠═╡ show_logs = false
remake(ladyprob, p = params_opt) |>
	x -> solve(x, reltol = 1e-9, abstol = 1e-9) |>
	x -> plot(x, yaxis = :log, ylims = (1, 1e4), label = ["Aphids" "Ladybug"], palette = :okabe_ito, title = "Evolution of an optimally bothered aphid population (& ladybugs)", size = (800, 500), margin = 3Plots.mm)

# ╔═╡ d85498b1-9180-4763-b79a-9fc5c326921e
md"""
We can see that a lower predation rate allows for the aphid population to flourish while the ladybug population steadily rises, until the ladybug population reaches a critical size near the end of the growing season. At this point there are so many ladybugs the aphid population starts collapsing ending with a strong crash just at the end of the growing season. Patience seems to be the key!
"""

# ╔═╡ 6619d108-791d-4cea-8690-005d1bf2b3e6
md"""
### Uncertainty assessment of ladybug behaviour 
"""

# ╔═╡ 79d65032-5b0e-412d-a69b-6aef7b1e51be
md"""
We've found a set of optimal parameters for ladybug behaviour, but that was for a deterministic case. We now want to get an idea of whether these parameters will truly result in high population sizes when considering a more realistic stochastic scenario. To achieve this, we perform a simple grid search over different values of the aphid predation rate $🍴🐞$ in the neighbourhood of the optimal value. We perform multiple runs per point to take into account the random nature of the outcome. We did not take the target aphid level $👌🦗$ into account as quickly playing around with different values showed very little change in the output (not shown for brevity).
"""

# ╔═╡ 48439a0a-cda3-4995-b69c-4d0ac8072178
function ladybug_performance(🍴🐞, jump_prob)
	sol = remake(jump_prob, p = [:🍴🐞 => 🍴🐞, :👌🦗 => res.minimizer[2]]) |>
		solve
	return sol[:🐞][end] # we're interested in the ladybug population level at the end of the growing season
end

# ╔═╡ 5f342652-3dc4-4560-a66d-1862be258e9b
begin
	Random.seed!(1337) # for reproducability
	num_repeats = 100
	jump_prob = DiscreteProblem(ladybug_rn, [], tspan, []) |>
		x -> JumpProblem(ladybug_rn, x)
	plot(title = "Final ladybug population sizes in function of aphid predation rate", xlabel = "Aphid predation rate", ylabel = "Ladybug population size", size = (800, 500), margin = 3Plots.mm)

	grid_points = 0.05:0.01:0.2
	for (grididx, 🍴🐞) in enumerate(grid_points)
		performances = [ladybug_performance(🍴🐞, jump_prob) for _ in 1:num_repeats] # do a few simulations with the same value for 🍴🐞
		boxplot!([grididx], performances, label = false, outliers = false, color = :orange)
	end
	xticks!(1:length(grid_points), string.(grid_points)) # boxplots always have a width of 1, so we need to fidget with the x-axis a little
	
end

# ╔═╡ 6605f38c-76d8-4a73-8134-c758ccf13085
md"""
From this figure we can see that ladybugs with a predation rate that is too low (<= 0.1) always go extinct by the end of the growing season, presumably because they die faster than they replenish the population. Around the optimal value we see a big peak in the expected and maximum final population sizes, steadily decreasing as the predation rate grows (with what seems another local optimum near 0.16). This confirms our previous finding that patient ladybugs will, on average, perform well with respect to their population size at the end of the growing season, though not quite as well as in a deterministic scenario.
"""

# ╔═╡ 60c63f1b-5a27-4539-a486-78c083457b0e
md"""
## Conclusion

In our project we made a simple model for the dynamics of a ladybug predating on aphids. We then optimized the predation strategy of the ladybugs in a deterministic scenario and found that a patient predation strategy can be very rewarding. We then checked whether this was also true in a stochastic scenario, and found that it was indeed the case.

Some following steps for this model would be to:
- Validate whether it can approximate real aphid-ladybug dynamics by trying to calibrate it on real data and change the model if required
- Introduce invasive ladybugs to the model and investigate their impact on the native ladybugs w.r.t. some of the model parameters
"""

# ╔═╡ 890afefc-f42b-4d74-b775-6dee5e5f0c2b
md"## Appendix"

# ╔═╡ 946cf45a-329b-42ca-9f77-a3c7bdc5611a
md"""
### References
"""

# ╔═╡ e59da868-a7f5-48f0-af1b-8a990f146f20
md"""
- Koch, R. L. (2003). The multicolored Asian lady beetle, Harmonia axyridis: a review of its biology, uses in biological control, and non-target impacts. Journal of insect Science, 3(1), 32.
- Brown, P. M. J., Adriaens, T., Bathon, H., Cuppen, J., Goldarazena, A., Hägg, T., ... & Roy, D. B. (2008). Harmonia axyridis in Europe: spread and distribution of a non-native coccinellid. From biological control to invasion: the ladybird Harmonia axyridis as a model species, 5-21.
"""

# ╔═╡ Cell order:
# ╟─14c7e803-c0ff-4211-8b60-2c9c246934dd
# ╟─6948149f-854e-4b3c-b51c-099dd221ab83
# ╟─89551690-500d-4e37-ae20-5beb71cc87ac
# ╟─aed58771-86c1-4928-a1ba-b7d2f503b188
# ╟─eaac6193-f4b6-4901-bd7b-1ff6025c739b
# ╟─98451eb8-8755-47fe-b8c6-ce2959256dd8
# ╟─ccfd3dbc-a9f0-4dfb-9a04-bdc69c95796b
# ╠═d5b50429-b4d7-4d1a-aa24-d7aa6dabba4e
# ╟─ec92e20d-bd34-4d2e-a828-281cdfff5184
# ╟─4273bd12-3df5-452d-a439-b2a612f8a138
# ╟─aa32c5e2-0f67-4910-a612-d958fa13169c
# ╟─9dcdbf69-0750-4cab-b174-0399095dc819
# ╠═7806ad6b-3a18-417c-98e0-a3657b9e8dd1
# ╟─d14a7e77-eb62-44dd-bc70-0853ccd41f08
# ╟─8959ff62-c0fd-42b8-b03a-f58439470fb1
# ╟─0bd9835b-cfcc-440b-829e-666fd4899411
# ╟─5f78a353-90d2-451e-ad9d-c55b9124f3ae
# ╟─61667718-20fc-40b3-87b8-a7546e2713fb
# ╟─f284fda5-aaf7-4121-a8f0-b996d501bec5
# ╟─2a0b0c1f-9510-4f71-9d65-b4fb7c854a98
# ╟─eacb84c4-94dc-48d0-95c4-05636650057a
# ╟─797a59c5-d120-4008-901d-14db799ab558
# ╠═eaa2c110-9ba9-4b28-bf41-13943f54a8e6
# ╠═ace97b1f-04ef-401a-93d8-4e6bffdb8985
# ╟─698811dd-81b0-4a17-98f4-29109912d85f
# ╠═88e28b17-7ed5-4da5-8310-860f1d1e8cff
# ╠═00d095f7-e90f-4924-b6f1-9703ff26987f
# ╟─3574468f-6996-4e6a-b86b-21f8ef1457a6
# ╟─d85498b1-9180-4763-b79a-9fc5c326921e
# ╟─6619d108-791d-4cea-8690-005d1bf2b3e6
# ╟─79d65032-5b0e-412d-a69b-6aef7b1e51be
# ╠═48439a0a-cda3-4995-b69c-4d0ac8072178
# ╟─5f342652-3dc4-4560-a66d-1862be258e9b
# ╟─6605f38c-76d8-4a73-8134-c758ccf13085
# ╟─60c63f1b-5a27-4539-a486-78c083457b0e
# ╟─890afefc-f42b-4d74-b775-6dee5e5f0c2b
# ╟─21357c48-f35d-11ee-23f8-2534bb1d82f4
# ╟─946cf45a-329b-42ca-9f77-a3c7bdc5611a
# ╟─e59da868-a7f5-48f0-af1b-8a990f146f20
