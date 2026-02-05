### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 06bdc430-b965-11ef-36a4-3d863afbaf6e
using Pkg; Pkg.activate("..")

# ╔═╡ ce9c8c34-3690-4241-b021-c08868157a55
using Turing, StatsPlots

# ╔═╡ ea5d107c-f171-401f-8972-6a928ab72907
using PlutoUI; TableOfContents()

# ╔═╡ eabed73e-19dc-4265-965a-cf762d630fb3
md"# Inference notebook #1: Intro"

# ╔═╡ daed8bc0-8a85-45ce-84ac-0d13ff1923f1
md"## Problem"

# ╔═╡ 6a79727c-5b5e-43cc-862f-7182c1ea878c
md"""
According to the [molecular clock hypothesis](https://en.wikipedia.org/wiki/Molecular_clock), the average amount of mutations in a gene $\bar{N}$ is proportional to how much time $t$ has passed, and identical for all species:

$\bar{N} = \alpha \, t \, .$

While this is a bit of an oversimplification, the concept has become an important tool in evolutionary biology to estimate how long ago species have diverged.
"""

# ╔═╡ 83aa4c69-c027-4ade-ac44-438d784b2b78
md"""
Consider the below figure of a small slice of the [tree of life](https://en.wikipedia.org/wiki/Tree_of_life_(biology)). Every animal represents a (fossilized) individual living during some point in evolution.
"""

# ╔═╡ 3ecd3ac7-2621-4ea6-8ff1-69962769d934
md"""
![Evolution example](https://raw.githubusercontent.com/Kermit-UGent/ModSim/2a369561ce842cf079d7660a36d0d9308739dc69/examples/ProbMod/figures/treeoflife.excalidraw.svg)
"""

# ╔═╡ 71bbd593-f2e4-40dd-b682-40e65305ebb3
md"""
We start at time 0 with a common ancestor of fish and terrestrial animals. 30 million years (Ma) later it diverges into ray-finned fish, which will give rise to most modern fish species, and lob-finned fish, which will give rise to i.e. mammals and reptiles. 

The ray-finned fish fossil is also one of the individuals for which we have DNA for its *cytochrome C* gene. The number represents that it has 25 mutations in this gene compared to the gene's sequence from our starting organism, the ancient bony fish fossil.
"""

# ╔═╡ 138538d4-f6b5-4b8e-af3a-273858cc463c
md"""
Taking into account all fossils, we can see that the number of mutations is roughly proportional with the time that has passed.
"""

# ╔═╡ f49c6d3f-5870-4473-b142-799fa84dbfb7
times = [30, 138, 375, 450]

# ╔═╡ 4ea82c1e-490f-4d88-b60c-5c2118409408
observed_mutations = [25, 94, 302, 335]

# ╔═╡ 26128162-e355-4316-8d74-291fbca194a6
scatter(times, observed_mutations, xlabel = "Time (My)", ylabel = "Number of mutations", legend = false, xlims = (0, 500), ylims = (0, 400))

# ╔═╡ d3fdfa6f-1ee2-48d0-bb6c-483a6fb4faf9
md"This leads us to our first question:"

# ╔═╡ d443d3e8-c1ee-467f-88b2-a0211d5b5cc0
md"""
!!! question
	- What are realistic values for *cytochrome C*'s mutation rate `α`?
"""

# ╔═╡ 5d83f336-372b-4493-baf6-8efffb663ff1
md"""
Figuring out the answer to the first question allows for some exciting follow-up questions! Consider for example that you find a new fossil of an ancient ancestor of the **seahorses**.
"""

# ╔═╡ 6666a388-1c45-4f6f-804f-eb3a260eae98
md"![Sharkmoment](https://raw.githubusercontent.com/Kermit-UGent/ModSim/2a369561ce842cf079d7660a36d0d9308739dc69/examples/ProbMod/figures/treeoflife2.excalidraw.svg)"

# ╔═╡ b0dbd481-ce36-4364-b995-b4fa8f36f76d
md"""
You don't know how old the fossil is, but you do find that the fossilized DNA contains **156** mutations in the _cytochrome c_ gene. How old should it be estimated as?
"""

# ╔═╡ 641a1a77-6fcb-4a6c-8f8f-3b171704efa8
md"Our second question is then:"

# ╔═╡ 442aa57f-79d8-4b22-93b2-466ac92c2b13
md"""
!!! question
	- What values are likely for the seahorse-ancestor fossil's age?
"""

# ╔═╡ 12f017e3-b7b8-408d-a677-52dd2ea900eb
md"## Explanation"

# ╔═╡ 7900fd13-cc69-41c5-972b-16d5b6d5452a
md"""### Making the model"""

# ╔═╡ 430927e8-7fb8-494c-b9da-d46f000c142f
md"""
We start again by defining a Turing model. Similar to the models of previous practical, it describes the *forward process*: how do you generate your observations (the amount of mutations $N$) based on your inputs (time that has passed $t$) and parameters (the mutation rate $α$)?

As discussed in the introduction, we assume a linear relationship between the **average** amount of mutations $\bar{N}$ and $t$:

$\bar{N} = \alpha \, t$
"""

# ╔═╡ 656d3299-e5d0-4799-84ae-da32f51bd92e
md"""
Since the accumulation of mutations is a random process, we can't expect the number of mutations $N$ to be exactly the predicted average $\bar{N}$. We can however expect it to be close to the predicted average. We can express this by saying the amount of mutations follows a distribution centered around the expected average:

$N \sim \text{Poisson}(\bar{N}) \, .$
"""

# ╔═╡ 6ec75944-0b62-43a3-80e9-d6cd2ebd4478
md"""
!!! note
	Can you explain why a Poisson distribution is a natural fit here?
"""

# ╔═╡ 0b5205c2-0eb4-4734-8470-055d29954f63
md"""
One problem: our model uses the gene's mutation rate `α`, but we don't know what it is. However, we do have some **prior** knowledge about mutation rates of essential genes (in general): they don't tend to be much larger than a few bp/My. We can encode this information by giving `α` the prior distribution `Exponential(2)`.
"""

# ╔═╡ 2034e550-84df-42ee-9f0d-2bf417563dd4
prior_alpha = Exponential(2)

# ╔═╡ 3221f614-ef9a-4e09-a7da-b649fff0ab61
plot(prior_alpha, title = "Prior belief of α", legend = false, xlabel = "α", ylabel = "Probability density")

# ╔═╡ 7a26cfec-b362-41c4-9b2e-1d040db8a5eb
md"""
We can sample some values from this prior distribution and use them to plot the *a priori* expected trends between time $t$ and average number of mutations $\bar{N}$:
"""

# ╔═╡ 1080294d-8da9-4326-a53b-afe158dc2ab9
begin
	scatter(times, observed_mutations, xlabel = "Time (My)",
		ylabel = "Number of mutations", label = false, xlims = (0, 500), 
		ylims = (0, 400), title = "A priori relationship between t and mean N"
	);
	for α in rand(prior_alpha, 1000)
		plot!(x -> α*x, color = :purple, alpha = 0.05, label = false);
	end
	plot!()
end

# ╔═╡ 624a0a67-c5c1-437c-89e1-e08e83003b41
md"""
We can see that some of the mutation rates from our prior distribution result in a relationship that matches the data well. However, our prior belief is *much too broad*: only a few of the lines are realistic for the data. Finding which of these lines from our prior belief are realistic for the data is essentially what inference does.
"""

# ╔═╡ 97257007-fbfa-4065-8788-869801ce3730
md"""
!!! note
	Why use `Exponential(2)` for the prior?

	Choosing a prior distribution is largely subjective and a big reason why some people are not fond of Bayesian modeling. There is no "one correct prior distribution". 

	However, different choices of reasonable priors often give very similar outcomes. Try running this notebook at the end with a different prior for `α`, such as `Exponential(10)` or `Uniform(0, 100)`. When are the results significantly different?
"""

# ╔═╡ 2a652734-bf47-405b-abb3-00e3d453bf47
md"""
**In summary**, we now have:
- a model that predicts our output $N$ from the input $t$:
$\bar{N} = α \, t \, .$
$N \sim \text{Poisson}(\bar{N}) \, .$
- a prior distribution for the parameter $α$.
"""

# ╔═╡ c52854b2-b2b6-4437-b65b-a1f912d054d5
md"We then translate it into a Turing model:"

# ╔═╡ 950c03a9-0056-4760-81f8-dd10d1e55cea
@model function mutations(ts)
	α ~ prior_alpha # prior distribution of parameter
	
	N = zeros(length(ts)) # output variable: in this model we have multiple values, so we need to preallocate a vector
	for i in 1:length(ts)
		N_average = α * ts[i]
		N[i] ~ Poisson(N_average)
	end
	
	return N
end

# ╔═╡ 31378eb3-51a5-4ad6-a713-7f77c7ceafcc
md"""
As you may notice, we have defined the mutation times $t$ as an input to the Turing model. This is not strictly necessary: you can also just hardcode the given values of $t$, variable `times`, in the Turing model, but this way you can easily define the model for different values of the input variables.

You simply instantiate the model with the correct values of $t$ as follows:
"""

# ╔═╡ 48f6b7dc-13aa-4057-8468-97db047773ba
mutation_model = mutations(times);

# ╔═╡ d4dbc185-6087-4862-b6e5-b5e4f51c2866
md"And can then generate random samples of the output as we are used to:"

# ╔═╡ 8b0bf05f-92a0-4ce7-8042-08c790088688
mutation_model() # random sample of N

# ╔═╡ 3e4998e7-4981-4945-8bf2-ddb5afcb43b1
chain = sample(mutation_model, Prior(), 2000);

# ╔═╡ 3421987d-1aab-4cab-bf83-dd3653715bce
α_sp = chain[:α]; # random sample of α

# ╔═╡ a9d54dfc-5337-412f-86b0-deeb5a0b6928
histogram(α_sp, title = "Sample of prior of α")

# ╔═╡ 425b4b6c-76b8-4676-8d0a-fd26711400d6
md"### Inference"

# ╔═╡ b8adbdd4-2642-4375-9979-0cb8f52c5bc8
md"""
The model so far has no extra information outside of our prior knowledge.
We can change this by **conditioning** the model on observed data as follows:
"""

# ╔═╡ 70f9e94d-a4e6-47d6-8d19-b60f7011d572
conditioned_model = mutation_model | (N = observed_mutations,)

# ╔═╡ 00c24cbb-88ba-49f2-9bf5-f44538fb2413
md"Note the syntax: we tell Turing that the value of the random variable `N` defined in our Turing model should be the values given by the variable `observed_mutations`."

# ╔═╡ 2d0c969d-03a2-4e4c-ace4-e439f81c771b
md"""
!!! danger
	Note the `,` at the end of `(N = observed_mutations,)`. This is important, as without it Julia thinks you simply put parentheses around a variable assignment and you'll get an error! See the below cell for an example.
"""

# ╔═╡ a35a43e2-e6b0-47ce-80b2-48148336274c
forgot_comma = mutation_model | (N = observed_mutations) 
	# errors because there is no `,` in the parentheses

# ╔═╡ c5f0dbb3-fba1-41f2-b7d2-740012603555
md"""
We can verify that for our conditioned model, the value of `N` has been set as constant: 
"""

# ╔═╡ d03cef36-3e82-4de4-89e7-af9f772edd8d
conditioned_model() # always returns `observed_mutations`

# ╔═╡ 371a48d5-daea-4d0b-968b-7e3056a65494
md"""
What we're after is our updated belief on the distribution of `α` given the observed data. We can do this by using the `sample` function on our model. We no longer use `Prior()` as second input, and instead choose one of the following sampling algorithms:
- `MH`: Metropolis-Hastings sampler
- `Gibbs`: Gibbs sampler
- `PG`: Particle Gibbs sampler
- `HMC`: Hamiltonian Monte Carlo sampler
- `NUTS`: No-U-Turn sampler

You can find more information about them in the corresponding Julia docs (see the `🔍Live Docs` in the bottom right corner). In practice, `NUTS` is often an excellent choice if all prior distributions are continuous and `PG` with 10-20 particles is a good default choice in all other cases. (`MH` and `Gibbs` also have their uses, but usually it takes more effort to make them work well.)
"""

# ╔═╡ 0d2c1359-434f-4f3d-8c04-c452c46d7ae8
mutation_chain = sample(conditioned_model, NUTS(), 2000)

# ╔═╡ 1c00437c-e2f3-44f6-b020-ca213e321239
md"It's always a good idea to check whether your sampling process has converged. You can do this by plotting the chain. It should look like a fuzzy caterpillar."

# ╔═╡ 4c79adff-0640-4e69-815d-ab94ebd9c937
plot(mutation_chain) # looks appropriately fuzzy!

# ╔═╡ f620d591-7982-4b67-9524-45cfac27436b
md"""
!!! note
	For an example of a non-converged chain, try using the `MH()` sampler instead of `NUTS()`. This sampling algorithm takes a lot of fiddling with its parameters (or a larger number of samples) for it to work well.
"""

# ╔═╡ 251a1b0e-7efc-4ce2-b0ea-48a5c13d2c63
md"""
The chain plot also shows the resulting **posterior distribution** of `α`. It is the prior distribution updated with the information contained in the data.
"""

# ╔═╡ 6ff29c57-aca3-4ebe-a206-e733e81bcc20
md"""
Taking the sampled values of the mutation rate from the chain and plotting a histogram will show us the exact same distribution as above. The one in the chain plot was simply smoothed to look continuous.
"""

# ╔═╡ b657217a-6ccc-4a41-b852-df4e39a7a10a
sp_alpha = mutation_chain[:α];

# ╔═╡ 9cf08616-d599-42a3-82e8-e8c98853c1d8
histogram(sp_alpha) # note the difference with the prior distribution!

# ╔═╡ e02c42dc-627c-4bc0-8ebb-fdb2b0f15b64
md"Plotting some sampled mutation rates from this distribution onto our data shows that they fit well:"

# ╔═╡ 87e70d5a-7a45-4a3e-b6c4-a894cc78621b
begin
	scatter(times, observed_mutations, xlabel = "Time (My)",
		ylabel = "Number of mutations", label = false, xlims = (0, 500), 
		ylims = (0, 400), title = "A posteriori relationship between t and mean N"
	);
	plot!([x -> αᵢ*x for αᵢ in sp_alpha[1:10:end]], color = :purple, opacity = 0.1, label = false)
end

# ╔═╡ 87059440-2919-4f6d-9d32-0df3ce75e2a2
md"And finally we can answer our first question:"

# ╔═╡ 644cff58-68b9-4d4b-8896-617fcacc39c5
mean(sp_alpha)

# ╔═╡ 5f2f5a78-ca90-4baf-8bf0-ff1fae88d785
sqrt(var(sp_alpha))

# ╔═╡ bfbf811f-7b2f-473b-b132-0c7e965c8b0d
md"""
 $α$ is ± normally distributed around 0.75 with a standard deviation of 0.025.
"""

# ╔═╡ d2d83d9a-0ed6-431b-8764-397c1bb019c2
md"### Seahorses (extra)"

# ╔═╡ 4e9dc370-6aca-40f7-807b-e85d412ab1a0
md"""
To answer how old the ancestral seahorse fossil is, we need to update the model a little.
So far the fossil ages were considered to be known exactly and given as input to the model `ts`. Since the seahorse fossil's age is unknown, we add a parameter for it called `fossil_age`. 

As prior knowledge we can use the fact that it must have evolved _after_ the ray-finned fish fossil (30 Ma after weird old fish), but _before_ modern seahorses (450 Ma after the bony fish fossil).
"""

# ╔═╡ fd15afe1-72d7-4663-b2b6-afa0dd219db8
@model function horsetations(ts)
	α ~ prior_alpha # prior distribution of parameter
	
	N = zeros(length(ts)) # output variable: in this model we have multiple values, so we need to preallocate a vector
	for i in eachindex(ts)
		N_average = α * ts[i]
		N[i] ~ Poisson(N_average)
	end
	
	fossil_age ~ Uniform(30, 450)
	horse_mutations ~ Poisson(α * fossil_age)
	
	return N
end

# ╔═╡ 43509549-f926-478b-a4da-995de443b3a7
md"Then we simply repeat model instantation, conditioning and sampling:"

# ╔═╡ dbb8ad46-4ac1-443e-b39a-89ddf938ede2
horse_model = horsetations(times)

# ╔═╡ 8b380751-0e55-4691-bab0-247fa1b7c510
horseditioned_model = horse_model | (N = observed_mutations, horse_mutations = 156);

# ╔═╡ e4eda27b-cff2-4ef9-bbdc-75ce2b19b10b
horse_chain = sample(horseditioned_model, NUTS(), 2000)

# ╔═╡ da82b71b-fdf5-4fba-8ef3-f2aac67d3494
md"And we have our posterior distribution of `fossil_age`! It seems like the seahorse ancestor lived about 200-220 million years after the bony fish fossil, or about 240 million years ago."

# ╔═╡ 7e5a39c7-cb66-4563-8f58-7245f88c9b85
histogram(horse_chain[:fossil_age])

# ╔═╡ 1e152790-a277-4dba-8785-d6a6120cc34f
md"## The essentials"

# ╔═╡ 0e044dc4-2fc4-4c41-9efa-e586fd138a69
md"This section again reitarates the essential code for this practical without much explanations."

# ╔═╡ e9d2f2ca-48ed-4cd8-bd0a-ff0ef640c8fb
let
	@model function mutations(ts)
		α ~ Exponential(2) # prior distribution of parameter
		
		N = zeros(length(ts)) # output variable: in this model we have multiple values, so we need to preallocate a vector
		for i in eachindex(ts)
			N_average = α * ts[i]
			N[i] ~ Poisson(N_average)
		end
		
		return N
	end

	mutation_model = mutations(times); # instantiate model
	conditioned_model = mutation_model | (N = observed_mutations,) 
		# condition model
		# mind the `,` after `observed_mutations`!
	mutation_chain = sample(conditioned_model, NUTS(), 2000);
		# run inference with an appropriate sampler
	α_sp = mutation_chain[:α] # get sample of posterior distribution of α
	histogram(α_sp, 
		title = "Posterior distribution of mutation rate α"
	)
end

# ╔═╡ Cell order:
# ╟─eabed73e-19dc-4265-965a-cf762d630fb3
# ╠═06bdc430-b965-11ef-36a4-3d863afbaf6e
# ╠═ce9c8c34-3690-4241-b021-c08868157a55
# ╠═ea5d107c-f171-401f-8972-6a928ab72907
# ╟─daed8bc0-8a85-45ce-84ac-0d13ff1923f1
# ╟─6a79727c-5b5e-43cc-862f-7182c1ea878c
# ╟─83aa4c69-c027-4ade-ac44-438d784b2b78
# ╟─3ecd3ac7-2621-4ea6-8ff1-69962769d934
# ╟─71bbd593-f2e4-40dd-b682-40e65305ebb3
# ╟─138538d4-f6b5-4b8e-af3a-273858cc463c
# ╠═f49c6d3f-5870-4473-b142-799fa84dbfb7
# ╠═4ea82c1e-490f-4d88-b60c-5c2118409408
# ╟─26128162-e355-4316-8d74-291fbca194a6
# ╟─d3fdfa6f-1ee2-48d0-bb6c-483a6fb4faf9
# ╟─d443d3e8-c1ee-467f-88b2-a0211d5b5cc0
# ╟─5d83f336-372b-4493-baf6-8efffb663ff1
# ╟─6666a388-1c45-4f6f-804f-eb3a260eae98
# ╟─b0dbd481-ce36-4364-b995-b4fa8f36f76d
# ╟─641a1a77-6fcb-4a6c-8f8f-3b171704efa8
# ╟─442aa57f-79d8-4b22-93b2-466ac92c2b13
# ╟─12f017e3-b7b8-408d-a677-52dd2ea900eb
# ╟─7900fd13-cc69-41c5-972b-16d5b6d5452a
# ╟─430927e8-7fb8-494c-b9da-d46f000c142f
# ╟─656d3299-e5d0-4799-84ae-da32f51bd92e
# ╟─6ec75944-0b62-43a3-80e9-d6cd2ebd4478
# ╟─0b5205c2-0eb4-4734-8470-055d29954f63
# ╠═2034e550-84df-42ee-9f0d-2bf417563dd4
# ╟─3221f614-ef9a-4e09-a7da-b649fff0ab61
# ╟─7a26cfec-b362-41c4-9b2e-1d040db8a5eb
# ╟─1080294d-8da9-4326-a53b-afe158dc2ab9
# ╟─624a0a67-c5c1-437c-89e1-e08e83003b41
# ╟─97257007-fbfa-4065-8788-869801ce3730
# ╟─2a652734-bf47-405b-abb3-00e3d453bf47
# ╟─c52854b2-b2b6-4437-b65b-a1f912d054d5
# ╠═950c03a9-0056-4760-81f8-dd10d1e55cea
# ╟─31378eb3-51a5-4ad6-a713-7f77c7ceafcc
# ╠═48f6b7dc-13aa-4057-8468-97db047773ba
# ╟─d4dbc185-6087-4862-b6e5-b5e4f51c2866
# ╠═8b0bf05f-92a0-4ce7-8042-08c790088688
# ╠═3e4998e7-4981-4945-8bf2-ddb5afcb43b1
# ╠═3421987d-1aab-4cab-bf83-dd3653715bce
# ╟─a9d54dfc-5337-412f-86b0-deeb5a0b6928
# ╟─425b4b6c-76b8-4676-8d0a-fd26711400d6
# ╟─b8adbdd4-2642-4375-9979-0cb8f52c5bc8
# ╠═70f9e94d-a4e6-47d6-8d19-b60f7011d572
# ╟─00c24cbb-88ba-49f2-9bf5-f44538fb2413
# ╟─2d0c969d-03a2-4e4c-ace4-e439f81c771b
# ╠═a35a43e2-e6b0-47ce-80b2-48148336274c
# ╟─c5f0dbb3-fba1-41f2-b7d2-740012603555
# ╠═d03cef36-3e82-4de4-89e7-af9f772edd8d
# ╟─371a48d5-daea-4d0b-968b-7e3056a65494
# ╠═0d2c1359-434f-4f3d-8c04-c452c46d7ae8
# ╟─1c00437c-e2f3-44f6-b020-ca213e321239
# ╠═4c79adff-0640-4e69-815d-ab94ebd9c937
# ╟─f620d591-7982-4b67-9524-45cfac27436b
# ╟─251a1b0e-7efc-4ce2-b0ea-48a5c13d2c63
# ╟─6ff29c57-aca3-4ebe-a206-e733e81bcc20
# ╠═b657217a-6ccc-4a41-b852-df4e39a7a10a
# ╠═9cf08616-d599-42a3-82e8-e8c98853c1d8
# ╟─e02c42dc-627c-4bc0-8ebb-fdb2b0f15b64
# ╟─87e70d5a-7a45-4a3e-b6c4-a894cc78621b
# ╟─87059440-2919-4f6d-9d32-0df3ce75e2a2
# ╠═644cff58-68b9-4d4b-8896-617fcacc39c5
# ╠═5f2f5a78-ca90-4baf-8bf0-ff1fae88d785
# ╟─bfbf811f-7b2f-473b-b132-0c7e965c8b0d
# ╟─d2d83d9a-0ed6-431b-8764-397c1bb019c2
# ╟─4e9dc370-6aca-40f7-807b-e85d412ab1a0
# ╠═fd15afe1-72d7-4663-b2b6-afa0dd219db8
# ╟─43509549-f926-478b-a4da-995de443b3a7
# ╠═dbb8ad46-4ac1-443e-b39a-89ddf938ede2
# ╠═8b380751-0e55-4691-bab0-247fa1b7c510
# ╠═e4eda27b-cff2-4ef9-bbdc-75ce2b19b10b
# ╟─da82b71b-fdf5-4fba-8ef3-f2aac67d3494
# ╠═7e5a39c7-cb66-4563-8f58-7245f88c9b85
# ╟─1e152790-a277-4dba-8785-d6a6120cc34f
# ╟─0e044dc4-2fc4-4c41-9efa-e586fd138a69
# ╠═e9d2f2ca-48ed-4cd8-bd0a-ff0ef640c8fb
