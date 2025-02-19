### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 06bdc430-b965-11ef-36a4-3d863afbaf6e
using Pkg; Pkg.activate("..")

# ╔═╡ ce9c8c34-3690-4241-b021-c08868157a55
using Turing, StatsPlots

# ╔═╡ eabed73e-19dc-4265-965a-cf762d630fb3
md"# Inference notebook #1: Intro"

# ╔═╡ daed8bc0-8a85-45ce-84ac-0d13ff1923f1
md"## Problem"

# ╔═╡ 6a79727c-5b5e-43cc-862f-7182c1ea878c
md"""
According to the [molecular clock hypothesis](https://en.wikipedia.org/wiki/Molecular_clock), the amount of mutations in a gene is proportional to how much time has passed, and identical for all species. While this is a bit of an oversimplification, the concept has become an important tool in evolutionary biology to estimate how long ago species have diverged.
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
scatter(times, observed_mutations, xlabel = "Time (My)", ylabel = "Number of mutations", legend = false, xlims = (0, 500))

# ╔═╡ 5d83f336-372b-4493-baf6-8efffb663ff1
md"""
Consider now that you find a new fossil of an ancient ancestor of the **seahorses**.
"""

# ╔═╡ 6666a388-1c45-4f6f-804f-eb3a260eae98
md"![Sharkmoment](https://raw.githubusercontent.com/Kermit-UGent/ModSim/2a369561ce842cf079d7660a36d0d9308739dc69/examples/ProbMod/figures/treeoflife2.excalidraw.svg)"

# ╔═╡ b0dbd481-ce36-4364-b995-b4fa8f36f76d
md"""
You don't know how old the fossil is, but you do find that the fossilized DNA contains **156** mutations in the _cytochrome c_ gene. How old should it be estimated as?
"""

# ╔═╡ 442aa57f-79d8-4b22-93b2-466ac92c2b13
md"""
!!! questions
	- What is *cytochrome C*'s mutation rate `α`?
	- What is the seahorse-ancestor fossil's age?
"""

# ╔═╡ 1e152790-a277-4dba-8785-d6a6120cc34f
md"## Copy-paste example"

# ╔═╡ 0e044dc4-2fc4-4c41-9efa-e586fd138a69
md"This section contains the essential code for this practical. A detailed explanation is given in the next section."

# ╔═╡ e9d2f2ca-48ed-4cd8-bd0a-ff0ef640c8fb
let
	@model function mutations(ts)
		α ~ Exponential(10)
		
		num_mutations = zeros(length(ts))
		for i in eachindex(ts)
			num_mutations[i] ~ Poisson(α * ts[i])
		end
		
		return num_mutations
	end

	mutation_model = mutations(times);
	conditioned_model = mutation_model | (num_mutations = observed_mutations,)
		# mind the `,` after `observed_mutations`!
	mutation_chain = sample(conditioned_model, NUTS(), 2000);
	α_sp = mutation_chain[:α]
	histogram(α_sp, 
		title = "Posterior distribution of mutation rate α"
	)
end

# ╔═╡ 12f017e3-b7b8-408d-a677-52dd2ea900eb
md"## Explanation"

# ╔═╡ 7900fd13-cc69-41c5-972b-16d5b6d5452a
md"""### Making the model"""

# ╔═╡ 430927e8-7fb8-494c-b9da-d46f000c142f
md"""
We start again by defining a Turing model. Similar to the models of previous practical, it describes the *forward process*: how do you generate your observations (the amount of mutations) based on your inputs (age of fossil) and parameters (the mutation rate)?

This may seem unintuitive, as we don't know the distribution of this gene's mutation rate `α`. However, we do have some **prior** knowledge about mutation rates of genes (in general): they don't tend to be much larger than a few bp/My. We can encode this information by giving `α` the prior distribution `Exponential(10)`.
"""

# ╔═╡ 2034e550-84df-42ee-9f0d-2bf417563dd4
prior_alpha = Exponential(10)

# ╔═╡ 3221f614-ef9a-4e09-a7da-b649fff0ab61
plot(prior_alpha, title = "Prior belief of α", legend = false, xlabel = "α", ylabel = "Probability density")

# ╔═╡ 97257007-fbfa-4065-8788-869801ce3730
md"""
!!! note
	Why use `Exponential(10)` for the prior and not `Exponential(1)`, or some other value?

	Choosing a prior distribution is largely subjective and a big reason why some people are not fond of Bayesian modeling. There is no "one correct prior distribution". 

	However, different choices of reasonable priors often give very similar outcomes. Try running this notebook at the end with a different prior for `α`, such as `Exponential(1)` or `Uniform(0, 100)`. When are the results significantly different?
"""

# ╔═╡ 2a652734-bf47-405b-abb3-00e3d453bf47
md"""
The rest of the model is pretty straightforward: if the mutation rate `α` is constant, the number of mutations after `t` million years should be about `α*t`. 

Since the accumulation of mutations is a random process, we can't expect the number of mutations to be exactly this number. Rather, we define it to follow a probability distribution centered around this number. A Poisson distribution is chosen as a good fit for count data.
"""

# ╔═╡ 950c03a9-0056-4760-81f8-dd10d1e55cea
@model function mutations(ts)
	α ~ prior_alpha
	
	num_mutations = zeros(length(ts))
	for i in eachindex(ts)
		num_mutations[i] ~ Poisson(α * ts[i])
	end
	
	return num_mutations
end

# ╔═╡ 31378eb3-51a5-4ad6-a713-7f77c7ceafcc
md"""
The model is instantiated with the correct inputs and can be used to generate samples as per usual.
"""

# ╔═╡ 48f6b7dc-13aa-4057-8468-97db047773ba
mutation_model = mutations(times);

# ╔═╡ 8b0bf05f-92a0-4ce7-8042-08c790088688
mutation_model() # random sample of num_mutations

# ╔═╡ 3e4998e7-4981-4945-8bf2-ddb5afcb43b1
chain = sample(mutation_model, Prior(), 2000);

# ╔═╡ a9d54dfc-5337-412f-86b0-deeb5a0b6928
histogram(chain[:α])

# ╔═╡ 425b4b6c-76b8-4676-8d0a-fd26711400d6
md"### Inference"

# ╔═╡ b8adbdd4-2642-4375-9979-0cb8f52c5bc8
md"""
The model so far has no extra information outside of our prior knowledge.
We can change this by **conditioning** the model on observed data as follows:
"""

# ╔═╡ 70f9e94d-a4e6-47d6-8d19-b60f7011d572
conditioned_model = mutation_model | (num_mutations = observed_mutations,)

# ╔═╡ 2d0c969d-03a2-4e4c-ace4-e439f81c771b
md"""
!!! danger
	Note the `,` at the end of `(num_mutations = observed_mutations,)`. This is important, as without it Julia thinks you simply put parentheses around a variable assignment and you'll get an error! See the below cell for an example.
"""

# ╔═╡ a35a43e2-e6b0-47ce-80b2-48148336274c
forgot_comma = mutation_model | (num_mutations = observed_mutations) 
	# errors because there is no `,` in the parentheses

# ╔═╡ c5f0dbb3-fba1-41f2-b7d2-740012603555
md"""
We can verify that for our conditioned model, the values of `num_mutations` has been set as constant: 
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

You can find more information about them in the corresponding Julia docs. In practice, `NUTS` is often an excellent choice if all variables are continuous and `PG` is a good default choice in all other cases. (`MH` and `Gibbs` also have their uses, but usually it takes more effort to make them work well.)
"""

# ╔═╡ 0d2c1359-434f-4f3d-8c04-c452c46d7ae8
mutation_chain = sample(conditioned_model, NUTS(), 2000)

# ╔═╡ 7441c82c-8aad-4255-92fe-14cd8ee93262
begin
	scatter(times, observed_mutations, xlabel = "Time (My)",
		ylabel = "Number of mutations", label = false, xlims = (0, 500), 
		title = "Predicted trend"
	);
	for α in mutation_chain[:α][1:10:end]
		plot!(x -> α*x, color = :purple, alpha = 0.05, label = false);
	end
	plot!()
end

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
Taking the sampled values of the mutation rate from the chain and plotting a histogram will show us the exact same distribution. The one in the chain plot was simply smoothed to look continuous.
"""

# ╔═╡ b657217a-6ccc-4a41-b852-df4e39a7a10a
sp_alpha = mutation_chain[:α];

# ╔═╡ 9cf08616-d599-42a3-82e8-e8c98853c1d8
histogram(sp_alpha)

# ╔═╡ e02c42dc-627c-4bc0-8ebb-fdb2b0f15b64
md"Plotting some sampled mutation rates from this distribution onto our data shows that they fit well:"

# ╔═╡ 87e70d5a-7a45-4a3e-b6c4-a894cc78621b
begin
	scatter(times, observed_mutations, xlabel = "Time (My)", ylabel = "Number of mutations", label = false, xlims = (0, 500))
	plot!([x -> αᵢ*x for αᵢ in sp_alpha[1:10:end]], color = :blue, opacity = 0.1, label = false)
end

# ╔═╡ 644cff58-68b9-4d4b-8896-617fcacc39c5
mean(sp_alpha)

# ╔═╡ 5f2f5a78-ca90-4baf-8bf0-ff1fae88d785
sqrt(var(sp_alpha))

# ╔═╡ bfbf811f-7b2f-473b-b132-0c7e965c8b0d
md"To answer our first question, α is ± normally distributed around 0.75 with a standard deviation of 0.025."

# ╔═╡ d2d83d9a-0ed6-431b-8764-397c1bb019c2
md"### Seahorses (extra)"

# ╔═╡ 4e9dc370-6aca-40f7-807b-e85d412ab1a0
md"""
To answer how old the ancestral seahorse fossil is, we need to update the model a little.
So far the fossil ages were considered to be known exactly and given as input to the model `ts`. Since the fossil's age is unknown, we add a parameter `fossil_age`. 

As prior knowledge we can use the fact that it must have evolved _after_ the ray-finned fish fossil (30 Ma after weird old fish), but _before_ modern seahorses (450 Ma after the bony fish fossil).
"""

# ╔═╡ fd15afe1-72d7-4663-b2b6-afa0dd219db8
@model function horsetations(ts)
	α ~ Exponential(10)

	num_mutations = zeros(length(ts))
	for i in eachindex(ts)
		num_mutations[i] ~ Poisson(α * ts[i])
	end
	
	fossil_age ~ Uniform(30, 450)
	horse_mutations ~ Poisson(α * fossil_age)
	
	return num_mutations
end

# ╔═╡ 43509549-f926-478b-a4da-995de443b3a7
md"Then we simply repeat model instantation, conditioning and sampling:"

# ╔═╡ dbb8ad46-4ac1-443e-b39a-89ddf938ede2
horse_model = horsetations(times)

# ╔═╡ 8b380751-0e55-4691-bab0-247fa1b7c510
horseditioned_model = horse_model | (num_mutations = observed_mutations, horse_mutations = 156);

# ╔═╡ e4eda27b-cff2-4ef9-bbdc-75ce2b19b10b
horse_chain = sample(horseditioned_model, NUTS(), 2000)

# ╔═╡ da82b71b-fdf5-4fba-8ef3-f2aac67d3494
md"And we have our posterior distribution of `fossil_age`! It seems like the seahorse ancestor lived about 200-220 million years after the bony fish fossil, or about 240 million years ago."

# ╔═╡ 7e5a39c7-cb66-4563-8f58-7245f88c9b85
histogram(horse_chain[:fossil_age])

# ╔═╡ Cell order:
# ╟─eabed73e-19dc-4265-965a-cf762d630fb3
# ╠═06bdc430-b965-11ef-36a4-3d863afbaf6e
# ╠═ce9c8c34-3690-4241-b021-c08868157a55
# ╟─daed8bc0-8a85-45ce-84ac-0d13ff1923f1
# ╟─6a79727c-5b5e-43cc-862f-7182c1ea878c
# ╟─83aa4c69-c027-4ade-ac44-438d784b2b78
# ╟─3ecd3ac7-2621-4ea6-8ff1-69962769d934
# ╟─71bbd593-f2e4-40dd-b682-40e65305ebb3
# ╟─138538d4-f6b5-4b8e-af3a-273858cc463c
# ╠═f49c6d3f-5870-4473-b142-799fa84dbfb7
# ╠═4ea82c1e-490f-4d88-b60c-5c2118409408
# ╟─26128162-e355-4316-8d74-291fbca194a6
# ╟─5d83f336-372b-4493-baf6-8efffb663ff1
# ╟─6666a388-1c45-4f6f-804f-eb3a260eae98
# ╟─b0dbd481-ce36-4364-b995-b4fa8f36f76d
# ╟─442aa57f-79d8-4b22-93b2-466ac92c2b13
# ╟─1e152790-a277-4dba-8785-d6a6120cc34f
# ╟─0e044dc4-2fc4-4c41-9efa-e586fd138a69
# ╠═e9d2f2ca-48ed-4cd8-bd0a-ff0ef640c8fb
# ╠═7441c82c-8aad-4255-92fe-14cd8ee93262
# ╟─12f017e3-b7b8-408d-a677-52dd2ea900eb
# ╟─7900fd13-cc69-41c5-972b-16d5b6d5452a
# ╟─430927e8-7fb8-494c-b9da-d46f000c142f
# ╠═2034e550-84df-42ee-9f0d-2bf417563dd4
# ╟─3221f614-ef9a-4e09-a7da-b649fff0ab61
# ╟─97257007-fbfa-4065-8788-869801ce3730
# ╟─2a652734-bf47-405b-abb3-00e3d453bf47
# ╠═950c03a9-0056-4760-81f8-dd10d1e55cea
# ╟─31378eb3-51a5-4ad6-a713-7f77c7ceafcc
# ╠═48f6b7dc-13aa-4057-8468-97db047773ba
# ╠═8b0bf05f-92a0-4ce7-8042-08c790088688
# ╠═3e4998e7-4981-4945-8bf2-ddb5afcb43b1
# ╠═a9d54dfc-5337-412f-86b0-deeb5a0b6928
# ╟─425b4b6c-76b8-4676-8d0a-fd26711400d6
# ╟─b8adbdd4-2642-4375-9979-0cb8f52c5bc8
# ╠═70f9e94d-a4e6-47d6-8d19-b60f7011d572
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
