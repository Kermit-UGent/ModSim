### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 39dfb960-d420-11ef-3bc9-d1015373ae44
using Pkg; Pkg.activate(".")

# ╔═╡ 8bfd8551-89d2-4c9f-8ebd-20e01b200e0e
using Turing, StatsPlots

# ╔═╡ f8f59613-e311-4c18-9a60-b09f04fde270
md"## Sampling"

# ╔═╡ dfa5b04a-225a-4bae-a749-74161788d9ae
md"""
Consider `X ~ Poisson(10)` and `Y ~ Poisson(X)`. By definition, $P(Y=y \mid X=x)$ is Poisson distributed. Is the unconditional PMF $P(Y = y)$ also Poisson distributed?

(Answer this question computationally by comparing the mean and the variance of $Y$, which should be equal for a Poisson distributed variable.)
"""

# ╔═╡ cc584408-e57e-4914-8008-b672f240d4ac
@model function poissonception()
	X ~ Poisson(10)
	Y ~ Poisson(X)

	return Y
end

# ╔═╡ 485eb36d-1a75-43ab-82c7-c5565160f7ff
sp_pc = [poissonception()() for i in 1:2000]

# ╔═╡ 73e59143-8f74-49d7-87a6-feebef4bccfd
histogram(sp_pc)

# ╔═╡ ede983e7-66f1-4451-b813-bede7619d251
mean(sp_pc)

# ╔═╡ 50cc7752-6ae3-4aad-ba61-f199e2a1e896
var(sp_pc) # E[Y] != var(Y) => not a poisson distribution

# ╔═╡ 4d6f284f-edd6-4de9-9d15-2ad9781d28a4
md"## Inference"

# ╔═╡ 01ba2358-b661-4fc4-9432-ef244757b0cb
md"""
There are two populations of fish living in the same pond. Consider the distributions of the lengths of both species: `X1 ~ Normal(90, 10)` and `X2 ~ Normal(70, 15)`.

You see fish of the following lengths:
`[94.0, 88.7, 89.6, 69.8, 52.8, 84.0, 89.3, 66.4, 95.1, 81.6]`

What is the ratio of each species in the pond?
"""

# ╔═╡ fa62bc2b-8551-4b0c-9e5d-6947f68c23e9
truedist = MixtureModel([Normal(90, 10), Normal(70, 15)], [0.84, 0.16])

# ╔═╡ dad5a9a5-12df-46a5-b7a5-dfd7c3816417
rand(truedist, 10) |> x -> round.(x, digits = 1) |> print

# ╔═╡ 36dad617-0ba0-432a-b94a-b8bae0fb9ed3
@model function fishmixture()
	fs1 ~ Uniform(0, 1) # fraction of species 1

	fishlens = zeros(10)
	for i in eachindex(fishlens)
		if rand() <= fs1
			fishlens[i] ~ Normal(90, 10)
		else
			fishlens[i] ~ Normal(70, 15)
		end
	end
end

# ╔═╡ 30150ff9-9d48-49cc-a803-c4e9838d8cb3
len_obs = [94.0, 88.7, 89.6, 69.8, 52.8, 84.0, 89.3, 66.4, 95.1, 81.6]

# ╔═╡ 9515dee8-71cb-4809-880a-4af9a75ce179
fishmodel = fishmixture() | (fishlens = len_obs,)

# ╔═╡ 078365e7-aeec-4139-a7b9-5ae0c8fae8ef
sample(fishmodel, MH(:fs1 => x -> Normal(x, 0.3)), 20000) |> plot

# ╔═╡ Cell order:
# ╠═39dfb960-d420-11ef-3bc9-d1015373ae44
# ╠═8bfd8551-89d2-4c9f-8ebd-20e01b200e0e
# ╟─f8f59613-e311-4c18-9a60-b09f04fde270
# ╟─dfa5b04a-225a-4bae-a749-74161788d9ae
# ╠═cc584408-e57e-4914-8008-b672f240d4ac
# ╠═485eb36d-1a75-43ab-82c7-c5565160f7ff
# ╠═73e59143-8f74-49d7-87a6-feebef4bccfd
# ╠═ede983e7-66f1-4451-b813-bede7619d251
# ╠═50cc7752-6ae3-4aad-ba61-f199e2a1e896
# ╟─4d6f284f-edd6-4de9-9d15-2ad9781d28a4
# ╟─01ba2358-b661-4fc4-9432-ef244757b0cb
# ╠═fa62bc2b-8551-4b0c-9e5d-6947f68c23e9
# ╠═dad5a9a5-12df-46a5-b7a5-dfd7c3816417
# ╠═36dad617-0ba0-432a-b94a-b8bae0fb9ed3
# ╠═30150ff9-9d48-49cc-a803-c4e9838d8cb3
# ╠═9515dee8-71cb-4809-880a-4af9a75ce179
# ╠═078365e7-aeec-4139-a7b9-5ae0c8fae8ef
