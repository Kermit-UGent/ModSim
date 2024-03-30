### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 27037ea0-cd76-11ee-3d7f-07c43460afea

begin
    using Pkg
	Pkg.activate("..")
    using Turing, Plots
end

# ╔═╡ dd5d4013-d775-4942-9ede-ab7674b8e3a1
md"""

Give the histogram of $X\sim$ binom(100, p) with a prior distribution for $p$ as
1.  $\delta(0.7)$ ($p$ is always 0.7),
2.  $p\sim$ Uniform(0.6, 0.8)
3.  $p\sim Normal(0.7, 0.1)$ (use TruncatedNormal to generate valid values for $p$)
4.  $p\sim$ Beta(7, 3)
5.  $p$ equal to 0.6 in 50% of the time and 0.8 elsewise

"""

# ╔═╡ a2ca08ee-9b33-45b2-ab3c-bb273dfd8158
@model function psample(pdist; n=50)
	p ~ pdist
	X ~ Binomial(n, p)
end

# ╔═╡ b06dd24d-e146-4803-8e70-f800219bf923
delta = Dirac(0.8)

# ╔═╡ b134c532-4439-45f3-ae50-217e101ca2c5
unif = Uniform(0.5, 0.9)

# ╔═╡ 5f72e7ff-6b04-4a98-8996-3c70accd82b6
plot(p->pdf(unif, p), 0, 1, xlab="p", ylab="f_p(p)", label="uniform")

# ╔═╡ daa8f129-d7bd-49f7-8adc-c21277cd45e7
norm = TruncatedNormal(0.7, 0.1, 0, 1)

# ╔═╡ 75d72d23-2eab-4ee1-b4d2-6d26c8f39611
plot(p->pdf(norm, p), 0, 1, xlab="p", ylab="f_p(p)", label="normal")

# ╔═╡ 2dca96c0-599d-48f7-a5f5-f45d93eb7103
beta = Beta(7, 3)

# ╔═╡ c0dbf5a3-edf3-49d5-af23-150f4f4e15e0
plot(p->pdf(beta, p), 0, 1, xlab="p", ylab="f_p(p)", label="beta")

# ╔═╡ 8846935e-ceea-46da-9f81-71449f8f4a68
mixture = MixtureModel([Dirac(0.6), Dirac(0.8)], [0.5, 0.5])

# ╔═╡ d62d94f2-7210-4039-97fb-fe0657166439
md"The means are all the same:"

# ╔═╡ ac87933d-d488-4fe9-ab97-bcad13223620
mean(delta)

# ╔═╡ 113d2b87-6a68-427d-995b-62c1b8f47e32
mean(unif)

# ╔═╡ 7dbfdd33-12c5-40b1-b77c-6a3e907dca32
mean(norm)

# ╔═╡ 79e068b4-2f09-4c58-a12b-f281beb0b310
mean(beta)

# ╔═╡ 64e62e77-7f0f-40a1-a253-16b95f1ced55
mean(mixture)

# ╔═╡ fceb1451-de7f-40da-be0b-dbf460d95b26
md"Let us make a function for the plots"

# ╔═╡ 1e5677be-cab4-4557-a828-d7e2f2833b9d
function plot_histogram(distribution; n_samples=10_000)
	sample = [rand(psample(distribution))[:X] for i in 1:10_000]
	return histogram(sample, label=string(distribution), xlab="x", legend = :outertop)
end

# ╔═╡ e2ca2876-a4dd-4171-9458-7ac48e7e42a3
plot_histogram(delta)

# ╔═╡ f76da90c-4285-41ea-b508-450aac6e85ec
plot_histogram(unif)

# ╔═╡ 0c42f962-2961-48dd-84bc-c1d9d1c53d96
plot_histogram(norm)

# ╔═╡ b3ec0fb1-46c4-4a8b-bac6-a2620daa2457
plot_histogram(beta)

# ╔═╡ ace3fed3-8362-486b-b811-647627c1cb02
plot_histogram(mixture)

# ╔═╡ Cell order:
# ╠═27037ea0-cd76-11ee-3d7f-07c43460afea
# ╟─dd5d4013-d775-4942-9ede-ab7674b8e3a1
# ╠═a2ca08ee-9b33-45b2-ab3c-bb273dfd8158
# ╠═b06dd24d-e146-4803-8e70-f800219bf923
# ╠═b134c532-4439-45f3-ae50-217e101ca2c5
# ╟─5f72e7ff-6b04-4a98-8996-3c70accd82b6
# ╠═daa8f129-d7bd-49f7-8adc-c21277cd45e7
# ╟─75d72d23-2eab-4ee1-b4d2-6d26c8f39611
# ╠═2dca96c0-599d-48f7-a5f5-f45d93eb7103
# ╟─c0dbf5a3-edf3-49d5-af23-150f4f4e15e0
# ╠═8846935e-ceea-46da-9f81-71449f8f4a68
# ╟─d62d94f2-7210-4039-97fb-fe0657166439
# ╠═ac87933d-d488-4fe9-ab97-bcad13223620
# ╠═113d2b87-6a68-427d-995b-62c1b8f47e32
# ╠═7dbfdd33-12c5-40b1-b77c-6a3e907dca32
# ╠═79e068b4-2f09-4c58-a12b-f281beb0b310
# ╠═64e62e77-7f0f-40a1-a253-16b95f1ced55
# ╠═fceb1451-de7f-40da-be0b-dbf460d95b26
# ╠═1e5677be-cab4-4557-a828-d7e2f2833b9d
# ╠═e2ca2876-a4dd-4171-9458-7ac48e7e42a3
# ╠═f76da90c-4285-41ea-b508-450aac6e85ec
# ╠═0c42f962-2961-48dd-84bc-c1d9d1c53d96
# ╠═b3ec0fb1-46c4-4a8b-bac6-a2620daa2457
# ╠═ace3fed3-8362-486b-b811-647627c1cb02
