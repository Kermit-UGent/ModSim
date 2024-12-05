### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ a2410616-5a17-403c-aa2f-dc93c2633c7f
using Pkg; Pkg.activate("../..");

# ╔═╡ 8307092d-368d-441d-8315-3dc312026534
using Turing, StatsPlots

# ╔═╡ aeb0aef0-b2ee-11ef-3cca-7f80b487ea17
md"# Sampling notebook #1: Intro"

# ╔═╡ 30957e05-85b8-4106-9635-82a5c11d9825
md"""
This notebook will guide you through the basics of sampling in Julia.
"""

# ╔═╡ 13636a8a-8c31-4397-8a7d-5cd7f899d7a5
md"""
To start off, load the required packages.
"""

# ╔═╡ 9f89e350-199f-4875-947b-61df653ffc19
md"## Problem"

# ╔═╡ 2e544855-16b8-4794-9ba7-70a1e7209dd2
md"""
To explain, let's go back to the circle throw example from the theory. 

The idea is simple: a circle with radius 1 has an area of π. If you throw darts at the unit square [-1, 1] x [-1, 1] with uniform probability, the probability of a dart landing inside that circle is the area of the circle over the area of the square.

```math
\begin{align}

A_{circle} &= π
\\ A_{square} &= (1 - (-1))^2 = 4
\\ P_{inside} &= \frac{A_{circle}}{A_{square}} = \frac{π}{4}

\end{align}
```
"""

# ╔═╡ 88b3087c-f78a-4df4-a9dc-5696bf4052d5
md"""
This means if we can estimate this probability, we can estimate π ezpz!

```math
\begin{align}

P_{inside} &= \frac{π}{4}
\\ π &= 4 \, P_{inside}

\end{align}
```
"""

# ╔═╡ 78f872c6-013b-4ad1-966f-9e0ce3288019
md"""
However simple the problem, this probability is not simple to calculate by hand:

```math
P_{inside} = \int_{0}^{1} \, \int_0^{\sqrt{1 - y^2}} dx \, dy
```
"""

# ╔═╡ bde6033f-15a8-4716-a41e-080f5d48e9d6
md"""
This, as one may guess, only gets worse for more complex problems. Which is why we use sampling instead!
"""

# ╔═╡ 1724e4e4-18d2-430d-a8b0-11f5891b09a3
md"""
!!! question
	Are you a real mathhead? Try computing the integral by hand! 
	
	Hint: there should be an inverse tangent function somewhere down the line.
"""

# ╔═╡ 5acc6791-7cb6-4a31-a28d-649e761329ee
md"## Copy-paste example"

# ╔═╡ f3a73486-485a-4c69-b3db-2153b5a06bd8
md"This section showcases the most essential code for the first practical. The next sections explain it in detail. 

The variable names are prefaced with `_` since Pluto doesn't let you use the same name twice."

# ╔═╡ 6b554b3e-3fc8-4a20-9aac-a2b351dd7f9e
_n_samples = 1000

# ╔═╡ 34925f20-0fc3-4bf3-8a8d-89a9596504c4
@model function _distances()
	x ~ Uniform(-1, 1)
	y ~ Uniform(-1, 1)
	dist = sqrt(x^2 + y^2)
	return dist
end

# ╔═╡ 0712f06a-039a-4735-8501-86d3fe8d9114
_dist_model = _distances();

# ╔═╡ c6f18791-892d-4915-b4e4-f540fb770104
_sp_dists = [_dist_model() for sample_idx in 1:_n_samples];

# ╔═╡ 5c7cf0d6-2d6b-487e-9609-546464f730a1
_sp_inside = _sp_dists .<= 1;

# ╔═╡ c52c7639-c400-4093-a1e2-a0473059766c
4 * mean(_sp_inside)

# ╔═╡ ccdcfa65-a02f-4110-a664-36090c3291d8
md"## Explanation"

# ╔═╡ 43e8f794-fec2-4ce4-9306-2fc8a9565343
md"### Defining the model"

# ╔═╡ ed788556-1627-4e9a-b901-e532272a8265
md"""
Turing models are defined as julia functions preceded by the `@model` macro.
Inside of them, you can define random variables with the "`var ~ Distribution(params)`" syntax, aside from doing the usual programming stuff. 
"""

# ╔═╡ 802e2769-613c-4ab6-b0e5-38173c220042
md"Our circle problem can be defined as follows:"

# ╔═╡ 9ea69ca8-328d-4d37-b3f5-40206353d91c
@model function distances()
    x ~ Uniform(-1, 1)
    y ~ Uniform(-1, 1)
    dist = sqrt(x^2 + y^2)
    return dist
end

# ╔═╡ b61eb0af-5f0c-4267-8c02-ab14acb31ece
md"Calling this function will return a Turing model:"

# ╔═╡ fbf3e659-2b59-4f97-8f88-dc783d384d75
dist_model = distances();

# ╔═╡ f88af70b-0af7-48f5-91dd-39d1bb42f409
md"### Solving the question"

# ╔═╡ 03840c86-2d22-4288-a202-93f4d141ec47
md"There's a number of things we can do with the model. The most simple is calling it, which will give the return value of the function after sampling a value for all random variables, here `x` and `y`."

# ╔═╡ df921cec-a4a4-40de-9537-194b52782411
dist_model()

# ╔═╡ 7d9c1044-98f5-4876-ba7e-a2d70de92597
md"We can use this to generate a large number of samples and make estimations about the probability"

# ╔═╡ 62952198-602c-49a3-81d5-890588b7262a
n_samples = 1_000

# ╔═╡ 1012cee3-e574-43d9-b6ed-c9aa8a5ec552
sp_dists = [dist_model() for sample_idx in 1:n_samples];

# ╔═╡ 4084cee7-1d2d-4fc0-8226-d9d416bb4eec
histogram(sp_dists, title = "Distances of points to origin", bins = 20, legend = nothing)

# ╔═╡ 321c50fc-36ff-438c-b98e-a14910714bca
md"Currently, we have samples of the distance to the origin. We can easily transform these to samples of being inside the circle or not, and subsequently estimate the desired probability."

# ╔═╡ 1df2ac6a-fa56-47e2-8c10-2a5d4da3a475
sp_inside = sp_dists .<= 1; # the circle has a radius of 1

# ╔═╡ 8c95747f-08dc-49e7-a9c7-9a5c5fef7a78
begin
	histogram(sp_dists[sp_inside], title = "Distances of points to origin", label = "Inside circle", bins = 15)
	histogram!(sp_dists[.!sp_inside], label = "Not inside circle", bins = 5)
end

# ╔═╡ 7d654226-04bd-4721-bc4b-26912fc87ba7
prob_inside = length(sp_dists[sp_inside]) / length(sp_dists)

# ╔═╡ f2bb0f76-836c-40eb-8fa9-f4fcf8c03671
mean(sp_inside) # shorter alternative

# ╔═╡ 73278ad6-00f2-4e57-8c17-4e0a47be5d57
4*prob_inside

# ╔═╡ 19162953-a1d6-4b3b-962e-947e36b032c5
md"""
!!! note
	How many samples do you need to get that beautiful `3.14` consistently? How about the yet even more charming `3.1415`?
"""

# ╔═╡ c79118af-a437-4979-a857-a6c84e1f789a
md"### Keeping track of the random variables"

# ╔═╡ 55206a66-fda0-4299-95f0-6dcc32287f3a
md"""
An alternate way to generate a number of samples from our model is to use the `sample` function. Rather than getting samples of the function's output, this returns sampled values of _all the random variables_.
"""

# ╔═╡ 8ef27c2f-31fe-4411-9d2f-a65041e07641
dist_chain = sample(dist_model, Prior(), n_samples)

# ╔═╡ ee6a2c1d-5522-4870-a0e1-664f8ccd5d8f
scatter(dist_chain[:x], dist_chain[:y], aspect_ratio = :equal, label = "Dart locations")

# ╔═╡ 563e0fa0-25b9-4b36-b213-c952199caa32
md"The sample values of a random variable can be acquired by indexing the resulting chain with the variable's name as a `Symbol` or `String`:"

# ╔═╡ 3994fc44-5fda-47dc-8cc6-e2550d9c237d
dist_chain[:x] # or dist_chain["x"]

# ╔═╡ eb4b6d56-1e41-4e75-914b-8e7444afb288
md"This can be useful for making plots, for example"

# ╔═╡ 3bfaf571-2aa8-4b13-8c48-e8434407597b
md"What if we want the function's return value too? We could calculate it based on our random variables by hand as `sqrt.(dist_chain[:x].^2 + dist_chain[:y].^2)`, or use the `generated_quantities` function."

# ╔═╡ 0dcc8d6f-3ab2-429d-bcf4-f022ea5d0124
sp_dists_alt = generated_quantities(dist_model, dist_chain);

# ╔═╡ 5a8adeb4-46cf-486f-adda-24661a79b2e9
sp_inside_alt = sp_dists_alt .<= 1;

# ╔═╡ 28c48f45-726e-4352-b3b2-416dbab9cb0d
scatter(dist_chain[:x], dist_chain[:y], aspect_ratio = :equal, groups = vec(sp_inside_alt), label = ["Outside of circle" "Inside of circle"]) 
# Note: The `sample` method returns matrices. For plotting, vectors are often preferred, which is why we convert `sp_inside_alt` to a vector here.

# ╔═╡ ab230de4-767f-43e8-8bdc-054234852715
md"### Dealing with many variables"

# ╔═╡ bf0904ce-fb72-4a39-a57b-1c99c8e9c82e
md"A common problem when defining the problem as a Turing model is many random variables being involved, often with the same distribution. Turing allows variables to be defined in a for-loop for this reason."

# ╔═╡ 4c54bccf-c674-41da-be15-edac57f96ee9
@model function distances_loop()
	coords = zeros(2)
	for i in eachindex(coords)
		coords[i] ~ Uniform(-1, 1)
	end

	dist = coords.^2 |> sum |> sqrt
end

# ╔═╡ c50b7bb9-57e9-4f86-a0fa-0a68488b2fa9
distloop_model = distances_loop();

# ╔═╡ 330531d8-c0cf-4d0f-84b3-42acf9e30b39
sp_loop = [distloop_model() for i in 1:n_samples];

# ╔═╡ 9d786f74-73de-4fe3-8842-fa394641ce29
4*mean(sp_loop .<= 1)

# ╔═╡ 6faecd7e-fa99-4616-86c5-2561985f98a1
md"### Working with Distributions"

# ╔═╡ 9621ceec-c532-4ef5-883d-da67af13bfa1
md"""
Under the hood, Turing makes use of Julia's `Distributions` package. Knowing some basic functionality of this package can be useful.
"""

# ╔═╡ 996e876a-8e2d-44d4-bd3a-7a8d4fd8c6ea
md"""
!!! note
	Turing automatically loads Distributions into the workspace, so "`using Distributions`" is not necessary when Turing has been loaded.
"""

# ╔═╡ e8e65be5-b6b2-4fd6-bfb3-c277433b7206
md"Considering the humble example of `X ~ Exponential(10)`, let's do some plotting, sampling and calculating."

# ╔═╡ 398e15bb-7cc4-4f61-919b-9edf790cacd3
Exponential(10) |> plot

# ╔═╡ 58a6ffb1-0fef-4edb-9fbe-de5b399b013d
spX = rand(Exponential(10), n_samples)

# ╔═╡ 6b3ac000-d31b-46ef-89bd-830c2cae8cdb
histogram(spX)

# ╔═╡ 84b9eb6b-d18a-4394-afff-3e49293aa1d4
pdf(Exponential(10), 0)

# ╔═╡ b577685b-00fb-4869-b930-1ec6cd163f7d
cdf(Exponential(10), 20)

# ╔═╡ 03ebb0a1-4fff-4d87-8d86-4d409d09363d
md"Just for fun, we can work out the circle example again without Turing."

# ╔═╡ 48fe2927-8c79-4482-8edc-21d8e0819619
begin
	
sp_dists_noturing = zeros(n_samples)
for i in 1:n_samples
    x = rand(Uniform(0, 1))
	y = rand(Uniform(0, 1))
    dist = x^2 + y^2
    sp_dists_noturing[i] = dist
end

sp_inside_noturing = sp_dists_noturing .<= 1
4 * mean(sp_inside_noturing)

end

# ╔═╡ Cell order:
# ╟─aeb0aef0-b2ee-11ef-3cca-7f80b487ea17
# ╟─30957e05-85b8-4106-9635-82a5c11d9825
# ╟─13636a8a-8c31-4397-8a7d-5cd7f899d7a5
# ╠═a2410616-5a17-403c-aa2f-dc93c2633c7f
# ╠═8307092d-368d-441d-8315-3dc312026534
# ╟─9f89e350-199f-4875-947b-61df653ffc19
# ╟─2e544855-16b8-4794-9ba7-70a1e7209dd2
# ╟─88b3087c-f78a-4df4-a9dc-5696bf4052d5
# ╟─78f872c6-013b-4ad1-966f-9e0ce3288019
# ╟─bde6033f-15a8-4716-a41e-080f5d48e9d6
# ╟─1724e4e4-18d2-430d-a8b0-11f5891b09a3
# ╟─5acc6791-7cb6-4a31-a28d-649e761329ee
# ╟─f3a73486-485a-4c69-b3db-2153b5a06bd8
# ╠═6b554b3e-3fc8-4a20-9aac-a2b351dd7f9e
# ╠═34925f20-0fc3-4bf3-8a8d-89a9596504c4
# ╠═0712f06a-039a-4735-8501-86d3fe8d9114
# ╠═c6f18791-892d-4915-b4e4-f540fb770104
# ╠═5c7cf0d6-2d6b-487e-9609-546464f730a1
# ╠═c52c7639-c400-4093-a1e2-a0473059766c
# ╟─ccdcfa65-a02f-4110-a664-36090c3291d8
# ╟─43e8f794-fec2-4ce4-9306-2fc8a9565343
# ╟─ed788556-1627-4e9a-b901-e532272a8265
# ╟─802e2769-613c-4ab6-b0e5-38173c220042
# ╠═9ea69ca8-328d-4d37-b3f5-40206353d91c
# ╟─b61eb0af-5f0c-4267-8c02-ab14acb31ece
# ╠═fbf3e659-2b59-4f97-8f88-dc783d384d75
# ╟─f88af70b-0af7-48f5-91dd-39d1bb42f409
# ╟─03840c86-2d22-4288-a202-93f4d141ec47
# ╠═df921cec-a4a4-40de-9537-194b52782411
# ╟─7d9c1044-98f5-4876-ba7e-a2d70de92597
# ╠═62952198-602c-49a3-81d5-890588b7262a
# ╠═1012cee3-e574-43d9-b6ed-c9aa8a5ec552
# ╟─4084cee7-1d2d-4fc0-8226-d9d416bb4eec
# ╟─321c50fc-36ff-438c-b98e-a14910714bca
# ╠═1df2ac6a-fa56-47e2-8c10-2a5d4da3a475
# ╟─8c95747f-08dc-49e7-a9c7-9a5c5fef7a78
# ╠═7d654226-04bd-4721-bc4b-26912fc87ba7
# ╠═f2bb0f76-836c-40eb-8fa9-f4fcf8c03671
# ╠═73278ad6-00f2-4e57-8c17-4e0a47be5d57
# ╟─19162953-a1d6-4b3b-962e-947e36b032c5
# ╟─c79118af-a437-4979-a857-a6c84e1f789a
# ╟─55206a66-fda0-4299-95f0-6dcc32287f3a
# ╠═8ef27c2f-31fe-4411-9d2f-a65041e07641
# ╟─ee6a2c1d-5522-4870-a0e1-664f8ccd5d8f
# ╟─563e0fa0-25b9-4b36-b213-c952199caa32
# ╠═3994fc44-5fda-47dc-8cc6-e2550d9c237d
# ╟─eb4b6d56-1e41-4e75-914b-8e7444afb288
# ╟─3bfaf571-2aa8-4b13-8c48-e8434407597b
# ╠═0dcc8d6f-3ab2-429d-bcf4-f022ea5d0124
# ╠═5a8adeb4-46cf-486f-adda-24661a79b2e9
# ╟─28c48f45-726e-4352-b3b2-416dbab9cb0d
# ╟─ab230de4-767f-43e8-8bdc-054234852715
# ╟─bf0904ce-fb72-4a39-a57b-1c99c8e9c82e
# ╠═4c54bccf-c674-41da-be15-edac57f96ee9
# ╠═c50b7bb9-57e9-4f86-a0fa-0a68488b2fa9
# ╠═330531d8-c0cf-4d0f-84b3-42acf9e30b39
# ╠═9d786f74-73de-4fe3-8842-fa394641ce29
# ╟─6faecd7e-fa99-4616-86c5-2561985f98a1
# ╟─9621ceec-c532-4ef5-883d-da67af13bfa1
# ╟─996e876a-8e2d-44d4-bd3a-7a8d4fd8c6ea
# ╟─e8e65be5-b6b2-4fd6-bfb3-c277433b7206
# ╠═398e15bb-7cc4-4f61-919b-9edf790cacd3
# ╠═58a6ffb1-0fef-4edb-9fbe-de5b399b013d
# ╠═6b3ac000-d31b-46ef-89bd-830c2cae8cdb
# ╠═84b9eb6b-d18a-4394-afff-3e49293aa1d4
# ╠═b577685b-00fb-4869-b930-1ec6cd163f7d
# ╟─03ebb0a1-4fff-4d87-8d86-4d409d09363d
# ╠═48fe2927-8c79-4482-8edc-21d8e0819619
