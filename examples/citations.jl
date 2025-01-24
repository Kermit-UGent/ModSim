### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 2fb807d8-f555-41cf-a374-8ff9cadf9533
using Pkg; Pkg.activate("..")

# ╔═╡ 6613efca-d42c-11ef-0e77-716f85ded3df
using Turing, PlutoUI, StatsPlots

# ╔═╡ e89d9c26-0ffe-48ca-ad39-93655d6868a4
md"Consider the following model:"

# ╔═╡ 19fcf183-322a-44c4-8c3c-987f667d882e
@model function citations()
	X ~ Normal(0, 1)
	Y ~ Normal(1, 2)
	Z ~ Poisson(exp(0.8X + 1.2Y + 0.2X*Y))
end

# ╔═╡ 8197a2a1-86fb-44eb-a09f-d4e8511167b2
md"Generate a sample of 10,000 observations of this distribution. 

Use this sample to create a histogram of `Z` and compute its quantiles. Additionally, compute the 99% quantile (this is the top-1% value of `Z`).

Next, compute the correlation between `X` and `Y` and draw a scatter plot. Does this match your observation?
"

# ╔═╡ aa73f4fa-23ca-41ff-b3cf-4bab0740d157
prior = sample(citations(), Prior(), 10_000)

# ╔═╡ e40b92ad-2a9a-414e-b8ee-cc658d505212
summarize(prior)

# ╔═╡ fba7efbd-8ad0-4fa3-aeb9-5a674e0295ce
histogram(prior[:Z])

# ╔═╡ a737c4c1-c27d-4800-9d0a-e71d4503025b
histogram(prior[:Z], yscale=:log10)

# ╔═╡ b1d4d0a6-6d1d-4cf2-b21e-51705c167b86
quantile(prior, q = [0.025, 0.25, 0.5, 0.75, 0.975, 0.99])

# ╔═╡ 4196765c-b2bd-4093-963b-82fb0ec2e343
cor(prior[:X], prior[:Y])

# ╔═╡ 5b4ef7a2-790e-46ba-89f1-136a1c9e7fdb
scatter(prior[:X], prior[:Y], xlab="X", ylab="Y")

# ╔═╡ b0bef20d-a10c-457b-ac8c-1e8caddc54e1
md"Now, make a conditional model `high_citation` in which you fix `Z = 100`. Use the No-U-Turn sampler to obtain 10,000 posterior samples.

Using this conditional sample, again compute the correlation between `X` and `Y` together with the scatter plot. What do you observe?
"

# ╔═╡ bb1681c5-ff86-497c-8cb8-f5e2bd638536
high_citation = citations() | (Z=100,)

# ╔═╡ 03f48a27-f0ce-46e7-989d-ab57cfc5babe
post = sample(high_citation, NUTS(), 10_000)

# ╔═╡ d5e6bd69-6a49-4d6a-8c37-9cda7fb97bb5
quantile(post)

# ╔═╡ 72577ba4-e20c-4eb1-af34-1d702957f1e3
cor(post[:X], post[:Y])

# ╔═╡ 0a64b09b-fd71-4628-8b7d-7df3c4b45a15
scatter(post[:X], post[:Y], xlabel="X", ylabel="Y")

# ╔═╡ bf8ac4d6-5a6e-49ff-b75e-6a0f7ad7ff4d
md"Now for some context, the model is represents the number of citations a journal article can attract after being published (`Z`). `X` represents a metric of **quality**, how well the study was done and whether everything was done correctly. `Y` represents the spectacularity or excitedness of the claims, i.e., whether its results would be of general interest. Re-interpret you results in light of this new information."

# ╔═╡ 830e9980-3f01-4f54-b29f-7247bb4531f4
plot(post)

# ╔═╡ Cell order:
# ╠═6613efca-d42c-11ef-0e77-716f85ded3df
# ╠═2fb807d8-f555-41cf-a374-8ff9cadf9533
# ╠═e89d9c26-0ffe-48ca-ad39-93655d6868a4
# ╠═19fcf183-322a-44c4-8c3c-987f667d882e
# ╠═8197a2a1-86fb-44eb-a09f-d4e8511167b2
# ╠═aa73f4fa-23ca-41ff-b3cf-4bab0740d157
# ╠═e40b92ad-2a9a-414e-b8ee-cc658d505212
# ╠═fba7efbd-8ad0-4fa3-aeb9-5a674e0295ce
# ╠═a737c4c1-c27d-4800-9d0a-e71d4503025b
# ╠═b1d4d0a6-6d1d-4cf2-b21e-51705c167b86
# ╠═4196765c-b2bd-4093-963b-82fb0ec2e343
# ╠═5b4ef7a2-790e-46ba-89f1-136a1c9e7fdb
# ╠═b0bef20d-a10c-457b-ac8c-1e8caddc54e1
# ╠═bb1681c5-ff86-497c-8cb8-f5e2bd638536
# ╠═03f48a27-f0ce-46e7-989d-ab57cfc5babe
# ╠═d5e6bd69-6a49-4d6a-8c37-9cda7fb97bb5
# ╠═72577ba4-e20c-4eb1-af34-1d702957f1e3
# ╠═0a64b09b-fd71-4628-8b7d-7df3c4b45a15
# ╠═bf8ac4d6-5a6e-49ff-b75e-6a0f7ad7ff4d
# ╠═830e9980-3f01-4f54-b29f-7247bb4531f4
