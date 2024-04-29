### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 3e70b82a-e4d3-4747-9679-aa5ab41d5b06
begin
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ ab31fa80-0606-11ef-06f5-bd74bc457ff5
using Plots, PlutoUI, LaTeXStrings, Latexify

# ╔═╡ a10a9f8e-4088-49e2-a497-b1d253d83ecc
using Optim, Turing, StatsPlots, StatsBase

# ╔═╡ 470ba643-ed55-4789-9384-c88982ca0c05
t = [0, 1, 5, 8, 10]

# ╔═╡ 03457e54-33ae-4cf8-a4a9-85dab26c1ef9
y = [1.2, 3, 4, 5, 5.4]

# ╔═╡ 0c6aebf8-b18f-4987-b13e-e687436af550
scatter(t, y)

# ╔═╡ 5dee9551-a275-4c92-8a35-9bf7d2a91e7d
logistic(t; r, K, y0=1) = K / (1 + (K-y0) / y0 * exp(-r*t))

# ╔═╡ 6984bd57-e95f-40b1-aecc-a864177bdb72
sq_loss(r, K) = sum(abs2, logistic.(t; r, K) .- y)

# ╔═╡ e4fd9a2f-92be-4d4f-a3b0-475848e7aed3
contourf(0.01:.01:5, 1:.01:10, sq_loss, color=:speed)

# ╔═╡ d43ee840-dfc5-4d2d-bb0f-14faf1573e73
@model function logistic_model(t, y)
	σ ~ InverseGamma()
	r ~ TriangularDist(0.1, 5, 2)
	K ~ TriangularDist(1, 10, 5)
	for i in 1:length(y)
		y[i] ~ Normal(logistic(t[i]; r, K), σ)
	end
end

# ╔═╡ 3b0ca57f-d7f3-4bbf-91b8-23d488c3d3cd
logistic_fit = logistic_model(t, y)

# ╔═╡ 7cff29b5-3478-49ae-abea-712cfbee762d
ll(r, K) = loglikelihood(logistic_fit, (;r, K, σ=0.5))

# ╔═╡ 8f08acd8-4480-418d-89dd-64f795be6347
contourf(0.01:.01:5, 1:.01:10, ll, color=cgrad(:speed, rev=true))

# ╔═╡ 9b82503e-48a0-4c78-b7e7-0102b2d1f9d5
ll(0.1, 2)

# ╔═╡ 1a283f5a-ad39-4f16-b2b2-ff531f1a9634
ml_log = optimize(logistic_fit, MLE(), BFGS())

# ╔═╡ 4acefe9c-8ea1-4491-9f03-0be80c5211c3
coeftable(ml_log)

# ╔═╡ e2eaca92-3b18-4448-967f-ec78d165bd68
lprior(r, K) = logprior(logistic_fit, (;r, K, σ=0.5))

# ╔═╡ cc104cdb-47d0-4b6b-adc4-103a6dad5421
contourf(0.01:.01:5, 1:.01:10, lprior, color=cgrad(:speed, rev=true))

# ╔═╡ f55335b3-c0ec-467f-9056-c5416b662088
lp(r, K) = ll(r, K) + lprior(r, K)

# ╔═╡ 99399d90-29cc-4648-871d-67adbebfec8a
contourf(0.01:.01:5, 1:.01:10, lp, color=cgrad(:speed, rev=true))

# ╔═╡ e7e7dd46-835a-4cf1-aec5-e5ae1d2b1b1b
map_log = optimize(logistic_fit, MAP(), BFGS())

# ╔═╡ 41c6b352-77f9-4480-be19-20db633c7a6e
map_log.values[:K]

# ╔═╡ c4652a3b-630d-4eb3-a5f0-63034f077730


# ╔═╡ b317a124-a044-48e0-af65-68be95593a6a


# ╔═╡ 087e8eb2-019e-4ffe-a520-636d4f13c01f


# ╔═╡ f5819e00-2e6a-4886-b363-f23ddc5abd72
coeftable(map_log)

# ╔═╡ 16b6935f-1a37-4b12-ac18-5803175555a5
chain = sample(logistic_fit, NUTS(), 10_000)

# ╔═╡ 255c791b-0012-4a04-b794-f9ffe1c62a26
summarize(chain)

# ╔═╡ 290d4176-5715-4b49-85d5-c90bde256b00
quantile(chain)

# ╔═╡ 80d28884-c66e-4384-ac44-d8768d8c810f
plot(chain)

# ╔═╡ 63ea89c5-df05-4651-b655-10d678fde222


# ╔═╡ dd7415cd-2646-4135-a389-baf8d6657fe1
md"## Appendix"

# ╔═╡ 2f28d23a-36a3-4831-a8df-ce08cfd7c44e
TableOfContents()

# ╔═╡ 55459b43-4e53-4f44-80d7-84c1f929fa52
plots = Dict()

# ╔═╡ a2be0a77-fca0-4d7f-af15-254faaa7e336
plots

# ╔═╡ Cell order:
# ╠═ab31fa80-0606-11ef-06f5-bd74bc457ff5
# ╠═a10a9f8e-4088-49e2-a497-b1d253d83ecc
# ╠═3e70b82a-e4d3-4747-9679-aa5ab41d5b06
# ╠═470ba643-ed55-4789-9384-c88982ca0c05
# ╠═03457e54-33ae-4cf8-a4a9-85dab26c1ef9
# ╠═0c6aebf8-b18f-4987-b13e-e687436af550
# ╠═5dee9551-a275-4c92-8a35-9bf7d2a91e7d
# ╠═6984bd57-e95f-40b1-aecc-a864177bdb72
# ╠═e4fd9a2f-92be-4d4f-a3b0-475848e7aed3
# ╠═d43ee840-dfc5-4d2d-bb0f-14faf1573e73
# ╠═3b0ca57f-d7f3-4bbf-91b8-23d488c3d3cd
# ╠═7cff29b5-3478-49ae-abea-712cfbee762d
# ╠═8f08acd8-4480-418d-89dd-64f795be6347
# ╠═9b82503e-48a0-4c78-b7e7-0102b2d1f9d5
# ╠═1a283f5a-ad39-4f16-b2b2-ff531f1a9634
# ╠═4acefe9c-8ea1-4491-9f03-0be80c5211c3
# ╠═e2eaca92-3b18-4448-967f-ec78d165bd68
# ╠═cc104cdb-47d0-4b6b-adc4-103a6dad5421
# ╠═f55335b3-c0ec-467f-9056-c5416b662088
# ╠═99399d90-29cc-4648-871d-67adbebfec8a
# ╠═e7e7dd46-835a-4cf1-aec5-e5ae1d2b1b1b
# ╠═41c6b352-77f9-4480-be19-20db633c7a6e
# ╠═c4652a3b-630d-4eb3-a5f0-63034f077730
# ╠═b317a124-a044-48e0-af65-68be95593a6a
# ╠═087e8eb2-019e-4ffe-a520-636d4f13c01f
# ╠═f5819e00-2e6a-4886-b363-f23ddc5abd72
# ╠═16b6935f-1a37-4b12-ac18-5803175555a5
# ╠═255c791b-0012-4a04-b794-f9ffe1c62a26
# ╠═290d4176-5715-4b49-85d5-c90bde256b00
# ╠═80d28884-c66e-4384-ac44-d8768d8c810f
# ╠═63ea89c5-df05-4651-b655-10d678fde222
# ╠═dd7415cd-2646-4135-a389-baf8d6657fe1
# ╠═2f28d23a-36a3-4831-a8df-ce08cfd7c44e
# ╠═55459b43-4e53-4f44-80d7-84c1f929fa52
# ╠═a2be0a77-fca0-4d7f-af15-254faaa7e336
