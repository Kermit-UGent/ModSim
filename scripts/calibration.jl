### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 3e70b82a-e4d3-4747-9679-aa5ab41d5b06
begin
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ ab31fa80-0606-11ef-06f5-bd74bc457ff5
using Plots, PlutoUI, LaTeXStrings, Latexify, LinearAlgebra, Random

# ╔═╡ a10a9f8e-4088-49e2-a497-b1d253d83ecc
using Optim, Turing, StatsPlots, StatsBase

# ╔═╡ 4b38ba22-8c3f-4c90-b42a-79063489641e
using Catalyst, DifferentialEquations

# ╔═╡ 0b6b3b7f-98e7-4922-82de-6a98453a627c
md"""

## Polynomial regression

Let us fit a polynomial regression model of the form

$$f(x) = \sum^{p-1}_{i=0} \beta_i \frac{x^i}{i!}\,.$$

"""

# ╔═╡ 1b055cac-6192-438f-8f3f-9891f5f06f49
polynom(θ) = x -> sum(θi * x^(i-1)/factorial(i-1) for (i, θi) in enumerate(θ))

# ╔═╡ ecceb6b5-74f4-4b5a-9c73-4ae1477b2c47
β = [2, 3, -6, -2, 3, 0, 0, 0, 0, 0]

# ╔═╡ 4ca12cc9-3bae-4535-a5c9-11cac5b2c89f
p = length(β)

# ╔═╡ ac3accfd-7edb-4bce-9891-d544210867e5
poly = polynom(β)

# ╔═╡ 6daf07bc-807b-4dc5-8492-fc6f107cb294
plot(poly, -5, 5)

# ╔═╡ 4f5f7c4f-ff5b-450c-901d-e7395f6623dd
xpoly = [-4, -3, -3.4, -2, -1, 0, 1, 1.2, 3.5, 4]

# ╔═╡ 81e20a6b-4309-464e-b57b-05a38558515f
σpoly = 3

# ╔═╡ bdf2debf-5bd7-4183-a428-109b01fc7488
λ = 0.05

# ╔═╡ 66527df1-bd8b-45d6-b499-b5d17db84c42
Apoly = [x^i/factorial(i) for x in xpoly, i in 0:p-1]

# ╔═╡ 0e5831ff-f542-40e1-9376-2abb420494b7
md"[mm?](https://chem.libretexts.org/Bookshelves/Biological_Chemistry/Supplemental_Modules_(Biological_Chemistry)/Enzymes/Enzymatic_Kinetics/Michaelis-Menten_Kinetics)
MM?"

# ╔═╡ 894bfbf5-26c0-43e2-8f5b-fef4c41e3438
# ╠═╡ disabled = true
#=╠═╡
data = [(14.4, 1.0), (15.9, 2.3), (17.3, 4.0), (18.7, 6.0), (20.1, 7.9), (21.6, 9.2), (23.1, 9.8), (24.6, 9.9)]
  ╠═╡ =#

# ╔═╡ 03457e54-33ae-4cf8-a4a9-85dab26c1ef9
#=╠═╡
t, y = first.(data), last.(data)
  ╠═╡ =#

# ╔═╡ 0c6aebf8-b18f-4987-b13e-e687436af550
#=╠═╡
scatter(t, y)
  ╠═╡ =#

# ╔═╡ 5dee9551-a275-4c92-8a35-9bf7d2a91e7d
logistic(t; r, K, y₀=.1) = K / (1 + (K-y₀) / y₀ * exp(-r*t))

# ╔═╡ 6984bd57-e95f-40b1-aecc-a864177bdb72
#=╠═╡
sq_loss(r, K) = sum(abs2, logistic.(t; r, K) .- y) |> sqrt
  ╠═╡ =#

# ╔═╡ e4fd9a2f-92be-4d4f-a3b0-475848e7aed3
#=╠═╡
contourf(0.01:.01:2, 1:.01:20, sq_loss, color=:speed)
  ╠═╡ =#

# ╔═╡ d43ee840-dfc5-4d2d-bb0f-14faf1573e73
@model function logistic_model(t, y; y₀=missing)
	σ ~ InverseGamma()
	r ~ Exponential(.1)
	K ~ Exponential(5)
	y₀ ~ Exponential(1)
	for i in 1:length(y)
		y[i] ~ Normal(logistic(t[i]; r, K, y₀), σ)
	end
end

# ╔═╡ 3b0ca57f-d7f3-4bbf-91b8-23d488c3d3cd
#=╠═╡
logistic_fit = logistic_model(t, y)
  ╠═╡ =#

# ╔═╡ 7cff29b5-3478-49ae-abea-712cfbee762d
#=╠═╡
ll(r, K) = loglikelihood(logistic_model(t, y; y₀=0.1), (;r, K, y₀=0.000213524, σ=0.5))
  ╠═╡ =#

# ╔═╡ 8f08acd8-4480-418d-89dd-64f795be6347
#=╠═╡
contourf(0.01:.01:2, 1:.01:20, ll, color=cgrad(:speed, rev=true))
  ╠═╡ =#

# ╔═╡ 9b82503e-48a0-4c78-b7e7-0102b2d1f9d5
#=╠═╡
ll(0.1, 2)
  ╠═╡ =#

# ╔═╡ 1a283f5a-ad39-4f16-b2b2-ff531f1a9634
#=╠═╡
ml_log = optimize(logistic_fit, MLE(), NelderMead())
  ╠═╡ =#

# ╔═╡ 4acefe9c-8ea1-4491-9f03-0be80c5211c3
#=╠═╡
coeftable(ml_log)
  ╠═╡ =#

# ╔═╡ e2eaca92-3b18-4448-967f-ec78d165bd68
#=╠═╡
lprior(r, K) = logprior(logistic_model(t, y; y₀=0.1), (;r, K, σ=0.5))
  ╠═╡ =#

# ╔═╡ cc104cdb-47d0-4b6b-adc4-103a6dad5421
#=╠═╡
contourf(0.01:.01:5, 1:.01:10, lprior, color=cgrad(:speed, rev=true))
  ╠═╡ =#

# ╔═╡ f55335b3-c0ec-467f-9056-c5416b662088
#=╠═╡
lp(r, K) = ll(r, K) + lprior(r, K)
  ╠═╡ =#

# ╔═╡ 99399d90-29cc-4648-871d-67adbebfec8a
#=╠═╡
contourf(0.01:.01:5, 1:.01:10, lp, color=cgrad(:speed, rev=true))
  ╠═╡ =#

# ╔═╡ e7e7dd46-835a-4cf1-aec5-e5ae1d2b1b1b
#=╠═╡
map_log = optimize(logistic_fit, MAP(), NelderMead())
  ╠═╡ =#

# ╔═╡ 41c6b352-77f9-4480-be19-20db633c7a6e
#=╠═╡
map_log.values[:K]
  ╠═╡ =#

# ╔═╡ f5819e00-2e6a-4886-b363-f23ddc5abd72
#=╠═╡
coeftable(map_log)
  ╠═╡ =#

# ╔═╡ 16b6935f-1a37-4b12-ac18-5803175555a5
#=╠═╡
chain = sample(logistic_fit, NUTS(), 10_000);
  ╠═╡ =#

# ╔═╡ 255c791b-0012-4a04-b794-f9ffe1c62a26
#=╠═╡
summarize(chain)
  ╠═╡ =#

# ╔═╡ 290d4176-5715-4b49-85d5-c90bde256b00
#=╠═╡
quantile(chain)
  ╠═╡ =#

# ╔═╡ 80d28884-c66e-4384-ac44-d8768d8c810f
#=╠═╡
plot(chain)
  ╠═╡ =#

# ╔═╡ 63ea89c5-df05-4651-b655-10d678fde222
@model function polynomial_regression(x, y, m)
	σm ~ InverseGamma(1)  # standard deviation of the error
	λ ~ InverseGamma(1)
	n = length(y)
	β ~ MultivariateNormal(m + 1, λ)
	x_stand = (x .- mean(x)) ./ std(x)
	for i in 1:n
		yv = 0.0
		for j in 0:m
			yv += β[j+1] * (x_stand[i])^j / factorial(j)
		end
		y[i] ~ Normal(yv, σm)
	end
end

# ╔═╡ 9b3e76fc-0c0b-4bd0-8af8-5e6053e24094
xp = 10rand(200) .- 5 |> sort!

# ╔═╡ 960814cb-611e-4057-8667-0c4de798a96e
yp = 3randn(200) .+ 4 .- 4xp .+ 8xp.^2

# ╔═╡ 9b317909-8997-46e2-809d-b005332c1022
sample(polynomial_regression(xp, yp, 5), NUTS(), 1000) |> summarize

# ╔═╡ bf80bb42-75c8-4064-952c-0613e3e8188d
optimize(polynomial_regression(xp, yp, 2), MLE(), NelderMead())

# ╔═╡ 5ff4ccf6-9bd8-4318-9de3-7752ed66f51a
optimize(polynomial_regression(xp, yp, 5), MAP()) |> coeftable

# ╔═╡ 63214fc8-2c39-4cbe-8aa7-da36e904e978
yeast = @reaction_network begin
	X * mm(G, μ, K), X + G => 2X
	m, X --> 0 
end

# ╔═╡ f241d041-34a1-4399-a724-8be968ed4677
yeast_ode = convert(ODESystem, yeast)

# ╔═╡ 0a283dff-c28e-4349-ad3f-a07b16bf274b
parameters(yeast_ode)

# ╔═╡ 0a461157-58f1-4478-9700-65bf5fb80ebe
prob_yeast = ODEProblem(yeast, [:X=>10.1, :G=>180], (0, 50), [:μ=>0.3, :K=>250, :m=>0.03])

# ╔═╡ 25ac2424-170c-463d-9a5e-c900f3ed8674
sol_yeast = solve(prob_yeast, Tsit5(), saveat=5)

# ╔═╡ 0e6ba65b-067c-42f4-90b7-8b8c387a2fa9
plot(sol_yeast)

# ╔═╡ 0c9c49af-59cc-4f2a-9de2-6095d379e18b
σ_X, σ_G = 7.4, 19

# ╔═╡ 018aa4ed-b4c2-4413-80c9-af4f34256865
length(sol_yeast)

# ╔═╡ a7105efd-a192-406e-949c-75c8d9e35da9
sol_yeast[:X,1]

# ╔═╡ 7e6dbcd1-be40-4315-87a9-129b7f8ce341
Xobs = sol_yeast[:X] .+ σ_X .* randn(length(sol_yeast))

# ╔═╡ 3fea2784-e662-4b01-896c-54bc913adf40
Gobs = sol_yeast[:G] .+ σ_G .* randn(length(sol_yeast))

# ╔═╡ 0f40a683-d4e7-4fc8-bf5e-18de8e3d6111
let
	scatter(sol_yeast.t, Xobs)
	scatter!(sol_yeast.t, Gobs)
end


# ╔═╡ 1ba6e0d0-3166-41ab-a106-df73e4579f39
tsteps = sol_yeast.t

# ╔═╡ 6185ec28-3286-40c8-bdcc-6d8e812bb7df
@model function yeast_inference(tsteps, X, G)
	σ_Xsq ~ InverseGamma()
	σ_Gsq ~ InverseGamma()
	μ ~ Uniform(0.01, 5)
	m ~ Uniform(0.01, 5)
	sol = solve(prob_yeast, Tsit5(), saveat=tsteps, p=[μ, 250.0, m])
	for i in 1:length(X)
		X[i] ~ Normal(sol[:X,i], sqrt(σ_Xsq))
		G[i] ~ Normal(sol[:G,i], sqrt(σ_Gsq))
	end
end

# ╔═╡ c4896649-a42c-4fd8-a8ef-ecc2eeee61d9
yeast_mod = yeast_inference(tsteps, Xobs, Gobs)

# ╔═╡ 95b92274-f1a2-4c4f-b058-d384211b9c1a
optimize(yeast_mod, MLE(), LBFGS()) |> coeftable

# ╔═╡ dd285391-6df6-47f7-a568-8786d2000f24
optimize(yeast_mod, MAP(), LBFGS()) |> coeftable

# ╔═╡ 16089b87-a760-4309-86b9-8300bb5a23cf
chain_yeast = sample(yeast_mod, NUTS(), MCMCSerial(), 5000, 5);

# ╔═╡ 05ed997f-e9b1-4f10-a186-78eac7f8aa8a
summarize(chain_yeast)

# ╔═╡ dd7415cd-2646-4135-a389-baf8d6657fe1
md"## Appendix"

# ╔═╡ 2f28d23a-36a3-4831-a8df-ce08cfd7c44e
TableOfContents()

# ╔═╡ 3af17761-3a77-4c83-a54b-b55e2a3eb708
rng = MersenneTwister(10)

# ╔═╡ 27eea16c-feff-4678-907e-f57fb7cf8c0e
ypoly = poly.(xpoly) + σpoly .* randn(rng,length(xpoly))

# ╔═╡ 6edf2f34-1e9e-4775-bd6a-624d003b6688
let
	scatter(xpoly, ypoly)
end

# ╔═╡ 70f8723f-1f85-4179-844d-0c4d0130d5c4
βls = (Apoly' * Apoly) \ (Apoly' * ypoly)

# ╔═╡ d194b539-f4c1-4ad5-a721-48b6cf73afed
βtik = (Apoly' * Apoly + λ * p * I) \ (Apoly' * ypoly)

# ╔═╡ 55459b43-4e53-4f44-80d7-84c1f929fa52
plots = Dict()

# ╔═╡ 181a4886-a515-4f06-9434-64b7387ee0d3
let
	p_ls = polynom(βls)
	p_tik = polynom(βtik)
	p = scatter(xpoly, ypoly, label="data", xlab=L"x")
	ylims!(-30, 25)
	plot!(poly, -5, 5, lw=2, label="f(x)", alpha=0.8, legend=:bottom)
	plot!(p_ls, -5, 5, lw=2, label="Penrose-Moore pseudo inverse", ls=:dash)
	plot!(p_tik, -5, 5, lw=2, label="Tikhonov inverse (λ=$λ)", ls=:dashdot)
	title!("Polynomial regression")
	plots["poly_regr"] = p
end

# ╔═╡ a2be0a77-fca0-4d7f-af15-254faaa7e336
plots

# ╔═╡ Cell order:
# ╠═ab31fa80-0606-11ef-06f5-bd74bc457ff5
# ╠═a10a9f8e-4088-49e2-a497-b1d253d83ecc
# ╠═4b38ba22-8c3f-4c90-b42a-79063489641e
# ╠═3e70b82a-e4d3-4747-9679-aa5ab41d5b06
# ╟─0b6b3b7f-98e7-4922-82de-6a98453a627c
# ╠═1b055cac-6192-438f-8f3f-9891f5f06f49
# ╠═ecceb6b5-74f4-4b5a-9c73-4ae1477b2c47
# ╠═4ca12cc9-3bae-4535-a5c9-11cac5b2c89f
# ╠═ac3accfd-7edb-4bce-9891-d544210867e5
# ╠═6daf07bc-807b-4dc5-8492-fc6f107cb294
# ╠═4f5f7c4f-ff5b-450c-901d-e7395f6623dd
# ╠═81e20a6b-4309-464e-b57b-05a38558515f
# ╠═bdf2debf-5bd7-4183-a428-109b01fc7488
# ╠═27eea16c-feff-4678-907e-f57fb7cf8c0e
# ╠═6edf2f34-1e9e-4775-bd6a-624d003b6688
# ╠═66527df1-bd8b-45d6-b499-b5d17db84c42
# ╠═70f8723f-1f85-4179-844d-0c4d0130d5c4
# ╠═d194b539-f4c1-4ad5-a721-48b6cf73afed
# ╠═181a4886-a515-4f06-9434-64b7387ee0d3
# ╠═0e5831ff-f542-40e1-9376-2abb420494b7
# ╠═894bfbf5-26c0-43e2-8f5b-fef4c41e3438
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
# ╠═f5819e00-2e6a-4886-b363-f23ddc5abd72
# ╠═16b6935f-1a37-4b12-ac18-5803175555a5
# ╠═255c791b-0012-4a04-b794-f9ffe1c62a26
# ╠═290d4176-5715-4b49-85d5-c90bde256b00
# ╠═80d28884-c66e-4384-ac44-d8768d8c810f
# ╠═63ea89c5-df05-4651-b655-10d678fde222
# ╠═9b3e76fc-0c0b-4bd0-8af8-5e6053e24094
# ╠═960814cb-611e-4057-8667-0c4de798a96e
# ╠═9b317909-8997-46e2-809d-b005332c1022
# ╠═bf80bb42-75c8-4064-952c-0613e3e8188d
# ╠═5ff4ccf6-9bd8-4318-9de3-7752ed66f51a
# ╠═63214fc8-2c39-4cbe-8aa7-da36e904e978
# ╠═f241d041-34a1-4399-a724-8be968ed4677
# ╠═0a283dff-c28e-4349-ad3f-a07b16bf274b
# ╠═0a461157-58f1-4478-9700-65bf5fb80ebe
# ╠═25ac2424-170c-463d-9a5e-c900f3ed8674
# ╠═0e6ba65b-067c-42f4-90b7-8b8c387a2fa9
# ╠═0c9c49af-59cc-4f2a-9de2-6095d379e18b
# ╠═018aa4ed-b4c2-4413-80c9-af4f34256865
# ╠═a7105efd-a192-406e-949c-75c8d9e35da9
# ╠═7e6dbcd1-be40-4315-87a9-129b7f8ce341
# ╠═3fea2784-e662-4b01-896c-54bc913adf40
# ╠═0f40a683-d4e7-4fc8-bf5e-18de8e3d6111
# ╠═1ba6e0d0-3166-41ab-a106-df73e4579f39
# ╠═6185ec28-3286-40c8-bdcc-6d8e812bb7df
# ╠═c4896649-a42c-4fd8-a8ef-ecc2eeee61d9
# ╠═95b92274-f1a2-4c4f-b058-d384211b9c1a
# ╠═dd285391-6df6-47f7-a568-8786d2000f24
# ╠═16089b87-a760-4309-86b9-8300bb5a23cf
# ╠═05ed997f-e9b1-4f10-a186-78eac7f8aa8a
# ╠═dd7415cd-2646-4135-a389-baf8d6657fe1
# ╠═2f28d23a-36a3-4831-a8df-ce08cfd7c44e
# ╠═3af17761-3a77-4c83-a54b-b55e2a3eb708
# ╠═55459b43-4e53-4f44-80d7-84c1f929fa52
# ╠═a2be0a77-fca0-4d7f-af15-254faaa7e336
