### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 3e70b82a-e4d3-4747-9679-aa5ab41d5b06
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	using Pkg
	Pkg.activate("..")
end
  ╠═╡ =#

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
# ╠═╡ disabled = true
#=╠═╡
p = length(β)
  ╠═╡ =#

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

# ╔═╡ 27eea16c-feff-4678-907e-f57fb7cf8c0e
ypoly = poly.(xpoly) + σpoly .* randn(MersenneTwister(2), length(xpoly))

# ╔═╡ 66527df1-bd8b-45d6-b499-b5d17db84c42
Apoly = [x^i/factorial(i) for x in xpoly, i in 0:9]

# ╔═╡ 70f8723f-1f85-4179-844d-0c4d0130d5c4
βls = (Apoly' * Apoly) \ (Apoly' * ypoly)

# ╔═╡ d194b539-f4c1-4ad5-a721-48b6cf73afed
βtik = (Apoly' * Apoly + λ * 10 * I) \ (Apoly' * ypoly)

# ╔═╡ e735ea0e-0c98-4707-bc5e-9c2bba04aa90
md"## Michaelis-Menten"

# ╔═╡ 5b80b912-e625-44a4-9653-8306189718c4
μmax_star, Ks_star, σmm = 17.0, 79, 1.2

# ╔═╡ 2f8746e7-2e4a-4e40-bd48-f69349f1259d
Xs = Cs = 5:5:180

# ╔═╡ 17252072-b71f-4c30-b108-8163dfd63fbf
n_mm = length(Cs)

# ╔═╡ bf757553-7650-4b42-a07e-41007a3c29db
μs = mm.(Cs, μmax_star, Ks_star)  # non-noisy concentrations

# ╔═╡ 234bebb6-0d05-41a0-871d-50cf752d17a2
@model function michaelis_menten(Xs, μs)
	sigmasq ~ InverseGamma()
	mu_max ~ Truncated(Normal(20, 10), 0, 50)
	Ks ~ Uniform(5, 500)
	for i in 1:length(Xs)
		μs[i] ~ Normal(mm(Xs[i], mu_max, Ks), sqrt(sigmasq))
	end
	return x -> mm(x, mu_max, Ks)
end

# ╔═╡ 4a336948-d46e-415f-bcaf-07c667da4089
Truncated(Normal(20, 10), 0, 50)

# ╔═╡ 07d91dc7-2db0-464e-bc61-9633f2bfe92b
@model function michaelis_menten_robust(Xs, μs)
	sigmasq ~ InverseGamma()
	mu_max ~ Truncated(Normal(20, 10), 0, 50)
	Ks ~ Uniform(5, 500)
	for i in 1:length(Xs)
		μs[i] ~ Laplace(mm(Xs[i], mu_max, Ks), sqrt(sigmasq))
	end
	return x -> mm(x, mu_max, Ks)
end

# ╔═╡ 0e5831ff-f542-40e1-9376-2abb420494b7
md"[mm?](https://chem.libretexts.org/Bookshelves/Biological_Chemistry/Supplemental_Modules_(Biological_Chemistry)/Enzymes/Enzymatic_Kinetics/Michaelis-Menten_Kinetics)
MM?"

# ╔═╡ 03373997-a9eb-4d5a-9ef5-8822c3f95285
md"## Poisson"

# ╔═╡ 70f8de3c-5c5f-4ff0-92eb-e6d2a06451e1
@bind λ_pois Slider(0:0.2:50, show_value=true, default=26)

# ╔═╡ cd3746d4-518d-476e-abf4-2ac92a2c1387
import ForwardDiff

# ╔═╡ 6dc06148-79bf-44bf-8a5c-b6221b55bce0
xs_tiny = [2, 4, 3, 4]

# ╔═╡ d733961a-094a-4724-9b16-ee669a559945
ll_pois2(λ) = sum(x->loglikelihood(Poisson(λ), x), xs_tiny)

# ╔═╡ 46551486-d8a3-4aca-a92d-8c7f1c4439c0
pois_map(λ) = ll_pois2(λ) - log(λ)

# ╔═╡ 0d6cc24c-6008-4bb7-8fb9-45223febc05b
λ_star = optimize(x->-pois_map(x[1]), [mean(xs_tiny)], Newton()).minimizer[1]

# ╔═╡ 1e634f77-286b-44ce-b768-72b060ee898b
md"Fisher information using the second-order derivative:"

# ╔═╡ 0082f110-22df-49d5-8958-29e9516e9ca8
fi = -ForwardDiff.derivative(l->ForwardDiff.derivative(pois_map, l), λ_star)

# ╔═╡ 362fe960-f247-4d7b-ab9d-01b67d1e647e
σ_λ = inv(√(fi))

# ╔═╡ 8feddb40-5cf8-4d76-a959-dd3071c20e27
post_Laplace = Normal(λ_star, σ_λ)

# ╔═╡ da1aa2a7-e876-4b80-9204-d3ea11bc0f44
md"## LV model"

# ╔═╡ 73e018a3-27f6-4631-a8c3-394644cb774e
let
	# Define Lotka-Volterra model.
function lotka_volterra(du, u, p, t)
    # Model parameters.
    α, β, γ, δ = p
    # Current state.
    x, y = u

    # Evaluate differential equations.
    du[1] = (α - β * y) * x # prey
    du[2] = (δ * x - γ) * y # predator

    return nothing
end

# Define initial-value problem.
u0 = [1.0, 1.0]
p = [1.5, 1.0, 3.0, 1.0]
tspan = (0.0, 10.0)
prob = ODEProblem(lotka_volterra, u0, tspan, p)

# Plot simulation.
plot(solve(prob, Tsit5()))
end

# ╔═╡ e79a72fc-9929-4faa-a76d-ad914fc731ab
function lotka_volterra!(du, u, p, t)
	# parameters
	alpha, beta, gamma, delta = p
	# states
	x, y = u[1], u[2]
	# differential equations
	du[1] = dxdt = alpha * x - beta * x * y
	du[2] = dydt = delta * x * y - gamma * y
end

# ╔═╡ 4c9f01b3-aef4-47ed-9633-1abda691b531
alpha = 1.5

# ╔═╡ c3008b5a-7138-48a5-b744-d8f3aa2ff2ad
beta = 1.0

# ╔═╡ db3f89f9-3bef-432f-b558-3d32c5e79abe
gamma = 3.0

# ╔═╡ eb132faa-fedb-4e18-868f-46783e503515
delta = 1.0

# ╔═╡ 2a25b088-7414-48a2-a702-8c87213b92d3
p = (alpha, beta, gamma, delta)

# ╔═╡ a05d6800-a168-4e44-bfcb-0809aedda09e
lv_prob = ODEProblem(lotka_volterra!, [1.0, 1.0], (0.0, 10.0), p)

# ╔═╡ b3c31cca-b472-4dec-b700-11cbff9465a9
plot(solve(lv_prob, Tsit5()), lw=2)

# ╔═╡ 14810fe9-e6fe-46c7-b6cf-9705b5974113
begin
	sol = solve(lv_prob, Tsit5(); saveat=0.1)
	odedata = Array(sol) + 0.8 * randn(size(Array(sol)))

	# Plot simulation and noisy observations.
	plot(sol; alpha=0.3)
	scatter!(sol.t, odedata'; color=[1 2], label="", marker=[:o :^])
end

# ╔═╡ 79473eed-7fba-43e7-80a7-51d011931578
@model function fitlv(data, lv_prob)
    # Prior distributions.
    σ ~ InverseGamma(2, 3)
    α ~ truncated(Normal(1.5, 0.5); lower=0.5, upper=2.5)
    β ~ truncated(Normal(1.2, 0.5); lower=0, upper=2)
    γ ~ truncated(Normal(3.0, 0.5); lower=1, upper=4)
    δ ~ truncated(Normal(1.0, 0.5); lower=0, upper=2)

    # Simulate Lotka-Volterra model. 
    p = [α, β, γ, δ]
    predicted = solve(lv_prob, Tsit5(); p=p, saveat=0.1)

    # Observations.
    for i in 1:length(predicted)
        data[:, i] ~ MvNormal(predicted[i], σ^2 * I)
    end

    return nothing
end


# ╔═╡ 3fc31430-b15d-4fa8-a326-5de13f380b5c
lv_model = fitlv(odedata, lv_prob)

# ╔═╡ 5f62cec0-3e8f-42f0-9cf6-e30239456803
lv_chain = sample(MersenneTwister(42), lv_model, NUTS(0.65), MCMCSerial(), 1000, 3; progress=false)

# ╔═╡ 6c6ce02e-1413-49f4-925b-0d07711690c0
summarize(lv_chain)

# ╔═╡ daae0473-49b2-4630-aa79-79e10dcb3120
@model function fitlv2(data::AbstractVector, prob)
    # Prior distributions.
    σ ~ InverseGamma(2, 3)
    α ~ truncated(Normal(1.5, 0.5); lower=0.5, upper=2.5)
    β ~ truncated(Normal(1.2, 0.5); lower=0, upper=2)
    γ ~ truncated(Normal(3.0, 0.5); lower=1, upper=4)
    δ ~ truncated(Normal(1.0, 0.5); lower=0, upper=2)

    # Simulate Lotka-Volterra model but save only the second state of the system (predators).
    p = [α, β, γ, δ]
    predicted = solve(prob, Tsit5(); p=p, saveat=0.1, save_idxs=2)

    # Observations of the predators.
    data ~ MvNormal(predicted.u, σ^2 * I)

    return nothing
end

# ╔═╡ 27ddc0aa-398c-439c-9a5c-a6153d081f51
model2 = fitlv2(odedata[2, :], lv_prob)

# ╔═╡ c7a19043-cf44-44cc-8317-5caff29269cc
count_data = rand.(Poisson.(Array(sol)))

# ╔═╡ 08da7b32-1e53-46da-bfeb-ea16b23dddf8
@model function fitlv3(data, prob)
    # Prior distributions.
    σ ~ InverseGamma(2, 3)
    α ~ truncated(Normal(1.5, 0.5); lower=0.5, upper=2.5)
    β ~ truncated(Normal(1.2, 0.5); lower=0, upper=2)
    γ ~ truncated(Normal(3.0, 0.5); lower=1, upper=4)
    δ ~ truncated(Normal(1.0, 0.5); lower=0, upper=2)

    # Simulate Lotka-Volterra model. 
    p = [α, β, γ, δ]
    predicted = solve(prob, Tsit5(); p=p, saveat=0.1)

    # Observations from a Poisson distribution
	n, m = size(data)
    for i in 1:n
		for j in 1:m
        	data[i,j] ~ Poisson(predicted[i,j])
		end
    end

    return nothing
end

# ╔═╡ 570db1f7-324c-41b8-b9eb-eea747ca8b9e
lv_model3 = fitlv3(count_data, lv_prob)

# ╔═╡ ec0dd911-b3e1-43c8-9fdc-8fffe49945d3
md"## Logistic model"

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
# ╠═╡ disabled = true
# ╠═╡ skip_as_script = true
#=╠═╡
logistic(t; r, K, y₀=.1) = K / (1 + (K-y₀) / y₀ * exp(-r*t))
  ╠═╡ =#

# ╔═╡ 6984bd57-e95f-40b1-aecc-a864177bdb72
#=╠═╡
sq_loss(r, K) = sum(abs2, logistic.(t; r, K) .- y) |> sqrt
  ╠═╡ =#

# ╔═╡ e4fd9a2f-92be-4d4f-a3b0-475848e7aed3
#=╠═╡
contourf(0.01:.01:2, 1:.01:20, sq_loss, color=:speed)
  ╠═╡ =#

# ╔═╡ d43ee840-dfc5-4d2d-bb0f-14faf1573e73
# ╠═╡ disabled = true
#=╠═╡
@model function logistic_model(t, y; y₀=missing)
	σ ~ InverseGamma()
	r ~ Exponential(.1)
	K ~ Exponential(5)
	y₀ ~ Exponential(1)
	for i in 1:length(y)
		y[i] ~ Normal(logistic(t[i]; r, K, y₀), σ)
	end
end
  ╠═╡ =#

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
# ╠═╡ disabled = true
#=╠═╡
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
  ╠═╡ =#

# ╔═╡ 9b3e76fc-0c0b-4bd0-8af8-5e6053e24094
# ╠═╡ disabled = true
#=╠═╡
xp = 10rand(200) .- 5 |> sort!
  ╠═╡ =#

# ╔═╡ 960814cb-611e-4057-8667-0c4de798a96e
# ╠═╡ disabled = true
#=╠═╡
yp = 3randn(200) .+ 4 .- 4xp .+ 8xp.^2
  ╠═╡ =#

# ╔═╡ 9b317909-8997-46e2-809d-b005332c1022
# ╠═╡ disabled = true
#=╠═╡
sample(polynomial_regression(xp, yp, 5), NUTS(), 1000) |> summarize
  ╠═╡ =#

# ╔═╡ bf80bb42-75c8-4064-952c-0613e3e8188d
# ╠═╡ disabled = true
#=╠═╡
optimize(polynomial_regression(xp, yp, 2), MLE(), NelderMead())
  ╠═╡ =#

# ╔═╡ 5ff4ccf6-9bd8-4318-9de3-7752ed66f51a
# ╠═╡ disabled = true
#=╠═╡
optimize(polynomial_regression(xp, yp, 5), MAP()) |> coeftable
  ╠═╡ =#

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
# ╠═╡ skip_as_script = true
#=╠═╡
optimize(yeast_mod, MLE(), LBFGS()) |> coeftable
  ╠═╡ =#

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
rng = MersenneTwister(6)

# ╔═╡ 2d82ae12-da64-426c-9ed1-f65480b88214
μ_noise = μs .+ rand(rng, Normal(0, σmm), n_mm)

# ╔═╡ f3fd99c7-94db-4cfa-a7f2-fe8d2bf33f8e
mm_model_normal = michaelis_menten(Xs, μ_noise)

# ╔═╡ a8de733e-8f92-443d-9a79-6dc4d17ca52b
plot(generated_quantities(mm_model_normal, rand(mm_model_normal)))

# ╔═╡ e2f28fc3-5cfd-4e5f-8e28-36a0e8bf43a2
# ╠═╡ skip_as_script = true
#=╠═╡
coeftable(optimize(mm_model_normal, MLE(), NelderMead()))
  ╠═╡ =#

# ╔═╡ b1e883f7-10bb-47b6-af26-0c278b15c0b6
chain_MM = sample(mm_model_normal, NUTS(), 10000)

# ╔═╡ 9a51f9f4-b403-4466-981f-1aaa74a48e5a
summarize(chain_MM)

# ╔═╡ 5ba1377a-b16e-4e70-8b72-c1477a75f4e8
mm_map = optimize(mm_model_normal, MAP(), NelderMead())

# ╔═╡ 0d262055-5220-42ae-b232-d96d4965667b
# ╠═╡ skip_as_script = true
#=╠═╡
coeftable(mm_map)
  ╠═╡ =#

# ╔═╡ ab9d0e48-6b7c-425b-a740-63bf65623a9a
Σ = informationmatrix(mm_map, expected=false) |> inv

# ╔═╡ ebb6d109-d4d2-4997-bf2b-6b7c12d47e3b
mu_max_map, Ks_map = mm_map.values[:mu_max], mm_map.values[:Ks]

# ╔═╡ da7b2976-fa09-405a-8c86-081984c0ac8c
function sq_loss_mm(μmax, Ks)
	L = 0.0
	for (X, μ) in zip(Xs, μ_noise)
		L += (mm(X, μmax, Ks) - μ)^2
	end
	return L
end

# ╔═╡ 1fda5580-0dcf-47f9-b03d-7e588d0e9c80
ls_sol = optimize(θ->sq_loss_mm(θ[1], θ[2]), [1.0, 10.0], NelderMead())

# ╔═╡ 6fdabf8f-3569-4b5e-a7d8-79c8cdd2fcfc
mu_max_ls, Ks_ls = ls_sol.minimizer

# ╔═╡ 794fec60-c448-4505-8190-b60e0119200f
begin
	μ_outliers = μs .+ rand(rng, Laplace(0, σmm), n_mm)
	μ_outliers[[ 5, 17, 15, 23, 27, 30]] .= 0
	μ_outliers
end

# ╔═╡ 7e8586ee-0d49-4455-8ad8-82a3994a52b4
scatter(Cs, μ_outliers)

# ╔═╡ 505b7f45-88af-4909-98e5-7f9311fafa01
xs = rand(rng, Poisson(λ_pois), 10)

# ╔═╡ 19221533-56dd-4d87-a00f-1b8f47160b00
ll_pois(λ) = sum(x->loglikelihood(Poisson(λ), x), xs)

# ╔═╡ 25fc4ca7-9e84-48f2-a4d4-4ebc0ad7183a
xs_large = rand(rng, Poisson(λ_pois), 50)

# ╔═╡ 37aac79a-ec46-4b6f-a9b8-abf16a5fe7b7
ll_pois_large(λ) = sum(x->loglikelihood(Poisson(λ), x), xs_large)

# ╔═╡ a882ce40-3cba-464e-9069-62731120fccc
chain2 = sample(rng, model2, NUTS(0.45), MCMCSerial(), 5000, 3; progress=false);

# ╔═╡ de189153-7f26-4225-9cc1-ddbbcfac7867
summarize(chain2)

# ╔═╡ 13372de2-77e8-40d9-b76a-7f379d0587e1
lv_chain3 = sample(rng, lv_model3, NUTS(), MCMCSerial(), 1000, 3; progress=false)

# ╔═╡ 0cbe46a7-790a-4948-9ca0-5d84575a3f72
neg(x) = -x

# ╔═╡ 55459b43-4e53-4f44-80d7-84c1f929fa52
plots = Dict()

# ╔═╡ 6edf2f34-1e9e-4775-bd6a-624d003b6688
let
	plots["poly_data"] = scatter(xpoly, ypoly)
end

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

# ╔═╡ 83bbcd7b-68e8-4f75-a007-400dd427913f
plots["MM_norm_data"] = scatter(Cs, μ_noise, xlab="X [mmol/L]", ylab="μ [mmol/s]", title="Michaelis-Menten data", label="data")

# ╔═╡ 2a154a97-2320-476b-8374-14cef9743c27
let
	p = contourf(0:0.1:50, 5:1:500, log ∘ sq_loss_mm, color=:speed, xlab=L"\mu_\max", ylab=L"K_s", title="Squared loss Michaelis-Menten")
	xlims!(0, 50)
	ylims!(5, 500)
	plots["MM_ls_loss"] = scatter!([mu_max_ls], [Ks_ls], label="minimizer")
end

# ╔═╡ 721c147c-1ecb-4898-a962-08abe2a764d4
let
	p = scatter(Cs, μ_noise, xlab="X [mmol/L]", ylab="μ [mmol/s]", label="data")
	title!(p, "MM least-squares")
	plot!(x->mm(x, mu_max_ls, Ks_ls), 0, 180, lw=2, label="LS fit")
	plots["MM_LS_fit"] = p
end

# ╔═╡ d5504f95-3fd4-4c02-9975-e201a9a09801
let
	ll(μmax, Ks) = -loglikelihood(mm_model_normal, (mu_max=μmax, Ks=Ks, sigmasq=1))
	p = contourf(0:0.1:50, 5:1:500, log ∘ ll, color=:speed, xlab=L"\mu_\max", ylab=L"K_s", title="Neg log-likelihood MM")
	xlims!(0, 50)
	ylims!(5, 500)
	plots["MM_ll"] = scatter!([mu_max_ls], [Ks_ls], label="minimizer")
end

# ╔═╡ 9cf24cc2-ac54-4239-a313-ee7ecd680f35
plots["MM_joint_post"] = marginalkde(chain_MM[:mu_max], chain_MM[:Ks], xlab=L"\mu_\max", ylab=L"K_s")

# ╔═╡ b69eb7d3-8660-44b2-992e-a2b9696d8e50
let
	p = scatter(Cs, μ_noise, xlab="X [mmol/L]", ylab="μ [mmol/s]", label="data")
	title!(p, "MM MAP + posterior")
	for mm_post in rand(generated_quantities(mm_model_normal, chain_MM), 200)
		plot!(mm_post, 0, 180, lw=0.1, color="#BBBBBB", alpha=0.5, label="")
	end
	plot!(x->mm(x, mu_max_map, Ks_map), 0, 180, lw=2, label="MAP fit", color=2)
	plot!([], lw=0.5, color="#BBBBBB", alpha=0.5, label="posterior sample")
	plots["MM_MAP_post"] = p
end

# ╔═╡ 12990770-b44f-44b7-b981-9a008d9028d1
plots["MM_outlier_data"] = scatter(Cs, μ_outliers, xlab="X [mmol/L]", ylab="μ [mmol/s]", title="Michaelis-Menten outliers", label="data")

# ╔═╡ b09cd3cb-c611-4679-ad6e-bc0eeddbf0f8
let
	mle = optimize(michaelis_menten(Xs, μ_outliers), MLE(), NelderMead())
	@show mu_max_n, Ks_n = mle.values[:mu_max], mle.values[:Ks]

	mle = optimize(michaelis_menten_robust(Xs, μ_outliers), MLE(), NelderMead())
	@show mu_max_r, Ks_r = mle.values[:mu_max], mle.values[:Ks]

	
	p = scatter(Cs, μ_outliers, xlab="X [mmol/L]", ylab="μ [mmol/s]", title="Michaelis-Menten outliers fit", label="data")
	plot!(x->mm(x, mu_max_n, Ks_n), 0, 180, lw=2, label="LS fit")
	plot!(x->mm(x, mu_max_r, Ks_r), 0, 180, lw=2, label="absolute error fit", ls=:dash)
	plots["MM_outliers_fit"] = p
end

# ╔═╡ d864db4d-e240-4b90-a8f7-0327eb9e18b0
let
	p = plot(ll_pois, 10, 50, label="log-likelihood 1 (small sample)", lw=2)
	plot!(ll_pois_large, 10, 50, label="log-likelihood 2 (large sample)", lw=2, ls=:dash)
	title!("Log-likelihoods different datasets")
	xlabel!(L"\lambda")
	ylabel!(L"\log L(\lambda)")
	plots["poisson_LL"] = p
end

# ╔═╡ 8b71bd3c-60b6-425c-87cf-64652d04c8b7
let
	p = plot(pois_map, 0.3, 8, lw=2, label="unnormalized log-posterior", legend=:bottom)
	plot!(x->pois_map(λ_star) - fi/2 * (λ_star-x)^2, 0, 8, label="second-order approx.", lw=2, ls=:dash)
	title!("Laplace approximation")
	xlabel!(L"\lambda")
	ylabel!(L"\log P(\lambda\mid D)")
	plots["Laplace_log"] = p
end

# ╔═╡ 6ed728a8-80a8-457a-8908-3c840b88efae
let
	# scaling factor posterior
	a = exp(logpdf(post_Laplace, λ_star) - pois_map(λ_star))
	p = plot(l -> a * exp(pois_map(l)), 0, 8, lw=2, label="posterior (scaled)", legend=:topright)
	plot!(x->pdf(post_Laplace, x), 0, 8, label="normal approx.", lw=2, ls=:dash)
	title!("Laplace approximation")
	xlabel!(L"\lambda")
	ylabel!(L"P(\lambda\mid D)")
	plots["Laplace"] = p
end

# ╔═╡ d30119bd-e036-45fa-a325-e480d7d334b9
plots["LV_diagnostic"] = plot(lv_chain)

# ╔═╡ 60022e21-b74a-41b7-a43d-6b14ed482112
let
	p = plot(; legend=false)
	posterior_samples = sample(lv_chain[[:α, :β, :γ, :δ]], 300; replace=false)
	for p in eachrow(Array(posterior_samples))
	    sol_p = solve(lv_prob, Tsit5(); p=p, saveat=0.1)
	    plot!(sol_p; alpha=0.1, color="#BBBBBB")
	end
	
	# Plot simulation and noisy observations.
	plot!(sol; color=[1 2], linewidth=2, ls=:auto)
	title!("Posteriors of Lotka-Volterra")
	scatter!(sol.t, odedata'; color=[1 2], marker=[:o :^])
	plots["LV_inf"] = p
end

# ╔═╡ 94cbbbd8-00f7-47fa-b814-0cbca085b886
let
	p = plot(; legend=false)
	posterior_samples = sample(chain2[[:α, :β, :γ, :δ]], 300; replace=false)
	for p in eachrow(Array(posterior_samples))
	    sol_p = solve(lv_prob, Tsit5(); p=p, saveat=0.1)
	    plot!(sol_p; alpha=0.1, color="#BBBBBB")
	end
	
	title!("Posteriors using only the predators")
		
	# Plot simulation and noisy observations.
	plot!(sol; color=[1 2], linewidth=2, ls=:auto)
	scatter!(sol.t, odedata'; color=[1 2], marker=[:o :^], alpha=[0.3 1])
	plots["LV_only_pred"] = p
end

# ╔═╡ b8f6255c-f0e7-49d0-a14e-2711df5c2419
plots["LV_counts_diag"] = plot(lv_chain3)

# ╔═╡ f845114c-b612-43c1-8a49-c3da8f52e278
let
	p = plot(; legend=false)
	posterior_samples = sample(lv_chain3[[:α, :β, :γ, :δ]], 300; replace=false)
	for p in eachrow(Array(posterior_samples))
	    sol_p = solve(lv_prob, Tsit5(); p=p, saveat=0.1)
	    plot!(sol_p; alpha=0.1, color="#BBBBBB")
	end
	
	# Plot simulation and noisy observations.
	plot!(sol; color=[1 2], linewidth=2, ls=:auto)
	yaxis!(0:maximum(count_data))
	title!("Posteriors of Lotka-Volterra (counts)")
	scatter!(sol.t, count_data'; color=[1 2], marker=[:o :^])
	plots["LV_counts"] = p
end

# ╔═╡ a2be0a77-fca0-4d7f-af15-254faaa7e336
plots

# ╔═╡ Cell order:
# ╠═ab31fa80-0606-11ef-06f5-bd74bc457ff5
# ╠═a10a9f8e-4088-49e2-a497-b1d253d83ecc
# ╠═4b38ba22-8c3f-4c90-b42a-79063489641e
# ╠═3e70b82a-e4d3-4747-9679-aa5ab41d5b06
# ╠═0b6b3b7f-98e7-4922-82de-6a98453a627c
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
# ╠═e735ea0e-0c98-4707-bc5e-9c2bba04aa90
# ╠═5b80b912-e625-44a4-9653-8306189718c4
# ╠═2f8746e7-2e4a-4e40-bd48-f69349f1259d
# ╠═17252072-b71f-4c30-b108-8163dfd63fbf
# ╠═bf757553-7650-4b42-a07e-41007a3c29db
# ╠═2d82ae12-da64-426c-9ed1-f65480b88214
# ╠═234bebb6-0d05-41a0-871d-50cf752d17a2
# ╠═4a336948-d46e-415f-bcaf-07c667da4089
# ╠═83bbcd7b-68e8-4f75-a007-400dd427913f
# ╠═f3fd99c7-94db-4cfa-a7f2-fe8d2bf33f8e
# ╠═a8de733e-8f92-443d-9a79-6dc4d17ca52b
# ╠═da7b2976-fa09-405a-8c86-081984c0ac8c
# ╠═1fda5580-0dcf-47f9-b03d-7e588d0e9c80
# ╠═6fdabf8f-3569-4b5e-a7d8-79c8cdd2fcfc
# ╠═2a154a97-2320-476b-8374-14cef9743c27
# ╠═721c147c-1ecb-4898-a962-08abe2a764d4
# ╠═d5504f95-3fd4-4c02-9975-e201a9a09801
# ╠═e2f28fc3-5cfd-4e5f-8e28-36a0e8bf43a2
# ╠═0d262055-5220-42ae-b232-d96d4965667b
# ╠═ab9d0e48-6b7c-425b-a740-63bf65623a9a
# ╠═b1e883f7-10bb-47b6-af26-0c278b15c0b6
# ╠═9cf24cc2-ac54-4239-a313-ee7ecd680f35
# ╠═5ba1377a-b16e-4e70-8b72-c1477a75f4e8
# ╠═ebb6d109-d4d2-4997-bf2b-6b7c12d47e3b
# ╠═b69eb7d3-8660-44b2-992e-a2b9696d8e50
# ╠═9a51f9f4-b403-4466-981f-1aaa74a48e5a
# ╠═794fec60-c448-4505-8190-b60e0119200f
# ╠═12990770-b44f-44b7-b981-9a008d9028d1
# ╠═07d91dc7-2db0-464e-bc61-9633f2bfe92b
# ╠═b09cd3cb-c611-4679-ad6e-bc0eeddbf0f8
# ╠═7e8586ee-0d49-4455-8ad8-82a3994a52b4
# ╠═0e5831ff-f542-40e1-9376-2abb420494b7
# ╠═03373997-a9eb-4d5a-9ef5-8822c3f95285
# ╠═70f8de3c-5c5f-4ff0-92eb-e6d2a06451e1
# ╠═505b7f45-88af-4909-98e5-7f9311fafa01
# ╠═25fc4ca7-9e84-48f2-a4d4-4ebc0ad7183a
# ╠═cd3746d4-518d-476e-abf4-2ac92a2c1387
# ╠═19221533-56dd-4d87-a00f-1b8f47160b00
# ╠═37aac79a-ec46-4b6f-a9b8-abf16a5fe7b7
# ╠═d864db4d-e240-4b90-a8f7-0327eb9e18b0
# ╠═6dc06148-79bf-44bf-8a5c-b6221b55bce0
# ╠═d733961a-094a-4724-9b16-ee669a559945
# ╠═46551486-d8a3-4aca-a92d-8c7f1c4439c0
# ╠═0d6cc24c-6008-4bb7-8fb9-45223febc05b
# ╟─1e634f77-286b-44ce-b768-72b060ee898b
# ╠═0082f110-22df-49d5-8958-29e9516e9ca8
# ╠═362fe960-f247-4d7b-ab9d-01b67d1e647e
# ╠═8feddb40-5cf8-4d76-a959-dd3071c20e27
# ╠═8b71bd3c-60b6-425c-87cf-64652d04c8b7
# ╠═6ed728a8-80a8-457a-8908-3c840b88efae
# ╠═da1aa2a7-e876-4b80-9204-d3ea11bc0f44
# ╠═73e018a3-27f6-4631-a8c3-394644cb774e
# ╠═e79a72fc-9929-4faa-a76d-ad914fc731ab
# ╠═4c9f01b3-aef4-47ed-9633-1abda691b531
# ╠═c3008b5a-7138-48a5-b744-d8f3aa2ff2ad
# ╠═db3f89f9-3bef-432f-b558-3d32c5e79abe
# ╠═eb132faa-fedb-4e18-868f-46783e503515
# ╠═2a25b088-7414-48a2-a702-8c87213b92d3
# ╠═a05d6800-a168-4e44-bfcb-0809aedda09e
# ╠═b3c31cca-b472-4dec-b700-11cbff9465a9
# ╠═14810fe9-e6fe-46c7-b6cf-9705b5974113
# ╠═79473eed-7fba-43e7-80a7-51d011931578
# ╠═3fc31430-b15d-4fa8-a326-5de13f380b5c
# ╠═5f62cec0-3e8f-42f0-9cf6-e30239456803
# ╠═6c6ce02e-1413-49f4-925b-0d07711690c0
# ╠═d30119bd-e036-45fa-a325-e480d7d334b9
# ╠═60022e21-b74a-41b7-a43d-6b14ed482112
# ╠═daae0473-49b2-4630-aa79-79e10dcb3120
# ╠═27ddc0aa-398c-439c-9a5c-a6153d081f51
# ╠═a882ce40-3cba-464e-9069-62731120fccc
# ╠═de189153-7f26-4225-9cc1-ddbbcfac7867
# ╠═94cbbbd8-00f7-47fa-b814-0cbca085b886
# ╠═c7a19043-cf44-44cc-8317-5caff29269cc
# ╠═08da7b32-1e53-46da-bfeb-ea16b23dddf8
# ╠═570db1f7-324c-41b8-b9eb-eea747ca8b9e
# ╠═13372de2-77e8-40d9-b76a-7f379d0587e1
# ╠═b8f6255c-f0e7-49d0-a14e-2711df5c2419
# ╠═f845114c-b612-43c1-8a49-c3da8f52e278
# ╠═ec0dd911-b3e1-43c8-9fdc-8fffe49945d3
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
# ╠═0cbe46a7-790a-4948-9ca0-5d84575a3f72
# ╠═55459b43-4e53-4f44-80d7-84c1f929fa52
# ╠═a2be0a77-fca0-4d7f-af15-254faaa7e336
