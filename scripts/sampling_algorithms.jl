### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 89971672-d3d8-11ee-0f21-13c70f67d6f2
begin
    using Pkg
    Pkg.activate("..")
    using Plots, PlutoUI, LaTeXStrings
	using Distributions, LinearAlgebra
end

# ╔═╡ 9d4d1812-0303-4fa7-a7e8-a10549b61fe6
using Zygote

# ╔═╡ 4ced4f92-1832-4c81-a7a6-85f87aab5375
#f(x, y; a=1, b=100) = (a - x)^2 + b * (y - x^2)^2
f(x, y; R=2) = exp(-abs((x^2 + y^2) - R^2))

# ╔═╡ 28475353-415b-49cf-9a1a-dd349f297b8c
p = MixtureModel([Normal(2, 0.4), Normal(4, 1.2)], [0.3, 0.7])

# ╔═╡ c1ebb67d-8cce-4a12-b13a-fcb823729993
plot(x->pdf(p, x), 0, 8)

# ╔═╡ da7fea34-4200-4ea2-a054-6f95d2313d47
f(xy) = f(xy...)

# ╔═╡ 0cbe1a45-7090-4a54-9873-8a8ef1c43dcd
p_donut = contourf(-4:0.01:4, -4:0.01:4, f, color=:speed, aspect_ratio=:equal
	, xlims=(-4, 4), ylims=(-4, 4), xlab=L"x", ylab=L"y")

# ╔═╡ df31ee7d-c8d8-4725-a1c2-ecf4bca8d61d
G = MultivariateNormal([0, 0], I)

# ╔═╡ c27bd239-25ef-428b-b267-773a025cb83b
M = 1.1f(0, 2) / pdf(G, [0, 2])

# ╔═╡ efa7ac57-04c9-4845-995e-f1798bf77aa0
function rejection_sampling(f, G; M, n=100)
	y = rand(G)
	samples = typeof(y)[]
	rejected = typeof(y)[]
	while length(samples) < n
		y = rand(G)
		α = f(y) / (M * pdf(G, y))
		if rand() < α
			push!(samples, y)
		else
			push!(rejected, y)
		end
	end
	return samples, rejected
end

# ╔═╡ c4710ca9-32a4-4afe-993c-9957af1c4956
samples, rejected = rejection_sampling(f, G; M, n=100)

# ╔═╡ 98e65ea8-a9e4-458f-8ce7-661fe0b71ea0
begin
	scatter!(deepcopy(p_donut), first.(samples), last.(samples), label="accepted", ms=2)
	scatter!(first.(rejected), last.(rejected), label="rejected", ms=2, alpha=0.2)
end

# ╔═╡ 7f979860-34de-40bf-8d8c-9fb937692185
acceptence_ratio = length(samples) / (length(samples) + length(rejected))

# ╔═╡ 3fbfbb95-223b-4ff3-a981-37cc69cd1f98
1 / M

# ╔═╡ ace4ea8d-9669-4aa4-a3fd-eb126697e7d9
function metropolings_hastings(f, g, x₀; n=100)
	samples = typeof(x₀)[]
	accepted = Bool[]
	xₜ = x₀
	while count(accepted) < n
		G = g(xₜ)
		x′ = rand(G)
		push!(samples, x′)
		α = f(x′) / f(xₜ)		
		if rand() ≤ α
			push!(accepted, true)
			xₜ = x′
		else
			push!(accepted, false)
		end
	end
	return samples, accepted
end


# ╔═╡ ed5e7bb7-1a0f-4c89-9b88-d9a960400bda
samples_MH, accepted_MH = metropolings_hastings(f, x->MultivariateNormal(x, 0.3I), randn(2); n=100)

# ╔═╡ d51d6929-a6c6-4b31-ba42-201a783a1869
begin
	scatter!(deepcopy(p_donut), first.(samples_MH[accepted_MH]), last.(samples_MH[accepted_MH]), label="accepted", ms=2)
	scatter!(first.(samples_MH[.!accepted_MH]), last.(samples_MH[.!accepted_MH]), label="rejected", ms=2, alpha=0.2)
	plot!(first.(samples_MH), last.(samples_MH), label="", alpha=0.5, lw=0.5, color="grey")
end

# ╔═╡ 429e5106-76f2-4318-ae8b-63fcce7f685a
function leapfrog(xₜ, pₜ, ϵ, ∇U)
	phalf = pₜ .- ϵ / 2 .* ∇U(xₜ)
	xₜ = xₜ .+ ϵ .* phalf
	pₜ = phalf .- ϵ /2 .* ∇U(xₜ)
	return xₜ, pₜ
end

# ╔═╡ eea393d0-5c87-4911-8063-2a1ce27f37b6
function hamiltonian_monte_carlo(f, x₀; σₚ=1, ϵ, L, n=100)
	U(x) = -log(f(x))
	∇U(x) = U'(x)
	samples = typeof(x₀)[]
	pₜ = randn(length(x₀)) .* σₚ
	xₜ = x₀
	for i in 1:n
		for _ in 1:L
			xₜ, pₜ = leapfrog(xₜ, pₜ, ϵ, ∇U)
		end
		pₜ .= randn(length(x₀)) .* σₚ
		push!(samples, xₜ)
	end
	return samples
end

# ╔═╡ 7a2dcc5c-b7b7-4e09-88b5-b38619a3e249
samples_hmc = hamiltonian_monte_carlo(f, randn(2), σₚ=1, ϵ=0.1, L=5, n=100)

# ╔═╡ d5450495-3d9a-46b1-a53d-c0a3d69a6f5f
begin
	scatter!(deepcopy(p_donut), first.(samples_hmc), last.(samples_hmc), label="HMC", ms=2)
	
	plot!(first.(samples_hmc), last.(samples_hmc), label="", alpha=0.5, lw=0.5, color="grey")
end

# ╔═╡ Cell order:
# ╠═89971672-d3d8-11ee-0f21-13c70f67d6f2
# ╠═4ced4f92-1832-4c81-a7a6-85f87aab5375
# ╠═28475353-415b-49cf-9a1a-dd349f297b8c
# ╠═c1ebb67d-8cce-4a12-b13a-fcb823729993
# ╠═da7fea34-4200-4ea2-a054-6f95d2313d47
# ╠═0cbe1a45-7090-4a54-9873-8a8ef1c43dcd
# ╠═df31ee7d-c8d8-4725-a1c2-ecf4bca8d61d
# ╠═c27bd239-25ef-428b-b267-773a025cb83b
# ╠═efa7ac57-04c9-4845-995e-f1798bf77aa0
# ╠═c4710ca9-32a4-4afe-993c-9957af1c4956
# ╠═98e65ea8-a9e4-458f-8ce7-661fe0b71ea0
# ╠═7f979860-34de-40bf-8d8c-9fb937692185
# ╠═3fbfbb95-223b-4ff3-a981-37cc69cd1f98
# ╠═ace4ea8d-9669-4aa4-a3fd-eb126697e7d9
# ╠═ed5e7bb7-1a0f-4c89-9b88-d9a960400bda
# ╠═d51d6929-a6c6-4b31-ba42-201a783a1869
# ╠═9d4d1812-0303-4fa7-a7e8-a10549b61fe6
# ╠═eea393d0-5c87-4911-8063-2a1ce27f37b6
# ╠═429e5106-76f2-4318-ae8b-63fcce7f685a
# ╠═7a2dcc5c-b7b7-4e09-88b5-b38619a3e249
# ╠═d5450495-3d9a-46b1-a53d-c0a3d69a6f5f
