### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 7e12cdbf-baf4-4ede-b6db-628544b1f211
begin
    using Pkg
	Pkg.activate("..")
	using Catalyst, OrdinaryDiffEq
end

# ╔═╡ 2fff8bf1-ee6d-4da2-b077-fe6075d580cd
infection = @reaction_network begin
	α * β, S + I --> 2S
	r * m, I --> D
	r * (1 - m), I --> R
end

# ╔═╡ e6e6034b-eeb8-4bde-a050-b4a0c4de2289
convert(ODESystem, infection)

# ╔═╡ 45bd99af-d739-4bf2-abbf-2a2731a0fb09
plankton = @reaction_network begin
	@parameters a b Y
	r * (1 - (F+Z)/K) * F, a * N + F => 2F
	α * F * Z, F + Z + b * N => Y * F  # to fix
	d, Z --> 0
	Nin * Q / V, 0 --> N
	Q / V, N --> 0
end

# ╔═╡ 5eef9443-b86a-49e8-98af-7ab453ef78fc
convert(ODESystem, plankton, combinatoric_ratelaws = false)

# ╔═╡ c8dc75f4-11d0-4660-bde7-062dc4afa74a
temperature = @reaction_network begin
	Qin / V1 / cp * (Tin - T1), 0 =>T1
	A1 * λ * (Ta - T1), T1 => 0 
	A2 * λ * (Ta - T2), T2 => 0 
	Qin / V2 / cp * (T1 - T2), 0 =>T2
	Ph / Cp / V2, 0 => T2
end

# ╔═╡ c54cd3ae-f97d-4e95-8692-d18d03883cc0
convert(ODESystem, temperature)

# ╔═╡ 8375fa94-7d78-4a59-a8a2-a1d7fa112fb6
crabs = @reaction_network begin
	n1, 0 --> C1
	0.01, (S1, S2) --> 0
	0.1, (C1, C2) --> (S1, S2)
	1/4, (S1, S2) --> (C2, C3)
	log(2), S3 --> 0  # harvesting
	c, S1 + C1 --> C1  # how to do cannibalism?
end

# ╔═╡ a140d269-7f21-4f69-a427-997e52fef43b
convert(ODESystem, crabs)

# ╔═╡ b1bfaf3d-4696-407d-9a63-2f2fc1a94ac5
yeast = @reaction_network begin
	@parameters k1 k2 k3 k4 k5 k6 k7 k8 k9 k10
	μ0 * X, k1 * G + k3 * O2 => X + k5 * CO2
	μr * X, k2 * G => X + k8 * E + k6 * CO2
	μe * X, k9 * E + k4 * O2 => X + k7 * CO2
	qm * X, k10 * G + O2 => k11 * CO2
	Gin * Qin, 0 --> G
	CO2out * V, CO2 --> 0
	O2in * V, 0 --> O2
end

# ╔═╡ ae4dd5f7-afc3-4b8c-a3bb-80b5dc451270
convert(ODESystem, yeast)

# ╔═╡ aef132cb-e305-4f3e-88a6-227daeb976b7


# ╔═╡ Cell order:
# ╠═7e12cdbf-baf4-4ede-b6db-628544b1f211
# ╠═2fff8bf1-ee6d-4da2-b077-fe6075d580cd
# ╠═e6e6034b-eeb8-4bde-a050-b4a0c4de2289
# ╠═45bd99af-d739-4bf2-abbf-2a2731a0fb09
# ╠═5eef9443-b86a-49e8-98af-7ab453ef78fc
# ╠═c8dc75f4-11d0-4660-bde7-062dc4afa74a
# ╠═c54cd3ae-f97d-4e95-8692-d18d03883cc0
# ╠═8375fa94-7d78-4a59-a8a2-a1d7fa112fb6
# ╠═a140d269-7f21-4f69-a427-997e52fef43b
# ╠═b1bfaf3d-4696-407d-9a63-2f2fc1a94ac5
# ╠═ae4dd5f7-afc3-4b8c-a3bb-80b5dc451270
# ╠═aef132cb-e305-4f3e-88a6-227daeb976b7
