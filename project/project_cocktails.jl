### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 1c4709fc-58a6-4a2c-857c-49a5442a180f
using DifferentialEquations, ModelingToolkit, Unitful

# ╔═╡ 21357c48-f35d-11ee-23f8-2534bb1d82f4
begin
	
	# make this cell invisible when you are finished
	title = "Cocktail shaking model"
	names = ["Michiel"]

	academic_year = "2023_20224"

	email_main_person = "mail@domain.be"

	using PlutoUI  # interactivity
	using Plots  # plotting
	TableOfContents()
end;

# ╔═╡ 14c7e803-c0ff-4211-8b60-2c9c246934dd
md"""
# $title

**$(join(names, ", ", " and "))**
"""

# ╔═╡ 6948149f-854e-4b3c-b51c-099dd221ab83
md"""
## Abstract

About 250 words about your project:
- (1-2 sentence) basic introduction to your topic, accessible to every bioengineering student
- (1-2 sentences) bit more specialized introduction
- (1-2 sentences) general goal of the project
- (2-3 sentences) short overview of how you built the model and what analysis you did
"""

# ╔═╡ 89551690-500d-4e37-ae20-5beb71cc87ac
md"""
## Model

general outline of the model + variables + parameters

For example, the metabolic rate $y$ as a function of the mass $m$ of an organism follows a power law.
"""

# ╔═╡ 7b1a87c5-0856-4968-ba2b-36da650cd0c8
begin
	@variables t #[unit = u"s"]   # mixing time in seconds
	@variables I(t)=200 #[unit = u"g"]  # amount of ice
	@variables L(t)=100 #[unit = u"L"]  # amount of liquid
	@variables T(t)=20 #[unit = u"°C"]  # cocktail temperature
	@variables S(t)=203 #[unit = u"g/L"]  # sugar concentration
	@variables Z(t)  
end;

# ╔═╡ 08b40947-7218-4885-b25c-982821ca0361
begin
	@constants Cmelt = 3.34e5 #[unit = u"J/kg"]
	@constants Cₚ = 4.187e3 #[unit = u"J/kg * K"]
end;

# ╔═╡ b1073c63-cb31-4160-98e1-ec1a0ca22d9c
10u"J/kg/K" * (10u"°C"-5u"°C")

# ╔═╡ 7110e56a-dfc1-4885-9eee-c59eeb37010f
Dₜ = Differential(t)

# ╔═╡ d560e96d-c6dc-4127-9616-1c186914cef5
melting = Dₜ(I) ~ -0.001 #* T #*I^0.66

# ╔═╡ 75b7ed6d-f27f-4f73-9c78-014b659323b8
#water_balance = Dₜ(I) + Dₜ(L) ~ 0 
water_balance = I+L ~ 100 

# ╔═╡ 3da54360-1f31-48a0-87cf-0d47132862bc
heat_balance = Cmelt * Dₜ(I) ~ Cₚ * L * Dₜ(T)

# ╔═╡ 3ccb5e2b-9ed3-486c-ac18-de62d2297e07
sugar_balance = expand_derivatives(Dₜ(S*L)) ~ 0
#sugar_balance = S*L ~ 10

# ╔═╡ aed58771-86c1-4928-a1ba-b7d2f503b188
md"Finally, the perceived sweetness depends on the sugar concentration and the temperature. Setting up such a relation is done in the field of psychophysics (e.g. [Steven's law](https://en.wikipedia.org/wiki/Stevens%27s_power_law)). Based on [this article](https://pubmed.ncbi.nlm.nih.gov/7100291/), we obtain a fairly simple emperical relation."

# ╔═╡ 98451eb8-8755-47fe-b8c6-ce2959256dd8
sweetness = Z ~ 14.0 * (S / 342.30)^(1.422 - 0.0146T)

# ╔═╡ aa32c5e2-0f67-4910-a612-d958fa13169c
@named sys_or = ODESystem([melting, water_balance,heat_balance,sugar_balance,sweetness])

# ╔═╡ 7f9c8c7d-90ef-496c-b907-300feba5edfc
simpsys = structural_simplify(sys, simplify=true)

# ╔═╡ 4595b9f6-3e9b-4160-8def-456b0c6c0caf
@mtkbuild sys = ODESystem([melting, water_balance,heat_balance,sugar_balance,sweetness])

# ╔═╡ 58d13e42-4edc-412c-a1e0-30e11fb2585e
states(simpsys)

# ╔═╡ 8cf34815-4ff7-4930-bc91-11081e1ab3f3
ODEProblem(simpsys, [S=>10, L=>0.3, S=>0.2, T=>21, I=>200], (0, 50))

# ╔═╡ f284fda5-aaf7-4121-a8f0-b996d501bec5
md"## Simulation and analysis"

# ╔═╡ 2a0b0c1f-9510-4f71-9d65-b4fb7c854a98
md"""
Explore your model
"""

# ╔═╡ 67380c07-a07c-4a0e-8654-80f27b951461
plot(y, 0.01, 1000, label="metabolic rate", xlab="mass (kg)")

# ╔═╡ 86307aaa-1349-444b-bc9b-e5c115727671


# ╔═╡ ef1b3843-a963-4a18-85c4-ca5a6e791d00


# ╔═╡ aa4ebdb8-27a1-493b-b95e-d2ac4c3e6d54


# ╔═╡ 60c63f1b-5a27-4539-a486-78c083457b0e
md"""
## Conclusion

A short conclusion of your analysis with a relection on how you would improve this model.
"""

# ╔═╡ 890afefc-f42b-4d74-b775-6dee5e5f0c2b
md"## Appendix"

# ╔═╡ 1f4b988a-6632-4589-8913-5de0574d94b3
perceived_sweetness(S, T) = 14.9 * (S / 342.30)^(1.42224 - 0.0146071T)

# ╔═╡ 797ce7d7-c49d-4a00-a8e3-058b6de4d3a9
heatmap(1:200, 4:40, perceived_sweetness, xlabel="Sugar concentration (g/L)", ylabel="temperature (°C)", title="perceived sweetness")

# ╔═╡ Cell order:
# ╟─14c7e803-c0ff-4211-8b60-2c9c246934dd
# ╠═6948149f-854e-4b3c-b51c-099dd221ab83
# ╠═1c4709fc-58a6-4a2c-857c-49a5442a180f
# ╠═89551690-500d-4e37-ae20-5beb71cc87ac
# ╠═7b1a87c5-0856-4968-ba2b-36da650cd0c8
# ╠═08b40947-7218-4885-b25c-982821ca0361
# ╠═b1073c63-cb31-4160-98e1-ec1a0ca22d9c
# ╠═7110e56a-dfc1-4885-9eee-c59eeb37010f
# ╠═d560e96d-c6dc-4127-9616-1c186914cef5
# ╠═75b7ed6d-f27f-4f73-9c78-014b659323b8
# ╠═3da54360-1f31-48a0-87cf-0d47132862bc
# ╠═3ccb5e2b-9ed3-486c-ac18-de62d2297e07
# ╠═aed58771-86c1-4928-a1ba-b7d2f503b188
# ╠═98451eb8-8755-47fe-b8c6-ce2959256dd8
# ╟─797ce7d7-c49d-4a00-a8e3-058b6de4d3a9
# ╠═aa32c5e2-0f67-4910-a612-d958fa13169c
# ╠═7f9c8c7d-90ef-496c-b907-300feba5edfc
# ╠═4595b9f6-3e9b-4160-8def-456b0c6c0caf
# ╠═58d13e42-4edc-412c-a1e0-30e11fb2585e
# ╠═8cf34815-4ff7-4930-bc91-11081e1ab3f3
# ╟─f284fda5-aaf7-4121-a8f0-b996d501bec5
# ╠═2a0b0c1f-9510-4f71-9d65-b4fb7c854a98
# ╠═67380c07-a07c-4a0e-8654-80f27b951461
# ╠═86307aaa-1349-444b-bc9b-e5c115727671
# ╠═ef1b3843-a963-4a18-85c4-ca5a6e791d00
# ╠═aa4ebdb8-27a1-493b-b95e-d2ac4c3e6d54
# ╠═60c63f1b-5a27-4539-a486-78c083457b0e
# ╟─890afefc-f42b-4d74-b775-6dee5e5f0c2b
# ╠═21357c48-f35d-11ee-23f8-2534bb1d82f4
# ╟─1f4b988a-6632-4589-8913-5de0574d94b3
