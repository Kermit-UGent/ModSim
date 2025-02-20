### A Pluto.jl notebook ###
# v0.20.4

#> [frontmatter]
#> order = "2"
#> title = "template (EN)"
#> date = "2025-02-07"
#> tags = ["project"]
#> description = "Project Template"
#> layout = "layout.jlhtml"
#> 
#>     [[frontmatter.author]]
#>     name = "Michiel Stock"

using Markdown
using InteractiveUtils

# ╔═╡ 21357c48-f35d-11ee-23f8-2534bb1d82f4
begin

    using Pkg
    Pkg.activate("../../pluto-deployment-environment")
	
	# make this cell invisible when you are finished
	title = "Our super cool project"
	names = ["Alice", "Bob", "Carol"]

	academic_year = "202x_202(x+1)"

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

# ╔═╡ aed58771-86c1-4928-a1ba-b7d2f503b188
# metabolic rate
y(m; a=0.75, C0=1) = C0 * m^a

# ╔═╡ 98451eb8-8755-47fe-b8c6-ce2959256dd8
# note that you likely will use a Catalyst model

# ╔═╡ 797ce7d7-c49d-4a00-a8e3-058b6de4d3a9


# ╔═╡ aa32c5e2-0f67-4910-a612-d958fa13169c


# ╔═╡ 7f9c8c7d-90ef-496c-b907-300feba5edfc


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

# ╔═╡ 23970226-3abe-4c37-8de6-b62c4acef55e
md"""
**Attribution**: one sentence about who did what mainly according to the [CRediT](https://en.wikipedia.org/wiki/Contributor_Roles_Taxonomy) (Contribution Roles Taxonomy) classification.
"""

# ╔═╡ 890afefc-f42b-4d74-b775-6dee5e5f0c2b
md"## Appendix"

# ╔═╡ Cell order:
# ╟─14c7e803-c0ff-4211-8b60-2c9c246934dd
# ╠═6948149f-854e-4b3c-b51c-099dd221ab83
# ╠═89551690-500d-4e37-ae20-5beb71cc87ac
# ╠═aed58771-86c1-4928-a1ba-b7d2f503b188
# ╠═98451eb8-8755-47fe-b8c6-ce2959256dd8
# ╠═797ce7d7-c49d-4a00-a8e3-058b6de4d3a9
# ╠═aa32c5e2-0f67-4910-a612-d958fa13169c
# ╠═7f9c8c7d-90ef-496c-b907-300feba5edfc
# ╟─f284fda5-aaf7-4121-a8f0-b996d501bec5
# ╠═2a0b0c1f-9510-4f71-9d65-b4fb7c854a98
# ╠═67380c07-a07c-4a0e-8654-80f27b951461
# ╠═86307aaa-1349-444b-bc9b-e5c115727671
# ╠═ef1b3843-a963-4a18-85c4-ca5a6e791d00
# ╠═aa4ebdb8-27a1-493b-b95e-d2ac4c3e6d54
# ╠═60c63f1b-5a27-4539-a486-78c083457b0e
# ╠═23970226-3abe-4c37-8de6-b62c4acef55e
# ╟─890afefc-f42b-4d74-b775-6dee5e5f0c2b
# ╠═21357c48-f35d-11ee-23f8-2534bb1d82f4
