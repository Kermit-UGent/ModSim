### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 57420478-d744-11ef-2ec5-735d47866238
import Pkg; Pkg.activate("..")

# ╔═╡ 6261bb1a-5f73-4c9d-8ca9-a52f51bd9c2e
using Plots, Turing, PlutoUI, StatsPlots

# ╔═╡ dc1fe99a-c3fc-4c3c-b149-5283dd581f6e
md"""
# Poppy inference 🌺

A farmer has a long field of 1 km. He knows that at some point, the soil type switches from sandy soil to clay-dominated soil. To estimate the point where the soil type switches, the farmer looks at the vegetation. Poppies are known to prefer sand-rich soils over clay. The farmer takes a look over the transect and notes the x-coordinate of where he found poppies sprouting.
"""

# ╔═╡ ed7a6888-f104-4513-91c0-6ab0f7c5ba65
poppy_locations = [
8.39358,
80.4798,
118.559,
177.407,
274.531,
309.858,
384.819,
388.152,
447.347,
544.646,
604.76,
641.367,
]

# ╔═╡ ec80dad4-20a7-4d29-a4d6-3ea4da8c7f43
scatter(poppy_locations, zeros(length(poppy_locations)), color="red", xlab="x [m]", xlims=[0,1000], yaxis=false, label="poppy location")

# ╔═╡ f0c2a3f9-f03c-482d-9df0-8f780e192557
md"""
1. Can you infer the location where the soil switches type using Bayesian reasoning? Logically, this point should be between after the last poppy and before 1000 m. Give the 95% credibility interval!
"""

# ╔═╡ 7c871370-5266-4b50-8971-c7e865990580
@model function poppies(locations)
    # sample boundary
	x_boundary ~ missing

	# number of poppies
	n = missing

	for i in 1:n
		# sample location
		locations[i] ~ missing
	end
end

# ╔═╡ d2843d5a-97d9-4b49-b55d-07a3ee3a432f
md"""

Additional questions (you need to modify the model):
2. The farmer is quite sure the transition spot is between 600 and 1000 meter, with its most likely value at 700 meter. Can you incorportate this information?
3. Can you estimate $\lambda$, the expected number of poppies per meter? Using a Bayesian approach, you can take TriangularDist(0, 0.1) as a prior. Hint: you need to treat the number of poppy locations `n` as a variable to condition on.

"""

# ╔═╡ 33d0742e-5bdb-4fa6-8fd5-01b7eee2b0ee
md"**Solution**"

# ╔═╡ 3a48a515-4040-43ab-a18c-376312061eeb
@model function poppies_sol(locations, n)
    # sample boundary
	
	#x_boundary ~ Uniform(maximum(locations), 1000)
	x_boundary ~ TriangularDist(600, 1000, 700)

	λ ~ TriangularDist(0, 0.1)

	n ~ Poisson(λ * x_boundary)

	for i in 1:n
		locations[i] ~ Uniform(0, x_boundary)
	end
end

# ╔═╡ aae1502f-f253-4891-8209-6583c5abe7ef
chain = sample(poppies_sol(poppy_locations, length(locations)), NUTS(), 10000)

# ╔═╡ d7affc89-0d96-411e-a890-c0c45374053b
quantile(chain)

# ╔═╡ a050822f-e2d1-4966-aec7-d3e95d5b9e37
plot(chain)

# ╔═╡ 621b0c2b-acbb-4a7d-bce7-aa37baf47fde
summarize(chain)

# ╔═╡ eb1a46dd-a68e-476a-a33d-20036d1cbbce
23/674

# ╔═╡ Cell order:
# ╠═6261bb1a-5f73-4c9d-8ca9-a52f51bd9c2e
# ╠═57420478-d744-11ef-2ec5-735d47866238
# ╠═dc1fe99a-c3fc-4c3c-b149-5283dd581f6e
# ╠═ed7a6888-f104-4513-91c0-6ab0f7c5ba65
# ╠═ec80dad4-20a7-4d29-a4d6-3ea4da8c7f43
# ╠═f0c2a3f9-f03c-482d-9df0-8f780e192557
# ╠═7c871370-5266-4b50-8971-c7e865990580
# ╠═d2843d5a-97d9-4b49-b55d-07a3ee3a432f
# ╠═33d0742e-5bdb-4fa6-8fd5-01b7eee2b0ee
# ╠═3a48a515-4040-43ab-a18c-376312061eeb
# ╠═aae1502f-f253-4891-8209-6583c5abe7ef
# ╠═d7affc89-0d96-411e-a890-c0c45374053b
# ╠═a050822f-e2d1-4966-aec7-d3e95d5b9e37
# ╠═621b0c2b-acbb-4a7d-bce7-aa37baf47fde
# ╠═eb1a46dd-a68e-476a-a33d-20036d1cbbce
