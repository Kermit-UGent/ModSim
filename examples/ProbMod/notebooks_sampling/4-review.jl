### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 882b942c-a5ae-445d-a946-6eb8cddc423e
using Pkg; Pkg.activate("..")

# ╔═╡ aa2737db-f89f-4e08-9abe-0a86e8880c19
using Turing, StatsPlots

# ╔═╡ 38e71480-b874-11ef-0388-155c2fbaeb60
md"# Review exercise: The Delta Works"

# ╔═╡ 20d494d3-e296-4ba5-b496-1cfd4add17ec
md"""
The Netherlands relies on a large number of constructions such as dams and dykes to prevent the country from being flooded. One important question when planning these constructions is: how high a wave should they be able to deal with?

The simplest answer to this question would be "as high as the highest wave that will hit these constructions". Starting from the probability distributions of the number and height of waves, this comes down to getting the probability distribution of the maximum height. Answering this question analytically is difficult enough that it warrants [its own field of mathematics](https://en.wikipedia.org/wiki/Extreme_value_theory), so we'll make use of sampling methods instead. 
"""

# ╔═╡ 8ada0eec-a9cf-479d-82ef-32ec2ed36ae6
md"""
Make the following assumptions:
- Every year has about 12 storm events (you can ignore waves that are not from a storm event).
- On a day with a storm event, 100 or so waves reach the shore.
- The height of a storm wave becomes less probable the larger it is.
- Storm waves have a 10/12 chance to be from a normal storm event with a mean wave height of 100 cm, or a 2/12 chance to be from a strong storm event with a mean wave height of 200 cm.
"""

# ╔═╡ 08c02317-3ea4-4dd4-a269-580281c1b6b1
md"""
!!! question
	What wave height should the Delta Works be able to deal with so that in the following 500 years, no wave higher than this value will occur with 99% probability?
"""

# ╔═╡ 726542f7-c7a7-416d-b47c-940338e295a2
md"""
!!! note
	This question can become computationally quite heavy. You can remedy this by replacing for-loops with the more efficient `filldist` function (see the docs!).
"""

# ╔═╡ 98deb640-4610-4d61-adb0-8d07b5a6bb1d
heightdist = MixtureModel([Exponential(100), Exponential(200)], [10/12, 2/12]);

# ╔═╡ 8e58dfaf-61c4-49f0-94e3-ecb55e349c80
rand(heightdist, 10000) |> histogram

# ╔═╡ 7a11c495-06ca-4a36-b594-0a720c25c42b
@model function waveheight(years)
	numstorms ~ Poisson(12 * years)
	numwaves ~ Poisson(100 * numstorms)
	
	wave_heights ~ filldist(heightdist, numwaves)

	return maximum(wave_heights)
end

# ╔═╡ 5ed65d57-3a1f-4c1e-8985-e03335b64427
wavemodel = waveheight(500)

# ╔═╡ b9352eaf-3af8-4f3b-91fa-088f7d866494
sp_wh = [wavemodel() for i in 1:200]

# ╔═╡ fb5b771e-7795-4d43-a570-0282fffb53f4
histogram(sp_wh, bins = 20)

# ╔═╡ a572a1d3-cf55-4b2b-9343-de690dcd69aa
quantile(sp_wh, 0.99)

# ╔═╡ 49266b7a-8a43-4f53-969b-e406e4adc41c
md"They should be able to deal with waves of ~30m."

# ╔═╡ Cell order:
# ╠═882b942c-a5ae-445d-a946-6eb8cddc423e
# ╠═aa2737db-f89f-4e08-9abe-0a86e8880c19
# ╟─38e71480-b874-11ef-0388-155c2fbaeb60
# ╟─20d494d3-e296-4ba5-b496-1cfd4add17ec
# ╟─8ada0eec-a9cf-479d-82ef-32ec2ed36ae6
# ╟─08c02317-3ea4-4dd4-a269-580281c1b6b1
# ╟─726542f7-c7a7-416d-b47c-940338e295a2
# ╠═98deb640-4610-4d61-adb0-8d07b5a6bb1d
# ╠═8e58dfaf-61c4-49f0-94e3-ecb55e349c80
# ╠═7a11c495-06ca-4a36-b594-0a720c25c42b
# ╠═5ed65d57-3a1f-4c1e-8985-e03335b64427
# ╠═b9352eaf-3af8-4f3b-91fa-088f7d866494
# ╠═fb5b771e-7795-4d43-a570-0282fffb53f4
# ╠═a572a1d3-cf55-4b2b-9343-de690dcd69aa
# ╟─49266b7a-8a43-4f53-969b-e406e4adc41c
