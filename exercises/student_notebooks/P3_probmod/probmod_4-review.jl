### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 882b942c-a5ae-445d-a946-6eb8cddc423e
using Pkg; Pkg.activate("..")

# ╔═╡ aa2737db-f89f-4e08-9abe-0a86e8880c19
using Turing, StatsPlots

# ╔═╡ 5be87389-1a97-481e-bec6-1d781f016077
md"# Review exercise: Buffon's needles"

# ╔═╡ a6435b94-f1af-4609-abfc-93d88730d023
md"""
A wise man once said: ["there is no greater joy than estimating π"](https://en.wikipedia.org/wiki/Approximations_of_%CF%80). Next to throwing darts at the unit square, another method to accomplish this is using [Buffon's needle problem](https://en.wikipedia.org/wiki/Buffon%27s_needle_problem).

The experiment is as follows: consider a floor with parallel lines all a distance of 1 away from eachother. Now drop a needle of length 1 (and width ~0) on the floor with a **random position and angle**. What is the probability $P_{cross}$ that the needle will cross one of the lines?
"""

# ╔═╡ df1758da-50fa-4476-a132-b952620110c6
md"The following image illustrates the problem (imagine $l$ = $t$ = 1) for two needles, where `a` crosses a line and `b` does not."

# ╔═╡ c2ba2183-957f-4b94-a22e-0c2f3dd957ad
md"""
![Buffon's needles](https://upload.wikimedia.org/wikipedia/commons/thumb/5/58/Buffon_needle.svg/1920px-Buffon_needle.svg.png)
"""

# ╔═╡ d71cf8dc-99e5-48e0-9abe-2242a6ccc30b
md"""
Using sampling magic, it's not difficult to make an estimate of this probability, $\hat{P}_{cross}$. Solving the problem analytically shows that the exact value is:
```math
P_{cross} = \frac{2}{\pi}
```

Therefore, our estimator for π is:
```math
\hat{π} = \frac{2}{\hat{P}_{cross}}
```
"""

# ╔═╡ 99747189-d201-4583-ac9f-6875b0b606f2
md"""
!!! question
	Estimate π using the Buffon's needle approximation.
"""

# ╔═╡ 6656cc86-db29-42ef-b612-c09252edfd49
md"""
!!! hint
	Assuming the lines are vertical, you only need to consider the **x-coordinates** of both ends of the needle.
"""

# ╔═╡ Cell order:
# ╠═882b942c-a5ae-445d-a946-6eb8cddc423e
# ╠═aa2737db-f89f-4e08-9abe-0a86e8880c19
# ╟─5be87389-1a97-481e-bec6-1d781f016077
# ╟─a6435b94-f1af-4609-abfc-93d88730d023
# ╟─df1758da-50fa-4476-a132-b952620110c6
# ╟─c2ba2183-957f-4b94-a22e-0c2f3dd957ad
# ╟─d71cf8dc-99e5-48e0-9abe-2242a6ccc30b
# ╟─99747189-d201-4583-ac9f-6875b0b606f2
# ╟─6656cc86-db29-42ef-b612-c09252edfd49
