### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 886c7932-da4b-45cc-ba73-8047389e4895
using Pkg; Pkg.activate("../..");

# ╔═╡ 80bc0e86-5ad3-4d61-9600-8dc05b86599d
using Turing, StatsPlots

# ╔═╡ f2327bc5-35ff-4d2d-b2dd-5a3964a417be
n_samples = 1000

# ╔═╡ 1ff166c0-b639-11ef-2143-ad11b6834d4c
md"# 1: Buffon's needles"

# ╔═╡ 5ec0a0ee-0c5f-437e-be85-162805af4731
md"""
A wise man once said: ["there is no greater joy than estimating π"](https://en.wikipedia.org/wiki/Approximations_of_%CF%80). Next to throwing darts at the unit square, another method to accomplish this is using [Buffon's needle problem](https://en.wikipedia.org/wiki/Buffon%27s_needle_problem).

The experiment is as follows: consider a floor with parallel lines all a distance of 1 away from eachother. Now drop a needle of length 1 (and width ~0) on the floor with a **random position and angle**. What is the probability $P_{cross}$ that the needle will cross one of the lines?
"""

# ╔═╡ 1d74c5dc-8cae-4c2d-9a17-d9a6c52a9169
md"The following image illustrates the problem (imagine $l$ = $t$ = 1) for two needles, where `a` crosses a line and `b` does not."

# ╔═╡ 4a67aba8-6015-4c95-9b5b-ca4fd2bed1bd
md"""
![Buffon's needles](https://upload.wikimedia.org/wikipedia/commons/thumb/5/58/Buffon_needle.svg/1920px-Buffon_needle.svg.png)
"""

# ╔═╡ 2382f49a-2e76-4e30-b1d0-eca298779d73
md"""
Using sampling magic, it's not difficult to make an estimate of this probability, $\hat{P}_{cross}$. Solving the problem analytically shows that the exact value is:
```math
P_{cross} = \frac{2}{\pi}
```

Therefore, our estimator for π is:
```math
π = \frac{2}{\hat{P}_{cross}}
```
"""

# ╔═╡ 2f277afb-eb92-44bf-b4fb-b994fefb7266
md"""
!!! question
	Estimate π using the Buffon's needle approximation.
"""

# ╔═╡ 870b3d09-de00-4b9c-a67f-a38a4cc289aa
@model function needles()
    x0 ~ Uniform(0, 1)
    theta ~ Uniform(0, pi)
    x_end = x0 + cos(theta)
    crossing = x_end < 0 || x_end > 1
    return crossing
end

# ╔═╡ af516e05-5f26-49f9-a1ca-d1f3d4f8ea99
needle_model = needles()

# ╔═╡ a50f848d-4d27-47e0-bfb6-76ef620660c0
sp_n = [needle_model() for i in 1:n_samples]

# ╔═╡ 4fb459cc-8449-4a3c-9164-7b89c3348c4e
mean(sp_n)

# ╔═╡ f2061e56-1de3-4be1-8c3b-5bcae649bf4d
pi_est = 2 / mean(sp_n)

# ╔═╡ 239deeda-98e2-4035-92a4-b489b272101f
md"# 2: The Delta Works"

# ╔═╡ 53a1cffb-1ae9-41f4-a7b7-2649e03a1c43
md"""
The Netherlands relies on a large number of constructions such as dams and dykes to prevent the country from being flooded. One important question when planning these constructions is: how high a wave should they be able to deal with?

The simplest answer to this question would be "as high as the highest wave that will hit these constructions". Starting from the probability distributions of the number and height of waves, this comes down to getting the probability distribution of the maximum height. Answering this question analytically is difficult enough that it warrants [its own field of mathematics](https://en.wikipedia.org/wiki/Extreme_value_theory), so we'll make use of sampling methods instead. 
"""

# ╔═╡ 0dd1245d-e445-4c96-8b83-634d71a11234
md"""
Make the following assumptions:
- The yearly amount of waves that reach the shore is about 100, with a similar variance
- The height of a wave becomes less probable the larger it is and has a mean of 100 cm
- Waves have a chance to be from a **storm** event, in which case the height distribution changes to have a mean of 200 cm.
- The chance of a storm event is 1% ± 0.1% .
"""

# ╔═╡ 3def2f49-aa48-4427-98b1-984a90b77f82
md"""
!!! question
	What wave height should the Delta Works be able to deal with so that in the following 500 years, no wave larger than this value will occur with 99% probability?
"""

# ╔═╡ fbc597ea-4e20-4c24-b4b1-6d0278493eae
md"""
!!! note
	Despite grossly underestimating the amount of yearly waves, this question can become computationally quite heavy. You can remedy this by replacing for-loops with the more efficient `filldist` function (see the docs!).
"""

# ╔═╡ cc93fa8a-06a3-4b9b-b80e-2d84714a58b3
plot(Poisson(100), title = "Distribution of yearly number of waves", legend = false)

# ╔═╡ 2ff42a70-6db7-4755-9392-14b6d2350c72
plot(Exponential(100), title = "Distribution of wave height", legend = false)

# ╔═╡ 266d32a3-ae67-4edc-88b9-70bc04543f0a
plot(Exponential(200), title = "Distribution of stormwave height", legend = false)

# ╔═╡ ed9eacb3-180e-4f35-b19d-4180c16e422f
@model function waveheight(years)
	yearly_nums ~ filldist(Poisson(100), years)
    total_wavenum = sum(yearly_nums)

	stormchance ~ Normal(0.01, 0.001)

	wavedist = MixtureModel([Exponential(100), Exponential(200)], [1-stormchance, stormchance])
	wave_heights ~ filldist(wavedist, total_wavenum)

	return maximum(wave_heights)
end

# ╔═╡ b8ab3318-460d-4fba-90e0-7da0df0c2a79
wavemodel = waveheight(500)

# ╔═╡ cc120f28-feca-4f20-bf7c-9da80ea02e63
sp_wh = [wavemodel() for i in 1:n_samples]

# ╔═╡ 6f7c3e7b-478d-4182-b206-9af830f2d9b7
histogram(sp_wh)

# ╔═╡ 1f7fd192-2fd6-49a2-9f35-ba67d6516c61
quantile(sp_wh, 0.99)

# ╔═╡ 72304beb-129a-4787-ac39-f85494348861
md"They should be able to deal with waves of ~22m."

# ╔═╡ Cell order:
# ╠═886c7932-da4b-45cc-ba73-8047389e4895
# ╠═80bc0e86-5ad3-4d61-9600-8dc05b86599d
# ╠═f2327bc5-35ff-4d2d-b2dd-5a3964a417be
# ╟─1ff166c0-b639-11ef-2143-ad11b6834d4c
# ╟─5ec0a0ee-0c5f-437e-be85-162805af4731
# ╟─1d74c5dc-8cae-4c2d-9a17-d9a6c52a9169
# ╟─4a67aba8-6015-4c95-9b5b-ca4fd2bed1bd
# ╟─2382f49a-2e76-4e30-b1d0-eca298779d73
# ╟─2f277afb-eb92-44bf-b4fb-b994fefb7266
# ╠═870b3d09-de00-4b9c-a67f-a38a4cc289aa
# ╠═af516e05-5f26-49f9-a1ca-d1f3d4f8ea99
# ╠═a50f848d-4d27-47e0-bfb6-76ef620660c0
# ╠═4fb459cc-8449-4a3c-9164-7b89c3348c4e
# ╠═f2061e56-1de3-4be1-8c3b-5bcae649bf4d
# ╟─239deeda-98e2-4035-92a4-b489b272101f
# ╟─53a1cffb-1ae9-41f4-a7b7-2649e03a1c43
# ╟─0dd1245d-e445-4c96-8b83-634d71a11234
# ╟─3def2f49-aa48-4427-98b1-984a90b77f82
# ╟─fbc597ea-4e20-4c24-b4b1-6d0278493eae
# ╟─cc93fa8a-06a3-4b9b-b80e-2d84714a58b3
# ╟─2ff42a70-6db7-4755-9392-14b6d2350c72
# ╟─266d32a3-ae67-4edc-88b9-70bc04543f0a
# ╠═ed9eacb3-180e-4f35-b19d-4180c16e422f
# ╠═b8ab3318-460d-4fba-90e0-7da0df0c2a79
# ╠═cc120f28-feca-4f20-bf7c-9da80ea02e63
# ╠═6f7c3e7b-478d-4182-b206-9af830f2d9b7
# ╠═1f7fd192-2fd6-49a2-9f35-ba67d6516c61
# ╟─72304beb-129a-4787-ac39-f85494348861
