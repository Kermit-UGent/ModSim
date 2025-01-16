### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 886c7932-da4b-45cc-ba73-8047389e4895
using Pkg; Pkg.activate("..");

# ╔═╡ 80bc0e86-5ad3-4d61-9600-8dc05b86599d
using Turing, StatsPlots

# ╔═╡ 52a38b60-178b-4a1d-ac32-e73fafd339f9
md"# Sampling notebook #3: Advanced"

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

# ╔═╡ e9ec0d51-784b-4e7a-a145-5aa456ef72cc
md"""
!!! hint
	Assuming the lines are vertical, you only need to consider the **x-coordinates** of both ends of the needle.
"""

# ╔═╡ 870b3d09-de00-4b9c-a67f-a38a4cc289aa
@model function needles()
    x0 ~ Uniform(0, 1)
    theta ~ Uniform(0, 2*pi)
    x_end = x0 + cos(theta)
    crossing = x_end < 0 || x_end > 1
    return crossing
end

# ╔═╡ af516e05-5f26-49f9-a1ca-d1f3d4f8ea99
needle_model = needles()

# ╔═╡ a50f848d-4d27-47e0-bfb6-76ef620660c0
sp_n = [needle_model() for i in 1:2000]

# ╔═╡ 4fb459cc-8449-4a3c-9164-7b89c3348c4e
mean(sp_n)

# ╔═╡ f2061e56-1de3-4be1-8c3b-5bcae649bf4d
pi_est = 2 / mean(sp_n)

# ╔═╡ c8941726-9e81-47fc-9b7e-cb3b5c0c61ca
md"# 2: Attraction"

# ╔═╡ 431023df-3724-4325-b0ac-96dbf5e4fd20
md"""
Following a course on electromagnetism will teach one that computing the net force between 2 arbitrary shapes can be a terrifying task. Tragedy has it then, that this is a very general problem with application from making fusion reactors to space travel. We can ease the pain by turning it into a sampling problem.

We'll start in a humble manner and simulate **the gravitational force between 2 cubes**. Both cubes are size 1. The first cube is in [0, 1] x [0, 1] x [0, 1], and the second cube in [1.1, 2.1] x [0, 1] x [0, 1], as shown in the figure below.
"""

# ╔═╡ 599ac984-ef1d-4c7a-8e87-9d4ddb1aa710
begin
	xe = [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0]
	ye = [0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1]
	ze = [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1]

	xe2 = [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0] .+ 1.1

	plot(xlims = (-0.5, 2.5), ylims = (-1, 2), zlims = (0, 3))
	plot!(xe, ye, ze; color = :blue, linewidth = 0.5, label = "cube 1")
	plot!(xe2, ye, ze; color = :orange, lw = 0.5, label = "cube 2")
end

# ╔═╡ bb8e4cd6-6106-4d3f-9e35-0dfd8b2c45f5
md"""
The gravitational force can be estimated by **randomly sampling a point from both cubes** and using a simplification of the formula for gravitational force between those points:

```math
F = \frac{1}{r^2}
```
"""

# ╔═╡ 058d2a08-7d40-46c5-8c5b-2b5ce9d2eb18
md"""
!!! questions
	- What is the estimated net force between the two cubes? Is this the same as if you had treated the cubes as point masses?
	- How many samples do you need to estimate this force reliably? Define a reliable estimator as one having a standard deviation of 10% of the mean value. Visualise the distribution of the estimator.
"""

# ╔═╡ 06d0e92f-1f07-41fd-b6ee-e94eb539627d
@model function cubeforce()
    x1 ~ Uniform(0, 1)
    y1 ~ Uniform(0, 1)
    z1 ~ Uniform(0, 1)

	x2 ~ Uniform(1.1, 2.1)
    y2 ~ Uniform(0, 1)
    z2 ~ Uniform(0, 1)

	r_squared = (x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2
	G = 1 # for simplicity
	F = G / r_squared
	
    return F
end

# ╔═╡ 43119060-bd35-4c4c-8831-d5c5df1d8dd5
cubemodel = cubeforce();

# ╔═╡ d5229ed1-8473-42e5-b1f9-98fb547d3e85
force_sp = [cubemodel() for _ in 1:2000];

# ╔═╡ 420120c3-1c8c-46e4-a7c0-0f2405dd159a
pointmass_force = 1 / 1.1^2

# ╔═╡ 1ce8f808-b917-4831-89b4-4c318f05724d
mean(force_sp)

# ╔═╡ a324a615-7066-4dbb-9410-95e1e1c46ac4
histogram(force_sp)

# ╔═╡ e38e8f51-90db-4136-84e2-f06cd03d502a
required_samples = 150 # manually change until standard deviation is about 0.1

# ╔═╡ a7e975a9-beb5-4169-9d22-3e970be2838b
averaged_force_sp = [mean([cubemodel() for _ in 1:required_samples]) for _ in 1:2000];

# ╔═╡ 7ff5cb25-47bd-432f-a404-359abaf44e3a
std(averaged_force_sp)

# ╔═╡ 99c7f6d6-eada-491b-b966-fdb1195fc111
histogram(averaged_force_sp)

# ╔═╡ Cell order:
# ╟─52a38b60-178b-4a1d-ac32-e73fafd339f9
# ╠═886c7932-da4b-45cc-ba73-8047389e4895
# ╠═80bc0e86-5ad3-4d61-9600-8dc05b86599d
# ╟─1ff166c0-b639-11ef-2143-ad11b6834d4c
# ╟─5ec0a0ee-0c5f-437e-be85-162805af4731
# ╟─1d74c5dc-8cae-4c2d-9a17-d9a6c52a9169
# ╟─4a67aba8-6015-4c95-9b5b-ca4fd2bed1bd
# ╟─2382f49a-2e76-4e30-b1d0-eca298779d73
# ╟─2f277afb-eb92-44bf-b4fb-b994fefb7266
# ╟─e9ec0d51-784b-4e7a-a145-5aa456ef72cc
# ╠═870b3d09-de00-4b9c-a67f-a38a4cc289aa
# ╠═af516e05-5f26-49f9-a1ca-d1f3d4f8ea99
# ╠═a50f848d-4d27-47e0-bfb6-76ef620660c0
# ╠═4fb459cc-8449-4a3c-9164-7b89c3348c4e
# ╠═f2061e56-1de3-4be1-8c3b-5bcae649bf4d
# ╟─c8941726-9e81-47fc-9b7e-cb3b5c0c61ca
# ╟─431023df-3724-4325-b0ac-96dbf5e4fd20
# ╟─599ac984-ef1d-4c7a-8e87-9d4ddb1aa710
# ╟─bb8e4cd6-6106-4d3f-9e35-0dfd8b2c45f5
# ╟─058d2a08-7d40-46c5-8c5b-2b5ce9d2eb18
# ╠═06d0e92f-1f07-41fd-b6ee-e94eb539627d
# ╠═43119060-bd35-4c4c-8831-d5c5df1d8dd5
# ╠═d5229ed1-8473-42e5-b1f9-98fb547d3e85
# ╠═420120c3-1c8c-46e4-a7c0-0f2405dd159a
# ╠═1ce8f808-b917-4831-89b4-4c318f05724d
# ╠═a324a615-7066-4dbb-9410-95e1e1c46ac4
# ╠═e38e8f51-90db-4136-84e2-f06cd03d502a
# ╠═a7e975a9-beb5-4169-9d22-3e970be2838b
# ╠═7ff5cb25-47bd-432f-a404-359abaf44e3a
# ╠═99c7f6d6-eada-491b-b966-fdb1195fc111
