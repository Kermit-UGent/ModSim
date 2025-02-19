### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 886c7932-da4b-45cc-ba73-8047389e4895
using Pkg; Pkg.activate("..");

# ╔═╡ 80bc0e86-5ad3-4d61-9600-8dc05b86599d
using Turing, StatsPlots, PlutoUI

# ╔═╡ 52a38b60-178b-4a1d-ac32-e73fafd339f9
md"# Sampling notebook #3: Advanced"

# ╔═╡ 116b840c-e766-4ff6-aafa-0977fb122992
md"## 1: Petridish peril"

# ╔═╡ 415b5ba8-3f6d-46ea-8f89-19fa7c0e74f9
md"""
Living the microbiology master thesis life, your mornings consist of inoculating petridishes with bacteria. Somewhere along the day, you need to split them. You want to do this **after** there's a decent amount of bacteria in the dish (>10\_000) but **before** they have overgrown the entire dish and start dying (<100\_000). This condition we call **splittable**.

You'd like to estimate how long after inoculation you should return to your bacteria so that they're most likely to be in a splittable state.
"""

# ╔═╡ 3691d6aa-c717-46e8-b8b3-f4aa56c9f761
md"""
Bacteria follow **logistic growth**, and you can use the following assumptions:
- The initial population size $P_0$ has a 75% chance of originating from a small droplet and a 25% chance for a big droplet
  - For small droplets, `P0` follows a `Poisson(10)`
  - For big droplets, `P0` follows a `Poisson(30)`
- The growth rate $r$ follows a `LogNormal(0.0, 0.3)`
- The growth capacity $K$ of the inoculated medium follows a `Normal(1e5, 1e4)`
"""

# ╔═╡ 749cacb9-b73d-470e-bf58-5550db5de7e0
md"""
!!! questions
	1. Plot the prior distribution of P0.
	2. What is the probability your bacteria are in a splittable state 8 hours after inoculation?
	3. Plot 100 of the sampled logistic growth curves from 0 to 12 hours.
"""

# ╔═╡ cec8b8e8-850e-4549-97fc-71eb35b8334b
md"### 1"

# ╔═╡ dbccef02-51e7-48dd-9dab-be44db064100
md"""
!!! tip
	A simple way of representing the distribution of P0 is through a mixture model. Look up the documentation of `MixtureModel` for how to make one in Turing.
"""

# ╔═╡ eb6dc4e7-e779-4bfc-b865-3defa3894181
dropletdist = MixtureModel([Poisson(10), Poisson(30)], [0.75, 0.25]);

# ╔═╡ 6517d682-b524-42e3-8a13-78ed5c1ce0dc
plot(dropletdist) # plot gives wrong result, the lines should stack

# ╔═╡ 60eabb68-7d65-4f78-97d9-b1e20c2dbd0a
rand(dropletdist, 10000) |> histogram 
	# plotting a histogram of random samples gives a better result

# ╔═╡ 107533fc-b300-4b8d-bea2-a3aa6a37938d
md"### 2"

# ╔═╡ 38c5b9cd-0baf-4a62-912e-4993614ddbf3
md"""
!!! tip
	You can `return` the logistic function estimated within the model and retrieve it using `generated_quantities` to make plotting easier later on.

	Remember: anonymous functions can be defined using `myfun = x -> ...`
"""

# ╔═╡ bde57599-1dca-41a4-94aa-498da72c2012
logistic(t, P0, r, K) =  K / (1 + (K - P0)/P0 * exp(-r*t))

# ╔═╡ 976cfb96-3196-4a61-bd5f-e4f5d24ba1e9
@model function petrigrowth()
	P0 ~ dropletdist
    r ~ LogNormal(0.0, 0.3)
	K ~ Normal(1e5, 1e4)

	logfun = t -> logistic(t, P0, r, K)
    return logfun
end

# ╔═╡ 4345b3dd-0731-4dd5-a319-627b4f91306e
petri_model = petrigrowth();

# ╔═╡ 9d747eef-a883-49b3-acb7-f0d077a2b902
chain_petri = sample(petri_model, Prior(), 2000);

# ╔═╡ ce8b53c3-0af1-4e9a-ad80-6fd7a2bf020b
logfuns = generated_quantities(petri_model, chain_petri);

# ╔═╡ e1d0c5d5-2b7a-4af4-a7ac-e1ca46e665c8
sp_petri = [logfun(8.0) for logfun in logfuns];

# ╔═╡ e3092915-a083-425e-8fa1-b7bf370abc8a
prob_splittable = mean((sp_petri .>= 1e4) .&& (sp_petri .<= 1e5))

# ╔═╡ e4f42f16-5cce-4fc2-aa01-8971f37c710e
md"### 3"

# ╔═╡ 14a8bc88-e10e-42d1-a4c0-7f3ebc5512a9
plot(logfuns[1:100], color = :violet, alpha = 0.3, label = false, xlims = (0, 12))

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
The gravitational force can be estimated by **randomly sampling a point from both cubes** and using the formula for gravitational force between those points, ignoring all constants:

```math
F = \frac{1}{r^2}
```
"""

# ╔═╡ 058d2a08-7d40-46c5-8c5b-2b5ce9d2eb18
md"""
!!! questions
	1. What is the estimated net force between the two cubes? Is this the same as if you had treated the cubes as point masses?
	1. How many samples do you need to estimate this force reliably? Define a reliable estimator as one having a standard deviation of 10% of the mean value. Visualise the distribution of the estimator.
"""

# ╔═╡ 3cd12379-5e7f-4f5d-a7bb-8b3a9d2e8497
md"### 1"

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

# ╔═╡ 1ce8f808-b917-4831-89b4-4c318f05724d
force_average = mean(force_sp)

# ╔═╡ 420120c3-1c8c-46e4-a7c0-0f2405dd159a
pointmass_force = 1 / 1.1^2

# ╔═╡ a324a615-7066-4dbb-9410-95e1e1c46ac4
histogram(force_sp)

# ╔═╡ 6d9209da-2c98-454e-88e6-32fe667efff7
md"### 2"

# ╔═╡ e38e8f51-90db-4136-84e2-f06cd03d502a
@bind required_samples Slider(10:10:200, default = 140, show_value = true)
	# manually change until standard deviation is about 0.1 (~10% of average force)

# ╔═╡ 788f60fc-b1a1-4635-a36a-f954738bcffc
estimator = mean([cubemodel() for _ in 1:required_samples])

# ╔═╡ a7e975a9-beb5-4169-9d22-3e970be2838b
estimator_sp = [mean([cubemodel() for _ in 1:required_samples]) for _ in 1:2000];

# ╔═╡ 7ff5cb25-47bd-432f-a404-359abaf44e3a
estimator_σ = std(estimator_sp)

# ╔═╡ 99c7f6d6-eada-491b-b966-fdb1195fc111
histogram(estimator_sp)

# ╔═╡ Cell order:
# ╟─52a38b60-178b-4a1d-ac32-e73fafd339f9
# ╠═886c7932-da4b-45cc-ba73-8047389e4895
# ╠═80bc0e86-5ad3-4d61-9600-8dc05b86599d
# ╟─116b840c-e766-4ff6-aafa-0977fb122992
# ╟─415b5ba8-3f6d-46ea-8f89-19fa7c0e74f9
# ╟─3691d6aa-c717-46e8-b8b3-f4aa56c9f761
# ╟─749cacb9-b73d-470e-bf58-5550db5de7e0
# ╟─cec8b8e8-850e-4549-97fc-71eb35b8334b
# ╟─dbccef02-51e7-48dd-9dab-be44db064100
# ╠═eb6dc4e7-e779-4bfc-b865-3defa3894181
# ╠═6517d682-b524-42e3-8a13-78ed5c1ce0dc
# ╠═60eabb68-7d65-4f78-97d9-b1e20c2dbd0a
# ╟─107533fc-b300-4b8d-bea2-a3aa6a37938d
# ╟─38c5b9cd-0baf-4a62-912e-4993614ddbf3
# ╠═bde57599-1dca-41a4-94aa-498da72c2012
# ╠═976cfb96-3196-4a61-bd5f-e4f5d24ba1e9
# ╠═4345b3dd-0731-4dd5-a319-627b4f91306e
# ╠═9d747eef-a883-49b3-acb7-f0d077a2b902
# ╠═ce8b53c3-0af1-4e9a-ad80-6fd7a2bf020b
# ╠═e1d0c5d5-2b7a-4af4-a7ac-e1ca46e665c8
# ╠═e3092915-a083-425e-8fa1-b7bf370abc8a
# ╟─e4f42f16-5cce-4fc2-aa01-8971f37c710e
# ╠═14a8bc88-e10e-42d1-a4c0-7f3ebc5512a9
# ╟─c8941726-9e81-47fc-9b7e-cb3b5c0c61ca
# ╟─431023df-3724-4325-b0ac-96dbf5e4fd20
# ╟─599ac984-ef1d-4c7a-8e87-9d4ddb1aa710
# ╟─bb8e4cd6-6106-4d3f-9e35-0dfd8b2c45f5
# ╟─058d2a08-7d40-46c5-8c5b-2b5ce9d2eb18
# ╟─3cd12379-5e7f-4f5d-a7bb-8b3a9d2e8497
# ╠═06d0e92f-1f07-41fd-b6ee-e94eb539627d
# ╠═43119060-bd35-4c4c-8831-d5c5df1d8dd5
# ╠═d5229ed1-8473-42e5-b1f9-98fb547d3e85
# ╠═1ce8f808-b917-4831-89b4-4c318f05724d
# ╠═420120c3-1c8c-46e4-a7c0-0f2405dd159a
# ╠═a324a615-7066-4dbb-9410-95e1e1c46ac4
# ╟─6d9209da-2c98-454e-88e6-32fe667efff7
# ╠═e38e8f51-90db-4136-84e2-f06cd03d502a
# ╠═788f60fc-b1a1-4635-a36a-f954738bcffc
# ╠═a7e975a9-beb5-4169-9d22-3e970be2838b
# ╠═7ff5cb25-47bd-432f-a404-359abaf44e3a
# ╠═99c7f6d6-eada-491b-b966-fdb1195fc111
