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
dropletdist = missing;

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
    r ~ missing
	K ~ missing

	logfun = missing
    return logfun
end

# ╔═╡ 4345b3dd-0731-4dd5-a319-627b4f91306e
petri_model = missing

# ╔═╡ 9d747eef-a883-49b3-acb7-f0d077a2b902
chain_petri = missing

# ╔═╡ ce8b53c3-0af1-4e9a-ad80-6fd7a2bf020b
logfuns = missing

# ╔═╡ e1d0c5d5-2b7a-4af4-a7ac-e1ca46e665c8
sp_petri = missing

# ╔═╡ e3092915-a083-425e-8fa1-b7bf370abc8a
prob_splittable = missing

# ╔═╡ e4f42f16-5cce-4fc2-aa01-8971f37c710e
md"### 3"

# ╔═╡ d6c193d7-9fe2-413b-801f-ebc33c772ee9
missing # plots

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

# ╔═╡ ceeeee76-e31c-4429-8ed0-e1c503433dbf
md"""
!!! questions
	1. What is the estimated net force between the two cubes? Is this the same as if you had treated the cubes as point masses?
	1. How many samples do you need to estimate this force reliably? Define a reliable estimator as one having a standard deviation of 10% of the mean value. Visualise the distribution of the estimator.
"""

# ╔═╡ 9be28327-d87f-4d50-bbcc-91e799f14dbf
md"### 1"

# ╔═╡ 06d0e92f-1f07-41fd-b6ee-e94eb539627d
@model function cubeforce()
    x1 ~ missing
    y1 ~ missing
    z1 ~ missing

	x2 ~ missing
    y2 ~ missing
    z2 ~ missing

	F = missing
	
    return F
end

# ╔═╡ 43119060-bd35-4c4c-8831-d5c5df1d8dd5
cubemodel = cubeforce();

# ╔═╡ d26ad946-2bb4-4383-860f-d601903ce1be
force_sp = missing

# ╔═╡ 7b117453-0875-4b41-a228-866c6c0a8208
force_average = missing

# ╔═╡ b28cfbae-2fac-4b38-b234-53f71e381bcd
pointmass_force = missing # (doesn't require Turing, only maths)

# ╔═╡ 304d6052-6d3d-487c-8c01-dab259276d6f
md"### 2"

# ╔═╡ e38e8f51-90db-4136-84e2-f06cd03d502a
@bind required_samples Slider(10:10:200, show_value = true)

# ╔═╡ b753fa82-28ab-46d1-9085-4c65301db046
estimator = mean([cubemodel() for _ in 1:required_samples]) 
	# a single estimation of the force given `required_samples` samples

# ╔═╡ 55cf127f-5534-4926-955a-487ea9553b70
estimator_sp = missing 
	# a sample of estimations given `required_samples` samples

# ╔═╡ f6a67c1e-3677-4e51-b6dd-5d11a5146ea7
estimator_σ = missing
	# standard deviation of the force estimator

# ╔═╡ 45f0813b-c4c9-4f13-8d66-1e58293c4422
missing # histogram

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
# ╠═d6c193d7-9fe2-413b-801f-ebc33c772ee9
# ╟─c8941726-9e81-47fc-9b7e-cb3b5c0c61ca
# ╟─431023df-3724-4325-b0ac-96dbf5e4fd20
# ╟─599ac984-ef1d-4c7a-8e87-9d4ddb1aa710
# ╟─bb8e4cd6-6106-4d3f-9e35-0dfd8b2c45f5
# ╟─ceeeee76-e31c-4429-8ed0-e1c503433dbf
# ╟─9be28327-d87f-4d50-bbcc-91e799f14dbf
# ╠═06d0e92f-1f07-41fd-b6ee-e94eb539627d
# ╠═43119060-bd35-4c4c-8831-d5c5df1d8dd5
# ╠═d26ad946-2bb4-4383-860f-d601903ce1be
# ╠═7b117453-0875-4b41-a228-866c6c0a8208
# ╠═b28cfbae-2fac-4b38-b234-53f71e381bcd
# ╟─304d6052-6d3d-487c-8c01-dab259276d6f
# ╠═e38e8f51-90db-4136-84e2-f06cd03d502a
# ╠═b753fa82-28ab-46d1-9085-4c65301db046
# ╠═55cf127f-5534-4926-955a-487ea9553b70
# ╠═f6a67c1e-3677-4e51-b6dd-5d11a5146ea7
# ╠═45f0813b-c4c9-4f13-8d66-1e58293c4422
