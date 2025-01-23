### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 94c6f31d-1a43-4221-b60c-1fa0ef8738b8
using Pkg; Pkg.activate("..")

# ╔═╡ 45bc5b66-c81b-4afb-8a7e-51aff9609c62
using Turing, StatsPlots

# ╔═╡ 41dc8060-cf5e-11ef-26f9-892577e77af0
md"# Inference notebook #2: Basics"

# ╔═╡ ed7c5547-bd82-4ca7-bf51-dd1c201d88af
md"## 1: Lights out"

# ╔═╡ a3c26839-2acc-49d8-98c4-3a7142dd6512
md"""
You use 4 of the same LED light in your room. Recently, you've had to replace one for the fourth time. In total, they lasted for 16000, 20000, 23000 and 41000 hours. They _should_ all have lasted 30000 - 50000 hours according to the manufacturer. You want to estimate the chance that the mean lifetime `μ` is smaller than advertised and the manufacturers are liars.

Make the following assumptions:
- The lifetime of a light follows an exponential distribution.
- The prior knowledge on the lifetime of the lights can be encoded as a `LogNormal(log(4e4), 0.3)` distribution.
"""

# ╔═╡ dac55632-6e17-479e-8090-5cf8eaa67dad
md"""
!!! questions
	- Plot the suggested prior distribution for the mean lifetime `μ`. What is the prior probability that `μ` < 30000?
	- What is the posterior probability that `μ` < 30000?
"""

# ╔═╡ 122b1b8f-8ba4-4d2d-9ffb-f640f04889bb
lights_prior = LogNormal(log(4e4), 0.3)

# ╔═╡ 93411707-e019-4348-9728-31cc0cae6040
plot(lights_prior)

# ╔═╡ a5bbb8ce-d134-4d7e-8ad7-54c2b06210ca
cdf(lights_prior, 30_000)

# ╔═╡ 126d3954-c77d-4c98-abe0-fd87d14e6265
@model function lights()
	μ ~ LogNormal(log(4e4), 0.3)
	lightlifes = zeros(4)
	for i in 1:4
		lightlifes[i] ~ Exponential(μ)
	end
end

# ╔═╡ 15fe25d9-5314-4056-8021-cf259ba27c94
lightmodel = lights() | (lightlifes = [16000, 20000, 23000, 41000],)

# ╔═╡ 2a142576-aec2-4311-9c0e-cb68665d59f6
lightschain = sample(lightmodel, NUTS(), 2000)

# ╔═╡ 0130c903-e3dd-415e-aa17-1b705f9e4ccc
plot(lightschain)

# ╔═╡ 11676f8a-d4fe-4999-9ea2-486fca461f65
mean(lightschain[:μ] .< 30_000)

# ╔═╡ ca8ec759-4a4f-4401-a747-aa904d3ee7b3
md"From the available data, there's about a 25% chance the manufacturers are liars."

# ╔═╡ dedbd315-dce1-47f0-bc49-47325cd5b170
md"## 2: Kitties"

# ╔═╡ db4af5b9-d424-458d-aeaa-f803dd06a222
md"""
Staying over at a friend's place, you wonder how much their two cats, Mochi and Momo, like you. In an attempt to find out, you're going to try to pet them a few times. You dub the success rate of them allowing it as their `affection` for you.

In true Bayesian fashion, you assign the following priors for the cats' `affection`:
- Mochi: `Beta(2, 3)`. She is very shy, so you expect her not to like you that much.
- Momo: `Beta(4, 2)`. Momo seems like a more sociable cat, so you have high expectations for him.

Finally, you collect data by trying to pet them. Mochi lets you pet her 3 out of 5 attempts, and Momo does so for 4 out of 11 attempts.
"""

# ╔═╡ d5f55a6e-c482-41be-ac32-2a3b57bb711e
md"""
!!! questions
	1. What is the expected `affection` for both cats?
	2. What is the chance Mochi likes you most?
	3. Visualise the prior and posterior `affection` for both cats. Where on the plots are the MLE and MAP estimates situated?
"""

# ╔═╡ 5642edf0-1187-4f46-ad4c-555c1f8256d3
md"### 1"

# ╔═╡ 9d129bf1-3527-420f-ac6e-296f9ddf50d2
@model function kitkat(a, b, n)
	affection ~ Beta(a, b)
	allowed_pets ~ Binomial(n, affection)
end

# ╔═╡ 7efedd08-ca9b-4902-b609-1f61d5aa5528
mochimodel = kitkat(2, 3, 5) | (allowed_pets = 4,)

# ╔═╡ 9d436e87-64fe-4913-8265-c7ad5081ca32
mochichain = sample(mochimodel, NUTS(), 2000)

# ╔═╡ 7f1cbb10-922a-49b0-a96c-d31992dba412
plot(mochichain)

# ╔═╡ c8276294-f16c-4d58-8f20-599834326e96
momomodel = kitkat(4, 2, 11) | (allowed_pets = 4,)

# ╔═╡ f4429d42-1b2e-462c-a5c6-a9b7e67f9832
momochain = sample(momomodel, NUTS(), 2000)

# ╔═╡ de595c61-12f2-44f2-815a-a54d70ccedbe
plot(momochain)

# ╔═╡ a90a2915-d260-481b-a8fa-c57f0654d85b
mean(mochichain[:affection])

# ╔═╡ 0ae66724-2188-40cd-8626-f757ca9f3d17
mean(momochain[:affection])

# ╔═╡ 7424fc21-d04b-423e-a344-11e04b81f803
md"### 2"

# ╔═╡ 22502264-7cad-4b0e-8af6-377b90e7eeec
mean(mochichain[:affection] .> momochain[:affection])

# ╔═╡ b54094be-b15e-4f62-a9d5-908ad85d5b09
md"### 3"

# ╔═╡ e03f656a-cd53-47ba-a77c-8c65a57c7778
begin
	plot(Beta(2, 3), title = "Mochi affection", label = "Prior")
	histogram!(mochichain[:affection], normalize = :pdf, alpha = 0.5, label = "Posterior")
	# the following 2 can also be answered without plotting
	vline!([median(mochichain[:affection])], label = "MAP", linewidth = 2)
	vline!([3/5], label = "MLE", linewidth = 2)
end

# ╔═╡ f2d68649-6012-49fe-a11e-7d735d28bb2d
begin
	plot(Beta(4, 2), title = "Momo affection", label = "Prior")
	histogram!(momochain[:affection], normalize = :pdf, alpha = 0.5, label = "Posterior")
	vline!([median(momochain[:affection])], label = "MAP", linewidth = 2)
	vline!([4/11], label = "MLE", linewidth = 2)
end

# ╔═╡ e9e09be2-62bc-4372-bdd4-3b78820188d7
md"## 3: Dice upon dice"

# ╔═╡ ff593160-9622-4994-b0bb-d277a9ae7aba
md"""
Your encounter with Carl the Chimera takes a tragic turn after failing to kill him in time. Rather than killing you, he **steals your magic** and sends you out of his lair laughing.

Now the only spell you have left to defend yourself is a crummy beginner spell: **Air-great-cubicuboctahedron**. You go to look up how much damage it deals again, but find someone spilled pepsi over the rulebook and the spell's information has become illegible.

All you remember is that it let you throw **3 dice**, they were either **4-sided, 6-sided or 8-sided**, and it usually dealt, like, **10 damage total**.
"""

# ╔═╡ 26d34e30-725d-4333-b887-ae72d242b88c
md"""
!!! question
	- Between the 3 possibilities for the amount of dice sides, which one is the most probable?
"""

# ╔═╡ 1a85b841-58bc-4ff0-ac05-b82dc6067103
md"""
!!! tip
	The `Categorical` and `Dirac` distributions can be useful for this question, though you can also solve it without them.
"""

# ╔═╡ 222f13eb-b7f2-470d-85be-2e1a2501cbfc
@model function dicing()
	num_sides ~ Categorical([0, 0, 0, 1/3, 0, 1/3, 0, 1/3]) 
		# or DiscreteUniform(4, 8)
	rolls = zeros(3)
	for i in eachindex(rolls)
		rolls[i] ~ DiscreteUniform(1, num_sides)
	end
	total ~ Dirac(sum(rolls))
end

# ╔═╡ 7e4aed1b-94cd-4282-bb4a-0014f1513441
dicemodel = dicing() | (total = 10,);

# ╔═╡ 3746762e-ae34-4a1c-897f-898f96d73b4c
dicechain = sample(dicemodel, PG(10), 2_000)

# ╔═╡ 6200e0eb-2044-450e-8da2-1f8fd0695cd6
plot(dicechain)

# ╔═╡ f0f33581-f31d-4d59-99be-0ffd04fe8af7
dicechain[:num_sides] |> median

# ╔═╡ cda36a5a-499b-4b07-9b09-cdd44a19c3c2
md"The dice being 6-sided is most probable."

# ╔═╡ 30088664-5157-4d99-8584-7a42d0acdfb8
md"## 4: Circleference"

# ╔═╡ 7dd2c189-79c0-4d29-9e17-9c24a78b5791
md"""
Given three (noisy) points $P_1=(x_1,y_1)$, $P_1=(x_2,y_2)$ and $P_3=(x_3,y_3)$, you want to infer the corresponding circle.

You can assume that the circle center can appear anywhere in the $[-20, 20]\times [-20, 20]$ square and the radius is between 0 and 50. Points are sampled randomly on the circle and have a slight amount of Gaussian noise ($\sigma=0.25$ works well).
"""

# ╔═╡ dcbf405c-786c-4226-b35c-dc718452bb61
md"""
!!! questions
	- Write a small probabilistic program that can infer the center and radius of the circle.
	- What does the inferred circle look like if you condition on only one or two of the circle points?
"""

# ╔═╡ 791ff4ed-c9d1-48e2-9dd0-4bf3979c6167
x1, y1 = 18.0, 2.1

# ╔═╡ cb189957-f9d4-480b-a492-92cfc2a8c2aa
x2, y2 = -7.3, 8.1

# ╔═╡ 81eca3b3-12d9-43d7-af14-e9aeb73f2471
x3, y3 = -13.0, -23.0

# ╔═╡ 28d28034-e999-4e34-b6b2-63c762094c59
begin
	
	function plotcircle!(p, R, xC, yC; dθ=0.01)
		θ = 0:dθ:2pi+0.1
		plot!(p, xC .+ R .* cos.(θ), yC .+ R .* sin.(θ), label="", alpha=0.5, color=:blue)
		return p
	end

	function plotsample(R=missing, xC=missing, yC=missing; kwargs...)
		p = plot(xlab="x", ylab="y", aspect_ratio=:equal;
				xlims=[-40, 40], ylims=[-40, 40], kwargs...)

		scatter!([x1], [y1], label="P1")
		scatter!([x2], [y2], label="P2")
		scatter!([x3], [y3], label="P3")
		ismissing(R) || plotcircle!(p, R, xC, yC; dθ=0.1)
		return p
	end

	function plotsample!(p, R=missing, xC=missing, yC=missing)
		scatter!([x1], [y1], label=false)
		scatter!([x2], [y2], label=false)
		scatter!([x3], [y3], label=false)
		ismissing(R) || plotcircle!(p, R, xC, yC; dθ=0.1)
	end

end

# ╔═╡ 4e96908a-4fc9-429d-bf37-7a569194a038
scatter([x1, x2, x3], [y1, y2, y3])

# ╔═╡ 84d27c98-9513-4ae3-8101-621c083a1b01
@model function circle(σ=0.25)
	# generate a circle center
	xC ~ Uniform(-20, 20)
	yC ~ Uniform(-20, 20)
	
	# generate a radius
	R ~ Uniform(0, 50)
	
	# three random points in polar coordinates
	θ1 ~ Uniform(0, 2pi)
	θ2 ~ Uniform(0, 2pi)
	θ3 ~ Uniform(0, 2pi)
	
	# P1
	x1 ~ Normal(xC + R * cos(θ1), σ)
	y1 ~ Normal(yC + R * sin(θ1), σ)
	# P2
	x2 ~ Normal(xC + R * cos(θ2), σ)
	y2 ~ Normal(yC + R * sin(θ2), σ)
	# P3
	x3 ~ Normal(xC + R * cos(θ3), σ)
	y3 ~ Normal(yC + R * sin(θ3), σ)
end

# ╔═╡ acc1d44a-f57e-4d2f-b33a-92cbddb67664
md"### All points"

# ╔═╡ 543c40d5-e8a7-492d-a0b8-e7e73e5953e2
circlemodel = circle() | (x1=x1, y1=y1, x2=x2, y2=y2, x3=x3, y3=y3);

# ╔═╡ 22173937-25a1-4ec0-877d-f9669653e43e
circlechain = sample(circlemodel, NUTS(), 2000)

# ╔═╡ 0b59a6bc-a886-49fd-a3f5-9b20e60882b9
plot(circlechain)

# ╔═╡ 0d04b0e8-3d1a-4281-b175-570148569ef2
begin
	p = plot(xlab="x", ylab="y", aspect_ratio=:equal; 
	xlims=[-40, 40], ylims=[-40, 40], title="Samples of P(circle|P1,P2,P3)")
	for i in 1:100
		plotsample!(p, circlechain[:R][i], circlechain[:xC][i], circlechain[:yC][i])
	end
	p
end

# ╔═╡ 0a140d2a-a24b-48df-9af8-7fa5d586a26f
md"### Some points"

# ╔═╡ b7649523-fd40-4fe3-8d86-fc2cb2c8c488
circle1 = circle() | (x1=x1, y1=y1)

# ╔═╡ b7ba8f91-f440-4166-af9e-11fcb5f1755f
circle2 = circle1 | (x2=x2, y2=y2)

# ╔═╡ 3b5f3623-3d00-4ba4-9182-5bddacc56567
chain1 = sample(circle1, NUTS(), 100);

# ╔═╡ c11452bc-f404-4e51-8f68-292f5b538c88
chain2 = sample(circle2, NUTS(), 100);

# ╔═╡ c1174a5c-a6a1-46a0-96be-6997e3201dfc
begin
	p1 = plot(xlab="x", ylab="y", aspect_ratio=:equal; 
	xlims=[-40, 40], ylims=[-40, 40], title="Samples of P(circle|P1)")
	for i in 1:100
		plotsample!(p1, chain1[:R][i], chain1[:xC][i], chain1[:yC][i])
	end
	p1
end

# ╔═╡ 793fc102-fc64-4041-a04e-e0b1b0741437
begin
	p2 = plot(xlab="x", ylab="y", aspect_ratio=:equal; 
	xlims=[-40, 40], ylims=[-40, 40], title="Samples of P(circle|P1,P2)")
	for i in 1:100
		plotsample!(p2, chain2[:R][i], chain2[:xC][i], chain2[:yC][i])
	end
	p2
end

# ╔═╡ Cell order:
# ╟─41dc8060-cf5e-11ef-26f9-892577e77af0
# ╠═94c6f31d-1a43-4221-b60c-1fa0ef8738b8
# ╠═45bc5b66-c81b-4afb-8a7e-51aff9609c62
# ╟─ed7c5547-bd82-4ca7-bf51-dd1c201d88af
# ╟─a3c26839-2acc-49d8-98c4-3a7142dd6512
# ╟─dac55632-6e17-479e-8090-5cf8eaa67dad
# ╠═122b1b8f-8ba4-4d2d-9ffb-f640f04889bb
# ╠═93411707-e019-4348-9728-31cc0cae6040
# ╠═a5bbb8ce-d134-4d7e-8ad7-54c2b06210ca
# ╠═126d3954-c77d-4c98-abe0-fd87d14e6265
# ╠═15fe25d9-5314-4056-8021-cf259ba27c94
# ╠═2a142576-aec2-4311-9c0e-cb68665d59f6
# ╠═0130c903-e3dd-415e-aa17-1b705f9e4ccc
# ╠═11676f8a-d4fe-4999-9ea2-486fca461f65
# ╟─ca8ec759-4a4f-4401-a747-aa904d3ee7b3
# ╟─dedbd315-dce1-47f0-bc49-47325cd5b170
# ╟─db4af5b9-d424-458d-aeaa-f803dd06a222
# ╟─d5f55a6e-c482-41be-ac32-2a3b57bb711e
# ╟─5642edf0-1187-4f46-ad4c-555c1f8256d3
# ╠═9d129bf1-3527-420f-ac6e-296f9ddf50d2
# ╠═7efedd08-ca9b-4902-b609-1f61d5aa5528
# ╠═9d436e87-64fe-4913-8265-c7ad5081ca32
# ╠═7f1cbb10-922a-49b0-a96c-d31992dba412
# ╠═c8276294-f16c-4d58-8f20-599834326e96
# ╠═f4429d42-1b2e-462c-a5c6-a9b7e67f9832
# ╠═de595c61-12f2-44f2-815a-a54d70ccedbe
# ╠═a90a2915-d260-481b-a8fa-c57f0654d85b
# ╠═0ae66724-2188-40cd-8626-f757ca9f3d17
# ╟─7424fc21-d04b-423e-a344-11e04b81f803
# ╠═22502264-7cad-4b0e-8af6-377b90e7eeec
# ╟─b54094be-b15e-4f62-a9d5-908ad85d5b09
# ╠═e03f656a-cd53-47ba-a77c-8c65a57c7778
# ╠═f2d68649-6012-49fe-a11e-7d735d28bb2d
# ╟─e9e09be2-62bc-4372-bdd4-3b78820188d7
# ╟─ff593160-9622-4994-b0bb-d277a9ae7aba
# ╟─26d34e30-725d-4333-b887-ae72d242b88c
# ╟─1a85b841-58bc-4ff0-ac05-b82dc6067103
# ╠═222f13eb-b7f2-470d-85be-2e1a2501cbfc
# ╠═7e4aed1b-94cd-4282-bb4a-0014f1513441
# ╠═3746762e-ae34-4a1c-897f-898f96d73b4c
# ╠═6200e0eb-2044-450e-8da2-1f8fd0695cd6
# ╠═f0f33581-f31d-4d59-99be-0ffd04fe8af7
# ╟─cda36a5a-499b-4b07-9b09-cdd44a19c3c2
# ╟─30088664-5157-4d99-8584-7a42d0acdfb8
# ╟─7dd2c189-79c0-4d29-9e17-9c24a78b5791
# ╟─dcbf405c-786c-4226-b35c-dc718452bb61
# ╟─28d28034-e999-4e34-b6b2-63c762094c59
# ╠═791ff4ed-c9d1-48e2-9dd0-4bf3979c6167
# ╠═cb189957-f9d4-480b-a492-92cfc2a8c2aa
# ╠═81eca3b3-12d9-43d7-af14-e9aeb73f2471
# ╠═4e96908a-4fc9-429d-bf37-7a569194a038
# ╠═84d27c98-9513-4ae3-8101-621c083a1b01
# ╟─acc1d44a-f57e-4d2f-b33a-92cbddb67664
# ╠═543c40d5-e8a7-492d-a0b8-e7e73e5953e2
# ╠═22173937-25a1-4ec0-877d-f9669653e43e
# ╠═0b59a6bc-a886-49fd-a3f5-9b20e60882b9
# ╠═0d04b0e8-3d1a-4281-b175-570148569ef2
# ╟─0a140d2a-a24b-48df-9af8-7fa5d586a26f
# ╠═b7649523-fd40-4fe3-8d86-fc2cb2c8c488
# ╠═b7ba8f91-f440-4166-af9e-11fcb5f1755f
# ╠═3b5f3623-3d00-4ba4-9182-5bddacc56567
# ╠═c11452bc-f404-4e51-8f68-292f5b538c88
# ╠═c1174a5c-a6a1-46a0-96be-6997e3201dfc
# ╠═793fc102-fc64-4041-a04e-e0b1b0741437
