### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 94c6f31d-1a43-4221-b60c-1fa0ef8738b8
using Pkg; Pkg.activate("..")

# ╔═╡ 45bc5b66-c81b-4afb-8a7e-51aff9609c62
using Turing, StatsPlots

# ╔═╡ 9daad9dc-b092-4845-a276-ed3a24249924
using PlutoUI; TableOfContents()

# ╔═╡ 41dc8060-cf5e-11ef-26f9-892577e77af0
md"# Inference notebook #2: Basics"

# ╔═╡ 299ba93b-0fc0-4bb3-9a2c-a571ce571f1b
md"## 1: Mole burrow"

# ╔═╡ 0d068aa6-f7c9-4a75-8382-cfbc903efc5b
md"""
Consider a mole's underground tunnel network of length `X` (in m). Now and then the mole makes a new molehill somewhere randomly above its tunnel, the locations of which we denote `Y`.
"""

# ╔═╡ 5d4e269d-bd5a-422e-9547-099e04b8eedf
md"""
                  Y1                                Y2   Y3
                  /\ 								/\   /\              <- molehills
    =================================================================    <- ground
                   |                                 |    |
             ------------------------------------------------            <- tunnel
             0                                              X
"""

# ╔═╡ 7abeb4d3-4648-4d8b-9383-de5213c8cc41
md"""
We're no mole experts, but for our prior information we can suppose a mole would not make a tunnel much longer than a few 100 m (probably a lot shorter).

We can formulate this as `X ~ Exponential(100)` and `Y ~ Uniform(0, X)`.

!!! questions
	1. Plot the prior of `X`. Is it diffuse or informative?
	1. Estimate `E[Y]`.
	1. Estimate `E[X|Y = 3]` and compare it with the prior expected value `E[X]`.
	1. Plot the histogram of `X` given `Y = 3.0`.
	1. You come back a day later and find even more molehills! You now measure the following values for `Y`: `[3.0, 1.5, 0.9, 5.7]`. How long do you estimate the tunnel given this extra information?
"""

# ╔═╡ 47aa2304-d312-40e9-b9c6-7c79a7d64de4
md"### 1: Prior plot"

# ╔═╡ 78c9ce62-e375-48ac-8083-55ee085c61de
X_prior = missing

# ╔═╡ c64b346a-059d-4bbe-b166-80116260e19b
missing # plot above

# ╔═╡ 78eb7779-f182-4419-b5d8-79a2f5c5d6da
md"### 2: Unconditional expected value of Y"

# ╔═╡ 3decb2ec-210a-4b2f-842d-6fd40dd3f77b
@model function mole()
	X ~ X_prior
	Y ~ missing
	return Y
end

# ╔═╡ 76dd814d-0d9b-4f7e-aff8-990da57d052b
molemodel = mole()

# ╔═╡ 61b8ef63-a613-4dcb-8f7a-936e0d59b862
Y_samples = missing

# ╔═╡ 4e8b9000-cb61-4f01-9ba9-17276ad0335e
E_Y = missing

# ╔═╡ 8777b133-7d7c-4a85-b89c-2f00093e9984
md"### 3: Conditional expected value of X"

# ╔═╡ 9695ee7e-359b-489c-962b-1bf84b052371
cond_mole = missing

# ╔═╡ 014538da-b5ef-41c8-b799-2c000b4c9134
molechain = missing

# ╔═╡ 977d0406-f154-4b87-b9c6-cafe4809a4a6
missing # plot molechain

# ╔═╡ af4811ce-c6e7-458a-8f89-0562355b3eb0
X_samplescondY = missing

# ╔═╡ fd540765-95e5-4071-81f7-e689b06cad0c
E_XcondY = missing

# ╔═╡ 5521933a-a42e-4a67-94b9-84eab52ddf07
E_X = missing

# ╔═╡ 9403c355-dce8-4aa9-80a3-824ea6193905
md"### 4: Conditional distribution of X"

# ╔═╡ bedd99ca-4285-4c3d-87e2-c22d414f8f07
# histogram

# ╔═╡ 21743f62-c7f5-4171-8277-0667edb71c56
md"### 5: Conditional distribution of X (with more data)"

# ╔═╡ b2b16c34-e12a-4bea-8098-313d01913bbf
@model function mole2()
	X ~ X_prior
	Ys = missing
	for i in missing
		Ys[i] ~ missing
	end
end

# ╔═╡ 26a245d2-c5b6-484d-bdf7-541948ff06ad
Y_obs = [3.0, 1.5, 0.9, 5.7]

# ╔═╡ dc6aadff-7cff-4a9a-a967-32cd322b2696
missing # (some lines of code)

# ╔═╡ ea2235f0-0cfc-4603-b027-d0a1bc6fa2c1
# I estimate the tunnel to be about (MISSING) long

# ╔═╡ 951d0913-1a52-4d5b-b5bb-168487e50ab2
md"## 2: Potatoes"

# ╔═╡ 49035a16-c419-4531-b157-a5ab357b44fe
md"""
Consider a number of potatoes `N` each with an average weight `W`. You weigh them together on an old balance to get an estimate of their total weight `T`.

Suppose the following priors: `N ~ Poisson(10)`, `W ~ Uniform(150, 250)`  and the following relationship between expected - and actual outcome: `T ~ Normal(N*W, 50)`.

!!! questions
	1. Plot a histogram of `N` given `T = 1200`.
	1. Estimate `P(N > 6, W > 175 | T = 1200)`.
	1. Estimate `P(N = 5 | T = 1200, W = 220)`.
"""

# ╔═╡ 5f8322ca-0ded-4750-8ad6-e66e8280daca
md"### 1: Conditional distribution of N"

# ╔═╡ ec479be9-0d8a-4c8d-97c9-6f64d861924c
@model function potatoes()
	N ~ missing
	W ~ missing
	T ~ missing
end

# ╔═╡ 438893d1-3107-493c-9299-d8813658a698
potato_model = missing

# ╔═╡ dd7d84b0-1f20-4478-96ba-b3f4e9df1d2c
potato_cond = missing

# ╔═╡ accca8d7-8822-49ae-8588-7259da82e37a
potato_chain = missing
	# be careful with your choice of sampling algorithm! are all priors continuous?

# ╔═╡ 519a952b-f040-4588-8acd-1cfaab88c4da
# plot chain

# ╔═╡ 57adc714-cdd1-46cf-a4ea-d0fb04924e1e
# histogram of N

# ╔═╡ 0fddec58-67d9-4dec-a72d-445675fa46a6
md"### 2: Probability"

# ╔═╡ b77d8796-9065-4c61-be21-3dc292b94492
p_potato1 = missing

# ╔═╡ 5bde987a-3f5c-4ff7-96a0-9e174f73fdfb
md"### 3: Probability (with more data)"

# ╔═╡ abffd02e-b917-4e9e-a94d-de197e1b0fc3
potato_cond2 = missing # alternative conditioned model now also given `W = 220`

# ╔═╡ 5b30d91a-7039-4d2a-bcf7-3a990622fd2a
pota2_chain = missing

# ╔═╡ 27939d84-d6e1-47b2-a939-5038c6bdcb48
# plot chain

# ╔═╡ 744cf31f-89c0-44e8-a6fb-0f4dcea04e5d
p_potato2 = missing

# ╔═╡ ed7c5547-bd82-4ca7-bf51-dd1c201d88af
md"## 3: Lights out"

# ╔═╡ a3c26839-2acc-49d8-98c4-3a7142dd6512
md"""
You use 4 of the same LED light in your room. Define:
-  $μ$ the **average** lifespan of your LED lights (in khr or 1000 hours) 
-  $Lᵢ$ the lifespan of the `i`-th LED light.

Assume that `μ ~ LogNormal(log(40), 0.5)`.
"""

# ╔═╡ dac55632-6e17-479e-8090-5cf8eaa67dad
md"""
!!! questions
	1. What is `E[μ]` (given no information about `Lᵢ`)?
	1. What is a sensible distribution for `Lᵢ`? (requires no code)
	1. What is `E[μ | L = [16, 20, 23, 41]]`?
	1. 🌟🌟🌟 (EXTRA DIFFICULT BONUS QUESTION): After 30 khr, two lights have died: one at 16 khr and one at 20 khr. The two other lights are still working. What is the expected value of `μ` given this information? 
"""

# ╔═╡ 7df97000-5cf0-4a29-b3c3-f867403f4318
md"### 1: Unconditional expected value"

# ╔═╡ 0fe475dd-411e-474e-9524-ef9fcceed7af
lights_prior = missing

# ╔═╡ db77b90b-33b5-4d70-8229-b4c7d74fa96c
E_mu = missing

# ╔═╡ 308cb1f5-e26b-43b2-be4d-9eabd68b3670
md"### 2: Choice of distribution"

# ╔═╡ 738b4097-56e1-4921-9195-2d535c78ae73
md"""
!!! note
	See theory p. 96 for an overview of elementary distributions.
"""

# ╔═╡ d797cb9b-430b-4d1d-8ded-21959a7cb3a9
# A sensible distribution for the lifespan of a LED-light is (missing)

# ╔═╡ cd37e5fc-ae9a-429d-844d-4b450b187b5e
md"### 3: Conditional expected value"

# ╔═╡ 126d3954-c77d-4c98-abe0-fd87d14e6265
@model function lights()
	μ ~ lights_prior
	lifespans = missing
	for i in missing
		lifespans[i] ~ missing
	end
end

# ╔═╡ 37bf6ca5-e42b-406e-94d6-bcd1e706e2bd
L_obs = [16, 20, 23, 41]

# ╔═╡ c49955c7-7082-4261-b4db-0056b5c637f1
E_mu_cond = missing

# ╔═╡ 9a168680-9d96-466f-ba75-122d6a391501
md"### 4 🌟🌟🌟: Working with censored data"

# ╔═╡ ab6331a3-676e-45e6-b5be-e9c4154ea071
md"""
!!! note
	This question is way above the exam's difficulty level. Don't feel bad if you can't find the answer right away!
"""

# ╔═╡ 1f35d962-a249-4be7-9a96-17eb83fca7d8
md"""
!!! hint
	You can model the number of lights that still work as a `Binomial` distribution, the success rate of which depends on `μ`.
"""

# ╔═╡ 8fc58fa4-b005-4f32-9eae-a8143582a1ae
@model function lights_censored(time_observed)
	missing
end

# ╔═╡ 410521d8-f767-4c4e-b19b-25cf31ec0f36
E_mu_cond🌟 = missing

# ╔═╡ 9dc0456b-7fd2-4120-8f9e-3de1984ff516
md"## 4: Fish"

# ╔═╡ 225cd579-1e0d-4680-8d3f-5a737a656eb8
md"""
There are two populations of fish living in the same pond. Let `fs1` be the fraction of fish belonging to species 1, `L1` the length of a fish of species 1 and `L2` the length of a fish of species 2.

Assume:
- You have no prior information about `fs1` except that it logically needs to be in `[0, 1]`.
- `L1 ~ Normal(90, 15)`.
- `L2 ~ Normal(60, 10)`.
"""

# ╔═╡ 3ffac5ca-3635-4aa1-bab7-7c28e7a801cb
md"""
!!! questions
	1. If `fs1 = 0.3`, what is the prior distribution of the lengths of **all** fish in the pond? Make a plot.
	1. Estimate `fs1` if you observe fish of the following lengths:
       `[94.0, 88.7, 89.6, 69.8, 52.8, 84.0, 89.3, 66.4, 95.1, 81.6]`.
	1. 🌟(BONUS QUESTION): What is the chance fish 4 belongs to species 1?
"""

# ╔═╡ b477b212-83a2-42f0-a616-52516e152d48
md"### 1: Prior distribution given `fs1`"

# ╔═╡ 23056f1e-128e-463e-a80f-56299397022e
md"""
!!! hint
	The distribution of fish lengths can be modelled as a `MixtureModel`.
"""

# ╔═╡ 3cd0c888-9f11-47f1-a293-b96aa80ea3b0
lengthdist = missing

# ╔═╡ e23756e0-7688-4e3b-ba8d-a827e621e862
missing # plot

# ╔═╡ 8e030a06-5104-4c5b-b1f2-f86464e66502
md"### 2: Conditional expected value"

# ╔═╡ 6ebf3a16-0e6a-491f-8280-b4327ed52cf0
len_obs = [94.0, 88.7, 89.6, 69.8, 52.8, 84.0, 89.3, 66.4, 95.1, 81.6]

# ╔═╡ b15de91a-fc48-4cdc-a35f-6453a9a59982
@model function fishmixture()
	fs1 ~ missing # fraction of species 1
	fishlendist = missing
	
	fishlens = missing
	for i in missing
		fishlens[i] ~ fishlendist
	end
end

# ╔═╡ 94a9e163-f7c0-4d22-8d37-5641bc49eb6b
fishmodel = missing

# ╔═╡ 40d460d2-1f2b-4d33-ba74-9a9d1e48f9ec
fishchain = missing # don't forget to plot to check convergence

# ╔═╡ 819ed364-4bc9-43ec-aa50-b965c8f1c826
fs1_est = missing

# ╔═╡ 952a941a-8703-45a3-aac1-a290a181e8c5
md"### 3🌟: Conditional expected value (spicy)"

# ╔═╡ 8629d049-b9fd-4e9e-9b55-401a3069e956
@model function fishmixture🌟()
	fs1 ~ missing # fraction of species 1
	
	fishlens = missing
	isspecies1 = missing
	for i in missing
		isspecies1[i] ~ missing
		if isspecies1[i] == 1.0
			fishlens[i] ~ missing
		else
			fishlens[i] ~ missing
		end
	end
end

# ╔═╡ 483cbe4d-64d3-4b16-b6fc-e97b21f174a6
p_fish4_is_species1 = missing

# ╔═╡ 30088664-5157-4d99-8584-7a42d0acdfb8
md"## 5: Circleference"

# ╔═╡ 7dd2c189-79c0-4d29-9e17-9c24a78b5791
md"""
Given three (noisy) points $P_1=(x_1,y_1)$, $P_1=(x_2,y_2)$ and $P_3=(x_3,y_3)$, you want to infer the corresponding circle.

You can assume that the circle center can appear anywhere in the $[-20, 20]\times [-20, 20]$ square and the radius is between 0 and 50. Points are sampled randomly on the circle and have a slight amount of Gaussian noise ($\sigma=0.25$ works well).
"""

# ╔═╡ dcbf405c-786c-4226-b35c-dc718452bb61
md"""
!!! questions
	1. Write a small probabilistic program that can infer the center and radius of the circle.
	1. What does the inferred circle look like if you condition on only one or two of the circle points?
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

# ╔═╡ 7eedb74d-eee1-4cf0-b2bf-5febf474edd2
md"### 1: All points"

# ╔═╡ 84d27c98-9513-4ae3-8101-621c083a1b01
@model function circle(σ=0.25)
	# generate a circle center
	xC ~ missing
	yC ~ missing
	
	# generate a radius
	R ~ missing
	
	# three random points in polar coordinates
	θ1 ~ missing
	θ2 ~ missing
	θ3 ~ missing
	
	# P1
	x1 ~ Normal(missing, σ)
	y1 ~ Normal(missing, σ)
	# P2
	x2 ~ Normal(missing, σ)
	y2 ~ Normal(missing, σ)
	# P3
	x3 ~ Normal(missing, σ)
	y3 ~ Normal(missing, σ)
end

# ╔═╡ 543c40d5-e8a7-492d-a0b8-e7e73e5953e2
circlemodel = circle() | (x1=x1, y1=y1, x2=x2, y2=y2, x3=x3, y3=y3);

# ╔═╡ 25ea29e7-b391-4bd1-bdbc-1957adb8c993
circlechain = missing

# ╔═╡ 0d04b0e8-3d1a-4281-b175-570148569ef2
begin
	p = plot(
		xlab="x", ylab="y", aspect_ratio=:equal, xlims=[-40, 40], ylims=[-40, 40],
		title="Samples of P(circle|P1,P2,P3)"
	)
	for i in 1:100
		plotsample!(p, circlechain[:R][i], circlechain[:xC][i], circlechain[:yC][i])
	end
	p
end

# ╔═╡ 0a140d2a-a24b-48df-9af8-7fa5d586a26f
md"### 2: Some points"

# ╔═╡ b7649523-fd40-4fe3-8d86-fc2cb2c8c488
circle1 = missing # given one point

# ╔═╡ b7ba8f91-f440-4166-af9e-11fcb5f1755f
circle2 = missing # given two points

# ╔═╡ 3b5f3623-3d00-4ba4-9182-5bddacc56567
chain1 = missing

# ╔═╡ c11452bc-f404-4e51-8f68-292f5b538c88
chain2 = missing

# ╔═╡ c1174a5c-a6a1-46a0-96be-6997e3201dfc
begin
	p1 = plot(
		xlab="x", ylab="y", aspect_ratio=:equal,
		xlims=[-40, 40], ylims=[-40, 40], title="Samples of P(circle|P1)"
	)
	for i in 1:100
		plotsample!(p1, chain1[:R][i], chain1[:xC][i], chain1[:yC][i])
	end
	p1
end

# ╔═╡ 793fc102-fc64-4041-a04e-e0b1b0741437
begin
	p2 = plot(
		xlab="x", ylab="y", aspect_ratio=:equal, 
		xlims=[-40, 40], ylims=[-40, 40], title="Samples of P(circle|P1,P2)"
	)
	for i in 1:100
		plotsample!(p2, chain2[:R][i], chain2[:xC][i], chain2[:yC][i])
	end
	p2
end

# ╔═╡ Cell order:
# ╟─41dc8060-cf5e-11ef-26f9-892577e77af0
# ╠═94c6f31d-1a43-4221-b60c-1fa0ef8738b8
# ╠═45bc5b66-c81b-4afb-8a7e-51aff9609c62
# ╠═9daad9dc-b092-4845-a276-ed3a24249924
# ╟─299ba93b-0fc0-4bb3-9a2c-a571ce571f1b
# ╟─0d068aa6-f7c9-4a75-8382-cfbc903efc5b
# ╟─5d4e269d-bd5a-422e-9547-099e04b8eedf
# ╟─7abeb4d3-4648-4d8b-9383-de5213c8cc41
# ╟─47aa2304-d312-40e9-b9c6-7c79a7d64de4
# ╠═78c9ce62-e375-48ac-8083-55ee085c61de
# ╠═c64b346a-059d-4bbe-b166-80116260e19b
# ╟─78eb7779-f182-4419-b5d8-79a2f5c5d6da
# ╠═3decb2ec-210a-4b2f-842d-6fd40dd3f77b
# ╠═76dd814d-0d9b-4f7e-aff8-990da57d052b
# ╠═61b8ef63-a613-4dcb-8f7a-936e0d59b862
# ╠═4e8b9000-cb61-4f01-9ba9-17276ad0335e
# ╟─8777b133-7d7c-4a85-b89c-2f00093e9984
# ╠═9695ee7e-359b-489c-962b-1bf84b052371
# ╠═014538da-b5ef-41c8-b799-2c000b4c9134
# ╠═977d0406-f154-4b87-b9c6-cafe4809a4a6
# ╠═af4811ce-c6e7-458a-8f89-0562355b3eb0
# ╠═fd540765-95e5-4071-81f7-e689b06cad0c
# ╠═5521933a-a42e-4a67-94b9-84eab52ddf07
# ╟─9403c355-dce8-4aa9-80a3-824ea6193905
# ╠═bedd99ca-4285-4c3d-87e2-c22d414f8f07
# ╟─21743f62-c7f5-4171-8277-0667edb71c56
# ╠═b2b16c34-e12a-4bea-8098-313d01913bbf
# ╠═26a245d2-c5b6-484d-bdf7-541948ff06ad
# ╠═dc6aadff-7cff-4a9a-a967-32cd322b2696
# ╠═ea2235f0-0cfc-4603-b027-d0a1bc6fa2c1
# ╟─951d0913-1a52-4d5b-b5bb-168487e50ab2
# ╟─49035a16-c419-4531-b157-a5ab357b44fe
# ╟─5f8322ca-0ded-4750-8ad6-e66e8280daca
# ╠═ec479be9-0d8a-4c8d-97c9-6f64d861924c
# ╠═438893d1-3107-493c-9299-d8813658a698
# ╠═dd7d84b0-1f20-4478-96ba-b3f4e9df1d2c
# ╠═accca8d7-8822-49ae-8588-7259da82e37a
# ╠═519a952b-f040-4588-8acd-1cfaab88c4da
# ╠═57adc714-cdd1-46cf-a4ea-d0fb04924e1e
# ╟─0fddec58-67d9-4dec-a72d-445675fa46a6
# ╠═b77d8796-9065-4c61-be21-3dc292b94492
# ╟─5bde987a-3f5c-4ff7-96a0-9e174f73fdfb
# ╠═abffd02e-b917-4e9e-a94d-de197e1b0fc3
# ╠═5b30d91a-7039-4d2a-bcf7-3a990622fd2a
# ╠═27939d84-d6e1-47b2-a939-5038c6bdcb48
# ╠═744cf31f-89c0-44e8-a6fb-0f4dcea04e5d
# ╟─ed7c5547-bd82-4ca7-bf51-dd1c201d88af
# ╟─a3c26839-2acc-49d8-98c4-3a7142dd6512
# ╟─dac55632-6e17-479e-8090-5cf8eaa67dad
# ╟─7df97000-5cf0-4a29-b3c3-f867403f4318
# ╠═0fe475dd-411e-474e-9524-ef9fcceed7af
# ╠═db77b90b-33b5-4d70-8229-b4c7d74fa96c
# ╟─308cb1f5-e26b-43b2-be4d-9eabd68b3670
# ╟─738b4097-56e1-4921-9195-2d535c78ae73
# ╠═d797cb9b-430b-4d1d-8ded-21959a7cb3a9
# ╟─cd37e5fc-ae9a-429d-844d-4b450b187b5e
# ╠═126d3954-c77d-4c98-abe0-fd87d14e6265
# ╠═37bf6ca5-e42b-406e-94d6-bcd1e706e2bd
# ╠═c49955c7-7082-4261-b4db-0056b5c637f1
# ╟─9a168680-9d96-466f-ba75-122d6a391501
# ╟─ab6331a3-676e-45e6-b5be-e9c4154ea071
# ╟─1f35d962-a249-4be7-9a96-17eb83fca7d8
# ╠═8fc58fa4-b005-4f32-9eae-a8143582a1ae
# ╠═410521d8-f767-4c4e-b19b-25cf31ec0f36
# ╟─9dc0456b-7fd2-4120-8f9e-3de1984ff516
# ╟─225cd579-1e0d-4680-8d3f-5a737a656eb8
# ╟─3ffac5ca-3635-4aa1-bab7-7c28e7a801cb
# ╟─b477b212-83a2-42f0-a616-52516e152d48
# ╟─23056f1e-128e-463e-a80f-56299397022e
# ╠═3cd0c888-9f11-47f1-a293-b96aa80ea3b0
# ╠═e23756e0-7688-4e3b-ba8d-a827e621e862
# ╟─8e030a06-5104-4c5b-b1f2-f86464e66502
# ╠═6ebf3a16-0e6a-491f-8280-b4327ed52cf0
# ╠═b15de91a-fc48-4cdc-a35f-6453a9a59982
# ╠═94a9e163-f7c0-4d22-8d37-5641bc49eb6b
# ╠═40d460d2-1f2b-4d33-ba74-9a9d1e48f9ec
# ╠═819ed364-4bc9-43ec-aa50-b965c8f1c826
# ╟─952a941a-8703-45a3-aac1-a290a181e8c5
# ╠═8629d049-b9fd-4e9e-9b55-401a3069e956
# ╠═483cbe4d-64d3-4b16-b6fc-e97b21f174a6
# ╟─30088664-5157-4d99-8584-7a42d0acdfb8
# ╟─7dd2c189-79c0-4d29-9e17-9c24a78b5791
# ╟─dcbf405c-786c-4226-b35c-dc718452bb61
# ╟─28d28034-e999-4e34-b6b2-63c762094c59
# ╠═791ff4ed-c9d1-48e2-9dd0-4bf3979c6167
# ╠═cb189957-f9d4-480b-a492-92cfc2a8c2aa
# ╠═81eca3b3-12d9-43d7-af14-e9aeb73f2471
# ╠═4e96908a-4fc9-429d-bf37-7a569194a038
# ╟─7eedb74d-eee1-4cf0-b2bf-5febf474edd2
# ╠═84d27c98-9513-4ae3-8101-621c083a1b01
# ╠═543c40d5-e8a7-492d-a0b8-e7e73e5953e2
# ╠═25ea29e7-b391-4bd1-bdbc-1957adb8c993
# ╠═0d04b0e8-3d1a-4281-b175-570148569ef2
# ╟─0a140d2a-a24b-48df-9af8-7fa5d586a26f
# ╠═b7649523-fd40-4fe3-8d86-fc2cb2c8c488
# ╠═b7ba8f91-f440-4166-af9e-11fcb5f1755f
# ╠═3b5f3623-3d00-4ba4-9182-5bddacc56567
# ╠═c11452bc-f404-4e51-8f68-292f5b538c88
# ╠═c1174a5c-a6a1-46a0-96be-6997e3201dfc
# ╠═793fc102-fc64-4041-a04e-e0b1b0741437
