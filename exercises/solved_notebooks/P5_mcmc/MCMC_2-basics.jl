### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 94c6f31d-1a43-4221-b60c-1fa0ef8738b8
using Pkg; Pkg.activate("..")

# ╔═╡ 45bc5b66-c81b-4afb-8a7e-51aff9609c62
using Turing, StatsPlots

# ╔═╡ 8f624627-5a65-441d-989b-3fd7397e341c
using PlutoUI; TableOfContents()

# ╔═╡ 41dc8060-cf5e-11ef-26f9-892577e77af0
md"# Inference notebook #2: Basics"

# ╔═╡ 299ba93b-0fc0-4bb3-9a2c-a571ce571f1b
md"## 1: Mole burrow"

# ╔═╡ 6a7be039-c7f5-4cc0-8f21-beb8ac20a181
md"""
Consider a mole's underground tunnel network of length `X` (in m). Now and then the mole makes a new molehill somewhere randomly above its tunnel, the locations of which we denote `Y`.
"""

# ╔═╡ a068aad0-59c2-4b45-aa0f-ce73bd184bdc
md"""
                  Y1                                Y2   Y3
                  /\ 								/\   /\              <- molehills
    =================================================================    <- ground
                   |                                 |    |
             ------------------------------------------------            <- tunnel
             0                                              X
"""

# ╔═╡ c72f1e5e-e8f6-41c3-bf75-ff5c6f36cd05
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

# ╔═╡ b8c69b3b-d9e7-4237-9cbf-95590b2d5770
md"### 1: Prior plot"

# ╔═╡ 5192d03d-b9d1-4f03-aaf1-f287e044e904
X_prior = Exponential(100)

# ╔═╡ 78c9ce62-e375-48ac-8083-55ee085c61de
plot(X_prior) # The prior only incorporates the knowledge that a mole's tunnel is probably less long than a few 100 ms - this is a diffuse to weakly informative prior

# ╔═╡ 78eb7779-f182-4419-b5d8-79a2f5c5d6da
md"### 2: Unconditional expected value of Y"

# ╔═╡ 3decb2ec-210a-4b2f-842d-6fd40dd3f77b
@model function mole()
	X ~ X_prior
	Y ~ Uniform(0, X)
	return Y
end

# ╔═╡ 76dd814d-0d9b-4f7e-aff8-990da57d052b
molemodel = mole()

# ╔═╡ b5bd5114-0a0c-4aa3-af00-4f3d05edc5e3
Y_samples = [molemodel() for i in 1:2000];

# ╔═╡ 4e8b9000-cb61-4f01-9ba9-17276ad0335e
E_Y = mean(Y_samples)

# ╔═╡ 8777b133-7d7c-4a85-b89c-2f00093e9984
md"### 3: Conditional expected value of X"

# ╔═╡ 9695ee7e-359b-489c-962b-1bf84b052371
cond_mole = molemodel | (Y = 3,);

# ╔═╡ 014538da-b5ef-41c8-b799-2c000b4c9134
molechain = sample(cond_mole, NUTS(), 2000);

# ╔═╡ b02dc714-e2bb-4ae2-acf9-c37a4389f953
plot(molechain)

# ╔═╡ b79dafee-4f6f-4301-9267-1b9798c125ea
X_samplescondY = molechain[:X];

# ╔═╡ fd540765-95e5-4071-81f7-e689b06cad0c
E_XcondY = mean(molechain[:X])

# ╔═╡ 5521933a-a42e-4a67-94b9-84eab52ddf07
E_X = mean(X_prior)

# ╔═╡ ca3e730c-a940-4c76-93eb-70ae4aa0e008
md"### 4: Conditional distribution of X"

# ╔═╡ 159caa6a-2ebf-44bd-87a5-b4ab8b085354
histogram(molechain[:X])

# ╔═╡ eda00c08-de49-4d2d-acc4-ba6e21ff0b11
md"### 5: Conditional distribution of X (with more data)"

# ╔═╡ b2b16c34-e12a-4bea-8098-313d01913bbf
@model function mole2()
	X ~ X_prior
	Ys = zeros(4)
	for i in 1:length(Ys)
		Ys[i] ~ Uniform(0, X)
	end
end

# ╔═╡ 1a78d96a-fe57-42a4-9785-003266117ddb
Y_obs = [3.0, 1.5, 0.9, 5.7]

# ╔═╡ 2f5c14e9-709a-40ec-a6cb-ddd870cc1a60
mole_cond2 = mole2() | (Ys = Y_obs,)

# ╔═╡ 7f90e18d-9af8-42fa-984b-aaa7c8c458b5
molechain2 = sample(mole_cond2, NUTS(), 2000);

# ╔═╡ caf0db69-a83f-433a-b848-7d7d8c2fa25e
plot(molechain2)

# ╔═╡ 3c1a508b-6d2c-4507-ab9f-752ee93709c9
histogram(molechain2[:X]) # it's probably 5-10m

# ╔═╡ 951d0913-1a52-4d5b-b5bb-168487e50ab2
md"## 2: Potatoes"

# ╔═╡ a2cdbd3a-9a53-4c34-949e-4604ff1874b4
md"""
Consider a number of potatoes `N` each with an average weight `W`. You weigh them together on an old balance to get an estimate of their total weight `T`.

Suppose the following priors: `N ~ Poisson(10)`, `W ~ Uniform(150, 250)`  and the following relationship between expected - and actual outcome: `T ~ Normal(N*W, 50)`.

!!! questions
	1. Plot a histogram of `N` given `T = 1200`.
	1. Estimate `P(N > 6, W > 175 | T = 1200)`.
	1. Estimate `P(N = 5 | T = 1200, W = 220)`.
"""

# ╔═╡ d90c08dc-5842-4c84-9dd3-62dc22a6aab8
md"### 1: Conditional distribution of N"

# ╔═╡ ec479be9-0d8a-4c8d-97c9-6f64d861924c
@model function potatoes()
	N ~ Poisson(10)
	W ~ Uniform(150, 250)
	T ~ Normal(N*W, 50)
end

# ╔═╡ 4e92d825-481f-430a-94a4-fbdf31b679eb
potato_model = potatoes()

# ╔═╡ 0ced04f8-322a-42f6-a1b7-1b30faf9024b
potato_cond = potato_model | (T = 1200,)

# ╔═╡ e48a2d96-6d7b-4861-b315-0ae580479eda
potato_chain = sample(potato_cond, PG(10), 2000);

# ╔═╡ ceb2b8fd-aa47-443b-8e5c-74a48eda1174
plot(potato_chain)

# ╔═╡ 8c1eb826-6023-4852-aa45-789f6d4b7051
histogram(potato_chain[:N])

# ╔═╡ 5098d57a-cf21-40c0-9091-5233439eb9c4
md"### 2: Probability"

# ╔═╡ 9dd6fa35-37d5-4ab4-ac2c-f4eefdabadd2
p_potato1 = mean(potato_chain[:N] .> 6 .&& potato_chain[:W] .> 175)

# ╔═╡ 989b755b-2275-4aae-a3e2-3b107a72fc0b
md"### 3: Probability (with more data)"

# ╔═╡ 6e3d4f8f-1766-49a6-9e4a-610671d8aa63
potato_cond2 = potato_model | (T = 1200, W = 220)

# ╔═╡ 4cd99dce-12df-4e3a-aa91-ff00c61544a8
pota2_chain = sample(potato_cond2, PG(10), 2000)

# ╔═╡ 48ea9893-16f7-457e-93b0-9f168a2f72f4
plot(pota2_chain)

# ╔═╡ c3d36685-e1dd-4734-bee9-f7aeebed3473
p_potato2 = mean(pota2_chain[:N] .== 5)

# ╔═╡ ed7c5547-bd82-4ca7-bf51-dd1c201d88af
md"## 3: Lights out"

# ╔═╡ ed124916-331d-43b1-be44-e78c250d79c5
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
lights_prior = LogNormal(log(40), 0.5)

# ╔═╡ 382db9a6-edcb-497b-9535-b576cd6badb0
plot(lights_prior) # not asked but a visualisation can always be useful

# ╔═╡ fce3e814-97a7-4e12-a7cc-b5b086dc8340
E_mu = mean(lights_prior)

# ╔═╡ ca403113-1c49-4819-9a54-5d1a20fb03cb
md"### 2: Choice of distribution"

# ╔═╡ 9c744027-38da-4393-bc08-1feb0a9d8c94
md"""
!!! note
	See theory p. 96 for an overview of elementary distributions.
"""

# ╔═╡ 490cbc67-3e0c-4c79-a83d-52cf493854a2
md"""
The exponential distribution is often used to model the waiting time for an event. This makes it a natural fit for a lamp's lifespan, which is the waiting time until it breaks. We know it needs to have a mean value of μ, so `Exponential(μ)` is a good choice.

One could also argue for a LogNormal distribution with mean μ or a Normal distribution with mean μ restricted to only the positive values. Both would need a large variance to reflect the lack of additional information outside of the mean lifespan.
"""

# ╔═╡ cd37e5fc-ae9a-429d-844d-4b450b187b5e
md"### 3: Conditional expected value"

# ╔═╡ 126d3954-c77d-4c98-abe0-fd87d14e6265
@model function lights()
	μ ~ lights_prior
	lifespans = zeros(4)
	for i in 1:length(lifespans)
		lifespans[i] ~ Exponential(μ)
	end
end

# ╔═╡ 13ed391e-a6d4-4c43-9b69-f60154bb129a
L_obs = [16, 20, 23, 41]

# ╔═╡ 15fe25d9-5314-4056-8021-cf259ba27c94
lightmodel = lights() | (lifespans = L_obs,)

# ╔═╡ 2a142576-aec2-4311-9c0e-cb68665d59f6
lightschain = sample(lightmodel, NUTS(), 2000)

# ╔═╡ 0130c903-e3dd-415e-aa17-1b705f9e4ccc
plot(lightschain)

# ╔═╡ 11676f8a-d4fe-4999-9ea2-486fca461f65
E_mu_cond = mean(lightschain[:μ])

# ╔═╡ bd7c310c-888e-4380-b261-da2525ce9c11
md"### 4 🌟🌟🌟: Working with censored data"

# ╔═╡ 4c7ead00-134b-41d0-89f5-49c7ca14e188
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
@model function lights_censored(time_observed, n_working)
	μ ~ lights_prior
	
	lifespans = zeros(2)
	for i in 1:length(lifespans)
		lifespans[i] ~ Exponential(μ)
	end
	p_stillworking = 1 - cdf(Exponential(μ), time_observed)
	n_working ~ Binomial(n_working + length(lifespans), p_stillworking) 
		# or simply `2 ~ Binomial(4, p_stillworking)`
end

# ╔═╡ fe2958a7-e9dd-4eca-979d-a80df12f8735
lightmodel_cens = lights_censored(30, 2) | (lifespans = [16, 20],)

# ╔═╡ 8abafb6a-dc83-422c-82d1-a721a0e1eca0
lightschain_cens = sample(lightmodel_cens, MH(), 10_000)

# ╔═╡ 5ba2886c-b2e1-49f4-90c6-549acc808f77
plot(lightschain_cens)

# ╔═╡ ed022247-5959-481f-a81e-41e5ee5a1448
E_mu_cond🌟 = mean(lightschain_cens[:μ])

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

# ╔═╡ de3151b7-e599-4f9d-9ee1-969830bac0e9
md"### 1: Prior distribution given `fs1`"

# ╔═╡ 23056f1e-128e-463e-a80f-56299397022e
md"""
!!! hint
	The distribution of fish lengths can be modelled as a `MixtureModel`.
"""

# ╔═╡ 3cd0c888-9f11-47f1-a293-b96aa80ea3b0
lengthdist = MixtureModel([Normal(90, 15), Normal(60, 10)], [0.3, 0.7])

# ╔═╡ afcf1260-db4d-4ff1-ac07-978161874e6c
histogram(rand(lengthdist, 10000))

# ╔═╡ 0017844e-2381-4f33-b11f-5711d3c6bda7
md"### 2: Conditional expected value"

# ╔═╡ 6ebf3a16-0e6a-491f-8280-b4327ed52cf0
len_obs = [94.0, 88.7, 89.6, 69.8, 52.8, 84.0, 89.3, 66.4, 95.1, 81.6]

# ╔═╡ b15de91a-fc48-4cdc-a35f-6453a9a59982
@model function fishmixture()
	fs1 ~ Uniform(0, 1) # fraction of species 1
	fishlendist = MixtureModel([Normal(90, 15), Normal(60, 10)], [fs1, 1-fs1])
	
	fishlens = zeros(10)
	for i in 1:length(fishlens)
		fishlens[i] ~ fishlendist
	end
end

# ╔═╡ 4f6bfecd-43a3-44c1-a20c-224c01b8469d
fishmodel = fishmixture() | (fishlens = len_obs,)

# ╔═╡ 41c1e279-4d0b-4447-a9ab-015863df8e91
fishchain = sample(fishmodel, NUTS(), 2000)

# ╔═╡ a786c349-34e3-4dda-b139-808259495753
plot(fishchain)

# ╔═╡ 4f5c0761-12ef-4997-9e40-050c30ec84ab
fs1_est = mean(fishchain[:fs1])

# ╔═╡ 952a941a-8703-45a3-aac1-a290a181e8c5
md"### 3🌟: Conditional expected value (spicy)"

# ╔═╡ 8629d049-b9fd-4e9e-9b55-401a3069e956
@model function fishmixture🌟()
	fs1 ~ Uniform(0, 1) # fraction of species 1
	
	fishlens = zeros(10)
	isspecies1 = zeros(10)
	for i in 1:length(fishlens)
		isspecies1[i] ~ Bernoulli(fs1)
		if isspecies1[i] == 1.0
			fishlens[i] ~ Normal(90, 15)
		else
			fishlens[i] ~ Normal(60, 10)
		end
	end
end

# ╔═╡ 75d4d482-2f30-4c21-be04-0b821635346f
fishmodel🌟 = fishmixture🌟() | (fishlens = len_obs,)

# ╔═╡ bdb0d902-3056-43fb-abb0-15f2494fcf9d
fishchain🌟 = sample(fishmodel🌟, PG(20), 2000)

# ╔═╡ 7372c616-05a8-428d-ab04-a7be9f65653d
plot(fishchain🌟)

# ╔═╡ 40cc67c1-8627-4f1f-b135-a9770f916b53
p_fish4_is_species1 = mean(fishchain🌟["isspecies1[4]"])

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
	xC ~ Uniform(-20, 20)
	yC ~ Uniform(-20, 20)
	
	# generate a radius
	R ~ Uniform(0, 50)
	
	# three random points in polar coordinates
	θ1 ~ Flat() 
		# `Uniform(0, 2*pi)` is also possible but can get the sampler stuck
		# at 0 or 2π!
	θ2 ~ Flat()
	θ3 ~ Flat()
	
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

# ╔═╡ 543c40d5-e8a7-492d-a0b8-e7e73e5953e2
circlemodel = circle() | (x1=x1, y1=y1, x2=x2, y2=y2, x3=x3, y3=y3);

# ╔═╡ 22173937-25a1-4ec0-877d-f9669653e43e
circlechain = sample(circlemodel, NUTS(), 2000)

# ╔═╡ 0b59a6bc-a886-49fd-a3f5-9b20e60882b9
plot(circlechain)

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
circle1 = circle() | (x1=x1, y1=y1)

# ╔═╡ b7ba8f91-f440-4166-af9e-11fcb5f1755f
circle2 = circle1 | (x2=x2, y2=y2)

# ╔═╡ 3b5f3623-3d00-4ba4-9182-5bddacc56567
chain1 = sample(circle1, NUTS(), 100);

# ╔═╡ c11452bc-f404-4e51-8f68-292f5b538c88
chain2 = sample(circle2, NUTS(), 100);

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
# ╠═8f624627-5a65-441d-989b-3fd7397e341c
# ╟─299ba93b-0fc0-4bb3-9a2c-a571ce571f1b
# ╟─6a7be039-c7f5-4cc0-8f21-beb8ac20a181
# ╟─a068aad0-59c2-4b45-aa0f-ce73bd184bdc
# ╟─c72f1e5e-e8f6-41c3-bf75-ff5c6f36cd05
# ╟─b8c69b3b-d9e7-4237-9cbf-95590b2d5770
# ╠═5192d03d-b9d1-4f03-aaf1-f287e044e904
# ╠═78c9ce62-e375-48ac-8083-55ee085c61de
# ╟─78eb7779-f182-4419-b5d8-79a2f5c5d6da
# ╠═3decb2ec-210a-4b2f-842d-6fd40dd3f77b
# ╠═76dd814d-0d9b-4f7e-aff8-990da57d052b
# ╠═b5bd5114-0a0c-4aa3-af00-4f3d05edc5e3
# ╠═4e8b9000-cb61-4f01-9ba9-17276ad0335e
# ╟─8777b133-7d7c-4a85-b89c-2f00093e9984
# ╠═9695ee7e-359b-489c-962b-1bf84b052371
# ╠═014538da-b5ef-41c8-b799-2c000b4c9134
# ╠═b02dc714-e2bb-4ae2-acf9-c37a4389f953
# ╠═b79dafee-4f6f-4301-9267-1b9798c125ea
# ╠═fd540765-95e5-4071-81f7-e689b06cad0c
# ╠═5521933a-a42e-4a67-94b9-84eab52ddf07
# ╟─ca3e730c-a940-4c76-93eb-70ae4aa0e008
# ╠═159caa6a-2ebf-44bd-87a5-b4ab8b085354
# ╟─eda00c08-de49-4d2d-acc4-ba6e21ff0b11
# ╠═b2b16c34-e12a-4bea-8098-313d01913bbf
# ╠═1a78d96a-fe57-42a4-9785-003266117ddb
# ╠═2f5c14e9-709a-40ec-a6cb-ddd870cc1a60
# ╠═7f90e18d-9af8-42fa-984b-aaa7c8c458b5
# ╠═caf0db69-a83f-433a-b848-7d7d8c2fa25e
# ╠═3c1a508b-6d2c-4507-ab9f-752ee93709c9
# ╟─951d0913-1a52-4d5b-b5bb-168487e50ab2
# ╟─a2cdbd3a-9a53-4c34-949e-4604ff1874b4
# ╟─d90c08dc-5842-4c84-9dd3-62dc22a6aab8
# ╠═ec479be9-0d8a-4c8d-97c9-6f64d861924c
# ╠═4e92d825-481f-430a-94a4-fbdf31b679eb
# ╠═0ced04f8-322a-42f6-a1b7-1b30faf9024b
# ╠═e48a2d96-6d7b-4861-b315-0ae580479eda
# ╠═ceb2b8fd-aa47-443b-8e5c-74a48eda1174
# ╠═8c1eb826-6023-4852-aa45-789f6d4b7051
# ╟─5098d57a-cf21-40c0-9091-5233439eb9c4
# ╠═9dd6fa35-37d5-4ab4-ac2c-f4eefdabadd2
# ╟─989b755b-2275-4aae-a3e2-3b107a72fc0b
# ╠═6e3d4f8f-1766-49a6-9e4a-610671d8aa63
# ╠═4cd99dce-12df-4e3a-aa91-ff00c61544a8
# ╠═48ea9893-16f7-457e-93b0-9f168a2f72f4
# ╠═c3d36685-e1dd-4734-bee9-f7aeebed3473
# ╟─ed7c5547-bd82-4ca7-bf51-dd1c201d88af
# ╟─ed124916-331d-43b1-be44-e78c250d79c5
# ╟─dac55632-6e17-479e-8090-5cf8eaa67dad
# ╟─7df97000-5cf0-4a29-b3c3-f867403f4318
# ╠═0fe475dd-411e-474e-9524-ef9fcceed7af
# ╠═382db9a6-edcb-497b-9535-b576cd6badb0
# ╠═fce3e814-97a7-4e12-a7cc-b5b086dc8340
# ╟─ca403113-1c49-4819-9a54-5d1a20fb03cb
# ╟─9c744027-38da-4393-bc08-1feb0a9d8c94
# ╟─490cbc67-3e0c-4c79-a83d-52cf493854a2
# ╟─cd37e5fc-ae9a-429d-844d-4b450b187b5e
# ╠═126d3954-c77d-4c98-abe0-fd87d14e6265
# ╠═13ed391e-a6d4-4c43-9b69-f60154bb129a
# ╠═15fe25d9-5314-4056-8021-cf259ba27c94
# ╠═2a142576-aec2-4311-9c0e-cb68665d59f6
# ╠═0130c903-e3dd-415e-aa17-1b705f9e4ccc
# ╠═11676f8a-d4fe-4999-9ea2-486fca461f65
# ╟─bd7c310c-888e-4380-b261-da2525ce9c11
# ╟─4c7ead00-134b-41d0-89f5-49c7ca14e188
# ╟─1f35d962-a249-4be7-9a96-17eb83fca7d8
# ╠═8fc58fa4-b005-4f32-9eae-a8143582a1ae
# ╠═fe2958a7-e9dd-4eca-979d-a80df12f8735
# ╠═8abafb6a-dc83-422c-82d1-a721a0e1eca0
# ╠═5ba2886c-b2e1-49f4-90c6-549acc808f77
# ╠═ed022247-5959-481f-a81e-41e5ee5a1448
# ╟─9dc0456b-7fd2-4120-8f9e-3de1984ff516
# ╟─225cd579-1e0d-4680-8d3f-5a737a656eb8
# ╟─3ffac5ca-3635-4aa1-bab7-7c28e7a801cb
# ╟─de3151b7-e599-4f9d-9ee1-969830bac0e9
# ╟─23056f1e-128e-463e-a80f-56299397022e
# ╠═3cd0c888-9f11-47f1-a293-b96aa80ea3b0
# ╠═afcf1260-db4d-4ff1-ac07-978161874e6c
# ╟─0017844e-2381-4f33-b11f-5711d3c6bda7
# ╠═6ebf3a16-0e6a-491f-8280-b4327ed52cf0
# ╠═b15de91a-fc48-4cdc-a35f-6453a9a59982
# ╠═4f6bfecd-43a3-44c1-a20c-224c01b8469d
# ╠═41c1e279-4d0b-4447-a9ab-015863df8e91
# ╠═a786c349-34e3-4dda-b139-808259495753
# ╠═4f5c0761-12ef-4997-9e40-050c30ec84ab
# ╟─952a941a-8703-45a3-aac1-a290a181e8c5
# ╠═8629d049-b9fd-4e9e-9b55-401a3069e956
# ╠═75d4d482-2f30-4c21-be04-0b821635346f
# ╠═bdb0d902-3056-43fb-abb0-15f2494fcf9d
# ╠═7372c616-05a8-428d-ab04-a7be9f65653d
# ╠═40cc67c1-8627-4f1f-b135-a9770f916b53
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
