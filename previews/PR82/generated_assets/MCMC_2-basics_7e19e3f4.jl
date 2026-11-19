### A Pluto.jl notebook ###
# v0.20.4

#> [frontmatter]
#> order = "21"
#> title = "4. MCMC basics"
#> date = "2025-08-06"
#> tags = ["exercises"]
#> description = "MCMC basics"
#> layout = "layout.jlhtml"
#> 
#>     [[frontmatter.author]]
#>     name = "Gauthier Vanhaelewyn"

using Markdown
using InteractiveUtils

# ╔═╡ 94c6f31d-1a43-4221-b60c-1fa0ef8738b8
using Pkg; Pkg.activate("../../pluto-deployment-environment")

# ╔═╡ 45bc5b66-c81b-4afb-8a7e-51aff9609c62
using Turing, StatsPlots

# ╔═╡ 41dc8060-cf5e-11ef-26f9-892577e77af0
md"# Inference notebook #2: Basics"

# ╔═╡ 299ba93b-0fc0-4bb3-9a2c-a571ce571f1b
md"## 1: Mole burrow"

# ╔═╡ 957386c5-775c-47f7-9a38-e1630e548689
md"""
Consider a mole's underground tunnel network of length `X` (in m). Now and then the mole makes a new molehill somewhere randomly above its tunnel, the locations of which we denote `Y`.

We can formulate this as `X ~ Exponential(100)` and `Y ~ Uniform(0, X)`.

!!! questions
	1. Plot the prior of `X`. Is it diffuse or informative?
	1. Estimate `E[Y]`.
	1. Estimate `E[X|Y = 3]` and compare it with the prior expected value `E[X]`.
	1. Plot the histogram of `X` given `Y = 3.0`.
	1. Plot the histogram of `X` again, but now given the following values for `Y`: `[3.0, 1.5, 0.9, 5.7]`.
"""

# ╔═╡ ca8057cc-6d70-4a88-88aa-c494af9659ef
md"""
**Questions:**
- How should I understand this?
- If the mole digs a tunnels that is X long, then at position Y there is one hole?
- Or, is there a hole at every position Y for an X meters long tunnel?
"""

# ╔═╡ 47aa2304-d312-40e9-b9c6-7c79a7d64de4
md"### 1"

# ╔═╡ 78c9ce62-e375-48ac-8083-55ee085c61de
plot(Exponential(100)) # The prior only incorporates the knowledge that a mole's tunnel is probably less long than a few km - this is a diffuse to weakly informative prior

# ╔═╡ 4ab8aaba-010c-4844-9cf6-627a746c2492
mijn_X = rand(Exponential(100))

# ╔═╡ 82289137-9c20-4d7a-97b2-56f417831409
mijn_Y = rand(Uniform(0, mijn_X))

# ╔═╡ 78eb7779-f182-4419-b5d8-79a2f5c5d6da
md"### 2"

# ╔═╡ 3decb2ec-210a-4b2f-842d-6fd40dd3f77b
@model function mole()
	X ~ Exponential(100)
	Y ~ Uniform(0, X)
	return Y
end

# ╔═╡ 76dd814d-0d9b-4f7e-aff8-990da57d052b
molemodel = mole()

# ╔═╡ 5c989060-40cd-4b27-b96f-3f98a85f5122
molemodel()   # returns Y

# ╔═╡ 4e8b9000-cb61-4f01-9ba9-17276ad0335e
E_Y = mean([molemodel() for i in 1:2000])

# ╔═╡ 8777b133-7d7c-4a85-b89c-2f00093e9984
md"### 3"

# ╔═╡ 9695ee7e-359b-489c-962b-1bf84b052371
cond_mole = molemodel | (Y = 3.0,);

# ╔═╡ 10deca95-9f88-45ff-a22f-30280776acef
md"""
**Questions:**
- Why doesn't this work for much higher values of Y? E.g. for Y = 8 or 10?
- Because if you sample `X = rand(Exponential(100))` and `Y = rand(Uniform(0, X))`, you usually get much larger numbers for Y than 3.
"""

# ╔═╡ 014538da-b5ef-41c8-b799-2c000b4c9134
molechain = sample(cond_mole, NUTS(), 2000)

# ╔═╡ b02dc714-e2bb-4ae2-acf9-c37a4389f953
plot(molechain)

# ╔═╡ fd540765-95e5-4071-81f7-e689b06cad0c
E_Xcond3 = mean(molechain[:X])

# ╔═╡ 5521933a-a42e-4a67-94b9-84eab52ddf07
E_X = mean(Exponential(100))

# ╔═╡ ca3e730c-a940-4c76-93eb-70ae4aa0e008
md"### 4"

# ╔═╡ 159caa6a-2ebf-44bd-87a5-b4ab8b085354
histogram(molechain[:X], normalized=:probability)

# ╔═╡ eda00c08-de49-4d2d-acc4-ba6e21ff0b11
md"### 5"

# ╔═╡ b2b16c34-e12a-4bea-8098-313d01913bbf
@model function mole2()
	X ~ Exponential(100)
	Ys = zeros(4)            # now we have 4 holes Ys (Ysample)
	for i in eachindex(Ys)
		Ys[i] ~ Uniform(0, X)
	end
end

# ╔═╡ 1a78d96a-fe57-42a4-9785-003266117ddb
Y_obs = [3.0, 1.5, 0.9, 5.7]

# ╔═╡ 2f5c14e9-709a-40ec-a6cb-ddd870cc1a60
mole_cond2 = mole2() | (Ys = Y_obs,)

# ╔═╡ 7f90e18d-9af8-42fa-984b-aaa7c8c458b5
molechain2 = sample(mole_cond2, NUTS(), 2000)

# ╔═╡ caf0db69-a83f-433a-b848-7d7d8c2fa25e
plot(molechain2)

# ╔═╡ 3c1a508b-6d2c-4507-ab9f-752ee93709c9
histogram(molechain2[:X], normalized=:probability)

# ╔═╡ 2d35375a-cd92-48d1-8779-761641e7b0af
mean(molechain2[:X])

# ╔═╡ 951d0913-1a52-4d5b-b5bb-168487e50ab2
md"## 2: Potatoes"

# ╔═╡ 49035a16-c419-4531-b157-a5ab357b44fe
md"""
Consider a number of potatoes `N` each with an average weight `W`. You weigh them together on an old balance to get an estimate of their total weight `T`.

We can formulate this as `N ~ Poisson(10)`, `W ~ Uniform(150, 250)`  and `T ~ Normal(N*W, 50)`.

!!! questions
	1. Plot a histogram of `N` given `T = 1200`.
	1. Estimate `P(N > 6, W > 175 | T = 1200)`.
	1. Estimate `P(N = 5 | T = 1200, W = 220)`.
"""

# ╔═╡ d5f544a4-ba51-4103-82ce-259ba20af11a
let
	N = rand(Poisson(10));
	W = rand(Uniform(150, 250));   # in grams
	T = rand(Normal(N*W, 50))      # in grams
end

# ╔═╡ 5f8322ca-0ded-4750-8ad6-e66e8280daca
md"### 1"

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
potato_chain = sample(potato_cond, PG(10), 2000)

# ╔═╡ ceb2b8fd-aa47-443b-8e5c-74a48eda1174
plot(potato_chain)

# ╔═╡ 8c1eb826-6023-4852-aa45-789f6d4b7051
histogram(potato_chain[:N], normalized=:probability)

# ╔═╡ 205efb15-e28a-418a-a1af-892afdf3c188
md"""
The above gives of the individual probabilities of having T=1200 with 5, 6, 7 and 8 potatoes.
"""

# ╔═╡ ad571d35-639c-4239-9e62-42952ddb48d0
# Probability of having T=1200 with N=7.
mean(potato_chain[:N] .== 7)

# ╔═╡ 0fddec58-67d9-4dec-a72d-445675fa46a6
md"### 2"

# ╔═╡ 9dd6fa35-37d5-4ab4-ac2c-f4eefdabadd2
p_potato1 = mean(potato_chain[:N] .> 6 .&& potato_chain[:W] .> 175)

# ╔═╡ 57ec05ca-bf37-4012-be9a-d83c1d9f0f8a
mean(potato_chain[:N] .<= 6 .|| potato_chain[:W] .<= 175)

# ╔═╡ 5bde987a-3f5c-4ff7-96a0-9e174f73fdfb
md"### 3"

# ╔═╡ 6e3d4f8f-1766-49a6-9e4a-610671d8aa63
potato_cond2 = potato_model | (T = 1200, W = 220,)
# potato_cond2 = potato_model | (T = 1100, W = 220,)

# ╔═╡ 4cd99dce-12df-4e3a-aa91-ff00c61544a8
pota2_chain = sample(potato_cond2, PG(10), 2000)

# ╔═╡ 48ea9893-16f7-457e-93b0-9f168a2f72f4
plot(pota2_chain)

# ╔═╡ c3d36685-e1dd-4734-bee9-f7aeebed3473
p_potato2 = mean(pota2_chain[:N] .== 5)

# ╔═╡ bdc40079-de7b-489f-9d9b-6e1e2130109f
mean(pota2_chain[:N] .== 4)

# ╔═╡ 22d009f0-390f-412d-9af8-33d169ac7371
mean(pota2_chain[:N] .== 6)

# ╔═╡ ed7c5547-bd82-4ca7-bf51-dd1c201d88af
md"## 3: Lights out"

# ╔═╡ a3c26839-2acc-49d8-98c4-3a7142dd6512
md"""
You use 4 of the same LED light in your room. Let `μ` be the **average** lifespan of your LED lights (in khr or 1000 hours) and `L`ᵢ the lifespan of the `i`-th LED light.

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
md"### 1"

# ╔═╡ 0fe475dd-411e-474e-9524-ef9fcceed7af
lights_prior = LogNormal(log(40), 0.5)

# ╔═╡ 382db9a6-edcb-497b-9535-b576cd6badb0
plot(lights_prior) # not asked but a visualisation can always be useful

# ╔═╡ fce3e814-97a7-4e12-a7cc-b5b086dc8340
E_mu = mean(lights_prior)

# ╔═╡ a7811e57-2820-4662-8c20-a4d2d7efd642
rand(lights_prior)

# ╔═╡ 308cb1f5-e26b-43b2-be4d-9eabd68b3670
md"### 2"

# ╔═╡ 490cbc67-3e0c-4c79-a83d-52cf493854a2
md"""
The exponential distribution is often used to model the waiting time for an event. This makes it a natural fit for a lamp's lifespan, which is the waiting time until it breaks. We know it needs to have a mean value of μ, so `Exponential(μ)` is a good choice.

One could also argue for a LogNormal distribution with mean μ or a Normal distribution with mean μ restricted to only the positive values. Both would need a large variance to reflect the lack of additional information outside of the mean lifespan.
"""

# ╔═╡ dbadf227-1399-45f7-ad48-df9d59e18343
# Average waiting time before it breaks is 45 khr. The longer we wait, the
# higher the chance that it breaks!
plot(Exponential(45))

# ╔═╡ 5d812594-733e-4b18-88dc-c28867cf0933
μ = 27

# ╔═╡ 0373b6a5-bb2c-47b9-a137-5a484b0c7bf2
plot(Normal(μ, sqrt(μ)))

# ╔═╡ ec06259f-67b5-4920-8635-35c52c080751
begin
	plot([LogNormal(log(μ), 1), Exponential(45)]; xlim=[0,400])
	vline!([mean(LogNormal(log(μ), 1)), mean(Exponential(45))])
end

# ╔═╡ a628bd3f-5d07-435e-80cc-b72662e13e7b
plot(cdf(Exponential(45), 0:400))

# ╔═╡ e2fd7500-277a-4d52-982c-5aea8421d756
cdf(Exponential(45), 100)   # probability that it breaks before 100 khr

# ╔═╡ 96ee6380-719e-47c1-a429-4fdd1931c9ce
1 - cdf(Exponential(45), 100)  # probability that it is still working after 100 khr

# ╔═╡ cd37e5fc-ae9a-429d-844d-4b450b187b5e
md"### 3"

# ╔═╡ 126d3954-c77d-4c98-abe0-fd87d14e6265
@model function lights()
	μ ~ lights_prior
	lifespans = zeros(4)
	for i in 1:length(lifespans)
		lifespans[i] ~ LogNormal(log(μ), 1) # Normal(μ, sqrt(μ)) #Exponential(μ) 
	end
end

# ╔═╡ 15fe25d9-5314-4056-8021-cf259ba27c94
lightmodel = lights() | (lifespans = [16, 20, 23, 41],)

# ╔═╡ 2a142576-aec2-4311-9c0e-cb68665d59f6
lightschain = sample(lightmodel, NUTS(), 2000)

# ╔═╡ 0130c903-e3dd-415e-aa17-1b705f9e4ccc
plot(lightschain)

# ╔═╡ 11676f8a-d4fe-4999-9ea2-486fca461f65
E_mu_cond = mean(lightschain[:μ])

# ╔═╡ f6bc2fad-24c7-48a7-9328-34b616fa106f
histogram(lightschain[:μ], normalize=true)

# ╔═╡ 347a4533-244f-4e93-b125-cdd56a5a08ad
begin
	plot([LogNormal(log(μ), 1), Exponential(45)]; xlim=[0,400])
	vline!([mean(LogNormal(log(μ), 1)), mean(Exponential(45))])
end

# ╔═╡ 9a168680-9d96-466f-ba75-122d6a391501
md"### 4 🌟🌟🌟"

# ╔═╡ 1f35d962-a249-4be7-9a96-17eb83fca7d8
md"""
!!! hint
	You can model the number of lights that still work as a `Binomial` distribution, the success rate of which depends on `μ`.
"""

# ╔═╡ 3ad29b60-375b-400b-b1e5-8918da6497ff
1-cdf(Exponential(40), 30)

# ╔═╡ 0cce9dee-43ac-42d3-9113-8b1503c7a73c
rand(Binomial(2 + 2, 0.47))

# ╔═╡ 8d501c17-58d2-43d3-92c2-0d21bc9113e3
md"""
After 30 khr, two lights have died: one at 16 khr and one at 20 khr. The two other lights are still working. What is the expected value of `μ` given this information?
"""

# ╔═╡ c626d3bf-4abe-4bdc-991a-363bee27b25b
md"""
Here you need to provide two arguments to the model function:
- How many lights still working? `n`
- At what time they are still working? `time_observed`

"""

# ╔═╡ 8fc58fa4-b005-4f32-9eae-a8143582a1ae
@model function lights_censored(n, time_observed)
	μ ~ lights_prior
	
	lifespans = zeros(2)
	for i in 1:length(lifespans)
		lifespans[i] ~ Exponential(μ)	
	end

	# Given the observation time, what is the probability that
	# a single light still works:
	p_stillworking = 1 - cdf(Exponential(μ), time_observed)
	# cdf(Exponential(μ), time_observed) is the  probability that it broke
	#     in [0, time_observed]
	# Number of lights still working with the above probability:
	n ~ Binomial(n + length(lifespans), p_stillworking)

	return (μ, lifespans, p_stillworking, n)
end

# ╔═╡ fe2958a7-e9dd-4eca-979d-a80df12f8735
lightmodel_cens = lights_censored(2, 30) | (lifespans = [16, 20],)

# ╔═╡ 93f6b680-aece-45d7-804e-a27073130d94
lightmodel_cens()

# ╔═╡ 8abafb6a-dc83-422c-82d1-a721a0e1eca0
lightschain_cens = sample(lightmodel_cens, NUTS(), 2000)

# ╔═╡ 5ba2886c-b2e1-49f4-90c6-549acc808f77
plot(lightschain_cens)

# ╔═╡ ed022247-5959-481f-a81e-41e5ee5a1448
E_mu_cond🌟 = mean(lightschain_cens[:μ])

# ╔═╡ c57046cd-5f53-43a2-9de4-9bfad27dcab9
@model function lights_censored2(time_observed)
	μ ~ lights_prior
	
	lifespans = zeros(2)
	for i in 1:length(lifespans)
		lifespans[i] ~ Exponential(μ)	
	end
	
	p_stillworking = 1 - cdf(Exponential(μ), time_observed)
	n ~ Binomial(4, p_stillworking)
	
	return (μ, lifespans, p_stillworking, n)
end

# ╔═╡ a2ebe90c-1006-4fc3-b8bd-bdf9c9e9fed5
plot(cdf(Exponential(10), 0:100))

# ╔═╡ 5b505f97-e776-4cc3-99f8-721dc241ae2a
lightmodel_cens2 = lights_censored2(30) | (lifespans = [16, 20], n = 2,)

# ╔═╡ deaf3485-9bc1-48da-809c-61327cce5ad9
lightmodel_cens2()

# ╔═╡ 56a91508-29fa-4e5f-8dcc-0bb9d67db1fe
lightschain_cens2 = sample(lightmodel_cens2, NUTS(), 2000)

# ╔═╡ 878e5eaf-2f1c-4c09-a9dd-47821e9295e5
plot(lightschain_cens2)

# ╔═╡ 7366d15b-9c03-4200-8720-6d879bb65476
mean(lightschain_cens2[:μ])

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
md"### 1"

# ╔═╡ 23056f1e-128e-463e-a80f-56299397022e
md"""
!!! hint
	The distribution of fish lengths can be modelled as a `MixtureModel`.
"""

# ╔═╡ 3cd0c888-9f11-47f1-a293-b96aa80ea3b0
lengthdist = MixtureModel([Normal(90, 15), Normal(60, 10)], [0.3, 0.7])

# ╔═╡ afcf1260-db4d-4ff1-ac07-978161874e6c
histogram(rand(lengthdist, 10000))

# ╔═╡ 8e030a06-5104-4c5b-b1f2-f86464e66502
md"### 2"

# ╔═╡ 6ebf3a16-0e6a-491f-8280-b4327ed52cf0
len_obs = [94.0, 88.7, 89.6, 69.8, 52.8, 84.0, 89.3, 66.4, 95.1, 81.6]

# ╔═╡ b15de91a-fc48-4cdc-a35f-6453a9a59982
@model function fishmixture()
	fs1 ~ Uniform(0, 1) # fraction of species 1
	# fish length distribution
	fishlendist = MixtureModel([Normal(90, 15), Normal(60, 10)], [fs1, 1-fs1])
	
	fishlens = zeros(10)    # fish lengths
	for i in eachindex(fishlens)
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
md"### 3🌟"

# ╔═╡ 8629d049-b9fd-4e9e-9b55-401a3069e956
@model function fishmixture🌟()
	fs1 ~ Uniform(0, 1)      # fraction of species 1 
	                         # or probability of belonging to species 1
	
	fishlens = zeros(10)     # samples with fish lengths
	isspecies1 = zeros(10)   # samples with 1's meaning belonging to species 1
	                         # samples with 0's meaning belonging to species 2
	for i in eachindex(fishlens)
		isspecies1[i] ~ Bernoulli(fs1) # samples belonging to species 1 or not
		if isspecies1[i] == 1.0        # if belongs to species 1
			fishlens[i] ~ Normal(90, 15)  # sample from distribution of species 1
		else
			fishlens[i] ~ Normal(60, 10)  # sample from distribution of species 1
		end
	end
end

# ╔═╡ 75d4d482-2f30-4c21-be04-0b821635346f
fishmodel🌟 = fishmixture🌟() | (fishlens = len_obs,)

# ╔═╡ bdb0d902-3056-43fb-abb0-15f2494fcf9d
fishchain🌟 = sample(fishmodel🌟, PG(20), 2000)

# ╔═╡ 7372c616-05a8-428d-ab04-a7be9f65653d
plot(fishchain🌟)

# ╔═╡ a60315f2-92f5-4f3b-b154-c3428816bfb5
# fish species 1 -> large fish
# fish species 2 -> small fish
#        [94.0, 88.7, 89.6, 69.8, 52.8, 84.0, 89.3, 66.4, 95.1, 81.6]
# fish:    1     2     3     4     5     6     7     8     9    10
# species: 1     1     1    1or2   2     1     1    1or2   1     1

# ╔═╡ 05c700e8-f24a-4c82-9d15-cb3450de9e2a
# We expect a very high chance here because 94.0 is a large fish
p_fish1_is_species1 = mean(fishchain🌟["isspecies1[1]"])

# ╔═╡ 40cc67c1-8627-4f1f-b135-a9770f916b53
# We expect medium chance here because 69.8 is between small and large
p_fish4_is_species1 = mean(fishchain🌟["isspecies1[4]"])

# ╔═╡ 8b9d7708-4a67-457a-9fe1-de12fb2a2de9
# We expect a very small chance here because 52.8 is a small fish
p_fish5_is_species1 = mean(fishchain🌟["isspecies1[5]"])

# ╔═╡ 5d8929fd-5aa4-4186-9fd0-e0a64ac06cbe
# We expect medium chance here because 66.4 is between small and large
p_fish8_is_species1 = mean(fishchain🌟["isspecies1[8]"])

# ╔═╡ 8e09d7f1-f8db-42a1-beb0-ed3afb043c19
mean(fishchain🌟[:fs1])   # should be the same as before

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
scatter([x1, x2, x3], [y1, y2, y3], aspect_ratio=:equal, xlim=[-40, 40], ylim=[-40, 40])

# ╔═╡ 7eedb74d-eee1-4cf0-b2bf-5febf474edd2
md"### 1"

# ╔═╡ abbb225e-fdcd-435e-9f81-c1e6fa1c8f5d
# Flat() -> p134 in syllabus
rand(Flat())

# ╔═╡ f1086281-ba2d-461d-bcfe-82e3d09aaae0
rand(2*π*Flat())

# ╔═╡ c4a2335d-5926-4da4-a46f-41da7b67fa01
rand(Uniform(0, 2*π))

# ╔═╡ 84d27c98-9513-4ae3-8101-621c083a1b01
@model function circle(σ=0.25)
	# generate a circle center
	xC ~ Uniform(-20, 20)
	yC ~ Uniform(-20, 20)
	
	# generate a radius
	R ~ Uniform(0, 50)
	
	# three random points in polar coordinates
	θ1 ~ 2*π*Flat() 
		# `Uniform(0, 2*pi)` is also possible but can get the sampler stuck
		# at 0 or 2π!
	θ2 ~ 2*π*Flat()
	θ3 ~ 2*π*Flat()
	
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

# ╔═╡ 121c291f-ec34-47cb-a9fc-01309b0f42c5
# Circle center at:
(mean(circlechain[:xC]), mean(circlechain[:yC]))

# ╔═╡ 7efa1868-08b9-49fb-b799-76f1c99d1a0e
mean.([circlechain[:xC], circlechain[:yC]])

# ╔═╡ 374d7cea-2dc6-4cb3-b929-ddaba1f4c6fb
# Radius:
mean(circlechain[:R])

# ╔═╡ 95d223a0-748b-4656-8777-e151c0dc8e93
mean(circlechain[:θ1])*180.0/π

# ╔═╡ 2668e03a-91eb-420f-af3f-b74d83126f76
mean(circlechain[:θ2])*180.0/π

# ╔═╡ 87be71c4-2bc3-432b-a977-f3d46424a68f
mean(circlechain[:θ3])*180.0/π

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
md"### 2"

# ╔═╡ b7649523-fd40-4fe3-8d86-fc2cb2c8c488
circle1 = circle() | (x1=x1, y1=y1)
# This can be any circle through (x1, y1) with centrer
# in [-20, 20]x[-20, 20] and radius between 0 and 50.

# ╔═╡ b7ba8f91-f440-4166-af9e-11fcb5f1755f
circle2 = circle1 | (x2=x2, y2=y2)
# This can be any circle through (x1, y1) and (x2, y2) with centrer
# in [-20, 20]x[-20, 20] and radius between 0 and 50.

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
# ╟─299ba93b-0fc0-4bb3-9a2c-a571ce571f1b
# ╟─957386c5-775c-47f7-9a38-e1630e548689
# ╟─ca8057cc-6d70-4a88-88aa-c494af9659ef
# ╟─47aa2304-d312-40e9-b9c6-7c79a7d64de4
# ╠═78c9ce62-e375-48ac-8083-55ee085c61de
# ╠═4ab8aaba-010c-4844-9cf6-627a746c2492
# ╠═82289137-9c20-4d7a-97b2-56f417831409
# ╟─78eb7779-f182-4419-b5d8-79a2f5c5d6da
# ╠═3decb2ec-210a-4b2f-842d-6fd40dd3f77b
# ╠═76dd814d-0d9b-4f7e-aff8-990da57d052b
# ╠═5c989060-40cd-4b27-b96f-3f98a85f5122
# ╠═4e8b9000-cb61-4f01-9ba9-17276ad0335e
# ╟─8777b133-7d7c-4a85-b89c-2f00093e9984
# ╠═9695ee7e-359b-489c-962b-1bf84b052371
# ╟─10deca95-9f88-45ff-a22f-30280776acef
# ╠═014538da-b5ef-41c8-b799-2c000b4c9134
# ╠═b02dc714-e2bb-4ae2-acf9-c37a4389f953
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
# ╠═2d35375a-cd92-48d1-8779-761641e7b0af
# ╟─951d0913-1a52-4d5b-b5bb-168487e50ab2
# ╟─49035a16-c419-4531-b157-a5ab357b44fe
# ╠═d5f544a4-ba51-4103-82ce-259ba20af11a
# ╟─5f8322ca-0ded-4750-8ad6-e66e8280daca
# ╠═ec479be9-0d8a-4c8d-97c9-6f64d861924c
# ╠═4e92d825-481f-430a-94a4-fbdf31b679eb
# ╠═0ced04f8-322a-42f6-a1b7-1b30faf9024b
# ╠═e48a2d96-6d7b-4861-b315-0ae580479eda
# ╠═ceb2b8fd-aa47-443b-8e5c-74a48eda1174
# ╠═8c1eb826-6023-4852-aa45-789f6d4b7051
# ╟─205efb15-e28a-418a-a1af-892afdf3c188
# ╠═ad571d35-639c-4239-9e62-42952ddb48d0
# ╟─0fddec58-67d9-4dec-a72d-445675fa46a6
# ╠═9dd6fa35-37d5-4ab4-ac2c-f4eefdabadd2
# ╠═57ec05ca-bf37-4012-be9a-d83c1d9f0f8a
# ╟─5bde987a-3f5c-4ff7-96a0-9e174f73fdfb
# ╠═6e3d4f8f-1766-49a6-9e4a-610671d8aa63
# ╠═4cd99dce-12df-4e3a-aa91-ff00c61544a8
# ╠═48ea9893-16f7-457e-93b0-9f168a2f72f4
# ╠═c3d36685-e1dd-4734-bee9-f7aeebed3473
# ╠═bdc40079-de7b-489f-9d9b-6e1e2130109f
# ╠═22d009f0-390f-412d-9af8-33d169ac7371
# ╟─ed7c5547-bd82-4ca7-bf51-dd1c201d88af
# ╟─a3c26839-2acc-49d8-98c4-3a7142dd6512
# ╟─dac55632-6e17-479e-8090-5cf8eaa67dad
# ╟─7df97000-5cf0-4a29-b3c3-f867403f4318
# ╠═0fe475dd-411e-474e-9524-ef9fcceed7af
# ╠═382db9a6-edcb-497b-9535-b576cd6badb0
# ╠═fce3e814-97a7-4e12-a7cc-b5b086dc8340
# ╠═a7811e57-2820-4662-8c20-a4d2d7efd642
# ╟─308cb1f5-e26b-43b2-be4d-9eabd68b3670
# ╟─490cbc67-3e0c-4c79-a83d-52cf493854a2
# ╠═dbadf227-1399-45f7-ad48-df9d59e18343
# ╠═5d812594-733e-4b18-88dc-c28867cf0933
# ╠═0373b6a5-bb2c-47b9-a137-5a484b0c7bf2
# ╠═ec06259f-67b5-4920-8635-35c52c080751
# ╠═a628bd3f-5d07-435e-80cc-b72662e13e7b
# ╠═e2fd7500-277a-4d52-982c-5aea8421d756
# ╠═96ee6380-719e-47c1-a429-4fdd1931c9ce
# ╟─cd37e5fc-ae9a-429d-844d-4b450b187b5e
# ╠═126d3954-c77d-4c98-abe0-fd87d14e6265
# ╠═15fe25d9-5314-4056-8021-cf259ba27c94
# ╠═2a142576-aec2-4311-9c0e-cb68665d59f6
# ╠═0130c903-e3dd-415e-aa17-1b705f9e4ccc
# ╠═11676f8a-d4fe-4999-9ea2-486fca461f65
# ╠═f6bc2fad-24c7-48a7-9328-34b616fa106f
# ╠═347a4533-244f-4e93-b125-cdd56a5a08ad
# ╟─9a168680-9d96-466f-ba75-122d6a391501
# ╟─1f35d962-a249-4be7-9a96-17eb83fca7d8
# ╠═3ad29b60-375b-400b-b1e5-8918da6497ff
# ╠═0cce9dee-43ac-42d3-9113-8b1503c7a73c
# ╟─8d501c17-58d2-43d3-92c2-0d21bc9113e3
# ╟─c626d3bf-4abe-4bdc-991a-363bee27b25b
# ╠═8fc58fa4-b005-4f32-9eae-a8143582a1ae
# ╠═fe2958a7-e9dd-4eca-979d-a80df12f8735
# ╠═93f6b680-aece-45d7-804e-a27073130d94
# ╠═8abafb6a-dc83-422c-82d1-a721a0e1eca0
# ╠═5ba2886c-b2e1-49f4-90c6-549acc808f77
# ╠═ed022247-5959-481f-a81e-41e5ee5a1448
# ╠═c57046cd-5f53-43a2-9de4-9bfad27dcab9
# ╠═a2ebe90c-1006-4fc3-b8bd-bdf9c9e9fed5
# ╠═5b505f97-e776-4cc3-99f8-721dc241ae2a
# ╠═deaf3485-9bc1-48da-809c-61327cce5ad9
# ╠═56a91508-29fa-4e5f-8dcc-0bb9d67db1fe
# ╠═878e5eaf-2f1c-4c09-a9dd-47821e9295e5
# ╠═7366d15b-9c03-4200-8720-6d879bb65476
# ╟─9dc0456b-7fd2-4120-8f9e-3de1984ff516
# ╟─225cd579-1e0d-4680-8d3f-5a737a656eb8
# ╟─3ffac5ca-3635-4aa1-bab7-7c28e7a801cb
# ╟─b477b212-83a2-42f0-a616-52516e152d48
# ╟─23056f1e-128e-463e-a80f-56299397022e
# ╠═3cd0c888-9f11-47f1-a293-b96aa80ea3b0
# ╠═afcf1260-db4d-4ff1-ac07-978161874e6c
# ╟─8e030a06-5104-4c5b-b1f2-f86464e66502
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
# ╠═a60315f2-92f5-4f3b-b154-c3428816bfb5
# ╠═05c700e8-f24a-4c82-9d15-cb3450de9e2a
# ╠═40cc67c1-8627-4f1f-b135-a9770f916b53
# ╠═8b9d7708-4a67-457a-9fe1-de12fb2a2de9
# ╠═5d8929fd-5aa4-4186-9fd0-e0a64ac06cbe
# ╠═8e09d7f1-f8db-42a1-beb0-ed3afb043c19
# ╟─30088664-5157-4d99-8584-7a42d0acdfb8
# ╟─7dd2c189-79c0-4d29-9e17-9c24a78b5791
# ╟─dcbf405c-786c-4226-b35c-dc718452bb61
# ╟─28d28034-e999-4e34-b6b2-63c762094c59
# ╠═791ff4ed-c9d1-48e2-9dd0-4bf3979c6167
# ╠═cb189957-f9d4-480b-a492-92cfc2a8c2aa
# ╠═81eca3b3-12d9-43d7-af14-e9aeb73f2471
# ╠═4e96908a-4fc9-429d-bf37-7a569194a038
# ╟─7eedb74d-eee1-4cf0-b2bf-5febf474edd2
# ╠═abbb225e-fdcd-435e-9f81-c1e6fa1c8f5d
# ╠═f1086281-ba2d-461d-bcfe-82e3d09aaae0
# ╠═c4a2335d-5926-4da4-a46f-41da7b67fa01
# ╠═84d27c98-9513-4ae3-8101-621c083a1b01
# ╠═543c40d5-e8a7-492d-a0b8-e7e73e5953e2
# ╠═22173937-25a1-4ec0-877d-f9669653e43e
# ╠═0b59a6bc-a886-49fd-a3f5-9b20e60882b9
# ╠═121c291f-ec34-47cb-a9fc-01309b0f42c5
# ╠═7efa1868-08b9-49fb-b799-76f1c99d1a0e
# ╠═374d7cea-2dc6-4cb3-b929-ddaba1f4c6fb
# ╠═95d223a0-748b-4656-8777-e151c0dc8e93
# ╠═2668e03a-91eb-420f-af3f-b74d83126f76
# ╠═87be71c4-2bc3-432b-a977-f3d46424a68f
# ╠═0d04b0e8-3d1a-4281-b175-570148569ef2
# ╟─0a140d2a-a24b-48df-9af8-7fa5d586a26f
# ╠═b7649523-fd40-4fe3-8d86-fc2cb2c8c488
# ╠═b7ba8f91-f440-4166-af9e-11fcb5f1755f
# ╠═3b5f3623-3d00-4ba4-9182-5bddacc56567
# ╠═c11452bc-f404-4e51-8f68-292f5b538c88
# ╠═c1174a5c-a6a1-46a0-96be-6997e3201dfc
# ╠═793fc102-fc64-4041-a04e-e0b1b0741437
