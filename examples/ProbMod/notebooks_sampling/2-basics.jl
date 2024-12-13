### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 7d4d4d20-b323-11ef-0926-b14785cb9ab5
using Pkg; Pkg.activate("..")

# ╔═╡ 4cfd4721-e29a-4270-8d15-021bcc966eb1
using Turing, StatsPlots

# ╔═╡ 0fc3fa95-0287-423c-ba5a-7e382202ff81
n_samples = 1_000

# ╔═╡ 9740ea64-cd4f-46b1-a741-02e392280601
md"""
# 1: Dice
"""

# ╔═╡ 187854bb-9e30-454d-9e03-cccf77aebb6b
md"You're playing a fun game of Caverns and Chimeras, and are facing off against the mighty Carl the Chimera. The fight is not going great and your next spell **needs to deal 50 or more damage** to slay the scary monster before it kills you. Spells deal **damage equal to the sum of the dice** they let you roll.

You can choose between your 2 mightiest spells:
- Watercube: lets you throw **4 dice with 20 sides** each
- Dirtprism: lets you throw **20 dice with 4 sides** each
"

# ╔═╡ 747e3c0a-357a-448a-b479-d0fcbe44a6c0
md"""
!!! questions
    - What are the probabilities of either attack doing the job? 
    - Plot a histogram of the damage of both attacks.
"""

# ╔═╡ eda95f45-083a-4e65-b57e-bd9890da1f9c
md"## 4d20"

# ╔═╡ 36477423-5628-4ed3-b54d-9a050557f6b7
@model function dice()
    roll1 ~ DiscreteUniform(1, 20)
	roll2 ~ DiscreteUniform(1, 20)
	roll3 ~ DiscreteUniform(1, 20)
	roll4 ~ DiscreteUniform(1, 20)
	# or use a loop

    dicesum = roll1 + roll2 + roll3 + roll4
    return dicesum
end

# ╔═╡ 02bc37e2-5cf6-404a-b3fe-b2120671adb2
dicemodel = dice()

# ╔═╡ 97021710-4b26-4a1c-a794-265c9559ce4d
sp_ds = [dicemodel() for i in 1:n_samples]

# ╔═╡ c0cd75bd-2a46-4f62-a59a-9bb8d47d45f2
mean(sp_ds .>= 50)

# ╔═╡ bc5f5d70-e269-45f9-867b-a00876ce8c40
histogram(sp_ds, bins = 30)

# ╔═╡ a1b933ac-5d1b-4800-a6e8-e942846b19d8
md"## 20d4"

# ╔═╡ 6b009ba3-83a8-4176-86d4-dd9f70ed29ec
@model function superdice()
    rolls = zeros(20)
	for i in eachindex(rolls)
		rolls[i] ~ DiscreteUniform(1, 4)
	end
    dicesum = sum(rolls)
    return dicesum
end

# ╔═╡ 83afb3c1-9e6a-4d18-b0c7-05ed0173df40
superdicemodel = superdice()

# ╔═╡ 700cc645-1123-4c8a-821e-1ae4fe8b0674
sp_sds = [superdicemodel() for i in 1:n_samples]

# ╔═╡ 98d04790-e86f-438a-9389-cb9697934bbf
mean(sp_sds .>= 50)

# ╔═╡ d5413573-5202-4fa2-94f6-e92a613d63aa
histogram(sp_sds, bins = 30)

# ╔═╡ 34f3014f-f4d4-43d1-b46f-bdca73aee33f
md"# 2: Super eggs"

# ╔═╡ 372436c4-262f-49b8-b1cf-626b043542bf
md"""
When a chicken lays an egg, there's a small chance it contains two egg yolks. This chance, as well as the number of eggs a chicken lays per year, go down as the chicken gets older. 
"""

# ╔═╡ 20111742-008a-44c3-8c27-62791cce3e1e
md"""
You can make the following assumptions
- The age $A$ of a random chicken (in years) is Poisson distributed with mean 2.
- The amount of eggs an $A$-year old chicken lays in a year $N$ is Poisson distributed with mean $300 - 20 \, A$.
- The probability of an $A$-year old chicken's egg having a double yolk $P$ is Normally distributed with mean $\dfrac{1}{800 + 100\,A}$ and a standard deviation of $10^{-4}$.
"""

# ╔═╡ 6e020801-983d-4ebc-a0e9-b5dd58f66c55
md"""
!!! questions
    - If someone hands you a random chicken, what is the probability it will lay 2 or more double eggs in a year? 
    - Compare the distributions of double eggs for 1-year old and 3-year old chickens.
"""

# ╔═╡ 2b3d930f-53d9-4869-9e4a-86a1a681b9d8
@model function eggs()
	A ~ Poisson(2)
	N ~ Poisson(300 - 20*A)
	P ~ Normal(1/(800+100*A), 1e-4)

	doubles ~ Binomial(N, P)
	return sum(doubles)
end

# ╔═╡ 834139c2-29bf-4882-bd73-1f89ad9f3803
egg_model = eggs();

# ╔═╡ ed67949b-c7b3-46c2-af99-ab64dbef4065
chain_egg = sample(egg_model, Prior(), n_samples)

# ╔═╡ 64d04ed9-6548-4251-8a4e-11b31e8f3143
sp_egg = generated_quantities(egg_model, chain_egg);

# ╔═╡ 00fd5517-5232-4d88-8ab0-3a0eb925eb3e
mean(sp_egg .>= 2)

# ╔═╡ 6fb5fe92-9f05-4567-8d6e-046528336d27
chicken_ages = chain_egg[:A];

# ╔═╡ fe87f9ee-a07c-4b4c-b89d-c3bb4ee7ee7f
histogram(sp_egg[chicken_ages .== 1], normalize = :probability)

# ╔═╡ 9d9b4091-cf44-46e6-b81b-74663c9c8dfe
histogram(sp_egg[chicken_ages .== 3], normalize = :probability)

# ╔═╡ af08eb05-51a6-49a9-9e0b-15c2f88c6273
md"# 3: Petridish peril"

# ╔═╡ 45dfe42b-2274-40bc-bc1c-9903cd285ea1
md"""
Living the microbiology master thesis life, your mornings consist of inoculating petridishes with bacteria. Somewhere along the day, you need to split them before they have overgrown the entire dish and start dying. You want to estimate the risk of overgrowth occuring in function of time so you can get out of bed as late as possible without too much risk of disaster.
"""

# ╔═╡ 6f0b6153-b0fc-4ac0-81f3-044bd1211010
md"""
Bacteria follow **logistic growth**, and you can use the following assumptions:
- The initial population size $P_0$ has a 75% chance of originating from a small droplet and a 25% chance for a big droplet
  - Small droplets follow a `TriangularDist(1, 20, 5)`
  - Big droplets follow a `TriangularDist(10, 50, 20)`
- The growth rate $r$ follows a `Normal(1.0, 0.2)`
- The growth capacity $K$ of the inoculated medium follows a `LogNormal(log(1e5), 0.3)`
"""

# ╔═╡ 0c0a8f60-8a1f-4445-81d9-08343e6235a6
md"""
!!! questions
	- Plot the prior distribution of P0.
	- What is the probability of a petridish having over $10^5$ bacteria after 7, 8 and 9 hours of incubating?
	- Plot 100 of the sampled logistic growth curves from 0 to 12 hours.
"""

# ╔═╡ 67001598-27eb-4813-8465-3ffab01d84f4
md"""
!!! note
	A simple way of representing the distribution of P0 is through a mixture model. Look up the documentation of `MixtureModel` for how to make one in Turing.
"""

# ╔═╡ 57982520-9627-4a8a-911f-0d26f8fbf5f2
logistic(t, P0, r, K) =  K / (1 + (K - P0)/P0 * exp(-r*t))

# ╔═╡ 78c87e03-9ffc-4b27-a9c7-78172e045c1a
@model function petrigrowth(t)
    smalldropdist = TriangularDist(1, 20, 5)
	bigdropdist = TriangularDist(10, 50, 20)
	P0 ~ MixtureModel([smalldropdist, bigdropdist], [0.75, 0.25])
    r ~ Normal(1.0, 0.2)
	K ~ LogNormal(log(1e5), 0.3)

    Pt = logistic(t, P0, r, K)
    return Pt
end

# ╔═╡ 517577f2-9dc2-4b0f-aa3f-fd18f7f6fac3
plot(MixtureModel([TriangularDist(1, 20, 5), TriangularDist(10, 50, 20)]))

# ╔═╡ 06dc3fc4-08f7-4bac-82d9-b020a533eab6
petri_model = petrigrowth(8);

# ╔═╡ 465bd995-0d0c-4c68-bfb1-69cc9156c5d4
for t in [7, 8, 9]
	sp_petri = [petrigrowth(t)() for sample in 1:n_samples]
	prob_overgrown = mean(sp_petri .>= 1e5)
	println(prob_overgrown)
end

# ╔═╡ 7da4a9f2-9df9-45da-a387-22804155d09f
chain_petri = sample(petri_model, Prior(), n_samples)

# ╔═╡ e05ae9b2-08cc-4e15-b4f2-4adf3ad30cf2
sp_petri = generated_quantities(petri_model, chain_petri);

# ╔═╡ 9bf2e9c1-c777-41df-856c-2cdb427e175a
logistfuns = [
	t -> logistic(t, P0, r, K) 
	for (P0, r, K) in zip(chain_petri[:P0], chain_petri[:r], chain_petri[:K])
];

# ╔═╡ de17d13e-65a8-459b-9aae-4e02a61a7fd8
@time plot(logistfuns[1:100], xlims = (0, 12), legend = false, color = :skyblue, alpha = 0.5)

# ╔═╡ ff06c070-50a2-43d0-9729-1c47e728ff52
md"# 4: Birthdays"

# ╔═╡ 6ac2238a-16fd-4a8d-b779-8627d87367ed
md"""
Sometimes, people are born on the same day of year.
"""

# ╔═╡ 01648616-bf50-4f66-82fc-eaae3de22a38
md"""
!!! question
	What is the probability 3+ students in a class of 150 share a birthday?
"""

# ╔═╡ da44d18c-8be3-446e-a5c2-905af545d2c6
md"""
!!! note
	You can solve this (among other possibilities) using either a for-loop and the `count_occurences` function given below, or the `Multinomial` distribution.
"""

# ╔═╡ 52cf545a-d7c7-41d8-ad89-617d2f8b3eb9
count_occurences(vec) = [sum(vec .== element) for element in unique(vec)]

# ╔═╡ 0a6ff75d-fbc7-48a1-924b-e16d2654749c
count_occurences([5, 364, 107, 364, 2, 99])

# ╔═╡ da48e557-e3a5-47a6-89e7-4168eb23cf9d
@model function birthdays(n_students)
	bdays = zeros(n_students)
	for bday_idx in eachindex(bdays)
		bdays[bday_idx] ~ DiscreteUniform(1, 365)
	end
    occurences = count_occurences(bdays)
	
	# or 
    # occurences ~ Multinomial(n_students, 365)

    return maximum(occurences)
end

# ╔═╡ 44407741-8c74-4c00-a040-897c4713e6d7
bday_model = birthdays(150)

# ╔═╡ b0443045-74a1-4b48-8f6c-a1b5e15eace4
sp_maxoccs = [bday_model() for i in 1:n_samples]

# ╔═╡ cac00880-15fa-483a-a09f-9b6d1219b0cf
mean(sp_maxoccs .>= 3)

# ╔═╡ 8d66c150-4502-4501-980e-b2ce0eb79221
md"""
Using the party trick #!
"""

# ╔═╡ 0e5dba08-2bf2-4fd3-8179-cfb64410318f
P_samebday = 1/365

# ╔═╡ e1b2c14f-1da4-4f8b-bf14-cf731e42d110
possible_events = binomial(150, 3)

# ╔═╡ aa6723e5-a502-4a6c-8471-a6609a3ed944
expected_events = possible_events * P_samebday^2

# ╔═╡ 51a7f215-beff-4da0-861e-03326198ae3a
1 - pdf(Poisson(expected_events), 0)

# ╔═╡ 4e1a6396-2c20-4921-bf50-8a0fe446ccce
1-exp(-expected_events)

# ╔═╡ Cell order:
# ╠═7d4d4d20-b323-11ef-0926-b14785cb9ab5
# ╠═4cfd4721-e29a-4270-8d15-021bcc966eb1
# ╠═0fc3fa95-0287-423c-ba5a-7e382202ff81
# ╟─9740ea64-cd4f-46b1-a741-02e392280601
# ╟─187854bb-9e30-454d-9e03-cccf77aebb6b
# ╟─747e3c0a-357a-448a-b479-d0fcbe44a6c0
# ╟─eda95f45-083a-4e65-b57e-bd9890da1f9c
# ╠═36477423-5628-4ed3-b54d-9a050557f6b7
# ╠═02bc37e2-5cf6-404a-b3fe-b2120671adb2
# ╠═97021710-4b26-4a1c-a794-265c9559ce4d
# ╠═c0cd75bd-2a46-4f62-a59a-9bb8d47d45f2
# ╠═bc5f5d70-e269-45f9-867b-a00876ce8c40
# ╟─a1b933ac-5d1b-4800-a6e8-e942846b19d8
# ╠═6b009ba3-83a8-4176-86d4-dd9f70ed29ec
# ╠═83afb3c1-9e6a-4d18-b0c7-05ed0173df40
# ╠═700cc645-1123-4c8a-821e-1ae4fe8b0674
# ╠═98d04790-e86f-438a-9389-cb9697934bbf
# ╠═d5413573-5202-4fa2-94f6-e92a613d63aa
# ╟─34f3014f-f4d4-43d1-b46f-bdca73aee33f
# ╟─372436c4-262f-49b8-b1cf-626b043542bf
# ╟─20111742-008a-44c3-8c27-62791cce3e1e
# ╟─6e020801-983d-4ebc-a0e9-b5dd58f66c55
# ╠═2b3d930f-53d9-4869-9e4a-86a1a681b9d8
# ╠═834139c2-29bf-4882-bd73-1f89ad9f3803
# ╠═ed67949b-c7b3-46c2-af99-ab64dbef4065
# ╠═64d04ed9-6548-4251-8a4e-11b31e8f3143
# ╠═00fd5517-5232-4d88-8ab0-3a0eb925eb3e
# ╠═6fb5fe92-9f05-4567-8d6e-046528336d27
# ╠═fe87f9ee-a07c-4b4c-b89d-c3bb4ee7ee7f
# ╠═9d9b4091-cf44-46e6-b81b-74663c9c8dfe
# ╟─af08eb05-51a6-49a9-9e0b-15c2f88c6273
# ╟─45dfe42b-2274-40bc-bc1c-9903cd285ea1
# ╟─6f0b6153-b0fc-4ac0-81f3-044bd1211010
# ╟─0c0a8f60-8a1f-4445-81d9-08343e6235a6
# ╟─67001598-27eb-4813-8465-3ffab01d84f4
# ╠═57982520-9627-4a8a-911f-0d26f8fbf5f2
# ╠═78c87e03-9ffc-4b27-a9c7-78172e045c1a
# ╠═517577f2-9dc2-4b0f-aa3f-fd18f7f6fac3
# ╠═06dc3fc4-08f7-4bac-82d9-b020a533eab6
# ╠═465bd995-0d0c-4c68-bfb1-69cc9156c5d4
# ╠═7da4a9f2-9df9-45da-a387-22804155d09f
# ╠═e05ae9b2-08cc-4e15-b4f2-4adf3ad30cf2
# ╠═9bf2e9c1-c777-41df-856c-2cdb427e175a
# ╠═de17d13e-65a8-459b-9aae-4e02a61a7fd8
# ╟─ff06c070-50a2-43d0-9729-1c47e728ff52
# ╟─6ac2238a-16fd-4a8d-b779-8627d87367ed
# ╟─01648616-bf50-4f66-82fc-eaae3de22a38
# ╟─da44d18c-8be3-446e-a5c2-905af545d2c6
# ╠═52cf545a-d7c7-41d8-ad89-617d2f8b3eb9
# ╠═0a6ff75d-fbc7-48a1-924b-e16d2654749c
# ╠═da48e557-e3a5-47a6-89e7-4168eb23cf9d
# ╠═44407741-8c74-4c00-a040-897c4713e6d7
# ╠═b0443045-74a1-4b48-8f6c-a1b5e15eace4
# ╠═cac00880-15fa-483a-a09f-9b6d1219b0cf
# ╟─8d66c150-4502-4501-980e-b2ce0eb79221
# ╠═0e5dba08-2bf2-4fd3-8179-cfb64410318f
# ╠═e1b2c14f-1da4-4f8b-bf14-cf731e42d110
# ╠═aa6723e5-a502-4a6c-8471-a6609a3ed944
# ╠═51a7f215-beff-4da0-861e-03326198ae3a
# ╠═4e1a6396-2c20-4921-bf50-8a0fe446ccce
