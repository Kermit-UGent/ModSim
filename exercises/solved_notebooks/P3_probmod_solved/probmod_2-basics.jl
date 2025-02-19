### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 7d4d4d20-b323-11ef-0926-b14785cb9ab5
using Pkg; Pkg.activate("..")

# ╔═╡ 4cfd4721-e29a-4270-8d15-021bcc966eb1
using Turing, StatsPlots

# ╔═╡ e4cb065e-12c6-4f1c-8497-1013fa9411d6
md"# Sampling notebook #2: Basics"

# ╔═╡ 7026f66f-9076-4aef-ada9-198450ef5da6
md"## 1: Double Poisson"

# ╔═╡ bb682b51-c3ac-4e31-9b79-4c13212d84e5
md"""
Let `X ∼ Poisson(10)` and `Y ~ Poisson(X)`.

1. Plot the (exact) distribution of `X` and use sampling (n = 10_000) to generate a histogram of `Y`.
2. Estimate the following probabilities:
    - `P(3 < Y ≤ 10)`.
    - `P(Y^2 > 100)`.
3. Consider `var(X|Y=15)` and `var(Y|X=15)`. 
    - Estimate them numerically.
    - One of the two has a simple analytical answer: which one is it, and what is its exact value?
"""

# ╔═╡ def27d26-2205-4a66-94f3-eddbc17483bf
md"### 1"

# ╔═╡ ff38df99-f843-414d-8e45-b46e06a65c22
@model function doublepoisson()
	X ~ Poisson(10)
	Y ~ Poisson(X)
	return Y
end

# ╔═╡ 7bdaeebd-eb98-44fb-891d-da6959323474
dpmodel = doublepoisson()

# ╔═╡ a9a0bf40-3b01-4da1-b881-55978a6fa80b
dpchain = sample(dpmodel, Prior(), 10_000)

# ╔═╡ f6d90b53-0f16-4e43-b77c-8eb05b8c6943
spY = dpchain[:Y];

# ╔═╡ 90835b5d-9705-43d0-b449-09cde57d5394
plot(Poisson(10))

# ╔═╡ 431d9be9-edfc-48fb-9df7-6a2b961fe4b4
histogram(spY)

# ╔═╡ 13989ba4-bcf8-4fdd-8aee-ab58c8905bc9
md"### 2"

# ╔═╡ 09a87037-9f4f-4dd9-bb6c-fe717db39ea6
probXY1 = mean(3 .< spY .<= 10)

# ╔═╡ 2e2f203a-34a6-4feb-8aa1-8a1c97a09165
probXY2 = mean(spY.^2 .> 100)

# ╔═╡ c243ca59-191d-4905-825e-6d7825a3c8a4
md"### 3"

# ╔═╡ 04cc5b42-7ab2-4055-a959-ba894c598f69
spX = dpchain[:X];

# ╔═╡ 597d97fb-cfbe-4cff-bfbb-57df9a2f42aa
varXcondY = var(spX[spY .== 15])

# ╔═╡ 6ca31280-a997-4e84-98bb-e7784d0c555c
varYcondX = var(spY[spX .== 15])

# ╔═╡ 9dcb122b-dd12-4e80-b391-f780e637d6fe
var(Poisson(15)) 
	# analytical answer of var(Y|X=15)
	# Given X, Y follows a Poisson distribution with mean X

# ╔═╡ ce7d57ed-4f31-4dcf-af3b-b37a2e2a9393
md"## 2: Combinations"

# ╔═╡ 087994ce-a26d-40c4-87eb-ef9f0ce7f1fb
md"""
Let `U ~ Uniform(0, 4)`, `V ∼ Normal(U, 1)` and `W ~ TriangularDist(0, 4, U)`. 
1. Use sampling (n = 10_000) to make a histogram of `|V − W|`.
2. Estimate `P(V > W)` and `P(V * W >= 10)`
3. Are `V` and `W` independent?
"""

# ╔═╡ 9e5cc347-c74f-46a3-9534-c5ad812844bf
md"### 1"

# ╔═╡ 47a43282-3892-4a9a-94b7-c359fa74e12b
@model function combinations()
	U ~ Uniform(0, 4)
	V ~ Normal(U, 1.0)
	W ~ TriangularDist(0, 4, U)
end

# ╔═╡ 890eaeb5-c63f-470c-85eb-d77cec68b740
comb_model = combinations()

# ╔═╡ 0d0297a3-60b6-47e5-b42b-bde99c09a16e
comb_chain = sample(comb_model, Prior(), 10_000);

# ╔═╡ c3cec28a-f472-46f5-9082-ddfdeb01ddf0
spV = comb_chain[:V];

# ╔═╡ 1911f031-d668-47df-995a-ca339ad4ae47
spW = comb_chain[:W];

# ╔═╡ 5c13b506-c05b-4258-95e7-93b648defa17
spVW = abs.(spV - spW)

# ╔═╡ d8ae7b0a-afee-41c0-af0c-f8b1d30086f9
histogram(spVW)

# ╔═╡ faa105b6-0700-4f4d-92fe-2bb72a4d6e44
md"### 2"

# ╔═╡ b4f156df-0b9c-495a-b95b-000c7f166243
probVW1 = mean(spV .> spW)

# ╔═╡ f7eb6faf-8884-4fea-9974-330d7fd555aa
probVW2 = mean(spV .* spW .>= 10)

# ╔═╡ 4decd959-aeb9-47d4-a381-14bbf4dbc5ab
md"### 3"

# ╔═╡ aa953baa-5105-49d7-82e6-94ca462624f7
md"""
!!! hint
	One way to disprove independence is showing that $E[V] \neq E[V \mid W \leq w]$ for any value $w$.
"""

# ╔═╡ b2999f6e-2dbd-44ab-a303-97234f665d33
mean(spV)

# ╔═╡ 5ea80b00-941a-4226-87cb-66da7ee7d76d
mean(spV[spW .<= 1]) # E[V] != E[V | W <= 1] => they are dependent

# ╔═╡ 9740ea64-cd4f-46b1-a741-02e392280601
md"""
## 3: Dice
"""

# ╔═╡ 187854bb-9e30-454d-9e03-cccf77aebb6b
md"You're playing a fun game of Caverns and Chimeras, and are facing off against the mighty Carl the Chimera. The fight is not going great and your next spell **needs to deal 50 or more damage** to slay the scary monster before it kills you. Spells deal **damage equal to the sum of the dice** they let you roll.

You can choose between your 2 mightiest spells:
- **Watercube**: lets you throw **4 dice with 20 sides** each.
- **Dirtprism**: lets you throw **20 dice with 4 sides** each.
"

# ╔═╡ 747e3c0a-357a-448a-b479-d0fcbe44a6c0
md"""
!!! questions
    - What are the probabilities of either attack doing the job? 
    - Plot a histogram of the damage of both attacks.
	- What is the probability of watercube dealing more damage than dirtprism?
"""

# ╔═╡ eda95f45-083a-4e65-b57e-bd9890da1f9c
md"### 4d20"

# ╔═╡ 36477423-5628-4ed3-b54d-9a050557f6b7
@model function watercube()
    roll1 ~ DiscreteUniform(1, 20)
	roll2 ~ DiscreteUniform(1, 20)
	roll3 ~ DiscreteUniform(1, 20)
	roll4 ~ DiscreteUniform(1, 20)
	# or use a loop

    dicesum = roll1 + roll2 + roll3 + roll4
    return dicesum
end

# ╔═╡ 02bc37e2-5cf6-404a-b3fe-b2120671adb2
watermodel = watercube()

# ╔═╡ 97021710-4b26-4a1c-a794-265c9559ce4d
sp_w = [watermodel() for i in 1:2000]

# ╔═╡ c0cd75bd-2a46-4f62-a59a-9bb8d47d45f2
p_watercube_kills = mean(sp_w .>= 50)

# ╔═╡ bc5f5d70-e269-45f9-867b-a00876ce8c40
histogram(sp_w, bins = 30)

# ╔═╡ a1b933ac-5d1b-4800-a6e8-e942846b19d8
md"### 20d4"

# ╔═╡ 6b009ba3-83a8-4176-86d4-dd9f70ed29ec
@model function dirtprism()
    rolls = zeros(20) # also possible to write out all 20 rolls by hand
	for i in eachindex(rolls)
		rolls[i] ~ DiscreteUniform(1, 4)
	end
    dicesum = sum(rolls)
    return dicesum
end

# ╔═╡ 83afb3c1-9e6a-4d18-b0c7-05ed0173df40
dirtmodel = dirtprism()

# ╔═╡ 700cc645-1123-4c8a-821e-1ae4fe8b0674
sp_d = [dirtmodel() for i in 1:2000]

# ╔═╡ 98d04790-e86f-438a-9389-cb9697934bbf
p_dirtprism_kills = mean(sp_d .>= 50)

# ╔═╡ d5413573-5202-4fa2-94f6-e92a613d63aa
histogram(sp_d, bins = 30)

# ╔═╡ 49790a8f-9f53-4ba7-9543-d6a879b520e0
md"### Comparison"

# ╔═╡ 6e098cb6-eddd-4924-a184-c2578dc28473
p_watercube_is_better = mean(sp_w .> sp_d)

# ╔═╡ 34f3014f-f4d4-43d1-b46f-bdca73aee33f
md"## 4: Super eggs"

# ╔═╡ 372436c4-262f-49b8-b1cf-626b043542bf
md"""
When a chicken lays an egg, there's a small chance it contains two egg yolks. This chance, as well as the number of eggs a chicken lays per year, go down as the chicken gets older. 
"""

# ╔═╡ 20111742-008a-44c3-8c27-62791cce3e1e
md"""
You can make the following assumptions
- The age $A$ of a random chicken (in years) is discrete and Uniformly distributed between 0 and 12.
- The number of eggs $N$ an $A$-year old chicken lays in a year is Poisson distributed with mean $300 - 20 \, A$.
- The probability of an $A$-year old chicken's egg having a double yolk $P$ is distributed as a `Beta(1, 800 + 100*A)`.
"""

# ╔═╡ 6e020801-983d-4ebc-a0e9-b5dd58f66c55
md"""
!!! questions
    - If someone hands you a random chicken, what is the probability it will lay 2 or more double eggs in a year? 
    - Compare the distributions of double eggs for 1-year old and 3-year old chickens.
"""

# ╔═╡ 2b3d930f-53d9-4869-9e4a-86a1a681b9d8
@model function eggs()
	A ~ DiscreteUniform(0, 12)
	N ~ Poisson(300 - 20*A)
	P ~ Beta(1, 800 + 100*A)

	doubles ~ Binomial(N, P)
	return sum(doubles)
end

# ╔═╡ 834139c2-29bf-4882-bd73-1f89ad9f3803
egg_model = eggs();

# ╔═╡ ed67949b-c7b3-46c2-af99-ab64dbef4065
chain_egg = sample(egg_model, Prior(), 2000)

# ╔═╡ 64d04ed9-6548-4251-8a4e-11b31e8f3143
sp_egg = generated_quantities(egg_model, chain_egg);

# ╔═╡ 00fd5517-5232-4d88-8ab0-3a0eb925eb3e
p_multiple_double_eggs = mean(sp_egg .>= 2) 
	# a ~2% chance for two or more double eggs in one year

# ╔═╡ 6fb5fe92-9f05-4567-8d6e-046528336d27
chicken_ages = chain_egg[:A];

# ╔═╡ fe87f9ee-a07c-4b4c-b89d-c3bb4ee7ee7f
histogram(sp_egg[chicken_ages .== 1], normalize = :probability)

# ╔═╡ 9d9b4091-cf44-46e6-b81b-74663c9c8dfe
histogram(sp_egg[chicken_ages .== 3], normalize = :probability)

# ╔═╡ ff06c070-50a2-43d0-9729-1c47e728ff52
md"## 5: Birthdays"

# ╔═╡ 6ac2238a-16fd-4a8d-b779-8627d87367ed
md"""
Sometimes, people are born on the same day of the year.
"""

# ╔═╡ 01648616-bf50-4f66-82fc-eaae3de22a38
md"""
!!! question
	What is the probability that, in a class of 150 students, 3 or more share a birthday?
"""

# ╔═╡ da44d18c-8be3-446e-a5c2-905af545d2c6
md"""
!!! note
	You can solve this (among other possibilities) using either a for-loop and the `count_occurences` function given below, or the `Multinomial` distribution.
"""

# ╔═╡ 52cf545a-d7c7-41d8-ad89-617d2f8b3eb9
count_occurences(vec) = [count(==(element), vec) for element in unique(vec)]

# ╔═╡ 0a6ff75d-fbc7-48a1-924b-e16d2654749c
count_occurences([5, 107, 364, 5, 5, 364]) # three 5's, one 107 and two 364's

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
sp_maxoccs = [bday_model() for i in 1:2000]

# ╔═╡ cac00880-15fa-483a-a09f-9b6d1219b0cf
mean(sp_maxoccs .>= 3)

# ╔═╡ 13002efe-15f1-4096-be4f-671432a8991e
md"""
!!! extra
	When calculating the probability of multiple rare events occuring within a given timespan, we can quickly estimate that probability as follows:
	- Calculate the probability that a single rare event occurs.
	- Calculate the expected amount of rare events in the given timespan, $E$.
	- The probability at least one of the rare events occurs is $1 - e^{-E}$.
"""

# ╔═╡ e59fe76a-b663-468e-8f0a-bce727f0fa27
p_triplebday = 1 - cdf(Binomial(150, 1/365), 2)
	# On one day you have 150 students (trials) with a 1/365 success rate for it to be their birthday. We can model the amount of birthdays on one day "B" with a Binomial distribution. Then P(B >= 3) = 1 - P(B <= 2) = 1 - cdf(Binomial(n, p), 2)

# ╔═╡ 805f44ee-11ca-4521-a689-1fce995b334a
E_triplebday = p_triplebday * 365 # The expected value of the amount of triple birthdays over 365 days

# ╔═╡ 40071507-bbc1-404b-8ed3-7333b5f1854c
1 - exp(-E_triplebday)

# ╔═╡ Cell order:
# ╟─e4cb065e-12c6-4f1c-8497-1013fa9411d6
# ╠═7d4d4d20-b323-11ef-0926-b14785cb9ab5
# ╠═4cfd4721-e29a-4270-8d15-021bcc966eb1
# ╟─7026f66f-9076-4aef-ada9-198450ef5da6
# ╟─bb682b51-c3ac-4e31-9b79-4c13212d84e5
# ╟─def27d26-2205-4a66-94f3-eddbc17483bf
# ╠═ff38df99-f843-414d-8e45-b46e06a65c22
# ╠═7bdaeebd-eb98-44fb-891d-da6959323474
# ╠═a9a0bf40-3b01-4da1-b881-55978a6fa80b
# ╠═f6d90b53-0f16-4e43-b77c-8eb05b8c6943
# ╠═90835b5d-9705-43d0-b449-09cde57d5394
# ╠═431d9be9-edfc-48fb-9df7-6a2b961fe4b4
# ╟─13989ba4-bcf8-4fdd-8aee-ab58c8905bc9
# ╠═09a87037-9f4f-4dd9-bb6c-fe717db39ea6
# ╠═2e2f203a-34a6-4feb-8aa1-8a1c97a09165
# ╟─c243ca59-191d-4905-825e-6d7825a3c8a4
# ╠═04cc5b42-7ab2-4055-a959-ba894c598f69
# ╠═597d97fb-cfbe-4cff-bfbb-57df9a2f42aa
# ╠═6ca31280-a997-4e84-98bb-e7784d0c555c
# ╠═9dcb122b-dd12-4e80-b391-f780e637d6fe
# ╟─ce7d57ed-4f31-4dcf-af3b-b37a2e2a9393
# ╟─087994ce-a26d-40c4-87eb-ef9f0ce7f1fb
# ╟─9e5cc347-c74f-46a3-9534-c5ad812844bf
# ╠═47a43282-3892-4a9a-94b7-c359fa74e12b
# ╠═890eaeb5-c63f-470c-85eb-d77cec68b740
# ╠═0d0297a3-60b6-47e5-b42b-bde99c09a16e
# ╠═c3cec28a-f472-46f5-9082-ddfdeb01ddf0
# ╠═1911f031-d668-47df-995a-ca339ad4ae47
# ╠═5c13b506-c05b-4258-95e7-93b648defa17
# ╠═d8ae7b0a-afee-41c0-af0c-f8b1d30086f9
# ╟─faa105b6-0700-4f4d-92fe-2bb72a4d6e44
# ╠═b4f156df-0b9c-495a-b95b-000c7f166243
# ╠═f7eb6faf-8884-4fea-9974-330d7fd555aa
# ╟─4decd959-aeb9-47d4-a381-14bbf4dbc5ab
# ╟─aa953baa-5105-49d7-82e6-94ca462624f7
# ╠═b2999f6e-2dbd-44ab-a303-97234f665d33
# ╠═5ea80b00-941a-4226-87cb-66da7ee7d76d
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
# ╟─49790a8f-9f53-4ba7-9543-d6a879b520e0
# ╠═6e098cb6-eddd-4924-a184-c2578dc28473
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
# ╟─13002efe-15f1-4096-be4f-671432a8991e
# ╠═e59fe76a-b663-468e-8f0a-bce727f0fa27
# ╠═805f44ee-11ca-4521-a689-1fce995b334a
# ╠═40071507-bbc1-404b-8ed3-7333b5f1854c
