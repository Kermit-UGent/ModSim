### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 7d4d4d20-b323-11ef-0926-b14785cb9ab5
using Pkg; Pkg.activate("..")

# ╔═╡ 4cfd4721-e29a-4270-8d15-021bcc966eb1
using Turing, StatsPlots

# ╔═╡ c88b8a8e-751c-43fd-b1d9-cf6ae2089a15
using PlutoUI; TableOfContents()

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
md"### 1: Plots"

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
Y_samples = dpchain[:Y];

# ╔═╡ 90835b5d-9705-43d0-b449-09cde57d5394
plot(Poisson(10))

# ╔═╡ 431d9be9-edfc-48fb-9df7-6a2b961fe4b4
histogram(Y_samples)

# ╔═╡ 13989ba4-bcf8-4fdd-8aee-ab58c8905bc9
md"### 2: Probabilities"

# ╔═╡ fff486bd-0d3e-4bee-9264-acf7cb50b847
md"""
!!! tip
	When comparing a vector of values to a single number, don't forget to use `.` to execute operations element-wise in Julia!
	
	✅ `Y_samples .< 1` compares every element of `Y_samples` to `1`
	
	❌ `Y_samples < 1` compares an entire vector with a single number → errors :(
"""

# ╔═╡ 09a87037-9f4f-4dd9-bb6c-fe717db39ea6
probXY1 = mean(3 .< Y_samples .<= 10)

# ╔═╡ 2e2f203a-34a6-4feb-8aa1-8a1c97a09165
probXY2 = mean(Y_samples.^2 .> 100)

# ╔═╡ c243ca59-191d-4905-825e-6d7825a3c8a4
md"### 3: Variances"

# ╔═╡ 601cd057-f066-4d28-9345-ad955a7037e1
md"""
!!! hint
	To create a sample of $X$ that is conditional on some value(s) of $Y$, you can start from a sample of $X$ and select only those elements for which the corresponding sample of $Y$ has the conditioned value(s).

	In other words, you'll need to index `X_samples` based on `Y_samples` (and vice versa for $\text{var}(Y ∣ X)$). 
"""

# ╔═╡ 04cc5b42-7ab2-4055-a959-ba894c598f69
X_samples = dpchain[:X];

# ╔═╡ 597d97fb-cfbe-4cff-bfbb-57df9a2f42aa
varXcondY = var(X_samples[Y_samples .== 15])

# ╔═╡ 6ca31280-a997-4e84-98bb-e7784d0c555c
varYcondX = var(Y_samples[X_samples .== 15])

# ╔═╡ 9dcb122b-dd12-4e80-b391-f780e637d6fe
var(Poisson(15)) 
	# analytical answer of var(Y|X=15)
	# Given X, Y follows a Poisson distribution with mean X

# ╔═╡ c976dd77-6a72-4acc-9b07-342cc01fc7ee
md"""
!!! warning
	Starting here, the entire structure of the answer will no longer be given. You are therefore expected to **add more code cells** yourself, using the `+` at the left in between two cells.
"""

# ╔═╡ 9740ea64-cd4f-46b1-a741-02e392280601
md"""
## 2: Dice
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
    1. What is the probability that Watercube does the job? Also plot a histogram of its damage. 
    1. Do the same for Dirtprism.
	1. What is the probability that watercube deals more damage than dirtprism?
"""

# ╔═╡ 0e16af61-c4d4-4574-85e2-508f3952b3ce
md"### 1: Watercube"

# ╔═╡ 13cdd6ee-04ec-4085-8083-3d52e88aa35c
md"""
!!! tip
	Consider the humble `DiscreteUniform` distribution. Not sure how it works? Open the **🔍 Live Docs** at the bottom right of the screen for more information
"""

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
watercube_samples = [watermodel() for i in 1:2000]

# ╔═╡ c0cd75bd-2a46-4f62-a59a-9bb8d47d45f2
p_watercube_kills = mean(watercube_samples .>= 50)

# ╔═╡ bc5f5d70-e269-45f9-867b-a00876ce8c40
histogram(watercube_samples, bins = 30)

# ╔═╡ 44a2cabc-aca9-4bb7-af9d-cd735c424d7b
md"### 2: Dirtprism"

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
dirtcube_samples = [dirtmodel() for i in 1:2000]

# ╔═╡ 98d04790-e86f-438a-9389-cb9697934bbf
p_dirtprism_kills = mean(dirtcube_samples .>= 50)

# ╔═╡ d5413573-5202-4fa2-94f6-e92a613d63aa
histogram(dirtcube_samples, bins = 30)

# ╔═╡ 49790a8f-9f53-4ba7-9543-d6a879b520e0
md"### 3: Comparison"

# ╔═╡ 6e098cb6-eddd-4924-a184-c2578dc28473
p_watercube_is_better = mean(watercube_samples .> dirtcube_samples)

# ╔═╡ 34f3014f-f4d4-43d1-b46f-bdca73aee33f
md"## 3: Super eggs"

# ╔═╡ 372436c4-262f-49b8-b1cf-626b043542bf
md"""
When a chicken lays an egg, there's a small chance it contains two egg yolks. This chance, as well as the number of eggs a chicken lays per year, goes down as the chicken gets older. 
"""

# ╔═╡ 20111742-008a-44c3-8c27-62791cce3e1e
md"""
You can make the following assumptions
- The age $A$ of a random chicken (in years) is discrete and Uniformly distributed between 0 and 10.
- The number of eggs $N$ an $A$-year old chicken lays in a year is Poisson distributed with mean $300 - 20 \, A$.
- The probability $P$ of an $A$-year old chicken's egg having a double yolk is distributed as a `Beta(1, 800 + 100*A)`.
"""

# ╔═╡ 6e020801-983d-4ebc-a0e9-b5dd58f66c55
md"""
!!! questions
    1. If someone hands you a random chicken, what is the probability it will lay 2 or more double eggs in a year? 
    1. Compare the distributions of double eggs for 1-year old and 5-year old chickens.
"""

# ╔═╡ 6ebe676b-aee6-459b-bcda-de0fe14418a3
md"### 1: Probability"

# ╔═╡ 34aad3d6-9deb-4298-b613-bf3a616f2b31
md"""
!!! tip
	In this exercise, the output variable (the number of double-yolked eggs) is **also a random variable**! In other words, it also follows some distribution.

	When considering what distribution, consider that each of the $N$ eggs represents a "trial" with a $P$ chance of success for a double yolk.
"""

# ╔═╡ 2b3d930f-53d9-4869-9e4a-86a1a681b9d8
@model function eggs()
	A ~ DiscreteUniform(0, 10)
	N ~ Poisson(300 - 20*A)
	P ~ Beta(1, 800 + 100*A)

	doubles ~ Binomial(N, P)
	return doubles
end

# ╔═╡ 834139c2-29bf-4882-bd73-1f89ad9f3803
egg_model = eggs();

# ╔═╡ ed67949b-c7b3-46c2-af99-ab64dbef4065
chain_egg = sample(egg_model, Prior(), 2000)

# ╔═╡ 64d04ed9-6548-4251-8a4e-11b31e8f3143
egg_samples = generated_quantities(egg_model, chain_egg);

# ╔═╡ 00fd5517-5232-4d88-8ab0-3a0eb925eb3e
p_multiple_double_eggs = mean(egg_samples .>= 2) 
	# a ~2% chance for two or more double eggs in one year

# ╔═╡ 80e3544e-86ec-469c-be44-77f94fabfc28
md"### 2: Histograms"

# ╔═╡ 6fb5fe92-9f05-4567-8d6e-046528336d27
chicken_ages = chain_egg[:A];

# ╔═╡ fe87f9ee-a07c-4b4c-b89d-c3bb4ee7ee7f
histogram(egg_samples[chicken_ages .== 1], normalize = :probability)

# ╔═╡ 9d9b4091-cf44-46e6-b81b-74663c9c8dfe
histogram(egg_samples[chicken_ages .== 5], normalize = :probability)

# ╔═╡ ff06c070-50a2-43d0-9729-1c47e728ff52
md"## 4: Birthdays"

# ╔═╡ 6ac2238a-16fd-4a8d-b779-8627d87367ed
md"""
Sometimes, people are born on the same day of the year.
"""

# ╔═╡ 01648616-bf50-4f66-82fc-eaae3de22a38
md"""
!!! question
	What is the probability that, in a class of 150 students, 3 or more share a birthday? Assume the probability for a person to be born is equal on every day of the year.
"""

# ╔═╡ 29438eef-a725-4873-9372-194f2f899f12
md"""
!!! tip
	You can solve this (among other possibilities) using either a for-loop and the `count_occurences` function given below, or the `Multinomial` distribution.
"""

# ╔═╡ 52cf545a-d7c7-41d8-ad89-617d2f8b3eb9
count_occurences(vec) = [count(==(element), vec) for element in unique(vec)]

# ╔═╡ 0a6ff75d-fbc7-48a1-924b-e16d2654749c
count_occurences([5, 107, 364, 5, 5, 364]) # three 5's, one 107 and two 364's

# ╔═╡ da48e557-e3a5-47a6-89e7-4168eb23cf9d
@model function birthdays(n_students)
	# number of students used as an input to the function so you can solve the question for any number of students - this was not asked but it's nice
	
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
bday_samples = [bday_model() for i in 1:2000]

# ╔═╡ cac00880-15fa-483a-a09f-9b6d1219b0cf
mean(bday_samples .>= 3)

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
# ╠═c88b8a8e-751c-43fd-b1d9-cf6ae2089a15
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
# ╟─fff486bd-0d3e-4bee-9264-acf7cb50b847
# ╠═09a87037-9f4f-4dd9-bb6c-fe717db39ea6
# ╠═2e2f203a-34a6-4feb-8aa1-8a1c97a09165
# ╟─c243ca59-191d-4905-825e-6d7825a3c8a4
# ╟─601cd057-f066-4d28-9345-ad955a7037e1
# ╠═04cc5b42-7ab2-4055-a959-ba894c598f69
# ╠═597d97fb-cfbe-4cff-bfbb-57df9a2f42aa
# ╠═6ca31280-a997-4e84-98bb-e7784d0c555c
# ╠═9dcb122b-dd12-4e80-b391-f780e637d6fe
# ╟─c976dd77-6a72-4acc-9b07-342cc01fc7ee
# ╟─9740ea64-cd4f-46b1-a741-02e392280601
# ╟─187854bb-9e30-454d-9e03-cccf77aebb6b
# ╟─747e3c0a-357a-448a-b479-d0fcbe44a6c0
# ╟─0e16af61-c4d4-4574-85e2-508f3952b3ce
# ╟─13cdd6ee-04ec-4085-8083-3d52e88aa35c
# ╠═36477423-5628-4ed3-b54d-9a050557f6b7
# ╠═02bc37e2-5cf6-404a-b3fe-b2120671adb2
# ╠═97021710-4b26-4a1c-a794-265c9559ce4d
# ╠═c0cd75bd-2a46-4f62-a59a-9bb8d47d45f2
# ╠═bc5f5d70-e269-45f9-867b-a00876ce8c40
# ╟─44a2cabc-aca9-4bb7-af9d-cd735c424d7b
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
# ╟─6ebe676b-aee6-459b-bcda-de0fe14418a3
# ╟─34aad3d6-9deb-4298-b613-bf3a616f2b31
# ╠═2b3d930f-53d9-4869-9e4a-86a1a681b9d8
# ╠═834139c2-29bf-4882-bd73-1f89ad9f3803
# ╠═ed67949b-c7b3-46c2-af99-ab64dbef4065
# ╠═64d04ed9-6548-4251-8a4e-11b31e8f3143
# ╠═00fd5517-5232-4d88-8ab0-3a0eb925eb3e
# ╟─80e3544e-86ec-469c-be44-77f94fabfc28
# ╠═6fb5fe92-9f05-4567-8d6e-046528336d27
# ╠═fe87f9ee-a07c-4b4c-b89d-c3bb4ee7ee7f
# ╠═9d9b4091-cf44-46e6-b81b-74663c9c8dfe
# ╟─ff06c070-50a2-43d0-9729-1c47e728ff52
# ╟─6ac2238a-16fd-4a8d-b779-8627d87367ed
# ╟─01648616-bf50-4f66-82fc-eaae3de22a38
# ╟─29438eef-a725-4873-9372-194f2f899f12
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
