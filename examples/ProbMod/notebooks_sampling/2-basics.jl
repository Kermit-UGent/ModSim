### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 7d4d4d20-b323-11ef-0926-b14785cb9ab5
using Pkg; Pkg.activate("..")

# ╔═╡ 4cfd4721-e29a-4270-8d15-021bcc966eb1
using Turing, StatsPlots

# ╔═╡ e4cb065e-12c6-4f1c-8497-1013fa9411d6
md"# Sampling notebook #2: Basics"

# ╔═╡ 9740ea64-cd4f-46b1-a741-02e392280601
md"""
## 1: Dice
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
mean(sp_w .>= 50)

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
mean(sp_d .>= 50)

# ╔═╡ d5413573-5202-4fa2-94f6-e92a613d63aa
histogram(sp_d, bins = 30)

# ╔═╡ 49790a8f-9f53-4ba7-9543-d6a879b520e0
md"### Comparison"

# ╔═╡ 6e098cb6-eddd-4924-a184-c2578dc28473
mean(sp_w .> sp_d)

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
- The amount of eggs $N$ an $A$-year old chicken lays in a year is Poisson distributed with mean $300 - 20 \, A$.
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
	A ~ Poisson(2)
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
Living the microbiology master thesis life, your mornings consist of inoculating petridishes with bacteria. Somewhere along the day, you need to split them. You want to do this **after** there's a decent amount of bacteria in the dish (>10\_000) but **before** they have overgrown the entire dish and start dying (<100\_000). This condition we call **splittable**.

You'd like to estimate how long after inoculation you should return to your bacteria so that they're most likely to be in a splittable state.
"""

# ╔═╡ 6f0b6153-b0fc-4ac0-81f3-044bd1211010
md"""
Bacteria follow **logistic growth**, and you can use the following assumptions:
- The initial population size $P_0$ has a 75% chance of originating from a small droplet and a 25% chance for a big droplet
  - For small droplets, `P0` follows a `Poisson(10)`
  - For big droplets, `P0` follows a `Poisson(30)`
- The growth rate $r$ follows a `LogNormal(0.0, 0.3)`
- The growth capacity $K$ of the inoculated medium follows a `Normal(1e5, 1e4)`
"""

# ╔═╡ 0c0a8f60-8a1f-4445-81d9-08343e6235a6
md"""
!!! questions
	1. Plot the prior distribution of P0.
	2. What is the probability your bacteria are in a splittable state 8 hours after inoculation?
	3. Plot 100 of the sampled logistic growth curves from 0 to 12 hours.
"""

# ╔═╡ 57982520-9627-4a8a-911f-0d26f8fbf5f2
logistic(t, P0, r, K) =  K / (1 + (K - P0)/P0 * exp(-r*t))

# ╔═╡ 8017759d-3d01-4195-b600-e2f15f18b34d
md"### 1"

# ╔═╡ 67001598-27eb-4813-8465-3ffab01d84f4
md"""
!!! tip
	A simple way of representing the distribution of P0 is through a mixture model. Look up the documentation of `MixtureModel` for how to make one in Turing.
"""

# ╔═╡ fb50a0cd-f0d7-42a2-ab5a-9f38ab133039
dropletdist = MixtureModel([Poisson(10), Poisson(30)], [0.75, 0.25]);

# ╔═╡ 517577f2-9dc2-4b0f-aa3f-fd18f7f6fac3
plot(dropletdist) # plot gives wrong result, the lines should stack

# ╔═╡ 1b06f0cb-77e4-4ee8-8a7d-731ecc36c72f
rand(dropletdist, 10000) |> histogram 
	# plotting a histogram of random samples gives a better result

# ╔═╡ b2ac8613-f14f-4187-bc60-f3ec884ef151
md"### 2"

# ╔═╡ 0dae6320-4678-41d1-94d8-bce20b7fe3c2
md"""
!!! tip
	You can `return` the logistic function estimated within the model and retrieve it using `generated_quantities` to make plotting easier later on.

	Remember: anonymous functions can be defined using `myfun(x) = ...`
"""

# ╔═╡ 78c87e03-9ffc-4b27-a9c7-78172e045c1a
@model function petrigrowth(t_obs)
	P0 ~ dropletdist
    r ~ LogNormal(0.0, 0.3)
	K ~ Normal(1e5, 1e4)

	logfun(t) = logistic(t, P0, r, K)
    Pt = logfun(t_obs)
    return logfun
end

# ╔═╡ 02b3485c-0b61-4427-ae50-a53f4537e20b
petri_model = petrigrowth(8.0);

# ╔═╡ 7da4a9f2-9df9-45da-a387-22804155d09f
chain_petri = sample(petri_model, Prior(), 2000);

# ╔═╡ 24e6daff-0708-4175-80b8-fb13595487c4
logfuns = generated_quantities(petri_model, chain_petri);

# ╔═╡ f492e7b0-7ee0-430a-8b31-6e82f8bf131d
sp_petri = [logfun(8.0) for logfun in logfuns]

# ╔═╡ 1018e455-12eb-4798-8fcf-3814b994ef51
prob_splittable = mean((sp_petri .>= 1e4) .&& (sp_petri .<= 1e5))

# ╔═╡ 05c8a2b3-91ce-47b2-8ef5-49c0b5db8aa0
md"### 3"

# ╔═╡ e4284ab8-869e-428f-ac6f-610f5b10f9ea
plot(logfuns[1:100], color = :violet, alpha = 0.3, label = false, xlims = (0, 12))

# ╔═╡ ff06c070-50a2-43d0-9729-1c47e728ff52
md"# 4: Birthdays"

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
sp_maxoccs = [bday_model() for i in 1:2000]

# ╔═╡ cac00880-15fa-483a-a09f-9b6d1219b0cf
mean(sp_maxoccs .>= 3)

# ╔═╡ Cell order:
# ╟─e4cb065e-12c6-4f1c-8497-1013fa9411d6
# ╠═7d4d4d20-b323-11ef-0926-b14785cb9ab5
# ╠═4cfd4721-e29a-4270-8d15-021bcc966eb1
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
# ╟─af08eb05-51a6-49a9-9e0b-15c2f88c6273
# ╟─45dfe42b-2274-40bc-bc1c-9903cd285ea1
# ╟─6f0b6153-b0fc-4ac0-81f3-044bd1211010
# ╟─0c0a8f60-8a1f-4445-81d9-08343e6235a6
# ╠═57982520-9627-4a8a-911f-0d26f8fbf5f2
# ╟─8017759d-3d01-4195-b600-e2f15f18b34d
# ╟─67001598-27eb-4813-8465-3ffab01d84f4
# ╠═fb50a0cd-f0d7-42a2-ab5a-9f38ab133039
# ╠═517577f2-9dc2-4b0f-aa3f-fd18f7f6fac3
# ╠═1b06f0cb-77e4-4ee8-8a7d-731ecc36c72f
# ╟─b2ac8613-f14f-4187-bc60-f3ec884ef151
# ╟─0dae6320-4678-41d1-94d8-bce20b7fe3c2
# ╠═78c87e03-9ffc-4b27-a9c7-78172e045c1a
# ╠═02b3485c-0b61-4427-ae50-a53f4537e20b
# ╠═7da4a9f2-9df9-45da-a387-22804155d09f
# ╠═24e6daff-0708-4175-80b8-fb13595487c4
# ╠═f492e7b0-7ee0-430a-8b31-6e82f8bf131d
# ╠═1018e455-12eb-4798-8fcf-3814b994ef51
# ╟─05c8a2b3-91ce-47b2-8ef5-49c0b5db8aa0
# ╠═e4284ab8-869e-428f-ac6f-610f5b10f9ea
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
