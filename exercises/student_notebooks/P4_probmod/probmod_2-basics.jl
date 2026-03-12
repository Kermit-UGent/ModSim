### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 7d4d4d20-b323-11ef-0926-b14785cb9ab5
using Pkg; Pkg.activate("..")

# ╔═╡ 4cfd4721-e29a-4270-8d15-021bcc966eb1
using Turing, StatsPlots

# ╔═╡ 73c2a5db-4019-4d91-b5f1-7ba378cb8c84
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
	X ~ missing
	Y ~ missing
	return Y
end

# ╔═╡ 42c18a70-efb3-436b-83bb-b586280d4a4e
dpmodel = doublepoisson();

# ╔═╡ b22a18d5-f70b-42a5-a09e-3de515148a6d
spY = missing

# ╔═╡ 34bd60df-272f-49d7-9346-fb4d125fe89b
histogram(spY)

# ╔═╡ 13989ba4-bcf8-4fdd-8aee-ab58c8905bc9
md"### 2: Probabilities"

# ╔═╡ aba42086-224f-44a0-b616-f8f651afdd18
md"""
!!! tip
	When comparing a vector of values to a single number, don't forget to use `.` to execute operations element-wise in Julia!
	
	✅ `spY .< 1` compares every element of `spY` to `1`
	
	❌ `spY < 1` compares an entire vector with a single number → errors :(
"""

# ╔═╡ 39e4f3eb-b7a9-4ece-846f-eb02dbd77860
probXY1 = missing

# ╔═╡ acca0e39-612f-43b8-9c25-74c143041978
probXY2 = missing

# ╔═╡ c243ca59-191d-4905-825e-6d7825a3c8a4
md"### 3: Variances"

# ╔═╡ 137a1727-08ab-4455-ad8e-88bbead80845
md"""
!!! hint
	To create a sample of $X$ that is conditional on some value(s) of $Y$, you can start from a sample of $X$ and select only those elements for which the corresponding sample of $Y$ has the conditioned value(s).

	In other words, you'll need to index `spX` based on `spY` (and vice versa for $\text{var}(Y ∣ X)$). 
"""

# ╔═╡ 263048e5-206c-499e-836c-bfe489ed9b74
spX = missing;

# ╔═╡ aedd0fe8-da3e-4463-b0cf-7c4f9a22db52
varXcondY = missing

# ╔═╡ 7ef53a87-e5df-4724-b896-3d1d46214c68
varYcondX = missing

# ╔═╡ ce9e2ce2-26e0-4f17-adf5-c922ba98239d
missing # analytical answer of (missing)

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

# ╔═╡ eda95f45-083a-4e65-b57e-bd9890da1f9c
md"### 1: Watercube"

# ╔═╡ 4fefa78b-746d-4e25-baee-d791eb22d930
md"""
!!! hint
	Consider the humble `DiscreteUniform` distribution. Not sure how it works? Open the **🔍 Live Docs** at the bottom right of the screen for more information
"""

# ╔═╡ 36477423-5628-4ed3-b54d-9a050557f6b7
@model function watercube()
    roll1 ~ missing
	roll2 ~ missing
	roll3 ~ missing
	roll4 ~ missing

    dicesum = roll1 + roll2 + roll3 + roll4
    return dicesum
end

# ╔═╡ 02bc37e2-5cf6-404a-b3fe-b2120671adb2
watermodel = watercube()

# ╔═╡ b5886255-5c1d-4d84-b7ee-6c690fa526dc
p_watercube_kills = missing

# ╔═╡ 84e06162-e3f1-4fd7-baf6-5095172413d2
missing # histogram

# ╔═╡ a1b933ac-5d1b-4800-a6e8-e942846b19d8
md"### 2: Dirtprism"

# ╔═╡ 6b009ba3-83a8-4176-86d4-dd9f70ed29ec
@model function dirtprism()
	# check the "For-loops for many variables" section from the intro notebook!
	
	dicesum = missing # consider the `sum` function
	return dicesum
end

# ╔═╡ 83afb3c1-9e6a-4d18-b0c7-05ed0173df40
dirtmodel = dirtprism()

# ╔═╡ d9152416-8a7b-480c-ba9f-7ab15404b7a6
p_dirtprism_kills = missing

# ╔═╡ 90e058d7-b3fe-4c42-a652-3c42bf9d851a
missing # histogram

# ╔═╡ 49790a8f-9f53-4ba7-9543-d6a879b520e0
md"### 3: Comparison"

# ╔═╡ 9271b6df-fa3e-4b79-8f0f-48a3b1287b42
p_watercube_is_better = missing

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

# ╔═╡ 98fcfcef-bf63-4eae-a325-ed4cef6d4fdd
md"### 1: Probability"

# ╔═╡ a0341046-c14a-494d-a9bd-a60c207c9e76
md"""
!!! hint
	In this exercise, the output variable (the number of double-yolked eggs) is **also a random variable**! In other words, it also follows some distribution.

	When considering what distribution, consider that each of the $N$ eggs represents a "trial" with a $P$ chance of success for a double yolk.
"""

# ╔═╡ 2b3d930f-53d9-4869-9e4a-86a1a681b9d8
@model function eggs()
	return missing
end

# ╔═╡ 49d274f7-5810-48b1-8954-22b6a0941a47
p_multiple_double_eggs = missing

# ╔═╡ 8dcdc5d7-07b5-4041-b955-485b6f830b75
md"### 2: Histograms"

# ╔═╡ cee7c02a-62db-4181-8321-b8bea8fb9339
missing # histogram 1

# ╔═╡ 6b3227c5-78aa-4031-b321-f938757d5ad8
missing # histogram 2

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

# ╔═╡ da44d18c-8be3-446e-a5c2-905af545d2c6
md"""
!!! tip
	You can solve this (among other possibilities) using either a for-loop and the `count_occurences` function given below, or the `Multinomial` distribution.
"""

# ╔═╡ 52cf545a-d7c7-41d8-ad89-617d2f8b3eb9
count_occurences(vec) = [count(==(element), vec) for element in unique(vec)]

# ╔═╡ 0a6ff75d-fbc7-48a1-924b-e16d2654749c
count_occurences([5, 107, 364, 5, 5, 364]) # three 5's, one 107 and two 364's

# ╔═╡ 2925bb92-32f8-415d-960f-b0a60a1037b8
@model function birthdays()
	missing
end

# ╔═╡ Cell order:
# ╟─e4cb065e-12c6-4f1c-8497-1013fa9411d6
# ╠═7d4d4d20-b323-11ef-0926-b14785cb9ab5
# ╠═4cfd4721-e29a-4270-8d15-021bcc966eb1
# ╠═73c2a5db-4019-4d91-b5f1-7ba378cb8c84
# ╟─7026f66f-9076-4aef-ada9-198450ef5da6
# ╟─bb682b51-c3ac-4e31-9b79-4c13212d84e5
# ╟─def27d26-2205-4a66-94f3-eddbc17483bf
# ╠═ff38df99-f843-414d-8e45-b46e06a65c22
# ╠═42c18a70-efb3-436b-83bb-b586280d4a4e
# ╠═b22a18d5-f70b-42a5-a09e-3de515148a6d
# ╠═34bd60df-272f-49d7-9346-fb4d125fe89b
# ╟─13989ba4-bcf8-4fdd-8aee-ab58c8905bc9
# ╟─aba42086-224f-44a0-b616-f8f651afdd18
# ╠═39e4f3eb-b7a9-4ece-846f-eb02dbd77860
# ╠═acca0e39-612f-43b8-9c25-74c143041978
# ╟─c243ca59-191d-4905-825e-6d7825a3c8a4
# ╟─137a1727-08ab-4455-ad8e-88bbead80845
# ╠═263048e5-206c-499e-836c-bfe489ed9b74
# ╠═aedd0fe8-da3e-4463-b0cf-7c4f9a22db52
# ╠═7ef53a87-e5df-4724-b896-3d1d46214c68
# ╠═ce9e2ce2-26e0-4f17-adf5-c922ba98239d
# ╟─9740ea64-cd4f-46b1-a741-02e392280601
# ╟─187854bb-9e30-454d-9e03-cccf77aebb6b
# ╟─747e3c0a-357a-448a-b479-d0fcbe44a6c0
# ╟─eda95f45-083a-4e65-b57e-bd9890da1f9c
# ╟─4fefa78b-746d-4e25-baee-d791eb22d930
# ╠═36477423-5628-4ed3-b54d-9a050557f6b7
# ╠═02bc37e2-5cf6-404a-b3fe-b2120671adb2
# ╠═b5886255-5c1d-4d84-b7ee-6c690fa526dc
# ╠═84e06162-e3f1-4fd7-baf6-5095172413d2
# ╟─a1b933ac-5d1b-4800-a6e8-e942846b19d8
# ╠═6b009ba3-83a8-4176-86d4-dd9f70ed29ec
# ╠═83afb3c1-9e6a-4d18-b0c7-05ed0173df40
# ╠═d9152416-8a7b-480c-ba9f-7ab15404b7a6
# ╠═90e058d7-b3fe-4c42-a652-3c42bf9d851a
# ╟─49790a8f-9f53-4ba7-9543-d6a879b520e0
# ╠═9271b6df-fa3e-4b79-8f0f-48a3b1287b42
# ╟─34f3014f-f4d4-43d1-b46f-bdca73aee33f
# ╟─372436c4-262f-49b8-b1cf-626b043542bf
# ╟─20111742-008a-44c3-8c27-62791cce3e1e
# ╟─6e020801-983d-4ebc-a0e9-b5dd58f66c55
# ╟─98fcfcef-bf63-4eae-a325-ed4cef6d4fdd
# ╟─a0341046-c14a-494d-a9bd-a60c207c9e76
# ╠═2b3d930f-53d9-4869-9e4a-86a1a681b9d8
# ╠═49d274f7-5810-48b1-8954-22b6a0941a47
# ╟─8dcdc5d7-07b5-4041-b955-485b6f830b75
# ╠═cee7c02a-62db-4181-8321-b8bea8fb9339
# ╠═6b3227c5-78aa-4031-b321-f938757d5ad8
# ╟─ff06c070-50a2-43d0-9729-1c47e728ff52
# ╟─6ac2238a-16fd-4a8d-b779-8627d87367ed
# ╟─01648616-bf50-4f66-82fc-eaae3de22a38
# ╟─da44d18c-8be3-446e-a5c2-905af545d2c6
# ╠═52cf545a-d7c7-41d8-ad89-617d2f8b3eb9
# ╠═0a6ff75d-fbc7-48a1-924b-e16d2654749c
# ╠═2925bb92-32f8-415d-960f-b0a60a1037b8
