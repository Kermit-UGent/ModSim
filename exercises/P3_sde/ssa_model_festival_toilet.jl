### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ a782b415-6f43-4cf2-a5d4-cff3aa31b187
using Pkg; Pkg.activate("..")

# ╔═╡ e84131d2-3ad2-11f0-02b8-2ffac00498cf
using Catalyst, Plots, DifferentialEquations

# ╔═╡ 0ab79c4e-0206-47e1-b057-3f84e4f66a17
using PlutoUI; TableOfContents()

# ╔═╡ 99c690a3-c1ed-4b41-8427-8b8b04fdcd0f
md"""
# Exercise: Modelling a Festival Toilet
"""

# ╔═╡ c7d127f5-d8b1-4f75-90db-262d70d21482
# https://www.linda.nl/lindanl-assets/uploads/2024/06/03140318/festivaltoilet-1800x1013.png
md"""
![Festival toiletten](https://www.linda.nl/lindanl-assets/uploads/2024/06/03140318/festivaltoilet-1800x1013.png)
"""

# ╔═╡ c8d03b21-b064-42e2-90b1-a7981b17a4a6
md"""
A good festival provides sufficient food and drink. This also means that there will be a lot of discharge. As a student worker, you are responsible for determining the toilet capacity by means of a simulation.

The festival has ten mobile unisex toilets and only one queue. Given that people are discreet and unpredictable, a Jump model seems appropriate.

These are the stocks in your model:
-  $L(t)$: the number of people in the queue (wacht_**L**_ijn) (initially zero);
-  $V(t)$: the number of free (_**V**_rije) toilets;
-  $B(t)$: the number of occupied (_**B**_ezette) toilets (initially zero).

You simulate minute by minute from 0 to 360 minutes (corresponding to 6 p.m. to 12 p.m.). You consider the following flows/processes:
- People are constantly joining the queue (zero-order process), at an average rate of 1.5 per minute.
- When a cubicle is free, one person from the queue enters and the cubicle is occupied. Only one person per cubicle is allowed. This happens quickly, as the kinetics have already been established.
- An occupied cubicle is occupied for an average of five minutes, then it is free again (and the person goes to do something else).

Use the following symbols for the parameters:
- `kₐ`: rate at which people arrive in the queue;
- `kₒ`: rate at which someone from the queue will occupy a free cubicle;
- `kₑ`: rate at which an occupied cubicle becomes available.
"""

# ╔═╡ 15c5c1e3-72d9-4beb-b637-f126400e8da1
md"""
The time interval used is:
"""

# ╔═╡ a30d6812-8337-4361-a16e-ee1bb86b2cb1
tspan = (0.0, 360.0)

# ╔═╡ 2bc030aa-4baa-4e75-90b7-c9fa1e7818ef
md"""
## Part 1: Simulating $L$, $V$ and $B$.
"""

# ╔═╡ f1e9eecb-7068-48fe-a0f2-bed45d41f673
md"""
Complete the model, create the (discrete) jump problem, solve it, and create a graph with the simulation.
"""

# ╔═╡ a9830054-02ee-44f6-a94c-5722959124d7
# toilet = @reaction_network begin
# 	@species L(t)=missing V(t)=missing B(t)=missing
# 	@parameters kₐ=missing kₒ=1/10 kₑ=missing
# 	missing
# 	kₒ, L + V --> B
# 	missing
# end

# ╔═╡ 9540af46-06b1-4c55-a541-47813d5d24d9
missing

# ╔═╡ 89c78b4a-e180-4ec3-a304-1873a39aeb8b
md"""
Create the (discrete) jump problem:
"""

# ╔═╡ a2e9c72a-a4a4-4633-851e-878440c82e36
dprob = missing;

# ╔═╡ b9bddaa3-4293-40a0-8b5b-a277d4fb9f5f
jprob = missing;

# ╔═╡ 08b146ff-fe8a-4780-939a-7917b0a0d08c
md"""
Solve the problem:
"""

# ╔═╡ 26c05ffa-c3c7-4b12-9fb9-7caaceeb9835
sol1 = missing;

# ╔═╡ 67f078fc-c3c2-4e3e-a29a-3200db58e4a3
md"""
Make the graph containing $L$, $V$ and $B$:
"""

# ╔═╡ be53c5e4-ac38-4fda-a19c-c02dd2be5989
missing

# ╔═╡ 0d1124a4-a2d2-4774-add8-4db92142154a
md"""
## Part 2: Maximum Length of the Queue
"""

# ╔═╡ cfaa00e1-4420-46f0-9ee9-03ae32a7a82b
md"""
Calculate the maximum length of the queue (wacht_**L**_ijn). Run your model `1000` times and save the maximum length each time in a vector `line_maxima` (via a `for`-loop). Create a histogram of those maximum lengths.
"""

# ╔═╡ dbd252b3-f8c9-4d9d-906f-463738291b3c
md"""
Calculate the maximum length of the queue $L$ using the solution `sol1` from Part 1:
"""

# ╔═╡ 5f4f76f6-03b2-4224-9969-83ef5be6ad12
missing

# ╔═╡ 2da3f077-1a67-4937-92db-e1baa7642903
md"""
Run your model 12 times and each time save the maximum length (via a `for`-loop) in a vector `line_maxima`.
"""

# ╔═╡ d80e6788-84cb-479f-9e0c-4e9affd86905
# begin
# 	line_maxima = []
# 	for i = 1:missing
# 		sol2 = missing
# 		missing
# 	end
# end

# ╔═╡ 97c72483-8cba-4317-a9ab-2761e153a4ac
# line_maxima

# ╔═╡ 4b0afbed-7325-48ec-b6d4-3bce0f7d5855
md"""
Make the histogram:
"""

# ╔═╡ 6abe6849-5bc0-474c-abbb-4ac45374aab8
missing

# ╔═╡ f3380a9b-7970-4084-b9c9-a6d0a473cf70
md"""
## Part 3: Introducing a Discrete Event
"""

# ╔═╡ 0b8bc24d-8e3d-432c-930e-96cc64a6bae4
md"""
Between 8 p.m. and 9 p.m. is a period that many people use to relieve themselves. During that period, an average of three people arrive per minute (and then returns to normal). Implement this with an event, simulate and create a plot.
"""

# ╔═╡ a808a52c-aef3-4f97-ad97-e6728a1d82db
md"""
Create the conditions, add the condition to the model, and recreate the (discrete) jump problem.
"""

# ╔═╡ 5d83a715-fc3c-46f0-80b8-6125d11403a2
rush_event = missing

# ╔═╡ 9d6ead94-0ac4-4a22-b091-551463151328
# @named toilet_with_rush = missing;

# ╔═╡ aa63d527-f671-4de2-beb6-16736579e18a
dprob3 = missing;

# ╔═╡ b72661ec-3b48-481f-94a2-bdcad402323c
jprob3 = missing;

# ╔═╡ 58ac9de1-8b5b-459b-aebc-77145bdf9901
md"""
Los het probleem op:
"""

# ╔═╡ a0a9c5ce-af7d-470f-98d7-dfb93ecf12e5
sol3 = missing;

# ╔═╡ 6a5126ee-23bc-4642-961f-a6ddb1320505
md"""
Create a graph showing $L$, $V$ and $B$:
"""

# ╔═╡ c8dae895-2315-4593-b975-3f3e24bc21ec
missing

# ╔═╡ a6fe82d9-7627-495b-878b-bc074897def1
md"""
## Part 4: Public Urination Events
"""

# ╔═╡ fccff7b3-8429-4494-a3c8-c6ae692eaac3
md"""
When the queue becomes too long (if there are more than 15 people in the queue), people start urinating in public (on average 1 per minute). These people disappear from the queue. Based on your model from 1), create a new model that simulates this. Create a new variable $W(t)$ (initially zero) that keeps track of the total number of public urination events (_**W**_ildplasevents). Plot the total number of public urination events over time.\
Hint: `>(L, 15)` is 1 when `L` is greater than 15, and zero in other cases. Alternatively `ifelse(L > 15, 1, 0)` is also 1 when `L` is greater than 15, and zero in other cases.
"""

# ╔═╡ 26c6c5e5-a7f8-4265-9423-c69b91b3d9e1
# toilet_wp = @reaction_network begin
# 	@species L(t)=missing V(t)=missing B(t)=missing W(t)=missing
# 	@parameters kₐ=missing kₒ=1/10 kₑ=missing
# 	missing
# 	kₒ, L + V --> B
# 	missing
# 	missing
# end

# ╔═╡ ee055acd-cd86-4e02-bf38-bf5779bd4f95
dprob_wp = missing;

# ╔═╡ b290acfd-52b1-45c4-aaaa-54ee200c5539
jprob_wp = missing;

# ╔═╡ 172782af-f395-4c23-b340-080a3ec90a20
sol_wp = missing;                  # voer enkele keren uit!

# ╔═╡ 0dc015e8-31d9-4e9d-950c-ac94d9ec523f
missing

# ╔═╡ Cell order:
# ╟─99c690a3-c1ed-4b41-8427-8b8b04fdcd0f
# ╟─c7d127f5-d8b1-4f75-90db-262d70d21482
# ╠═a782b415-6f43-4cf2-a5d4-cff3aa31b187
# ╠═e84131d2-3ad2-11f0-02b8-2ffac00498cf
# ╠═0ab79c4e-0206-47e1-b057-3f84e4f66a17
# ╟─c8d03b21-b064-42e2-90b1-a7981b17a4a6
# ╟─15c5c1e3-72d9-4beb-b637-f126400e8da1
# ╠═a30d6812-8337-4361-a16e-ee1bb86b2cb1
# ╟─2bc030aa-4baa-4e75-90b7-c9fa1e7818ef
# ╟─f1e9eecb-7068-48fe-a0f2-bed45d41f673
# ╠═a9830054-02ee-44f6-a94c-5722959124d7
# ╠═9540af46-06b1-4c55-a541-47813d5d24d9
# ╟─89c78b4a-e180-4ec3-a304-1873a39aeb8b
# ╠═a2e9c72a-a4a4-4633-851e-878440c82e36
# ╠═b9bddaa3-4293-40a0-8b5b-a277d4fb9f5f
# ╟─08b146ff-fe8a-4780-939a-7917b0a0d08c
# ╠═26c05ffa-c3c7-4b12-9fb9-7caaceeb9835
# ╟─67f078fc-c3c2-4e3e-a29a-3200db58e4a3
# ╠═be53c5e4-ac38-4fda-a19c-c02dd2be5989
# ╟─0d1124a4-a2d2-4774-add8-4db92142154a
# ╟─cfaa00e1-4420-46f0-9ee9-03ae32a7a82b
# ╟─dbd252b3-f8c9-4d9d-906f-463738291b3c
# ╠═5f4f76f6-03b2-4224-9969-83ef5be6ad12
# ╟─2da3f077-1a67-4937-92db-e1baa7642903
# ╠═d80e6788-84cb-479f-9e0c-4e9affd86905
# ╠═97c72483-8cba-4317-a9ab-2761e153a4ac
# ╟─4b0afbed-7325-48ec-b6d4-3bce0f7d5855
# ╠═6abe6849-5bc0-474c-abbb-4ac45374aab8
# ╟─f3380a9b-7970-4084-b9c9-a6d0a473cf70
# ╟─0b8bc24d-8e3d-432c-930e-96cc64a6bae4
# ╟─a808a52c-aef3-4f97-ad97-e6728a1d82db
# ╠═5d83a715-fc3c-46f0-80b8-6125d11403a2
# ╠═9d6ead94-0ac4-4a22-b091-551463151328
# ╠═aa63d527-f671-4de2-beb6-16736579e18a
# ╠═b72661ec-3b48-481f-94a2-bdcad402323c
# ╟─58ac9de1-8b5b-459b-aebc-77145bdf9901
# ╠═a0a9c5ce-af7d-470f-98d7-dfb93ecf12e5
# ╟─6a5126ee-23bc-4642-961f-a6ddb1320505
# ╠═c8dae895-2315-4593-b975-3f3e24bc21ec
# ╟─a6fe82d9-7627-495b-878b-bc074897def1
# ╟─fccff7b3-8429-4494-a3c8-c6ae692eaac3
# ╠═26c6c5e5-a7f8-4265-9423-c69b91b3d9e1
# ╠═ee055acd-cd86-4e02-bf38-bf5779bd4f95
# ╠═b290acfd-52b1-45c4-aaaa-54ee200c5539
# ╠═172782af-f395-4c23-b340-080a3ec90a20
# ╠═0dc015e8-31d9-4e9d-950c-ac94d9ec523f
