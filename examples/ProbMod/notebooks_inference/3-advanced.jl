### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 75581580-2fb2-4112-b397-2b775eb64630
using Pkg; Pkg.activate("..")

# ╔═╡ e07a1ae5-43b7-4c12-831d-43e1738eeac0
using Turing, StatsPlots

# ╔═╡ f84d9259-69c0-4165-8bc0-d924fef18182
md"# Inference notebook #3: Advanced"

# ╔═╡ a0ba1731-b915-4f50-a600-34c0c47c7e2e
n_samples = 1000

# ╔═╡ aa2b6263-7d13-4683-bc92-25663ed02604
md"## Fast and furious cars"

# ╔═╡ f3510387-1ab6-4aa2-bed9-9c8297a8b3c5
md"""
The CEO of an electrical car company has put you in charge of programming the new **self-driving cars**. To start off, we simply want to use the car's coordinates to predict its trajectory and see if it's going to coincide with a wall.

Make the following assumptions:
- The car moves in a straight line.
- The slope of the car's trajectory is 0.0 ± 1.0.
- There is noise on the car's y-coordinates. The deviation between the car's observed y-coordinates and real y-coordinates is expected to be 0.5
- For simplicity, assume the x-coordinates are exactly know.
"""

# ╔═╡ 103de21d-9d21-453f-842f-edc8fda0a274
md"""
!!! question
	Given the data below, what is the probability the car will crash into the wall?
"""

# ╔═╡ 8d8c91bd-71ff-4010-b8cf-2af0eb9034b6
xs = [0, 1, 2, 3, 4, 5, 6, 7, 8] # assumed to be known exactly

# ╔═╡ 58f173d2-3659-417e-b943-914228fe88d5
ys_obs = [5.8, 6.4, 4.4, 3.9, 4.5, 3.3, 4.7, 5.6, 5.4] # noisy

# ╔═╡ 9e87ebb9-8229-4d93-a081-59ba1e49c55d
begin
	scatter(xs, ys_obs, xlims = (0, 10), ylims = (0, 10), label = "Estimated car positions")
	xs_wall = [10.0, 10.0]
    ys_wall = [-100.0, 5.0]
    plot!(xs_wall, ys_wall, label = "THE WALL", linewidth = 10, color = :red);
end

# ╔═╡ b629b597-840b-46ca-b98c-867b164922ef
@model function boink(xs)
    y0 ~ Uniform(0, 10)
    slope ~ Normal(0, 1)

	ys_obs = zeros(length(xs))
	for i in eachindex(xs)
		y = slope * xs[i] + y0
        ys_obs[i] ~ Normal(y, 0.5)
    end
end

# ╔═╡ 7303f8cb-2a22-4c94-8767-e5816e9a76b5
crash_model = boink(xs) | (ys_obs = ys_obs,)

# ╔═╡ 1c59425a-29e0-4b7b-9691-781ea0f40a95
crash_chain = sample(crash_model, NUTS(), 2000)

# ╔═╡ 98e837ce-b0f6-459b-b30c-f3201369244d
plot(crash_chain)

# ╔═╡ 2dec98ac-bb7d-4260-8d0a-38dbdcc86c00
predicted_lines = [
	x -> x*crash_chain[:slope][i] + crash_chain[:y0][i]
	for i in eachindex(crash_chain[:y0])
]

# ╔═╡ 10a262d3-c335-4067-8ee9-072cc62ca9f6
begin
	scatter(xs, ys_obs, xlims = (0, 10), ylims = (0, 10))
	plot!(predicted_lines, color = :blue, alpha = 0.03, label = false)
end

# ╔═╡ 1ab2f14e-13ea-4fbd-8a0a-eb033f47cc1a
sp_crash = [pline(10) <= 5 for pline in predicted_lines]

# ╔═╡ 9a2f9d11-7962-4167-a271-ecdd67c775af
mean(sp_crash)

# ╔═╡ 0b5b0684-0c71-48cb-aef4-658d5e20c8fa
md"There's about an 88% chance the car is heading towards the wall!"

# ╔═╡ 79145c9d-c6e3-4a0c-88ba-1ee7c37d2b6f
md"### Tokyo drift"

# ╔═╡ cec0c2a9-bba7-43b1-b1a1-6cf0528d70db
md"""
As you share your first results, someone notifies you that, actually, the car **has already learned how to drift**. This means that the car's trajectory may actually not be just a single line, but rather two (connected) lines.

Extend your previous model so that the trajectory of the car consists of **two** straight lines, where the car follows line 1 before the lines intersect, and line 2 afterwards.
"""

# ╔═╡ 6f388078-0918-4667-82e2-ffbed70350cc
md"""
!!! questions
	- Does the new model fit the data better?
	- What is the chance of wall-collision according to the new model?
"""

# ╔═╡ 89002e16-e1b8-403a-ae9c-446f7e9c9957
@model function drifting(xs)
    y0 ~ Uniform(0, 10)
    slope ~ Normal(0, 1)
	
	x_drift ~ Uniform(4, 10)
	y_dr = y0 + x_drift*slope
	slope_dr ~ Normal(0, 1)
	y0_dr ~ Normal(y_dr - (slope_dr*x_drift), 0.1)
	
	ys_obs = zeros(length(xs))
	for i in eachindex(xs)
		if xs[i] <= x_drift
            y = y0 + slope*xs[i]
        else
            y = y0_dr + slope_dr*xs[i]
        end
        ys_obs[i] ~ Normal(y, 0.5)
    end
end

# ╔═╡ 1c0f631c-0679-46d5-a802-7b10390a1d68
drift_model = drifting(xs) | (ys_obs = ys_obs,)

# ╔═╡ e1ab69ac-de9d-4a20-a2ba-48a57eaf3451
drift_chain = sample(drift_model, NUTS(), 2000)

# ╔═╡ a0d56497-ea1c-4d8c-ab30-db9702064ab9
plot(drift_chain)

# ╔═╡ e7a95eb3-b2b8-42dc-8602-8ec82d70434d
predrift_lines = [
	x -> x * drift_chain[:slope][i] + drift_chain[:y0][i]
	for i in eachindex(drift_chain[:y0])
];

# ╔═╡ fca024ed-7b43-466d-a3d2-63faf168eeae
drifting_lines = [
	x -> x * drift_chain[:slope_dr][i] + drift_chain[:y0_dr][i]
	for i in eachindex(drift_chain[:y0_dr])
];

# ╔═╡ 406d6c5d-4973-411c-9273-9e8f92f92b88
begin
	scatter(xs, ys_obs, xlims = (0, 10), ylims = (0, 10))
	plot!(predrift_lines, color = :blue, alpha = 0.01, label = false)
	plot!(drifting_lines, color = :blue, alpha = 0.01, label = false)
end

# ╔═╡ 8700c2cd-d241-4de9-9941-5e74efbcbd6d
md"Comparing visually, the new predicted trajectory lines up much better with the data than the old trajectory. It seems safe to assume the car actually did turn."

# ╔═╡ 290b3ab5-909f-4153-a7b6-6274413b515c
sp_driftingcrash = [pline(10) <= 5 for pline in drifting_lines]

# ╔═╡ 8ee83453-46ed-4244-b838-b1c808554f45
mean(sp_driftingcrash)

# ╔═╡ 779309c5-c2fb-40c6-87f1-dfede49cc9eb
md"Under our new model, the car crashing seems very unlikely."

# ╔═╡ 34bf815a-3bd3-49b5-b06d-d4297ca213a8
md"## Petridish peril (inference edition)"

# ╔═╡ 78571fb5-c44f-4e7f-afde-4436db6c945b
md"""
We continue with the "petridish peril" question from the previous practical. 

You've made a model to predict bacterial population levels at certain timepoints based on your knowledge of how the species in question grows. You'd now like to update the model with information about the specific strain you're using, so you inoculate a few petri dishes and count them after a short incubation period.

Incorporate the following information into the model to make it more accurate:
- The population levels after 5 hours of incubating were 6000, 21000 and 4000.
- You expect the amount of bacteria you count to be normally distributed around the actual number with a standard deviation of 1000.
"""

# ╔═╡ 2a84472c-cb6f-4607-9b98-c88cc2744e3d
md"""
!!! questions
	- Now taking into account the measurements, what are the chances of your petridish being in a splittable state after 7, 8 and 9 hours?
	- Visualise the updated growth curves.
"""

# ╔═╡ 4b94e5a6-1ed6-4095-ab18-06b76d4fec99
logistic(t, P0, r, K) =  K / (1 + (K - P0)/P0 * exp(-r*t))

# ╔═╡ f080f708-a457-40a3-936c-b82d5159975d
dropletdist = MixtureModel([Uniform(1, 30), Uniform(20, 50)], [0.75, 0.25]);

# ╔═╡ 7b4aef69-10ec-4935-b7fd-4c1d49aa9b3d
@model function petrigrowth(t)
	P0 ~ dropletdist
    r ~ Normal(1.0, 0.2)
	K ~ LogNormal(log(1e5), 0.3)

	Pt = logistic(t, P0, r, K)
    P_obs ~ Normal(Pt, 1000)
    return Pt
end

# ╔═╡ d0cdebbe-3339-4ebc-9aed-cfdfbd8da3c3
P_obs = [6_000, 21_000, 4_000]

# ╔═╡ a5b5a65e-2bac-49fa-b4ca-2fcaf64e4ada
petrimodel = petrigrowth(5) | (P_obs = P_obs,)

# ╔═╡ b579ccc6-993e-451b-a137-6fff6b630b49
petri_chain = sample(petrimodel, NUTS(), n_samples)

# ╔═╡ 9633f07a-d583-4680-b3cc-f5701540968f
plot(petri_chain)

# ╔═╡ 7d0e5f0b-2cdf-4946-a186-f70774e363bb
P0s, rs, Ks = petri_chain[:P0], petri_chain[:r], petri_chain[:K]

# ╔═╡ 21895fd3-0b4a-48ad-9e9e-04fbac17bcbb
for t in [7, 8, 9]
	sp_petri = [logistic(t, P0s[i], rs[i], Ks[i]) for i in eachindex(P0s)]
	prob_splittable = mean((sp_petri .>= 1e4) .&& (sp_petri .<= 1e5))
	println("After $t hours the probability is ", prob_splittable)
end

# ╔═╡ 7445b378-8ee3-4646-9d93-3283c02c20b7
logistfuns = [
	t -> logistic(t, P0, r, K) 
	for (P0, r, K) in zip(petri_chain[:P0], petri_chain[:r], petri_chain[:K])
];

# ╔═╡ 3d9b73ef-2aae-4c19-bddc-4081927ec92d
plot(logistfuns[1:100], xlims = (0, 12), legend = false, color = :skyblue, alpha = 0.5)

# ╔═╡ Cell order:
# ╟─f84d9259-69c0-4165-8bc0-d924fef18182
# ╠═75581580-2fb2-4112-b397-2b775eb64630
# ╠═e07a1ae5-43b7-4c12-831d-43e1738eeac0
# ╠═a0ba1731-b915-4f50-a600-34c0c47c7e2e
# ╟─aa2b6263-7d13-4683-bc92-25663ed02604
# ╟─f3510387-1ab6-4aa2-bed9-9c8297a8b3c5
# ╟─103de21d-9d21-453f-842f-edc8fda0a274
# ╠═8d8c91bd-71ff-4010-b8cf-2af0eb9034b6
# ╠═58f173d2-3659-417e-b943-914228fe88d5
# ╟─9e87ebb9-8229-4d93-a081-59ba1e49c55d
# ╠═b629b597-840b-46ca-b98c-867b164922ef
# ╠═7303f8cb-2a22-4c94-8767-e5816e9a76b5
# ╠═1c59425a-29e0-4b7b-9691-781ea0f40a95
# ╠═98e837ce-b0f6-459b-b30c-f3201369244d
# ╠═2dec98ac-bb7d-4260-8d0a-38dbdcc86c00
# ╠═10a262d3-c335-4067-8ee9-072cc62ca9f6
# ╠═1ab2f14e-13ea-4fbd-8a0a-eb033f47cc1a
# ╠═9a2f9d11-7962-4167-a271-ecdd67c775af
# ╟─0b5b0684-0c71-48cb-aef4-658d5e20c8fa
# ╟─79145c9d-c6e3-4a0c-88ba-1ee7c37d2b6f
# ╟─cec0c2a9-bba7-43b1-b1a1-6cf0528d70db
# ╟─6f388078-0918-4667-82e2-ffbed70350cc
# ╠═89002e16-e1b8-403a-ae9c-446f7e9c9957
# ╠═1c0f631c-0679-46d5-a802-7b10390a1d68
# ╠═e1ab69ac-de9d-4a20-a2ba-48a57eaf3451
# ╠═a0d56497-ea1c-4d8c-ab30-db9702064ab9
# ╠═e7a95eb3-b2b8-42dc-8602-8ec82d70434d
# ╠═fca024ed-7b43-466d-a3d2-63faf168eeae
# ╟─406d6c5d-4973-411c-9273-9e8f92f92b88
# ╟─8700c2cd-d241-4de9-9941-5e74efbcbd6d
# ╠═290b3ab5-909f-4153-a7b6-6274413b515c
# ╠═8ee83453-46ed-4244-b838-b1c808554f45
# ╟─779309c5-c2fb-40c6-87f1-dfede49cc9eb
# ╟─34bf815a-3bd3-49b5-b06d-d4297ca213a8
# ╟─78571fb5-c44f-4e7f-afde-4436db6c945b
# ╟─2a84472c-cb6f-4607-9b98-c88cc2744e3d
# ╠═4b94e5a6-1ed6-4095-ab18-06b76d4fec99
# ╠═f080f708-a457-40a3-936c-b82d5159975d
# ╠═7b4aef69-10ec-4935-b7fd-4c1d49aa9b3d
# ╠═d0cdebbe-3339-4ebc-9aed-cfdfbd8da3c3
# ╠═a5b5a65e-2bac-49fa-b4ca-2fcaf64e4ada
# ╠═b579ccc6-993e-451b-a137-6fff6b630b49
# ╠═9633f07a-d583-4680-b3cc-f5701540968f
# ╠═7d0e5f0b-2cdf-4946-a186-f70774e363bb
# ╠═21895fd3-0b4a-48ad-9e9e-04fbac17bcbb
# ╠═7445b378-8ee3-4646-9d93-3283c02c20b7
# ╠═3d9b73ef-2aae-4c19-bddc-4081927ec92d
