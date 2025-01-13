### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ ef127ffc-1e2e-4c30-945b-d6cded4d6515
using Pkg; Pkg.activate("..")

# ╔═╡ deb237f9-f0ff-4dfc-8627-43053ccdeb20
using Turing, StatsPlots

# ╔═╡ 2d709e7b-3350-47ec-a3ed-8aa47ad5a8c2
function generate_data(n_wasps = 10;
	minbound = 0, maxbound = 1000, spoil_answer = false
	)
	
	x_n, y_n = rand(DiscreteUniform(minbound, maxbound), 2)
	xs, ys = [rand(DiscreteUniform(minbound, maxbound), n_wasps) for _ in 1:2]
	v_wasps = rand(Uniform(5, 10), n_wasps)
	ts = [2*sqrt((x-x_n)^2 + (y-y_n)^2)/v_wasp for (x, y, v_wasp) in zip(xs, ys, v_wasps)]
	spoil_answer && println("The true location of the nest is: $x_n, $y_n")

	return xs, ys, ts
end;

# ╔═╡ af96af94-d969-4e9a-93f4-e205f8b7f576
md"# Review exercise: Hornet nests"

# ╔═╡ 07f05baa-1336-4e1a-9cbc-5f4506e5b34a
md"""
In recent years, the Asian giant hornet (_Vespa mandarinia_) has become an invasive species in a number of countries, including Belgium. Since they become aggressive when people get close to their nests, the nests often need to be removed when they appear in residential areas. Finding the nests, however, can be a difficult task: the hornets can go hunting over a km from their nest.
"""

# ╔═╡ 90aee617-84d5-4915-afb4-e7d422cfb4b7
md"""
One method for finding the nest is to set up a feeder station, mark any hornets gathering food, and record how long it takes for them to fly back to their nest with it and return for more. Making an estimate of their flight speed, the return time can be used to infer the distance of that location to the nest. Repeated measurements in other locations gives enough information for a triangulation of sorts.
"""

# ╔═╡ 082a4cc0-b151-4f8b-87da-278484218e6e
md"""
![The Asian giant hornet](https://upload.wikimedia.org/wikipedia/commons/thumb/1/19/Vespa_mandarinia_japonica1.jpg/1280px-Vespa_mandarinia_japonica1.jpg)
*The Asian giant hornet (credit: Picture by KENPEI on Wikipedia)*
"""

# ╔═╡ 484b56ec-57b2-496d-a5cf-1b1c0da97c58
md"""
!!! question
	Consider below the coordinates of marked hornets, as well as their return times. Where is the hornet nest located? You may assume the nest is somewhere within the plot's boundaries.
"""

# ╔═╡ c1850e64-9e2c-46fb-b7cd-22e8af81d3aa
xs, ys, ts = generate_data(spoil_answer = true);

# ╔═╡ 28cb3363-6856-4164-b60e-36ec2e88ed56
scatter(xs, ys, label = "wasp locations", marker_z = ts, title = "Locations of wasps colored by return time", xlims = (0, 1000), ylims = (0, 1000))

# ╔═╡ ac0496d3-aa62-4c07-8233-16fe67c52b24
@model function horenaars(xs, ys, ts)
    x_nest ~ Uniform(0, 1000)
    y_nest ~ Uniform(0, 1000)
    v_wasp ~ Gamma(8)

    dists = @. sqrt((xs - x_nest)^2 + (ys - y_nest)^2)
    for i in eachindex(ts)
        ts[i] ~ Normal(2*dists[i] / v_wasp, 10)
    end
end

# ╔═╡ 88e95234-4e43-4ead-bd0f-f32f9fc109c2
n_samples = 1000

# ╔═╡ c5a94b3e-a4b8-42e6-a8c7-79e942212852
chain = sample(horenaars(xs, ys, ts), NUTS(), n_samples)

# ╔═╡ 1c84ab96-b3db-494e-860e-80e2d178da43
plot(chain)

# ╔═╡ 0c6be1fa-bff5-495b-90e0-90dabe03105a
begin
	scatter(chain[:x_nest], chain[:y_nest], opacity = 0.1, color = :blue, label = "Estimated nest locations", xlims = (0, 1000), ylims = (0, 1000));
	scatter!(xs, ys, color = :orange, label = "wasp locations", marker_z = ts)
end

# ╔═╡ Cell order:
# ╠═ef127ffc-1e2e-4c30-945b-d6cded4d6515
# ╠═deb237f9-f0ff-4dfc-8627-43053ccdeb20
# ╟─2d709e7b-3350-47ec-a3ed-8aa47ad5a8c2
# ╟─af96af94-d969-4e9a-93f4-e205f8b7f576
# ╟─07f05baa-1336-4e1a-9cbc-5f4506e5b34a
# ╟─90aee617-84d5-4915-afb4-e7d422cfb4b7
# ╟─082a4cc0-b151-4f8b-87da-278484218e6e
# ╟─484b56ec-57b2-496d-a5cf-1b1c0da97c58
# ╠═c1850e64-9e2c-46fb-b7cd-22e8af81d3aa
# ╟─28cb3363-6856-4164-b60e-36ec2e88ed56
# ╠═ac0496d3-aa62-4c07-8233-16fe67c52b24
# ╠═88e95234-4e43-4ead-bd0f-f32f9fc109c2
# ╠═c5a94b3e-a4b8-42e6-a8c7-79e942212852
# ╠═1c84ab96-b3db-494e-860e-80e2d178da43
# ╟─0c6be1fa-bff5-495b-90e0-90dabe03105a
