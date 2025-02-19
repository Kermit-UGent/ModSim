### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ ef127ffc-1e2e-4c30-945b-d6cded4d6515
using Pkg; Pkg.activate("..")

# ╔═╡ deb237f9-f0ff-4dfc-8627-43053ccdeb20
using Turing, StatsPlots

# ╔═╡ 2d709e7b-3350-47ec-a3ed-8aa47ad5a8c2
function generate_data(n_wasps = 10; minbound = 0, maxbound = 1000)
	
	x_n, y_n = rand(DiscreteUniform(minbound, maxbound), 2)
	xs, ys = [rand(DiscreteUniform(minbound, maxbound), n_wasps) for _ in 1:2]
	v_wasps = rand(Uniform(5, 10), n_wasps)
	ts = [2*sqrt((x-x_n)^2 + (y-y_n)^2)/v_wasp for (x, y, v_wasp) in zip(xs, ys, v_wasps)]

	return xs, ys, ts, [x_n, y_n]
end;

# ╔═╡ af96af94-d969-4e9a-93f4-e205f8b7f576
md"# Review exercise: Hornet nests"

# ╔═╡ 07f05baa-1336-4e1a-9cbc-5f4506e5b34a
md"""
In recent years, the Asian giant hornet (_Vespa mandarinia_) has become an invasive species in a number of countries, including Belgium. Since they become aggressive when people get close to their nests, the nests often need to be removed when they appear in residential areas. Finding the nests, however, can be a difficult task: the hornets can go hunting over a kilometer from their nest.
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

# ╔═╡ bfb49d79-772d-4066-bd96-14143a1b5eeb
md"""
Consider below the coordinates of marked hornets, as well as their return times.
"""

# ╔═╡ c1850e64-9e2c-46fb-b7cd-22e8af81d3aa
xs, ys, ts, true_location = generate_data();

# ╔═╡ 28cb3363-6856-4164-b60e-36ec2e88ed56
scatter(
	xs, ys, label = "wasp locations", marker_z = ts,
	title = "Locations of wasps colored by return time",
	xlims = (0, 1000), ylims = (0, 1000)
)

# ╔═╡ 484b56ec-57b2-496d-a5cf-1b1c0da97c58
md"""
!!! question
	Where is the hornet nest located? You may assume the nest is somewhere within the plot's boundaries.
"""

# ╔═╡ 66862953-aeb2-4310-90bc-673260f0ecc0
x_nest_sp = missing # vector with possible values of the nest's x-coordinate

# ╔═╡ f522cbcb-112b-4d56-b860-5aee44b1aa2f
y_nest_sp = missing # vector with possible values of the nest's y-coordinate

# ╔═╡ 511f27c3-7cdb-426f-b794-85897ff44135
begin
	scatter(x_nest_sp, y_nest_sp, opacity = 0.1, color = :blue, label = "Estimated nest locations", xlims = (0, 1000), ylims = (0, 1000));
	scatter!(xs, ys, color = :orange, label = "wasp locations", marker_z = ts)
	scatter!(true_location[1:1], true_location[2:2], color = RGB(0, 1, 0), label = "True nest location")
end

# ╔═╡ Cell order:
# ╠═ef127ffc-1e2e-4c30-945b-d6cded4d6515
# ╠═deb237f9-f0ff-4dfc-8627-43053ccdeb20
# ╟─2d709e7b-3350-47ec-a3ed-8aa47ad5a8c2
# ╟─af96af94-d969-4e9a-93f4-e205f8b7f576
# ╟─07f05baa-1336-4e1a-9cbc-5f4506e5b34a
# ╟─90aee617-84d5-4915-afb4-e7d422cfb4b7
# ╟─082a4cc0-b151-4f8b-87da-278484218e6e
# ╟─bfb49d79-772d-4066-bd96-14143a1b5eeb
# ╠═c1850e64-9e2c-46fb-b7cd-22e8af81d3aa
# ╟─28cb3363-6856-4164-b60e-36ec2e88ed56
# ╟─484b56ec-57b2-496d-a5cf-1b1c0da97c58
# ╠═66862953-aeb2-4310-90bc-673260f0ecc0
# ╠═f522cbcb-112b-4d56-b860-5aee44b1aa2f
# ╠═511f27c3-7cdb-426f-b794-85897ff44135
