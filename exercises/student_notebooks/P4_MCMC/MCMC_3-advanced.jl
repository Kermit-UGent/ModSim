### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 75581580-2fb2-4112-b397-2b775eb64630
using Pkg; Pkg.activate("..")

# ╔═╡ e07a1ae5-43b7-4c12-831d-43e1738eeac0
using Turing, StatsPlots

# ╔═╡ d12cdd62-f90a-4ba8-8610-5e86e922e881
using PlutoUI

# ╔═╡ f84d9259-69c0-4165-8bc0-d924fef18182
md"# Inference notebook #3: Advanced"

# ╔═╡ aa2b6263-7d13-4683-bc92-25663ed02604
md"## GPS"

# ╔═╡ f3510387-1ab6-4aa2-bed9-9c8297a8b3c5
md"""
GPS systems need to decide what road a car is following based on noisy positional data. We consider here a simplified example. 

At some known timepoints `ts`, we get noisy observations on the car's vertical position `ys_obs` (imagine it as a lattitude of sorts). There are two parallel roads (lines) the car can actually be on, which both have a constant vertical position. If the car is on road 1, then `y = 0`. If it is on road 2, then `y = 1`. The problem is visualized below.
"""

# ╔═╡ 90f185f3-91d4-4ee7-8de3-b76251c0c169
ts = 1:10;

# ╔═╡ 599ae818-fb28-4803-b720-9fcf2bb162b7
ys_obs = [0.6, 0.0, 0.8, -0.7, -0.5, 0.2, 1.0, 1.2, 1.8, 1.1];

# ╔═╡ b2da9dc3-006c-49d0-81a7-6375a2dc6872
begin
	p_cardata = scatter(ts, ys_obs, label = "Observed car positions",
		xlabel = "Time", ylabel = "Vertical position"
	)
	hline!([0.0], color = :orange, label = "Road 1", linewidth = 2)
	hline!([1.0], color = :blue, label = "Road 2", linewidth = 2)
end

# ╔═╡ 46d3214c-e724-4b98-bc1f-5b9913d2b14a
md"""
At some point `t_switch` ∈ [0, 10], the car switches from lane 1 to lane 2. We can describe the model as follows:
- If `t <= t_switch`, then `y ~ Normal(0.0, σ)`,
- If `t > t_switch`, then `y ~ Normal(1.0, σ)`,
with `σ` some (small) noise parameter.
"""

# ╔═╡ 0bdac14c-db71-47b5-95d5-d83fb1e68cab
md"""
Below is a plot showing the car's trajectory for some value of `t_switch`. You can adjust the slider to change this guess value.
"""

# ╔═╡ 9e3db8b6-5466-48b6-9ca9-0d24515b9de3
@bind switchtime Slider(0:0.1:10, default = 5.0, show_value = true)

# ╔═╡ cc0d6b32-6f12-46b8-915f-8a471317c35e
begin
	plot!(
		deepcopy(p_cardata), [0.0, switchtime, switchtime, 10.0],
		[0.0, 0.0, 1.0, 1.0], label = "Car trajectory", color = :black, 
		linewidth = 2, xticks = ([0, 10, switchtime], ["0", "10", "t_switch"])
	)
end

# ╔═╡ 96fc1412-83fc-4f67-8995-065b163d739c
md"""
!!! question
	Infer the posterior probability of `t_switch` given the data.
"""

# ╔═╡ 70fd793a-3ce3-446f-adc1-65ce0a68e48a
@model function cars(ts)
	t_switch ~ missing
	σ ~ missing
	
	for pointidx in 1:length(ts)
		missing
	end
end

# ╔═╡ fb202b86-c058-4f03-9062-ab282c71d5c4
missing # histogram of `t_switch`

# ╔═╡ 34bf815a-3bd3-49b5-b06d-d4297ca213a8
md"## Petridish peril (inference edition)"

# ╔═╡ 78571fb5-c44f-4e7f-afde-4436db6c945b
md"""
We continue with the "petridish peril" question from the previous practical. 

You've made a model to predict bacterial population levels at certain timepoints based on your knowledge of how the species in question grows. You'd now like to update the model with information about the specific strain you're using, so you inoculate a petri dish and count the number of bacteria after a short incubation period.

Incorporate the following information into the model to make it more accurate:
- The population level after 5 hours of incubating was 21000.
- You expect the number of bacteria you count to be Poisson distributed around the actual number.
"""

# ╔═╡ 2a84472c-cb6f-4607-9b98-c88cc2744e3d
md"""
!!! questions
	1. Now taking into account the measurement, what are the chances of your petridish being in a splittable state after 8 hours?
	1. Visualise the updated growth curves.
	1. 🌟(BONUS): The prior for P0 being discrete doesn't allow for the use of a continuous sampler. Change the prior with a sufficiently similar continuous one to fix this. How does this affect the results?
"""

# ╔═╡ 34120240-6fac-42b8-b5c8-bc3271927248
md"""
!!! tip
	Just like in the previous version of the question, `return`ing the estimated logistic function can be useful.
"""

# ╔═╡ 4b94e5a6-1ed6-4095-ab18-06b76d4fec99
logistic(t, P0, r, K) =  K / (1 + (K - P0)/P0 * exp(-r*t))

# ╔═╡ 234b2c87-fe91-4be4-bf0c-a6f20dcc38fe
md"### 1"

# ╔═╡ f080f708-a457-40a3-936c-b82d5159975d
dropletdist = MixtureModel([Poisson(10), Poisson(30)], [0.75, 0.25]);

# ╔═╡ 7b4aef69-10ec-4935-b7fd-4c1d49aa9b3d
@model function petrigrowth()
	P0 ~ dropletdist
    r ~ LogNormal(0.0, 0.3)
	K ~ Normal(1e5, 1e4)

	logfun = t -> logistic(t, P0, r, K)
	Pt = missing
    P_obs ~ missing
	
    return logfun
end

# ╔═╡ 3e98c640-4bb0-4b5d-bae0-769133a599a7
prob_splittable = missing

# ╔═╡ 0ea95b67-d4da-4c5a-ad2e-05024ad074a3
md"### 2"

# ╔═╡ 47b15322-0b8c-48f1-b123-271ebae92655
missing # plot

# ╔═╡ 9ab88be4-4cf8-4747-ac36-3f1b82899be0
md"### 3🌟"

# ╔═╡ 113d8311-7bdc-461c-b077-920e23b33d39
dropletdist🌟 = missing

# ╔═╡ 4e730df9-f619-464a-b8a3-57448132404b
begin
	plot(dropletdist, label = ["Original prior" ""], color = :orange)
	plot!(dropletdist🌟, label = ["Continuous alternative" ""], color = :blue)
		# With mixture models it takes some fiddling to make the labels look nice - don't worry about this, it's not important for the course
end

# ╔═╡ d9b50958-20ae-4085-80d6-19420c7d89df
@model function petrigrowth🌟()
	P0 ~ dropletdist🌟
    r ~ LogNormal(0.0, 0.3)
	K ~ Normal(1e5, 1e4)

	logfun = t -> logistic(t, P0, r, K)
	Pt = missing
    P_obs ~ missing
	
    return logfun
end

# ╔═╡ 0704ed4d-5b3e-401d-a496-98782cb20b09
prob_splittable🌟 = missing

# ╔═╡ Cell order:
# ╟─f84d9259-69c0-4165-8bc0-d924fef18182
# ╠═75581580-2fb2-4112-b397-2b775eb64630
# ╠═e07a1ae5-43b7-4c12-831d-43e1738eeac0
# ╠═d12cdd62-f90a-4ba8-8610-5e86e922e881
# ╟─aa2b6263-7d13-4683-bc92-25663ed02604
# ╟─f3510387-1ab6-4aa2-bed9-9c8297a8b3c5
# ╠═90f185f3-91d4-4ee7-8de3-b76251c0c169
# ╠═599ae818-fb28-4803-b720-9fcf2bb162b7
# ╟─b2da9dc3-006c-49d0-81a7-6375a2dc6872
# ╟─46d3214c-e724-4b98-bc1f-5b9913d2b14a
# ╟─0bdac14c-db71-47b5-95d5-d83fb1e68cab
# ╟─9e3db8b6-5466-48b6-9ca9-0d24515b9de3
# ╟─cc0d6b32-6f12-46b8-915f-8a471317c35e
# ╟─96fc1412-83fc-4f67-8995-065b163d739c
# ╠═70fd793a-3ce3-446f-adc1-65ce0a68e48a
# ╠═fb202b86-c058-4f03-9062-ab282c71d5c4
# ╟─34bf815a-3bd3-49b5-b06d-d4297ca213a8
# ╟─78571fb5-c44f-4e7f-afde-4436db6c945b
# ╟─2a84472c-cb6f-4607-9b98-c88cc2744e3d
# ╟─34120240-6fac-42b8-b5c8-bc3271927248
# ╠═4b94e5a6-1ed6-4095-ab18-06b76d4fec99
# ╟─234b2c87-fe91-4be4-bf0c-a6f20dcc38fe
# ╠═f080f708-a457-40a3-936c-b82d5159975d
# ╠═7b4aef69-10ec-4935-b7fd-4c1d49aa9b3d
# ╠═3e98c640-4bb0-4b5d-bae0-769133a599a7
# ╟─0ea95b67-d4da-4c5a-ad2e-05024ad074a3
# ╠═47b15322-0b8c-48f1-b123-271ebae92655
# ╟─9ab88be4-4cf8-4747-ac36-3f1b82899be0
# ╠═113d8311-7bdc-461c-b077-920e23b33d39
# ╠═4e730df9-f619-464a-b8a3-57448132404b
# ╠═d9b50958-20ae-4085-80d6-19420c7d89df
# ╠═0704ed4d-5b3e-401d-a496-98782cb20b09
