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

# ╔═╡ 309035dd-5653-48a6-a53d-817e743279fa
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 6b342f14-e7d5-11ef-1ea0-77ceb0d78f32
using Markdown

# ╔═╡ 7bfcc024-7d7f-4c0b-b918-8f7626b10974
using InteractiveUtils

# ╔═╡ 71140c81-af29-4857-8020-4f94c8bd64b3
using Catalyst, DifferentialEquations, Plots, Distributions

# ╔═╡ 284f5847-9c15-41f3-a595-1e12a22df69f
using PlutoUI; TableOfContents()

# ╔═╡ 1f975552-b0b8-4830-8dcc-214574d4fc38
md"""
# Exercise: Modeling a simple Bike Share System
"""

# ╔═╡ d2f32eab-0b35-4794-9219-5bcbb4c069c5
md"""
Imagine a bike share system for students traveling between Olin College and Wellesley College, which are about three miles apart in eastern Massachusetts. Suppose the system contains 12 bikes and two bike racks, one at Olin and one at Wellesley, each with the capacity to hold 12 bikes. As students arrive, check out a bike, and ride to the other campus, the number of bikes in each location changes.
Initially there are 10 bikes at Olin and, hence, 2 bikes at Wellesley. For this simple model, we will aslo assume that the changes in the number of bikes at both locations is instantenuously. The rate at which a bike is moved from Olin to Wellesley is denoted as $p_1$ ($min^{-1}$); the rate at which a bike is moved from Wellesley to Olin is denoted as $p_2$ ($min^{-1}$). Both processes are zeroth-order and we want to see the evolution of bikes during $1\,h = 60\,min$.

This is a discreet and stochastic problem and you need to solve it with SSA.
"""

# ╔═╡ 016842c9-9479-4061-a27e-9dc006121f23
md"""
Create a *reaction network object* model for the aforementioned problem in order to simulate the evolution of the number of bikes at Olin ($O$) and Wellesley ($W$) with time. Name it `bike_sharing`.
"""

# ╔═╡ 6c97bf81-ef32-45a4-aa7c-c8c26ba2d2c3
# bike_sharing = @reaction_network begin
# 	missing
# 	...
# end
bike_sharing = @reaction_network begin
    @species O(t)=10 W(t)=2
	p₁, O => W
	p₂, W => O
end

# ╔═╡ 9c7ab7fb-7380-41a3-85ea-714478ade218
md"""
Convert the system to a symbolic differential equation model and verify, by analyzing the differential equation, that your model is correctly implemented.
"""

# ╔═╡ 1536fe23-0f8d-4b86-98d2-076248b35954
# osys = missing
osys = convert(ODESystem, bike_sharing)

# ╔═╡ e6e2ff5c-38eb-4ba3-b430-c9031483a0a5
md"""
Initialize a vector `u0` with the initial conditions:
"""

# ╔═╡ ab6af765-1cde-4da8-bbc1-a5fab391db54
# u0 = missing
u0 = [:O => 10, :W => 2]

# ╔═╡ 378878a0-5c09-4eb0-ac43-1031014ff12a
md"""
Set the timespan for the simulation:
"""

# ╔═╡ 3ae98e83-7beb-4597-89be-80c813d4349b
tspan = (0.0, 60.0)

# ╔═╡ 988f79c0-9c7b-4752-a7f2-d4473ad73ce6
md"""
Create a slider for the variable `p₁` in the range of $0.0$ and $1.0$ with a step of $0.1$. Take a default value of $0.0$.
"""

# ╔═╡ 0d8f53f8-0a14-4ac6-bd0c-2190d4db0909
# @bind p₁ missing
@bind p₁ Slider(0.0:0.1:1, default=0.0, show_value=true)

# ╔═╡ 08d43ac8-a973-4d7b-baf7-4c37e54cfe24
md"
Initialize vector `parms` with parameter values, `p₁` is the slider value and assign a constant value of `0.3` to `p₂`.
"

# ╔═╡ e20e4dd8-bdbb-4005-af68-6bf7e4ec130e
# parms = missing
parms = [:p₁=>p₁, :p₂=>0.30]

# ╔═╡ 238e1120-34af-4d57-8efa-aa80ab28a874
md"""
Create a DiscreteProblem and store it in `dprob`:
"""

# ╔═╡ d4c45709-70c9-4ba0-8fb8-6b600473723d
# dprob = missing
dprob = DiscreteProblem(bike_sharing, u0, tspan, parms)

# ╔═╡ d06fb076-76e4-4248-a940-96804ea68833
md"""
Create a JumpProblem and store it in `jdprob`. Use the simulation method `Direct()` and an additional option `save_positions=(false, false)`. The latter prohibits to save the values just before and after the jump event (later, when solving the problem we will namely use `saveat=1.0`).
"""

# ╔═╡ 7644adf4-d992-48b1-b40a-12fdf30f6cb5
jdprob = JumpProblem(bike_sharing, dprob, Direct(), save_positions=(false, false));
# https://docs.sciml.ai/JumpProcesses/dev/jump_solve/#JumpProcesses.jl

# ╔═╡ 59d2d3e1-354b-4444-b8a1-16ad8ea2ba94
md"""
If we would now solve the problem, you might for example encounter negative values for the number of bikes at Olin and, hence, a value larger than 12 at  Wellesley. Of course, this makes no sense. Therefore, we will create the so-called `condition` function that needs to invoke the so-called `affect!` function at each jump event in order to check on $O$ and $W$ and setting them to valid values if necessary.
"""

# ╔═╡ 53e15d08-7d1d-4bc2-9fd3-4bc6fbb9de84
md"""
Create the `condition` function.
"""

# ╔═╡ 1511d269-9706-41b1-b8e2-f85ebcedc2d8
# function condition(u, t, integrator)
# 	missing
# end
function condition(u, t, integrator)
	true
end

# ╔═╡ fd9b4521-30d2-4af3-a080-40e9dbf27aa4
md"""
Create the `affect` function.

Hints:
- If the number of bikes at Olin is larger than 12 (this implies that the number of bikes at Wellesley is negative), then the number of bikes at Olin should be set to 12 and the number of bikes at Wellesley should be set to 0.
- If the number of bikes at Olin is negative (this implies that the number of bikes at Wellesley is larger than 12), then the number of bikes at Olin should be set to 0 and the number of bikes at Wellesley should be set to 12.
- The number of bikes at Olin is accessed through `integrator.u[1]` and the number of bikes at Wellesley is accessed through `integrator.u[2]`.
"""

# ╔═╡ 3cac94d6-14ae-42f2-b0b6-f47d87cdb518
# function affect!(integrator)
# 	if missing
# 		missing
# 		missing
# 	end
# 	if missing
# 		missing
# 		missing
# 	end
# end
function affect!(integrator)
	if integrator.u[1] > 12
		integrator.u[1] = 12
		integrator.u[2] = 0
	end
	if integrator.u[1] < 0
		integrator.u[1] = 0
		integrator.u[2] = 12
	end
end

# ╔═╡ b2797fd2-b6e1-4bb0-af23-26931fe8be69
md"""
Create the discrete callback function. Store it in `cb`. Again use the option `save_positions=(false,false)`.
"""

# ╔═╡ f947c2d9-9123-422d-8972-157717c85b3c
# cb = missing
cb = DiscreteCallback(condition, affect!, save_positions=(false,false));

# ╔═╡ 74708270-b1ec-48c7-af32-3b970b92c706
md"""
Solve the problem and store it in `jdsol`. Use the `SSAStepper()` stepping algorithm, the option `saveat=1.0` and the callback function.
"""

# ╔═╡ 2b00df5d-994e-47a1-8068-c93ce3f1a618
jdsol = solve(jdprob, SSAStepper(), saveat=1.0, callback=cb);

# ╔═╡ 9d06c31e-3525-4889-a1de-3fe02413c7d8
md"""
Plot the solution.
"""

# ╔═╡ 9a90f800-3669-4831-b50b-c5405bbb9a03
plot(jdsol, ylim=(0, 12))

# ╔═╡ a554fd16-aa3d-48ca-8de6-5582725c27d8
md"""
Analyse the results. See what happens when you:
- run the notebook cell with the `solve` function repeatly
- change the value of `p₁` using the slider
"""

# ╔═╡ d6872046-b5ef-4c2d-a9bb-2418f57f715d
md"""
!!! question "Question"
	From what value of `p₁` do you start to get empty bike racks at Olin?
"""

# ╔═╡ 747e20c4-b06b-4e78-a09a-55053cf42bf4
md"- Answer: missing"

# ╔═╡ 92181028-60fc-4830-afba-2380ac91455d
md"""
You can inspect the actual number of bike values at Olin by using `jdsol[:O]`:
"""

# ╔═╡ f8942b10-773a-4b22-baad-8004fba8bd34
# missing
jdsol[:O]

# ╔═╡ 43d41284-053d-4dfe-8d5b-96be70c0495c
md"""
If you want to have a `true` boolean value on positions where the vector value is zero (and `false` on non-zero values), then you would compare `jdsol[:O]` element wise with `0`. In Julia, if you want to do element wise operations with/on vectors, you always need to place a dot (`.`) in front of the operator, like for example `.==`.

Compare in that way `jdsol[:O]` with `0`:
"""

# ╔═╡ a999ae2a-7567-41e7-9c0c-e94fad6f5d46
# missing
jdsol[:O] .== 0

# ╔═╡ 9ebb5b44-04d7-4b89-acdb-e40a245703d2
md"""
Furthermore, if you want the count the number of `true` values in the latter (hence, the zero element values), you can simply use the function `count(...)`. Count the number of zeros:
"""

# ╔═╡ 049de8d5-b221-452b-b2c4-9bc1e0c17f48
# missing
count(jdsol[:O] .== 0)

# ╔═╡ a73a2853-1f48-4179-9771-083794d3f137
md"""
Using the aforementioned way to count zeros in a vector, we will now count the zeros for a range of $p$ values. Because of the stochastic behaviour of the system, for each $p$ values we will count the zeros for a $1000$ simulations and then storing only the average value.

To introduce a new value for $p_1$ you can take a deepcopy of the problem and remake the problem like this:
- `jdprob_re = remake(deepcopy(jdprob); p=[:p₁=>p_val])`
and then solving the problem and store it in `jdsol_re`.

In the layout below, `mean_zero_counts` while contain the final mean values of the averaged numbers of zeros from a 1000 simulations using a specific $p$ value, `zero_counts_p_val` will contain the actual number of zeros for a 1000 simulations using a specific $p$ value.

Use the layout below to fill in `mean_zero_counts`.
"""

# ╔═╡ 682e9120-0e1c-4dfa-9ec6-66bb0a3f4374
# begin
# 	p_values = 0.0:0.1:1.0  # different p-values
# 	mean_zero_counts = []   # vector to store the corresponding mean zero values
# 	for p_val = p_values    # p_val will be each of the p_values
# 		zero_counts_p_val = []  # vector to store the zeros for the 1000 simulations
# 		for i = 1:1000          # do a 1000 simulations
# 			# take a deepcopy and remake the problem for the specific p-value
# 			jdprob_re = missing
# 			# solve the problem
# 			jdsol_re = missing
# 			# append the number of zeros to zero_counts_p_val
# 			append!(..., ...)
# 		end
# 		# append the mean number of zeros to mean_zero_counts
# 		append!(..., ...)
# 	end
# end
begin
	p_values = 0.0:0.1:1.0  # different p-values
	mean_zero_counts = []   # vector to store the corresponding mean zero values
	for p_val = p_values    # p_val will be each of the p_values
		zero_counts_p_val = []  # vector to store the zeros for the 1000 simulations
		for i = 1:1000          # do a 1000 simulation
			# take a deepcopy and remake the problem for the specific p-value
			jdprob_re = remake(deepcopy(jdprob); p=[:p₁=>p_val])
			# solve the problem
			jdsol_re = solve(jdprob_re, SSAStepper(), saveat=1.0, callback=cb);
			# append the number of zeros to zero_counts_p_val
			append!(zero_counts_p_val, count(jdsol_re[:O] .== 0))
		end
		# append the mean number of zeros to mean_zero_counts
		append!(mean_zero_counts, mean(zero_counts_p_val))
	end
end

# ╔═╡ 705d3fcb-20b6-4481-a304-1d3ccd623674
md"""
Have a look at the mean zero counts by typing `mean_zero_counts`:
"""

# ╔═╡ 5968317a-6c07-4655-8137-6702656bb3b4
# missing
mean_zero_counts

# ╔═╡ ff9370d8-3395-4382-9f51-afa11748319e
md"""
Plot the mean zero counts as a function of the $p$-values.
"""

# ╔═╡ 48be49d0-0b60-44f3-8152-1ca917a4232e
plot(p_values, mean_zero_counts)

# ╔═╡ d6452915-bdf0-48f0-8c7d-3df83c7bce72
md"""
!!! question "Question"
	- From what value of $p$ do the empty number of bike racks at Olin clearly begin to rise?
	- Reflect on this, does this make sense? Hint: change the value of $p_2$ and observe what happens.
"""

# ╔═╡ Cell order:
# ╠═6b342f14-e7d5-11ef-1ea0-77ceb0d78f32
# ╠═7bfcc024-7d7f-4c0b-b918-8f7626b10974
# ╠═309035dd-5653-48a6-a53d-817e743279fa
# ╠═71140c81-af29-4857-8020-4f94c8bd64b3
# ╠═284f5847-9c15-41f3-a595-1e12a22df69f
# ╟─1f975552-b0b8-4830-8dcc-214574d4fc38
# ╟─d2f32eab-0b35-4794-9219-5bcbb4c069c5
# ╟─016842c9-9479-4061-a27e-9dc006121f23
# ╠═6c97bf81-ef32-45a4-aa7c-c8c26ba2d2c3
# ╟─9c7ab7fb-7380-41a3-85ea-714478ade218
# ╠═1536fe23-0f8d-4b86-98d2-076248b35954
# ╟─e6e2ff5c-38eb-4ba3-b430-c9031483a0a5
# ╠═ab6af765-1cde-4da8-bbc1-a5fab391db54
# ╟─378878a0-5c09-4eb0-ac43-1031014ff12a
# ╠═3ae98e83-7beb-4597-89be-80c813d4349b
# ╟─988f79c0-9c7b-4752-a7f2-d4473ad73ce6
# ╠═0d8f53f8-0a14-4ac6-bd0c-2190d4db0909
# ╟─08d43ac8-a973-4d7b-baf7-4c37e54cfe24
# ╠═e20e4dd8-bdbb-4005-af68-6bf7e4ec130e
# ╟─238e1120-34af-4d57-8efa-aa80ab28a874
# ╠═d4c45709-70c9-4ba0-8fb8-6b600473723d
# ╟─d06fb076-76e4-4248-a940-96804ea68833
# ╠═7644adf4-d992-48b1-b40a-12fdf30f6cb5
# ╟─59d2d3e1-354b-4444-b8a1-16ad8ea2ba94
# ╟─53e15d08-7d1d-4bc2-9fd3-4bc6fbb9de84
# ╠═1511d269-9706-41b1-b8e2-f85ebcedc2d8
# ╟─fd9b4521-30d2-4af3-a080-40e9dbf27aa4
# ╠═3cac94d6-14ae-42f2-b0b6-f47d87cdb518
# ╟─b2797fd2-b6e1-4bb0-af23-26931fe8be69
# ╠═f947c2d9-9123-422d-8972-157717c85b3c
# ╟─74708270-b1ec-48c7-af32-3b970b92c706
# ╠═2b00df5d-994e-47a1-8068-c93ce3f1a618
# ╟─9d06c31e-3525-4889-a1de-3fe02413c7d8
# ╠═9a90f800-3669-4831-b50b-c5405bbb9a03
# ╟─a554fd16-aa3d-48ca-8de6-5582725c27d8
# ╟─d6872046-b5ef-4c2d-a9bb-2418f57f715d
# ╟─747e20c4-b06b-4e78-a09a-55053cf42bf4
# ╟─92181028-60fc-4830-afba-2380ac91455d
# ╠═f8942b10-773a-4b22-baad-8004fba8bd34
# ╟─43d41284-053d-4dfe-8d5b-96be70c0495c
# ╠═a999ae2a-7567-41e7-9c0c-e94fad6f5d46
# ╟─9ebb5b44-04d7-4b89-acdb-e40a245703d2
# ╠═049de8d5-b221-452b-b2c4-9bc1e0c17f48
# ╟─a73a2853-1f48-4179-9771-083794d3f137
# ╠═682e9120-0e1c-4dfa-9ec6-66bb0a3f4374
# ╟─705d3fcb-20b6-4481-a304-1d3ccd623674
# ╠═5968317a-6c07-4655-8137-6702656bb3b4
# ╟─ff9370d8-3395-4382-9f51-afa11748319e
# ╠═48be49d0-0b60-44f3-8152-1ca917a4232e
# ╟─d6452915-bdf0-48f0-8c7d-3df83c7bce72
