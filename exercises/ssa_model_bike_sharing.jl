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
using Catalyst

# ╔═╡ cc1e8c62-c65e-4e59-8c1b-bc0026dd8e06
using DifferentialEquations, Plots

# ╔═╡ 284f5847-9c15-41f3-a595-1e12a22df69f
using PlutoUI; TableOfContents()

# ╔═╡ eb8a1d29-3da2-489e-8e37-5fbaecc3d973
using Distributions

# ╔═╡ 6c97bf81-ef32-45a4-aa7c-c8c26ba2d2c3
bike_sharing = @reaction_network begin
    @species O(t)=10 W(t)=2
    @parameters p₁=0.50 p₂=0.33
	p₁, O => W
	p₂, W => O
end

# ╔═╡ 1536fe23-0f8d-4b86-98d2-076248b35954
osys = convert(ODESystem, bike_sharing)

# ╔═╡ ab6af765-1cde-4da8-bbc1-a5fab391db54
u0 = [:O => 10, :W => 2]

# ╔═╡ 3ae98e83-7beb-4597-89be-80c813d4349b
tspan = (0.0, 60.0)

# ╔═╡ 0d8f53f8-0a14-4ac6-bd0c-2190d4db0909
@bind p₁ Slider(0.0:0.1:1, default=0.0, show_value=true)

# ╔═╡ e20e4dd8-bdbb-4005-af68-6bf7e4ec130e
parms = [:p₁=>p₁, :p₂=>0.20]

# ╔═╡ d4c45709-70c9-4ba0-8fb8-6b600473723d
dprob = DiscreteProblem(bike_sharing, u0, tspan, parms)

# ╔═╡ 7644adf4-d992-48b1-b40a-12fdf30f6cb5
jdprob = JumpProblem(bike_sharing, dprob, Direct(), save_positions=(false, false));
# https://docs.sciml.ai/JumpProcesses/dev/jump_solve/#JumpProcesses.jl

# ╔═╡ 1511d269-9706-41b1-b8e2-f85ebcedc2d8
function condition(u, t, integrator)
	true
end

# ╔═╡ 3cac94d6-14ae-42f2-b0b6-f47d87cdb518
function affect!(integrator)
	if integrator.u[1] > 12 && integrator.u[2] < 0
		integrator.u[1] = 12
		integrator.u[2] = 0
	end
	if integrator.u[1] < 0 && integrator.u[2] > 12
		integrator.u[1] = 0
		integrator.u[2] = 12
	end
end

# ╔═╡ f947c2d9-9123-422d-8972-157717c85b3c
cb = DiscreteCallback(condition, affect!, save_positions=(false,false));

# ╔═╡ 2b00df5d-994e-47a1-8068-c93ce3f1a618
jdsol = solve(jdprob, SSAStepper(), saveat=1.0, callback=cb);

# ╔═╡ 9a90f800-3669-4831-b50b-c5405bbb9a03
plot(jdsol, ylim=(0, 12))

# ╔═╡ 24af8260-493f-4123-b0dc-928fbadaae16
jdsol.t

# ╔═╡ f8942b10-773a-4b22-baad-8004fba8bd34
jdsol[:O]

# ╔═╡ 049de8d5-b221-452b-b2c4-9bc1e0c17f48
count(jdsol[:O] .== 0)

# ╔═╡ 51d8f821-4aa6-4de8-9072-acee63d4f7ec
jdprob.ps[:p₁]

# ╔═╡ 682e9120-0e1c-4dfa-9ec6-66bb0a3f4374
begin
	p_values = 0.0:0.05:1.0
	mean_zero_counts = []
	for p_val = p_values
		zero_counts_p_val = []
		for i = 1:1000
			jdprob_re = remake(deepcopy(jdprob); p=[:p₁=>p_val])
			jdsol_re = solve(jdprob_re, SSAStepper(), saveat=1.0, callback=cb);
			append!(zero_counts_p_val, count(jdsol_re[:O] .== 0))
		end
		append!(mean_zero_counts, mean(zero_counts_p_val))
	end
end

# ╔═╡ 5968317a-6c07-4655-8137-6702656bb3b4
mean_zero_counts

# ╔═╡ 48be49d0-0b60-44f3-8152-1ca917a4232e
plot(p_values, mean_zero_counts)

# ╔═╡ Cell order:
# ╠═6b342f14-e7d5-11ef-1ea0-77ceb0d78f32
# ╠═7bfcc024-7d7f-4c0b-b918-8f7626b10974
# ╠═309035dd-5653-48a6-a53d-817e743279fa
# ╠═71140c81-af29-4857-8020-4f94c8bd64b3
# ╠═cc1e8c62-c65e-4e59-8c1b-bc0026dd8e06
# ╠═284f5847-9c15-41f3-a595-1e12a22df69f
# ╠═eb8a1d29-3da2-489e-8e37-5fbaecc3d973
# ╠═6c97bf81-ef32-45a4-aa7c-c8c26ba2d2c3
# ╠═1536fe23-0f8d-4b86-98d2-076248b35954
# ╠═ab6af765-1cde-4da8-bbc1-a5fab391db54
# ╠═3ae98e83-7beb-4597-89be-80c813d4349b
# ╠═0d8f53f8-0a14-4ac6-bd0c-2190d4db0909
# ╠═e20e4dd8-bdbb-4005-af68-6bf7e4ec130e
# ╠═d4c45709-70c9-4ba0-8fb8-6b600473723d
# ╠═7644adf4-d992-48b1-b40a-12fdf30f6cb5
# ╠═1511d269-9706-41b1-b8e2-f85ebcedc2d8
# ╠═3cac94d6-14ae-42f2-b0b6-f47d87cdb518
# ╠═f947c2d9-9123-422d-8972-157717c85b3c
# ╠═2b00df5d-994e-47a1-8068-c93ce3f1a618
# ╠═9a90f800-3669-4831-b50b-c5405bbb9a03
# ╠═24af8260-493f-4123-b0dc-928fbadaae16
# ╠═f8942b10-773a-4b22-baad-8004fba8bd34
# ╠═049de8d5-b221-452b-b2c4-9bc1e0c17f48
# ╠═51d8f821-4aa6-4de8-9072-acee63d4f7ec
# ╠═682e9120-0e1c-4dfa-9ec6-66bb0a3f4374
# ╠═5968317a-6c07-4655-8137-6702656bb3b4
# ╠═48be49d0-0b60-44f3-8152-1ca917a4232e
