### A Pluto.jl notebook ###
# v0.20.3

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

# ╔═╡ 99ff6d3d-0a9c-485c-87b0-d361080fd7ca
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	using Pkg
	Pkg.activate("..")
end
  ╠═╡ =#

# ╔═╡ 1e42a500-f18b-11ee-1d1b-c54e26589084
using PlutoUI, LaTeXStrings, Plots, LinearAlgebra

# ╔═╡ 6f264016-64c4-4680-a0d3-6f161b7a6d67
using Distributions

# ╔═╡ 279d873f-292e-4590-b4d0-0b933bcf325e
using Turing

# ╔═╡ fb45cbb5-1117-4201-ab57-a6c2cdf96599
using Optim

# ╔═╡ ce7c6458-edad-49a6-9742-2377eb2a1071
f_rosenbrock((x1, x2); a=1, b=100) = (a-x1)^2 + b * (x2-x1^2)^2

# ╔═╡ d4af97db-51c2-43a6-9ffc-d8269af53a71
f_rosenbrock(x1,x2) = f_rosenbrock((x1, x2))

# ╔═╡ c2a0d824-6e0e-4594-9755-3521deac9051
f = x->f_rosenbrock(x)

# ╔═╡ b7ea4943-a8e1-40ea-90bd-9220aa0eb901
md"## Grid search and random search"

# ╔═╡ 3927ad66-4111-494c-9cb7-9c659ec57151
x_grid = argmin(f_rosenbrock, ((x1, x2) for x1 in range(-2, 2, length=10)
							for x2 in range(-1, 3, length=10)))

# ╔═╡ 71a0c6ef-c1a5-4ca3-8114-e496bc24813e
distr = product_distribution(Uniform(-2,2), Uniform(-1, 3))

# ╔═╡ 1e5398ef-8666-4371-8fbd-78e5b35b3253
x_rand = argmin(f, (rand(distr) for i in 1:100))

# ╔═╡ ad7a5769-f325-4b6e-8a5e-04a65478dc0f
@model function sampling_dist()
	x ~ Uniform(-2, 2)
	y ~ Uniform(-1, 3)
end

# ╔═╡ 3308c8df-93a0-4eea-b45d-34859d5b2c0d
argmin(f, (rand(sampling_dist()) for i in 1:100))

# ╔═╡ edb750a6-291a-42bb-9cba-b289f8ba801a
@model function sampling_dist2()
	x ~ Normal(0, 1)
	y ~ Normal(1, 1)
end

# ╔═╡ d7f8e367-4404-4d3f-98ef-a5fcb9436566
argmin(f, (rand(sampling_dist2()) for i in 1:100))

# ╔═╡ caa25351-585b-4d4d-a51f-2a7ee68f6bbe
x0 = [-1.0, 2.0]

# ╔═╡ f06fbf86-05fe-44d5-9dc2-e17bb0347597
opts = Optim.Options(store_trace=true, iterations=10_000, extended_trace=true);

# ╔═╡ b0ded35f-3a77-479f-abf8-1f98b9c0bb7a
sol_GD = optimize(f_rosenbrock, x0, GradientDescent(), opts, autodiff = :forward)

# ╔═╡ 42c40178-ea40-44ce-a416-aefc1ae03662
sol_GD.minimizer

# ╔═╡ 2adcab36-3436-443b-a11a-c25e2dac19c1
# ╠═╡ skip_as_script = true
#=╠═╡
p = sol_GD.trace[1]
  ╠═╡ =#

# ╔═╡ b2aaff23-c1d4-4f86-978f-32e592ac7072
# ╠═╡ skip_as_script = true
#=╠═╡
p.value
  ╠═╡ =#

# ╔═╡ 9d34850a-4069-42d5-abdf-4f8e7199a947
# ╠═╡ skip_as_script = true
#=╠═╡
p.metadata["x"]
  ╠═╡ =#

# ╔═╡ af2940b8-991b-4231-887a-ae48c60212d6
[p.metadata["x"] for p in sol_GD.trace]

# ╔═╡ f5c307dc-7610-4f6a-abda-191443d3901c
sol_GD.time_run

# ╔═╡ eb974079-2187-4d2d-b2c7-19830d5ce6ee
sol_newton = optimize(f, x0, Newton(), opts, autodiff = :forward)

# ╔═╡ 0640276b-c6d8-46cf-a928-79a425306319
sol_newton.time_run

# ╔═╡ f23d0a8c-5045-4071-be2d-389f1a796484


# ╔═╡ 49718e41-61ea-4934-918b-ef03a259d609
@bind θ Slider(0:0.1:(2π), default=deg2rad(67))

# ╔═╡ d6e95454-b2ea-453a-aad2-75ce39f1a6e3
θ

# ╔═╡ 3b285e28-d047-499a-99d8-0dde5c80cab2
R = [cos(θ) -sin(θ); 
	sin(θ) cos(θ)]

# ╔═╡ 414a8764-2310-4a5a-90ba-ed7b126fb4f8
rot = x -> R * x 

# ╔═╡ 082e119f-5a2b-44f4-8357-55f518c0a64d
frot(x, y) = f_rosenbrock(rot([x, y]))

# ╔═╡ 4392fa7e-3e56-4944-aee2-66150a02e5f4
frot((x,y)) = frot(x, y)

# ╔═╡ b7dda6c0-1b52-4f30-a9b8-427d6768d793
optimize(frot, rot(x0), Newton(), opts, autodiff=:forward)

# ╔═╡ 376db0fa-cfc6-48c6-a65d-f58e3afb24f7
rot(x0)

# ╔═╡ f4e05e34-8aba-4610-8a2d-fcf6752e0f1d


# ╔═╡ da830667-9ef4-4b88-a07b-9ec50dd4e815
frescaled(x, y) = f_rosenbrock(10x, y/10)

# ╔═╡ 3556cd4c-eb43-48f8-84a6-80246683a8bc
frescaled((x, y)) = frescaled(x, y)

# ╔═╡ 4047e168-0906-4704-a191-1d4d6a0adba4
optimize(frescaled, [x0[1]/10, x0[2]*10], Newton(), opts)

# ╔═╡ 0267e113-71d9-4b26-bed6-41ceae37c1ce
optimize(frescaled, [x0[1]/10, x0[2]*10], GradientDescent(), opts)

# ╔═╡ f053529a-2726-48c6-a848-f6bc220ee1f0
sol_LBFGS = optimize(f, x0, LBFGS(), opts, autodiff=:forward)

# ╔═╡ 7fd70049-b91c-41d1-bf34-688d3461c9df
sol_LBFGS.trace

# ╔═╡ 4e4e6c03-f2dd-4f7f-8c00-82d47286d019
sol_LBFGS.time_run

# ╔═╡ b60d91f5-f635-4dc0-95fd-a35fbac97978
sol_NM = optimize(f, x0, NelderMead(), opts)

# ╔═╡ 63231b16-b6f5-41be-8039-4dd3b254c31e
sol_NM.time_run

# ╔═╡ 171dab60-921a-40e3-b215-6ba712a780eb
sol_NM.trace

# ╔═╡ a1f5e697-3363-4406-8dd4-8dd66c76f23c
sol_PSO = optimize(f, x0, ParticleSwarm(lower=[-2.0, -1.0],
						upper=[2.0, 3.0], n_particles=20), opts)

# ╔═╡ ee5f554b-ead2-40bb-bfe8-b1cc7053a070
#=╠═╡
typeof(p)
  ╠═╡ =#

# ╔═╡ 6bbee6f5-ae16-428e-86bd-67b7b35ed4c5
@model function g(x)
    m ~ Normal(0, 1)
	b ~ Exponential( 10)
    x ~ Normal(exp(-m), b)
end

# ╔═╡ 523a9f90-4a02-4716-abb5-fda48c948816
optimize(g(10), MLE())

# ╔═╡ cf1a8701-ec6a-4dfa-b72f-64d4168f332e
optimize(g(10), MAP())

# ╔═╡ a8fcc4bb-0360-402a-8992-5c5883138b11
md"# Appendix"

# ╔═╡ 17e6402f-2bc0-41b1-ab62-54b8edd6ffa8
plots = Dict()

# ╔═╡ 90125dda-fa19-415c-a8d7-eb37bc1dbf14
begin
	 p_ros = contourf(-2:0.01:2, -1:0.01:3, log10 ∘ f_rosenbrock, color= cgrad(:speed), xlab=L"x_1", ylab=L"x_2", title="Rosenbrock function", xlims=[-2,2], ylims=[-1, 3])
	plots["rosenbrock"] = scatter!([1], [1], label="", color="red")
end

# ╔═╡ 1adee9b4-7cbe-4d14-8004-8cf54613231d
begin
	 p_ros2 = contour(-2:0.01:2, -1:0.01:3, log10 ∘ f_rosenbrock, color= cgrad(:speed), xlab=L"x_1", ylab=L"x_2", title="Rosenbrock function", xlims=[-2,2], ylims=[-1, 3])
	plots["rosenbrock_light"] = scatter!([1], [1], label="", color="red")
end

# ╔═╡ bd20d253-9238-4260-ae87-159dc34e10d9
let
	p = contour(-2:0.01:2, -1:0.01:3, log10 ∘ f_rosenbrock, color= cgrad(:speed), xlab=L"x_1", ylab=L"x_2", title="Rosenbrock function", xlims=[-2,2], ylims=[-1, 3])
	scatter!([1], [1], label=L"x^\star", color="red")
	grid = [(x, y) for x in range(-2, 2, length=10)
							for y in range(-1, 3, length=10)]
	plots["grid_search"] = scatter!(first.(grid), last.(grid), title="Grid search",
			color="orange", label="", ms=3, m=:x)
end

# ╔═╡ adb796f0-8296-4451-9086-8ba731f80608
let
	p = contour(-2:0.01:2, -1:0.01:3, log10 ∘ f_rosenbrock, color= cgrad(:speed), xlab=L"x_1", ylab=L"x_2", title="Rosenbrock function", xlims=[-2,2], ylims=[-1, 3])
	scatter!([1], [1], label=L"x^\star", color="red")
	grid = [rand(distr) for i in 1:100]
	plots["random_search"] = scatter!(first.(grid), last.(grid), title="Random search (uniform)",
			color="orange", label="", ms=3, m=:x)
end

# ╔═╡ 6391bbce-bdad-4833-82a9-3cde82fa0952
let
	contour(-2:0.01:2, -1:0.01:3, log10 ∘ f_rosenbrock, color= cgrad(:speed), xlab=L"x_1", ylab=L"x_2", title="Rosenbrock function", xlims=[-2,2], ylims=[-1, 3])
	scatter!([1], [1], label="", color="red")
	grid = [rand(sampling_dist2()) for i in 1:100]
	plots["random_search_norm"] = scatter!(first.(grid), last.(grid), title="Random search (non-uniform)",
			color="orange", label="", ms=3, m=:x)
end

# ╔═╡ 462774d1-c1c5-48d7-8f25-2c091ad86943
let
	p = contour(-2:0.01:2, -1:0.01:3, log10 ∘ f_rosenbrock, color= cgrad(:speed), xlab=L"x", ylab=L"y", title="Stratified sampling", xlims=[-2,2], ylims=[-1, 3])
	vline!(range(-2, 2, length=7), ls=:dash, color=:blue, alpha=0.7, label="")
	hline!(range(-1, 3, length=7), ls=:dash, color=:blue, alpha=0.7, label="")
	points = []
	for a in range(-2, length=6, step=2/3)
		for b in range(-1, length=6, step=2/3)
			x = 2rand()/3 + a
			y = 2rand()/3 + b
			push!(points, (x, y))
		end
	end
	plots["stratified_sampling"] = scatter!(first.(points), last.(points), color="orange", label="", ms=3, m=:x)
end

# ╔═╡ ee76b80d-64b2-4e6a-8384-d3179c996d55
let
	sol = sol_GD
	x1, x2 = sol.minimizer
	p = contour(-2:0.01:2, -1:0.01:3, log10 ∘ f_rosenbrock, color= cgrad(:speed), xlab=L"x_1", ylab=L"x_2", title="Gradient descent", xlims=[-2,2], ylims=[-1, 3])
	scatter!([1], [1], label=L"x^\star", color="red")
	scatter!([x1], [x2], label="minimizer", color=:gold)
	path = [p.metadata["x"] for p in sol.trace]
	plot!(first.(path), last.(path), lw=2, color=:blue, alpha=0.7, label="path")
	plots["Rosenbrock_GD"] = p
end

# ╔═╡ 320939b3-617e-43f8-a915-7b01637cf3d0
let
	sol = sol_newton
	x1, x2 = sol.minimizer
	p = contour(-2:0.01:2, -1:0.01:3, log10 ∘ f_rosenbrock, color= cgrad(:speed), xlab=L"x_1", ylab=L"x_2", title="Newton's method", xlims=[-2,2], ylims=[-1, 3])
	scatter!([1], [1], label=L"x^\star", color="red")
	scatter!([x1], [x2], label="minimizer", color=:gold)
	path = [p.metadata["x"] for p in sol.trace]
	plot!(first.(path), last.(path), lw=2, color=:orange, alpha=0.7, label="path")
	plots["Rosenbrock_Newton"] = p
end

# ╔═╡ 4eaa94d1-00ff-4baa-91ff-0bd380f30fb1
let
	quadr(x) = x' * A * x / 2
	A = [5 -2; -2 1]
	x = [-1., 3.]
	x1, x2 = x
	p = contourf(-5:.01:5, -5:.01:5, (x,y)->quadr([x,y]), color=:speed, aspect_ratio=:equal, xlims=(-5,5), ylims=(-5,5))
	title!("GD vs. Newton step")
	scatter!([0], [0], label=L"x^\star", color=:red)
	dx_gd = - A * x
	dx_gd .*= 2/norm(dx_gd)
	dx_newton = - x
	dx_newton .*= 2/norm(dx_newton) 
	quiver!([x1], [x2], quiver=([dx_gd[1]], [dx_gd[2]]), label=L"\delta \mathbf{x}", lw=2, color=:blue2)
	quiver!([x1], [x2], quiver=([dx_newton[1]], [dx_newton[2]]), lw=2, ls=:dash, color=:red)
	xlabel!(L"x_1")
	ylabel!(L"x_2")
	title!("Negative gradient vs\nNewton step")
	plots["steps"] = p
end

# ╔═╡ f8e0b77b-1d8b-4ae8-be97-1a233c8a94d2
let
	sol = optimize(frot, rot(x0), Newton(), opts)
	x1, x2 = sol.minimizer
	p = contour(-5:0.01:5, -3:0.01:3, log10 ∘ frot, color= cgrad(:speed), xlab=L"x_1", ylab=L"x_2", title="Newton's method", xlims=[-5,5], ylims=[-3, 3])
	x1m, x2m = rot([1, 1])
	scatter!([x1m], [x2m], label=L"x^\star", color="red")
	scatter!([x1], [x2], label="minimizer", color=:gold)
	path = [p.metadata["x"] for p in sol.trace]
	plot!(first.(path), last.(path), lw=2, color=:orange, alpha=0.7, label="path")
	plots["Rosenbrock_Newton_rot"] = p
end

# ╔═╡ 99bd0581-596b-4390-8248-d0de51caceac
let
	sol = optimize(frescaled, [x0[1]/10, x0[2]*10], Newton(), opts)
	x1, x2 = sol.minimizer
	p = contour(-.2:0.001:.2, -10:0.1:30, log10 ∘ frescaled, color= cgrad(:speed), xlab=L"x_1", ylab=L"x_2", title="Newton's method (rescaled)", xlims=[-.2,.2], ylims=[-10, 30])
	scatter!([.1], [10], label=L"x^\star", color="red")
	scatter!([x1], [x2], label="minimizer", color=:gold)
	path = [p.metadata["x"] for p in sol.trace]
	plot!(first.(path), last.(path), lw=2, color=:orange, alpha=0.7, label="path")
	plots["Rosenbrock_Newton_rescaled"] = p
end

# ╔═╡ 60b92e5e-64bd-489b-8b87-873a5db3d7b1
let
	sol = optimize(frescaled, [x0[1]/10, x0[2]*10], GradientDescent(), opts)
	x1, x2 = sol.minimizer
	p = contour(-.2:0.001:.2, -10:0.1:30, log10 ∘ frescaled, color= cgrad(:speed), xlab=L"x_1", ylab=L"x_2", title="Gradient descent (rescaled)", xlims=[-.2,.2], ylims=[-10, 30])
	scatter!([.1], [10], label=L"x^\star", color="red")
	scatter!([x1], [x2], label="minimizer", color=:gold)
	path = [p.metadata["x"] for p in sol.trace]
	plot!(first.(path), last.(path), lw=2, color=:blue, alpha=0.7, label="path")
	plots["Rosenbrock_GD_rescaled"] = p
end

# ╔═╡ 401783fe-d82f-4c3c-84ad-ef982aa34f7d
let
	sol = sol_LBFGS
	x1, x2 = sol.minimizer
	p = contour(-2:0.01:2, -1:0.01:3, log10 ∘ f_rosenbrock, color= cgrad(:speed), xlab=L"x_1", ylab=L"x_2", title="L-BFGS", xlims=[-2,2], ylims=[-1, 3])
	scatter!([1], [1], label=L"x^\star", color="red")
	scatter!([x1], [x2], label="minimizer", color=:gold)
	path = [p.metadata["x"] for p in sol.trace]
	plot!(first.(path), last.(path), lw=2, color=:red, alpha=0.7, label="path")
	plots["Rosenbrock_LBGS"] = p
end

# ╔═╡ d8ca4713-5063-47f4-8763-ae72694dac0b
let
	p = plot(title="Convergence rates\nRosenbrock function", xscale=:log10, yscale=:log10, xlab=L"k+1", ylab=L"f(x^{(k)})-f(x^\star)", ls=:auto)
	f_gd = [p.value for p in sol_GD.trace]
	plot!(f_gd, label="gradient descent", lw=2, color=:blue, ls=:solid)
	f_newton = [p.value for p in sol_newton.trace]
	plot!(f_newton, label="Newton's method", lw=2, color=:orange, ls=:dash)
	f_bfgs = [p.value for p in sol_LBFGS.trace]
	plot!(f_bfgs, label="L-BFGS", lw=2, color=:red, ls=:dot)
	plots["convergence_FO"] = p
end

# ╔═╡ f3baf67a-2085-4601-99ee-9dde9ac44dec
let
	sol = sol_NM
	x1, x2 = sol.minimizer
	p = contour(-2:0.01:2, -1:0.01:3, log10 ∘ f_rosenbrock, color= cgrad(:speed), xlab=L"x_1", ylab=L"x_2", title="Nelder-Mead", xlims=[-2,2], ylims=[-1, 3])
	scatter!([1], [1], label=L"x^\star", color="red")
	scatter!([x1], [x2], label="minimizer", color=:gold)
	path = [p.metadata["centroid"] for p in sol.trace]
	scatter!(first.(path), last.(path), color="orange", ms=3, m=:x, label="evaluations")
	plots["Rosenbrock_NM"] = p
end

# ╔═╡ ef8c9e9e-974f-4970-b48a-6cde8ccd586b
plots

# ╔═╡ Cell order:
# ╠═1e42a500-f18b-11ee-1d1b-c54e26589084
# ╠═99ff6d3d-0a9c-485c-87b0-d361080fd7ca
# ╠═ce7c6458-edad-49a6-9742-2377eb2a1071
# ╠═d4af97db-51c2-43a6-9ffc-d8269af53a71
# ╠═c2a0d824-6e0e-4594-9755-3521deac9051
# ╠═90125dda-fa19-415c-a8d7-eb37bc1dbf14
# ╠═1adee9b4-7cbe-4d14-8004-8cf54613231d
# ╠═b7ea4943-a8e1-40ea-90bd-9220aa0eb901
# ╠═3927ad66-4111-494c-9cb7-9c659ec57151
# ╠═bd20d253-9238-4260-ae87-159dc34e10d9
# ╠═6f264016-64c4-4680-a0d3-6f161b7a6d67
# ╠═71a0c6ef-c1a5-4ca3-8114-e496bc24813e
# ╠═1e5398ef-8666-4371-8fbd-78e5b35b3253
# ╠═adb796f0-8296-4451-9086-8ba731f80608
# ╠═279d873f-292e-4590-b4d0-0b933bcf325e
# ╠═ad7a5769-f325-4b6e-8a5e-04a65478dc0f
# ╠═3308c8df-93a0-4eea-b45d-34859d5b2c0d
# ╠═edb750a6-291a-42bb-9cba-b289f8ba801a
# ╠═d7f8e367-4404-4d3f-98ef-a5fcb9436566
# ╠═6391bbce-bdad-4833-82a9-3cde82fa0952
# ╠═462774d1-c1c5-48d7-8f25-2c091ad86943
# ╠═fb45cbb5-1117-4201-ab57-a6c2cdf96599
# ╠═caa25351-585b-4d4d-a51f-2a7ee68f6bbe
# ╠═f06fbf86-05fe-44d5-9dc2-e17bb0347597
# ╠═b0ded35f-3a77-479f-abf8-1f98b9c0bb7a
# ╠═42c40178-ea40-44ce-a416-aefc1ae03662
# ╠═2adcab36-3436-443b-a11a-c25e2dac19c1
# ╠═b2aaff23-c1d4-4f86-978f-32e592ac7072
# ╠═9d34850a-4069-42d5-abdf-4f8e7199a947
# ╠═af2940b8-991b-4231-887a-ae48c60212d6
# ╠═ee76b80d-64b2-4e6a-8384-d3179c996d55
# ╠═f5c307dc-7610-4f6a-abda-191443d3901c
# ╠═eb974079-2187-4d2d-b2c7-19830d5ce6ee
# ╠═0640276b-c6d8-46cf-a928-79a425306319
# ╠═320939b3-617e-43f8-a915-7b01637cf3d0
# ╠═f23d0a8c-5045-4071-be2d-389f1a796484
# ╠═49718e41-61ea-4934-918b-ef03a259d609
# ╠═d6e95454-b2ea-453a-aad2-75ce39f1a6e3
# ╠═3b285e28-d047-499a-99d8-0dde5c80cab2
# ╠═414a8764-2310-4a5a-90ba-ed7b126fb4f8
# ╠═082e119f-5a2b-44f4-8357-55f518c0a64d
# ╠═4392fa7e-3e56-4944-aee2-66150a02e5f4
# ╠═b7dda6c0-1b52-4f30-a9b8-427d6768d793
# ╠═376db0fa-cfc6-48c6-a65d-f58e3afb24f7
# ╠═f4e05e34-8aba-4610-8a2d-fcf6752e0f1d
# ╠═4eaa94d1-00ff-4baa-91ff-0bd380f30fb1
# ╠═f8e0b77b-1d8b-4ae8-be97-1a233c8a94d2
# ╠═da830667-9ef4-4b88-a07b-9ec50dd4e815
# ╠═3556cd4c-eb43-48f8-84a6-80246683a8bc
# ╠═4047e168-0906-4704-a191-1d4d6a0adba4
# ╠═99bd0581-596b-4390-8248-d0de51caceac
# ╠═0267e113-71d9-4b26-bed6-41ceae37c1ce
# ╠═60b92e5e-64bd-489b-8b87-873a5db3d7b1
# ╠═f053529a-2726-48c6-a848-f6bc220ee1f0
# ╠═7fd70049-b91c-41d1-bf34-688d3461c9df
# ╠═401783fe-d82f-4c3c-84ad-ef982aa34f7d
# ╠═4e4e6c03-f2dd-4f7f-8c00-82d47286d019
# ╠═d8ca4713-5063-47f4-8763-ae72694dac0b
# ╠═b60d91f5-f635-4dc0-95fd-a35fbac97978
# ╠═63231b16-b6f5-41be-8039-4dd3b254c31e
# ╠═171dab60-921a-40e3-b215-6ba712a780eb
# ╠═f3baf67a-2085-4601-99ee-9dde9ac44dec
# ╠═a1f5e697-3363-4406-8dd4-8dd66c76f23c
# ╠═ee5f554b-ead2-40bb-bfe8-b1cc7053a070
# ╠═6bbee6f5-ae16-428e-86bd-67b7b35ed4c5
# ╠═523a9f90-4a02-4716-abb5-fda48c948816
# ╠═cf1a8701-ec6a-4dfa-b72f-64d4168f332e
# ╠═a8fcc4bb-0360-402a-8992-5c5883138b11
# ╠═17e6402f-2bc0-41b1-ab62-54b8edd6ffa8
# ╠═ef8c9e9e-974f-4970-b48a-6cde8ccd586b
