### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 14176146-cb6c-11ee-008f-5d711daf184d
begin
    using Pkg
	Pkg.activate("..")
    using Catalyst, Plots,  DifferentialEquations, PlutoUI
end


# ╔═╡ e7f76855-48a4-4477-91c7-b9c27d234a37
bdp = @reaction_network begin
	@parameters η1 η2 η3 η4 
  b, 0 --> B  # money enters on the bank
  f, B --> S  # invest
  μ, S --> 2S  # stocks or index funds
  r, S --> B  # return back to the bank
end

# ╔═╡ e379260c-80b7-4855-b7aa-69e813ab90e6
@variables T

# ╔═╡ fc9dea51-def4-4885-b392-a7fc85ac9b66
u₀ = [:S => 10.0, :B=>0.0]

# ╔═╡ 2475e3ec-7380-450d-bb57-2a7421637f29
tspan = (0., 5.)

# ╔═╡ c503365a-13e0-465d-b804-50f37e0199e0
convert(SDESystem, bdp)

# ╔═╡ 3cfd8147-1bd2-49bc-b806-b3acb784e8c8
@bind f Slider(0:0.01:10, show_value=true)

# ╔═╡ 3411d1a1-c0a3-42e1-b571-488bb5c029a3
p = (:μ => log(1.2), :r=>0.0, :b => 0, :f=>0f,
			:η1 =>0.0, :η2 => 0.0, :η3 => 2.175, :η4 => 0.)

# ╔═╡ 07969ca9-688a-4a45-a60e-8934e37c1194
sprob = SDEProblem(bdp, u₀, tspan, p, noise_scaling = @parameters η1 η2 η3 η4 )

# ╔═╡ 0b31a796-9f29-481e-a1a1-2aa83e5e4c1b
sol = solve(sprob, LambaEM(), tstops = range(0., 15, step = 1e-2))

# ╔═╡ a5697d4d-446b-43f6-abf1-c7e78a9d599e
plot(sol)

# ╔═╡ 7d30c7a2-a32b-4417-af0a-1627412380b4
eprob = EnsembleProblem(sprob)

# ╔═╡ c627d93f-1ee5-4978-9102-52eb3efc0061
esol = solve(eprob; trajectories=50, tstops = range(0., 4, step = 1e-2))

# ╔═╡ abf0e148-b623-4a03-b65c-af9127c6ddda
summ = EnsembleSummary(esol)

# ╔═╡ cff5ed3f-3d6a-41fd-99fc-d4e54e5f6e4d
plot(summ, idx=2)

# ╔═╡ 65d7ab2b-d1ae-4a51-9b9c-fea2e53501b7
plot(esol)

# ╔═╡ bcb92576-178f-475b-81d4-e9cf8235bbee


# ╔═╡ 58047756-417f-4d40-be84-4970c0588e48
let 

	rn_2 = @reaction_network begin
    (k1,k2), X1 <--> X2
end
u0 = [:X1 => 10.0, :X2 => 10.0]
tspan = (0.0, 10.0)
p_2 = [:k1 => 1.0, :k2 => 1.0]

sprob_2 = SDEProblem(rn_2, u0, tspan, p_2)
plot(solve(sprob_2))
end

# ╔═╡ b5d7cef4-aced-4f5b-bf2c-fae018a2f417
let

function f(u, p, t)
	du = zeros(2)
	du[1] = 1. + 0.2 * u[2] - p * u[1]
	du[2] = u[2] * log(1.08) - p * u[2] + 0.2 * u[1]
	return du
end

function g(u, p, t)
	du = zeros(2)
	du[1] = 0
	du[2] = 1 * u[2]
	return du
end
	α = 1
β = 1
u₀ = [0.0, 0.0]

prob = SDEProblem(f, g, u₀, tspan)

	plot(solve(prob))
end

# ╔═╡ a67fd056-8921-4ec6-971d-101b39e276d0


# ╔═╡ 6ea78c27-1e60-4d26-9231-e168fdf2206e


# ╔═╡ 842b4e73-ff8c-40e2-bd27-f1c76cf156c9
let
	α = 1
β = 1
u₀ = 1 / 2
f(u, p, t) = log(1.08) * u + 1
g(u, p, t) = 0 * u
dt = 1 // 2^(4)
tspan = (0.0, 5.0)
prob = SDEProblem(f, g, u₀, (0.0, 5.0))

	plot(solve(prob))
end

# ╔═╡ edd31323-6975-4e02-bdb9-06c69fcaee09


# ╔═╡ 68527b74-28a4-445a-bfba-40b4752c87d5


# ╔═╡ e4d97ea3-59f4-4827-9c13-be6bea819d80


# ╔═╡ 2bc62d6b-9613-4f13-9104-1d26aafd99e8


# ╔═╡ ba011284-d2d3-4117-a545-42d6945c1c7b


# ╔═╡ b937c06d-fbc2-4897-ab4b-9d1ebb7b3192


# ╔═╡ Cell order:
# ╠═14176146-cb6c-11ee-008f-5d711daf184d
# ╠═e7f76855-48a4-4477-91c7-b9c27d234a37
# ╠═3411d1a1-c0a3-42e1-b571-488bb5c029a3
# ╠═e379260c-80b7-4855-b7aa-69e813ab90e6
# ╠═fc9dea51-def4-4885-b392-a7fc85ac9b66
# ╠═2475e3ec-7380-450d-bb57-2a7421637f29
# ╠═c503365a-13e0-465d-b804-50f37e0199e0
# ╠═07969ca9-688a-4a45-a60e-8934e37c1194
# ╠═0b31a796-9f29-481e-a1a1-2aa83e5e4c1b
# ╠═a5697d4d-446b-43f6-abf1-c7e78a9d599e
# ╠═3cfd8147-1bd2-49bc-b806-b3acb784e8c8
# ╠═7d30c7a2-a32b-4417-af0a-1627412380b4
# ╠═c627d93f-1ee5-4978-9102-52eb3efc0061
# ╠═abf0e148-b623-4a03-b65c-af9127c6ddda
# ╠═cff5ed3f-3d6a-41fd-99fc-d4e54e5f6e4d
# ╠═65d7ab2b-d1ae-4a51-9b9c-fea2e53501b7
# ╠═bcb92576-178f-475b-81d4-e9cf8235bbee
# ╠═58047756-417f-4d40-be84-4970c0588e48
# ╠═b5d7cef4-aced-4f5b-bf2c-fae018a2f417
# ╠═a67fd056-8921-4ec6-971d-101b39e276d0
# ╠═6ea78c27-1e60-4d26-9231-e168fdf2206e
# ╠═842b4e73-ff8c-40e2-bd27-f1c76cf156c9
# ╠═edd31323-6975-4e02-bdb9-06c69fcaee09
# ╠═68527b74-28a4-445a-bfba-40b4752c87d5
# ╠═e4d97ea3-59f4-4827-9c13-be6bea819d80
# ╠═2bc62d6b-9613-4f13-9104-1d26aafd99e8
# ╠═ba011284-d2d3-4117-a545-42d6945c1c7b
# ╠═b937c06d-fbc2-4897-ab4b-9d1ebb7b3192
