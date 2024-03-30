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

# ╔═╡ 50686ef4-e76b-11ee-0815-cf5064e3ccc0
begin
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ bf38a285-6107-4c51-a8d4-8f2098e92300
using Catalyst, DifferentialEquations, Plots, PlutoUI

# ╔═╡ 24f6f3e1-7966-42d5-96bc-4e0e8a4ed03a
rs = @reaction_network begin
	kb*A^order, A => B
end

# ╔═╡ 273208e6-962a-4d74-9aa1-ab042678f4ea
Graph(rs)

# ╔═╡ 32297ce9-956e-4b28-b872-56dd25d853d3
complexgraph(rs)

# ╔═╡ 79f198d1-5688-4739-9728-7042709788f9
convert(ODESystem, rs)

# ╔═╡ a9bd6434-d95e-4851-a9a8-965a2fe41167
@bind order Slider(0.0:1:2.5, default=1.0)

# ╔═╡ 39e88920-97bb-4b43-a1de-3543b762366e
@bind xmax Slider(10:5:100, default=10)

# ╔═╡ df4ec301-683e-43c4-8782-13e13ca58c55
let
	p = (:kb => 0.2, :order => order)
	tspan = (0., 100.)
	u0 = [:A => 10., :B => 0.]
	# solve ODEs
	oprob = ODEProblem(rs, u0, tspan, p)
	# NOTE: callback is used to avoid negative concentrations
	osol  = solve(oprob, Tsit5(), callback = PositiveDomain([1.0, 1.0]; save = false, abstol = 1E-9))
	plot(osol)
	plot!(xlim=(0, xmax))
end

# ╔═╡ Cell order:
# ╠═50686ef4-e76b-11ee-0815-cf5064e3ccc0
# ╠═bf38a285-6107-4c51-a8d4-8f2098e92300
# ╠═24f6f3e1-7966-42d5-96bc-4e0e8a4ed03a
# ╠═273208e6-962a-4d74-9aa1-ab042678f4ea
# ╠═32297ce9-956e-4b28-b872-56dd25d853d3
# ╠═79f198d1-5688-4739-9728-7042709788f9
# ╠═a9bd6434-d95e-4851-a9a8-965a2fe41167
# ╠═39e88920-97bb-4b43-a1de-3543b762366e
# ╠═df4ec301-683e-43c4-8782-13e13ca58c55
