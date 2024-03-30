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

# ╔═╡ e719b7c4-c7f9-4d50-9a87-0662020bc2a1
begin
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ f88a1d39-4304-41b3-82b3-2534412b39eb
using Catalyst, DifferentialEquations, Plots, PlutoUI

# ╔═╡ a9620163-a5b4-4fb6-85aa-a3a1293b3571
@bind α Slider(0.0004:0.0001:0.02, default=0.01)

# ╔═╡ 557289b4-a45d-448d-9c24-97dc15932393
@bind β Slider(0.001:0.001:0.05, default=0.01)

# ╔═╡ 10362f30-e606-11ee-0464-53a08feae6cf
begin
	rs = @reaction_network begin
	    α, S + I --> 2I
	    β, I --> R
	end
	# p = [:α => .1/100, :β => .01]
	p = [:α => α, :β => β]
	tspan = (0.0,500.0)
	u0 = [:S => 99.0, :I => 1.0, :R => 0.0]
	
	# Solve ODEs.
	oprob = ODEProblem(rs, u0, tspan, p)
	osol = solve(oprob)
	
	# Solve Jumps.
	dprob = DiscreteProblem(rs, u0, tspan, p)
	jprob = JumpProblem(rs, dprob, Direct())
	jsol = solve(jprob, SSAStepper())
	
	plot(plot(osol; title = "Reaction Rate Equation ODEs", legend=:outertopright),
	     plot(jsol; title = "Gillespie Jump Simulation", legend=:outertopright);
	     layout = (2, 1))
end

# ╔═╡ Cell order:
# ╠═e719b7c4-c7f9-4d50-9a87-0662020bc2a1
# ╠═f88a1d39-4304-41b3-82b3-2534412b39eb
# ╠═a9620163-a5b4-4fb6-85aa-a3a1293b3571
# ╠═557289b4-a45d-448d-9c24-97dc15932393
# ╠═10362f30-e606-11ee-0464-53a08feae6cf
