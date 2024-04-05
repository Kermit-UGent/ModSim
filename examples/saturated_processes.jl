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

# ╔═╡ 54edfaf1-32f9-48de-a7af-769d9e37539b
begin
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ 94b38aa5-a494-4175-a20c-c2406733dd74
using Catalyst, DifferentialEquations, Plots, PlutoUI

# ╔═╡ 9b62ea0d-e88c-41df-997f-a93a2a4c0dc5
md"# Michaelis-Menten"

# ╔═╡ 0cd772f6-bd4f-4a60-b30b-442c9261dbab
begin
	rs_MM = @reaction_network begin
	  c1, S + E --> SE
	  c2, SE --> S + E
	  c3, SE --> P + E
	end
end

# ╔═╡ cf478fd2-b318-4db8-bc88-6778a48cd74a
convert(ODESystem, rs_MM)

# ╔═╡ ef819284-c262-4a55-9064-933845a7aac0
complexgraph(rs_MM)

# ╔═╡ 09e61efe-e608-11ee-31cb-cfff58d95d0d
let
	p = (:c1 => 0.00166, :c2 => 0.0001, :c3 => 0.1)
	tspan = (0., 100.)
	u0 = [:S => 301., :E => 100., :SE => 0., :P => 0.]
	
	# solve ODEs
	oprob = ODEProblem(rs_MM, u0, tspan, p)
	osol  = solve(oprob, Tsit5())
	
	plot(osol; title = "Reaction Rate Equation ODEs")
end

# ╔═╡ 4d340b68-2455-4abc-971a-13a464731005
md"# Monod"

# ╔═╡ 036fab2d-b5a7-4039-a50f-9164f07cdeff
rs_Monod = @reaction_network begin
	@species X(t) S(t)
	@parameters μmax Ks Y
	μmax / (Ks + S) * X * S, (1/Y)*S => X
end

# ╔═╡ 32a9fc59-73b7-413f-bb41-a6140f84d16a
Graph(rs_Monod)

# ╔═╡ be4fe8d2-36d2-4de0-bacd-a0339ffc27b2
convert(ODESystem, rs_Monod)

# ╔═╡ 936d9d89-c319-4be0-966e-073dd68860cd
Graph(rs_Monod)

# ╔═╡ e769dea8-e336-4e4c-875e-3d6f7d7325dd
md"""
### Maximum Specific Growth Rate ($\mu\_{\\max}$)

*   For bacteria, $\mu\_{\\max}$ typically ranges from 0.10 to 1.0 per hour. Fast-growing bacteria in optimal conditions can have higher rates.
*   For yeast and fungi, the range is generally wider, from 0.03 to 0.40 per hour, as they usually grow slower than bacteria.

### Half-Saturation Constant ($K\_s$)

*   The $K\_s$ value is indicative of the affinity of the microorganism for the substrate: lower values suggest a high affinity. For common substrates like glucose, $K\_s$ might range from 10 to 100 mg/L for various bacteria and yeasts.
*   It's important to note that $K\_s$ can vary widely depending on the substrate and the microorganism. For example, for more complex substrates or those that are less readily metabolized, $K\_s$ can be significantly higher.

### Yield Coefficient ($Y$)

*   The yield coefficient, indicating the amount of biomass produced per unit of substrate consumed, often ranges from 0.30 to 0.60 biomass/g substrate for many bacterial growth systems.
*   The specific value of Y depends on the efficiency of substrate utilization and the metabolic pathways involved. For example, aerobic conditions tend to result in higher yield coefficients compared to anaerobic conditions due to the more efficient energy extraction from the substrate.
"""

# ╔═╡ 548c71ce-c661-468c-baac-554901b74dc5
@bind μmax Slider(0.01:0.001:0.1, default=0.1)

# ╔═╡ c6133ecb-2324-484d-8025-c2e9d8a85c3c
@bind Ks Slider(0.1:0.1:100., default=10.)

# ╔═╡ 18eb4e7d-4d73-42b7-a210-34386134e365
@bind Y Slider(0.3:0.01:0.6, default=0.3)

# ╔═╡ 68f571da-3f78-482c-9d62-90fce7924356
let
	# time and reaction rate in hours!
	p = (:μmax => μmax, :Ks => Ks, :Y => Y)
	tspan = (0., 100.)
	u0 = [:S => 100., :X => 10.]
	
	# solve ODEs
	oprob = ODEProblem(rs_Monod, u0, tspan, p)
	# NOTE: callback is used to avoid negative concentrations
	osol  = solve(oprob, Tsit5(), callback = PositiveDomain([1.0, 1.0]; save = false, abstol = 1E-9))
	
	plot(osol; title = "Reaction Rate Equation ODEs")
end

# ╔═╡ Cell order:
# ╠═54edfaf1-32f9-48de-a7af-769d9e37539b
# ╠═94b38aa5-a494-4175-a20c-c2406733dd74
# ╟─9b62ea0d-e88c-41df-997f-a93a2a4c0dc5
# ╠═0cd772f6-bd4f-4a60-b30b-442c9261dbab
# ╠═cf478fd2-b318-4db8-bc88-6778a48cd74a
# ╠═32a9fc59-73b7-413f-bb41-a6140f84d16a
# ╠═ef819284-c262-4a55-9064-933845a7aac0
# ╠═09e61efe-e608-11ee-31cb-cfff58d95d0d
# ╟─4d340b68-2455-4abc-971a-13a464731005
# ╠═036fab2d-b5a7-4039-a50f-9164f07cdeff
# ╠═be4fe8d2-36d2-4de0-bacd-a0339ffc27b2
# ╠═936d9d89-c319-4be0-966e-073dd68860cd
# ╟─e769dea8-e336-4e4c-875e-3d6f7d7325dd
# ╠═548c71ce-c661-468c-baac-554901b74dc5
# ╠═c6133ecb-2324-484d-8025-c2e9d8a85c3c
# ╠═18eb4e7d-4d73-42b7-a210-34386134e365
# ╠═68f571da-3f78-482c-9d62-90fce7924356
