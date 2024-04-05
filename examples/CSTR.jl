### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 4ad03253-636a-450c-863f-b7c8f4367bf5
begin
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ ad5fede2-e600-11ee-374e-85534f76407e
using Catalyst, DifferentialEquations, Plots

# ╔═╡ 077fef83-15c6-48cd-9e9f-4ee2644a0c2b
begin
	# generate a reaction network resembling 3 CSTRs
	# with A connected to the input, C to the output, and B some dead volume tank
	# check the call to `complexgraph(rs)` below, that will make it clear
	rs = @reaction_network begin
	    # note that using @parameters and @species is optional, this is the explicit declaration
	    # of the parameters and species in the reaction network, they can also be inferred
	    # automatically.
	    @parameters C_in=1. t0=1. a2b=1. a2c=1. b2c=1. c2b=1. C_out=1.
	    @species A(t) B(t) C(t)
	    # C_in, ∅ → A  # this is a birth rate, from an external source to A, or A is generated de novo
	    C_in*(sign(t-t0)+1.)/2., ∅ → A # same as line above, but C_in only starts after t0
	    a2b, A → B # transfer from A to B with speed a2b
	    a2c, A → C # transfer from A to C with speed a2c
	    (b2c, c2b), B ↔ C
	    C_out, C → ∅
	end
end

# ╔═╡ ca063487-960b-4cfc-a0f3-f15f8c7d0279
rs

# ╔═╡ 8149786f-5c8f-4150-a355-7eaf190e7239
Graph(rs)

# ╔═╡ c2834401-6e2d-4d57-8dab-4e57cc849b7f
complexgraph(rs)

# ╔═╡ 6f418969-14f1-4211-849d-a24e741b99ee
# Matrix-Vector Reaction Rate Equation Representation
netstoichmat(rs)

# ╔═╡ 77a50173-311f-4d33-8b88-611638c29f23
begin
	p = (
	    :C_in  => 1., 
	    :t0    => 12., 
	    :a2b   => 0.,
	    :a2c   => 0.7,
	    :b2c   => 0.4,
	    :c2b   => 0.1,
	    :C_out => 0.7
	) 
	tspan = (0., 30.)
	# start at concentration = 0
	u0 = [
	    :A => 0.5, 
	    :B => 0.5,
	    :C => 0.5
	]  
	
	# solve ODEs: the reaction network is transformed into an ODEsystem under the hood
	oprob = ODEProblem(rs, u0, tspan, p)
	osol  = solve(oprob, Tsit5())
	plot(osol)
end

# ╔═╡ 76a5d5a5-6e9d-4d83-9ca3-b144fc7a8b39


# ╔═╡ 6dc0a4f8-44fd-478f-88b7-e8786d4c8077
# get the species in the reaction network
species(rs)

# ╔═╡ cbc53770-b875-4124-9ebb-239f29d94a73
# ODEs satisfied by x(t) are then
# \frac{d\mathbf{x}}{dt} = N \mathbf{v}(\mathbf{x}(t),t),
# N is a constant M by K matrix with the net stoichiometric coefficient of species m in reaction k. 
N = netstoichmat(rs)

# ╔═╡ 9b290fdd-9df6-4e84-927d-c4a204d452e3
# the different reactions ν:
rxs = reactions(rs)

# ╔═╡ 3c2c6aaa-dc4a-41db-99aa-281c9f90163c
ν = oderatelaw.(rxs)

# ╔═╡ Cell order:
# ╠═4ad03253-636a-450c-863f-b7c8f4367bf5
# ╠═ad5fede2-e600-11ee-374e-85534f76407e
# ╠═077fef83-15c6-48cd-9e9f-4ee2644a0c2b
# ╠═ca063487-960b-4cfc-a0f3-f15f8c7d0279
# ╠═8149786f-5c8f-4150-a355-7eaf190e7239
# ╠═c2834401-6e2d-4d57-8dab-4e57cc849b7f
# ╠═6f418969-14f1-4211-849d-a24e741b99ee
# ╠═77a50173-311f-4d33-8b88-611638c29f23
# ╠═76a5d5a5-6e9d-4d83-9ca3-b144fc7a8b39
# ╠═6dc0a4f8-44fd-478f-88b7-e8786d4c8077
# ╠═cbc53770-b875-4124-9ebb-239f29d94a73
# ╠═9b290fdd-9df6-4e84-927d-c4a204d452e3
# ╠═3c2c6aaa-dc4a-41db-99aa-281c9f90163c
