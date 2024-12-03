### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 8e495fb4-b14e-11ef-005d-d3d50fec748d
begin
	import Pkg
	Pkg.activate("..")
end

# ╔═╡ bfed4c8f-2bbe-4a52-a426-6a66a0bc3e32
using Catalyst, Plots, DifferentialEquations

# ╔═╡ 4bc0f5d9-24d1-4149-bec0-0e4b504a99cc
md"""# Toilets

Five toilets are connected to a shared water reservoir of 100L. Each toilet's tank slowly fills from the big reservoir $T$, with a rate of $rT$. When a toilet's tank reaches a capacity of 10 liters, it empties its tank, and its content is transferred to $T$. Create a simulation of this system, starting with random initial water levels in each toilet tank. Run the simulation and observe the flushing patterns of the toilets over time. Does a pattern emerge in the way the toilets flush?
"""

# ╔═╡ 61f59d02-2e70-4e68-a328-4a090967e20d
toilets = @reaction_network begin
	r, T --> (V1, V2, V3, V4, V5)
end

# ╔═╡ 17fb9569-021d-445a-a0c7-6e1d9fca97a8
@unpack V1, V2, V3, V4, V5, T = toilets;

# ╔═╡ 2e92b4d3-3040-4b4f-a1b7-d90abbe19ab8
v = 10

# ╔═╡ 4773ce79-6a4a-42a6-898b-84bf118f2f40
r = 0.2

# ╔═╡ 0b3f9241-daf8-418c-920f-3bc3ae18bf70
flushing = [[V1 ~ v] => [T~T+v, V1~0],
			[V2 ~ v] => [T~T+v, V2~0],
			[V3 ~ v] => [T~T+v, V3~0],
			[V4 ~ v] => [T~T+v, V4~0],
			[V5 ~ v] => [T~T+v, V5~0],]

# ╔═╡ 0fe8f520-6545-496c-88ca-f756a7ce5e99
V0 = 10rand(5)

# ╔═╡ d8f1990c-2830-4785-96f0-202409c54d0d
T0 = 50

# ╔═╡ 5e0c47c4-f157-4371-a81d-294786980318
begin 
	@named toilets_flush = ReactionSystem(equations(toilets), continuous_events=flushing)
	toilets_flush = complete(toilets_flush)
end

# ╔═╡ d0adedb4-9045-43de-abc1-90a89dc55cb7
oprob = ODEProblem(toilets_flush, [T=>T0, V1=>V0[1], V2=>V0[2],V3=>V0[3],V4=>V0[4],V5=>V0[5],], (0.0, 10.0), [:r=>r])

# ╔═╡ 1572af6e-f2e3-4a45-9e78-359f1d976c31
sol = solve(oprob, Tsit5())

# ╔═╡ 360f4aa0-6a93-4ecb-9dfa-a43e3c9954c6
plot(sol, lw=2)

# ╔═╡ Cell order:
# ╠═bfed4c8f-2bbe-4a52-a426-6a66a0bc3e32
# ╠═8e495fb4-b14e-11ef-005d-d3d50fec748d
# ╟─4bc0f5d9-24d1-4149-bec0-0e4b504a99cc
# ╠═61f59d02-2e70-4e68-a328-4a090967e20d
# ╠═17fb9569-021d-445a-a0c7-6e1d9fca97a8
# ╠═2e92b4d3-3040-4b4f-a1b7-d90abbe19ab8
# ╠═4773ce79-6a4a-42a6-898b-84bf118f2f40
# ╠═0b3f9241-daf8-418c-920f-3bc3ae18bf70
# ╠═0fe8f520-6545-496c-88ca-f756a7ce5e99
# ╠═d8f1990c-2830-4785-96f0-202409c54d0d
# ╠═5e0c47c4-f157-4371-a81d-294786980318
# ╠═d0adedb4-9045-43de-abc1-90a89dc55cb7
# ╠═1572af6e-f2e3-4a45-9e78-359f1d976c31
# ╠═360f4aa0-6a93-4ecb-9dfa-a43e3c9954c6
