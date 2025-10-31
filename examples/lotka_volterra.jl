### A Pluto.jl notebook ###
# v0.20.10

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 1d4c01b0-0ed4-4d7b-be00-722f43076f29
using Pkg; Pkg.activate("..")

# ╔═╡ ac1e6158-d41f-11ef-07a3-bbf7014dbbcd
using Catalyst, OrdinaryDiffEq, Plots

# ╔═╡ b29ecd80-0e11-4501-a3bd-fbb778630072
using PlutoUI

# ╔═╡ f1d5bd29-2142-415a-84e2-2ab4e21cd0db
osys = convert(ODESystem, lotka_volterra)

# ╔═╡ df247423-a290-440e-a87b-020584db645e
# ╠═╡ disabled = true
#=╠═╡
oprob = ODEProblem(lotka_volterra, [:🐰=>10.0, :🦊 => 1.0],
						(0.0, 50.0), [:α=>0.8, :β=>0.2, :γ=>0.6], combinatoric_ratelaws=false)
  ╠═╡ =#

# ╔═╡ a5acd96a-b759-4177-8097-a5f73ddda483
oprob = ODEProblem(complete(osys), [:🐰=>10.0, :🦊 => 1.0],
						(0.0, 50.0), [:α=>0.8, :β=>0.2, :γ=>0.6, :δ=>0.2])

# ╔═╡ 9bfd05e1-6a66-40e9-ad2d-0b2911095647
sol = solve(oprob, RK4())

# ╔═╡ 79212598-3858-43e5-af34-4059fea38246
plot(sol)

# ╔═╡ b051487d-12fc-4a38-b580-3379ddb1521d
plot(sol, idxs=(:🐰,:🦊), xlab="rabbits", ylab="foxes", arrow=true)

# ╔═╡ 10880546-4f88-430c-8453-8a18fec6b611
md"Can you adapt the model to include a carrying capacity of 20 for the rabbits?"

# ╔═╡ 5449da57-bd43-4fa4-8519-ff0c91cf8001
@bind showsol CheckBox()

# ╔═╡ f4e99c48-acfc-4acc-a70f-c4d4d33c2d72
showsol && md"""
```
lotka_volterra = @reaction_network begin
	α * (1 - 🐰/20), 🐰 --> 2🐰
	β, 🦊 --> ∅
	γ, 🐰 + 🦊 --> 🦊 + 0.1🦊
end
```
"""

# ╔═╡ 4dd925bf-d490-4b0a-8363-db376cbad00b
lotka_volterra = @reaction_network begin
	α, 🐰 --> 2🐰         # reproduction prey
	β, 🐰 + 🦊 --> 🦊     # predation
	γ, 🦊 --> ∅           # mortality predator
	δ * 🐰, 🦊  --> 2🦊   # predation growths
end

# ╔═╡ 2604b8c5-3ed7-4be5-bb5d-20ec5c2598ba
# ╠═╡ disabled = true
#=╠═╡
lotka_volterra = @reaction_network begin
	α, 🐰 --> 2🐰   
	β, 🐰 + 🦊 --> 🦊 + δ * 🦊
	γ, 🦊 --> ∅
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═1d4c01b0-0ed4-4d7b-be00-722f43076f29
# ╠═ac1e6158-d41f-11ef-07a3-bbf7014dbbcd
# ╠═b29ecd80-0e11-4501-a3bd-fbb778630072
# ╠═4dd925bf-d490-4b0a-8363-db376cbad00b
# ╠═2604b8c5-3ed7-4be5-bb5d-20ec5c2598ba
# ╠═f1d5bd29-2142-415a-84e2-2ab4e21cd0db
# ╠═df247423-a290-440e-a87b-020584db645e
# ╠═a5acd96a-b759-4177-8097-a5f73ddda483
# ╠═9bfd05e1-6a66-40e9-ad2d-0b2911095647
# ╠═79212598-3858-43e5-af34-4059fea38246
# ╠═b051487d-12fc-4a38-b580-3379ddb1521d
# ╟─10880546-4f88-430c-8453-8a18fec6b611
# ╟─5449da57-bd43-4fa4-8519-ff0c91cf8001
# ╟─f4e99c48-acfc-4acc-a70f-c4d4d33c2d72
