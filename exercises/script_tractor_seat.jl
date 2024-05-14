#======================================================================================
=======================================================================================
  Tractor seat
=======================================================================================
======================================================================================#

#=
https://docs.sciml.ai/Catalyst/stable/introduction_to_catalyst/catalyst_for_new_julia_users/
https://docs.sciml.ai/Catalyst/stable/catalyst_functionality/dsl_description/
=#

using Catalyst
using DifferentialEquations
using Plots


function heaviside(t, t0, v0, v1)
    (v1 - v0) .* (t .> t0) .+ v0
end

#======================================================================================
  Solving using Catalyst
======================================================================================#

Fin(t) = heaviside(t, 1, 0, -100.0)
rn = @reaction_network begin
    x2, 0 --> x1
    2*ζ*ω, x2 --> 0
    ω^2*x1/x2, x2 --> 0
    K*ω^2*Fin(t), 0 --> x2
end

osys  = convert(ODESystem, rn)
println(osys)
# Differential(t)(x1(t)) ~ x2(t)
# Differential(t)(x2(t)) ~ -x1(t)*(ω^2) - 2x2(t)*ζ*ω - 100.0*K*(t > 0)*(ω^2)

β = 920.0; k = 22600; m = 80.0
u0 = [:x1 => 0.0, :x2 => 0.0]
params = [:ζ => (1/2)*β/sqrt(k*m), :ω => sqrt(k/m), :K => 1/k]
tspan = (0.0, 3.0)  # the time interval to solve on
oprob = ODEProblem(rn, u0, tspan, params)
sol = solve(oprob, abstol=1e-9, reltol=1e-9)

# convert 92-element Vector{Vector{Float64}} into 92×2 Matrix{Float64}
solmat = mapreduce(permutedims, vcat, sol.u)

t = sol.t
x = solmat[:,1]

plot(t, x, linecolor=:blue, linewidth=2, xaxis="time [s]", yaxis="movement [m]", label="x")


#======================================================================================
  Solving using DifferentialEquations
======================================================================================#

Fin(t) = heaviside(t, 1, 0, -100.0)
function tractorseat!(du, u, p, t)
    ζ, ω, K = p
    du[1] = u[2]
    du[2] = -2*ζ*ω*u[2] - ω^2*u[1] + K*ω^2*Fin(t)
end

β = 920.0; k = 22600; m = 80.0
u0 = [0.0, 0.0]
params = ((1/2)*β/sqrt(k*m), sqrt(k/m), 1/k)
tspan = (0.0, 3.0)
oprob = ODEProblem(tractorseat!, u0, tspan, params)

sol = solve(oprob, abstol=1e-9, reltol=1e-9)

# convert 92-element Vector{Vector{Float64}} into 92×2 Matrix{Float64}
solmat = mapreduce(permutedims, vcat, sol.u)

t = sol.t
x = solmat[:,1]

plot(t, x, linecolor=:blue, linewidth=2, xaxis="time [s]", yaxis="movement [m]", label="x")


# Influence of the parameter ω


# Influence of the parameter ζ


# Influence of the parameter K


