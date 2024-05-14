#======================================================================================
=======================================================================================
  Sludge system - nitrification of ammonia
=======================================================================================
======================================================================================#

#=
https://docs.sciml.ai/Catalyst/stable/introduction_to_catalyst/catalyst_for_new_julia_users/
https://docs.sciml.ai/Catalyst/stable/catalyst_functionality/dsl_description/
https://www.matecdev.com/posts/julia-plotting-multiple-plots.html
=#

using Catalyst
using DifferentialEquations
using Plots
# using Symbolics

#======================================================================================
  Define the step function
======================================================================================#

# Step function, go from v0 to v1 at t = t0
function heaviside(t, t0, v0, v1)
    (v1 - v0) .* (t .> t0) .+ v0
end

#======================================================================================
  Solving using Catalyst
======================================================================================#

A₀(t) = heaviside(t, 0, 0.0, 0.03)

# Define all reactions
sludge_system = @reaction_network begin
    F₀/V₁*A₀(t), 0 --> A₁
    μmax*N₁/(K+A₁), A₁ --> Y*N₁
    F₀*(R+1)/V₁, A₁ --> (V₁/V₂)*A₂
    F₀*(R+1)/V₁, N₁ --> (V₁/V₂)*N₂
    F₀*(R+1)/C/V₂, N₂ --> (R/(R+1)*V₂/V₁*C)*N₁
    F₀*(R+1)/C/V₂, A₂ --> (R/(R+1)*V₂/V₁*C)*A₁
    F₀*(R+1)*(1-1/C)/V₂, A₂ --> 0
end

# Convert to an ODE system
osys  = convert(ODESystem, sludge_system, combinatoric_ratelaws=false)

# Check out the diff. eqns. equations
equations(osys)
# 4-element Vector{Equation}:
#  Differential(t)(A₁(t)) ~ (-F₀*(1 + R)*A₁(t)) / V₁ + (-N₁(t)*A₁(t)*μmax) / (K + A₁(t)) + (A₀*F₀) / V₁ + (F₀*R*A₂(t)) / V₁
#  Differential(t)(N₁(t)) ~ (F₀*R*N₂(t)) / V₁ + (Y*N₁(t)*A₁(t)*μmax) / (K + A₁(t)) + (-F₀*(1 + R)*N₁(t)) / V₁
#  Differential(t)(A₂(t)) ~ (-F₀*(1 + R)*A₂(t)*(1 + -1 / C)) / V₂ + (F₀*(1 + R)*A₁(t)) / V₂ + (-F₀*(1 + R)*A₂(t)) / (C*V₂)
#  Differential(t)(N₂(t)) ~ (F₀*(1 + R)*N₁(t)) / V₂ + (-F₀*(1 + R)*N₂(t)) / (C*V₂)
# Check out the states
states(osys)
# 4-element Vector{SymbolicUtils.BasicSymbolic{Real}}:
#  A₁(t)
#  N₁(t)
#  A₂(t)
#  N₂(t)
# Check out the parameters
parameters(osys)
# 9-element Vector{SymbolicUtils.BasicSymbolic{Real}}:
#  F₀
#  V₁
#  A₀
#  μmax
#  K
#  Y
#  R
#  V₂
#  C

R = 0.20     # recycle factor
C = 1.30     # concentration factor

Y  = 0.70    # yield coefficient [kg/kg]

V₁ = 400.0   # volume reactor [m^3]
V₂ = 400.0   # volume settler [m^3]

μmax = 0.08  # Maximum growth rate of the nitrifying organisms [1/h]
K = 0.02      # Saturation constant of nitrifying organisms [kg/m^3]

# A0 = 0.03;
# r1 = 0
# r2 = 0.03
# tswitch = 0.0
F₀ = 5

# Initial conditions
A₁₀ = 0.0     # Initial concentration of ammonia in reactor (1) [kg/m^3]
N₁₀ = 0.01    # Initial concentration of nitrifiers in reactor (1) [kg/m^3]
A₂₀ = 0.0     # Initial concentration of ammonia in settler (2) [kg/m^3]
N₂₀ = 0.02    # Initial concentration of nitrifiers in settler (2) [kg/m^3]

u₀ = [:A₁ => A₁₀, :N₁ => N₁₀, :A₂ => A₂₀, :N₂ => N₂₀]
params = [:R => R, :C => C, :Y => Y, :V₁ => V₁, :V₂ => V₂, :μmax => μmax, :K => K, :F₀ => F₀]
tspan = (0.0, 600.0)

oprob = ODEProblem(sludge_system, u₀, tspan, params, combinatoric_ratelaws=false)

sol = solve(oprob, Tsit5(), saveat=1)

plot(sol)
