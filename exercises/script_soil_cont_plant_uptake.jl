
using Catalyst
using DifferentialEquations, Plots

#======================================================================================
  Solving using Catalyst
======================================================================================#

# Define all reactions
soil_cont_plant_uptake = @reaction_network begin
    Cin, ∅ --> C
    k₁, C --> ∅
    k₂, C + P --> 2P
    # k₂, C --> P            # graph is not so interesting
    k₃, P --> ∅
end

parameters(soil_cont_plant_uptake)
species(soil_cont_plant_uptake)

osys  = convert(ODESystem, soil_cont_plant_uptake)

# Check out the diff. eqns. equations
equations(osys)
# Differential(t)(C(t)) ~ Cin - k₁*C(t) - k₂*C(t)*P(t)
# Differential(t)(P(t)) ~ -k₃*P(t) + k₂*C(t)*P(t)

Cin = 0.06
k₁ = 4.1e-3    # 1e-3 to 1e-1 per day
k₂ = 1.9e-2    # 1e-3 to 1e-1 per day
k₃ = 2.2e-2    # 1e-3 to 1e-1 per day

C₀ = 0.001     # mg/kg
P₀ = 0.001

u₀ = [:C => C₀, :P => P₀]
# Set time span for simulation
tspan = (0.0, 400)
# Set parameters for ODE problem
params = [:Cin => Cin, :k₁ => k₁, :k₂ => k₂, :k₃ => k₃]
# Create ODE problem
# oprob = ODEProblem(soil_cont_plant_uptake, u₀, tspan, params, combinatoric_ratelaws=false)
oprob = ODEProblem(soil_cont_plant_uptake, u₀, tspan, params)

osol = solve(oprob)

plot(osol)