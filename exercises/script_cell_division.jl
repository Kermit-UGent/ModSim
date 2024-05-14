using Catalyst
using DifferentialEquations, Plots

# Define all reactions
interaction_mpf_cycline = @reaction_network begin
    α*C, ∅ --> M
    β*C*M^2, ∅ --> M
    γ/(1+M), M --> ∅
    δ, ∅ --> C
    -M, ∅ --> C
end

species(interaction_mpf_cycline)
# M(t)
# C(t)

parameters(interaction_mpf_cycline)
# α
# β
# γ
# δ

# Convert to an ODE system
osys  = convert(ODESystem, interaction_mpf_cycline)

# Check out the diff. eqns. equations
equations(osys)
# Differential(t)(M(t)) ~ (-M(t)*γ) / (1 + M(t)) + C(t)*α + C(t)*(M(t)^2)*β
# Differential(t)(C(t)) ~ -M(t) + δ

# Set initial condition for ODE problem
u₀ = [:M => 2.0, :C => 0.1]
# Set time span for simulation
tspan = (0.0, 20.0)
# Set parameters for ODE problem
params = [:α => 2.0, :β => 1.0, :γ => 10.0, :δ => 1.0]
# Create ODE problem
oprob = ODEProblem(interaction_mpf_cycline, u₀, tspan, params)

osol = solve(oprob, Tsit5(), saveat=0.01)
plot(osol)