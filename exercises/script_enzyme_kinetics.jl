using Catalyst
using DifferentialEquations, Plots

# https://en.wikipedia.org/wiki/Michaelis%E2%80%93Menten_kinetics

#======================================================================================
  Solving using Catalyst
======================================================================================#

# Define all reactions
enzyme_kin = @reaction_network begin
    k₁, S + E --> ES
    k₂, ES --> S + E
    k₃, ES --> P + E
end

parameters(enzyme_kin)
# k₁
# k₂
# k₃
species(enzyme_kin)
# S(t)
# E(t)
# ES(t)
# P(t)

osys  = convert(ODESystem, enzyme_kin)

# Check out the diff. eqns. equations
equations(osys)
# Differential(t)(S(t)) ~ k₂*ES(t) - k₁*E(t)*S(t)
# Differential(t)(E(t)) ~ k₂*ES(t) + k₃*ES(t) - k₁*E(t)*S(t)
# Differential(t)(ES(t)) ~ -k₂*ES(t) - k₃*ES(t) + k₁*E(t)*S(t)
# Differential(t)(P(t)) ~ k₃*ES(t)

# Set initial condition for ODE problem
u₀ = [:S => 0.05, :E => 0.02, :ES => 0.0, :P => 0.0]
# Set time span for simulation
tspan = (0.0, 500.0)
# Set parameters for ODE problem
params = [:k₁ => 1.2, :k₂ => 1.4, :k₃ => 0.08]
# params = [:k₁ => 0.41, :k₂ => 0.4, :k₃ => 0.39, :Km => 0.5]
# Create ODE problem
# oprob = ODEProblem(enzyme_kin, u₀, tspan, params, combinatoric_ratelaws=false)
oprob = ODEProblem(enzyme_kin, u₀, tspan, params)

# u_guess = [:S => 0.0, :E => 0.18, :ES => 0.0, :P => 0.05]
# Seq, Eeq, ESeq, Peq = solve(SteadyStateProblem(ODEProblem(enzyme_kin, u_guess, tspan, params, combinatoric_ratelaws=false)))

osol = solve(oprob, Tsit5(), saveat=0.1)
plot(osol)

osol[end]