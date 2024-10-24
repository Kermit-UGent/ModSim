
# https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/#callbacks

using Catalyst
using DifferentialEquations, Plots

#=
In a simple birth-death model for mice, the birth rate of mice represents the
rate at which new individuals are added to the population through reproduction.
This rate is influenced by factors such as the number of reproductive females,
their fecundity, and the frequency of reproduction cycles. Conversely,
the death rate reflects the rate at which individuals are removed from
the population due to mortality factors such as predation, disease,
and environmental stressors. Together, these rates interact dynamically
to shape the population dynamics of mice in their natural habitat.
Denote the number of mice by X, de birth rate by b, and the death rate by d.
=#

birth_death = @reaction_network begin
    b, 0 --> X
    d, X --> 0
end

species(birth_death)
# X(t)
parameters(birth_death)
# b
# d
reactions(birth_death)
# b, ∅ --> X
# d, X --> ∅

osys  = convert(ODESystem, birth_death)   # convert rn to a ODESystem
equations(osys)
# Differential(t)(X(t)) ~ b - d*X(t)
unknowns(osys)
# X(t)
parameters(osys)
# b
# d

u0 = [:X => 2.0]

tspan = (0.0, 520.0)
params = [:b => 3.0, :d => 0.015]

dprob = DiscreteProblem(birth_death, u0, tspan, params)

jdprob = JumpProblem(birth_death, dprob, Direct())

jdsol = solve(jdprob, SSAStepper())

plot(jdsol)

# Events in time or states don't seem to work with Jump processes... (SSA)