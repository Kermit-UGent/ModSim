
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

u0 = [:X => 2.0]
tspan = (0.0, 520.0)
params = [:b => 3.0, :d => 0.015]

oprob = ODEProblem(birth_death, u0, tspan, params)
osol = solve(oprob, Tsit5(), saveat=1.0)

plot(osol)


#=
Suppose that at t = 250 the death rate of the mice population increases by 30 % due
to a new predator species in the area.
=#

condition2 = [250.0]
function affect2!(integrator)
    integrator.p[2] *= (1 + 30/100)      # Qin is the 3rd parameter !!!
end
cb2 = PresetTimeCallback(condition2, affect2!)
osol2 = solve(deepcopy(oprob), Tsit5(), saveat=1.0, callback=cb2)
plot(osol2)


#=
In overcrowded environments, female mice may experience decreased reproductive success
due to factors such as limited access to nesting sites and increased competition for mates.
Suppose that when the mice population reaches 150, the birth rate decreases by 20 %.
=#

function condition2(u, t, integrator)
    u[1] - 150
end
function affect!(integrator)
    integrator.p[1] *= (1 - 20/100)
end
# function affect!(integrator)
#     integrator.p[1] *= (1 - 20/100)
#     integrator.p[2] *= (1 + 20/100)
# end
ps_cb = ContinuousCallback(condition2, affect!)
osol = solve(deepcopy(oprob), Tsit5(), saveat=1.0, callback=ps_cb)
plot(osol)
