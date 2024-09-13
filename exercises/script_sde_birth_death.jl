
using Catalyst
using DifferentialEquations, Plots

bd_ode_process = @reaction_network begin
    (p,d), ∅ <--> X
end

u0 = [:X => 1]
tspan = (0.0, 200.0)
params = [:p => 0.05, :d => 0.01]

oprob = ODEProblem(bd_ode_process, u0, tspan, params)
osol = solve(oprob)
plot(osol)


# https://docs.sciml.ai/Catalyst/stable/catalyst_applications/advanced_simulations/

bd_sde_process = @reaction_network begin
    @parameters η1 η2
    # (p,d), ∅ <--> X
    k1, 0 --> X
    k2, X --> 0
end

u0 = [:X => 1]
tspan = (0.0, 200.0)
params = [:k1 => 0.05, :k2 => 0.01, :η1 => 0.2, :η2 => 0.2]
sprob = SDEProblem(bd_sde_process, u0, tspan, params; noise_scaling = @parameters η1 η2)
# ssol = solve(sprob, STrapezoid())
ssol = solve(sprob, LambaEM())
# ssol = solve(sprob, EM(), dt=0.1)
plot(ssol)

