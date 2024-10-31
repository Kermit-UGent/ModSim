#=
Stochastic Differential Equations (SDEs) are mathematical equations used to model
systems influenced by random noise. They extend ordinary differential equations (ODEs)
by incorporating terms that represent stochastic processes, typically in the form of
a Wiener process or Brownian motion. SDEs are widely used in various fields, such as
physics, biology, finance, and engineering, to describe the evolution of systems
under uncertainty or with inherent randomness.
=#

using Catalyst
using DifferentialEquations, Plots

infection_sde_model = @reaction_network begin
	@parameters η=40
	@default_noise_scaling η
	α * β, S + I --> 2I, [noise_scaling = 60]
	r * m, I --> D
	r * (1 - m), I --> R
end

species(infection_sde_model)
parameters(infection_sde_model)

osys = convert(ODESystem, infection_sde_model)
equations(osys)
# Differential(t)(S(t)) ~ -S(t)*I(t)*α*β
# Differential(t)(I(t)) ~ (-1 + m)*r*I(t) - m*r*I(t) + S(t)*I(t)*α*β
# Differential(t)(D(t)) ~ m*r*I(t)
# Differential(t)(R(t)) ~ (1 - m)*r*I(t)

u0 = [:S => 9999000, :I => 1000, :D => 0, :R => 0]
tspan = (0.0, 90.0)
params = [:α => 0.08, :β => 1.0e-6, :r => 0.2, :m => 0.4]


# https://docs.sciml.ai/DiffEqDocs/stable/tutorials/sde_example/
# If we allow for some randomness in some of the coefficients of a differential
# equation we often obtain a more realistic mathematical model of the situation.

# https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/
# https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/#Full-List-of-Methods

sprob = SDEProblem(infection_sde_model, u0, tspan, params)

# EM - The Euler-Maruyama method. Fixed time step only.
# ssol = solve(sprob, EM(), dt=0.1)

# LambaEM - A modified Euler-Maruyama method with adaptive time stepping
ssol = solve(sprob, EM(), dt=0.1)
# ssol = solve(sprob, LambaEM(), tstops=range(0.0, step=0.1, length=901))

plot(ssol)


# https://stackoverflow.com/questions/69049991/simulating-a-reflecting-boundary-sdeproblem
function condition(u,t,integrator)
	true
end
function affect!(integrator)
	for i = 1:4
		if integrator.u[i] > 10000000
			integrator.u[i] = 10000000
		end
		if integrator.u[i] < 0
			integrator.u[i] = 0
		end
	end
end
# https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/#discrete_callback
cb = DiscreteCallback(condition,affect!;save_positions=(false,true))
# save_positions: Boolean tuple for whether to save before and after the affect!.
#     This saving will occur just before and after the event, only at event times,
#     and does not depend on options like saveat,
ssol = solve(sprob, EM(), dt=0.1, callback=cb, save_everystep=true)
plot(ssol)

esprob = EnsembleProblem(sprob)
essol = solve(esprob, EM(), dt=0.1, callback=cb, save_everystep=true, trajectories=100)
plot(essol)



# dprob = DiscreteProblem(infection_model, u0, tspan, params)
# jprob = JumpProblem(infection_model, dprob, Direct())
# dsol = solve(jprob, SSAStepper())
# plot(dsol)

# ensemblejprob = EnsembleProblem(jprob)
# edsol = solve(ensemblejprob, SSAStepper(), EnsembleThreads(), trajectories = 2)
# plot(edsol)