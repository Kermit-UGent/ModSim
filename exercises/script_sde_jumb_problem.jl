#=
https://docs.juliahub.com/DiffEqJump/Lhfzo/8.6.3/tutorials/discrete_stochastic_example/

=#

using Catalyst
using DifferentialEquations, Plots

infection_model = @reaction_network begin
	α * β, S + I --> 2I
	r * m, I --> D
	r * (1 - m), I --> R
end

species(infection_model)

u0 = [:S => 50, :I => 1, :D => 0, :R => 0]
tspan = (0.0, 60.0)
params = [:α => 0.15, :β => 0.1, :r => 0.2, :m => 0.6]

# https://docs.sciml.ai/DiffEqDocs/stable/tutorials/sde_example/
# If we allow for some randomness in some of the coefficients of a differential
# equation we often obtain a more realistic mathematical model of the situation.
sprob = SDEProblem(infection_model, u0, tspan, params)
ssol = solve(sprob, EM(), dt=0.1)
plot(ssol)
# https://stackoverflow.com/questions/69049991/simulating-a-reflecting-boundary-sdeproblem
function condition(u,t,integrator)
	true
end
function affect!(integrator)
	for i = 1:4
		if integrator.u[i] > 50
			integrator.u[i] = 50
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
ssol = solve(sprob, EM(), dt=0.1, callback=cb, save_everystep=false)
plot(ssol)


dprob = DiscreteProblem(infection_model, u0, tspan, params)
jprob = JumpProblem(infection_model, dprob, Direct())
dsol = solve(jprob, SSAStepper())
plot(dsol)

# ensemblejprob = EnsembleProblem(jprob)
# edsol = solve(ensemblejprob, SSAStepper(), EnsembleThreads(), trajectories = 2)
# plot(edsol)

#=

https://docs.juliahub.com/DiffEqJump/Lhfzo/8.6.3/tutorials/discrete_stochastic_example/

Imagine a bike share system for students traveling between Olin College and Wellesley
College, which are about 5 km apart in eastern Massachusetts.
Suppose the system contains 12 bikes and 2 bike racks, one at Olin and one at Wellesley,
each with the capacity to hold 12 bikes. As students arrive, check out a bike, and
ride to the other campus, the number of bikes in each location changes.

Here Direct() indicates that we will determine the random times and types of reactions
using Gillespie's Direct stochastic simulation algorithm (SSA), also known as Doob's
method or Kinetic Monte Carlo. 

We will use SSAStepper() to handle time-stepping the Direct method from jump to jump.

=#

bikeshare_model = @reaction_network begin
	p1, Bo --> Bw
	p2, Bw --> Bo
end

osys  = convert(ODESystem, bikeshare_model)
equations(osys)

u0 = [:Bo => 10, :Bw => 2]
tspan = (0.0, 30.0)
params = [:p1 => 0.50, :p2 => 0.20]

dprob = DiscreteProblem(bikesharing_model, u0, tspan, params)
jprob = JumpProblem(bikesharing_model, dprob, Direct())

djsol = solve(jprob, SSAStepper())
plot(djsol, ylims=(0, 12))
# dsol.t[end-5:end]


bikesharing_model = @reaction_network begin
	p1, Bo --> Tow
	r, Tow --> Bw
	p2, Bw --> Two
	r, Two --> Bo
end

species(bikesharing_model)

osys  = convert(ODESystem, bikesharing_model)
equations(osys)

u0 = [:Bo => 12, :Bw => 0, :Tow => 0, :Two => 0]
tspan = (0.0, 30.0)
params = [:p1 => 0.5, :p2 => 0.50, :r => 0.1]

dprob = DiscreteProblem(bikesharing_model, u0, tspan, params)
jprob = JumpProblem(bikesharing_model, dprob, Direct())

dsol = solve(jprob, SSAStepper())
# dsol[:Bo]
# dsol[:Bw]
# plot(dsol, ylims=(0, 12))
begin
	plot(dsol.t, dsol[:Bo], ylims=(0, 12))
	plot!(dsol.t, dsol[:Bw], ylims=(0, 12))
end
