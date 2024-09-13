#======================================================================================
=======================================================================================
  Fermenter with first-order kinetics
=======================================================================================
======================================================================================#

#=
https://docs.sciml.ai/Catalyst/stable/introduction_to_catalyst/catalyst_for_new_julia_users/
https://docs.sciml.ai/Catalyst/stable/catalyst_functionality/dsl_description/
https://www.matecdev.com/posts/julia-plotting-multiple-plots.html
=#

using Catalyst
using DifferentialEquations, Plots


#======================================================================================
  ODE problem
======================================================================================#

# Define all reactions
fermenter_secondorder = @reaction_network begin
    k*X, S --> Y*X            # Y*X is created from one S at a rate β
    Q/V*Sin, 0 --> S       # S is created at a rate Q/V*Sin(t)
    Q/V, (S, X) --> 0                  # S degraded at a rate Q/V*S, X degraded at a rate Q/V*X
end

parameters(fermenter_secondorder)
# k
# Y
# Q
# V
# Sin

# Convert to an ODE system
osys = convert(ODESystem, fermenter_secondorder)

# Check out the diff. eqns. equations
equations(osys)
# Differential(t)(S(t)) ~ (Q*Sin) / V + (-Q*S(t)) / V - k*X(t)*S(t)
# Differential(t)(X(t)) ~ (-Q*X(t)) / V + Y*k*X(t)*S(t)

# Check out the states
states(osys)
#  S(t)
#  X(t)
# Check out the parameters
parameters(osys)
# k
# Y
# Q
# V
# Sin



# Set initial condition for ODE problem
u0 = [:S => 0.0, :X => 0.1]
# Set time span for simulation
tspan = (0.0, 120.0)  # the time interval to solve on
# Set parameters for ODE problem
params = [:k => 0.2, :Y => 0.76, :Q => 2, :V => 40, :Sin => 2.2]
# Create ODE problem
oprob = ODEProblem(fermenter_secondorder, u0, tspan, params, combinatoric_ratelaws=false)
osol = solve(oprob, Tsit5(), saveat=0.5)

plot(osol)

# #  Finding steady state
# osol[:S][end]
# osol[:X][end]
# u_guess = [:S => osol[:S][end], :X => osol[:X][end]]
# Seq, Xeq = solve(SteadyStateProblem(ODEProblem(fermenter_secondorder, u_guess, tspan, params, combinatoric_ratelaws=false)))
# # 0.36549707602339176
# # 1.3942222222222225


#======================================================================================
  SDE problem
======================================================================================#

fermenter_sde_secondorder = @reaction_network begin
  @parameters η1 η2 η3 η4
  k*X, S --> Y*X            # Y*X is created from one S at a rate β
  Q/V*Sin, 0 --> S       # S is created at a rate Q/V*Sin(t)
  Q/V, (S, X) --> 0                  # S degraded at a rate Q/V*S, X degraded at a rate Q/V*X
end

osys = convert(ODESystem, fermenter_sde_secondorder)
equations(osys)
# Differential(t)(S(t)) ~ (Q*Sin) / V + (-Q*S(t)) / V - k*X(t)*S(t)
# Differential(t)(X(t)) ~ (-Q*X(t)) / V + Y*k*X(t)*S(t)

# Set initial condition for ODE problem
u0 = [:S => 0.0, :X => 0.1]
# Set time span for simulation
tspan = (0.0, 120.0)  # the time interval to solve on
# Set parameters for ODE problem
params = [:k => 0.2, :Y => 0.76, :Q => 2, :V => 40, :Sin => 2.2, :η1 => 0.2, :η2 => 0.1, :η3 => 0.0, :η4 => 0.0]
# Create ODE problem
sprob = SDEProblem(fermenter_sde_secondorder, u0, tspan, params, combinatoric_ratelaws=false; noise_scaling = @parameters η1 η2 η3 η4)
ssol = solve(sprob, EM(), dt=0.1)
plot(ssol)

esprob = EnsembleProblem(sprob)
essol = solve(esprob, EM(), dt=0.1; trajectories=100)
plot(essol, ylim=(0.0,2.0); linealpha = 0.5)
