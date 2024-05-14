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
  Solving using Catalyst
======================================================================================#

# Define all reactions
fermenter_firstorder = @reaction_network begin
    β, S --> Y*X            # Y*X is created from one S at a rate β
    Q/V*Sin, 0 --> S       # S is created at a rate Q/V*Sin(t)
    Q/V, (S, X) --> 0                  # S degraded at a rate Q/V*S, X degraded at a rate Q/V*X
end

parameters(fermenter_firstorder)
# β
# Y
# Q
# V
# Sin

# Convert to an ODE system
osys = convert(ODESystem, fermenter_firstorder)

# Check out the diff. eqns. equations
equations(osys)
# Differential(t)(S(t)) ~ (Q*Sin) / V + (-Q*S(t)) / V - S(t)*β
# Differential(t)(X(t)) ~ (-Q*X(t)) / V + Y*S(t)*β
# Check out the states
states(osys)
#  S(t)
#  X(t)
# Check out the parameters
parameters(osys)
# β
# Y
# Q
# V
# Sin



# Set initial condition for ODE problem
u0 = [:S => 0.0, :X => 0.1]
# Set time span for simulation
tspan = (0.0, 120.0)  # the time interval to solve on
# Set parameters for ODE problem
params = [:β => 0.98, :Y => 0.80, :Q => 2, :V => 40, :Sin => 2.2]
# Create ODE problem
oprob = ODEProblem(fermenter_firstorder, u0, tspan, params, combinatoric_ratelaws=false)
osol = solve(oprob, Tsit5(), saveat=0.5)

plot(osol)

#  Finding steady state
osol[:S][end]
osol[:X][end]
u_guess = [:S => osol[:S][end], :X => osol[:X][end]]
Seq, Xeq = solve(SteadyStateProblem(ODEProblem(fermenter_firstorder, u_guess, tspan, params, combinatoric_ratelaws=false)))
# 0.10679611650485438
# 1.6745631067961169

# Use steady state conditions as initial conditions
u0 = [:S => Seq, :X => Xeq]
# Create ODE problem
oprob = ODEProblem(fermenter_firstorder, u0, tspan, params, combinatoric_ratelaws=false)

condition = [20.0]
function affect!(integrator)
    integrator.p[5] = 3.4    # Sin is the 5th parameter !!!
end      
cb = PresetTimeCallback(condition, affect!)
osol = solve(deepcopy(oprob), Tsit5(), saveat=0.5, callback=cb)
plot(osol)

#  Finding new steady state
osol[:S][end]
osol[:X][end]
u_guess = [:S => osol[:S][end], :X => osol[:X][end]]
params = [:β => 0.98, :Y => 0.80, :Q => 2, :V => 40, :Sin => 3.4]
Seq, Xeq = solve(SteadyStateProblem(ODEProblem(fermenter_firstorder, u_guess, tspan, params, combinatoric_ratelaws=false)))
# 0.16504854368932037
# 2.5879611650485437
