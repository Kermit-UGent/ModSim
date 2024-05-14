#======================================================================================
=======================================================================================
  Fermenter with monod kinetics
=======================================================================================
======================================================================================#

using Catalyst
using DifferentialEquations, Plots

#======================================================================================
  Solving using Catalyst
======================================================================================#

# f1(S, X, Q) = (Sin-S)*Q/V-1/Y*mu_max*S/(Ks+S)*X
# f2(S, X, Q) = -X*Q/V+mu_max*S/(Ks+S)*X

# Define all reactions
fermenter_monod = @reaction_network begin
    X * mm(S, μmax, Ks), S --> Y*X     # Y*X is created from one S at a rate X * mm(S, μmax, Ks)
    Q/V, (S, X) --> 0                  # S and X are degraded at a rate Q/V*S
    Q/V*Sin, 0 --> S                   # S is created at a rate Q/V*Sin 
end

parameters(fermenter_monod)
# μmax
# Ks
# Y
# Q
# V
# Sin

# Convert to an ODE system
osys  = convert(ODESystem, fermenter_monod)

# Check out the diff. eqns. equations
equations(osys)
#  Differential(t)(S(t)) ~ (Q*Sin) / V + (-Q*S(t)) / V - X(t)*S(t)*Catalyst.mm(S(t), μmax, Ks)
#  Differential(t)(X(t)) ~ (-Q*X(t)) / V + Y*X(t)*S(t)*Catalyst.mm(S(t), μmax, Ks)

# Check out the states
states(osys)
#  S(t)
#  X(t)

# Check out the parameters
parameters(osys)
#  μmax
#  Ks
#  Y
#  Q
#  V
#  Sin

S₀ = 0.02
X₀ = 5.0e-5

μmax = 0.39
Ks = 0.48
Y = 0.77
Q = 2
V = 40
Sin = 2.2

# Set initial condition for ODE problem
u₀ = [:S => S₀, :X => X₀]
# Set time span for simulation
tspan = (0.0, 100)
# Set parameters for ODE problem
params = [:μmax => μmax, :Ks => Ks, :Y => Y, :Q => Q, :V => V, :Sin => Sin]
# Create ODE problem
oprob = ODEProblem(fermenter_monod, u₀, tspan, params)

osol = solve(oprob, Tsit5(), saveat=0.5)
plot(osol)


################################################################
## Sensitivity analysis
################################################################

using SciMLSensitivity

oprob_sens = ODEForwardSensitivityProblem(oprob.f, [S₀, X₀], tspan, [μmax, Ks, Y, Q, V, Sin])
# oprob.f        : the ODE function
# [0.02, 5.0e-5] : initial conditions, the order is [S, X], see species(fermenter_monod)
# tspan          : the time span
# [0.39, 0.48, 0.77, 2, 40, 2.2] : the parameter values, the order is [μmax, Ks, Y, Q, V Sin], see parameters(fermenter_monod)

osol_sens = solve(oprob_sens, Tsit5(), saveat=0.5)
u, dp = extract_local_sensitivities(osol_sens)

sol_S = u[1,:]     # solution for S
sol_X = u[2,:]     # solution for X

sens_μmax = dp[1]  # sensitivity for S and X on μmax  
sens_Ks   = dp[2]  # sensitivity for S and X on Ks
sens_Sin  = dp[6]  # sensitivity for S and X on Sin


# Absolute sensitivity of S and X to μmax
plot(osol_sens.t, sens_μmax', title="Absolute sensitivity of S and X to μmax", label=["dS/dμmax" "dX/dμmax"])

# Absolute sensitivity of S and X to Ks
plot(osol_sens.t, sens_Ks', title="Absolute sensitivity of S and X to Ks", label=["dS/dKs" "dX/dKs"])

# Absolute sensitivity of S and X to Sin
plot(osol_sens.t, sens_Sin', title="Absolute sensitivity for S and X on Ks", label=["dS/dSin" "dX/dSin"])

# Absolute sensitivity of S to μmax and Ks and Sin
plot(osol_sens.t, sens_μmax[1,:], title="Absolute sensitivity of S and X to Sin", label="dS/dμmax")
plot!(osol_sens.t, sens_Ks[1,:], label="dS/dKs")
plot!(osol_sens.t, sens_Sin[1,:], label="dS/dSin")

# Absolute sensitivity of X to μmax and Ks and Sin
plot(osol_sens.t, sens_μmax[2,:], title="Absolute sensitivity of X to μmax and Ks and Sin", label="dX/dμmax")
plot!(osol_sens.t, sens_Ks[2,:], label="dX/dKs")
plot!(osol_sens.t, sens_Sin[2,:], label="dX/dSin")


# Relative sensitivity of S to μmax and Ks and Sin
plot(osol_sens.t, sens_μmax[1,:].*μmax, title="Relative sensitivity of S to μmax and Ks and Sin", label="dS/dμmax * μmax")
plot!(osol_sens.t, sens_Ks[1,:].*Ks, label="dS/dKs * Ks")
plot!(osol_sens.t, sens_Sin[1,:].*Sin, label="dS/dSin * Sin")

# Relative sensitivity of X to μmax and Ks and Sin
plot(osol_sens.t, sens_μmax[2,:].*μmax, title="Relative sensitivity for X on μmax and Ks", label="dX/dμmax * μmax")
plot!(osol_sens.t, sens_Ks[2,:].*Ks, label="dX/dKs * Ks")
plot!(osol_sens.t, sens_Sin[2,:].*Sin, label="dX/dSin * Sin")


# Total relative sensitivity of S and X to μmax and Ks and Sin
plot(osol_sens.t, sens_μmax[1,:].*μmax./sol_S, title="Total relative sensitivity for S and X on μmax and Ks", label="dS/dμmax * μmax / S")
plot!(osol_sens.t, sens_Ks[1,:].*Ks./sol_S, label="dS/dKs * Ks / S")
plot!(osol_sens.t, sens_Sin[1,:].*Sin./sol_S, label="dS/dSin * Sin / S")
plot!(osol_sens.t, sens_μmax[2,:].*μmax./sol_X, label="dX/dμmax * μmax / X")
plot!(osol_sens.t, sens_Ks[2,:].*Ks./sol_X, label="dX/dKs * Ks / X")
plot!(osol_sens.t, sens_Sin[2,:].*Sin./sol_X, label="dX/dSin * Sin / X")


sum( sqrt.( sum(( sens_μmax'./u'.*0.39 ).^2, dims=2) ./ 2 ) ) / 201
sum( sqrt.( sum(( sens_Ks'./u'.*0.48 ).^2, dims=2) ./ 2 ) ) / 201


################################################################
## Calibration - parameter estimation
################################################################

using Distributions
x = rand(Truncated(Normal(0, 1), -3, 3), 10)
x = rand(Normal(0, 1), 10)

u₀ = [:S => 0.02, :X => 5.0e-5]
tspan = (0.0, 100)
params = [:μmax => 0.39, :Ks => 0.48, :Y => 0.77, :Q => 2, :V => 40, :Sin => 2.2]
oprob = ODEProblem(fermenter_monod, u₀, tspan, params)
osol = solve(oprob, Tsit5(), saveat=0.5)
# plot(osol)
size(osol, 2)

Ssim = osol[:S]
Xsim = osol[:X]

maximum(Ssim)

sigma_S = 0.14
sigma_X = 0.057

ratio = 0.8
Sdta = Ssim + rand(Normal(0, sigma_S), size(osol, 2)) .* ((1 - ratio) .+ ratio.*Ssim./maximum(Ssim))
Xdta = Xsim + rand(Normal(0, sigma_X), size(osol, 2)) .* ((1 - ratio) .+ ratio.*Xsim./maximum(Xsim))

scatter(osol.t, Sdta)
scatter!(osol.t, Xdta)
