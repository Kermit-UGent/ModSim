using Catalyst


chemostat = @reaction_network Chemostat begin
    X * mm(S, μmax, Ks), S => 0.8X   # growth
    c * q / V, ∅ --> S               # inflow substrate
    q / V, (S, X) --> ∅              # outflow
end

mm

species(chemostat)
#  S(t)
#  X(t)

parameters(chemostat)
#  μmax
#  Ks
#  c
#  q
#  V

ode = convert(ODESystem, chemostat, combinatoric_ratelaws=false)
equations(ode)
#  Differential(t)(S(t)) ~ (-q*S(t)) / V + (c*q) / V - X(t)*Catalyst.mm(S(t), μmax, Ks)
#  Differential(t)(X(t)) ~ (-q*X(t)) / V + 0.8X(t)*Catalyst.mm(S(t), μmax, Ks)
states(ode)
#  S(t)
#  X(t)
parameters(ode)
#  μmax
#  Ks
#  c
#  q
#  V

################################################################
## Simulating the system
################################################################

using DifferentialEquations, Plots

u₀map = [:S => 0.0, :X => 0.1]
# :S => 0.0
# :X => 0.1

pmap  = [:μmax => 1.7, :Ks => 29, :c=>50, :q => 10, :V => 100]
#  :μmax => 1.7
#    :Ks => 29.0
#     :q => 10.0
#     :V => 100.0
#     :c => 50.0

tspan = (0.0, 24.0)

oprob = ODEProblem(chemostat, u₀map, tspan, pmap, combinatoric_ratelaws=false)

osol = solve(oprob, Tsit5(), saveat=0.5)

plot(osol)



################################################################
## Finding steady state
################################################################

ss_guess = [:X => osol[:X][end], :S => osol[:S][end]]  # last value as guess
# 2-element Vector{Pair{Symbol, Float64}}:
#  :X => 34.31236861749778
#  :S => 2.5849812423175864

Seq, Xeq = solve(SteadyStateProblem(ODEProblem(chemostat, ss_guess, tspan, pmap, combinatoric_ratelaws=false)))
# retcode: Success
# u: 2-element Vector{Float64}:
#   2.3015873015873014
#  38.15873015873016

Seq, Xeq
# (2.3015873015873014, 38.15873015873016)

################################################################
## Parameter uncertainty
################################################################

using Measurements

# ± with \pm
pmap_uncertainty  = [:μmax => 1.7±0.2, :Ks => 29±3, :c=>50, :q => 10, :V => 100]
# :μmax => 1.7 ± 0.2
# :Ks => 29.0 ± 3.0
#  :c => 50.0 ± 0.0
#  :q => 10.0 ± 0.0
#  :V => 100.0 ± 0.0

oprob_uncertainty = ODEProblem(chemostat, u₀map, tspan, pmap_uncertainty, combinatoric_ratelaws=false)

osol_uncertainty = solve(oprob_uncertainty, Tsit5(), saveat=0.5)

plot(osol_uncertainty)


################################################################
## Sensitivity analysis
################################################################

# import Pkg; Pkg.add("SciMLSensitivity")
using SciMLSensitivity

oprob_sens = ODEForwardSensitivityProblem(oprob.f, [0.0, 0.1], tspan, [1.7, 29, 50, 10, 100])
# oprob.f    : the ODE function
# [0.0, 0.1] : initial conditions, the order is [S, X], see species(chemostat)
# tspan      : the time span
# [1.7, 29, 50, 10, 100] : the parameter values, the order is [μmax, Ks, c, q, V], see parameters(chemostat)

osol_sens = solve(oprob_sens, saveat=0.1)
# osol_sens.t : the vector with the time points
size(osol_sens.t)
# osol_sens.u : vector of vector with solution and sensitivities
size(osol_sens.u)

u, dp = extract_local_sensitivities(osol_sens)

# u  : solution for S and X of the ODE problem
size(u)
# (2, 241)
# u[1,:] : solution for S
# u[2,:] : solution for X
plot(osol_sens.t, u', title="Solution S and X", label=["S" "X"])
#  Note that there is no comma in label=["S" "X"] !!!

# dp : sensitivities of S and X (cf. dS/dt and dX/dt) with respect to all parameters
size(dp)
# (5,)          # because there are 5 parameters
sens_c = dp[3]  # c is the 3rd parameter in [μmax, Ks, c, q, V]
size(sens_c)
# (2, 241)    # two row because we have two variables (S and X)
# sens_c[1,:]     # sensitivity for S on c
# sens_c[2,:]     # sensitivity for X on c

plot(osol_sens.t, sens_c', title="Sensitivity for S and X on the ingoing concentration c", label=["dS/dt" "dX/dt"])


################################################################
## Calibration
################################################################


Xmeas =  [0.2,  0.52,  1.12,  2.5,  5.52,  11.98,  24.6,  45.15,  66.66,  77.16]
tmeas = 2:2:20
scatter(tmeas, Xmeas, label="biomass", xlabel="t")

# Xtmp = Xmeas .+ 2
# abs2.(Xtmp - Xmeas)
# sum(abs2, Xtmp - Xmeas)  # same as: sum((Xtmp-Xmeas).^2)

function loss(μmax, Ks)
    u₀map = [:X => 0.1, :S => 100]
    pmap  = [:μmax => μmax, :Ks => Ks, :q => 0, :V =>1, :c=>50]
    tspan = (0.0, 24.0)
    oprob = ODEProblem(chemostat, u₀map, tspan, pmap, combinatoric_ratelaws=false)
    sol = solve(oprob, Tsit5(), saveat=tmeas)
    Xsim = sol[:X]
    mse = sum(abs2, Xsim - Xmeas)
    return log(mse)  # we take a logarithm because there is a huge difference in scale
end

contourf(0.1:0.1:3, 1:0.5:100, loss, xlab="μmax", ylab="Ks", title="log(MSE)")

(minimum_val, (μmax_best, Ks_best)) = minimum((loss(μmax, Ks), (μmax, Ks)) for μmax in 0.1:0.1:4 for Ks in 1:0.5:100)
# (-6.001829196352962, (0.9, 78.0))

u₀map = [:X => 0.1, :S => 100]
pmap  = [:μmax => μmax_best, :Ks => Ks_best, :q => 0, :V =>1, :c=>50]
tspan = (0.0, 24.0)
oprob = ODEProblem(chemostat, u₀map, tspan, pmap, combinatoric_ratelaws=false)
osol = solve(oprob, Tsit5(), saveat=tmeas)
plot(osol)
scatter!(tmeas, Xmeas, label="biomass", xlabel="t")


# Using an optimization function

function Jtheta(thetas)           # function return to be minimized
    u₀map = [:X => 0.1, :S => 100]
    pmap  = [:μmax => thetas[1], :Ks => thetas[2], :q => 0, :V => 1, :c => 50]
    tspan = (0.0, 24.0)
    oprob = ODEProblem(chemostat, u₀map, tspan, pmap, combinatoric_ratelaws=false)
    sol = solve(oprob, Tsit5(), saveat=tmeas)
    Xsim = sol[:X]
    mse = sum(abs2, Xsim - Xmeas)
    return log(mse)  # we take a logarithm because there is a huge difference in scale
end

using Optim

results = optimize(Jtheta, [0.9, 78.0], NelderMead())
Optim.minimizer(results)
# 0.8994536987092343
# 77.9150642441882