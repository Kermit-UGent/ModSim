
using Catalyst
using DifferentialEquations, Plots

#======================================================================================
  Solving using Catalyst
======================================================================================#

# Define all reactions
water_evap_infil = @reaction_network begin
    I, ∅ --> W
    -O, ∅ --> W
    k₁, W --> ∅
    k₂, W --> G
    k₂, G --> W
end

parameters(water_evap_infil)
# I
# O
# k₁
# k₂
species(water_evap_infil)
# W(t)
# G(t)

osys  = convert(ODESystem, water_evap_infil)

# Check out the diff. eqns. equations
equations(osys)
# Differential(t)(W(t)) ~ I - O - k₁*W(t) + k₂*G(t) - k₂*W(t)
# Differential(t)(G(t)) ~ -k₂*G(t) + k₂*W(t)

I = 2.7
O = 20.0
k₁ = 0.4    # 1e-3 to 1e-1 per day
k₂ = 1.0    # 1e-3 to 1e-1 per day

W₀ = 6.75
G₀ = 6.75   # 5.5

u₀ = [:W => W₀, :G => G₀]
# Set time span for simulation
tspan = (0.0, 20.0)
# Set parameters for ODE problem
params = [:I => I, :O => O, :k₁ => k₁, :k₂ => k₂]

# u_guess = [:W => 6, :G => 6]
# Seq, Xeq = solve(SteadyStateProblem(ODEProblem(water_evap_infil, u_guess, tspan, params, combinatoric_ratelaws=false)))
# 6.75
# 6.75

# Create ODE problem
oprob = ODEProblem(water_evap_infil, u₀, tspan, params, combinatoric_ratelaws=false)

# osol = solve(oprob)
# plot(osol)


function condition(u, t, integrator)
  u[1]
end
function affect!(integrator)
  integrator.p[2] = 0.0
end
ps_cb = ContinuousCallback(condition, affect!)
osol = solve(deepcopy(oprob), Tsit5(), saveat=0.1, callback=ps_cb)
plot(osol)



################################################################
## Sensitivity analysis
################################################################

using SciMLSensitivity

oprob_sens = ODEForwardSensitivityProblem(oprob.f, [W₀, G₀], tspan, [I, O, k₁, k₂])

osol_sens = solve(oprob_sens, Tsit5(), saveat=0.1)
u, dp = extract_local_sensitivities(osol_sens)

sens_I = dp[1]
sens_O = dp[2]
sens_k1 = dp[3]
sens_k2 = dp[4]

# Absolute sensitivity for W and G on I
plot(osol_sens.t, sens_I', title="Absolute sensitivity for W and G on I", label=["dW/dI" "dG/dI"])
plot(osol_sens.t, sens_O', title="Absolute sensitivity for W and G on O", label=["dW/dO" "dG/dO"])
plot(osol_sens.t, sens_k1', title="Absolute sensitivity for W and G on k₁", label=["dW/dk₁" "dG/dk₁"])
plot(osol_sens.t, sens_k2', title="Absolute sensitivity for W and G on k₂", label=["dW/dk₂" "dG/dk₂"])