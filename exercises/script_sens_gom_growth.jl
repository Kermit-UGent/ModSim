using Catalyst
using DifferentialEquations, Plots


################################################################
## Modelling Gompertz growth
################################################################

growth_mod_gom = @reaction_network begin
    -μ, W --> ∅
    D*log(W), W --> ∅
end

parameters(growth_mod_gom)
# μ
# D

osys_gom  = convert(ODESystem, growth_mod_gom)

# Check out the diff. eqns. equations
equations(osys_gom)
# Differential(t)(W(t)) ~ W(t)*μ - D*W(t)*log(W(t))

# Set initial condition for ODE problem
u₀_gom = [:W => 2.0]
# Set time span for simulation
tspan = (0.0, 100.0)
# Set parameters for ODE problem
params_gom = [:μ => 0.09, :D => 0.04]
# Create ODE problem
oprob_gom = ODEProblem(growth_mod_gom, u₀_gom, tspan, params_gom)

osol_gom = solve(oprob_gom, Tsit5(), saveat=0.5)
plot(osol_gom)


################################################################
## Sensitivity analysis
################################################################

using SciMLSensitivity

oprob_sens_gom = ODEForwardSensitivityProblem(oprob_gom.f, [2.0], tspan, [0.09, 0.04])

osol_sens_gom = solve(oprob_sens_gom, Tsit5(), saveat=0.5)
u_gom, dp_gom = extract_local_sensitivities(osol_sens_gom)

size(u_gom)
size(dp_gom)

W_sol_gom = u_gom[1,:]     # solution for W

sens_μ_gom = dp_gom[1]'  # sensitivity for W on μ
sens_D_gom = dp_gom[2]'  # sensitivity for W on Wf

# Absolute sensitivity of W to μ
plot(osol_sens_gom.t, sens_μ_gom, title="Absolute sensitivity of W to μ", label=["dW/dμ"], xlabel="Time (day)")

# Absolute sensitivity of W to D
plot(osol_sens_gom.t, sens_D_gom, title="Absolute sensitivity of W to D", label=["dW/dD"], xlabel="Time (day)")


################################################################
## Uncertainty analysis
################################################################

using Measurements

params_uncert_gom = [:μ => 0.09±0.01, :D => 0.040]
params_uncert_gom = [:μ => 0.09, :D => 0.040±0.002]
params_uncert_gom = [:μ => 0.09±0.01, :D => 0.040±0.002]
# Create ODE problem
oprob_uncert_gom = ODEProblem(growth_mod_gom, u₀_gom, tspan, params_uncert_gom)

osol_uncert_gom = solve(oprob_uncert_gom, Tsit5(), saveat=2.0)
plot(osol_uncert_gom)



################################################################
## Calibration
################################################################

W_meas = [2.3, 4.5, 6.6, 7.6, 9.0, 9.1, 9.4]
t_meas = [0, 20, 29, 41, 50, 65, 72]
scatter(t_meas, W_meas, label="Growth", xlabel="t", xlims=(0, 80), ylims=(0, 10))

function Jtheta_gom(thetas)           # function return to be minimized
    u₀ = [:W => thetas[1]]
    params  = [:μ => thetas[2], :D => thetas[3]]
    oprob = ODEProblem(growth_mod_gom, u₀, tspan, params, combinatoric_ratelaws=false)
    sol = solve(oprob, Tsit5(), saveat=t_meas)
    W_sol = sol[:W]
    J = (1/W_sigma^2)*sum(abs2, W_sol - W_meas)
    return J
end

# W₀ = 2.0
# μ = 0.0964
# D = 0.041

using Optim

results_gom = optimize(Jtheta_gom, [2.0, 0.09, 0.04], NelderMead())
opt_gom = Optim.minimizer(results_gom)
# 2.0424460825820874
# 0.09620691573096914
# 0.04108675572344885
Optim.minimum(results_gom)    # J-value / mse-value
# 3.420229732575734

u₀_opt_gom = [:W => opt_gom[1]]
params_opt_gom = [:μ => opt_gom[2], :D => opt_gom[3]]
oprob_opt_gom = ODEProblem(growth_mod_gom, u₀_opt_gom, tspan, params_opt_gom)
osol_opt_gom = solve(oprob_opt_gom, Tsit5(), saveat=0.5)
plot(osol_opt_gom, label="Gompertz growth", xlabel="t", xlims=(0, 80), ylims=(0, 10))
scatter!(t_meas, W_meas, label="Yield")

