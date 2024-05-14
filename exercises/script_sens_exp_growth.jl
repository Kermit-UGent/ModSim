using Catalyst
using DifferentialEquations, Plots


################################################################
## Modelling Exponential growth
################################################################

growth_mod_exp = @reaction_network begin
    μ*Wf, ∅ --> W
    μ, W --> ∅
end

parameters(growth_mod_exp)
# μ
# Wf

osys_exp = convert(ODESystem, growth_mod_exp)

# Check out the diff. eqns. equations
equations(osys_exp)
#  Differential(t)(W(t)) ~ Wf*μ - W(t)*μ

# W₀ = 2.0
# μ = 0.02
# Wf = 10.0

# Set initial condition for ODE problem
u₀_exp = [:W => 2.0]
# Set time span for simulation
tspan = (0.0, 100.0)
# Set parameters for ODE problem
params_exp = [:μ => 0.02, :Wf => 10.0]
# Create ODE problem
oprob_exp = ODEProblem(growth_mod_exp, u₀_exp, tspan, params_exp)

osol_exp = solve(oprob_exp, Tsit5(), saveat=0.5)
plot(osol_exp)


################################################################
## Sensitivity analysis
################################################################

using SciMLSensitivity

oprob_sens_exp = ODEForwardSensitivityProblem(oprob_exp.f, [2.0], tspan, [0.02, 10.0])

osol_sens_exp = solve(oprob_sens_exp, Tsit5(), saveat=0.5)
u_exp, dp_exp = extract_local_sensitivities(osol_sens_exp)

size(u_exp)
size(dp_exp)

W_sol_exp = u_exp[1,:]     # solution for W

sens_μ_exp  = dp_exp[1]'  # sensitivity of W to μ
sens_Wf_exp = dp_exp[2]'  # sensitivity of W to Wf

# Absolute sensitivity for W on μ
plot(osol_sens_exp.t, sens_μ_exp, title="Absolute sensitivity of W to μ", label=["dW/dμ"], xlabel="Time (day)")

# Absolute sensitivity for W on Wf
plot(osol_sens_exp.t, sens_Wf_exp, title="Absolute sensitivity of W to Wf", label=["dW/dWf"], xlabel="Time (day)")


################################################################
## Uncertainty analysis
################################################################

using Measurements

params_uncert_exp = [:μ => 0.02±0.01, :Wf => 10.0]
params_uncert_exp = [:μ => 0.02, :Wf => 10.0±1.5]
params_uncert_exp = [:μ => 0.02±0.01, :Wf => 10.0±1.5]
# Create ODE problem
oprob_uncert_exp = ODEProblem(growth_mod_exp, u₀_exp, tspan, params_uncert_exp)

osol_uncert_exp = solve(oprob_uncert_exp, Tsit5(), saveat=2.0)
plot(osol_uncert_exp)


################################################################
## Calibration
################################################################

W_meas = [2.3, 4.5, 6.6, 7.6, 9.0, 9.1, 9.4]
t_meas = [0, 20, 29, 41, 50, 65, 72]
scatter(t_meas, W_meas, label="Growth", xlabel="t", xlims=(0, 80), ylims=(0, 10))

function Jtheta_exp(thetas)           # function return to be minimized
    u₀ = [:W => thetas[1]]
    params  = [:μ => thetas[2], :Wf => thetas[3]]
    oprob = ODEProblem(growth_mod_exp, u₀, tspan, params, combinatoric_ratelaws=false)
    sol = solve(oprob, Tsit5(), saveat=t_meas)
    W_sol = sol[:W]
    J = (1/W_sigma^2)*sum(abs2, W_sol - W_meas)
    return J
end

using Optim

results_exp = optimize(Jtheta_exp, [2.0, 0.02, 10.0], NelderMead())
opt_exp = Optim.minimizer(results_exp)
# 2.089662889475103
# 0.01853062477032477
# 12.358676119868713
Optim.minimum(results_exp)    # J-value / mse-value
# 5.175755636174577

u₀_opt_exp = [:W => opt_exp[1]]
params_opt_exp = [:μ => opt_exp[2], :Wf => opt_exp[3]]
oprob_opt_exp = ODEProblem(growth_mod_exp, u₀_opt_exp, tspan, params_opt_exp)
osol_opt_exp = solve(oprob_opt_exp, Tsit5(), saveat=0.5)
plot(osol_opt_exp, label="Exponential growth", xlabel="t", xlims=(0, 80), ylims=(0, 10))
scatter!(t_meas, W_meas, label="Yield")

