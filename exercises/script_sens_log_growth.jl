using Catalyst
using DifferentialEquations, Plots


################################################################
## Modelling Logistic growth
################################################################

# growth_mod_log = @reaction_network begin
#     -μ*(1 - W/Wf), W --> ∅
# end

growth_mod_log = @reaction_network begin
    μ*W, ∅ --> W
    μ/Wf*W, W --> ∅
end

osys_log  = convert(ODESystem, growth_mod_log)
equations(osys_log)

parameters(growth_mod_log)
# μ
# Wf

osys_log  = convert(ODESystem, growth_mod_log)

# Check out the diff. eqns. equations
equations(osys_log)
# Differential(t)(W(t)) ~ (-(W(t)^2)*μ) / Wf + W(t)*μ

# 2.0 0.07 10.0

# Set initial condition for ODE problem
u₀_log = [:W => 2.0]
# Set time span for simulation
tspan = (0.0, 100.0)
# Set parameters for ODE problem
params_log = [:μ => 0.07, :Wf => 10.0]
# Create ODE problem
oprob_log = ODEProblem(growth_mod_log, u₀_log, tspan, params_log)

osol_log = solve(oprob_log, Tsit5(), saveat=0.5)
plot(osol_log)


################################################################
## Sensitivity analysis
################################################################

using SciMLSensitivity

# 2.0 0.065 9.8

oprob_sens_log = ODEForwardSensitivityProblem(oprob_log.f, [2.0], tspan, [0.07, 10.0])

osol_sens_log = solve(oprob_sens_log, Tsit5(), saveat=0.5)
u_log, dp_log = extract_local_sensitivities(osol_sens_log)

size(u_log)
size(dp_log)

W_sol_log = u_log[1,:]     # solution for W

sens_μ_log  = dp_log[1]'  # sensitivity for W on μ
sens_Wf_log = dp_log[2]'  # sensitivity for W on Wf

sens_μ_log[5]

size(osol_sens_log.t)
size(sens_μ_log)

maximum(sens_μ_log)
max_val, index = findmax(sens_μ_log)
osol_sens_log.t[index]

findmax(sens_Wf_log)

# Absolute sensitivity for W on μ
plot(osol_sens_log.t, sens_μ_log, title="Absolute sensitivity of W to μ", label=["dW/dμ"], xlabel="Time (day)")

# Absolute sensitivity for W on Wf
plot(osol_sens_log.t, sens_Wf_log, title="Absolute sensitivity of W to Wf", label=["dW/dWf"], xlabel="Time (day)")


################################################################
## Uncertainty analysis
################################################################

using Measurements

params_uncert_log = [:μ => 0.07±0.02, :Wf => 10.0]
params_uncert_log = [:μ => 0.07, :Wf => 10.0±1.5]
params_uncert_log = [:μ => 0.07±0.02, :Wf => 10.0±1.5]
# Create ODE problem
oprob_uncert_log = ODEProblem(growth_mod_log, u₀_log, tspan, params_uncert_log)

osol_uncert_log = solve(oprob_uncert_log, Tsit5(), saveat=2.0)
plot(osol_uncert_log)



################################################################
## Calibration
################################################################

W_meas = [2.3, 4.5, 6.6, 7.6, 9.0, 9.1, 9.4]
t_meas = [0, 20, 29, 41, 50, 65, 72]
W_sigma = 0.5
scatter(t_meas, W_meas, title="Grass growth data",
                        label="Yield",
                        xlabel="t",
                        xlims=(0, 80),
                        ylims=(0, 10))

function Jtheta_log(thetas)           # function return to be minimized
    u₀ = [:W => thetas[1]]
    params  = [:μ => thetas[2], :Wf => thetas[3]]
    oprob = ODEProblem(growth_mod_log, u₀, tspan, params, combinatoric_ratelaws=false)
    sol = solve(oprob, Tsit5(), saveat=t_meas)
    W_sol = sol[:W]
    J = (1/W_sigma^2)*sum(abs2, W_sol - W_meas)
    return J
end

# 2.0 0.065 9.8

using Optim

results_log = optimize(Jtheta_log, [2.0, 0.07, 10.0], NelderMead())
opt_log = Optim.minimizer(results_log)
# 2.049546098598697
# 0.06665793744212963
# 9.735504889362218
Optim.minimum(results_log)    # J-value / mse-value
# 2.2056801547033844

u₀_opt_log = [:W => opt_log[1]]
params_opt_log = [:μ => opt_log[2], :Wf => opt_log[3]]
oprob_opt_log = ODEProblem(growth_mod_log, u₀_opt_log, tspan, params_opt_log)
osol_opt_log = solve(oprob_opt_log, Tsit5(), saveat=0.5)
plot(osol_opt_log, label="Logistic growth", xlabel="t", xlims=(0, 80), ylims=(0, 10))
scatter!(t_meas, W_meas, label="Yield")


# osol_opt_log = solve(oprob_opt_log, Tsit5(), saveat=t_meas)
# sqrt(sum(abs2, osol_opt_log[:W] .- W_meas) / (length(W_meas)-1))
# # 0.3031556218757122

# using Statistics
# std(osol_opt_log[:W] .- W_meas)
# # 0.30269436851836035