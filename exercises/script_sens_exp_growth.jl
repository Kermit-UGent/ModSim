using Catalyst
using DifferentialEquations, Plots


################################################################
## Modelling Exponential growth
################################################################

growth_exp = @reaction_network begin
	@species W(t)=2.0
	@parameters μ=0.02 Wf=10.0
    μ*Wf, ∅ --> W
    μ, W --> ∅
end

parameters(growth_exp)
# μ
# Wf

osys_exp = convert(ODESystem, growth_exp)

# Check out the diff. eqns. equations
equations(osys_exp)
#  Differential(t)(W(t)) ~ Wf*μ - W(t)*μ

# W₀ = 2.0
# μ = 0.02
# Wf = 10.0

# Set initial condition for ODE problem
u0_exp = [:W => 2.0]
# Set time span for simulation
tspan = (0.0, 100.0)
# Set parameters for ODE problem
params_exp = [:μ => 0.02, :Wf => 10.0]
# Create ODE problem
oprob_exp = ODEProblem(growth_exp, u0_exp, tspan, params_exp)

osol_exp = solve(oprob_exp, Tsit5(), saveat=0.5)
plot(osol_exp)


################################################################
## Sensitivity analysis
################################################################

using ForwardDiff

function growth_exp_sim(params)
	μ, Wf = params
    u0_exp = [:W => 2.0]
    tspan = (0.0, 100.0)
    oprob_exp = ODEProblem(growth_exp, u0_exp, tspan, [:μ=>μ, :Wf=>Wf])
	osol_exp = solve(oprob_exp, Tsit5(), saveat=0.5)
	return osol_exp
end

growth_exp_sim([0.02, 10])

growth_exp_sim_W(params) = growth_exp_sim(params)[:W]

t_vals = 0:0.5:100.0

W_exp = growth_exp_sim_W([0.02, 10.0])

sens_W_exp = ForwardDiff.jacobian(growth_exp_sim_W, [0.02, 10.0])

sens_W_on_μ_exp = sens_W_exp[:,1]
sens_W_on_Wf_exp = sens_W_exp[:,2]

sens_W_on_μ_rel_exp = sens_W_on_μ_exp .* 0.02 ./ W_exp
sens_W_on_Wf_rel_exp = sens_W_on_Wf_exp .* 10.0 ./ W_exp

plot(t_vals, [sens_W_on_μ_rel_exp, sens_W_on_Wf_rel_exp], title="Normalized sensitivities", label=["W on μ" "W on Wf"], xlabel="Time (day)")


# ################################################################
# ## Sensitivity analysis
# ################################################################

# using SciMLSensitivity

# oprob_sens_exp = ODEForwardSensitivityProblem(oprob_exp.f, [2.0], tspan, [0.02, 10.0])

# osol_sens_exp = solve(oprob_sens_exp, Tsit5(), saveat=0.5)
# u_exp, dp_exp = extract_local_sensitivities(osol_sens_exp)

# size(u_exp)
# size(dp_exp)

# W_sol_exp = u_exp[1,:]     # solution for W

# sens_μ_exp  = dp_exp[1]'  # sensitivity of W to μ
# sens_Wf_exp = dp_exp[2]'  # sensitivity of W to Wf

# # Absolute sensitivity for W on μ
# plot(osol_sens_exp.t, sens_μ_exp, title="Absolute sensitivity of W to μ", label=["dW/dμ"], xlabel="Time (day)")

# # Absolute sensitivity for W on Wf
# plot(osol_sens_exp.t, sens_Wf_exp, title="Absolute sensitivity of W to Wf", label=["dW/dWf"], xlabel="Time (day)")


################################################################
## Uncertainty analysis
################################################################

using Measurements

params_uncert_exp = [:μ => 0.02±0.01, :Wf => 10.0]
params_uncert_exp = [:μ => 0.02, :Wf => 10.0±1.5]
params_uncert_exp = [:μ => 0.02±0.01, :Wf => 10.0±1.5]
# Create ODE problem
oprob_uncert_exp = ODEProblem(growth_exp, u₀_exp, tspan, params_uncert_exp)

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
    oprob = ODEProblem(growth_exp, u₀, tspan, params, combinatoric_ratelaws=false)
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
oprob_opt_exp = ODEProblem(growth_exp, u₀_opt_exp, tspan, params_opt_exp)
osol_opt_exp = solve(oprob_opt_exp, Tsit5(), saveat=0.5)
plot(osol_opt_exp, label="Exponential growth", xlabel="t", xlims=(0, 80), ylims=(0, 10))
scatter!(t_meas, W_meas, label="Yield")


################################################################
## Calibration 2
################################################################

using Optim, Turing, StatsPlots
using LinearAlgebra
using StatsBase

W_meas = [1.87, 2.45, 3.72, 4.32, 5.28, 7.01, 6.83, 8.62, 9.45, 10.31, 10.56, 11.72, 11.05, 11.53, 11.39, 11.7, 11.15, 11.49, 12.04, 11.95, 11.68]
# W_meas = [1.97, 2.14, 3.16, 4.14, 4.69, 5.48, 6.54, 7.31, 7.45, 8.79, 8.04, 9.77, 9.78, 9.36, 9.62, 9.81, 10.26, 9.22, 9.67, 9.77, 9.74]
t_meas = 0:5:100
scatter(t_meas, W_meas, title="Grass growth data",
    label="Yield", xlabel="t", xlims=(0, 100), ylims=(0, 14))

W_meas[1]    # 1.87
W_meas[end]  # 11.68
(W_meas[end] - W_meas[1]) / t_meas[end]   # 0.0981


##########################################################################
# Setting up the Turing model
##########################################################################

u₀_exp = [:W => 2.0]
tspan = (0.0, 100.0)
params_exp = [:μ => 0.02, :Wf => 10.0]
oprob_exp = ODEProblem(growth_exp, u₀_exp, tspan, params_exp)


@model function growth_gom_inference(t_meas, W)
    σ_W ~ InverseGamma()
    W₀ ~ Uniform(0, 10)
    μ ~ Uniform(0, 1)
    Wf ~ Uniform(0, 100)
    # IMPORTANT:
    # - If you want to overwrite the initial condition, use
    #   remake(oprob_log; u0=[W₀]) instead of just oprob_log.
    # - If you want to overwrite the parameters, place all
    #   parameters in the option p=[...] in the same order
    #   as in the reaction network model.
    osol_exp = solve(remake(oprob_exp; u0=[W₀]), Tsit5(), saveat=t_meas, p=[μ, Wf])
    W ~ MvNormal(osol_exp[:W], σ_W^2 * I)
end

growth_gom_inf = growth_gom_inference(t_meas, W_meas)


##########################################################################
# MLE = Maximum Likelihood Estimation
##########################################################################

results_gom_mle = optimize(growth_gom_inf, MLE(), NelderMead())
# ModeResult with maximized lp of -18.25
# [0.576955037418865, 0.9583148988579651, 0.02804914844427458, 13.022008899256887]
results_gom_mle |> coeftable
# ──────────────────────────────────────────────────────────────────────────
#           Coef.  Std. Error         z      Pr(>|z|)  Lower 95%   Upper 95%
# ──────────────────────────────────────────────────────────────────────────
# σ_W   0.576955   0.0890281    6.48059  9.13625e-11    0.402463   0.751447
# W₀    0.958315   0.410975     2.33181  0.0197107      0.152819   1.76381
# μ     0.0280491  0.00308942   9.07911  1.09469e-19    0.021994   0.0341043
# Wf   13.022      0.452884    28.7535   8.17814e-182  12.1344    13.9096
# ──────────────────────────────────────────────────────────────────────────

coefnames(results_gom_mle)
#  :σ_W
#  :W₀
#  :μ
#  :Wf

coef(results_gom_mle)
# 4-element Named Vector{Float64}
# A   │ 
# ────┼──────────
# σ_W │  0.576955
# W₀  │  0.958315
# μ   │ 0.0280491
# Wf  │    13.022

stderror(results_gom_mle)
# 4-element Named Vector{Float64}
# A   │ 
# ────┼───────────
# σ_W │  0.0890281
# W₀  │   0.410975
# μ   │ 0.00308942
# Wf  │   0.452884

coef(results_gom_mle)[:σ_W]
coef(results_gom_mle)[:W₀]
coef(results_gom_mle)[:μ]
coef(results_gom_mle)[:Wf]

W₀_opt_exp = coef(results_gom_mle)[:W₀]
μ_opt_exp = coef(results_gom_mle)[:μ]
Wf_opt_exp = coef(results_gom_mle)[:Wf]

u₀_opt_exp = [:W => W₀_opt_exp]
tspan = (0.0, 100.0)
params_opt_exp = [:μ => μ_opt_exp, :Wf => Wf_opt_exp]
oprob_opt_exp = ODEProblem(growth_exp, u₀_opt_exp, tspan, params_opt_exp)
osol_opt_exp = solve(oprob_opt_exp, Tsit5(), saveat=0.5)
plot(osol_opt_exp, label="Exponential growth", xlabel="t", xlims=(0, 100), ylims=(0, 14))
scatter!(t_meas, W_meas, label="Yield")


##########################################################################
# MAP = Maximum A Posterior
##########################################################################

results_gom_map = optimize(growth_gom_inf, MAP(), NelderMead())
# ModeResult with maximized lp of -19.68
# [0.5753665525787777, 1.1151426631391854, 0.027287620330604608, 13.093154662763009]
results_gom_map |> coeftable
# ──────────────────────────────────────────────────────────────────────────
#           Coef.  Std. Error         z      Pr(>|z|)   Lower 95%  Upper 95%
# ──────────────────────────────────────────────────────────────────────────
# σ_W   0.575367   0.0869705    6.61565  3.69923e-11    0.404907    0.745826
# W₀    1.11514    0.377315     2.95547  0.00312193     0.375619    1.85467
# μ     0.0272876  0.00293547   9.29584  1.46049e-20    0.0215342   0.033041
# Wf   13.0932     0.455917    28.7183   2.25457e-181  12.1996     13.9867
# ──────────────────────────────────────────────────────────────────────────


##########################################################################
# NUTS = No-U-Turn Sampler
##########################################################################

iterations = 100;
results_gom_nuts = sample(growth_gom_inf, NUTS(), iterations)
summarize(results_gom_nuts)
# parameters      mean       std      mcse   ess_bulk   ess_tail      rhat   ess_per_sec 
# Symbol   Float64   Float64   Float64    Float64    Float64   Float64       Float64 
#    σ_W    0.6107    0.0926    0.0325     9.1883    23.7338    1.0755        1.0051
#     W₀    1.1870    0.3920    0.2992     1.9982    31.6824    1.5206        0.2186
#      μ    0.0272    0.0025    0.0016     2.2759    39.2638    1.4013        0.2490
#     Wf   13.1114    0.3441    0.2356     2.0701    21.0246    1.5079        0.2264

iterations = 500;
results_gom_nuts = sample(growth_gom_inf, NUTS(), iterations)
summarize(results_gom_nuts)
# parameters      mean       std      mcse   ess_bulk   ess_tail      rhat   ess_per_sec 
# Symbol   Float64   Float64   Float64    Float64    Float64   Float64       Float64 
#    σ_W    0.6771    0.1339    0.0117   127.5052   143.3003    1.0019        8.5666
#     W₀    1.1532    0.4832    0.0387   156.4447   180.4347    1.0136       10.5109
#      μ    0.0273    0.0036    0.0003   192.2302   217.8047    1.0030       12.9152
#     Wf   13.1411    0.5885    0.0456   194.7689   193.9734    1.0005       13.0858