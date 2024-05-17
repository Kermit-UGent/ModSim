using Catalyst
using DifferentialEquations, Plots

#=
https://discourse.julialang.org/t/get-mle-parameters-e-g-p-value-confidence-intervals-from-a-turing-model-using-optim-jl/101433/3
=#


################################################################
## Modelling Logistic growth
################################################################

# growth_mod_log = @reaction_network begin
#     -μ*(1 - W/Wf), W --> ∅
# end

growth_log = @reaction_network begin
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
oprob_log = ODEProblem(growth_log, u₀_log, tspan, params_log)
# oprob_log = remake(oprob_log; u0=[4.0])

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
oprob_uncert_log = ODEProblem(growth_log, u₀_log, tspan, params_uncert_log)

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
    oprob = ODEProblem(growth_log, u₀, tspan, params, combinatoric_ratelaws=false)
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
oprob_opt_log = ODEProblem(growth_log, u₀_opt_log, tspan, params_opt_log)
osol_opt_log = solve(oprob_opt_log, Tsit5(), saveat=0.5)
plot(osol_opt_log, label="Logistic growth", xlabel="t", xlims=(0, 80), ylims=(0, 10))
scatter!(t_meas, W_meas, label="Yield")


# osol_opt_log = solve(oprob_opt_log, Tsit5(), saveat=t_meas)
# sqrt(sum(abs2, osol_opt_log[:W] .- W_meas) / (length(W_meas)-1))
# # 0.3031556218757122

# using Statistics
# std(osol_opt_log[:W] .- W_meas)
# # 0.30269436851836035


################################################################
## Calibration 2
################################################################

using Optim, Turing, StatsPlots
using LinearAlgebra
using StatsBase

##########################################################################
# Creating measured data
##########################################################################

u₀_log = [:W => 1.6]
tspan = (0.0, 100.0)
params_log = [:μ => 0.084, :Wf => 11.6]
oprob_log = ODEProblem(growth_log, u₀_log, tspan, params_log)
osol_log = solve(oprob_log, Tsit5(), saveat=5.0)
plot(osol_log)
σ_W = 0.35
Wobs = abs.(osol_log[:W] .+ σ_W .* randn(length(osol_log)))
scatter(osol_log.t, Wobs)

println(round.(Wobs, digits=2))

W_meas = [1.87, 2.45, 3.72, 4.32, 5.28, 7.01, 6.83, 8.62, 9.45, 10.31, 10.56, 11.72, 11.05, 11.53, 11.39, 11.7, 11.15, 11.49, 12.04, 11.95, 11.68]
# W_meas = [1.97, 2.14, 3.16, 4.14, 4.69, 5.48, 6.54, 7.31, 7.45, 8.79, 8.04, 9.77, 9.78, 9.36, 9.62, 9.81, 10.26, 9.22, 9.67, 9.77, 9.74]
t_meas = 0:5:100
scatter(t_meas, W_meas, title="Grass growth data",
    label="Yield", xlabel="t", xlims=(0, 100), ylims=(0, 14))

# W_meas = [2.3, 4.5, 6.6, 7.6, 9.0, 9.1, 9.4]
# t_meas = [0, 20, 29, 41, 50, 65, 72]

W_meas[1]    # 1.87
W_meas[end]  # 11.68
(W_meas[end] - W_meas[1]) / t_meas[end]   # 0.0981

##########################################################################
# Setting up the Turing model
##########################################################################

@model function growth_log_inference(t_meas, W)
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
    osol_log = solve(remake(oprob_log; u0=[W₀]), Tsit5(), saveat=t_meas, p=[μ, Wf])
    W ~ MvNormal(osol_log[:W], σ_W^2 * I)
end

growth_log_inf = growth_log_inference(t_meas, W_meas)


##########################################################################
# MLE = Maximum Likelihood Estimation
##########################################################################

results_log_mle = optimize(growth_log_inf, MLE(), NelderMead())
# ModeResult with maximized lp of -6.20
# [0.3250847442926941, 1.8475435183171234, 0.07641755219669163, 11.835357113482539]
results_log_mle |> coeftable
# ──────────────────────────────────────────────────────────────────────────
#           Coef.  Std. Error         z     Pr(>|z|)   Lower 95%   Upper 95%
# ──────────────────────────────────────────────────────────────────────────
# σ_W   0.325085   0.0501555    6.48154  9.07939e-11   0.226782    0.423388
# W₀    1.84754    0.144674    12.7704   2.4003e-37    1.56399     2.1311
# μ     0.0764176  0.00385156  19.8407   1.32662e-87   0.0688686   0.0839665
# Wf   11.8354     0.127587    92.763    0.0          11.5853     12.0854
# ──────────────────────────────────────────────────────────────────────────

coefnames(results_log_mle)
#  :σ_W
#  :W₀
#  :μ
#  :Wf

coef(results_log_mle)
# 4-element Named Vector{Float64}
# A   │ 
# ────┼──────────
# σ_W │  0.325085
# W₀  │   1.84754
# μ   │ 0.0764176
# Wf  │   11.8354

coef(results_log_mle)[:σ_W]
coef(results_log_mle)[:W₀]
coef(results_log_mle)[:μ]
coef(results_log_mle)[:Wf]

stderror(results_log_mle)
# 4-element Named Vector{Float64}
# A   │ 
# ────┼───────────
# σ_W │  0.0501555
# W₀  │   0.144674
# μ   │ 0.00385156
# Wf  │   0.127587


W₀_opt = coef(results_log_mle)[:W₀]
μ_opt = coef(results_log_mle)[:μ]
Wf_opt = coef(results_log_mle)[:Wf]


u₀_opt_log = [:W => W₀_opt]
tspan = (0.0, 100.0)
params_opt_log = [:μ => μ_opt, :Wf => Wf_opt]
oprob_opt_log = ODEProblem(growth_mod_log, u₀_opt_log, tspan, params_opt_log)
osol_opt_log = solve(oprob_opt_log, Tsit5(), saveat=0.5)
plot(osol_opt_log, label="Logistic growth", xlabel="t", xlims=(0, 100), ylims=(0, 14))
scatter!(t_meas, W_meas, label="Yield")


##########################################################################
# MAP = Maximum A Posterior
##########################################################################

results_log_map = optimize(growth_log_inf, MAP(), NelderMead())
# ModeResult with maximized lp of -9.57
# [0.33313095709980506, 1.8506166079656807, 0.07634774217311349, 11.83633424737603]
results_log_map |> coeftable
# ──────────────────────────────────────────────────────────────────────────
#           Coef.  Std. Error         z     Pr(>|z|)   Lower 95%   Upper 95%
# ──────────────────────────────────────────────────────────────────────────
# σ_W   0.333131   0.0507983    6.55792  5.45629e-11   0.233568    0.432694
# W₀    1.85062    0.146758    12.61     1.86054e-36   1.56298     2.13826
# μ     0.0763477  0.00391111  19.5207   7.31648e-85   0.0686821   0.0840134
# Wf   11.8363     0.130578    90.6455   0.0          11.5804     12.0923
# ──────────────────────────────────────────────────────────────────────────


##########################################################################
# NUTS = No-U-Turn Sampler
##########################################################################

iterations = 100;
results_log_nuts = sample(growth_log_inf, NUTS(), iterations)
summarize(results_log_nuts)
# parameters      mean       std      mcse   ess_bulk   ess_tail      rhat   ess_per_sec 
# Symbol   Float64   Float64   Float64    Float64    Float64   Float64       Float64 
#    σ_W    0.3865    0.0797    0.0192    22.1543    21.0246    1.0401        2.7958
#     W₀    1.8196    0.1601    0.0295    31.4596    24.3582    1.0285        3.9702
#      μ    0.0774    0.0050    0.0010    30.2795    31.7840    1.0258        3.8212
#     Wf   11.8245    0.1341    0.0236    32.1964    58.9391    1.0227        4.0632

iterations = 500;
results_log_nuts = sample(growth_log_inf, NUTS(), iterations)
summarize(results_log_nuts)
# parameters      mean       std      mcse   ess_bulk   ess_tail      rhat   ess_per_sec 
# Symbol   Float64   Float64   Float64    Float64    Float64   Float64       Float64 
#    σ_W    0.3803    0.0680    0.0044   222.5733   147.6067    0.9983       29.8476
#     W₀    1.8623    0.1690    0.0092   330.4091   314.9978    1.0024       44.3086
#      μ    0.0762    0.0043    0.0003   269.5818   333.2991    1.0069       36.1515
#     Wf   11.8442    0.1519    0.0137   123.0974   202.3017    1.0192       16.5076

