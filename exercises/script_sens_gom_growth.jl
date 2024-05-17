using Catalyst
using DifferentialEquations, Plots


################################################################
## Modelling Gompertz growth
################################################################

growth_gom = @reaction_network begin
    -μ, W --> ∅
    D*log(W), W --> ∅
end

parameters(growth_gom)
# μ
# D

osys_gom  = convert(ODESystem, growth_gom)

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
oprob_gom = ODEProblem(growth_gom, u₀_gom, tspan, params_gom)

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
oprob_uncert_gom = ODEProblem(growth_gom, u₀_gom, tspan, params_uncert_gom)

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
    oprob = ODEProblem(growth_gom, u₀, tspan, params, combinatoric_ratelaws=false)
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
oprob_opt_gom = ODEProblem(growth_gom, u₀_opt_gom, tspan, params_opt_gom)
osol_opt_gom = solve(oprob_opt_gom, Tsit5(), saveat=0.5)
plot(osol_opt_gom, label="Gompertz growth", xlabel="t", xlims=(0, 80), ylims=(0, 10))
scatter!(t_meas, W_meas, label="Yield")



################################################################
## Calibration 2
################################################################

using Optim, Turing, StatsPlots
using LinearAlgebra
using StatsBase

W_meas = [1.87, 2.45, 3.72, 4.32, 5.28, 7.01, 6.83, 8.62, 9.45, 10.31, 10.56, 11.72, 11.05, 11.53, 11.39, 11.7, 11.15, 11.49, 12.04, 11.95, 11.68]
t_meas = 0:5:100
scatter(t_meas, W_meas, title="Grass growth data",
    label="Yield", xlabel="t", xlims=(0, 100), ylims=(0, 14))

##########################################################################
# Setting up the Turing model
##########################################################################

u₀_gom = [:W => 2.0]
params_gom = [:μ => 0.09, :D => 0.04]
oprob_gom = ODEProblem(growth_gom, u₀_gom, tspan, params_gom)



@model function growth_gom_inference(t_meas, W)
    σ_W ~ InverseGamma()
    W₀ ~ Uniform(0, 10)
    μ ~ Uniform(0, 2)
    D ~ Uniform(0, 1)
    osol_gom = solve(remake(oprob_gom; u0=[W₀]), Tsit5(), saveat=t_meas, p=[μ, D])
    W ~ MvNormal(osol_gom[:W], σ_W^2 * I)
end

growth_gom_inf = growth_gom_inference(t_meas, W_meas)


##########################################################################
# MLE = Maximum Likelihood Estimation
##########################################################################

results_gom_mle = optimize(growth_gom_inf, MLE(), NelderMead())
# ModeResult with maximized lp of -10.34
# [0.39598609421047243, 1.4367305292920083, 0.1305278175841566, 0.052317705283294276]
results_gom_mle |> coeftable
# ───────────────────────────────────────────────────────────────────────
#          Coef.  Std. Error         z     Pr(>|z|)  Lower 95%  Upper 95%
# ───────────────────────────────────────────────────────────────────────
# σ_W  0.395986   0.0611004    6.4809   9.11748e-11   0.276231  0.515741
# W₀   1.43673    0.203328     7.06606  1.5939e-12    1.03821   1.83525
# μ    0.130528   0.00779095  16.7538   5.31524e-63   0.115258  0.145798
# D    0.0523177  0.00334683  15.632    4.40485e-55   0.045758  0.0588774
# ───────────────────────────────────────────────────────────────────────

coefnames(results_gom_mle)
# 4-element Vector{Symbol}:
#  :σ_W
#  :W₀
#  :μ
#  :D

coef(results_gom_mle)
# 4-element Named Vector{Float64}
# A   │ 
# ────┼──────────
# σ_W │  0.395986
# W₀  │   1.43673
# μ   │  0.130528
# D   │ 0.0523177

stderror(results_gom_mle)
# 4-element Named Vector{Float64}
# A   │ 
# ────┼───────────
# σ_W │  0.0611004
# W₀  │   0.203328
# μ   │ 0.00779095
# D   │ 0.00334683

coef(results_gom_mle)[:σ_W]
coef(results_gom_mle)[:W₀]
coef(results_gom_mle)[:μ]
coef(results_gom_mle)[:D]

W₀_opt_gom = coef(results_gom_mle)[:W₀]
μ_opt_gom = coef(results_gom_mle)[:μ]
D_opt_gom = coef(results_gom_mle)[:D]

u₀_opt_gom = [:W => W₀_opt_gom]
tspan = (0.0, 100.0)
params_opt_gom = [:μ => μ_opt_gom, :D => D_opt_gom]
oprob_opt_gom = ODEProblem(growth_gom, u₀_opt_gom, tspan, params_opt_gom)
osol_opt_gom = solve(oprob_opt_gom, Tsit5(), saveat=0.5)
plot(osol_opt_gom, label="Gompertz growth", xlabel="t", xlims=(0, 100), ylims=(0, 14))
scatter!(t_meas, W_meas, label="Yield")


##########################################################################
# MAP = Maximum A Posterior
##########################################################################

results_gom_map = optimize(growth_gom_inf, MAP(), NelderMead())
# ModeResult with maximized lp of -14.01
# [0.40074571995723124, 1.4367182028274699, 0.13052775236058195, 0.05231763826052395]
results_gom_map |> coeftable
# ───────────────────────────────────────────────────────────────────────
#          Coef.  Std. Error         z     Pr(>|z|)  Lower 95%  Upper 95%
# ───────────────────────────────────────────────────────────────────────
# σ_W  0.400746   0.0607573    6.59584  4.22844e-11  0.281664   0.519828
# W₀   1.43672    0.205771     6.98213  2.90744e-12  1.03341    1.84002
# μ    0.130528   0.00788457  16.5548   1.47733e-61  0.115074   0.145981
# D    0.0523176  0.00338704  15.4464   7.97348e-54  0.0456792  0.0589561
# ───────────────────────────────────────────────────────────────────────

##########################################################################
# NUTS = No-U-Turn Sampler
##########################################################################

iterations = 500;
results_gom_nuts = sample(growth_gom_inf, NUTS(), iterations)
summarize(results_gom_nuts)
# parameters      mean       std      mcse   ess_bulk   ess_tail      rhat   ess_per_sec 
# Symbol   Float64   Float64   Float64    Float64    Float64   Float64       Float64 
#    σ_W    0.4524    0.0756    0.0054   193.4209   173.1813    1.0050      130.6898
#     W₀    1.4270    0.2275    0.0213   114.7213   126.2865    1.0202       77.5144
#      μ    0.1317    0.0092    0.0010    90.4793   152.7701    1.0170       61.1347
#      D    0.0528    0.0039    0.0004    86.4987   160.6863    1.0180       58.4451

iterations = 1000;
results_gom_nuts = sample(growth_gom_inf, NUTS(), iterations)
summarize(results_gom_nuts)
# parameters      mean       std      mcse   ess_bulk   ess_tail      rhat   ess_per_sec 
# Symbol   Float64   Float64   Float64    Float64    Float64   Float64       Float64 
#    σ_W    0.4620    0.0796    0.0047   274.0115   312.3791    1.0032      119.4470
#     W₀    1.4270    0.2374    0.0163   206.3392   410.0721    1.0011       89.9473
#      μ    0.1316    0.0091    0.0006   209.7804   329.0615    1.0011       91.4474
#      D    0.0528    0.0039    0.0003   203.4118   369.4011    1.0017       88.6712