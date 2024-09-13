
using Catalyst
using DifferentialEquations, Plots

#======================================================================================
  Solving using Catalyst
======================================================================================#

irrigation_mod = @reaction_network begin
    k/Smax, S1 --> S2
    v, 2S2 --> ∅
    R * (1 - S1res / Smax), ∅ --> S1
    R/Smax, S1 --> ∅
end

parameters(irrigation_mod)
# k
# Smax
# v
# R
# S1res

# Convert to an ODE system
osys  = convert(ODESystem, irrigation_mod)

# Check out the diff. eqns. equations
equations(osys)
# Differential(t)(S1(t)) ~ (-k*S1(t)) / Smax + (-R*S1(t)) / Smax + R*(1 + (-S1res) / Smax)
# Differential(t)(S2(t)) ~ (k*S1(t)) / Smax - v*(S2(t)^2)
# Check out the states
states(osys)
# S1(t)
# S2(t)
# Check out the parameters
parameters(osys)
# k
# Smax
# v
# R
# S1res

S1res = 10
S2res = 15
S10 = 40 - S1res
S20 = 40 - S2res
tstep = 60
R = 5
Smax = 150
k = 3
v = 1e-3

# Set initial conditions
u0 = [:S1 => S10, :S2 => S20]
u0 = [:S1 => 0.0, :S2 => 0.0]
u0 = [:S1 => 140.0, :S2 => 135.0]
# Set time span for simulation
tspan = (0.0, 150.0)  # the time interval to solve on
# Set parameters for ODE problem
params = [:k => 3.0, :Smax => 150.0, :v => 1e-3, :R => 5.0, :S1res => 10.0]
# Create ODE problem
oprob = ODEProblem(irrigation_mod, u0, tspan, params)
# osol = solve(oprob, Tsit5(), saveat=0.5)
# plot(osol)
# osol[:S1]
# osol[:S2]

condition = [tstep]
function affect!(integrator)
    integrator.p[4] += 5.0      # R is the 4th parameter !!!
end
cb = PresetTimeCallback(condition, affect!)
osol = solve(deepcopy(oprob), Tsit5(), saveat=10.0, callback=cb)
plot(osol)


osol.t
osol[:S1]
osol[:S2]

plot(osol.t, (osol[:S1] .+ S1res) ./ Smax .* R, xaxis="time [s]", label="runoff")
plot!(osol.t, (osol[:S2].^2) .* v, label="outflow")
plot!(osol.t, k .* osol[:S1] ./ Smax, label="percolation")


################################################################
## Calibration
################################################################

# # u0 = [:S1 => 0.0, :S2 => 0.0]
# u0 = [:S1 => 140.0, :S2 => 135.0]
# tspan = (0.0, 150.0)  # the time interval to solve on
# params = [:k => 4.2, :Smax => 143.0, :v => 1e-3, :R => 5.0, :S1res => 10.0]
# oprob = ODEProblem(irrigation_mod, u0, tspan, params)
# osol = solve(oprob, Tsit5(), saveat=10.0)
# plot(osol)

# S1sim = osol[:S1]
# S2sim = osol[:S2]

# maximum(S1sim)
# maximum(S2sim)

# S1_sigma = 4.0
# S2_sigma = 6.0

# using Distributions
# ratio = 0.8
# S1dta = S1sim + rand(Normal(0, S1_sigma), size(osol, 2)) .* ((1 - ratio) .+ ratio.*S1sim./maximum(S1sim))
# S2dta = S2sim + rand(Normal(0, S2_sigma), size(osol, 2)) .* ((1 - ratio) .+ ratio.*S2sim./maximum(S2sim))

# begin
#     scatter(osol.t, S1dta)
#     scatter!(osol.t, S2dta)
# end

# println(round.(S1dta, digits=2))
# println(round.(S2dta, digits=2))

# S1_meas1 = round.(S1dta, digits=2)
# S2_meas1 = round.(S2dta, digits=2)

# S1_meas2 = round.(S1dta, digits=2)
# S2_meas2 = round.(S2dta, digits=2)

S1_meas1 = [0.2, 35.94, 52.49, 66.86, 60.66, 67.81, 73.22, 71.31, 72.94, 64.08, 70.11, 68.53, 70.54, 63.63, 67.39, 62.84]
S2_meas1 = [0.63, 6.2, 17.67, 22.96, 35.41, 44.08, 43.5, 53.34, 47.57, 47.77, 43.96, 52.22, 46.67, 46.74, 46.46, 39.92]

S1_meas2 = [137.96, 106.15, 90.15, 84.64, 76.15, 75.73, 73.32, 68.48, 70.06, 69.36, 70.91, 72.13, 76.25, 74.34, 74.93, 71.58]
S2_meas2 = [124.08, 80.14, 60.15, 50.12, 49.66, 47.78, 46.56, 48.41, 42.7, 43.72, 49.03, 51.91, 48.24, 46.14, 51.22, 43.78]

t_meas = 0:10:150

S1_sigma = 4.0
S2_sigma = 6.0

begin
  scatter(t_meas, S1_meas1, label="S1 meas", color=:blue, title="Experment 1")
  scatter!(t_meas, S2_meas1, label="S2 meas", color=:red, ylims=(0, 150))
end

begin
  scatter(t_meas, S1_meas2, label="S1 meas", color=:blue, title="Experment 2")
  scatter!(t_meas, S2_meas2, label="S2 meas", color=:red, ylims=(0, 150))
end

tspan = (0.0, 150.0)

function Jtheta_irrigation(thetas)           # function return to be minimized
  params = [:k => thetas[1], :Smax => thetas[2], :v => 1e-3, :R => 5, :S1res => 10]
  # Experiment 1
  u0 = [0.0, 0.0]
  oprob = ODEProblem(irrigation_mod, u0, tspan, params)
  osol = solve(oprob, Tsit5(), saveat=t_meas)
  S1_sol = osol[:S1]
  S2_sol = osol[:S2]
  J1 = (1/S1_sigma^2)*sum(abs2, S1_sol - S1_meas1) + (1/S2_sigma^2)*sum(abs2, S2_sol - S2_meas1)
  # Experiment 2
  u0 = [140.0, 135.0]
  oprob = ODEProblem(irrigation_mod, u0, tspan, params)
  osol = solve(oprob, Tsit5(), saveat=t_meas)
  S1_sol = osol[:S1]
  S2_sol = osol[:S2]
  J2 = (1/S1_sigma^2)*sum(abs2, S1_sol - S1_meas2) + (1/S2_sigma^2)*sum(abs2, S2_sol - S2_meas2)
  return J1+J2
end

thetas_init = [3.0, 150.0]
Jtheta_irrigation(thetas_init)

using Optim

results = optimize(Jtheta_irrigation, thetas_init, NelderMead())
thetas_opt = Optim.minimizer(results)
Optim.minimum(results)    # J-value / mse-value

params = [:k => thetas_opt[1], :Smax => thetas_opt[2], :v => 1e-3, :R => 5, :S1res => 10]

# Experiment 1
u0₁ = [0.0, 0.0]
oprob1 = ODEProblem(irrigation_mod, u0₁, tspan, params)
osol1 = solve(oprob1, Tsit5(), saveat=1.0)

begin
  plot(osol1, labels=["S1 sim" "S2 sim"], xlabel="t")   # , xlims=(0, 80), ylims=(0, 10)
  scatter!(t_meas, S1_meas1, label="S1 meas", color=:blue, title="Experment 1")
  scatter!(t_meas, S2_meas1, label="S2 meas", color=:red, ylims=(0, 150))
end


# Experiment 2
u0₂ = [140.0, 135.0]
oprob2 = ODEProblem(irrigation_mod, u0₂, tspan, params)
osol2 = solve(oprob2, Tsit5(), saveat=1.0)

begin
  plot(osol2, labels=["S1 sim" "S2 sim"], xlabel="t")   # , xlims=(0, 80), ylims=(0, 10)
  scatter!(t_meas, S1_meas2, label="S1 meas", color=:blue, title="Experment 2")
  scatter!(t_meas, S2_meas2, label="S2 meas", color=:red, ylims=(0, 150))
end


################################################################
## Calibration with Turing
################################################################

using Optim, Turing, StatsPlots
using LinearAlgebra
using StatsBase


S1_meas1 = [0.2, 35.94, 52.49, 66.86, 60.66, 67.81, 73.22, 71.31, 72.94, 64.08, 70.11, 68.53, 70.54, 63.63, 67.39, 62.84]
S2_meas1 = [0.63, 6.2, 17.67, 22.96, 35.41, 44.08, 43.5, 53.34, 47.57, 47.77, 43.96, 52.22, 46.67, 46.74, 46.46, 39.92]

S1_meas2 = [137.96, 106.15, 90.15, 84.64, 76.15, 75.73, 73.32, 68.48, 70.06, 69.36, 70.91, 72.13, 76.25, 74.34, 74.93, 71.58]
S2_meas2 = [124.08, 80.14, 60.15, 50.12, 49.66, 47.78, 46.56, 48.41, 42.7, 43.72, 49.03, 51.91, 48.24, 46.14, 51.22, 43.78]

t_meas = 0.0:10.0:150.0

u01 = [:S1 => 0.0, :S2 => 0.0]
# u02 = [:S1 => 140.0, :S2 => 135.0]
tspan = (0.0, 150.0)  # the time interval to solve on
params = [:k => 3.0, :Smax => 150.0, :v => 1e-3, :R => 5.0, :S1res => 10.0]
oprob = ODEProblem(irrigation_mod, u01, tspan, params)

@model function irrigation_inference(t_meas, S1_meas1, S2_meas1, S1_meas2, S2_meas2)
	σ_S1 ~ InverseGamma()
	σ_S2 ~ InverseGamma()
  k ~ Uniform(0, 10)
	Smax ~ Uniform(100, 200)
	osol1 = solve(remake(oprob; u0=[0.0, 0.0]), Tsit5(), saveat=t_meas, p=[k, Smax, 1e-3, 5, 10])
	S1_meas1 ~ MvNormal(osol1[:S1], σ_S1^2 * I)
	S2_meas1 ~ MvNormal(osol1[:S2], σ_S2^2 * I)
	osol2 = solve(remake(oprob; u0=[140.0, 135.0]), Tsit5(), saveat=t_meas, p=[k, Smax, 1e-3, 5, 10])
	S1_meas2 ~ MvNormal(osol2[:S1], σ_S1^2 * I)
	S2_meas2 ~ MvNormal(osol2[:S2], σ_S2^2 * I)
end

irrigation_inf = irrigation_inference(t_meas, S1_meas1, S2_meas1, S1_meas2, S2_meas2)


##########################################################################
# MLE = Maximum Likelihood Estimation
##########################################################################

results_mle = optimize(irrigation_inf, MLE(), NelderMead())
# ModeResult with maximized lp of -172.11
# [3.3732684275881772, 3.761249262218508, 4.387686310921089, 142.78702446781438]
results_mle |> coeftable
# ────────────────────────────────────────────────────────────────────────
#           Coef.  Std. Error         z     Pr(>|z|)  Lower 95%  Upper 95%
# ────────────────────────────────────────────────────────────────────────
# σ_S1    3.37327    0.421682   7.99955  1.24871e-15    2.54679    4.19975
# σ_S2    3.76125    0.470151   8.00009  1.24325e-15    2.83977    4.68273
# k       4.38769    0.23363   18.7805   1.09019e-78    3.92978    4.84559
# Smax  142.787      3.53529   40.3891   0.0          135.858    149.716
# ────────────────────────────────────────────────────────────────────────

coefnames(results_mle)
# 4-element Vector{Symbol}:
#  :σ_S1
#  :σ_S2
#  :k
#  :Smax

coef(results_mle)
# 4-element Named Vector{Float64}
# A    │ 
# ─────┼────────
# σ_S1 │ 3.37327
# σ_S2 │ 3.76125
# k    │ 4.38769
# Smax │ 142.787

stderror(results_mle)
# 4-element Named Vector{Float64}
# A    │ 
# ─────┼─────────
# σ_S1 │ 0.421682
# σ_S2 │ 0.470151
# k    │  0.23363
# Smax │  3.53529

k_opt = coef(results_mle)[:k]
Smax_opt = coef(results_mle)[:Smax]

tspan = (0.0, 150.0)  # the time interval to solve on
params_opt = [:k => k_opt, :Smax => Smax_opt, :v => 1e-3, :R => 5.0, :S1res => 10.0]

u01 = [:S1 => 0.0, :S2 => 0.0]
oprob1_opt = ODEProblem(irrigation_mod, u01, tspan, params_opt)
osol1_opt = solve(oprob1_opt, Tsit5(), saveat=0.5)

begin
  plot(osol1_opt, labels=["S1 sim" "S2 sim"], xlabel="t")   # , xlims=(0, 80), ylims=(0, 10)
  scatter!(t_meas, S1_meas1, label="S1 meas1", color=:blue, title="Experment 1")
  scatter!(t_meas, S2_meas1, label="S2 meas1", color=:red)
end

u02 = [:S1 => 140.0, :S2 => 135.0]
oprob2_opt = ODEProblem(irrigation_mod, u02, tspan, params_opt)
osol2_opt = solve(oprob2_opt, Tsit5(), saveat=0.5)

begin
  plot(osol2_opt, labels=["S1 sim" "S2 sim"], xlabel="t")   # , xlims=(0, 80), ylims=(0, 10)
  scatter!(t_meas, S1_meas2, label="S1 meas2", color=:blue, title="Experment 2")
  scatter!(t_meas, S2_meas2, label="S2 meas2", color=:red)
end


##########################################################################
# MAP = Maximum A Posterior
##########################################################################

results_map = optimize(irrigation_inf, MAP(), NelderMead())
# ModeResult with maximized lp of -184.62
# [3.2872658080787143, 3.6636698380367383, 4.387702230545071, 142.78723499931294]
results_map |> coeftable
# ────────────────────────────────────────────────────────────────────────
#           Coef.  Std. Error         z     Pr(>|z|)  Lower 95%  Upper 95%
# ────────────────────────────────────────────────────────────────────────
# σ_S1    3.28727    0.399552   8.22738  1.91354e-16    2.50416    4.07037
# σ_S2    3.66367    0.445172   8.22979  1.87548e-16    2.79115    4.53619
# k       4.3877     0.22757   19.2807   7.80715e-83    3.94167    4.83373
# Smax  142.787      3.44378   41.4623   0.0          136.038    149.537
# ────────────────────────────────────────────────────────────────────────


##########################################################################
# NUTS = No-U-Turn Sampler
##########################################################################

iterations = 5000;
chain_irrigation_nuts = sample(irrigation_inf, NUTS(), iterations)    # MCMCThreads()

# ┌ Info: Found initial step size
# └   ϵ = 0.0125
# Chains MCMC chain (5000×16×1 Array{Float64, 3}):
# Iterations        = 1001:1:6000
# Number of chains  = 1
# Samples per chain = 5000
# Wall duration     = 14.28 seconds
# Compute duration  = 14.28 seconds
# parameters        = σ_S1, σ_S2, k, Smax
# internals         = lp, n_steps, is_accept, acceptance_rate, log_density, hamiltonian_energy, hamiltonian_energy_error, max_hamiltonian_energy_error, tree_depth, numerical_error, step_size, nom_step_size
# Summary Statistics
#   parameters       mean       std      mcse    ess_bulk    ess_tail      rhat   ess_per_sec 
#       Symbol    Float64   Float64   Float64     Float64     Float64   Float64       Float64 
#         σ_S1     3.4809    0.4470    0.0080   3286.0019   2698.0620    0.9999      230.1605
#         σ_S2     3.8827    0.4982    0.0088   3588.6160   2739.9655    1.0000      251.3565
#            k     4.4152    0.2479    0.0048   2634.9634   2260.3591    1.0003      184.5600
#         Smax   143.1711    3.7492    0.0726   2718.0773   2218.6103    1.0001      190.3815
# Quantiles
#   parameters       2.5%      25.0%      50.0%      75.0%      97.5% 
#       Symbol    Float64    Float64    Float64    Float64    Float64 
#         σ_S1     2.7161     3.1734     3.4354     3.7520     4.4421
#         σ_S2     3.0606     3.5128     3.8315     4.1862     4.9797
#            k     3.9489     4.2414     4.4082     4.5799     4.9117
#         Smax   136.0869   140.5947   143.0800   145.5909   150.8676

summarize(chain_irrigation_nuts)
# parameters       mean       std      mcse    ess_bulk    ess_tail      rhat   ess_per_sec 
# Symbol    Float64   Float64   Float64     Float64     Float64   Float64       Float64 
#   σ_S1     3.4809    0.4470    0.0080   3286.0019   2698.0620    0.9999      230.1605
#   σ_S2     3.8827    0.4982    0.0088   3588.6160   2739.9655    1.0000      251.3565
#      k     4.4152    0.2479    0.0048   2634.9634   2260.3591    1.0003      184.5600
#   Smax   143.1711    3.7492    0.0726   2718.0773   2218.6103    1.0001      190.3815

mean(chain_irrigation_nuts)
# Mean
#   parameters       mean 
#       Symbol    Float64 
#         σ_S1     3.4809
#         σ_S2     3.8827
#            k     4.4152
#         Smax   143.1711

histogram(chain_irrigation_nuts[:k], bins=50)
histogram(chain_irrigation_nuts[:Smax], bins=50)

mean(chain_irrigation_nuts[:k])
std(chain_irrigation_nuts[:k])

mean(chain_irrigation_nuts[:Smax])
std(chain_irrigation_nuts[:Smax])
