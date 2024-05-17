#======================================================================================
=======================================================================================
  Fermenter with monod kinetics
=======================================================================================
======================================================================================#

#=
https://docs.sciml.ai/Catalyst/stable/introduction_to_catalyst/catalyst_for_new_julia_users/
https://docs.sciml.ai/Catalyst/stable/catalyst_functionality/dsl_description/
https://www.matecdev.com/posts/julia-plotting-multiple-plots.html
=#

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
# 2-element Vector{Equation}:
#  Differential(t)(S(t)) ~ (Q*Sin) / V + (-Q*S(t)) / V - X(t)*S(t)*Catalyst.mm(S(t), μmax, Ks)
#  Differential(t)(X(t)) ~ (-Q*X(t)) / V + Y*X(t)*S(t)*Catalyst.mm(S(t), μmax, Ks)
# Check out the states
states(osys)
# 2-element Vector{SymbolicUtils.BasicSymbolic{Real}}:
#  S(t)
#  X(t)
# Check out the parameters
parameters(osys)
# 6-element Vector{SymbolicUtils.BasicSymbolic{Real}}:
#  μmax
#  Ks
#  Y
#  Q
#  V
#  Sin

# Set initial condition for ODE problem
u₀ = [:S => 0.0, :X => 0.01]
# Set time span for simulation
tspan = (0.0, 200)
# Set parameters for ODE problem
params = [:μmax => 0.30, :Ks => 0.15, :Y => 0.80, :Q => 2, :V => 40, :Sin => 2.2]
# Create ODE problem
oprob = ODEProblem(fermenter_monod, u₀, tspan, params)



osol1 = solve(oprob, Tsit5(), saveat=0.5)
plot(osol1)

#  Finding steady state
osol1[:S][end]
osol1[:X][end]
u_guess = [:S => osol1[:S][end], :X => osol1[:X][end]]
Seq, Xeq = solve(SteadyStateProblem(ODEProblem(fermenter_monod, u_guess, tspan, params, combinatoric_ratelaws=false)))
# 0.3093512042040855
# 1.5125190366367316


# Double the Sin at t=100
condition2 = [100.0]
function affect2!(integrator)
    integrator.p[6] = 2.8
end
cb2 = PresetTimeCallback(condition2, affect2!)
osol2 = solve(deepcopy(oprob), Tsit5(), saveat=0.5, callback=cb2)
plot(osol2)

#  Finding new steady state
osol2[:S][end]
osol2[:X][end]
u_guess = [:S => osol2[:S][end], :X => osol2[:X][end]]
params_mod = [:Y => 0.80, :Q => 2.0, :V => 40, :Sin => 2.8, :μmax => 0.30, :Ks => 0.15]
Seq, Xeq = solve(SteadyStateProblem(ODEProblem(fermenter_monod, u_guess, tspan, params_mod, combinatoric_ratelaws=false)))
# 0.30935120420397816
# 1.9925190366368177


# Double the flow Q at t=100
condition3 = [100.0]
function affect3!(integrator)
    integrator.p[4] *= 2.0
end
cb3 = PresetTimeCallback(condition3, affect3!)
osol3 = solve(deepcopy(oprob), Tsit5(), saveat=0.5, callback=cb3)
plot(osol3)

#  Finding new steady state
osol3[:S][end]
osol3[:X][end]
u_guess = [:S => osol3[:S][end], :X => osol3[:X][end]]
params = [:Y => 0.80, :Q => 4.0, :V => 40, :Sin => 2.2, :μmax => 0.30, :Ks => 0.15]
Seq, Xeq = solve(SteadyStateProblem(ODEProblem(fermenter_monod, u_guess, tspan, params, combinatoric_ratelaws=false)))
# 0.5337604031627765
# 1.3329916774697788


################################################################
## Sensitivity analysis
################################################################

using SciMLSensitivity

# u₀ = [:S => 0.0, :X => 0.01]
# # Set time span for simulation
# tspan = (0.0, 200)
# # Set parameters for ODE problem
# params = [:μmax => 0.30, :Ks => 0.15, :Y => 0.80, :Q => 2, :V => 40, :Sin => 2.2]

tspan = (0.0, 100)

oprob_sens = ODEForwardSensitivityProblem(oprob.f, [0.0, 0.01], tspan, [0.3, 0.15, 0.80, 2.0, 40.0, 2.2])
# oprob.f        : the ODE function
# [0.02, 5.0e-5] : initial conditions, the order is [S, X], see species(fermenter_monod)
# tspan          : the time span
# [0.39, 0.48, 0.77, 2, 40, 2.2] : the parameter values, the order is [μmax, Ks, Y, Q, V Sin], see parameters(fermenter_monod)

osol_sens = solve(oprob_sens, Tsit5(), saveat=0.5)
u, dp = extract_local_sensitivities(osol_sens)

sol_S = u[1,:]     # solution for S
sol_X = u[2,:]     # solution for X

sens_μmax = dp[1]'  # sensitivity for S and X on μmax  
sens_Ks   = dp[2]'  # sensitivity for S and X on Ks
sens_Sin  = dp[6]'  # sensitivity for S and X on Sin


# Absolute sensitivity of S and X to μmax
plot(osol_sens.t, sens_μmax, title="Absolute sensitivity of S and X to μmax", label=["dS/dμmax" "dX/dμmax"])

# Absolute sensitivity of S and X to Ks
plot(osol_sens.t, sens_Ks, title="Absolute sensitivity of S and X to Ks", label=["dS/dKs" "dX/dKs"])

# Absolute sensitivity of S and X to Sin
plot(osol_sens.t, sens_Sin, title="Absolute sensitivity for S and X on Sin", label=["dS/dSin" "dX/dSin"])


# trsens_μmax = sens_μmax.*0.3./[sol_S sol_X]
# trsens_Ks   = sens_Ks.*0.15./[sol_S sol_X]
# trsens_Sin  = sens_Sin.*2.2./[sol_S sol_X]

# # Total relative sensitivity of S and X to μmax
# plot(osol_sens.t, trsens_μmax, title="Total relative sensitivity of S and X to μmax", label=["dS/dμmax" "dX/dμmax"])

# # Total relative sensitivity of S and X to Ks
# plot(osol_sens.t, trsens_Ks, title="Total relative sensitivity of S and X to Ks", label=["dS/dKs" "dX/dKs"])

# # Total relative sensitivity of S and X to Sin
# plot(osol_sens.t, trsens_Sin, title="Total relative sensitivity for S and X on Ks", label=["dS/dSin" "dX/dSin"])


################################################################
## Parameter uncertainty
################################################################

using Measurements

# ± with \pm
params_uncert = [:μmax => 0.30, :Ks => 0.15, :Y => 0.80, :Q => 2, :V => 40, :Sin => 2.2±0.4]
params_uncert = [:μmax => 0.30, :Ks => 0.15, :Y => 0.80, :Q => 2, :V => 40±1, :Sin => 2.2]
params_uncert = [:μmax => 0.30, :Ks => 0.15, :Y => 0.80, :Q => 2±0.2, :V => 40, :Sin => 2.2]
params_uncert = [:μmax => 0.30, :Ks => 0.15±0.03, :Y => 0.80, :Q => 2, :V => 40, :Sin => 2.2]
params_uncert = [:μmax => 0.30±0.06, :Ks => 0.15, :Y => 0.80, :Q => 2, :V => 40, :Sin => 2.2]
params_uncert = [:μmax => 0.30±0.06, :Ks => 0.15±0.03, :Y => 0.80, :Q => 2, :V => 40, :Sin => 2.2±0.4]

oprob_uncert = ODEProblem(fermenter_monod, u₀, tspan, params_uncert, combinatoric_ratelaws=false)

osol_uncert = solve(oprob_uncert, Tsit5(), saveat=2.0)

plot(osol_uncert)


################################################################
## Calibration
################################################################

using Distributions
x = rand(Truncated(Normal(0, 1), -3, 3), 10)
x = rand(Normal(0, 1), 10)

# u₀ = [:S => 0.0, :X => 0.01]
# tspan = (0.0, 100)
# params = [:μmax => 0.26, :Ks => 0.19, :Y => 0.72, :Q => 2, :V => 40, :Sin => 2.2]
# oprob = ODEProblem(fermenter_monod, u₀, tspan, params)
# osol = solve(oprob, Tsit5(), saveat=2)
# # plot(osol)
# size(osol, 2)

# Ssim = osol[:S]
# Xsim = osol[:X]

# maximum(Ssim)

# sigma_S = 0.05
# sigma_X = 0.10

# ratio = 0.8
# Sdta = Ssim + rand(Normal(0, sigma_S), size(osol, 2)) .* ((1 - ratio) .+ ratio.*Ssim./maximum(Ssim))
# Xdta = Xsim + rand(Normal(0, sigma_X), size(osol, 2)) .* ((1 - ratio) .+ ratio.*Xsim./maximum(Xsim))

# scatter(osol.t, Sdta)
# scatter!(osol.t, Xdta)

# S_meas = Sdta
# X_meas = Xdta
# t_meas = osol.t

# print(S_meas)
# S_meas = [0.001748527898167188, 0.20428991317610207, 0.42308290211345584, 0.5812773896467612, 0.6760955310938963, 0.8730785271065923, 0.9502601158945785, 1.131196692751299, 1.228020843813872, 1.2428339237402242, 1.3409788577268795, 1.3855474603964635, 1.4387872770278027, 1.442784786201463, 1.450970270838476, 1.4814660520467369, 1.2270769217315953, 1.119240885364341, 1.0068863615158918, 0.8232574692493287, 0.8005116019872506, 0.6401634395227176, 0.5764656652440903, 0.5142843390157951, 0.48473334602595114, 0.46653570762816676, 0.45822303517525076, 0.4238518030171025, 0.4153083110486855, 0.43834474729111095, 0.4634342293342907, 0.37470678318422046, 0.45145216531521837, 0.4126438618131415, 0.4196912946692012, 0.40807048202868285, 0.38166341437318296, 0.46286939554611733, 0.39997029595499917, 0.3667422851420893, 0.3754225717634865, 0.4164556699035795, 0.381149816287548, 0.4083011083028045, 0.3823032605431713, 0.36564420671184467, 0.3782627236152731, 0.4155201254906501, 0.3751484002874023, 0.3861453303962444, 0.3931825702538002]
# print(X_meas)
# X_meas = [0.03540919867413992, 0.013675712660366962, 0.007226934232513934, 0.010107563084744696, 0.0036780743364808197, 0.03870270476467011, 0.07416649472147932, 0.003674806536739371, -0.003046359317726847, 0.01317026282226353, 0.042985358002208436, 0.060052865902412045, 0.0780362051764657, 0.10046926924661961, 0.18808852185040334, 0.18449433987203973, 0.3433540452298016, 0.5056554789903555, 0.5513529677789504, 0.6221808322430916, 0.8131884825991468, 0.9671402897884847, 0.832205439079949, 1.0512156552889849, 1.088611197626013, 1.1766034356773647, 1.2280158439777207, 1.1661520146082969, 1.2269064196887367, 1.1109293107027678, 1.0439927716992385, 1.1463662707389317, 1.2511189528458733, 1.269791213962837, 1.290139900627747, 1.1923591193123506, 1.3594548592488953, 1.2831852961253911, 1.0973447130065035, 1.215967277395141, 1.2304130840739567, 1.250198933272166, 1.2270270845206828, 1.0700501615028046, 1.2058227150414467, 1.4043418567038835, 1.3112939081244053, 1.3776685598728324, 1.2634137570689783, 1.330686413825863, 1.2408185373970535]

# S_meas = round.(S_meas, digits=3)
# X_meas = round.(X_meas, digits=3)

# print(S_meas)
# print(X_meas)

S_meas = [0.002, 0.204, 0.423, 0.581, 0.676, 0.873, 0.95, 1.131, 1.228, 1.243, 1.341, 1.386, 1.439, 1.443, 1.451, 1.481, 1.227, 1.119, 1.007, 0.823, 0.801, 0.64, 0.576, 0.514, 0.485, 0.467, 0.458, 0.424, 0.415, 0.438, 0.463, 0.375, 0.451, 0.413, 0.42, 0.408, 0.382, 0.463, 0.4, 0.367, 0.375, 0.416, 0.381, 0.408, 0.382, 0.366, 0.378, 0.416, 0.375, 0.386, 0.393]
X_meas = [0.035, 0.014, 0.007, 0.01, 0.004, 0.039, 0.074, 0.004, 0.000, 0.013, 0.043, 0.06, 0.078, 0.1, 0.188, 0.184, 0.343, 0.506, 0.551, 0.622, 0.813, 0.967, 0.832, 1.051, 1.089, 1.177, 1.228, 1.166, 1.227, 1.111, 1.044, 1.146, 1.251, 1.27, 1.29, 1.192, 1.359, 1.283, 1.097, 1.216, 1.23, 1.25, 1.227, 1.07, 1.206, 1.404, 1.311, 1.378, 1.263, 1.331, 1.241]
t_meas = 0.0:2.0:100.0

scatter(t_meas, S_meas, color=:blue)
scatter!(t_meas, X_meas, color=:red)

S_sigma = 0.05
X_sigma = 0.10

u₀ = [:S => 0.0, :X => 0.01]
tspan = (0.0, 100)

function Jtheta_fermenter_monod(thetas)           # function return to be minimized
    params = [:μmax => thetas[1], :Ks => thetas[2], :Y => thetas[3], :Q => 2, :V => 40, :Sin => 2.2]
    oprob = ODEProblem(fermenter_monod, u₀, tspan, params, combinatoric_ratelaws=false)
    sol = solve(oprob, Tsit5(), saveat=t_meas)
    S_sol = sol[:S]
    X_sol = sol[:X]
    J = (1/S_sigma^2)*sum(abs2, S_sol - S_meas) + (1/X_sigma^2)*sum(abs2, X_sol - X_meas)
    return J
end

thetas_init = [0.30, 0.15, 0.80]
Jtheta_fermenter_monod(thetas_init)

using Optim

results = optimize(Jtheta_fermenter_monod, thetas_init, NelderMead())
thetas_opt = Optim.minimizer(results)
Optim.minimum(results)    # J-value / mse-value

params = [:μmax => thetas_opt[1], :Ks => thetas_opt[2], :Y => thetas_opt[3], :Q => 2, :V => 40, :Sin => 2.2]
oprob = ODEProblem(fermenter_monod, u₀, tspan, params, combinatoric_ratelaws=false)
osol = solve(oprob, Tsit5(), saveat=t_meas)

plot(osol, labels=["S sim" "X sim"], xlabel="t")   # , xlims=(0, 80), ylims=(0, 10)
scatter!(t_meas, S_meas, label="S meas", color=:blue)
scatter!(t_meas, X_meas, label="X meas", color=:red)


################################################################
## Calibration 2
################################################################

using Optim, Turing, StatsPlots
using LinearAlgebra
using StatsBase

u₀ = [:S => 0.0, :X => 0.01]
tspan = (0.0, 100)
params = [:μmax => 0.30, :Ks => 0.15, :Y => 0.80, :Q => 2, :V => 40, :Sin => 2.2]
oprob = ODEProblem(fermenter_monod, u₀, tspan, params)

osol = solve(oprob, Tsit5(), saveat=5)
plot(osol)

σ_S, σ_X = 0.02, 0.04

N = length(osol)

# Sobs = abs.(osol[:S] .+ σ_S .* randn(length(osol)))
# Xobs = abs.(osol[:X] .+ σ_X .* randn(length(osol)))

let
	scatter(osol.t, Sobs)
	scatter!(osol.t, Xobs)
end

println(round.(Sobs, digits=4))
println(round.(Xobs, digits=4))
S_meas = [0.0133, 0.5252, 0.864, 1.1218, 1.2865, 1.2405, 0.8459, 0.5288, 0.4028, 0.393, 0.325, 0.3536, 0.3346, 0.2992, 0.2944, 0.3495, 0.3236, 0.3259, 0.3147, 0.3216, 0.3021]
X_meas = [0.0673, 0.0283, 0.0219, 0.0897, 0.1282, 0.2227, 0.6763, 0.9794, 1.201, 1.18, 1.3214, 1.4004, 1.4115, 1.3927, 1.4711, 1.5209, 1.4962, 1.5648, 1.4552, 1.4781, 1.5051]
t_meas = 0:5:100
# tsteps = osol.t

@model function biomass_inference(t_meas, S, X)
	σ_S ~ InverseGamma()
	σ_X ~ InverseGamma()
  μmax ~ LogNormal()
	Ks ~ LogNormal()
  # IMPORTANT: place all parameters in the option p=[...] in the same order
  # as in the reaction network model.
	osol = solve(oprob, Tsit5(), saveat=t_meas, p=[μmax, Ks, 0.80, 2, 40, 2.2])
	S ~ MvNormal(osol[:S], σ_S^2 * I)
	X ~ MvNormal(osol[:X], σ_X^2 * I)
end

biomass_mod = biomass_inference(t_meas, S_meas, X_meas)

#=
Information about MLE and MAP:
https://towardsdatascience.com/mle-vs-map-a989f423ae5c
https://roger010620.medium.com/maximum-likelihood-estimation-mle-and-maximum-a-posterior-map-in-machine-learning-923e7bc743b2
=#


##########################################################################
# MLE = Maximum Likelihood Estimation
##########################################################################

results_mle = optimize(biomass_mod, MLE(), NelderMead())
# ModeResult with maximized lp of 89.32
# [0.018749286626397788, 0.04439754303721986, 0.3063384375324556, 0.17030119660821036]
results_mle |> coeftable
# ────────────────────────────────────────────────────────────────────────
#           Coef.  Std. Error         z     Pr(>|z|)  Lower 95%  Upper 95%
# ────────────────────────────────────────────────────────────────────────
# σ_S   0.0187493  0.00290906   6.44514  1.15493e-10  0.0130476  0.0244509
# σ_X   0.0443975  0.00688693   6.44664  1.1436e-10   0.0308994  0.0578957
# μmax  0.306338   0.00593128  51.648    0.0          0.294713   0.317964
# Ks    0.170301   0.0186784    9.11754  7.68472e-20  0.133692   0.20691
# ────────────────────────────────────────────────────────────────────────


##########################################################################
# MAP = Maximum A Posterior
##########################################################################

results_map = optimize(biomass_mod, MAP(), NelderMead())
# ModeResult with maximized lp of 50.36
# [0.05005592910567663, 0.06922016251253287, 0.3070229211180945, 0.17563313150439827]
results_map |> coeftable
# ─────────────────────────────────────────────────────────────────────────
#           Coef.  Std. Error         z      Pr(>|z|)  Lower 95%  Upper 95%
# ─────────────────────────────────────────────────────────────────────────
# σ_S   0.0500559  0.00986625   5.07345  3.90666e-7    0.0307184  0.0693934
# σ_X   0.0692202  0.0123952    5.58444  2.34455e-8    0.044926   0.0935143
# μmax  0.307023   0.0141115   21.757    5.93219e-105  0.279365   0.334681
# Ks    0.175633   0.0445641    3.94114  8.10961e-5    0.0882892  0.262977
# ─────────────────────────────────────────────────────────────────────────


##########################################################################
# NUTS = No-U-Turn Sampler
##########################################################################

iterations = 5000;
chain_biomass_nuts = sample(biomass_mod, NUTS(), iterations)    # MCMCThreads()
# ┌ Info: Found initial step size
# └   ϵ = 0.4
# Chains MCMC chain (5000×16×1 Array{Float64, 3}):

# Iterations        = 1001:1:6000
# Number of chains  = 1
# Samples per chain = 5000
# Wall duration     = 17.43 seconds
# Compute duration  = 17.43 seconds
# parameters        = σ_S, σ_X, μmax, Ks
# internals         = lp, n_steps, is_accept, acceptance_rate, log_density, hamiltonian_energy, hamiltonian_energy_error, max_hamiltonian_energy_error, tree_depth, numerical_error, step_size, nom_step_size
# Summary Statistics
#   parameters      mean       std      mcse    ess_bulk    ess_tail      rhat   ess_per_sec 
#       Symbol   Float64   Float64   Float64     Float64     Float64   Float64       Float64 
#          σ_S    0.0582    0.0129    0.0002   3392.6984   2741.4693    1.0000      194.6024
#          σ_X    0.0764    0.0149    0.0003   3760.6838   2805.2152    1.0006      215.7097
#         μmax    0.3132    0.0171    0.0003   2724.8752   2708.5298    0.9999      156.2966
#           Ks    0.1962    0.0548    0.0011   2514.5185   2573.0328    1.0003      144.2307
# Quantiles
#   parameters      2.5%     25.0%     50.0%     75.0%     97.5% 
#       Symbol   Float64   Float64   Float64   Float64   Float64 
#          σ_S    0.0390    0.0491    0.0562    0.0653    0.0883
#          σ_X    0.0528    0.0659    0.0745    0.0847    0.1110
#         μmax    0.2849    0.3012    0.3114    0.3230    0.3525
#           Ks    0.1085    0.1574    0.1903    0.2273    0.3189

summarize(chain_biomass_nuts)
# parameters      mean       std      mcse    ess_bulk    ess_tail      rhat   ess_per_sec 
# Symbol   Float64   Float64   Float64     Float64     Float64   Float64       Float64 
#    σ_S    0.0582    0.0129    0.0002   3392.6984   2741.4693    1.0000      194.6024
#    σ_X    0.0764    0.0149    0.0003   3760.6838   2805.2152    1.0006      215.7097
#   μmax    0.3132    0.0171    0.0003   2724.8752   2708.5298    0.9999      156.2966
#     Ks    0.1962    0.0548    0.0011   2514.5185   2573.0328    1.0003      144.2307
