
using Catalyst
using DifferentialEquations, Plots

#======================================================================================
  Solving using Catalyst
======================================================================================#

irrigation_mod = @reaction_network begin
    k/Smax, S₁ --> S₂
    v, 2S₂ --> ∅
    R * (1 - S₁res / Smax), ∅ --> S₁
    R/Smax, S₁ --> ∅
end

parameters(irrigation_mod)
# k
# Smax
# v
# R
# S₁res

# Convert to an ODE system
osys  = convert(ODESystem, irrigation_mod)

# Check out the diff. eqns. equations
equations(osys)
# Differential(t)(S₁(t)) ~ (-k*S₁(t)) / Smax + (-R*S₁(t)) / Smax + R*(1 + (-S₁res) / Smax)
# Differential(t)(S₂(t)) ~ (k*S₁(t)) / Smax - v*(S₂(t)^2)
# Check out the states
states(osys)
# S₁(t)
# S₂(t)
# Check out the parameters
parameters(osys)
# k
# Smax
# v
# R
# S₁res

S₁res = 10
S₂res = 15
S₁₀ = 40 - S₁res
S₂₀ = 40 - S₂res
tstep = 60
R = 5
Smax = 150
k = 3
v = 1e-3

# Set initial conditions
u₀ = [:S₁ => S₁₀, :S₂ => S₂₀]
u₀ = [:S₁ => 0.0, :S₂ => 0.0]
u₀ = [:S₁ => 140.0, :S₂ => 135.0]
# Set time span for simulation
tspan = (0.0, 150.0)  # the time interval to solve on
# Set parameters for ODE problem
params = [:k => 3.0, :Smax => 150.0, :v => 1e-3, :R => 5.0, :S₁res => 10.0]
# Create ODE problem
oprob = ODEProblem(irrigation_mod, u₀, tspan, params)
# osol = solve(oprob, Tsit5(), saveat=0.5)
# plot(osol)
# osol[:S₁]
# osol[:S₂]

condition = [tstep]
function affect!(integrator)
    integrator.p[4] += 5.0      # R is the 4th parameter !!!
end
cb = PresetTimeCallback(condition, affect!)
osol = solve(deepcopy(oprob), Tsit5(), saveat=10.0, callback=cb)
plot(osol)


osol.t
osol[:S₁]
osol[:S₂]

plot(osol.t, (osol[:S₁] .+ S₁res) ./ Smax .* R, xaxis="time [s]", label="runoff")
plot!(osol.t, (osol[:S₂].^2) .* v, label="outflow")
plot!(osol.t, k .* osol[:S₁] ./ Smax, label="percolation")


################################################################
## Calibration
################################################################

# # u₀ = [:S₁ => 0.0, :S₂ => 0.0]
# u₀ = [:S₁ => 140.0, :S₂ => 135.0]
# tspan = (0.0, 150.0)  # the time interval to solve on
# params = [:k => 4.2, :Smax => 143.0, :v => 1e-3, :R => 5.0, :S₁res => 10.0]
# oprob = ODEProblem(irrigation_mod, u₀, tspan, params)
# osol = solve(oprob, Tsit5(), saveat=10.0)
# plot(osol)

# S1sim = osol[:S₁]
# S2sim = osol[:S₂]

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
  params = [:k => thetas[1], :Smax => thetas[2], :v => 1e-3, :R => 5, :S₁res => 10]
  # Experiment 1
  u₀ = [0.0, 0.0]
  oprob = ODEProblem(irrigation_mod, u₀, tspan, params)
  osol = solve(oprob, Tsit5(), saveat=t_meas)
  S1_sol = osol[:S₁]
  S2_sol = osol[:S₂]
  J1 = (1/S1_sigma^2)*sum(abs2, S1_sol - S1_meas1) + (1/S2_sigma^2)*sum(abs2, S2_sol - S2_meas1)
  # Experiment 2
  u₀ = [140.0, 135.0]
  oprob = ODEProblem(irrigation_mod, u₀, tspan, params)
  osol = solve(oprob, Tsit5(), saveat=t_meas)
  S1_sol = osol[:S₁]
  S2_sol = osol[:S₂]
  J2 = (1/S1_sigma^2)*sum(abs2, S1_sol - S1_meas2) + (1/S2_sigma^2)*sum(abs2, S2_sol - S2_meas2)
  return J1+J2
end

thetas_init = [3.0, 150.0]
Jtheta_irrigation(thetas_init)

using Optim

results = optimize(Jtheta_irrigation, thetas_init, NelderMead())
thetas_opt = Optim.minimizer(results)
Optim.minimum(results)    # J-value / mse-value

params = [:k => thetas_opt[1], :Smax => thetas_opt[2], :v => 1e-3, :R => 5, :S₁res => 10]

# Experiment 1
u₀₁ = [0.0, 0.0]
oprob1 = ODEProblem(irrigation_mod, u₀₁, tspan, params)
osol1 = solve(oprob1, Tsit5(), saveat=1.0)

begin
  plot(osol1, labels=["S1 sim" "S2 sim"], xlabel="t")   # , xlims=(0, 80), ylims=(0, 10)
  scatter!(t_meas, S1_meas1, label="S1 meas", color=:blue, title="Experment 1")
  scatter!(t_meas, S2_meas1, label="S2 meas", color=:red, ylims=(0, 150))
end


# Experiment 2
u₀₂ = [140.0, 135.0]
oprob2 = ODEProblem(irrigation_mod, u₀₂, tspan, params)
osol2 = solve(oprob2, Tsit5(), saveat=1.0)

begin
  plot(osol2, labels=["S1 sim" "S2 sim"], xlabel="t")   # , xlims=(0, 80), ylims=(0, 10)
  scatter!(t_meas, S1_meas2, label="S1 meas", color=:blue, title="Experment 2")
  scatter!(t_meas, S2_meas2, label="S2 meas", color=:red, ylims=(0, 150))
end
