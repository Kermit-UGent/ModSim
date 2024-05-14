#======================================================================================
=======================================================================================
  Water tower
=======================================================================================
======================================================================================#

#=
https://docs.sciml.ai/Catalyst/stable/introduction_to_catalyst/catalyst_for_new_julia_users/
https://docs.sciml.ai/Catalyst/stable/catalyst_functionality/dsl_description/
https://docs.sciml.ai/DiffEqDocs/stable/getting_started/
=#

using Catalyst
using DifferentialEquations
using Plots



#======================================================================================
  Solving using Catalyst
======================================================================================#

water_tower = @reaction_network begin
  # reaction rate, substrate --> product
  1/τ, y --> 0           # y is degraded at rate b
  K/τ*Qin, 0 --> y       # y is created at rate K/τ*Qin
end

species(water_tower)
# y(t)
numspecies(water_tower)
# 1
parameters(water_tower)
# τ
# K
# Qin
numparams(water_tower)
# 3
reactions(water_tower)
# 1 / τ, y --> ∅
# (K*Qin) / τ, ∅ --> y
numreactions(water_tower)
# 2

osys  = convert(ODESystem, water_tower)   # convert rn to a ODESystem
equations(osys)
# Differential(t)(y(t)) ~ (-y(t)) / τ + (K*Qin) / τ
states(osys)
# y(t)
parameters(osys)
# τ
# K
# Qin

u0 = [:y => 0.0]
tspan = (0.0, 3600.0)  # the time interval to solve on
params = [:τ => 300.0, :K => 1.5, :Qin => 2.2]
oprob = ODEProblem(water_tower, u0, tspan, params)
osol = solve(oprob)
plot(osol)


condition = [1800.0]
affect!(integrator) = integrator.p[3] = 1.0      # Qin is the 3rd parameter !!!
ps_cb = PresetTimeCallback(condition, affect!)
osol = solve(deepcopy(oprob), Tsit5(), saveat=1.0, callback=ps_cb)
plot(osol)


function heaviside(t, t0, v0, v1)
    (v1 - v0) .* (t .> t0) .+ v0
end

Qin(t) = heaviside(t, 1800, 1.6, 2.0)

#======================================================================================
  Solving using Catalyst
======================================================================================#

rn = @reaction_network begin
    # reaction rate, substrate --> product
    1/τ, y --> 0           # y is degraded at rate b
    K/τ*Qin(t), 0 --> y    # y is created at rate K/τ*Qin(t)
end

# Graph(rn)

osys  = convert(ODESystem, rn)   # convert rn to a ODESystem
println(osys)
# Differential(t)(y(t)) ~ (-y(t)) / τ + (K*(1.6 + 0.3999999999999999(t > 1800))) / τ

u0 = [:y => 2.5]
tspan = (0.0, 3600.0)  # the time interval to solve on
params = [:τ => 300.0, :K => 1.5]
oprob = ODEProblem(rn, u0, tspan, params)
sol = solve(oprob, abstol=1e-9, reltol=1e-9)

plot(sol, linecolor=:red, linewidth=2, xaxis="time [s]", yaxis="", label="pressure [bar]")
plot!(sol.t, Qin(sol.t), linecolor=:blue, ls=:dash, linewidth=2)


# Influence of parameter τ

Qin(t) = heaviside(t, 0, 0.0, 1.0)
rn = @reaction_network begin
  1/τ, y --> 0
  K/τ*Qin(t), 0 --> y
end
u0 = [:y => 0.0]
attributes = [(300.0, :blue, :solid, "K=1.5, τ=300"),
             (225.0, :red, :dashdot, "K=1.5, τ=225"),
             (375.0, :black, :dash, "K=1.5, τ=375")]
P = plot(xaxis="time (s)", yaxis="step response")
for item in attributes
    τ=item[1]; lc=item[2]; ls=item[3]; label=item[4]
    params = [:τ => τ, :K => 1.5]
    oprob = ODEProblem(rn, u0, tspan, params)
    sol = solve(oprob, abstol=1e-9, reltol=1e-9)
    plot!(P, sol, lw=2, lc=lc, ls=ls, label=label)
end
plot!(P, sol.t, Qin(sol.t), lw=2, lc=:black, ls=:dot, label="step input")
display(P)

# Influence of parameter K

Qin(t) = heaviside(t, 0, 0.0, 1.0)
rn = @reaction_network begin
  1/τ, y --> 0
  K/τ*Qin(t), 0 --> y
end
u0 = [:y => 0.0]
attributes = [(1.0, :blue, :solid, "K=1.0, τ=300"),
             (1.5, :red, :dashdot, "K=1.5, τ=300"),
             (2.0, :black, :dash, "K=2.0, τ=300")]
P = plot(xaxis="time (s)", yaxis="step response")
for item in attributes
    K=item[1]; lc=item[2]; ls=item[3]; label=item[4]
    params = [:τ => 300, :K => K]
    oprob = ODEProblem(rn, u0, tspan, params)
    sol = solve(oprob, abstol=1e-9, reltol=1e-9)
    plot!(P, sol, lw=2, lc=lc, ls=ls, label=label)
end
plot!(P, sol.t, Qin(sol.t), lw=2, lc=:black, ls=:dot, label="step input")
display(P)


#======================================================================================
  Solving using DifferentialEquations
======================================================================================#

Qin(t) = heaviside(t, 1800, 1.6, 2.0)
function watertower!(du, u, p, t)
    τ, K = p
    du[1] = K/τ*Qin(t) - 1/τ*u[1]
end

u0 = [2.5]
params = (300.0, 1.5)
tspan = (0, 3600)
oprob = ODEProblem(watertower!, u0, tspan, params)
sol = solve(oprob, abstol=1e-9, reltol=1e-9)

plot(sol, linecolor=:red, linewidth=2, xaxis="time [s]", yaxis="", label="pressure [bar]")
plot!(sol.t, Qin(sol.t), linecolor=:blue, ls=:dash, linewidth=2)


# Influence of parameter τ

Qin(t) = heaviside(t, 0, 0.0, 1.0)
u0 = [0.0]
attributes = [(300.0, :blue, :solid, "K=1.5, τ=300"),
             (225.0, :red, :dashdot, "K=1.5, τ=225"),
             (375.0, :black, :dash, "K=1.5, τ=375")]
P = plot(xaxis="time (s)", yaxis="step response")
for item in attributes
    τ=item[1]; lc=item[2]; ls=item[3]; label=item[4]
    params = (τ, 1.5)
    oprob = ODEProblem(watertower!, u0, tspan, params)
    sol = solve(oprob, abstol=1e-9, reltol=1e-9)
    plot!(P, sol, lw=2, lc=lc, ls=ls, label=label)
end
plot!(P, sol.t, Qin(sol.t), lw=2, lc=:black, ls=:dot, label="step input")
display(P)

# Influence of parameter K

Qin(t) = heaviside(t, 0, 0.0, 1.0)
u0 = [0.0]
attributes = [(1.0, :blue, :solid, "K=1.0, τ=300"),
             (1.5, :red, :dashdot, "K=1.5, τ=300"),
             (2.0, :black, :dash, "K=2.0, τ=300")]
P = plot(xaxis="time (s)", yaxis="step response")
for item in attributes
    K=item[1]; lc=item[2]; ls=item[3]; label=item[4]
    params = (300.0, K)
    oprob = ODEProblem(watertower!, u0, tspan, params)
    sol = solve(oprob, abstol=1e-9, reltol=1e-9)
    plot!(P, sol, lw=2, lc=lc, ls=ls, label=label)
end
plot!(P, sol.t, Qin(sol.t), lw=2, lc=:black, ls=:dot, label="step input")
display(P)
