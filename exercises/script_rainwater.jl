#======================================================================================
=======================================================================================
  Tractor seat
=======================================================================================
======================================================================================#

using DifferentialEquations
using Plots

function heaviside(t, t0, v0, v1)
    (v1 - v0) .* (t .> t0) .+ v0
end

#======================================================================================
  Solving using DifferentialEquations
======================================================================================#

#======================================================================================
  Simulation - Filling
======================================================================================#

function rainwater_fill!(du, u, p, t)
    A, c = p
    du[1] = Q(t) / A
end

t1 = 8; t2 = 18; Qin = 1.4
Q(t) = heaviside(t, t1, 0, Qin) + heaviside(t, t2, 0, -1.4)
c = 2
A = 5
h0 = 0

u0 = [h0]
params = A, c
tspan = (0.0, 24.0)
oprob = ODEProblem(rainwater_fill!, u0, tspan, params)

sol = solve(oprob, abstol=1e-9, reltol=1e-9)

solmat = mapreduce(permutedims, vcat, sol.u)

t = sol.t
h = solmat[:,1]

P = plot(t, h, linecolor=:blue, linewidth=2, xaxis="time [s]", yaxis="height [m]", label="h")
plot!(P, t, Q(t), linecolor=:black, lw=2, ls=:dot, label="Q_in")


#======================================================================================
  Simulation - Draining
======================================================================================#

function rainwater_drain!(du, u, p, t)
    A, c = p
    du[1] = -c / A * u[1]
end

c = 2
A = 5
h0 = 3

u0 = [h0]
params = A, c
tspan = (0.0, 10.0)
oprob = ODEProblem(rainwater_drain!, u0, tspan, params)

sol = solve(oprob, abstol=1e-9, reltol=1e-9)

solmat = mapreduce(permutedims, vcat, sol.u)

t = sol.t
h = solmat[:,1]

plot(t, h, linecolor=:blue, linewidth=2, xaxis="time [s]", yaxis="height [m]", label="h")
