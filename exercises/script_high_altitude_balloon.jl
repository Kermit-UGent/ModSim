#======================================================================================
=======================================================================================
  High Altitude Balloon
=======================================================================================
======================================================================================#

using DifferentialEquations
using Plots

#======================================================================================
  Solving using DifferentialEquations
======================================================================================#

function bullet!(du, u, p, t)
    ρ₀, p₀, g, Cd, R, m = p
    y = u[1]
    v = u[2]
    ydot = v
    vdot = -g - (ρ₀*π*R^2*Cd)/(2*m)*v^2*exp(-(ρ₀*g*y)/p₀)
    du[1] = ydot
    du[2] = vdot
end

ρ₀  = 1.293      # [kg/m^3]
p₀ = 101300     # [Pa]
g = 9.81        # [m/s^2]
Cd = 0.2;       # drag coefficient
R = 10.6/2*1e-3 # radius of bullet
m = 0.026       # mass bullet [kg]

y₀ = 0.0
v₀ = 960.0

u0 = [y₀, v₀]
params = ρ₀, p₀, g, Cd, R, m
tspan = (0.0, 30.0)
oprob = ODEProblem(bullet!, u0, tspan, params)

sol = solve(oprob, abstol=1e-9, reltol=1e-9)

# convert 92-element Vector{Vector{Float64}} into 92×2 Matrix{Float64}
solmat = mapreduce(permutedims, vcat, sol.u)

t = sol.t
y = solmat[:,1]
v = solmat[:,2]

Py = plot(t, y, linecolor=:black, linewidth=1.5, xaxis="time [s]", yaxis="altitude [m]", label="y")
Pv = plot(t, v, linecolor=:blue, linewidth=1.5, xaxis="time [s]", yaxis="velocity [m/s]", label="v")
plot(Py, Pv, layout=(1, 2))



#======================================================================================
  Solving using ModelingToolkit
======================================================================================#

using ModelingToolkit

@parameters ρ₀, p₀, g, Cd, R, m
@variables t y(t) v(t)
D = Differential(t)

eqs = [D(y) ~ v
       D(v) ~ -g - (ρ₀*π*R^2*Cd)/(2*m)*v^2*exp(-(ρ₀*g*y)/p₀)]

@named osys = ODESystem(eqs, t, [y, v], [ρ₀, p₀, g, Cd, R, m])

u0 = [y => 0.0, v => 960.0]
params = [ρ₀ => 1.293, p₀ => 101300, g => 9.81, Cd => 0.2, R => 10.6/2*1e-3, m => 0.026]
tspan = (0.0, 30.0)

oprob = ODEProblem(osys, u0, tspan, params)
sol = solve(oprob, abstol=1e-9, reltol=1e-9)

solmat = mapreduce(permutedims, vcat, sol.u)

t = sol.t
y = solmat[:,1]
v = solmat[:,2]

Py = plot(t, y, linecolor=:black, linewidth=1.5, xaxis="time [s]", yaxis="altitude [m]", label="y")
Pv = plot(t, v, linecolor=:blue, linewidth=1.5, xaxis="time [s]", yaxis="velocity [m/s]", label="v")
plot(Py, Pv, layout=(1, 2))
