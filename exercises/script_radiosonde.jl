#======================================================================================
=======================================================================================
  Radiosonde Model
=======================================================================================
======================================================================================#

using DifferentialEquations
using Plots

#======================================================================================
  Solving using DifferentialEquations
======================================================================================#

#======================================================================================
  Model functions
======================================================================================#

function radiosonde_rise!(du, u, p, t)
    ρ₀, p₀, R₀, Rprime, m, g, Cw, Cwprime = p

    y = u[1]; v = u[2]

    ydot = v
    vdot = (4 * ρ₀ * π * R₀^3 * g) / (3 * m) * exp(ρ₀ * g * y / 2 / p₀) - g - (ρ₀ * π * R₀ * Cw) / (2 * m) * v^2

    du[1] = ydot; du[2] = vdot
end


function radiosonde_fall!(du, u, p, t)
    ρ₀, p₀, R₀, Rprime, m, g, Cw, Cwprime = p

    y = u[1]; v = u[2]

    ydot = v
    vdot = - g - sign(v) * (ρ₀ * π * Rprime^2 * Cwprime) / (2 * m) * v^2 * exp(-ρ₀ * g * y / p₀)

    du[1] = ydot; du[2] = vdot
end


#======================================================================================
  Inputs, Outputs, States, Parameters, Constants
======================================================================================#

# Inputs: None
# Outputs: y and v
# States: y and v
# Parameters: R₀, Rprime, m, Cw, Cwprime
# Constants: ρ₀, p₀, g

#======================================================================================
  Parameter and constant values
======================================================================================#

ρ₀ = 1.293     # air density at sea level [kg/m3]
p₀ = 101300    # pressure at sea level [Pa]
R₀ = 1.0       # radius of balloon at sea level [m]
Rprime = 0.75  # radius of parachute [m]
m = 4.0        # mass of radiosonde [kg]
g = 9.81       # gravitational acceleration [m/s2]
Cw = 0.5       # drag/friction coefficient for balloon [-]
Cwprime = 1.0  # drag/friction coefficient for parachute [-]

#======================================================================================
  Part I, with balloon - Simulation in the range [0 s, 10 s]
======================================================================================#

# initial values
y₀ = 0.0; v₀ = 0.0

u0 = [y₀, v₀]
params = ρ₀, p₀, R₀, Rprime, m, g, Cw, Cwprime
tspan = (0.0, 10.0)
oprob = ODEProblem(radiosonde_rise!, u0, tspan, params)

sol = solve(oprob, abstol=1e-9, reltol=1e-9, saveat=0.01)

solmat = mapreduce(permutedims, vcat, sol.u)

t = sol.t
y = solmat[:,1]
v = solmat[:,2]

Py = plot(t, y, linecolor=:black, linewidth=1.5, xaxis="time [s]", yaxis="altitude [m]", label="y")
Pv = plot(t, v, linecolor=:blue, linewidth=1.5, xaxis="time [s]", yaxis="velocity [m/s]", label="v")
plot(Py, Pv, layout=(1, 2))


#======================================================================================
  Part I, with balloon - Simulation in the range [0 s, 3600 s]
======================================================================================#

# initial values
y₀ = 0.0; v₀ = 0.0

u0 = [y₀, v₀]
params = ρ₀, p₀, R₀, Rprime, m, g, Cw, Cwprime
tspan = (0.0, 3600.0)
oprob = ODEProblem(radiosonde_rise!, u0, tspan, params)

sol = solve(oprob, abstol=1e-9, reltol=1e-9, saveat=1)

solmat = mapreduce(permutedims, vcat, sol.u)

t = sol.t
y = solmat[:,1]
v = solmat[:,2]

Py = plot(t, y, linecolor=:black, linewidth=1.5, xaxis="time [s]", yaxis="altitude [m]", label="y")
Pv = plot(t, v, linecolor=:blue, linewidth=1.5, xaxis="time [s]", yaxis="velocity [m/s]", label="v")
plot(Py, Pv, layout=(1, 2))


# What is the time needed for the sonde for reaching an altitude of 30 km,
# and what is its velocity at that point?

i = findfirst(>(30000), y)    # 3556
t[i]      # 3555.0 s
y[i]      # 30012.563523298468 m
v[i]      # 17.431599242694382 m/s


#======================================================================================
  Part II, with parachute - Simulation in the range [0 s, 20 s]
======================================================================================#

# initial values
y₀ = 30000.0; v₀ = 20.0

u0 = [y₀, v₀]
params = ρ₀, p₀, R₀, Rprime, m, g, Cw, Cwprime
tspan = (0.0, 20.0)
oprob = ODEProblem(radiosonde_fall!, u0, tspan, params)

sol = solve(oprob, abstol=1e-9, reltol=1e-9, saveat=0.01)

solmat = mapreduce(permutedims, vcat, sol.u)

t = sol.t
y = solmat[:,1]
v = solmat[:,2]

Py = plot(t, y, linecolor=:black, linewidth=1.5, xaxis="time [s]", yaxis="altitude [m]", label="y")
Pv = plot(t, v, linecolor=:blue, linewidth=1.5, xaxis="time [s]", yaxis="velocity [m/s]", label="v")
plot(Py, Pv, layout=(1, 2))


#======================================================================================
  Part II, with parachute - Simulation in the range [0 s, 3600 s]
======================================================================================#

# initial values
y₀ = 30000.0; v₀ = 20.0

u0 = [y₀, v₀]
params = ρ₀, p₀, R₀, Rprime, m, g, Cw, Cwprime
tspan = (0.0, 3600.0)
oprob = ODEProblem(radiosonde_fall!, u0, tspan, params)

sol = solve(oprob, abstol=1e-9, reltol=1e-9, saveat=1)

solmat = mapreduce(permutedims, vcat, sol.u)

t = sol.t
y = solmat[:,1]
v = solmat[:,2]

Py = plot(t, y, linecolor=:black, linewidth=1.5, xaxis="time [s]", yaxis="altitude [m]", label="y")
Pv = plot(t, v, linecolor=:blue, linewidth=1.5, xaxis="time [s]", yaxis="velocity [m/s]", label="v")
plot(Py, Pv, layout=(1, 2))


i = findfirst(<(0.0), y)    # 2314
t[i]      # 2313.0 s
y[i]      # -4.774964905901471 m
v[i]      # -5.859516605186464 m/s
