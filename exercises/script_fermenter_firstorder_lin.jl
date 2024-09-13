#======================================================================================
=======================================================================================
  Fermenter with first-order kinetics
=======================================================================================
======================================================================================#

#=
https://docs.sciml.ai/Catalyst/stable/introduction_to_catalyst/catalyst_for_new_julia_users/
https://docs.sciml.ai/Catalyst/stable/catalyst_functionality/dsl_description/
https://www.matecdev.com/posts/julia-plotting-multiple-plots.html
=#

using Catalyst
using DifferentialEquations
using Plots
using Symbolics


@variables V Q β Y S_in S X t
Dt = Differential(t)

f1 = Q/V*(S_in - S) - β/Y*S
f2 = -Q/V*X + β*S
g1 = S + X
g2 = S
g3 = X
g4 = S_in - S

A = Symbolics.jacobian([f1, f2], [S, X])
# 2×2 Matrix{Num}:
#  (-β) / Y + (-Q) / V         0
#                    β  (-Q) / V
B = Symbolics.jacobian([f1, f2], [S_in])
# 2×1 Matrix{Num}:
#  Q / V
#      0
C = Symbolics.jacobian([g1, g2, g3, g4], [S, X])
# 3×2 Matrix{Num}:
#  1  1
#  1  0
#  0  1
# -1  0
D = Symbolics.jacobian([g1, g2, g3, g4], [S_in])
# 4×1 Matrix{Num}:
#  0
#  0
#  0
#  1

[Differential(t)(S); Differential(t)(X)] ~ A*[S;X] + B*[S_in]


Differential(t)(S) ~ Q/V*(S_in - S) - β/Y*S

#======================================================================================
  Define the step function
======================================================================================#

# Step function, go from v0 to v1 at t = t0
function heaviside(t, t0, v0, v1)
    (v1 - v0) .* (t .> t0) .+ v0
end

#======================================================================================
  Calculating the operating point
======================================================================================#

# Defining symboloc variables
@variables V Q β Y S_in Sw Xw

# Setting up equations for finding the operating point (Sw, Xw)
eq1 = 0 ~ Q/V*(S_in - Sw) - β/Y*Sw
eq2 = 0 ~ -Q/V*Xw + β*Sw
# Solve for Sw and Xw
optn = Symbolics.solve_for([eq1, eq2],[Sw, Xw])
# There is one operating point here, define the values as initial conditions
S0 = substitute(optn[1], [V => 20, Q => 2, β => 0.03, Y => 0.67, S_in => 0.02])
X0 = substitute(optn[2], [V => 20, Q => 2, β => 0.03, Y => 0.67, S_in => 0.02])

#======================================================================================
  Solving using Catalyst
======================================================================================#

# Sin goes from 0.02 to 0.04 at t = 10
Sin(t) = heaviside(t, 10, 0.02, 0.04)

# Define all reactions
fermenter_firstorder = @reaction_network begin
    Q/V, S --> 0           # S is degraded at a rate Q/V*S
    Q/V*Sin(t), 0 --> S    # S is created at a rate Q/V*Sin(t)
    β/Y, S --> 0           # S is degraded at a rate β/Y*S
    β*S, 0 --> X           # X is created at a rate β*S
    Q/V, X --> 0
end

# Convert to an ODE system
osys  = convert(ODESystem, fermenter_firstorder)

# Check out the diff. eqns. equations
equations(osys)
# 2-element Vector{Equation}:
#  Differential(t)(S(t)) ~ (-S(t)*β) / Y + (Q*Sin) / V + (-Q*S(t)) / V
#  Differential(t)(X(t)) ~ (-Q*X(t)) / V + S(t)*β
# Check out the states
states(osys)
# 2-element Vector{SymbolicUtils.BasicSymbolic{Real}}:
#  S(t)
#  X(t)
# Check out the parameters
parameters(osys)
# 5-element Vector{SymbolicUtils.BasicSymbolic{Real}}:
#  Q
#  V
#  β
#  Y

# Assign values to parameters
V = 20; Q = 2; β = 0.03; Y = 0.67;
# Initial concentraation for Sin
Sin0 = Sin(0)
# Set the operating point (Sw, Xw) when Sin = Sin0 (0.02) to the
# initial condition (S0, X0)
S0 = (Q/V)/(Q/V + β/Y)*Sin0
X0 = β*V/Q*S0
# Set initial condition for ODE problem
u0 = [:S => S0, :X => X0]
# Set parameters for ODE problem
params = [:V => V, :Q => Q, :β => β, :Y => Y]
# Set time span for simulation
tspan = (0.0, 100.0)  # the time interval to solve on
# Create ODE problem
oprob = ODEProblem(fermenter_firstorder, u0, tspan, params)
# Solve ODE problem
sol = solve(oprob, abstol=1e-9, reltol=1e-9)
# sol consists of:
# sol.t   ->  vector of with time values
# sol.u   ->  vector of vector with [S, X] values
#                             [S1, X1]
#                             [S2, X2]
#                             ...
#                             [Sn, Xn]

# Convert sol.u to a matrix with a column for S and a column for X
solmat = mapreduce(permutedims, vcat, sol.u)
# Select the column for S, select the column for X
S = solmat[:,1]; X = solmat[:,2]
# Define t as the time vector
t = sol.t

P1 = plot(t, Sin(t), linecolor=:green, lw=2, ls=:dash, xaxis="time [h]", label="Sin")
plot!(P1, t, S+X, linecolor=:blue, lw=2, ls=:solid, label="S+X")
P2 = plot(t, Sin(t), linecolor=:green, lw=2, ls=:dash, xaxis="time [h]", label="Sin")
plot!(P2, t, S, linecolor=:red, lw=2, ls=:dashdot, label="S")
plot!(P2, t, X, linecolor=:blue, lw=2, ls=:solid, label="X")
P3 = plot(t, Sin(t), linecolor=:green, lw=2, ls=:dash, xaxis="time [h]", label="Sin")
plot!(P3, t, Sin(t)-S, linecolor=:blue, lw=2, ls=:solid, label="Sin-S")
plot(P1, P2, P3, layout=(3, 1))


