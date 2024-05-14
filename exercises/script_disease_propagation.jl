#======================================================================================
=======================================================================================
  Disease propagation
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


#======================================================================================
  Solving using Catalyst
======================================================================================#

# Ontmoeting van S en I (S + I) levert zal S omzetten in I, dus
# dit levert 2I in het totaal (S + I --> 2I) aan een tempo α * β
# I wordt omgezet in D aan een tempo m*r
# I wordt omgezet in R aan een tempo (1-m)*r
infection = @reaction_network begin
	α * β, S + I --> 2I
	r * m, I --> D
	r * (1 - m), I --> R
end

# Alternatief:
# S wordt omgezet in I aan een tempo α*β*I
# I wordt omgezet in D aan een tempo m*r
# I wordt omgezet in R aan een tempo (1-m)*r
# infection = @reaction_network begin
#   α*β*I, S --> I
#   m*r, I --> D
#   (1-m)*r, I --> R
# end

# Alternative:
# infection = @reaction_network begin
#     α*β*I, S --> I
#     (m*r, (1-m)*r), I --> (D, R)
# end

# Convert to an ODE system
osys  = convert(ODESystem, infection)
# println(osys)
# Check the differential equations
equations(osys)
# 4-element Vector{Equation}:
#  Differential(t)(S(t)) ~ -S(t)*I(t)*α*β
#  Differential(t)(I(t)) ~ (-1 + m)*r*I(t) - m*r*I(t) + S(t)*I(t)*α*β
#  Differential(t)(D(t)) ~ m*r*I(t)
#  Differential(t)(R(t)) ~ (1 - m)*r*I(t)
# Check the states
states(osys)
# 4-element Vector{SymbolicUtils.BasicSymbolic{Real}}:
#  S(t)
#  I(t)
#  D(t)
#  R(t)
# Check the parameters
parameters(osys)
# 4-element Vector{SymbolicUtils.BasicSymbolic{Real}}:
#  α
#  β
#  r
#  m
# Check the independend variable
osys.iv
# t

# Assigne values to the parameters
α = 0.05; β = 1e-6; r = 0.1; m = 0.6;
# Total number of individuals
N = 1e7;
# Assign values to the initial conditions
I0 = 1000; S0 = N - I0; D0 = 0; R0 = 0;
# Set initial condition for the ODE problem
u0 = [:S => S0, :I => I0, :D => D0, :R => R0]
# 4-element Vector{Pair{Symbol, Float64}}:
#  :S => 9.999e6
#  :I => 1000.0
#  :D => 0.0
#  :R => 0.0
# Set parameters for the ODE problem
params = [:α => α, :β => β, :r => r, :m => m]
# 4-element Vector{Pair{Symbol, Float64}}:
#  :α => 0.05
#  :β => 1.0e-6
#  :r => 0.5
#  :m => 0.6
# Set time span for simulation
tspan = (0.0, 150.0)  # the time interval to solve on
# Create ODE problem
oprob = ODEProblem(infection, u0, tspan, params)
# Solve ODE problem
sol = solve(oprob, abstol=1e-9, reltol=1e-9)
# sol consists of:
# sol.t   ->  vector of with time values
# sol.u   ->  vector of vector with [S, I, R, D] values for eacht time instant
#                             [S1, I1, D1, R1]
#                             [S2, I2, D2, R2]
#                             ...
#                             [Sn, In, Dn, Rn]

plot(sol)


# Convert sol.u to a matrix with columns for S, I, D and R
solmat = mapreduce(permutedims, vcat, sol.u)
# Select the column for S, I, D and R
S = solmat[:,1]; I = solmat[:,2]; D = solmat[:,3]; R = solmat[:,4];
# Define t as the time vector
t = sol.t

plot(t, S, linecolor=:blue, lw=2, ls=:solid, xaxis="time [s]", label="S")
plot!(t, I, linecolor=:red, lw=2, ls=:dashdot, xaxis="time [s]", label="I")
plot!(t, D, linecolor=:black, lw=2, ls=:dash, xaxis="time [s]", label="D")
plot!(t, R, linecolor=:lightblue, lw=2, ls=:dot, xaxis="time [s]", label="R")


# Influence of parameter r

u0 = [:S => S0, :I => I0, :D => D0, :R => R0]  # initial values, unmodified    
# Each item in attibutes contains:
#   (new r-value, line color, line style, label)
attributes = [(0.1, :blue, :solid, "r=0.1"),
             (0.2, :red, :dash, "r=0.2"),
             (0.5, :green, :dashdot, "r=0.5")]
PS = plot(xaxis="time (s)", yaxis="S")    # Create plot for S-values
PI = plot(xaxis="time (s)", yaxis="I")    # Create plot for I-values
PD = plot(xaxis="time (s)", yaxis="D")    # Create plot for D-values
PR = plot(xaxis="time (s)", yaxis="R")    # Create plot for R-values
# Iterate over attributes in order to use the new r-values
for item in attributes
    r=item[1]; lc=item[2]; ls=item[3]; label=item[4]     # get r, lc, ls and label
    params = [:α => α, :β => β, :r => r, :m => m]        # put new r in there
    oprob = ODEProblem(infection, u0, tspan, params)     # create ODE problem
    sol = solve(oprob, abstol=1e-9, reltol=1e-9)         # solve problem
    solmat = mapreduce(permutedims, vcat, sol.u)         # convert to matrix
    S = solmat[:,1]; I = solmat[:,2]; D = solmat[:,3]; R = solmat[:,4]; t = sol.t;
    plot!(PS, t, S, lw=2, lc=lc, ls=ls, label=label)     # add S vector to plot
    plot!(PI, t, I, lw=2, lc=lc, ls=ls, label=label)     # add I vector to plot
    plot!(PD, t, D, lw=2, lc=lc, ls=ls, label=label)     # add D vector to plot
    plot!(PR, t, R, lw=2, lc=lc, ls=ls, label=label)     # add R vector to plot
end
plot(PS, PI, PD, PR, layout=(2,2))                       # put plots in a 2x2 layout
