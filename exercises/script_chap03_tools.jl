
using Catalyst
using DifferentialEquations, Plots


####################################################################################
# Combustion model
####################################################################################

function radius!(du, u, p, t)
    du[1] = u[1]^2 - u[1]^3
    return du
end

u0 = [1/100]
u0 = [0.999]
tspan = (0.0, 200.0)
params = []
oprob_radius = ODEProblem(radius!, u0, tspan, params)
osol_radius = solve(oprob_radius, Tsit5())
plot(osol_radius)
osol_radius = solve(oprob_radius, Rosenbrock23())
plot(osol_radius)


####################################################################################
# ???
####################################################################################

some_process = @reaction_network begin
    0.04, A --> B
    3.0e7, 2B --> C + B
    10000.0, B + C --> A + C
end

species(some_process)

osys_some_process = convert(ODESystem, some_process)
equations(osys_some_process)
# Differential(t)(A(t)) ~ -0.04A(t) + 10000.0B(t)*C(t)
# Differential(t)(B(t)) ~ 0.04A(t) - 1.5e7(B(t)^2) - 10000.0B(t)*C(t)
# Differential(t)(C(t)) ~ 1.5e7(B(t)^2)

oprob_some_process = ODEProblem(some_process, [:A => 1.0, :B => 0.0, :C => 0.0], (0.0, 1000.0), [])
osol_some_process = solve(oprob_some_process, Tsit5())
plot(osol_some_process)


####################################################################################
# Example 15 (Bouncing ball)
####################################################################################

function ball!(du, u, g, t)
    y, v = u
    du[1] = v
    du[2] = -g
    return du
end

function condition(u, t, integrator)
    u[1]  # check when u[1] (i.e. x) == 0
end
function affect!(integrator)
    # nearly elastic collision
    integrator.u[2] = -0.9integrator.u[2]
end

# create the callback
cb = ContinuousCallback(condition, affect!)

u0 = [50.0, 0.0]
tspan = (0.0, 15.0)
g = 9.81
oprob_ball = ODEProblem(ball!, u0, tspan, g)
osol_ball = solve(deepcopy(oprob_ball), Tsit5(), callback=cb)
plot(osol_ball)


# Eénmaal botsen
stop_condition = [0]
function condition2(u, t, integrator)
    u[1] - 50*stop_condition[1]  # check when u[1] (i.e. x) == 0
end
function affect2!(integrator, param=stop_condition)
    # nearly elastic collision
    integrator.u[2] = -0.9integrator.u[2]
    param[1] = 1
end

# create the callback
cb2 = ContinuousCallback(condition2, affect2!)

u0 = [50.0, 0.0]
tspan = (0.0, 15.0)
g = 9.81
oprob_ball = ODEProblem(ball!, u0, tspan, g)
osol_ball = solve(deepcopy(oprob_ball), Tsit5(), callback=cb2)
plot(osol_ball)


####################################################################################
# Example 16 (Bioreactor dosing)
####################################################################################

bacterial_growth = @reaction_network begin
    @species X(t)=10 G(t)=8
    @parameters r=0.2 m=0.8
    r, X + G --> 2X
    m, X --> 0
end

species(bacterial_growth)
# X(t)
# G(t)

osys_bacterial_growth = convert(ODESystem, bacterial_growth)
equations(osys_bacterial_growth)
# Differential(t)(X(t)) ~ -m*X(t) + r*X(t)*G(t)
# Differential(t)(G(t)) ~ -r*X(t)*G(t)

dosetimes = 5:5:20
function affect!(integrator)
    integrator.u[2] += 10
end
cb = PresetTimeCallback(dosetimes, affect!)

oprob_bg = ODEProblem(bacterial_growth, [], (0, 20), [])
osol_bg = solve(oprob_bg, Tsit5(), callback=cb)
plot(osol_bg, lw=2, title="Dosed bioreactor")


####################################################################################
# Example 17 (Noisy competition)
####################################################################################

two_species_competition2 = @reaction_network begin
    @species A(t)=2.0 B(t)=2.5
    @parameters K=100 r=0.5 d=0.1 η1 η2 η3 η4
    r*(1-(A+B)/K), A --> 2A
    r*(1-(A+B)/K), B --> 2B
    d, (A, B) --> 0
end

reactions(two_species_competition2)
# r*(1 + (-B(t) - A(t)) / K), A --> 2*A
# r*(1 + (-B(t) - A(t)) / K), B --> 2*B
# d, A --> ∅
# d, B --> ∅

parameters(two_species_competition2)
# η1
# η2
# η3
# η4
# r
# K
# d

osys_two_species_competition2 = convert(ODESystem, two_species_competition2)
equations(osys_two_species_competition2)
# Differential(t)(A(t)) ~ -d*A(t) + r*A(t)*(1 + (-B(t) - A(t)) / K)
# Differential(t)(B(t)) ~ -d*B(t) + r*B(t)*(1 + (-B(t) - A(t)) / K)

tspan = (0.0, 50.0)
params = [:η1 => 1.0, :η2 => 1.0, :η3 => 1.0, :η4 => 1.0]
sprob_tsc2 = SDEProblem(two_species_competition2, [], tspan, params; noise_scaling = @parameters η1 η2 η3 η4)
ssol_tsc2 = solve(sprob_tsc2, STrapezoid())
ssol_tsc2 = solve(sprob_tsc2, LambaEM(), tstops=range(0.0, step=0.5, length=501))
ssol_tsc2 = solve(sprob_tsc2, LambaEM())
ssol_tsc2 = solve(sprob_tsc2, EM(), dt=0.1)
plot(ssol_tsc2)


####################################################################################
# Discrete SIR Model
####################################################################################

sir2 = @reaction_network begin
    @species S(t)=50 I(t)=5 R(t)=0
    @parameters η1 η2 η3
    beta, S + I --> 2I
    gamma, I --> R
end

osys_sir2 = convert(ODESystem, sir2)
equations(osys_sir2)
#  Differential(t)(S(t)) ~ -beta*S(t)*I(t)
#  Differential(t)(I(t)) ~ -gamma*I(t) + beta*S(t)*I(t)
#  Differential(t)(R(t)) ~ gamma*I(t)
# oprob_sir = ODEProblem(sir, [], (0, 120), [:beta => 0.005, :gamma => 0.05])
# oprob_sir = ODEProblem(sir, [], (0, 30), [:beta => 0.03, :gamma => 0.3])
# osol_sir = solve(oprob_sir, Rosenbrock23())
# plot(osol_sir)
# oprob_sir = ODEProblem(sir, [], (0, 30), [:beta => 0.01, :gamma => 0.3])
# osol_sir = solve(oprob_sir, Rosenbrock23())
# plot(osol_sir)
