#=
CATALYST.JL CHEAT SHEET

AUTHOR: Michiel Stock (michiel.stock@ugent.be)
=#

using Catalyst
using Plots, DifferentialEquations

# defining a reaction system

mm = @reaction_network begin
    (kB, kD), S + E <--> ES  # reversible binding
    kP, ES --> P + E         # conversion of substrate by enzyme
end

species(mm)  # check the species

parameters(mm)  # check the parameters

equations(mm)  # check the equationd

@unpack S = mm  # extract parameter / variable

# convert ReactionSystem in an ODE system
osys = convert(ODESystem, mm)

# SIMULATION

# define paramters and intial values
u0map = [:S => 10.0, :E => 0.1, :P => 0.0, :ES => 0]
pmap = [:kB => 0.5, :kD => 0.1, :kP => 2.2]

# alternative, just as good!
u0map = [mm.S => 10.0, mm.E => 0.1, mm.P => 0.0, mm.ES => 0]
pmap = [mm.kB => 0.5, mm.kD => 0.1, mm.kP => 2.2]

tspan = (0.0, 100.)

oprob = ODEProblem(mm, u0map, tspan, pmap)

sol = solve(oprob, Tsit5())

plot(sol)  # plot all variables

plot(sol, idxs=[:S, :P])  # plot only S and P

# alternative, unpack variables
@unpack P, S = mm  # extract product and substrate
plot(sol, idxs=[P, S])

plot(sol, idxs=P/(S + P), title="Fraction of substrate converted")

# DEFAULT OPTIONS

mm = @reaction_network begin
    @species S(t)=10. E(t)=0.1 P(t)=0 ES(t)=0
    @parameters kB=0.5 kD=0.1 kP=2.2
    (kB, kD), S + E <--> ES
    kP, ES --> P + E
end

# no need to specify initial states and parameters
oprob = ODEProblem(mm, [], tspan)

plot(solve(oprob, Tsit5()))

# overwriting defauls
oprob = ODEProblem(mm, [:E=>0.2], tspan, [:kP=>1.5])

# adding annotation to the parameters
mm = @reaction_network begin
    @species S(t)=10. [description="substrate"] P(t)=0 [description="product"]
    @parameters kB=0.5 [description="binding rate"] kD=0.1 [description="dissociated rate rate"] kP=2.2 [description="conversion rate"]
    (kB, kD), S + E <--> ES
    kP, ES --> P + E
end

# OBESRVABLES

mm = @reaction_network begin
    @observables begin
        Etot ~ E + ES
    end
    (kB, kD), S + E <--> ES
    kP, ES --> P + E
end

observed(mm)

oprob = ODEProblem(mm, u0map, tspan, pmap)

sol = solve(oprob, Tsit5())

plot(sol, idxs=:Etot, ylims=(0,1))

# REMOVING CONSERVED QUANTITIES

equations(osys)

mm_conserved = convert(ODESystem, mm, remove_conserved = true)

# enzyme (E + ES) is conserved
equations(mm_conserved)

# still in observables
observed(mm_conserved)


# EVENTS

oprob = ODEProblem(mm, u0map, tspan, pmap)

sol = solve(oprob, Tsit5())

plot(sol)

## DISCRETE EVENTS (TIME-BAESD)

dilluting = [20.0] => [mm.E ~ mm.E/2, mm.S ~ mm.S/2,
                            mm.ES ~ mm.ES / 2, mm.P ~ mm.P/2]

@named mm_dillute = ReactionSystem(equations(mm), discrete_events=dilluting)

mm_dillute = complete(mm_dillute)

oprob = ODEProblem(mm_dillute, u0map, tspan, pmap)

sol = solve(oprob, Tsit5(), tstops=20)

plot(sol)

## CONTINUOUS EVENTS

substrate_feeding = [mm.S ~ 5] => [mm.S ~ mm.S + 5]

@named mm_fed = ReactionSystem(equations(mm), continuous_events=substrate_feeding)

mm_fed = complete(mm_fed)

oprob = ODEProblem(mm_fed, u0map, tspan, pmap)

sol = solve(oprob, Tsit5(), tstops=20)

plot(sol)

# STEADY STATE

mm_continuous = @reaction_network begin
    (kB, kD), S + E <--> ES  # reversible binding
    kP, ES --> P + E         # conversion of substrate by enzyme
    1, 0 --> S  # constant inflow
    1, (S, P) --> 0  # constant outflow
end

ssprob = SteadyStateProblem(mm_continuous, u0map, pmap)

sol = solve(ssprob)

# JUMP SIMULATION

# starting with 500 substrate molecules and 5 enzyme molecules
u0map_int = [:S => 500, :E => 5, :P => 0, :ES => 0]

jsys = JumpInputs(mm, u0map_int, tspan, pmap)
jprob = JumpProblem(jsys)

jsol = solve(jprob)
plot(jsol)


# STOCHASTIC DIFFERENTIAL EQUATIONS

sys = @reaction_network begin
    (k1, k2), A <--> B
end

u0map = [:A => 10., :B=> 200.]
pmap = [:k1 => 2, :k2=> 5]


sprob = SDEProblem(sys, u0map, (0., 20.), pmap)

sol = solve(sprob, STrapezoid(), dt=0.02)
plot(sol)

sys = @reaction_network begin
    @parameters η
    @default_noise_scaling η
    (k1, k2), A <--> B
end

pmap = [:k1 => 2, :k2=> 5, :η=>0.1]

sprob = SDEProblem(sys, u0map, (0., 20.), pmap)
sol = solve(sprob, STrapezoid(), dt=0.02)
plot(sol)

sys = @reaction_network begin
    @parameters η
    @default_noise_scaling η
    k1, A --> B, [noise_scaling = 0.0]
    k2, B --> A, [noise_scaling = η]
end

pmap = [:k1 => 2, :k2=> 5, :η=>1]
sprob = SDEProblem(sys, u0map, (0., 20.), pmap)
sol = solve(sprob, STrapezoid(), dt=0.02)
plot(sol)