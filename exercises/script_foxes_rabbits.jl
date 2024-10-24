
using Catalyst
using DifferentialEquations, Plots

#=
Rabbits live on some secluded territory. Their maximum growth rate coefficient
is gr [T^-1], and their population capacity is Rm [#rabbits]. The rabbits die
of age or sickness with a dying rate coefficient dr [T^-1].
At t=0, foxes intrude the territory and stay there. The foxes exclusively feed
themselves with the rabbits. They hunt the rabits at a rate proportional to the
number of foxes (proportionality factor is hr [T^-1 (#foxes)^-1]). The population 
of foxes grows at a rate proportional to the number of rabbits (proportionality
factor is fr [T^-1 (#rabbits)^-1]).

The foxes die of age or sickness with a dying rate coefficient dr [T^-1].

The initial number of rabbits on the territory is 89, the initial number of
foxes intruding the territory is 2.

1) Make simulations of the evolution of rabbits and foxes as a ODE
problem in the time interval [0, 10].

2) Make simulations of the evolution of rabbits and foxes as a discrete (jump) 
problem in the time interval [0, 10].

Assume the following parameter values:

gr=18.4 Rm=120 dr=2.0 hr=1.4
df=1.0 gf=0.05
=#

foxes_rabbits_rn = @reaction_network begin
    @species R(t)=89 F(t)=2
    @parameters gr=18.4 Rm=120 dr=2.0 hr=1.4 df=1.0 gf=0.05
    gr*(1 - R/Rm), R --> 2R          # natural population growth of the rabbits
    dr, R --> 0                      # deaths by age or sickness of the rabbits
    hr*F, R --> 0                    # hunting of rabbits by the foxes
    gf*R, F --> 2F                   # gainings of the foxes by hunting rabbits
    df, F --> 0                      # deaths by age or sickness of the foxes
end

# species(foxes_rabbits_rn)
# parameters(foxes_rabbits_rn)
# reactions(foxes_rabbits_rn)

osys  = convert(ODESystem, foxes_rabbits_rn)   # convert rn to a ODESystem
equations(osys)
# Differential(t)(R(t)) ~ -dr*R(t) + gr*(1 + (-R(t)) / Rm)*R(t) - hr*R(t)*F(t)
# Differential(t)(F(t)) ~ -df*F(t) + gf*R(t)*F(t)

u0 = [:R => 89, :F => 2]
tspan = (0.0, 5.0)
params = [:gr=>18.4, :Rm=>120, :dr=>2.0, :hr=>1.4, :gf=>0.05, :df=>1.0]

# 1) ODE problem

oprob = ODEProblem(foxes_rabbits_rn, [], tspan, [])
osol = solve(oprob, Tsit5(), saveat=0.05)
plot(osol)

# 2) Discrete (jump) problem

dprob = DiscreteProblem(foxes_rabbits_rn, u0, tspan, params)
jdprob = JumpProblem(foxes_rabbits_rn, dprob, Direct())

jdsol = solve(jdprob, SSAStepper())
plot(jdsol)

