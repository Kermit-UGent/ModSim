using Catalyst
using DifferentialEquations, Plots

#=
The garbage is damage casual for aging, which we will call X. To be concrete,
for mammals we assume that X is the total number of senescent cells in the body.
These senescent cells secrete factors that cause chronic inflammation
and reduced regeneration, leading to disease and decline. Later, when
we discuss other organisms, X will represent other forms of damage.

What is cellular senescence? Senescent cells are unique in that they
eventually stop multiplying but don't die off when they should. They
instead remain and continue to release chemicals that can trigger
inflammation.

Natural killer cells, also known as NK cells, are a type of cytotoxic
lymphocyte critical to the innate immune system.

De rommel is schade die voortkomt uit veroudering, die we X zullen noemen.
Om het concreet te maken, voor zoogdieren nemen we aan dat X het totale aantal
senescente cellen in het lichaam is. Deze senescente cellen scheiden
factoren af ​​die chronische ontstekingen veroorzaken en verminderde
regeneratie, wat leidt tot ziekte en achteruitgang. Later, wanneer
we andere organismen bespreken, zal X andere vormen van schade
vertegenwoordigen.

Het lichaam kan de meeste cellen die beschadigd raken repareren, maar
met een klein deel daarvan lukt dit niet en hun aantal neemt elk jaar
toe. De cel deelt dan niet verder, maar blijft wel in het lichaam zitten
en doet daar minder goed zijn werk. Dit is het proces dat we senescence
noemen.
=#

senescent_cells_rn = @reaction_network begin
    @species X(t)=0.0
    @parameters μ=0.00558 β=0.4464 κ=1.116 η=0.1
    @default_noise_scaling η
    μ*t, 0 --> X
    # β/(κ+X), X --> 0
    mm(X, β, κ), X => 0, [noise_scaling = 0.5]
end

osys  = convert(ODESystem, senescent_cells_rn)   # convert rn to a ODESystem
equations(osys)

u0 = [:X => 0.0]
tspan = (0.0, 120.0)
params = [:μ=>0.00558, :β=>0.4464, :κ=>1.116]

# 0.00186*3
# 0.1488*3

# 1) ODE problem

oprob = ODEProblem(senescent_cells_rn, u0, tspan, params)
osol = solve(oprob, Tsit5(), saveat=0.05)
plot(osol)


# 2) SDE problem

sprob = SDEProblem(senescent_cells_rn, u0, tspan, params)
ssol = solve(sprob, EM(), dt=0.1)
plot(ssol, ylim=(0, 6))

esprob = EnsembleProblem(sprob)
essol = solve(esprob, EM(), dt=0.1, save_everystep=true, trajectories=100)
plot(essol, ylim=(0, 6))
#=
Interpret the results. Ask yourself the following question:
1. Suppose that $5$ trillion senescent cells is about the
maximum a human body can bear. What is the (approximate)
corresponding range of ages?
=#