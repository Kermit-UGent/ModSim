# https://docs.sciml.ai/DiffEqDocs/dev/features/ensemble/


using Catalyst
using DifferentialEquations, Plots

infection_model = @reaction_network begin
	@parameters α=0.08 β=1.0e-6 r=0.2 m=0.4
	@species S(t)=9999000 I(t)=1000 D(t)=0 R(t)=0
	α * β, S + I --> 2I
	r * m, I --> D
	r * (1 - m), I --> R
end

species(infection_model)
# S(t)
# I(t)
# D(t)
# R(t)
parameters(infection_model)
# α
# β
# r
# m
reactions(infection_model)
# α*β, S + I --> 2*I
# m*r, I --> D
# (1 - m)*r, I --> R

osys  = convert(ODESystem, infection_model)
equations(osys)
unknowns(osys)
parameters(osys)

# u0 = [:S => 50, :I => 1, :D => 0, :R => 0]
# params = [:α => 0.15, :β => 0.1, :r => 0.2, :m => 0.6]
u0 = [:S => 9999000, :I => 1000, :D => 0, :R => 0]
tspan = (0.0, 90.0)
params = [:α => 0.08, :β => 1.0e-6, :r => 0.2, :m => 0.4]

oprob = ODEProblem(infection_model, u0, tspan, params)
# oprob = ODEProblem(infection_model, [], tspan, [])

osol = solve(oprob, Tsit5(), saveat=0.5)

plot(osol)
round((osol.u[end][3] / 1e7)*100, digits=1)    # 39.2 % of the population deceased
i = findfirst(==(2), osol.t)      # 5
osol.u[i]
# 4-element Vector{Float64}:
#     9.995907834374713e6
#  3318.927328621464
#   309.29531866588485
#   463.9429779988272


# Influence of r

# Influence of the duration of infection 1/r
# Check the effect of the average period in which infected persons are contagious if they
# are on average either 10, 5, 2 days or 1 day contagious. Try to interpret the results yourself.
# Ask yourself the following questions: (1) What are the trends in the results obtained?
# (2) How can this be explained from the model structure? (3) Is this consistent with what
# you intuitively have in mind?

params = [:α => 0.08, :β => 1.0e-6, :r => 0.1, :m => 0.4]
oprob = ODEProblem(infection_model, u0, tspan, params)
osol = solve(oprob, Tsit5(), saveat=0.5)
plot(osol)

params = [:α => 0.08, :β => 1.0e-6, :r => 0.2, :m => 0.4]
oprob = ODEProblem(infection_model, u0, tspan, params)
osol = solve(oprob, Tsit5(), saveat=0.5)
plot(osol)

params = [:α => 0.08, :β => 1.0e-6, :r => 0.5, :m => 0.4]
oprob = ODEProblem(infection_model, u0, tspan, params)
osol = solve(oprob, Tsit5(), saveat=0.5)
plot(osol)

params = [:α => 0.08, :β => 1.0e-6, :r => 1.0, :m => 0.4]
oprob = ODEProblem(infection_model, u0, tspan, params)
osol = solve(oprob, Tsit5(), saveat=0.5)
plot(osol)


# Suppose that regulations are such that on day 14 people need to reduce their contacts by 50%

params = [:α => 0.08, :β => 1.0e-6, :r => 0.2, :m => 0.4]
oprob = ODEProblem(infection_model, u0, tspan, params)
condition2 = [14.0]
function affect2!(integrator)
	integrator.ps[:β] = 0.5e-6     # β 
end
cb2 = PresetTimeCallback(condition2, affect2!)
osol2 = solve(deepcopy(oprob), Tsit5(), saveat=0.5, callback=cb2)
plot(osol2)
round((osol2.u[end][3] / 1e7)*100, digits=2)

# # On day 2 there are about 3318 people infected. Suppose that 3300 of them are put into isolation.

# params = [:α => 0.08, :β => 1.0e-6, :r => 0.2, :m => 0.4]
# oprob = ODEProblem(infection_model, u0, tspan, params)
# condition = [2.0]
# affect!(integrator) = integrator.u[2] -= 3300      # I is the 2-nd state variable
# # affect!(integrator) = integrator[:I] -= 3300
# ps_cb = PresetTimeCallback(condition, affect!)
# osol = solve(deepcopy(oprob), Tsit5(), saveat=0.5, callback=ps_cb)
# plot(osol)

# When there are 1e6 people infected, suppose that we can put 0.999e6 of them into isolation.
proceed_with_condition = [1]
function condition3(u, t, integrator)
	u[2] - 1.0e6*proceed_with_condition[1]
end
function affect3!(integrator, param=proceed_with_condition)
	integrator.u[2] -= 0.999e6
	param[1] = 0
end
cb3 = ContinuousCallback(condition3, affect3!)
osol3 = solve(deepcopy(oprob), Tsit5(), saveat=0.1, callback=cb3)
plot(osol3)
osol3.u[end]

osol3.u[end][3] + 0.999e6*0.4

###################################################################################
## EXERCISES
###################################################################################

# Exercise: Influence of α

# Evaluate the effect of a decreasing risk of infection after contact with an infected person,
# i.e. r = 0.2, β = 0.1 and α = 8%, 12%, 16% or 20%. Try to interpret the results yourself.
# Ask yourself the following questions: (1) What are the trends in the obtained results?
# (2) How can this be explained from the model structure? (3) Is this consistent with what
# you intuitively have in mind?

# params = [:α => 0.10, :β => 0.1, :r => 0.2, :m => 0.6]
params_ex1 = [:α => 0.08, :β => 1.0e-6, :r => 0.2, :m => 0.4]
oprob_ex1 = ODEProblem(infection_model, u0, tspan, params_ex1)
osol_ex1 = solve(oprob_ex1, Tsit5(), saveat=0.5)
plot(osol_ex1)

# params = [:α => 0.20, :β => 0.1, :r => 0.2, :m => 0.6]
params_ex1 = [:α => 0.12, :β => 1.0e-6, :r => 0.2, :m => 0.4]
oprob_ex1 = ODEProblem(infection_model, u0, tspan, params_ex1)
osol_ex1 = solve(oprob_ex1, Tsit5(), saveat=0.5)
plot(osol_ex1)

# params = [:α => 0.30, :β => 0.1, :r => 0.2, :m => 0.6]
params_ex1 = [:α => 0.16, :β => 1.0e-6, :r => 0.2, :m => 0.4]
oprob_ex1 = ODEProblem(infection_model, u0, tspan, params_ex1)
osol_ex1 = solve(oprob_ex1, Tsit5(), saveat=0.5)
plot(osol)

# params = [:α => 0.40, :β => 0.1, :r => 0.2, :m => 0.6]
params_ex1 = [:α => 0.20, :β => 1.0e-6, :r => 0.2, :m => 0.4]
oprob_ex1 = ODEProblem(infection_model, u0, tspan, params_ex1)
osol_ex1 = solve(oprob_ex1, Tsit5(), saveat=0.5)
plot(osol_ex1)


#----------------------------------------------------------------------------------
# Exercise: Administration of medicinal products

# Scientists have developed a medicine that heals sick people and makes them immune to
# the disease. After administering medication, the infection duration is reduced to two
# days. All treated patients heal and acquire immunity to the virus. The model will have
# to be extended with two additional parameters.
# - b: the fraction of infected persons undergoing treatment
# - rb : the rate at which the infected persons treated are no longer contagious (day)

# Administering the drug to a fraction of the infected affects three differential equations:
# those of I, D and R.
# - The fraction of infected persons treated (b) has a reduced infection duration.
# - The fraction of infected persons not receiving treatment (1 − b) still has the
#   same duration of infection.
# - Mortality m only affects the group of sick people who were not given any
#   medication.
# - All subjects treated acquire resistance
# - A fraction of the untreated also heals
# Check the effect on the epidemic when 0%, 25%, 50%, 75%, 100% of sick individuals are
# treated (with I0 = 1000, r = 0.2, α = 0.08, rb = 0.5). Interpret the figures obtained for
# this purpose:
# 1) Why does the peak in the number of infected residents shift to the right when the value
#    of b increases?
# 2) Why does the number of resistant residents first rise when the value of b increases and
#    then fall when the value of b continues to increase?

# Set-up a reaction model. Base yourself on infection_model
infection_med = @reaction_network begin
	@parameters α=0.08 β=1.0e-6 b=0.0 m=0.4 r=0.2 rb=0.5
	@species S(t)=9999000 I(t)=1000 D(t)=0 R(t)=0
	α * β, S + I --> 2I
	(1 - b) * m * r, I --> D
	((1 - b) * (1 - m) * r, b * rb), I --> R
    # b * rb, I --> R
end

parameters(infection_med)
# α
# β
# b
# m
# r
# rb

# Check out the differential equations
species(infection_med)
osys_ex2  = convert(ODESystem, infection_med)
equations(osys_ex2)
# 4-element Vector{Equation}:
#  Differential(t)(S(t)) ~ -S(t)*I(t)*α*β
#  Differential(t)(I(t)) ~ -b*rb*I(t) + (-1 + b)*(1 - m)*r*I(t) + (-1 + b)*m*r*I(t) + S(t)*I(t)*α*β
#  Differential(t)(D(t)) ~ (1 - b)*m*r*I(t)
#  Differential(t)(R(t)) ~ b*rb*I(t) + (1 - b)*(1 - m)*r*I(t)
unknowns(osys_ex2)
parameters(osys_ex2)

# u0 = [:S => 50, :I => 1, :D => 0, :R => 0]
# tspan = (0.0, 60.0)
# params = [:α => 0.15, :β => 0.1, :r => 0.2, :m => 0.6, :rb => 0.5, :b => 0.0]

u0 = [:S => 9999000, :I => 1000, :D => 0, :R => 0]
tspan = (0.0, 90.0)

params_ex2 = [:α => 0.08, :β => 1.0e-6, :b => 0.0, :m => 0.4, :r => 0.2, :rb => 0.5]

oprob_ex2 = ODEProblem(infection_med, u0, tspan, params_ex2)
osol_ex2 = solve(oprob_ex2, Tsit5(), saveat=0.5)
plot(osol_ex2)

params = [:α => 0.08, :β => 1.0e-6, :r => 0.2, :m => 0.4, :rb => 0.5, :b => 0.25]
oprob = ODEProblem(infection_med, u0, tspan, params)
osol = solve(oprob, Tsit5(), saveat=0.5)
plot(osol)

params = [:α => 0.08, :β => 1.0e-6, :r => 0.2, :m => 0.4, :rb => 0.5, :b => 0.50]
oprob = ODEProblem(infection_med, u0, tspan, params)
osol = solve(oprob, Tsit5(), saveat=0.5)
plot(osol)

params = [:α => 0.08, :β => 1.0e-6, :r => 0.2, :m => 0.4, :rb => 0.5, :b => 0.75]
oprob = ODEProblem(infection_med, u0, tspan, params)
osol = solve(oprob, Tsit5(), saveat=0.5)
plot(osol)

params = [:α => 0.08, :β => 1.0e-6, :r => 0.2, :m => 0.4, :rb => 0.5, :b => 1.00]
oprob = ODEProblem(infection_med, u0, tspan, params)
osol = solve(oprob, Tsit5(), saveat=0.5)
plot(osol)


#----------------------------------------------------------------------------------
# Exercise: Adding vaccination to the model

# Scientists have developed a vaccine that makes healthy people immediately immune to
# the disease.
# - Vaccination affects several differential equations:
# - Susceptible individuals are vaccinated at a rate of p vac (with unit 1/day).
# These persons can therefore no longer be infected.
# - The vaccinated persons become resistant.

infection_med_vac = @reaction_network begin
	@parameters α=0.08 β=1.0e-6 b=0.0 m=0.4 r=0.2 rb=0.5 pvac=0.00
	@species S(t)=9999000 I(t)=1000 D(t)=0 R(t)=0
	α * β, S + I --> 2I
	(1 - b) * m * r, I --> D
	(1 - b) * (1 - m) * r, I --> R
    b * rb, I --> R
    pvac, S --> R
end

parameters(infection_med_vac)
# α
# β
# b
# m
# r
# rb
# pvac

# Check out the differential equations
species(infection_med_vac)
osys  = convert(ODESystem, infection_med_vac)
equations(osys)
# 4-element Vector{Equation}:
#  Differential(t)(S(t)) ~ -pvac*S(t) - S(t)*I(t)*α*β
#  Differential(t)(I(t)) ~ -b*rb*I(t) + (-1 + b)*(1 - m)*r*I(t) + (-1 + b)*m*r*I(t) + S(t)*I(t)*α*β
#  Differential(t)(D(t)) ~ (1 - b)*m*r*I(t)
#  Differential(t)(R(t)) ~ pvac*S(t) + b*rb*I(t) + (1 - b)*(1 - m)*r*I(t)
unknowns(osys)
parameters(osys)

u0 = [:S => 9999000, :I => 1000, :D => 0, :R => 0]
tspan = (0.0, 90.0)

params = [:α => 0.08, :β => 1.0e-6, :r => 0.2, :m => 0.4, :rb => 0.5, :b => 0.25, :pvac => 0.00]

oprob = ODEProblem(infection_med_vac, u0, tspan, params)
osol = solve(oprob, Tsit5(), saveat=0.5)
plot(osol)
osol.u[end]

# https://docs.sciml.ai/Catalyst/stable/catalyst_applications/advanced_simulations/#advanced_simulations_callbacks

# Vaccination only from day 7 on.

condition4 = [7.0]
function affect4!(integrator)
	integrator.ps[:pvac] = 0.05     # pvac is the 7-th parameter !!!
end
ps_cb = PresetTimeCallback(condition4, affect4!)
osol = solve(deepcopy(oprob), Tsit5(), saveat=0.5, callback=ps_cb)
plot(osol)


