####################################################################################

using Catalyst

rn = @reaction_network begin
    k1, A + B --> C
    k2, C --> B
    k3, A + A --> B + C
end

species(rn)
# A(t)
# B(t)
# C(t)

parameters(rn)
# k1
# k2
# k3

reactions(rn)
# k1, A + B --> C
# k2, C --> B
# k3, 2*A --> B + C

osys = convert(ODESystem, rn)

equations(osys)
# Differential(t)(A(t)) ~ -k1*B(t)*A(t) - k3*(A(t)^2)
# Differential(t)(B(t)) ~ k2*C(t) - k1*B(t)*A(t) + (1//2)*k3*(A(t)^2)
# Differential(t)(C(t)) ~ -k2*C(t) + k1*B(t)*A(t) + (1//2)*k3*(A(t)^2)


####################################################################################
# Example 4 (Hydrogen and nitric oxide reaction)
####################################################################################

using DifferentialEquations, Plots

reactionsys = @reaction_network begin
    (r1, r1), 2NO <--> N2O2
    r2, N2O2 + H2 --> N2O + H2O
    r3, N2O + H2 --> N2 + H2O
end

osys_example_04 = convert(ODESystem, reactionsys)

equations(osys_example_04)
# Differential(t)(NO(t)) ~ 2r1*N2O2(t) - r1*(NO(t)^2)
# Differential(t)(N2O2(t)) ~ -r1*N2O2(t) + (1//2)*r1*(NO(t)^2) - r2*H2(t)*N2O2(t)
# Differential(t)(H2(t)) ~ -r2*H2(t)*N2O2(t) - r3*N2O(t)*H2(t)
# Differential(t)(N2O(t)) ~ r2*H2(t)*N2O2(t) - r3*N2O(t)*H2(t)
# Differential(t)(H2O(t)) ~ r2*H2(t)*N2O2(t) + r3*N2O(t)*H2(t)
# Differential(t)(N2(t)) ~ r3*N2O(t)*H2(t)

# define an ODE problem
reactsyst_prob = ODEProblem(reactionsys,
        [:NO=>5.2, :H2O=>0, :H2=>5.1, :N2O2=>0, :N2O=>0, :N2=>0],  # initial values
        (0.0, 60.0),                                               # time interval
        [:r1=>1e2, :r2=>0.1, :r3=>50])                             # par. values

# solve the system
sol04 = solve(reactsyst_prob, Rosenbrock23())

# plot the solution
plot(sol04.t, sol04[:NO], color=:blue, linestyle=:dash, xlims=(0, 60), label="NO", xlabel="t", ylabel="Concentration [mol/L]")
plot!(sol04.t, sol04[:N2O2], color=:red, linestyle=:dashdot, label="N2O2")
plot!(sol04.t, sol04[:H2], color=:green, linestyle=:dashdotdot, label="H2")
plot!(sol04.t, sol04[:N2O], color=:darkblue, linestyle=:dot, label="N2O")

plot(sol04)


####################################################################################
# Example 5 (Hydrogen and nitric oxide reaction (bis))
####################################################################################

reactionsys2 = @reaction_network begin
    r * NO^2 * H2, 2H2 + 2NO => 0
end

osys_example_05 = convert(ODESystem, reactionsys2)

equations(osys_example_05)
# Differential(t)(H2(t)) ~ -2r*(NO(t)^2)*H2(t)
# Differential(t)(NO(t)) ~ -2r*(NO(t)^2)*H2(t)

# define an ODE problem
reactsyst2_prob = ODEProblem(reactionsys2,
        [:NO=>5.2, :H2=>5.1],  # initial values
        (0.0, 60.0),           # time interval
        [:r=>0.1])                             # par. values

# solve the system
sol05 = solve(reactsyst2_prob, Rosenbrock23())

plot(sol05)


####################################################################################
# Example 6 (Water tank (bis))
####################################################################################

tank = @reaction_network begin
    @species V(t)=0
    @parameters q=10 A=0.5^2*pi*200 r=10
    q, 0 --> V       # incoming water
    r / A, V --> 0   # emptying
end

osys_tank = convert(ODESystem, tank)

equations(osys_tank)
# Differential(t)(V(t)) ~ q + (-r*V(t)) / A

oprob_tank = ODEProblem(tank, [], (0.0, 50.0), [])

osol_tank = solve(oprob_tank, Rosenbrock23())
plot(osol_tank)


####################################################################################
# Example 7 (Bacterial growth)
####################################################################################

# A bacteria B (expressed in colony forming units, CFU) that grows with a doubling time of one hour.
growth1 = @reaction_network begin
    @species B(t)=1
    @parameters r=log(2)
    r, B --> 2B
end

osys_growth1 = convert(ODESystem, growth1)
equations(osys_growth1)
# Differential(t)(B(t)) ~ r*B(t)

oprob_growth1 = ODEProblem(growth1, [], (0, 20), [])
osol_growth1 = solve(oprob_growth1, Rosenbrock23())
plot(osol_growth1, title="Exponential growth")


# Assume that the bacteria growth follows logistic growth
growth2 = @reaction_network begin
    @species B(t)=1
    @parameters r=log(2) K=1e3
    r*(1-B/K), B --> 2B
end

reactions(growth2)
#  r*(1 + (-B(t)) / K), B --> 2*B

osys_growth2 = convert(ODESystem, growth2)
equations(osys_growth2)
# Differential(t)(B(t)) ~ r*B(t)*(1 + (-B(t)) / K)

oprob_growth2 = ODEProblem(growth2, [], (0, 20), [])
osol_growth2 = solve(oprob_growth2, Rosenbrock23())
plot(osol_growth2, title="Logistic growth")


# Equivalent with growth2
growth3 = @reaction_network begin
    @species B(t)=1
    @parameters r=log(2) K=1e3
    r, B --> 2B
    r/K, B + B --> 0
end

reactions(growth3)
# r, B --> 2*B
# r / K, 2*B --> ∅

osys_growth3 = convert(ODESystem, growth3)
equations(osys_growth3)
# Differential(t)(B(t)) ~ (-r*(B(t)^2)) / K + r*B(t)

oprob_growth3 = ODEProblem(growth3, [], (0, 20), [])
osol_growth3 = solve(oprob_growth3, Rosenbrock23())
plot(osol_growth3, title="Logistic growth")


####################################################################################
# Example 8 (Lotka-Volterra)
####################################################################################

lotka_volterra = @reaction_network begin
    alpha, x --> 2x        # reproduction prey
    beta, y --> 0          # mortality pred
    gamma, x + y --> 1.1y  # predation
end

reactions(lotka_volterra)

osys_lv = convert(ODESystem, lotka_volterra, combinatoric_ratelaws=false)
equations(osys_lv)
# Differential(t)(x(t)) ~ alpha*x(t) - gamma*x(t)*y(t)
# Differential(t)(y(t)) ~ -beta*y(t) + 0.10000000000000009gamma*x(t)*y(t)

oprob_lv = ODEProblem(lotka_volterra, [:x => 1.0, :y => 1.0], (0.0, 50.0), [:alpha => 1.5, :beta => 0.5, :gamma => 0.5], combinatoric_ratelaws=false)
osol_lv = solve(oprob_lv, Rosenbrock23())
plot(osol_lv)


####################################################################################
# Zeroth-order processes
####################################################################################

# In Catalyst, production is zeroth-order by default.
zop = @reaction_network begin
    r, 0 --> y
end

osys_zop = convert(ODESystem, zop)
equations(osys_zop)
# Differential(t)(y(t)) ~ r


# Degradation is assumed to be first-order, so to specific a process where y
# is removed at a constant rate, one uses r, y => 0
zod = @reaction_network begin
    r, y => 0
end

osys_zod = convert(ODESystem, zod)
equations(osys_zod)
# Differential(t)(y(t)) ~ -r

oprob_zod = ODEProblem(zod, [:y => 10.0], (0, 20), [:r => 1])
osol_zod = solve(oprob_zod, Rosenbrock23())
plot(osol_zod)


# One has to be cautious when having zeroth-order removal of a
# species, as these states usually cannot attain negative values.
zod2 = @reaction_network begin
    ifelse(y>0, r, 0), y => 0
end

osys_zod2 = convert(ODESystem, zod2)
equations(osys_zod2)
# Differential(t)(y(t)) ~ -ifelse(y(t) > 0, r, 0)

oprob_zod2 = ODEProblem(zod2, [:y => 10.0], (0, 20), [:r => 1])
osol_zod2 = solve(oprob_zod2, Rosenbrock23())
plot(osol_zod2, xlims=(0, 20))


####################################################################################
# First-order processes
####################################################################################

# First order production
fop = @reaction_network begin
    r, y --> 2y
    # r, y => 2y
end
osys_fop = convert(ODESystem, fop)
equations(osys_fop)
# Differential(t)(y(t)) ~ r*y(t)
oprob_fop = ODEProblem(fop, [:y => 1.0], (0, 20), [:r => 0.1])
osol_fop = solve(oprob_fop, Rosenbrock23())
plot(osol_fop, xlims=(0, 20), ylims=(0, 10), title="First order production")


# First order degradation
fod = @reaction_network begin
    r, y --> 0
end
osys_fod = convert(ODESystem, fod)
equations(osys_fod)
#  Differential(t)(y(t)) ~ -r*y(t)
oprob_fod = ODEProblem(fod, [:y => 10.0], (0, 20), [:r => 0.1])
osol_fod = solve(oprob_fod, Rosenbrock23())
plot(osol_fod, xlims=(0, 20), ylims=(0, 10), title="First order degradation")


# First order kinetics
fok = @reaction_network begin
    r, A --> B
end
osys_fok = convert(ODESystem, fok)
equations(osys_fok)
# Differential(t)(y(t)) ~ r*y(t)
oprob_fok = ODEProblem(fok, [:A => 10.0, :B => 1.0], (0, 20), [:r => 0.1])
osol_fok = solve(oprob_fok, Rosenbrock23())
plot(osol_fok, xlims=(0, 20), ylims=(0, 10), title="First order kinetics")


####################################################################################
# Second-order processes
####################################################################################

# Second order production
sop = @reaction_network begin
    r, y + y --> 4y
    # r, 2y --> 4y      # equivalent
end
osys_sop = convert(ODESystem, sop)
equations(osys_sop)
#  Differential(t)(y(t)) ~ r*(y(t)^2)
oprob_sop = ODEProblem(sop, [:y => 1.0], (0, 20), [:r => 0.045])
osol_sop = solve(oprob_sop, Rosenbrock23())
plot(osol_sop, xlims=(0, 20), ylims=(0, 10), title="Second order production")
# Asymptote: (r*y0)^(-1)
(0.045*1.0)^(-1)


# Second order production
sop_abc = @reaction_network begin
    r, A + B --> C
end
osys_sop_abc = convert(ODESystem, sop_abc)
equations(osys_sop_abc)
# Differential(t)(A(t)) ~ -r*B(t)*A(t)
# Differential(t)(B(t)) ~ -r*B(t)*A(t)
# Differential(t)(C(t)) ~ r*B(t)*A(t)
oprob_sop_abc = ODEProblem(sop_abc, [:A => 10.0, :B => 8.0, :C => 0.5], (0, 20), [:r => 0.05])
osol_sop_abc = solve(oprob_sop_abc, Rosenbrock23())
plot(osol_sop_abc, xlims=(0, 20), ylims=(0, 10), title="Second order production A + B --> C")


# Second order production
sop_aac = @reaction_network begin
    r, A + A --> C
end
osys_sop_aac = convert(ODESystem, sop_aac)
equations(osys_sop_aac)
# Differential(t)(A(t)) ~ -r*(A(t)^2)
# Differential(t)(C(t)) ~ (1//2)*r*(A(t)^2)
oprob_sop_aac = ODEProblem(sop_aac, [:A => 10.0, :C => 0.5], (0, 20), [:r => 0.05])
osol_sop_aac = solve(oprob_sop_aac, Rosenbrock23())
plot(osol_sop_aac, xlims=(0, 20), ylims=(0, 10), title="Second order production A + A --> C")


####################################################################################
# Saturated processes
####################################################################################

# Michaelis-Menten 
# In biochemistry, Michaelis-Menten (MM) kinetics describes the relationship between
# reaction rate and substrate concentration in enzyme-catalyzed reactions.
# It is:
#     - first-order at low substrate concentrations and
#     - zeroth-order at high substrate concentrations,
# showing the limited capacity of the catalyst.
mm_kin = @reaction_network begin
    (kp1, km1), E + S <--> ES
    k2, ES --> E + P
end
osys_mm_kin = convert(ODESystem, mm_kin)
equations(osys_mm_kin)
# Differential(t)(E(t)) ~ k2*ES(t) + km1*ES(t) - kp1*E(t)*S(t)
# Differential(t)(S(t)) ~ km1*ES(t) - kp1*E(t)*S(t)
# Differential(t)(ES(t)) ~ -k2*ES(t) - km1*ES(t) + kp1*E(t)*S(t)
# Differential(t)(P(t)) ~ k2*ES(t)
species(mm_kin)
# E(t)
# S(t)
# ES(t)
# P(t)
oprob_mm_kin = ODEProblem(mm_kin, [:E => 10.0, :S => 8.0, :ES => 0.0, :P => 0.0], (0, 40), [:kp1 => 0.1, :km1 => 0.1, :k2 => 0.2])
osol_mm_kin = solve(oprob_mm_kin, Rosenbrock23())
plot(osol_mm_kin, xlims=(0, 40), ylims=(0, 12), title="MM kinetics")
# P grows at a rate determined by the Michaelis-Menten equation!!!
# Et = E + ES = 10.0 + 0.0 = 10.0
# (k2 + km1)/kp1 = (0.2 + 0.1)/0.1 = 3.0   this will be Ks
# k2*Et = 0.2*10.0 = 2.0                   this will be vmax
# mm(X, v, K) = v*X/(X + K)
S = 0:0.1:5
vmax = 1.0
Ks = 1.0
plot(S, mm.(S, vmax, Ks), title="Michael-Menten kinetics", xlabel="[S]", ylabel="v", label="reaction rate", ylims=(0, 1))
# Ks is the slope in the beginning of the plot v versus S
# vmax is the horizontal asymptote

plot(S, mm.(S, 1.0, 1.0), title="Michael-Menten kinetics parameters", xlabel="[S]", ylabel="v", label="vmax=1, Ks=1", ylims=(0, 2), color=:blue, linestyle=:solid)
plot!(S, mm.(S, 2.0, 1.0), label="vmax=2, Ks=1", color=:red, linestyle=:dashdot)
plot!(S, mm.(S, 1.0, 2.0), label="vmax=1, Ks=2", color=:green, linestyle=:dashdotdot)
plot!(S, mm.(S, 1.0, 0.5), label="vmax=1, Ks=0.5", color=:blue, linestyle=:dot)


# Reverse Michaelis-Menten
S = 0:0.1:5
vmax = 1.0
Ks = 1.0
plot(S, mmr.(S, vmax, Ks), title="Repressing Michael-Menten kinetics", xlabel="[S]", ylabel="v", label="reaction rate", ylims=(0, 1))


# Hill-equation
L = 0:0.1:5
vmax = 1.0
Kl = 1.0
n = 5
plot(S, hill.(L, vmax, Kl, n), title="Hill equation", xlabel="[L]", ylabel="v", label="vmax=1.0, Kl=1.0, n=5", ylims=(0, 1))


# Reverse Hill-equation
L = 0:0.1:5
vmax = 1.0
Kl = 1.0
n = 5
plot(S, hillr.(L, vmax, Kl, n), title="Reverse Hill equation", xlabel="[L]", ylabel="v", label="vmax=1.0, Kl=1.0, n=5", ylims=(0, 1))


####################################################################################
# Example 9 (Repression of gene expression)
####################################################################################

repressor = @reaction_network begin
    @species R(t)=0 mRNA(t)=0
    hillr(R, vtranscr, Ki, 4), 0 --> mRNA
    vtransl, mRNA --> mRNA + R
    d, mRNA --> 0
    r, R --> 0
end
species(repressor)
# R(t)
# mRNA(t)
parameters(repressor)
# vtranscr
# Ki
# vtransl
# d
# r
osys_repressor = convert(ODESystem, repressor)
equations(osys_repressor)
# Differential(t)(R(t)) ~ -r*R(t) + vtransl*mRNA(t)
# Differential(t)(mRNA(t)) ~ Catalyst.hillr(R(t), vtranscr, Ki, 4) - d*mRNA(t)
oprob_repressor = ODEProblem(repressor, [], (0, 60), [:vtranscr => 0.5, :Ki => 1.0, :vtransl => 0.1, :d => 0.6, :r => 0.05])
osol_repressor = solve(oprob_repressor, Rosenbrock23())
plot(osol_repressor)


####################################################################################
# Logistic growth
####################################################################################

r = 1.0; K = 100.0
P = 0:1:120
plot(P, r.*(1 .- P./K).*P, title="Logistic growth", xlabel="P", ylabel="growth", label="r=1, K=100")

P0 = 2
t = 0:0.1:10
plot(t, K ./ (1 .+ (K-P0)/P0.*exp.(-r.*t)), title="Logistic growth equation", xlabel="t", ylabel="P", label="P")


####################################################################################
# Example 10 (Two species with competition)
####################################################################################

two_species_competition = @reaction_network begin
    @species A(t)=1.0 B(t)=5.0
    @parameters ra=0.2 rb=0.15 K=100 m=0.1
    ra*(1-(A+B)/K), A --> 2A
    rb*(1-(A+B)/K), B --> 2B
    m, (A, B) --> 0
end

osys_two_species_competition = convert(ODESystem, two_species_competition)
equations(osys_two_species_competition)
# Differential(t)(A(t)) ~ -m*A(t) + ra*A(t)*(1 + (-B(t) - A(t)) / K)
# Differential(t)(B(t)) ~ -m*B(t) + rb*B(t)*(1 + (-B(t) - A(t)) / K)
oprob_two_species_competition = ODEProblem(two_species_competition, [:B => 1.0], (0, 100), [])
osol_two_species_competition = solve(oprob_two_species_competition, Rosenbrock23())
plot(osol_two_species_competition)

oprob_two_species_competition = ODEProblem(two_species_competition, [], (0, 100), [])
osol_two_species_competition = solve(oprob_two_species_competition, Rosenbrock23())
plot(osol_two_species_competition)

oprob_two_species_competition = ODEProblem(two_species_competition, [:B => 1.0], (0, 100), [:rb => 0.199])
osol_two_species_competition = solve(oprob_two_species_competition, Rosenbrock23())
plot(osol_two_species_competition)


####################################################################################
# Example 11 (Pharmacokinetics)
####################################################################################

# When a drug is taken, it is first present in the gut (G), before it is
# absorbed and present in the blood (B), before it can enter the tissue (T)
# of interest.
pharkin = @reaction_network begin
    @species B(t)=0 T(t)=0
    k1, G --> B
    k2, B --> T
end

osys_pharkin = convert(ODESystem, pharkin)
equations(osys_pharkin)
# Differential(t)(B(t)) ~ -k1*B(t) + k1*G(t)
# Differential(t)(T(t)) ~ k1*B(t)
# Differential(t)(G(t)) ~ -k1*G(t)
oprob_pharkin = ODEProblem(pharkin, [:G => 5.0], (0, 50), [:k1 => 0.5, :k2 => 0.1])
osol_pharkin = solve(oprob_pharkin, Rosenbrock23())
plot(osol_pharkin)


####################################################################################
# Example 12 (Reactor with a dead zone)
####################################################################################

reactor2 = @reaction_network begin
    @species Ad(t)=0 B(t)=0
    @parameters V=100 k=0.3 r=0.05 q=1 c=.1
    # amount of A entering the bulk
    q*c, 0 --> Am
    # amount of A leaving the bulk
    q/(V*(1-f)), (Am, B) --> 0
    r, Am --> B # reaction
    # exchange bulk-dead zone
    (k/(V*(1-f)), k/(V*f)), Am <--> Ad
end

osys_reactor2 = convert(ODESystem, reactor2)
equations(osys_reactor2)
# Differential(t)(Ad(t)) ~ (k*Am(t)) / (V*(1 - f)) + (-k*Ad(t)) / (V*f)
# Differential(t)(B(t)) ~ r*Am(t)
# Differential(t)(Am(t)) ~ (k*Ad(t)) / (V*f) + (-q*Am(t)) / (V*(1 - f)) + (-k*Am(t)) / (V*(1 - f)) + c*q - r*Am(t)
oprob_reactor2 = ODEProblem(reactor2, [:Am => 0.0], (0, 500), [:f => 0.2])
osol_reactor2 = solve(oprob_reactor2, Rosenbrock23())
plot(osol_reactor2)


####################################################################################
# Example 13 (The SIR model)
####################################################################################

sir = @reaction_network begin
    @species S(t)=50 I(t)=5 R(t)=0
    beta, S + I --> 2I
    gamma, I --> R
end

osys_sir = convert(ODESystem, sir)
equations(osys_sir)
#  Differential(t)(S(t)) ~ -beta*S(t)*I(t)
#  Differential(t)(I(t)) ~ -gamma*I(t) + beta*S(t)*I(t)
#  Differential(t)(R(t)) ~ gamma*I(t)
# oprob_sir = ODEProblem(sir, [], (0, 120), [:beta => 0.005, :gamma => 0.05])
oprob_sir = ODEProblem(sir, [], (0, 30), [:beta => 0.03, :gamma => 0.3])
osol_sir = solve(oprob_sir, Rosenbrock23())
plot(osol_sir)
oprob_sir = ODEProblem(sir, [], (0, 30), [:beta => 0.01, :gamma => 0.3])
osol_sir = solve(oprob_sir, Rosenbrock23())
plot(osol_sir)


####################################################################################
# Example 14 (Butterfly model)
####################################################################################

butterfly = @reaction_network begin
    @species E(t)=10 C(t)=0 P(t)=0 B(t)=0
    @parameters f=1.5 s=0.05 m=0.05
    s, E --> C
    s, C --> P
    s, P --> B
    m, (E, C, P, B) --> 0
    f, B --> B + E
end

osys_butterfly = convert(ODESystem, butterfly)
equations(osys_butterfly)
# Differential(t)(E(t)) ~ f*B(t) - m*E(t) - s*E(t)
# Differential(t)(C(t)) ~ -m*C(t) + s*E(t) - s*C(t)
# Differential(t)(P(t)) ~ -m*P(t) + s*C(t) - s*P(t)
# Differential(t)(B(t)) ~ -m*B(t) + s*P(t)
oprob_butterfly = ODEProblem(butterfly, [], (0, 100), [])
osol_butterfly = solve(oprob_butterfly, Rosenbrock23())
plot(osol_butterfly)


####################################################################################
# Exercise 2
####################################################################################

fermentation = @reaction_network begin
    @species Y(t)=0.05 G(t)=150.0 E(t)=0.0
    @parameters mumax=0.1 KsG=5.0 KsE=30.0
    mm(Y, mumax, KsG), G + Y --> 2Y
    mm(E, mumax, KsE), G + E --> 2E
    2*50/100, 0 --> G
    2/100, Y --> 0
    2/100, E --> 0
end

osys_fermentation = convert(ODESystem, fermentation)
equations(osys_fermentation)
oprob_fermentation = ODEProblem(fermentation, [], (0, 50), [])
osol_fermentation = solve(oprob_fermentation, Rosenbrock23())
plot(osol_fermentation)


####################################################################################
# Exercise 3
####################################################################################

lotka_volterra2 = @reaction_network begin
    alpha, x --> 2x        # reproduction prey
    beta, y --> 0          # mortality pred
    gamma*(x*y)^k, x + y --> 1.1y  # predation
end

osys_lv2 = convert(ODESystem, lotka_volterra2, combinatoric_ratelaws=false)
equations(osys_lv2)
# Differential(t)(x(t)) ~ alpha*x(t) - gamma*x(t)*y(t)
# Differential(t)(y(t)) ~ -beta*y(t) + 0.10000000000000009gamma*x(t)*y(t)

oprob_lv2 = ODEProblem(lotka_volterra2, [:x => 1.0, :y => 1.0], (0.0, 50.0), [:alpha => 1.5, :beta => 0.5, :gamma => 0.5, :k => 0.5], combinatoric_ratelaws=false)
osol_lv2 = solve(oprob_lv2, Rosenbrock23())
plot(osol_lv2)


####################################################################################
# Exercise 4
####################################################################################

# https://www.palomar.edu/anthro/mendel/mendel_2.htm

zygotes = @reaction_network begin
    @parameters f=1.0 r=1/100
    @species AA(t)=100*f Aa(t)=0 aa(t)=100*(1-f)
    r, AA + AA --> AA + AA
    r, AA + Aa --> AA + Aa
    r, AA + aa --> Aa + Aa
    r, Aa + Aa --> 0.25*AA + 0.50*Aa + 0.25*aa
    r, aa + Aa --> Aa + aa
    r, aa + aa --> aa + aa
end

species(zygotes)
reactions(zygotes)

osys_zygotes = convert(ODESystem, zygotes, combinatoric_ratelaws=false)
equations(osys_zygotes)
# Differential(t)(AA(t)) ~ 0.25r*(Aa(t)^2.0) - r*AA(t)*aa(t)
# Differential(t)(Aa(t)) ~ -1.5r*(Aa(t)^2.0) + 2r*AA(t)*aa(t)
# Differential(t)(aa(t)) ~ 0.25r*(Aa(t)^2.0) - r*AA(t)*aa(t)

oprob_zygotes = ODEProblem(zygotes, [], (0.0, 100.0), [:f => 0.7], combinatoric_ratelaws=false)
osol_zygotes = solve(oprob_zygotes, Rosenbrock23())
plot(osol_zygotes)
