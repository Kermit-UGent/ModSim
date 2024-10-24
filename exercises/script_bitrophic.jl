
using Catalyst
using DifferentialEquations, Plots

# % parameter definitions
# theta = 0.2;      % [/d]
# k = 4000;         % [kg/ha]
# f = 0.001;        % [ha/(kg d)]
# phi = 0.2;        % efficiency ratio
# mu = 0.1;         % mortality ratio [/d]
# %mu = 0.07 -> for figure in solution part
# p = 3;            % insecticed factor

# % initial conditions
# C0 = 100;         % initial amount of crops [kg/ha]
# A0 = 0.5;         % initial amount of insects [kg/ha]

# % simulation
# simtime = 200;     % simulation time [d]
# sampletime = 0.25;


# Define all reactions
bitrophic_model = @reaction_network begin
    θ*(1-C/k), C --> 2C
    f, C + A --> (1 + ϕ)*A
    (1+p)*μ, A --> ∅
end

species(bitrophic_model)
# C(t)
# A(t)

parameters(bitrophic_model)
# θ
# k
# f
# ϕ
# p
# μ

# Convert to an ODE system
osys  = convert(ODESystem, bitrophic_model)

# Check out the diff. eqns. equations
equations(osys)
# Differential(t)(C(t)) ~ -f*A(t)*C(t) + C(t)*(1 + (-C(t)) / k)*θ
# Differential(t)(A(t)) ~ (-1 - p)*A(t)*μ + f*A(t)*C(t)*ϕ

# Set initial condition for ODE problem
u0 = [:C => 100, :A => 0.5]
# Set time span for simulation
tspan = (0.0, 200.0)
# Set parameters for ODE problem
params = [:θ => 0.2, :k => 4000.0, :f => 0.001, :ϕ => 0.2, :p => 3, :μ => 0.1]
# Create ODE problem
oprob = ODEProblem(bitrophic_model, u0, tspan, params)

osol = solve(oprob, Tsit5(), saveat=0.5)
plot(osol)

ss_guess = [:C => osol[:C][end], :A => osol[:A][end]]  # last value as guess
# :C => 1999.9428643956621
# :A => 99.93555955847543
Ceq, Aeq = solve(SteadyStateProblem(ODEProblem(bitrophic_model, ss_guess, tspan, params, combinatoric_ratelaws=false)))
# 1999.9999999999998
# 100.00000000000006

# θ, k, f, ϕ, p, μ

################################################################
## Sensitivity analysis
################################################################

using ForwardDiff

function bitrophic_model_sim(params)
	θ, ϕ, p = params
    u0 = [:C => 100, :A => 0.5]
    tspan = (0.0, 200.0)
    params = [:θ => θ, :k => 4000.0, :f => 0.001, :ϕ => ϕ, :p => p, :μ => 0.1]
    oprob = ODEProblem(bitrophic_model, u0, tspan, params, combinatoric_ratelaws=false)
	osol = solve(oprob, Tsit5(), saveat=0.5)
	return osol
end

bitrophic_model_sim_C(params) = bitrophic_model_sim(params)[:C]
bitrophic_model_sim_A(params) = bitrophic_model_sim(params)[:A]

t_vals = 0:0.5:200.0
θ = 0.2
ϕ = 0.2
p = 3

C_sim = bitrophic_model_sim_C([θ, ϕ, p])
A_sim = bitrophic_model_sim_A([θ, ϕ, p])

sens_C = ForwardDiff.jacobian(bitrophic_model_sim_C, [θ, ϕ, p])
sens_A = ForwardDiff.jacobian(bitrophic_model_sim_A, [θ, ϕ, p])

sens_C_on_θ = sens_C[:, 1]
sens_C_on_ϕ = sens_C[:, 2]
sens_C_on_p = sens_C[:, 3]

sens_A_on_θ = sens_A[:, 1]
sens_A_on_ϕ = sens_A[:, 2]
sens_A_on_p = sens_A[:, 3]

sens_C_on_θ_rel = sens_C_on_θ .* θ ./ C_sim
sens_C_on_ϕ_rel = sens_C_on_ϕ .* ϕ ./ C_sim
sens_C_on_p_rel = sens_C_on_p .* p ./ C_sim

sens_A_on_θ_rel = sens_A_on_θ .* θ ./ A_sim
sens_A_on_ϕ_rel = sens_A_on_ϕ .* ϕ ./ A_sim
sens_A_on_p_rel = sens_A_on_p .* p ./ A_sim

plot(t_vals, [sens_C_on_θ_rel, sens_A_on_θ_rel], title="Normaized sensitivities", label=["C on θ" "A on θ"], xlabel="Time (days)")
sens_C_on_θ_rel[end]
sens_A_on_θ_rel[end]
sens_A_on_θ_rel[end] / sens_C_on_θ_rel[end]
"""
In steady state, the maximum growth rate θ has very little effect on C,
because if θ would incease then more crops C would grow but then
they would at the same time be eaten up by the insects. You can also
see that θ has a positive effect on A, actually a factor of about 140
larger than on C.
"""

plot(t_vals, [sens_C_on_ϕ_rel, sens_A_on_ϕ_rel], title="Normaized sensitivities", label=["C on ϕ" "A on ϕ"], xlabel="Time (days)")
"""
The parameter ϕ can be considered as the yield factor for the insects A.
In steady state, hence, its influence is positive, and the more insects,
the less the number of crops. That is the reason why the effect of ϕ is
negative for C.
"""

plot(t_vals, [sens_C_on_p_rel, sens_A_on_p_rel], title="Normaized sensitivities", label=["C on p" "A on p"], xlabel="Time (days)")
"""
The parameter p is related to the amount of pesticides used.
The more pesticides, the more insects will die and the more
crops will be able to grow. Hence, p has an positive influence
on C.
"""

plot(t_vals, [sens_C_on_θ_rel, sens_C_on_ϕ_rel, sens_C_on_p_rel], title="Normaized sensitivities", label=["C on θ" "C on ϕ" "C on p"], xlabel="Time (days)")

51

# ################################################################
# ## Sensitivity analysis
# ################################################################

# using SciMLSensitivity

# oprob_sens = ODEForwardSensitivityProblem(oprob.f, [100, 0.5], tspan, [0.2, 4000.0, 0.001, 0.2, 3, 0.1])
# # oprob.f        : the ODE function
# # [100, 0.5]     : initial conditions, the order is [C, A], see species(bitrophic_model)
# # tspan          : the time span
# # [0.2, 4000.0, 0.001, 0.2, 3, 0.07] : the parameter values, the order is [θ, k, f, ϕ, p, μ], see parameters(bitrophic_model)

# osol_sens = solve(oprob_sens, Tsit5(), saveat=0.5)
# u, dp = extract_local_sensitivities(osol_sens)

# sol_C = u[1,:]     # solution for C
# sol_A = u[2,:]     # solution for A

# sens_theta = dp[1]'  # Select sensitivity of C and A to θ
# sens_phi   = dp[4]'  # Select sensitivity of C and A to ϕ
# sens_p     = dp[5]'  # Select sensitivity of C and A to p

# # Plot sensitivity of C and A to θ
# plot(osol_sens.t, sens_theta, title="Absolute sensitivity of C and A to θ", label=["dC/dθ" "dA/dθ"])
# sens_theta[end, :]

# # Plot sensitivity of C and A to ϕ
# plot(osol_sens.t, sens_phi, title="Absolute sensitivity of C and A to ϕ", label=["dC/dϕ" "dA/dϕ"])
# sens_phi[end, :]

# # Plot sensitivity of C and A to p
# plot(osol_sens.t, sens_p, title="Absolute sensitivity of C and A to p", label=["dC/dp" "dA/dp"])
# sens_p[end, :]

# trsens_theta = sens_theta.*0.2./[sol_C sol_A]
# trsens_phi = sens_phi.*0.2./[sol_C sol_A]
# trsens_p = sens_p.*3.0./[sol_C sol_A]

# # Absolute sensitivity of C and A to θ
# plot(osol_sens.t, trsens_theta, title="Absolute sensitivity of C and A to θ", label=["dC/dp * θ / C" "dA/dp * θ / A"])
# trsens_theta[end, :]

# # Total relative sensitivity of C and A to ϕ
# plot(osol_sens.t, trsens_phi, title="Total relative sensitivity of C and A to ϕ", label=["dC/dϕ * ϕ / C" "dA/dϕ * ϕ / A"])
# trsens_phi[end, :]

# # Total relative sensitivity of C and A to p
# plot(osol_sens.t, sens_p.*3.0./u', title="Total relative sensitivity of C and A to p", label=["dC/dp * p / C" "dA/dp * p / A"])
# trsens_p[end, :]


################################################################
## Parameter uncertainty
################################################################

using Measurements

# ± with \pm
params_uncert = [:θ => 0.20±0.02, :k => 4000.0, :f => 0.001, :ϕ => 0.2±0.02, :p => 3.0±0.2, :μ => 0.1]
params_uncert = [:θ => 0.20±0.02, :k => 4000.0, :f => 0.001, :ϕ => 0.2, :p => 3.0, :μ => 0.1]
params_uncert = [:θ => 0.20, :k => 4000.0, :f => 0.001, :ϕ => 0.2±0.02, :p => 3.0, :μ => 0.1]
params_uncert = [:θ => 0.20, :k => 4000.0, :f => 0.001, :ϕ => 0.2, :p => 3.0±0.2, :μ => 0.1]

oprob_uncert = ODEProblem(bitrophic_model, u₀, tspan, params_uncert, combinatoric_ratelaws=false)

osol_uncert = solve(oprob_uncert, Tsit5(), saveat=2.0)

plot(osol_uncert)


