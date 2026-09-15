using Catalyst, GlobalSensitivity, Statistics, DifferentialEquations, Plots

# use the SIR model with parameters β and γ for the analysis.
sir = @reaction_network begin
	@species S(t)=50 I(t)=5 R(t)=0
	β, S + I --> 2I
	γ, I --> R
end
# convert the reaction network to an ODE system
sir_sys = convert(ODESystem, sir, combinatoric_ratelaws=false)

prob = ODEProblem(
    sir_sys,
    [],
    [0, 100], # range for time t
    [0.005, .01] # β, γ
)

# Now, let's create a function that takes in a parameter set 
# and calculates the maximum of the infected population and 
# the average of the susceptible population for those parameter values.
# To do this, we will make use of the remake function, which 
# creates a new ODEProblem, and uses the p keyword argument 
# to set the new parameters:
f1 = function (p)
    # β, γ = p
    prob1 = remake(prob; p = p)
    sol = solve(prob1, Tsit5(); saveat = 0.1)
    [mean(sol[1, :]), maximum(sol[2, :])]
end

p_range = [[1e-5, 1], [1e-5, 1]] # parameter ranges

# compute the GSA using Sobol
m = gsa(f1, Sobol(), p_range, samples = 1000)

m.S1 # first order effects
m.ST # total effects

m.S1[1, :] # first order effect of β and γ on the mean susceptible population
m.S1[2, :] # first order effect of β and γ on the maximum infected population

scatter(m.S1[1, :], ["β", "γ"], label = "Mean susceptible population (first order effects)", ylabel = "Parameter", xlabel = "Effects")
scatter!(m.ST[1, :], ["β", "γ"], label = "Mean susceptible population (total effects)")
scatter!(m.S1[2, :], ["β", "γ"], label = "Maximum infected population (first order effects)")
scatter!(m.ST[2, :], ["β", "γ"], label = "Maximum infected population (total effects)")
savefig("sobol.pdf")

# compute the GSA using Morris
m = gsa(f1, Morris(num_trajectory=50), p_range)

m.means
m.variances

# TODO check if abs is needed? by definition we should not have negative values
scatter([abs(m.means[1, 1])], [m.variances[1, 1]], label = "beta, Mean susceptible population", xlabel = "μ*", ylabel = "σ")
scatter!([abs(m.means[1, 2])], [m.variances[1, 2]], label = "gamma, Mean susceptible population")
scatter!([abs(m.means[2, 1])], [m.variances[2, 1]], label = "beta, Maximum infected population")
scatter!([abs(m.means[2, 2])], [m.variances[2, 2]], label = "gamma, Maximum infected population")
savefig("morris.pdf")

# Morris paths visualisation
# https://uqpyproject.readthedocs.io/en/latest/sensitivity/morris.html
# unit square
plot([0, 1], [0, 0], color=:black)
plot!([0, 1], [1, 1], color=:black)
plot!([0, 0], [0, 1], color=:black)
plot!([1, 1], [0, 1], color=:black)
# lines
plot!([0, 0.2], [0.2, 0.2], color="CornflowerBlue", lw=3)
plot!([0.2, 0.2], [0.2, 0.8], color="CornflowerBlue", lw=3)
plot!([0.2, 0.8], [0.8, 0.8], color="CornflowerBlue", lw=3)

plot!([0.5, 0.5], [0.9, 0.4], color="DarkOrange", lw=3)
plot!([0.5, 0.7], [0.4, 0.4], color="DarkOrange", lw=3)
plot!([0.7, 0.7], [0.4, 0.2], color="DarkOrange", lw=3)

plot!([0.3, 0.9], [0.6, 0.6], color="LightPink", lw=3)
plot!([0.9, 0.9], [0.6, 0.1], color="LightPink", lw=3)
plot!(xlabel = "X₁", ylabel = "X₂", xlims = (-0.2, 1.2), ylims = (-.2, 1.2), size = (500, 500), legend=false, xticks = 0:1, yticks = 0:1)

savefig("randomtrajectories.pdf")

# a plot showcasing uncertainty
plot()
plot!(x -> sin(x), 0, 2pi, label = "sin(x)", lw = 2)
plot!(x -> sin(x)+0.3x, 0, 2pi, label = "upper bound", lw = 2)
plot!(x -> sin(x)-0.3x, 0, 2pi, label = "lower bound", lw = 2)
plot!(x -> sin(x)-0.3x, 0, 2pi, label = "error", lw = 0, fillrange=x -> sin(x)+0.3x, fillalpha=0.1, fillcolor=:red, z_order=:back)
plot!(xlabel="time", ylabel="response")

savefig("uncertainplot.pdf")
