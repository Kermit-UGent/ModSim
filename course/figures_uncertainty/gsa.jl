using GlobalSensitivity, Plots, Statistics, DifferentialEquations, Catalyst

# code for the examples in the chapter uncertainty

sir = @reaction_network begin
	@species S(t)=50 I(t)=5 R(t)=0
	β, S + I --> 2I
	γ, I --> R
end

sir_sys = convert(ODESystem, sir, combinatoric_ratelaws=false)

prob = ODEProblem(sir_sys, [], [0, 100], [0.005, .01])
plot(solve(prob), lw=2)

# Now, let's create a function that takes in a parameter set and calculates the maximum of 
# the predator population and the average of the prey population for those parameter values.
# To do this, we will make use of the remake function, which creates a new ODEProblem, and 
# use the p keyword argument to set the new parameters:
f1 = function (p)
    # β, γ = p
    prob1 = remake(prob; p = p)
    sol = solve(prob1, Tsit5())
    [mean(sol[1, :]), maximum(sol[2, :])]
end

p_range = [[1e-5, 1], [1e-5, 1]] # parameter ranges
m = gsa(f1, Morris(total_num_trajectory = 1000, num_trajectory = 150), p_range)

m.means
m.variances

scatter(
    m.means[1, :], m.variances[1, :], series_annotations = [:β, :γ], color = :gray)

m = gsa(
    f1, Sobol(),
    p_range, samples = 1000)

m.S1, m.ST
