using Distributions
using StatsPlots
using Turing

p_prior_values = rand(Uniform(0, 1), 10000)

histogram(p_prior_values, normalized=true, bins=10,
    xlabel="p value", ylabel="frequency (normalized)", legend=false,
    alpha=0.6, color="red", size=(600, 350))

histogram(p_prior_values, bins=10,
    xlabel="p value", ylabel="frequency (absolute)", legend=false,
    alpha=0.6, color="red", size=(600, 350))