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

#=
All processes in which we have two possible outcomes -heads or tails
in our case-, and some probability p of success -probability of heads-,
these are called Bernoulli trials.
=#

N = 100
p = 0.5

# Geeft the PDF weer van het aantal HEADS indien we N keer opgooien.

bar(Binomial(N,p), xlim=(0, N), label=false,
    xlabel="Succeses", ylabel="Probability", title="Binomial distribution", color="green",
    alpha=0.8, size=(600, 350)) 

#=
Bernoulli(p)  is een PDF

The Bernoulli distribution, named after Swiss mathematician Jacob Bernoulli,
is the discrete probability distribution of a random variable which takes the
value 1 with probability p and the value 0 with
probability q = 1 − p.
=#

rand(Bernoulli(p), 10)       # returns a vector of length 10 with random 0's and 1's
println(rand(Bernoulli(p), 100))
# Bool[0, 1, 1, 0, 0, 1, 0, 1, 1, 1]


######################################################################
# Prior distribution of p is Uniform
######################################################################

@model coinflip(y) = begin
    # Our prior belief about the probability of heads in a coin.
    p ~ Uniform(0, 1)

    # The number of observations.
    N = length(y)
    for n in 1:N
        # Heads or tails of a coin are drawn from a Bernoulli distribution.
        y[n] ~ Bernoulli(p)
    end
end

outcome = [0, 1, 1, 0, 1, 0, 0, 1, 1, 1]
# outcome = [0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1]

# Settings of the Hamiltonian Monte Carlo (HMC) sampler.
iterations = 1000
ϵ = 0.05
τ = 10

# Start sampling.
chain = sample(coinflip(outcome[1]), HMC(ϵ, τ), iterations, progress=false)

histogram(chain[:p], legend=false, xlabel="p", ylabel="Probability", title="Posterior distribution of p after getting tails", alpha=0.7, xlim=(0, 1), size=(600, 350))

samples = []
for i in 2:10   
    global chain_
    chain_ = sample(coinflip(outcome[1:i]), HMC(ϵ, τ), iterations, progress=false)
    push!(samples, chain_[:p])
end

plots = [histogram(samples[i], normalized=true, legend=false, bins=10, title="$(i+1) outcomes", titlefont = font(8)) for i in 1:9];
plot(plots...)


######################################################################
# Prior distribution of p is Beta(2, 2)
######################################################################

plot(Beta(2,2), legend=false, xlabel="p", ylabel="Probability", title="Prior distribution for p", fill=(0, .5,:dodgerblue), ylim=(0,2), size=(600, 350))


@model coinflip_beta_prior(y) = begin
    # Our prior belief about the probability of heads in a coin.
    p ~ Beta(2, 2)

    # The number of observations.
    N = length(y)
    for n in 1:N
        # Heads or tails of a coin are drawn from a Bernoulli distribution.
        y[n] ~ Bernoulli(p)
    end
end

samples_beta_prior = []
for i in 2:10   
    global chain__
    chain__ = sample(coinflip_beta_prior(outcome[1:i]), HMC(ϵ, τ), iterations, progress=false)
    push!(samples_beta_prior, chain__[:p])
end


plots = [histogram(samples_beta_prior[i], normalized=true, legend=false, bins=10, title="$(i+1) outcomes", titlefont = font(10), color="red", alpha=0.6) for i in 1:9];
plot(plots...)