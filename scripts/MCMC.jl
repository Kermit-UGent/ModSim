### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 103e5ba0-cfdc-11ee-13b1-cf53dfdd9a3b
# ╠═╡ skip_as_script = true
#=╠═╡
begin
    using Pkg
	Pkg.activate("..")	
end
  ╠═╡ =#

# ╔═╡ d778ce04-8df4-42ef-95f1-9cf4880e0420
using Turing, StatsPlots, Distributions

# ╔═╡ 0350aa5b-d105-4dfa-a454-59873672b3a0
using Plots, PlutoUI, LaTeXStrings, LinearAlgebra, Random

# ╔═╡ e87dd9a4-f4ed-46d7-9f9d-dae0cadb7d05
md"""
# Bayesian reasoning and advanced sampling methods

In the previous chapter, we explored how one can build joint probability distributions and generate samples from them. This chapter explores inference: we will fix one or multiple variables, which represents observing them, and will sample from the conditional distribution. 
"""

# ╔═╡ 0f6cf892-b218-4740-a298-43f47e51acae
md"""
## Introduction to Bayesian reasoning

In the previous chapter, we have seen how to build complex distributions from simple ones. We stressed that our probabilistic programming framework is extremely flexible and allows us to write general code to compute things or even include sophisticated models, such as the ordinary differential equations we will use in subsequent chapters. In effect, we have built a simulator that can generate data $\mathcal{D}$ from one or more inputs, parameters, or hidden variables $\theta$. The simulator is stochastic, so it can be seen as a probability distribution. We will denote the simulator as $P(\mathcal{D}\mid \theta)$, i.e., we convert parameters and inputs into data. The problem we want to study in this chapter is that of *inference*: given that we observe data, what can we say about the parameters, i.e., can we find $P(\theta\mid \mathcal{D})$? We have already seen how this distribution can be attained using Bayes' theorem:

$$P(\theta\mid\mathcal{D}) = \frac{P(\mathcal{D}\mid\theta)\,P(\theta)}{P(\mathcal{D})}\,,$$

of which each component represents an important concept:
-  $P(\mathcal{D} \mid \theta)$ is called the *likelihood distribution*, it represents the probability of the data, given the parameters;
-   $P(\theta)$ is the *prior distribution* representing our beliefs and knowledge about the parameters before we observe the data;
-   $P(\mathcal{D})$ is the model *evidence*, representing how probable this model or simulator can generate this data over all possible values of the parameters;
-  $P(\theta\mid \mathcal{D})$ is the desired *posterior distribution*, representing our beliefs about the values of the parameters, given that we have observed data.

So, Bayes' theorem for data in words is

$$\text{posterior} = \frac{\text{likelihood}\times\text{prior}}{\text{evidence}}\,.$$

As noted earlier, the model evidence is not obvious to compute in general, as we have to integrate over all possible values of $\theta$, which may be high-dimensional to get this. However, given that we have observed the data, $P(\mathcal{D})$ can be seen as a constant that ensures that the posterior is a normalized distribution. So, we have

$$P(\mathcal{D}\mid \theta) \propto P(\theta\mid \mathcal{D})\,P(\theta)\,.$$

To get the *maximum a posteriori probability* value of $\theta$, we have to solve the following optimization problem:

$$\theta^\star_{MAP}=\text{arg max}_{\theta}\, P(\theta\mid \mathcal{D})\,P(\theta)\,.$$

The MAP estimate balances the prior on $\theta$ with the likelihood of the data. The final estimate of the parameters or state will be something in between both, depending on their relative strength.

Compare the MAP estimate with the *maximum likelihood* estimator for $\theta$ you have seen in statistics (e.g., to derive the estimator of the mean of a normal distribution), which is

$$\theta^\star_{ML}=\text{arg max}_{\theta}\, P(\mathcal{D} \mid \theta)\,.$$

It only depends on the likelihood, not on any prior beliefs on $\theta$. 

Computations are usually performed after log-transforming the distributions. This is because high-dimensional distributions can represent many datasets, so any one point, the evaluation value will be extremely small. Furthermore, many probability distributions have very friendly forms if you take their logarithm and discard constants that don't depend on the state variables; for example, the log-PDF of a normal distribution is given by

$$-\frac{(x-\mu)^2}{\sigma^2} + \text{cst}\,,$$
where the constant term does not depend on $x$ but merely ensures normalization. As the logarithm is a monotonically increasing transformation, it does not change optimization problems: minimizing the log-likelihood gives the same estimator as minimizing the likelihood and ditto for the posterior. Note that this is the reason why maximum-likelihood estimation is equivalent to least-squares under the assumption of homoskedastic normally distributed noise.

The MAP and ML are point estimators; they represent only a single value. In statistics, one often represents an estimation by a $(1-\alpha)$ (e.g., 95%) *confidence interval*. These are two bounds $[l, u]$ of which the data indicates it is very likely to contain the true parameter.

> The proper interpretation under frequentist statistics is that if similar data is collected many times and one uses suitable methods to create the confidence intervals, the true, unknown parameter $\theta$ would fall in these intervals $(1-\alpha)$% of the times.

Thinks are much more straightforward using Bayesian reasoning: we have the full posterior distribution $P(\theta\mid \mathcal{D})$, which contains all information about $\theta$ given the data, the model and our initial beliefs. We can use quantiles to generate an interval, for example of 95%. This directly represents our belief what $\theta$ could be. In Bayesian statistics, these types of intervals are called *credibility intervals*. We can ask questions like whether it is likely that $\theta>0.8$ etc! 

You might wonder how we obtained the prior $P(\theta)$ in the first place. The prior is the main reason why some statisticians dislike Bayesian statistics, which is that they are inherently subjective. We can identify different types of priors based on the amount of information that it contains:
• **Informative prior:** This prior incorporates existing knowledge or beliefs about the parameter, potentially biasing the posterior towards specific values. For example, if we know a historical germination rate for seeds, we could use that information to inform the prior distribution for the current experiment.
• **Weakly informative prior:** This prior acknowledges some prior knowledge but avoids strong biases. It might specify a range or a general shape for the parameter distribution, allowing the data to significantly influence the posterior.
• **Diffuse or uninformative prior:** This prior represents minimal or no prior knowledge about the parameter. It is often chosen as a uniform distribution across the possible parameter values, letting the data dictate the posterior distribution entirely.

Let us revise the pepper seed germination example from the previous chapter to make Bayesian reasoning more concrete. Here, we used a Binomial distribution Binom(10, $p$) to model the number of seeds that could germinate (the likelihood). Suppose we want to infer the germination probability $p$ given that eight of the ten seeds germinated. We can plot this likelihood function by plotting the PMF $p_X(k)=P(X=k)$  of this distribution for values of $p$ keeping $k$ fixed to 10. Note that our prior (the distribution of $p$ over different strains) was given by a Beta(8, 3) distribution, which represents a bump at around $p\approx 0.7272$ with a considerate margin for higher or lower values. When we look at the likelihood of $p$ given that $X=8$, we see a sharp peak at $p=0.8$. The data-generating model assumes hence tell us that $p=0.8$ (the maximum-likelihood estimation for $p$) is the most likely, with higher or lower values being considerably less likely. It is almost impossible that eight seeds have germinated given that $p < 0.5$. The posterior (footnote of these together) is a distribution that lies between both distributions and has a peak at $p^\star_\text{MAP}=0.79$. As we see, the posterior is more conservative!
"""

# ╔═╡ eb37dd4b-3af0-4534-92c2-e495b83024be
md" Number of seeds germinated : $(@bind k Slider(0:10, default=9, show_value=true))"

# ╔═╡ 29928abf-f010-438e-9b60-e6673a51ddf5
k

# ╔═╡ 768ea152-00a5-4409-86eb-9a1c554eb3b2
mean(Beta(8, 3))  

# ╔═╡ 959cd079-59a4-42f5-a868-4cc0675c694b
seed_likelihood = p -> pdf(Binomial(10, p), k)

# ╔═╡ 39cb8a88-9bfe-4b51-8f2f-b89b45dccc95
seed_prior = p -> pdf(Beta(8, 3), p)

# ╔═╡ 667f2069-fd32-4c3e-aebb-b75a38bf84b9
seed_posterior = p -> pdf(Beta(8+k, 3+10-k), p)

# ╔═╡ 472cea3d-12a8-449a-8e87-6c1b5aaa1fd7
md"In this case, we can compute the exact posterior because the beta distribution and the binomial distribibition are *conjugate*, meaning that a binomial likelihood on a beta prior results in a beta posterior with updated parameters. Though conjugated priors are important in Bayesian statistitics, we won't make use of them as we will use a sampling approach." 

# ╔═╡ 3ef47022-f9e3-434e-87a3-20b1870a7c24
p_ML = argmax(seed_likelihood, 0:0.01:1)

# ╔═╡ 6004b26a-8f10-492d-bdd2-1d317d8bb394
p_MAP = argmax(seed_posterior, 0:0.01:1)

# ╔═╡ 2ab01769-42cd-44e1-8d4c-2c8b477c732a
md"Let us consider a bit of a more of a classical statistical example. Suppose we have a vector $\mathbf{y}$ with values i.i.d. distrubuted from a normal distribution $y_i\sim$Normal($\mu$, $\sigma$), with unknown mean and standard deviation."

# ╔═╡ 4d6d6c9e-c2a8-4c7d-9d18-89d2e9d0c213
y = [7.2, 8.3, 5.4, 9.8, 7.9]

# ╔═╡ ba06e8a5-c142-4963-97b4-783ed4c706bf
mean(y), std(y)  # sample mean and standard deviation

# ╔═╡ f44abcc5-9010-4d10-b722-ef8147c8fd66
md"""Your basic statistics course would course have estimators for $\mu$ and $\sigma$, the sample mean and the square root of the sample variance, respectively. Let us use a Bayesian reasoning and place a prior on these two unknown parameters:

$$\mu \sim \text{Normal}(0, 20)$$

$$\sigma \sim \text{InverseGamma}(1/10)$$

which are chosen to be "reasonable". 
"""

# ╔═╡ 0c29af96-949d-47b3-869f-274ed39c9d25
md"The full, hierarchical model can be implemented as a Turing model."

# ╔═╡ 2cbe07e6-8cf0-4dc4-8aea-f56059ffa367
md"Note that this is now a distribution over 2 + $n$ variables as $\mathbf{y}$ can be of arbirary length. By giving $y$ as an argument, we fix the variables, though the parameters $\mu$ and $\sigma$ remain random variables."

# ╔═╡ afb7bdca-a8de-4342-9212-5e81fa47963a
md"The prior over $\mu$ and $\sigma$ is a product distribution of a normal (centered around 0 and large standard deviation) and an inverse Gamma."

# ╔═╡ 1fbd13c1-cc93-4e5b-98f1-ace195d20f97
md"The likelihood of $\mu$ and $\sigma$ can be obtainded from the PDF of a normal distirbution over the parameters:

$$L(\mu, \sigma\mid\mathbf{y}) = \prod_{i=1}^n \frac{1}{\sqrt{2\pi}\sigma}e^{-\frac{(y_i-\mu)^2}{2\sigma^2}}\,.$$"

# ╔═╡ 213465af-f980-4848-9e5a-0cb5429f4af8
md"The posterior is the product of the prior and the likelihood."

# ╔═╡ 7f2a7d88-442f-431e-8e09-42082bda8d24
md"The posterior is slightly shifted compared to the likelihood, though both are quite close! We see that the posterio distribution is centred around the sample mean and standard deviation. The posterior identifies which parameters of the normal distribution have likely given rise to the data. At sight, the mean $\mu$ is likely between 6 and 9, while the standard deviation is likely to be between 1 and 2. Using numerical methods, we can identify the MAP estimator of the parameters and credibility intervals. Inference in two dimensions is hence still quite tractable. However, in what follows, we will use sampling to get observations from the posterior, which we can use to learn everything we want to know about the distributions and hence our parameters or variables of interest."

# ╔═╡ 985c1d3e-bc02-4dc5-9a28-d0ec9edf5cc9
md"""
## Rejection sampling

> A group of little kids are playing at the beach. They have drawn a circle in the sand and throw little pebbles, counting which fell into its boundaries.

Given some constraints, rejection sampling is a simple framework to sample from reasonably complex distributions. We have already encountered rejection sampling, without calling it as such, in our Monte-Carlos algorithm for estimating $\pi$. Here, we discuss the univariate version. The generalization to multivariate distributions is trivial.

As we have argued earlier, we often know the probability density or mass function up to the normalization constant $Z$, which requires integration:

$$p(x) = \frac{1}{Z}\tilde{p}(x)\,,$$

where $\tilde{p}(x)$ is the unnormalized PDF, from which we could theoretically get the normalization constant as

$$Z = \int_{-\infty}^\infty\tilde{p}(x)\mathrm{d}x\,.$$

For example, remember that when we obtain the unnormalized posterior distribution by multiplying the prior with the likelihood with the hard-to-compute evidence. The normalization is constant is merely a scaling factor. 

Rejection sampling uses a *proposal distribution* $q(x)$ from which it is easy to sample, e.g., a uniform or normal distribution. Importantly, the proposal distribution and the target distribution need to have the same support. Next, we need to have a constant $M$, such that $Mq(x)\ge \tilde{p}(x)$ for all $x$, meaning that we rescale the proposal distribution such that it lies above the target distribution. Rejection sampling is done in three steps using two random samples:

1. generate a sample from the proposal distribution: $x \sim q(x)$;
2. compute the *acceptance probability* $\alpha=\frac{\tilde{p}(x)}{Mq(x)}$;
3. generate a uniformly-distributed number: $u\sim$Unif(0,1):
   - if $u\le \alpha$, then **accept** $x$ as a sample of $p(x)$
   - if $u>\alpha$, **reject** $x$ and run the algorithm.

The acceptance probability, $\alpha = \frac{\tilde{p}(x)}{Mq(x)}$, ensures that the accepted samples are drawn from the target distribution $p(x)$. When a sample $x$ is generated from the proposal distribution $q(x)$ and an acceptance probability $\alpha$ is computed, it effectively scales the probability of accepting $x$ based on how likely it is under the target distribution relative to the proposal distribution. Rejection sampling is an exact method; accepted samples follow the target distribution.

We might fail at generating a sample (rejection), so we must rerun it until it accepts a sample. The probability of accepting a sample is 

$$P(\text{accept}) = \int_{-\infty}^\infty\alpha q(x)\mathrm{d}x=\int_{-\infty}^\infty\frac{\tilde{p}(x)}{Mq(x)} q(x)\mathrm{d}x=\frac{Z}{M}\,,$$

meaning that we need on average $M/Z$ throws to generate a sample from $p(x)$ (this follows a geometric distribution). It should be clear that very low acceptance rates are the main potential reasons why rejection sampling won't scale for certain problems. A low acceptance rate could be due to:
- a poor match between $q(x)$ and $p(z)$ (try to find a proposal distribution that closely resembles the shape of the target distribution);
- the constant $M$ being too high (try to find as small a value of $M$ as possible).

One modification to deal with low acceptance rates is *adaptive rejection sampling*, which constructs a sampling distribution $q(x)$ on the fly based on the generated samples. For example, you could make $q(x)$ piecewise linear and, with every draw, update these pieces to better match the target distribution.
"""

# ╔═╡ a0a1deec-75b6-4061-aa32-ad5e940d4123
md"M : $(@bind M Slider(2.1:0.2:10, show_value=true))"

# ╔═╡ d91cea07-fb14-4b93-8ad1-31fb8f3afe8d
md" n throws : $(@bind n_rejection_sampling Slider(10:10:5000, default=100, show_value=true))"

# ╔═╡ d669d181-dc13-4921-9aca-88acae1197eb
md"""
## Markov chain Monte Carlo methods

Rejection sampling and similar methods scale poorly to problems with three dimensions and more. Here, we discuss a much more general and powerful framework called Markov chain Monte Carlo (MCMC). These methods have their origins in physics and it was since the 1980's that they significantly impacted statistical modeling. In contrast to the earlier seen sampling methods, MCMC generates samples by performing a biased random walk in the sampling space, which converges (hopefully!) to the target distribution. These methods produces 'chains' of samples that are not independent but posses a certain memory that slowly decays.  When making inference, one has to be cautious the the effective sample size using MCMC can be much lower than when the samples would have been i.i.d. For this reason, the samples that are obtained by MCMC methods are often called *pseudosamples*.
"""

# ╔═╡ 9c48f6a4-db8c-4e15-8486-4dd9b281a80d
md"""
### Brief refresher of Markov chains

Before delving into the specifics of the MCMC, let us recap the basics of elementary Markov chains. A Markov chain is a mathematical model that describes sequences of states (usually changing in time) where the probability distribution of the next state only depends on the next states. Markov chains are ubiquitously used in science and technology for modeling sequences: text, DNA sequences, the weather, stock markets, etc. The changes in states are called *transitions*, and the changes going from one state to another are called the transition probability, which is stored in a *transition matrix*. For example, suppose we model the weather where there are only two states: "sunny" and "rainy". We can model weather transitions from one day to the next as:

$$T = \begin{bmatrix}0.9&0.1\\0.5 & 0.5\end{bmatrix}\,,$$
in which $T_{ij}$ represents the probability of state $i$ transitioning into state $j$. Here, it means that when it is sunny, we have a 90% chance that the next day will also be sunny and 10% that the next day will be rainy (first row). If we represent the state probability vector at time $t$ as $\mathbf{q}_t$, then we can update these probabilities as

$$\mathbf{q}_{t+1} = T^\intercal\mathbf{q}_t\,.$$
For example, if at $t=0$, it rains, i.e., $\mathbf{q}_0 = [0,1]^\intercal$ (we are 100% sure), then the next day, it has 50% of raining $\mathbf{q}_1 = [0.5,0.5]'$ and the next day, the state probabilities are $\mathbf{q}_2 = [0.5\times 0.9+0.5\times 0.5, 0.5\times 0.1+0.5\times0.5]^\intercal = [0.7, 0.3]^\intercal$. We see that our initial 'sharp' belief about the weather changes into a vector of state probabilities.

"""

# ╔═╡ e7bd6347-ac8a-4460-b0f7-4692968205f7
Tweather = [0.9 0.1;
			0.5 0.5]

# ╔═╡ 807fe053-bf07-446b-aac6-e8de716c5e2a
q0 = [0,1]

# ╔═╡ 012ff9ea-08ae-46fd-afab-49341f35f57a
q1 = Tweather' * q0  # distribution after one time step

# ╔═╡ cf858c95-d2b9-4f53-af1a-d5a511f531ce
q2 = Tweather' * q1  # distribution after two time steps

# ╔═╡ d77c1db1-925f-4134-a860-8a1683b3adee
md"""
We see that our initial 'sharp' belief about the weather changes into a vector of state probabilities. There are two ways to interpret these state probability vectors:

1. **population level**: In this interpretation, the state probability vector represents the fraction of the population or system occupying each state at a particular time step. When we consider a large number of random points in time (or far-away) locations, 50% of locations with rain will have rain the next day, and 50% will have sunny weather the next day.
2. **individual level**: Here, we consider a single individual (sometimes called a *particle* in terms of MCMC) over a long period of time and model the change in probability of being in different states. Here, a person who sees rain can compute the probability of seeing rain the next day and the day after that and so on.

 A natural question would be what the weather would likely be after a very long time, e.g., 100 days. We can compute this directly as:

$$\mathbf{q}_{100} = (T^\intercal)^{100}\mathbf{q}_0\,,$$

which evaluates as $\mathbf{q}_{100}\approx [5/6, 1/6]$.
"""

# ╔═╡ cc92d39d-830f-451c-81b4-de00f74545bf
q100 = Tweather'^100 * q0  # stationary distribution

# ╔═╡ c37ea17f-57b5-4871-a4d9-5aa5e609c846
md"""
Over long times, this chain seems to converge to a *stationary distribution*, often denoted as $\boldsymbol{\pi}$. In the population interpretation, this would mean that, from many regions for which the Markov chain can describe the weather transition, about 83% will have sunny weather in the long run. For an individual region, it means that if you wait long enough, there is an 83% chance of having a sun. Note that here, it would not matter with what initial state vector one starts with. We always converge to $\boldsymbol{\pi}$. For a stationary distribution, it should hold that

$$\boldsymbol{\pi} = T^\intercal \boldsymbol{\pi}\,,$$

meaning that the state probability vector is 'in equilibrium':the state probabilities do not change anymore. We can find the stationary distribution by normalizing the eigenvector of $T^\intercal$ that corresponds to an eigenvalue of 1. A Markov chain does not necessarily have a stationary distribution, for example, if it is periodic. 
"""

# ╔═╡ 38a6bd5f-78ab-43be-a010-6fe4c209b1a7
md"Consider a second example, in which we have ten states of which only transitions of $i$ to $i-1$, $i$, or $i+1$ are allowed (this is called a birth-death process). For example, the transition matrix is given by"

# ╔═╡ eab53e8a-b7cd-4c97-83f1-bcfd9c59cc4f
begin
	N = 10
	rng = MersenneTwister(11)
	
	p_remain = rand(rng, 0.1:.1:0.9, N)
	T = [i==j ? p_remain[i] : 0.0 for i in 1:N, j in 1:N]
	for i in 1:N-1
		T[i,i+1] = 1 - p_remain[i]
	end
	T[end,1] = 1 - p_remain[end]
	
	#T = Tridiagonal(0.3ones(9), 0.1ones(10), 0.6ones(9))
	#T ./= sum(T, dims=2)
end;

# ╔═╡ 7915a972-df48-462c-bcc7-2ca94a98e441
T

# ╔═╡ 53799772-c8a1-4055-a01e-42f147ffc253
md"Starting from a single state, we can model how the state vector evolves to a the stationary distribution."

# ╔═╡ b3c45c85-d075-452a-a827-b7d1e00a17a1
md"Start from state $(@bind i0 Select(1:10, default=4))"

# ╔═╡ fcfa5203-73df-40f6-a6c9-eda531054b38
p₀ = [i==i0 ? 1.0 : 0.0 for i in 1:size(T,2)]

# ╔═╡ b2cbf023-a0ea-4714-b30b-b331a081c7fa
π_mc = (T')^100 * p₀

# ╔═╡ 45e848aa-3e6a-47b5-a0cc-fa09dd9a4ae3
md"Note that, except for states 1 and 10, the probability of moving to a higher state is twice that of moving to a lower state. Remaining in the same state occurs with a probability of 0.1. When we start in the first state, we see that the state probability vector quickly settles to an equilibrium distribution where higher states are much more likely than lower states. The time needed to approximately reach the stationary distribution is called the *mixing time*."

# ╔═╡ ae080b48-87f8-4af6-b9d1-224080d4b787
md"""
*Ergodicity* is a crucial property for MCMC (Markov chain Monte Carlo) methods to function effectively. It means that a long chain will eventually visit all the states. In MCMC, we design Markov chains where the target distribution (the distribution we want to sample from) becomes the stationary distribution. This means the chain will eventually spend most of its time in regions with high probability under the target distribution.

Ergodicity ensures two key aspects for successful sampling:

1. **Recurrence**: The chain can revisit any state from any starting point after a finite number of steps. This guarantees the chain explores the entire sample space of the target distribution, preventing it from getting stuck in specific regions.
2. **Aperiodicity**: The chain doesn't get trapped in cycles of states. It can move freely between all possible states without getting stuck in a repeating pattern.

These properties allow MCMC to converge to the target distribution over time. Without ergodicity, the chain might get stuck in specific areas that don't represent the actual distribution, leading to biased samples. In essence, ergodicity ensures the chain effectively explores the entire landscape of the target distribution, enabling MCMC to generate representative samples.
"""

# ╔═╡ 0a4f1446-e49e-4281-be44-52d3ac4a7854
md"""
### Metropolis-Hastings

The *Metropolis-Hastings (MH) algorithm* works similarly to rejection sampling; we again use a proposal distribution from which we can easily sample. These samples can again be accepted or rejected based on information from the (unnormalized) target distribution. The main difference is that rejection sampling generates each sample independently, while the MH algorithm generates a path where each new candidate sample from the proposal distribution is based on the previously accepted sample. This way, the MH algorithm generates a random walk in the space one wants to sample. This random walk is constructed so that, in the long run, the probability of being in a particular region matches the target density or mass function.

We will denote the samples that are generated using indices $\mathbf{x}_1, \mathbf{x}_2, \mathbf{x}_3,\ldots$ and a candidate sample (which may or may not end up in the chain) is denoted as $\mathbf{z}^\star$. Our proposal distribution generates a new sample based on the previous sample at step $t$:

$$\mathbf{x}^\star \sim q(\mathbf{x}\mid \mathbf{x}_t)\,.$$

In the original formulation, the proposal distribution is symmetric (i.e., $q(\mathbf{x}\mid \mathbf{x}')=q(\mathbf{x}'\mid \mathbf{x})$), so this is simply referred to as the Metropolis algorithm. For example, as a proposal distribution, you could use a normal distribution centered around $\mathbf{x}_t$. When generating a candidate sample, it is accepted with an acceptance probability

$$\alpha(\mathbf{x}^\star, \mathbf{x}_t) = \min\left(1, \frac{\tilde{p}(\mathbf{x^\star})}{\tilde{p}(\mathbf{x}_t)}\right)\,.$$

In practice, one generates a random number $u$ over the unit interval $(0,1)$ and accepts the sample if $\alpha(\mathbf{x}^\star, \mathbf{x}_t)>u$ (note that you can ignore the minimum when the candidate has a higher probability than the previous sample, we always accept). So, if the candidate is accepted, we set $\mathbf{z}_{t+1}=\mathbf{z}^\star$; otherwise the sample is discarded and $\mathbf{z}_{t+1}=\mathbf{z}_t$ and another candidate is drawn from the same proposal distribution. In practice, however, only a single copy is kept, often with an integer weighting factor recording how often that state appears. Again, we emphasise that $\mathbf{x}_1, \mathbf{x}_2, \mathbf{x}_3,\ldots$ is not an independent sample, however, $\mathbf{x}_t$ tends to $p(\mathbf{x})$ when $t\rightarrow \infty$. Successive samples will be highly correlated. If you want independent samples, you can either:
- only retain every $M$-th sample in a long chain (e.g., every 100 samples);
- run several independent, shorter chains and take the final sample of each. 

It is not trivial to say which strategy is best in general. A chain will need a certain burn-in period to match the target distribution well, so a single, long chain might be effective. However, a set of independent chains could quickly explore different regions of our distribution, in addition to being able to run at several cores of your computer or cluster at the same time.

The reason why this algorithm works can be seen by seeing this as a continuous Markov chain with transition function $T(\mathbf{x}_{t+1}, \mathbf{x}_{t})$, for which a stationary distribution should satisfy:

$$p(\mathbf{x})=\sum_{\mathbf{x}'}T(\mathbf{x}', \mathbf{x})p(\mathbf{x}')\,.$$

A sufficient condition for ensuring that $p(\mathbf{x})$ is a stationary distribution if the transition probabilities satisfy *detailed balance* , defined by

$$p(\mathbf{x})T(\mathbf{x}, \mathbf{x}') = p(\mathbf{x}')T(\mathbf{x}', \mathbf{x})\,.$$

You might already see, that if we choose $T(\mathbf{x}, \mathbf{x}')$ to be symmetric, we can obtain the acceptance probability from choosing an transition probability

$$T(\mathbf{x}, \mathbf{x}') = \alpha(\mathbf{x}', \mathbf{x}) q(\mathbf{x}', \mathbf{x})\,.$$

From this, we can see that the suggested probability satisfies the detailed balance. 

Following this reasoning, the general MH algorithm with non-symmetric proposal distributions uses a slightly different acceptance probability:
$$\alpha(\mathbf{x}^\star, \mathbf{x}_t) = \min\left(1, \frac{\tilde{p}(\mathbf{x^\star})q(\mathbf{x}_t\mid \mathbf{x}^\star)}{\tilde{p}(\mathbf{x}_t)q(\mathbf{x}^\star\mid \mathbf{x}_t)}\right)\,,$$
The Metropolis algorithm is a special case. In Turing, you can specify the proposal distribution, though the default uses the prior distribution.
"""

# ╔═╡ b9587b83-1521-407a-a466-9566d022d402
function metropolis_hastings(p, x₀, q; n=100)
	samples = [x₀]
	accepted = [false]
	xₜ = x₀
	while sum(accepted) < n
		x′ = rand(q(xₜ))
		push!(samples, x′)
		α = min(1.0, p(x′) / p(xₜ))
		if rand() < α
			xₜ = x′
			push!(accepted, true)
		else
			push!(accepted, false)
		end
	end
	return samples, accepted
end		

# ╔═╡ 5d8894c1-25d7-41e5-a99c-bca8119851fe
q_MH = x -> Normal(x, 1)

# ╔═╡ 094b4405-6205-4b85-a311-c79cf6028498
md"Instead of using a homebrewn MH algorithm, it might now be the ideal time to illustrate the built-in Turing sampler on our `uncertain_normal` example. Sampling from the posterior distribution can be done using the `sample` function with the appropriate distribution, algorithm and sample size. Here, for didactive purposes, we fixed the seed in the optional first argument."

# ╔═╡ e1cb4a40-4a8a-4ff9-9406-f99aee6fdd57
md"By deafult, the MH algorithm uses the prior diistribution. You can set this to whichever distribution you like. We can use the function `summarize` to look at some summary statistics of our chains."

# ╔═╡ 56f7dfd2-2c81-403a-9c11-0d93c10f604e
md"For now, we focuss on the just the mean and standard deviation. The average and standard deviation of $\mu$ and $\sigma$ our chain more or less matches the what we expect based on the earlier plot. We can also use the function `quantile` to see these basic quantile distributions."

# ╔═╡ a8d1b7e7-03df-4d05-be5f-758f98bb1c26
md"Plotting the values of the chain, together with empirical data (can be done using just `plot(chain)`) shows that the MH chain does not generate a lote of unique values."

# ╔═╡ f97d0f4c-8b3c-4ea7-bcd2-6019ab831734
md" Most of the candidates are rejected! This is because our prior is far too broad as a candidate distribution. Let us create a more clever MH version."

# ╔═╡ 2cf4022e-1b68-4a69-b6d0-669c3e5ac230
my_MH = MH(:σ=>s->InverseGamma(s),
		   :μ=>m->Normal(m, 1/2))

# ╔═╡ 0c652efb-0edb-4a4e-a14c-7d5ac6eb5ffa
md"This looks much better! Our proposal is closer to our canidates, so we can stick closer to high-density regions. Likely using this chain, we can make sensible inferences about that parameters of the model! 

If we plot the chain on the posterior, we see a good fit."

# ╔═╡ 69f2c623-632f-4df4-8566-a48574c0f9fb
md"""
### Gibbs sampling

Gibbs sampling can be seen as a special case of the Metropolis-Hasting's algorithm. The idea behind Gibbs sampling is quite simple: though sampling from an unnormalized distribution $p(x_1, x_2, \ldots, x_n)$ is hard, it might be easier to sample each $x_i$ one at the time and conditioning on the previous samples of the other variables. For example, suppose we want to sample from a distribution with three variables, $p(x_1, x_2,x_3)$ and at step $t$, we have the sample $(x_{1,t}, x_{2,t}, x_{3,t}$). We cycle through the variables in turn, replacing $x_{1,t}$ by sampling

$$x_{1,t+1}\sim p(x_1\mid x_{2,t}, x_{3,t})\,.$$

Next, we sample the second variable as:

$$x_{2,t+1}\sim p(x_2\mid x_{1,t+1}, x_{3,t})\,,$$

and, finally:

$$x_{3,t+1}\sim p(x_2\mid x_{1,t+1}, x_{2,t+1})\,.$$

Then, we can cycle again through the variables to obtain a new sample $\mathbf{x}_{t+2}$ at infinitum. It is quite easy to see that this is a special case of the Metropolis-Hastings algorithm with for variable $x_k$ the proposal distribution $q_k(\mathbf{x}^\star\mid\mathbf{x}) = p(x_k\mid \mathbf{x}_{\setminus k})$ ($\mathbf{x}_{\setminus k}$ is the probability vector $\mathbf{x}$ without ${x}_{k}$) with as acceptance probability:

$$\alpha(\mathbf{x}^\star, \mathbf{x}_t) = \min\left(1, \frac{{p}(\mathbf{x^\star})q_k(\mathbf{x}_t\mid\mathbf{x}^\star)}{{p}(\mathbf{x}_t)q_k(\mathbf{x}^\star\mid\mathbf{x}_t)}\right)=\min\left(1, \frac{{p}({x}_k^\star\mid \mathbf{x}_{\setminus k}^\star)p(\mathbf{x}_{\setminus k}^\star)p(x_{k,t}\mid \mathbf{x}_{\setminus k}^\star)}{{p}({x}_{k,t}\mid \mathbf{x}_{\setminus k, t})p(\mathbf{x}_{\setminus k,t})p(x_k^\star\mid \mathbf{x}_{\setminus k, t})}\right)=1\,,$$
using $\mathbf{x}^\star_{\setminus k}=\mathbf{x}_{\setminus k,t}$. The Gibbs steps are hence always accepted.

In the types of models that we are building, Gibbs sampling is easy to use, as our distributions (Bayesian networks) are constructed from simple distributions with efficient sampling routines.

"""

# ╔═╡ 1a24999d-1ab7-4220-8e22-b1ac33ad03ac
md"As an example, consider sampling from a bivariate normal distribution, where we slice from one axis at a time."

# ╔═╡ 5e6c6135-f4f9-447c-8f7a-a49ec25e414a
md"""
σ₁: $(@bind σ₁ Slider(0.1:0.2:2, show_value=true, default=1))

σ₂: $(@bind σ₂ Slider(0.1:0.2:2, show_value=true, default=2))

ρ: $(@bind ρ Slider(-0.95:0.05:0.95, show_value=true, default=0.9))

n : $(@bind n_gibbs Slider(5:5:100, show_value=true, default=25))
"""

# ╔═╡ 3fca01f3-ac6e-4fdc-b8ca-48040ca562f3
μ = [1, 0]

# ╔═╡ d8ff625a-6204-4579-a8f5-9650544c7222
@model function uncertain_normal(y=missing)
	μ ~ Normal(0, 20)
	σ ~ InverseGamma(1/10)
	for i in 1:length(y)
		y[i] ~ Normal(μ, σ)
	end
end

# ╔═╡ b04ce3ab-c5c0-4a13-8d82-05198532b22d
uncertain_normal(y)

# ╔═╡ 7cfc4306-746b-434a-8160-52341c1e6519
norm_prior(μ, σ) = exp(logprior(uncertain_normal(y), (μ=μ, σ=σ)))

# ╔═╡ fb35033f-3d11-4ab1-bed8-2c22396064f8
norm_likelihood(μ, σ) = exp(loglikelihood(uncertain_normal(y), (μ=μ, σ=σ)))

# ╔═╡ e3e51306-385a-461b-94e8-e2a18b5250ab
norm_posterior(μ, σ) = norm_prior(μ, σ) * norm_likelihood(μ, σ)

# ╔═╡ 63af963c-a304-4142-8f23-a56a5fc14e76
chain_MH = sample(MersenneTwister(1), uncertain_normal(y), MH(), 1_000)

# ╔═╡ a9e11333-aab9-4e63-aee4-6bceb6b771bc
summarize(chain_MH)

# ╔═╡ 6309cc65-59d4-4df2-b756-30c6e2f877cd
quantile(chain_MH)

# ╔═╡ 1b92c00e-8b2a-4cc8-b6f2-23d91cf4f989
chain_MH2 = sample(MersenneTwister(1), uncertain_normal(y), my_MH, 1_000);

# ╔═╡ 37986301-c9d2-42e3-8b73-7c833bf7ae17
summarize(chain_MH2)

# ╔═╡ 138bf97c-9295-440d-99f4-24f7fe58b774
Σ = [σ₁^2 σ₁*σ₂*ρ;
	 σ₁*σ₂*ρ σ₂^2]

# ╔═╡ 87b984ee-2daa-4cca-b5eb-60c4aa3b2fb5
mvn = MultivariateNormal(μ, Σ)

# ╔═╡ 06e83a23-fa3a-49f7-8f43-d46f3285e943
chain_Gibbs = sample(uncertain_normal(y), Gibbs(MH(:μ), MH(:σ)), 1000)

# ╔═╡ a9bd2d76-16e6-4084-b65d-0645ec94aa9c
sample(uncertain_normal(y), Gibbs(HMC(0.2, 3, :μ), PG(20, :σ)), 1000)

# ╔═╡ 9b4d67dc-9c75-45ee-9328-d229c527fdc6
summarize(chain_Gibbs)

# ╔═╡ 1ba9efb2-cd95-4097-b8f6-3a770205c935
md"Gibbs sampling is especially powerful because it can recombine samplers for different variables."

# ╔═╡ e5d604e0-44c4-4d10-8faa-7ade5f5fbfff
md"""
### Hamiltonian Monte-Carlo

The previous MCMC methods could only take relatively small steps in the state space. In the worst case, they exhibit random walk behavior where the distance traveled through the state space grows with the square root of the number of steps (which we called slow in the previous chapter). They are also "blind" in the sense that they do not use a lot of information about the distribution. Hamiltonian Monte Carlo (HMC) methods are more sophisticated, as they transform the sampling problem into a dynamical system using ideas from physics. Solvers similar to those that are used to solve differential equations can then be used to generate the chain.

The main idea behind constructing the Markov chain is using concepts of classical mechanics to determine the trajectory of a particle, determined by the probability distribution. This particle has a state vector $\mathbf{x}$ (the position in state space) and a momentum vector $\mathbf{p}$, representing its velocity in state space. We interpret the (unnormalized) probability distribution as a potential function according to:

$$U(\mathbf{x}) = -\log \tilde{p}(\mathbf{x})\,,$$

where we see that regions with higher probability correspond to those with lower potential energy, the unknown normalization constant boils down to a constant shift of this potential and is, hence, of no importance. The kinetic energy of the particle is given by $\frac{1}{2}\mathbf{p}^\intercal\mathbf{p}$. One can combine the potential and kinetic energy into the *Hamiltonian function*:
$$H(\mathbf{x}, \mathbf{p}) = U(\mathbf{x}) +\frac{1}{2} \mathbf{p}^T\mathbf{p}\,.$$

This Hamiltonian completely determines the dynamics of the particle according to the relations:

$$\frac{\partial x_i}{\partial t} = \frac{\partial H}{\partial p_i}$$

and

$$\frac{\partial p_i}{\partial t} = -\frac{\partial H}{\partial x_i}\,.$$

It is easy to see that these equations of motion preserve the Hamiltonian ($\frac{\mathrm{d}H}{\mathrm{d}t}=0$). An essential property of Hamiltonian dynamical systems is that they preserve volume in state space. This is known as *Liouville's theorem*., meaning that if you consider a region given the space of variables $(\mathbf{x},\mathbf{p})$, the shape of this region will change, but the volume will remain the same. 

Steps in the state space can be taken by numerical integration, for which Liouville's theorem still holds. One set of integration schemes is the *leapfrog integration*, given by

$$\mathbf{p}(t+\epsilon/2) = \mathbf{p}(t) - \frac{\epsilon}{2}\nabla U (\mathbf{x})$$

$$\mathbf{x}(t+\epsilon) = \mathbf{x}(t) + \epsilon\, \mathbf{p}(t-\epsilon /2)$$

$$\mathbf{p}(t+\epsilon) = \mathbf{p}(t+\epsilon/2) - \frac{\epsilon}{2}\nabla U (\mathbf{x})$$

where $\epsilon$ is the step size. Note that we have here used functional dependence of $t$ to distinguish continuous integration steps with discrete steps in the Markov chain ($\mathbf{x}(t)$ vs $\mathbf{x}_t$). One takes several steps, denoted by integer $L$ leapfrog, to generate a new sample in the chain. Afterward, one usually resamples the momentum vector $\mathbf{p}$ to obtain an ergodic sampling scheme. To recap, the steps of the HMC are

1. sample an initial state vector $\mathbf{x}_0$ from the prior and generate a random momentum vector $\mathbf{p}_0$;
2. repeat $T$ times:
   - perform $L$ leapfrog steps using step size $\epsilon$ on $(\mathbf{x}_t,\mathbf{p}_t)$ to obtain $(\mathbf{x}_{t+1},\mathbf{p}_{t+1})$
   - resample $\mathbf{p}_{t+1}$

"""

# ╔═╡ 4019f495-ce19-44e4-ad79-18349fa6ec18
chain_HMC = sample(MersenneTwister(1), uncertain_normal(y), HMC(0.1, 5), 1000);

# ╔═╡ afbf92d5-c72e-4041-a265-50d762df58f7
summarize(chain_HMC)

# ╔═╡ bec15afc-3f02-4c7d-b9e2-01f496e31689
summarize(chain_HMC)

# ╔═╡ 786c8319-806c-461d-af08-f082d6f12339
md"You may have noticed that HMC is more computationally demanding than the previously seen methods. We need to perform several integration steps to generate a single pseudo-sample (for which we need to keep track of momentum variables in addition to state variables) and compute the gradient of the log-probability density function. This gradient is usually computed automatically using automatic differentiation, which generally does not require too much computational resources. It also implies that the basic version of HMC can only be used for real-valued random variables, as one has to be able to compute derivatives. Given the same computational budget, one will obtain much fewer pseudo-samples using HMC than with MH. However, the samples generated by HMC will often be of a higher quality and cover the state space much more thoroughly, so it is likely worth the trade-off.
"

# ╔═╡ 81863647-168f-46b1-87f1-cc1169b7c6ff
md"""
Making HMC behave well requires carefully tuning the parameters $\epsilon$ and $L$. For this reason, the [No-U-Turn Sampler](https://arxiv.org/pdf/1111.4246.pdf) (NUTS). This automatically detects when sufficient leapfrog steps are taken by randomly moving forward and backward in time until it detects a "U-turn", i.e., it backtracks in the state space.
"""

# ╔═╡ aa022156-efb4-4b2b-ac22-9b4228ee2578
chain_NUTS = sample(MersenneTwister(1), uncertain_normal(y), NUTS(), 1000);

# ╔═╡ 5f41a9e3-70e0-40d5-bdff-d90555130817
summarize(chain_NUTS)

# ╔═╡ e6201f93-d5fc-4f33-9e3c-f7d8eb0146b5
md"""
## MCMC in practice

Sampling from probability distributions is more complex than numerically solving differential equations. For both, using existing, tested, high-quality software is strongly advised. Picking the right sampling algorithm and making it work well requires insight into the model, choosing hood proposal distributions, and some experience.

Rather than looking at the number of pseudo samples generated, it usually makes more sense to consider the [[effective number of samples]]. This metric is corrected for autocorrelation, i.e., samples close in a chain are not completely independent. There is no definitive answer to how many samples one needs. If the goal is pinpointing the posterior mean, a couple hundred could suffice. When one wants to characterize the exact shape of the posterior or analyse the tails, for example, 1% or 99% quantiles, many more might be required.

Every chain needs some warm-up time to reach the stationary distribution. For this reason, the first fraction of the chain is usually discarded. One often runs several Markov chains in parallel, often on different computer threads or different nodes on a computer cluster. This might be valuable for diagnostic purposes, to see if the chains converge. However, note that every one of these chains will likely need a burn-in time, so a large part of the computational efforts will be wasted. Some authors advise trying several short chains for debugging and a very long chain for the final sampling.

To check if a chain is converging nicely, trace plots are usually the most informative. Remember, a trace plot shows the value of the variable throughout the steps of the chain. An ideal trace plot should look somewhat like a hairy caterpillar, with fluctuations around a mean and no systematic trends. When the chain is behaving badly, expect to see sharp peaks, flat lines where the chain is stuck and systematic trends.

> When you have a computational problem, there is often a problem with your model.

The best way to have a chain that works well is to have a model that describes the data well. Here, there is ideally a single peak that corresponds to the optimal parameter configurations. When you see that the chain is not progressing well, this can usually be improved by adding some informative priors, even if they are very weak ones. This will often tamper with erratic behaviour.
"""

# ╔═╡ 1b2a38dc-0caf-4c86-9910-8574fe10c48f
@model function diffuse_prior(y1, y2)
	μ ~ Turing.Flat()
	σ ~ Turing.FlatPos(0.0)
	y1 ~ Normal(μ, σ)
	y2 ~ Normal(μ, σ)
end

# ╔═╡ 9e1c4bc8-882f-45bf-b437-38c369867dac
y1, y2 = 7.8, 9.8

# ╔═╡ 15c7c5ed-5bcf-4336-aecb-22fd3ef38d05
chain_diff = sample(diffuse_prior(y1, y2), NUTS(), 10_000)

# ╔═╡ cde76a39-51a0-4398-ba2c-81170e395a07
@model function weak_prior(y1, y2)
	μ ~ Normal(0, 1000)
	σ ~ Exponential(10) #Uniform(0, 100)
	y1 ~ Normal(μ, σ)
	y2 ~ Normal(μ, σ)
end

# ╔═╡ 42a81a70-f8a1-4a11-a475-c4adddefde73
chain_weak = sample(weak_prior(y1, y2), NUTS(), 10_000)

# ╔═╡ 930902f1-5d0d-41a9-8c37-07e0d913aa6e
quantile(chain_weak)

# ╔═╡ 96f05537-6fb3-45ce-8d75-28efdfa03279
@model function donut(R=5, σ=1.5)
	θ ~ Uniform(0.0, 2pi)
	x ~ Normal(R * cos(θ), σ)
	y ~ Normal(R * sin(θ), σ)
end

# ╔═╡ 3708eb56-53e4-43b8-8af4-e1f6e302a0cf
let
	samples_donut = [rand(donut()) for i in 1:1000]
	x = [s[:x] for s in samples_donut]
	y = [s[:y] for s in samples_donut]
	scatter(x, y, aspect_ratio=:equal, label="", xlab=L"x", ylab=L"y", alpha=0.75)
end

# ╔═╡ 07da284a-8f21-4088-b936-244d11f517dd
x_nuts = sample(donut(), HMC(0.1, 10), 10000)

# ╔═╡ 9d14bb2f-2c8e-4b45-9775-e168a11005ff
begin
	pdf_donut(x, y) = exp(())
	#heatmap(-15:0.1:15, 0.1:0.02:10, pdf_donut, color=:speed)
	scatter(x_nuts[:x], x_nuts[:y], aspect_ratio=:equal, label="", xlab=L"x", ylab=L"y", alpha=0.25)
end

# ╔═╡ 5ce17f7f-7e93-4c43-9f75-b3b677af9fd1
@model function trending(yl, ym, yr)
	σa ~ Gamma()
	σn ~ Gamma()
	a ~ Normal(0, σa)
	b ~ Normal(0, σa)
	yl ~ Normal(-a + b, σn)
	ym ~ Normal(b, σn)
	yr ~ Normal(a + b, σn)
end

# ╔═╡ 71c210f8-34ad-469a-89fa-dc99efb6b72a
quantile(sample(trending(-2.3, 0.3, -2.8), NUTS(), 10_000))

# ╔═╡ d5db3698-8546-4606-8346-85475619a3e0
quantile

# ╔═╡ 6a178a79-ae85-4496-abef-11bf62d64da8
md"# Appendix 🐉"

# ╔═╡ d9d417f3-a427-409e-9090-c76c8226be04
TableOfContents()

# ╔═╡ 0d45d6f9-f27c-4d8e-b117-26132418be55
dist2pdf(distribution) = x -> pdf(distribution, x)

# ╔═╡ 037b15f6-0858-41fa-a8fd-c6cc8c8d8e7a
function seed_bayesian_plot(k; a=8, b=3, n=10, title="Inference of p given $k of $n seed germinated", labels=true)

	seed_likelihood(p) = pdf(Binomial(n, p), k)
	seed_prior(p) = pdf(Beta(a, b), p)
	seed_posterior(p) = pdf(Beta(a+k, b+n-k), p)

	p_ML = argmax(seed_likelihood, 0:0.01:1)
	p_MAP = argmax(seed_posterior, 0:0.01:1)
		
	plot_seed = plot(seed_prior, 0, 1, lw=2, label=labels ? "prior (Beta($a, $b))" : "", xlab=L"p", color="green"; title)
	plot!(p->n * seed_likelihood(p), 0, 1, lw=2, label=labels ? "likelihood (rescaled)" : "", color="blue", ls=:dash)
	plot!(seed_posterior, 0, 1, lw=2, label=labels ? "posterior" : "", color="orange", ls=:dashdot)
	vline!([p_ML], lw=2, alpha=0.5, label= labels ? "p* (maximum likelihood)" : "", ls=:dash)
	vline!([p_MAP], lw=2, alpha=0.5, label=labels ? "p* (maximum posterior likelihood)" : "", ls=:dot)
	return plot_seed
end

# ╔═╡ df004387-ea5f-42b1-9409-7b37d1a5031a
md"Distributions for rejection sampling"

# ╔═╡ 01fb5888-a01a-4dd8-af08-bb70894a99b8
p_univar = MixtureModel([Normal(2, 0.4), Normal(4, 1.2)], [0.25, 0.75])

# ╔═╡ 66c3e374-389c-4be9-bf16-341565d4b4c9
p = dist2pdf(p_univar)

# ╔═╡ 730e5cbb-1787-44b3-9bb4-32dba7035ea8
@model function seed_germ(X=missing, n=10)
	p ~ Beta(8, 3)
	X ~ Binomial(n, p)
end

# ╔═╡ 91cfb385-5617-4129-b57b-69aff3621ce5
x_mh, acc_mh = metropolis_hastings(p, 4.0, q_MH, n=200)

# ╔═╡ f578aafc-318f-4dc4-a574-0843057b1204
q = Truncated(Normal(2.5, 3), -1, 10)

# ╔═╡ efefaff1-e6c7-46fa-bfa7-b4e6dec67641
plots = Dict()

# ╔═╡ cc904bc9-7f2c-4943-90ff-e0753aea08b8
let
	prior_informative = Normal(5, 0.5) |> dist2pdf
	prior_weakly = TriangularDist(1, 8, 5) |> dist2pdf
	prior_diffuse = Uniform(0, 10) |> dist2pdf

	p = plot(prior_informative, 0, 12, lw=2, label="informative prior", xlab=L"\theta", ylab=L"f_\theta(\theta)")
	plot!(prior_weakly, 0, 12, label="weakly informative prior", lw=2, ls=:auto)
	plot!(prior_diffuse, 0, 12, label="diffuse informative prior", lw=2, ls=:auto)
	title!("Different types of priors")
	plots["priors"] = p
end

# ╔═╡ 0bc4cdff-f86d-4174-8aac-d2ca61717f96
plots["seeds_likelihood"] = plot(seed_likelihood, 0, 1, lw=2, label="likelihood", xlab=L"p", color="blue", ls=:dash)

# ╔═╡ 08d33254-5013-4e4f-ae0b-1f5b6f3a2fb2
plots["seeds_prior"] = plot(seed_prior, 0, 1, lw=2, label="prior", xlab=L"p", color="green")

# ╔═╡ 217edb24-a555-4deb-ab58-2969d8564a7b
plots["seeds_posterior"] = plot(seed_posterior, 0, 1, lw=2, label="posterior", xlab=L"p", color="orange", ls=:dashdot)

# ╔═╡ 34a5bb1a-8919-4e41-89a5-fdd4e48bec26
plots["seeds_bayes"] = seed_bayesian_plot(k)

# ╔═╡ 92cd9bfb-8a58-4b73-ab02-5ac4447499d0
plots["seeds_10"] = seed_bayesian_plot(10, labels=false)

# ╔═╡ fb486bb6-b462-40f5-8189-2b46c5817fc1
plots["seeds_2"] = seed_bayesian_plot(2, labels=false)

# ╔═╡ 780b3e4e-7959-401f-aca5-5dad90276128
plots["seeds_weak_prior"] = seed_bayesian_plot(k, a=1, b=1, labels=false)

# ╔═╡ ba15b522-e4e8-4de7-882e-a0d14da2fac1
plots["seeds_strong_prior"] = seed_bayesian_plot(k, a=80, b=30, labels=false)

# ╔═╡ 9609ee3c-2a5e-450f-8822-d616552cae09
plots["seeds_large_dataset"] = seed_bayesian_plot(10k, n=100, labels=false)

# ╔═╡ 3eac18d0-f755-42aa-bbcc-9714cc64dc5a
plots["norm_muprior"] = plot(dist2pdf(Normal(0, 20)), -50, 50, xlabel=L"\mu", label="Normal(0, 20)", title="Prior for μ", lw=2)

# ╔═╡ cd8bf990-4f36-4a56-83c2-c58217448b1e
plots["norm_sigmaprior"] = plot(dist2pdf(InverseGamma(.1)), 0, 20, xlabel=L"\sigma", label="InverseGamma(1)", title="Prior for σ", lw=2)

# ╔═╡ 45e175bc-0eaa-4280-aaf6-c8b318a6c12d
plots["norm_priors"] = plot(plots["norm_muprior"], plots["norm_sigmaprior"],size=(800, 400))

# ╔═╡ a09eda2e-0bb7-46c9-8d60-2b4447257772
plots["norm_multiprior"] = heatmap(-15:0.1:15, 0.1:0.02:5, norm_prior, color=:speed, ylab=L"\sigma", xlab=L"\mu", title="Prior")

# ╔═╡ 3b821e04-a27a-489f-924e-e82a08693371
plots["norm_likelihood"] = heatmap(-15:0.1:15, 0.1:0.02:5, norm_likelihood, color=:speed, ylab=L"\sigma", xlab=L"\mu", title="Likelihood")

# ╔═╡ b4eef7a8-a61c-471d-baf0-1d0b5fb997e6
plots["norm_posterior"] = plot_post = heatmap(-15:0.1:15, 0.1:0.02:5, norm_posterior, color=:speed, ylab=L"\sigma", xlab=L"\mu", title="Posterior")

# ╔═╡ 2c81f5d1-2b8b-4fb5-b810-d5fcfdc8a5d6
plots["norm_posterior"] = plot_post;

# ╔═╡ fae76f2c-caae-49ff-9dfe-b67cd53d8dea
begin
    x_prop = rand(q, n_rejection_sampling)
	α_prop = pdf.(Ref(p_univar), x_prop) ./ (M .* pdf.(Ref(q), x_prop))
	u_rj = rand(n_rejection_sampling) 
	acc = u_rj .≤ α_prop
	u_rj .*= M .* pdf.(Ref(q), x_prop)

	x_rejection_sampling = x_prop[acc]
	acceptance_probability = length(x_rejection_sampling) / n_rejection_sampling
	
	prs = plot(x->pdf(p_univar, x), -1, 10, lw=2, label="target p(x)", xlabel=L"x", title="Rejection sampling\np(accept)=$(round(acceptance_probability, digits=2))",
	legend=:outerbottom)
	plot!(x->M*pdf(q, x), -1, 10, lw=2, label="proposal distribiton M * q(x)")
	plot!(zero, -1:0.02:10, fillrange=x->pdf(p_univar, x), fillalpha=0.3, lw=2, label="acceptance region", linealpha=0)
	plot!(x->pdf(p_univar, x), -1, 10, fillrange=x->M*pdf(q, x), fillalpha=0.3, lw=2, label="rejection region", linealpha=0)
	scatter!(x_prop[acc], u_rj[acc], label="accepted", ms=3, alpha=0.8)
	scatter!(x_prop[.!acc], u_rj[.!acc], label="rejected", ms=3, alpha=0.8)
	plots["rejection_sampling"] = prs
end

# ╔═╡ 9db5d774-8e6d-4bd1-a85a-d0e0705197ae
x_rejection_sampling

# ╔═╡ 525c0e2e-6353-4002-a486-bbdde0911c38
length(x_rejection_sampling)  # number of accepted samples

# ╔═╡ d3b6cf36-a2d1-40fa-952c-f20d1468eb38
acceptance_probability

# ╔═╡ ff813a28-5ab7-4b49-8b05-b4b6a73b94a5
plots["MC_ss"] = bar(1:10, π_mc, xlab="state", ylab="PMF", label="", title="Stationary distribution of the cycle", xticks=1:10)

# ╔═╡ bfd58734-b396-4424-9343-4562932f70e3
plots["MC_ss_conv"] = groupedbar(0:30, vcat([p₀' * T^i for i in 0:30]...), bar_position = :stack, xlab = L"t", ylab = "fraction in state i", xticks=0:30, label=reshape(["state $i" for i in 1:size(T,2)], 1, :), title="State evolution of the cycle", legend=:outertopright)

# ╔═╡ ecfc2307-6926-42c1-b08d-6ff93502c5cd
let
	mc = ifelse.(acc_mh, :blue, :orange)
	ind = 0:length(x_mh)-1
	
	p = plot(x->pdf(p_univar, x), -1, 10, lw=2, label="", xlabel=L"x", title="Metropolis-Hastings", ylab="PDF")
	
	plot!(twinx(), x_mh, ind, lw=0.5, label="", ylab="sample number", color=:green, ylim=(0, 300))
	scatter!(twinx(), x_mh[acc_mh], ind[acc_mh], label="accepted", color=:blue, ms=2, ylim=(0, 300))
	scatter!(twinx(), x_mh[.!acc_mh], ind[.!acc_mh], label="rejected", color=:orange, ms=2, ylim=(0, 300))
	plots["MH"] = p
end

# ╔═╡ 11fc7beb-2a11-4acb-8545-9842ca783682
let
	p = histogram(x_mh[acc_mh], normalize=true, label="pseudosamples MH", xlab=L"x")
	plot!(p, 0, 10, label="PDF target", lw=2, title="Histogram Metropolis-Hastings")
	plots["MH_hist"] = p
end

# ╔═╡ 6227e4ee-6f89-4e45-b1d6-a087aa1bad40
plots["MH_diagn"] = plot(chain_MH)

# ╔═╡ 35bc2b8f-716d-4a85-8e60-7242709448c3
plots["MH_diagn_tuned"] = plot(chain_MH2)

# ╔═╡ ab09925c-4689-4f87-be44-f8eca12a5c9c
let
	p = contour(-15:0.1:15, 0.1:0.02:5, norm_posterior, color=:speed, ylab=L"\sigma", xlab=L"\mu", title="Metropolis-Hasting sampling", xlim=(-15,15), ylim=(0.1, 5))
	plot!(chain_MH2[:μ], chain_MH2[:σ], label="", alpha=0.7)
	scatter!(chain_MH2[:μ], chain_MH2[:σ], ms=0.2, label="chain")
	plots["MH_sampling"] =  p
end

# ╔═╡ 76b0886a-fc91-40ef-89d5-dcf5d261907f
let	
	X1 = Normal(μ[1], σ₁)
	X1cond(x2) = Normal(μ[1] + Σ[1,2]*inv(Σ[2,2])*(x2-μ[2]), 
						√(Σ[1,1] - Σ[2,1]*inv(Σ[2,2])*Σ[1,2]))
	X2cond(x1) = Normal(μ[2] + Σ[1,2]*inv(Σ[1,1])*(x1-μ[1]), 
					√(Σ[2,2] - Σ[1,2]*inv(Σ[1,1])*Σ[2,1]))

	x1 = rand(X1)
	x2 = rand(X2cond(x1))
	(x1, x2)

	x1_gibbs = [x1]
	x2_gibbs = [x2]
	for t in 1:n_gibbs
		x2 = last(x2_gibbs)
		x1 = rand(X1cond(x2))
		push!(x1_gibbs, x1)
		push!(x2_gibbs, x2)
		x2 = rand(X2cond(x1))
		push!(x1_gibbs, x1)
		push!(x2_gibbs, x2)
	end
	p = contour(-5:0.05:5, -5:0.05:5, (x,y)->pdf(mvn, [x,y]), color=:speed, xlab=L"x", ylab=L"y")
	plot!(x1_gibbs, x2_gibbs, color="blue", alpha=0.7, label="Gibbs chain", lw=1.3)
	scatter!(x1_gibbs[2:2:end], x2_gibbs[2:2:end], label="sample", ms=2, alpha=0.7)
	#scatter!([last(x1_gibbs)], [last(x2_gibbs)])
	#plot!(x->pdf(X1cond(last(x2_gibbs)), x)-5, -4, 4, label="X1|x2")
	plots["Gibbs"] = p
	
end

# ╔═╡ 444a89d1-a8c1-4df8-945d-6e8dd32ae204
plots["Gibbs_diagn"] = plot(chain_Gibbs)

# ╔═╡ 4d8fd482-11cb-4f2e-8c47-2428ee7e66c6
let
	p = contour(-15:0.1:15, 0.1:0.02:5, norm_posterior, color=:speed, ylab=L"\sigma", xlab=L"\mu", title="Gibbs sampling", xlim=(-15,15), ylim=(0.1, 5))
	plot!(chain_Gibbs[:μ], chain_Gibbs[:σ], label="", alpha=0.7)
	scatter!(chain_Gibbs[:μ], chain_Gibbs[:σ], ms=0.2, label="chain")
	plots["Gibbs_sampling"] = p
end

# ╔═╡ 99d638e6-1fae-4020-af2f-28703cc67377
plots["HMC_diagn"] = plot(chain_HMC)

# ╔═╡ 7f729aa1-78c0-4d2e-9ae5-c27502be22e8
let
	p = contour(-15:0.1:15, 0.1:0.02:5, norm_posterior, color=:speed, ylab=L"\sigma", xlab=L"\mu", title="HMC sampling")
	plot!(chain_HMC[:μ], chain_HMC[:σ], label="", alpha=0.7)
	scatter!(chain_HMC[:μ], chain_HMC[:σ], ms=0.2, label="chain")
	plots["HMC_sampling"] = p
end

# ╔═╡ 575ad310-8cbd-4b05-aee0-d84fb55ad6ed
plots["NUTS_diagn"] = plot(chain_NUTS)

# ╔═╡ 4a86bf84-434a-4ef3-a9c6-4c72df5068e0
let
	p = contour(-15:0.1:15, 0.1:0.02:5, norm_posterior, color=:speed, ylab=L"\sigma", xlab=L"\mu", title="NUTS sampling")
	plot!(chain_NUTS[:μ], chain_NUTS[:σ], label="", alpha=0.7)
	scatter!(chain_NUTS[:μ], chain_NUTS[:σ], ms=0.2, label="chain")
	plots["NUTS_sampling"] = p
end

# ╔═╡ 8de88cd0-e326-4926-bfe9-776e3516d6e2
plots["uninformative_prior"] = plot(chain_diff)

# ╔═╡ 0969b685-e9c5-44d5-8698-a141dd01dca7
plots["weak_diffusive_prior"] = plot(chain_weak)

# ╔═╡ 262879f7-9441-4305-899a-b7c2881958fe
plots

# ╔═╡ d2414490-d2a5-4862-81e8-8ae8ffc5213a
length(plots)

# ╔═╡ Cell order:
# ╠═103e5ba0-cfdc-11ee-13b1-cf53dfdd9a3b
# ╠═d778ce04-8df4-42ef-95f1-9cf4880e0420
# ╠═0350aa5b-d105-4dfa-a454-59873672b3a0
# ╟─e87dd9a4-f4ed-46d7-9f9d-dae0cadb7d05
# ╟─0f6cf892-b218-4740-a298-43f47e51acae
# ╟─cc904bc9-7f2c-4943-90ff-e0753aea08b8
# ╟─eb37dd4b-3af0-4534-92c2-e495b83024be
# ╠═29928abf-f010-438e-9b60-e6673a51ddf5
# ╠═730e5cbb-1787-44b3-9bb4-32dba7035ea8
# ╠═768ea152-00a5-4409-86eb-9a1c554eb3b2
# ╠═959cd079-59a4-42f5-a868-4cc0675c694b
# ╠═39cb8a88-9bfe-4b51-8f2f-b89b45dccc95
# ╠═667f2069-fd32-4c3e-aebb-b75a38bf84b9
# ╟─472cea3d-12a8-449a-8e87-6c1b5aaa1fd7
# ╟─0bc4cdff-f86d-4174-8aac-d2ca61717f96
# ╠═08d33254-5013-4e4f-ae0b-1f5b6f3a2fb2
# ╠═217edb24-a555-4deb-ab58-2969d8564a7b
# ╠═3ef47022-f9e3-434e-87a3-20b1870a7c24
# ╠═6004b26a-8f10-492d-bdd2-1d317d8bb394
# ╠═34a5bb1a-8919-4e41-89a5-fdd4e48bec26
# ╠═92cd9bfb-8a58-4b73-ab02-5ac4447499d0
# ╠═fb486bb6-b462-40f5-8189-2b46c5817fc1
# ╠═780b3e4e-7959-401f-aca5-5dad90276128
# ╠═ba15b522-e4e8-4de7-882e-a0d14da2fac1
# ╠═9609ee3c-2a5e-450f-8822-d616552cae09
# ╟─2ab01769-42cd-44e1-8d4c-2c8b477c732a
# ╠═4d6d6c9e-c2a8-4c7d-9d18-89d2e9d0c213
# ╠═ba06e8a5-c142-4963-97b4-783ed4c706bf
# ╟─f44abcc5-9010-4d10-b722-ef8147c8fd66
# ╠═3eac18d0-f755-42aa-bbcc-9714cc64dc5a
# ╠═cd8bf990-4f36-4a56-83c2-c58217448b1e
# ╟─45e175bc-0eaa-4280-aaf6-c8b318a6c12d
# ╟─0c29af96-949d-47b3-869f-274ed39c9d25
# ╠═d8ff625a-6204-4579-a8f5-9650544c7222
# ╟─2cbe07e6-8cf0-4dc4-8aea-f56059ffa367
# ╠═b04ce3ab-c5c0-4a13-8d82-05198532b22d
# ╟─afb7bdca-a8de-4342-9212-5e81fa47963a
# ╠═7cfc4306-746b-434a-8160-52341c1e6519
# ╠═a09eda2e-0bb7-46c9-8d60-2b4447257772
# ╟─1fbd13c1-cc93-4e5b-98f1-ace195d20f97
# ╠═fb35033f-3d11-4ab1-bed8-2c22396064f8
# ╠═3b821e04-a27a-489f-924e-e82a08693371
# ╟─213465af-f980-4848-9e5a-0cb5429f4af8
# ╠═e3e51306-385a-461b-94e8-e2a18b5250ab
# ╟─b4eef7a8-a61c-471d-baf0-1d0b5fb997e6
# ╟─2c81f5d1-2b8b-4fb5-b810-d5fcfdc8a5d6
# ╟─7f2a7d88-442f-431e-8e09-42082bda8d24
# ╟─985c1d3e-bc02-4dc5-9a28-d0ec9edf5cc9
# ╟─a0a1deec-75b6-4061-aa32-ad5e940d4123
# ╟─d91cea07-fb14-4b93-8ad1-31fb8f3afe8d
# ╟─fae76f2c-caae-49ff-9dfe-b67cd53d8dea
# ╠═9db5d774-8e6d-4bd1-a85a-d0e0705197ae
# ╠═525c0e2e-6353-4002-a486-bbdde0911c38
# ╠═d3b6cf36-a2d1-40fa-952c-f20d1468eb38
# ╟─d669d181-dc13-4921-9aca-88acae1197eb
# ╟─9c48f6a4-db8c-4e15-8486-4dd9b281a80d
# ╠═e7bd6347-ac8a-4460-b0f7-4692968205f7
# ╠═807fe053-bf07-446b-aac6-e8de716c5e2a
# ╠═012ff9ea-08ae-46fd-afab-49341f35f57a
# ╠═cf858c95-d2b9-4f53-af1a-d5a511f531ce
# ╟─d77c1db1-925f-4134-a860-8a1683b3adee
# ╠═cc92d39d-830f-451c-81b4-de00f74545bf
# ╟─c37ea17f-57b5-4871-a4d9-5aa5e609c846
# ╟─38a6bd5f-78ab-43be-a010-6fe4c209b1a7
# ╟─eab53e8a-b7cd-4c97-83f1-bcfd9c59cc4f
# ╠═7915a972-df48-462c-bcc7-2ca94a98e441
# ╟─53799772-c8a1-4055-a01e-42f147ffc253
# ╟─fcfa5203-73df-40f6-a6c9-eda531054b38
# ╠═b2cbf023-a0ea-4714-b30b-b331a081c7fa
# ╟─ff813a28-5ab7-4b49-8b05-b4b6a73b94a5
# ╟─b3c45c85-d075-452a-a827-b7d1e00a17a1
# ╠═bfd58734-b396-4424-9343-4562932f70e3
# ╟─45e848aa-3e6a-47b5-a0cc-fa09dd9a4ae3
# ╟─ae080b48-87f8-4af6-b9d1-224080d4b787
# ╟─0a4f1446-e49e-4281-be44-52d3ac4a7854
# ╠═b9587b83-1521-407a-a466-9566d022d402
# ╠═5d8894c1-25d7-41e5-a99c-bca8119851fe
# ╠═91cfb385-5617-4129-b57b-69aff3621ce5
# ╟─ecfc2307-6926-42c1-b08d-6ff93502c5cd
# ╠═11fc7beb-2a11-4acb-8545-9842ca783682
# ╟─094b4405-6205-4b85-a311-c79cf6028498
# ╠═63af963c-a304-4142-8f23-a56a5fc14e76
# ╟─e1cb4a40-4a8a-4ff9-9406-f99aee6fdd57
# ╠═a9e11333-aab9-4e63-aee4-6bceb6b771bc
# ╟─56f7dfd2-2c81-403a-9c11-0d93c10f604e
# ╠═6309cc65-59d4-4df2-b756-30c6e2f877cd
# ╟─a8d1b7e7-03df-4d05-be5f-758f98bb1c26
# ╠═6227e4ee-6f89-4e45-b1d6-a087aa1bad40
# ╟─f97d0f4c-8b3c-4ea7-bcd2-6019ab831734
# ╠═2cf4022e-1b68-4a69-b6d0-669c3e5ac230
# ╠═1b92c00e-8b2a-4cc8-b6f2-23d91cf4f989
# ╠═37986301-c9d2-42e3-8b73-7c833bf7ae17
# ╠═35bc2b8f-716d-4a85-8e60-7242709448c3
# ╟─0c652efb-0edb-4a4e-a14c-7d5ac6eb5ffa
# ╠═ab09925c-4689-4f87-be44-f8eca12a5c9c
# ╟─69f2c623-632f-4df4-8566-a48574c0f9fb
# ╟─1a24999d-1ab7-4220-8e22-b1ac33ad03ac
# ╟─5e6c6135-f4f9-447c-8f7a-a49ec25e414a
# ╠═3fca01f3-ac6e-4fdc-b8ca-48040ca562f3
# ╟─138bf97c-9295-440d-99f4-24f7fe58b774
# ╠═87b984ee-2daa-4cca-b5eb-60c4aa3b2fb5
# ╟─76b0886a-fc91-40ef-89d5-dcf5d261907f
# ╠═06e83a23-fa3a-49f7-8f43-d46f3285e943
# ╠═a9bd2d76-16e6-4084-b65d-0645ec94aa9c
# ╠═444a89d1-a8c1-4df8-945d-6e8dd32ae204
# ╠═9b4d67dc-9c75-45ee-9328-d229c527fdc6
# ╠═4d8fd482-11cb-4f2e-8c47-2428ee7e66c6
# ╟─1ba9efb2-cd95-4097-b8f6-3a770205c935
# ╟─e5d604e0-44c4-4d10-8faa-7ade5f5fbfff
# ╠═4019f495-ce19-44e4-ad79-18349fa6ec18
# ╠═afbf92d5-c72e-4041-a265-50d762df58f7
# ╠═99d638e6-1fae-4020-af2f-28703cc67377
# ╠═bec15afc-3f02-4c7d-b9e2-01f496e31689
# ╠═7f729aa1-78c0-4d2e-9ae5-c27502be22e8
# ╟─786c8319-806c-461d-af08-f082d6f12339
# ╟─81863647-168f-46b1-87f1-cc1169b7c6ff
# ╠═aa022156-efb4-4b2b-ac22-9b4228ee2578
# ╠═5f41a9e3-70e0-40d5-bdff-d90555130817
# ╠═575ad310-8cbd-4b05-aee0-d84fb55ad6ed
# ╟─4a86bf84-434a-4ef3-a9c6-4c72df5068e0
# ╟─e6201f93-d5fc-4f33-9e3c-f7d8eb0146b5
# ╠═1b2a38dc-0caf-4c86-9910-8574fe10c48f
# ╠═9e1c4bc8-882f-45bf-b437-38c369867dac
# ╠═15c7c5ed-5bcf-4336-aecb-22fd3ef38d05
# ╠═8de88cd0-e326-4926-bfe9-776e3516d6e2
# ╠═cde76a39-51a0-4398-ba2c-81170e395a07
# ╠═42a81a70-f8a1-4a11-a475-c4adddefde73
# ╠═0969b685-e9c5-44d5-8698-a141dd01dca7
# ╠═930902f1-5d0d-41a9-8c37-07e0d913aa6e
# ╠═96f05537-6fb3-45ce-8d75-28efdfa03279
# ╠═3708eb56-53e4-43b8-8af4-e1f6e302a0cf
# ╠═07da284a-8f21-4088-b936-244d11f517dd
# ╠═9d14bb2f-2c8e-4b45-9775-e168a11005ff
# ╠═5ce17f7f-7e93-4c43-9f75-b3b677af9fd1
# ╠═71c210f8-34ad-469a-89fa-dc99efb6b72a
# ╠═d5db3698-8546-4606-8346-85475619a3e0
# ╠═6a178a79-ae85-4496-abef-11bf62d64da8
# ╠═d9d417f3-a427-409e-9090-c76c8226be04
# ╠═0d45d6f9-f27c-4d8e-b117-26132418be55
# ╠═037b15f6-0858-41fa-a8fd-c6cc8c8d8e7a
# ╟─df004387-ea5f-42b1-9409-7b37d1a5031a
# ╠═01fb5888-a01a-4dd8-af08-bb70894a99b8
# ╠═66c3e374-389c-4be9-bf16-341565d4b4c9
# ╠═f578aafc-318f-4dc4-a574-0843057b1204
# ╠═efefaff1-e6c7-46fa-bfa7-b4e6dec67641
# ╠═262879f7-9441-4305-899a-b7c2881958fe
# ╠═d2414490-d2a5-4862-81e8-8ae8ffc5213a
