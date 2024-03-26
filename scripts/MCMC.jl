### A Pluto.jl notebook ###
# v0.19.40

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
using Turing, StatsPlots, Distributions

# ╔═╡ 2c4d5454-b734-4ecb-8f71-d38496951307
using Plots, PlutoUI, LaTeXStrings, LinearAlgebra, Random

# ╔═╡ e87dd9a4-f4ed-46d7-9f9d-dae0cadb7d05
md"""
# Bayesian reasoning and advanced sampling methods

ADD INTRODUCTION HERE

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
	T = Tridiagonal(0.3ones(9), 0.1ones(10), 0.6ones(9))
	T ./= sum(T, dims=2)
end;

# ╔═╡ 7915a972-df48-462c-bcc7-2ca94a98e441
T

# ╔═╡ 53799772-c8a1-4055-a01e-42f147ffc253
md"Starting from a single state, we can model how the state vector evolves to a the stationary distribution."

# ╔═╡ b3c45c85-d075-452a-a827-b7d1e00a17a1
md"Start from state $(@bind i0 Select(1:10))"

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

# ╔═╡ b4eef7a8-a61c-471d-baf0-1d0b5fb997e6
plot_post = heatmap(-15:0.1:15, 0.1:0.02:5, norm_posterior, color=:speed, ylab=L"\sigma", xlab=L"\mu", title="Posterior")

# ╔═╡ 63af963c-a304-4142-8f23-a56a5fc14e76
chain_MH = sample(MersenneTwister(1), uncertain_normal(y), MH(), 1_000)

# ╔═╡ a9e11333-aab9-4e63-aee4-6bceb6b771bc
summarize(chain_MH)

# ╔═╡ 6309cc65-59d4-4df2-b756-30c6e2f877cd
quantile(chain_MH)

# ╔═╡ 1b92c00e-8b2a-4cc8-b6f2-23d91cf4f989
chain_MH2 = sample(MersenneTwister(1), uncertain_normal(y), my_MH, 1_000);

# ╔═╡ 138bf97c-9295-440d-99f4-24f7fe58b774
Σ = [σ₁^2 σ₁*σ₂*ρ;
	 σ₁*σ₂*ρ σ₂^2]

# ╔═╡ 87b984ee-2daa-4cca-b5eb-60c4aa3b2fb5
mvn = MultivariateNormal(μ, Σ)

# ╔═╡ 06e83a23-fa3a-49f7-8f43-d46f3285e943
chain_Gibbs = sample(MersenneTwister(1), uncertain_normal(y), Gibbs(MH(:μ), MH(:σ)), 1000)

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

# ╔═╡ e6201f93-d5fc-4f33-9e3c-f7d8eb0146b5
md"""
## Analysing chains and diagnostics

TO BE COMPLETED
"""

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
	plot!(p->n * seed_likelihood(p), 0, 1, lw=2, label=labels ? "likelihood (rescaled)" : "", color="blue")
	plot!(seed_posterior, 0, 1, lw=2, label=labels ? "posterior" : "", color="orange")
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

# ╔═╡ 34cf42a4-b38d-4a94-89ee-6fb85c44acbd
q = TruncatedNormal(2.5, 3, -1, 10)

# ╔═╡ efefaff1-e6c7-46fa-bfa7-b4e6dec67641
plots = Dict()

# ╔═╡ 0bc4cdff-f86d-4174-8aac-d2ca61717f96
plots["seeds_likelihood"] = plot(seed_likelihood, 0, 1, lw=2, label="likelihood", xlab=L"p", color="blue")

# ╔═╡ 08d33254-5013-4e4f-ae0b-1f5b6f3a2fb2
plots["seeds_prior"] = plot(seed_prior, 0, 1, lw=2, label="prior", xlab=L"p", color="green")

# ╔═╡ 217edb24-a555-4deb-ab58-2969d8564a7b
plots["seeds_posterior"] = plot(seed_posterior, 0, 1, lw=2, label="posterior", xlab=L"p", color="orange")

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
plots["norm_muprior"] = plot(dist2pdf(Normal(0, 20)), -50, 50, xlabel=L"\mu", label="Normal(0, 20)", title="prior for μ", lw=2)

# ╔═╡ cd8bf990-4f36-4a56-83c2-c58217448b1e
plots["norm_sigmaprior"] = plot(dist2pdf(InverseGamma(.1)), 0, 20, xlabel=L"\sigma", label="InverseGamma(1)", title="prior for σ", lw=2)

# ╔═╡ a09eda2e-0bb7-46c9-8d60-2b4447257772
plots["norm_multiprior"] = heatmap(-15:0.1:15, 0.1:0.02:5, norm_prior, color=:speed, ylab=L"\sigma", xlab=L"\mu", title="Prior")

# ╔═╡ 3b821e04-a27a-489f-924e-e82a08693371
plots["norm_likelihood"] = heatmap(-15:0.1:15, 0.1:0.02:5, norm_likelihood, color=:speed, ylab=L"\sigma", xlab=L"\mu", title="Likelihood")

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
	
	prs = plot(x->pdf(p_univar, x), -1, 10, lw=2, label="target p(x)", xlabel=L"x", title="Rejection sampling\np(accept)=$(round(acceptance_probability, digits=2))")
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
plots["MC_ss"] = bar(π_mc, xlab="state", ylab="PMF", label="", title="Stationary distribution birth-death process")

# ╔═╡ bfd58734-b396-4424-9343-4562932f70e3
plots["MC_ss_conv"] = groupedbar(0:30, vcat([p₀' * T^i for i in 0:30]...), bar_position = :stack, xlab = L"t", ylab = "fraction in state i", xticks=0:30, label=reshape(["state $i" for i in 1:size(T,2)], 1, :), title="State evolution in a birth-death process")

# ╔═╡ ecfc2307-6926-42c1-b08d-6ff93502c5cd
let
	mc = ifelse.(acc_mh, :blue, :orange)
	ind = 0:length(x_mh)-1
	
	p = plot(x->pdf(p_univar, x), -1, 10, lw=2, label="", xlabel=L"x", title="Metropolis-Hastings", ylab="PDF")
	
	plot!(twinx(), x_mh, ind, lw=0.5, label="", ylab="sample number")
	scatter!(twinx(), x_mh[acc_mh], ind[acc_mh], label="accepted", color=:blue, ms=2)
	scatter!(twinx(), x_mh[.!acc_mh], ind[.!acc_mh], label="rejected", color=:orange, ms=2)
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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
StatsPlots = "f3b207a7-027a-5e70-b257-86293d7955fd"
Turing = "fce5fe82-541a-59a6-adf8-730c64b5f9a0"

[compat]
Distributions = "~0.25.107"
LaTeXStrings = "~1.3.1"
Plots = "~1.39.0"
PlutoUI = "~0.7.58"
StatsPlots = "~0.15.7"
Turing = "~0.30.5"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.0"
manifest_format = "2.0"
project_hash = "de341c8265cb7d85b638d1805df494129ecbe989"

[[deps.ADTypes]]
git-tree-sha1 = "41c37aa88889c171f1300ceac1313c06e891d245"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "0.2.6"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractMCMC]]
deps = ["BangBang", "ConsoleProgressMonitor", "Distributed", "FillArrays", "LogDensityProblems", "Logging", "LoggingExtras", "ProgressLogging", "Random", "StatsBase", "TerminalLoggers", "Transducers"]
git-tree-sha1 = "b0489adc45a7c8cf0d8e2ddf764f89c1c3decebd"
uuid = "80f14c24-f653-4e6a-9b94-39d6b0f70001"
version = "5.2.0"

[[deps.AbstractPPL]]
deps = ["AbstractMCMC", "DensityInterface", "Random", "Setfield", "SparseArrays"]
git-tree-sha1 = "917ad8da4becae82028aba80b7e25197f0c76dd1"
uuid = "7a57a42e-76ec-4ea3-a279-07e840d6d9cf"
version = "0.7.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "c278dfab760520b8bb7e9511b968bf4ba38b7acc"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.3"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "LinearAlgebra", "MacroTools", "Test"]
git-tree-sha1 = "cb96992f1bec110ad211b7e410e57ddf7944c16f"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.35"

    [deps.Accessors.extensions]
    AccessorsAxisKeysExt = "AxisKeys"
    AccessorsIntervalSetsExt = "IntervalSets"
    AccessorsStaticArraysExt = "StaticArrays"
    AccessorsStructArraysExt = "StructArrays"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    Requires = "ae029012-a4dd-5104-9daa-d747884805df"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "0fb305e0253fd4e833d486914367a2ee2c2e78d0"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.0.1"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AdvancedHMC]]
deps = ["AbstractMCMC", "ArgCheck", "DocStringExtensions", "InplaceOps", "LinearAlgebra", "LogDensityProblems", "LogDensityProblemsAD", "ProgressMeter", "Random", "Requires", "Setfield", "SimpleUnPack", "Statistics", "StatsBase", "StatsFuns"]
git-tree-sha1 = "dfa0e3508fc3df81d28624b328f3b937c1df8bc2"
uuid = "0bf59076-c3b1-5ca4-86bd-e02cd72cde3d"
version = "0.6.1"

    [deps.AdvancedHMC.extensions]
    AdvancedHMCCUDAExt = "CUDA"
    AdvancedHMCMCMCChainsExt = "MCMCChains"
    AdvancedHMCOrdinaryDiffEqExt = "OrdinaryDiffEq"

    [deps.AdvancedHMC.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    MCMCChains = "c7f686f2-ff18-58e9-bc7b-31028e88f75d"
    OrdinaryDiffEq = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"

[[deps.AdvancedMH]]
deps = ["AbstractMCMC", "Distributions", "FillArrays", "LinearAlgebra", "LogDensityProblems", "Random", "Requires"]
git-tree-sha1 = "16589dbdd36c782ff01700908e962b303474f641"
uuid = "5b7e9947-ddc0-4b3f-9b55-0d8042f74170"
version = "0.8.1"
weakdeps = ["DiffResults", "ForwardDiff", "MCMCChains", "StructArrays"]

    [deps.AdvancedMH.extensions]
    AdvancedMHForwardDiffExt = ["DiffResults", "ForwardDiff"]
    AdvancedMHMCMCChainsExt = "MCMCChains"
    AdvancedMHStructArraysExt = "StructArrays"

[[deps.AdvancedPS]]
deps = ["AbstractMCMC", "Distributions", "Random", "Random123", "Requires", "StatsFuns"]
git-tree-sha1 = "672f7ce648e06f93fceefde463c5855d77b6915a"
uuid = "576499cb-2369-40b2-a588-c64705576edc"
version = "0.5.4"
weakdeps = ["Libtask"]

    [deps.AdvancedPS.extensions]
    AdvancedPSLibtaskExt = "Libtask"

[[deps.AdvancedVI]]
deps = ["Bijectors", "Distributions", "DistributionsAD", "DocStringExtensions", "ForwardDiff", "LinearAlgebra", "ProgressMeter", "Random", "Requires", "StatsBase", "StatsFuns", "Tracker"]
git-tree-sha1 = "1f919a9c59cf3dfc68b64c22c453a2e356fca473"
uuid = "b5ca4192-6429-45e5-a2d9-87aec30a685c"
version = "0.2.4"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra", "Logging"]
git-tree-sha1 = "9b9b347613394885fd1c8c7729bfc60528faa436"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.5.4"

[[deps.Arpack_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "5ba6c757e8feccf03a1554dfaf3e26b3cfc7fd5e"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.1+1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "c5aeb516a84459e0318a02507d2261edad97eb75"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.7.1"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Atomix]]
deps = ["UnsafeAtomics"]
git-tree-sha1 = "c06a868224ecba914baa6942988e2f2aade419be"
uuid = "a9b6321e-bd34-4604-b9c9-b65b8de01458"
version = "0.1.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.BangBang]]
deps = ["Compat", "ConstructionBase", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables"]
git-tree-sha1 = "7aa7ad1682f3d5754e3491bb59b8103cae28e3a3"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.40"

    [deps.BangBang.extensions]
    BangBangChainRulesCoreExt = "ChainRulesCore"
    BangBangDataFramesExt = "DataFrames"
    BangBangStaticArraysExt = "StaticArrays"
    BangBangStructArraysExt = "StructArrays"
    BangBangTypedTablesExt = "TypedTables"

    [deps.BangBang.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    TypedTables = "9d95f2ec-7b3d-5a63-8d20-e2491e220bb9"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.Bijectors]]
deps = ["ArgCheck", "ChainRules", "ChainRulesCore", "ChangesOfVariables", "Compat", "Distributions", "Functors", "InverseFunctions", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "MappedArrays", "Random", "Reexport", "Requires", "Roots", "SparseArrays", "Statistics"]
git-tree-sha1 = "199dc2c4151db557549a0ad8888ce1a60337ff42"
uuid = "76274a88-744f-5084-9051-94815aaf08c4"
version = "0.13.8"

    [deps.Bijectors.extensions]
    BijectorsDistributionsADExt = "DistributionsAD"
    BijectorsForwardDiffExt = "ForwardDiff"
    BijectorsLazyArraysExt = "LazyArrays"
    BijectorsReverseDiffExt = "ReverseDiff"
    BijectorsTrackerExt = "Tracker"
    BijectorsZygoteExt = "Zygote"

    [deps.Bijectors.weakdeps]
    DistributionsAD = "ced4e74d-a319-5a8a-b0ac-84af2272839c"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    LazyArrays = "5078a376-72f3-5289-bfd5-ec5146d43c02"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.BitFlags]]
git-tree-sha1 = "2dc09997850d68179b69dafb58ae806167a32b1b"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.8"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRules]]
deps = ["Adapt", "ChainRulesCore", "Compat", "Distributed", "GPUArraysCore", "IrrationalConstants", "LinearAlgebra", "Random", "RealDot", "SparseArrays", "SparseInverseSubset", "Statistics", "StructArrays", "SuiteSparse"]
git-tree-sha1 = "4e42872be98fa3343c4f8458cbda8c5c6a6fa97c"
uuid = "082447d4-558c-5d27-93f4-14fc19e9eca2"
version = "1.63.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "aef70bb349b20aa81a82a19704c3ef339d4ee494"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.22.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.ChangesOfVariables]]
deps = ["LinearAlgebra", "Test"]
git-tree-sha1 = "2fba81a302a7be671aefe194f0525ef231104e7f"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.8"
weakdeps = ["InverseFunctions"]

    [deps.ChangesOfVariables.extensions]
    ChangesOfVariablesInverseFunctionsExt = "InverseFunctions"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "Random", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "9ebb045901e9bbf58767a9f34ff89831ed711aae"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.15.7"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "59939d8a997469ee05c4b4944560a820f9ba0d73"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.4"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "c955881e3c981181362ae4088b35995446298b80"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.14.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.2+0"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "9c4708e3ed2b799e6124b5673a712dda0b596a9b"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.3.1"

[[deps.ConsoleProgressMonitor]]
deps = ["Logging", "ProgressMeter"]
git-tree-sha1 = "3ab7b2136722890b9af903859afcf457fa3059e8"
uuid = "88cd18e8-d9cc-4ea6-8889-5259c0d15c8b"
version = "0.1.2"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c53fc348ca4d40d7b371e71fd52251839080cbc9"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.4"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1fb174f0d48fe7d142e1109a10636bc1d14f5ac2"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.17"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "66c4c81f259586e8f002eacebc177e1fb06363b0"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.11"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "7c302d7a5fec5214eb8a5a4c466dcf7a51fcf169"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.107"
weakdeps = ["ChainRulesCore", "DensityInterface", "Test"]

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

[[deps.DistributionsAD]]
deps = ["Adapt", "ChainRules", "ChainRulesCore", "Compat", "Distributions", "FillArrays", "LinearAlgebra", "PDMats", "Random", "Requires", "SpecialFunctions", "StaticArrays", "StatsFuns", "ZygoteRules"]
git-tree-sha1 = "060a19f3f879773399a7011676eb273ccc265241"
uuid = "ced4e74d-a319-5a8a-b0ac-84af2272839c"
version = "0.6.54"

    [deps.DistributionsAD.extensions]
    DistributionsADForwardDiffExt = "ForwardDiff"
    DistributionsADLazyArraysExt = "LazyArrays"
    DistributionsADReverseDiffExt = "ReverseDiff"
    DistributionsADTrackerExt = "Tracker"

    [deps.DistributionsAD.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    LazyArrays = "5078a376-72f3-5289-bfd5-ec5146d43c02"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.DynamicPPL]]
deps = ["ADTypes", "AbstractMCMC", "AbstractPPL", "BangBang", "Bijectors", "Compat", "ConstructionBase", "Distributions", "DocStringExtensions", "LinearAlgebra", "LogDensityProblems", "LogDensityProblemsAD", "MacroTools", "OrderedCollections", "Random", "Requires", "Setfield", "Test"]
git-tree-sha1 = "60a3a231813a89bd796a639d5ce389386e89f4e3"
uuid = "366bfd00-2699-11ea-058f-f148b4cae6d8"
version = "0.24.7"

    [deps.DynamicPPL.extensions]
    DynamicPPLChainRulesCoreExt = ["ChainRulesCore"]
    DynamicPPLEnzymeCoreExt = ["EnzymeCore"]
    DynamicPPLForwardDiffExt = ["ForwardDiff"]
    DynamicPPLMCMCChainsExt = ["MCMCChains"]
    DynamicPPLReverseDiffExt = ["ReverseDiff"]
    DynamicPPLZygoteRulesExt = ["ZygoteRules"]

    [deps.DynamicPPL.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    MCMCChains = "c7f686f2-ff18-58e9-bc7b-31028e88f75d"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    ZygoteRules = "700de1a5-db45-46bc-99cf-38207098b444"

[[deps.EllipticalSliceSampling]]
deps = ["AbstractMCMC", "ArrayInterface", "Distributions", "Random", "Statistics"]
git-tree-sha1 = "e611b7fdfbfb5b18d5e98776c30daede41b44542"
uuid = "cad2338a-1db2-11e9-3401-43bc07c9ede2"
version = "2.0.0"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "4820348781ae578893311153d69049a93d05f39d"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random"]
git-tree-sha1 = "5b93957f6dcd33fc343044af3d48c215be2562f1"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.9.3"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Logging", "Printf"]
git-tree-sha1 = "fb409abab2caf118986fc597ba84b50cbaf00b87"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.3"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Functors]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "166c544477f97bbadc7179ede1c1868e0e9b426b"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.4.7"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "ff38ba61beff76b8f4acad8ab0c97ef73bb670cb"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.9+0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "ec632f177c0d990e64d955ccc1b8c04c485a0950"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.6"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "27442171f28c952804dede8ff72828a96f2bfc1f"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.10"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "025d171a2847f616becc0f84c8dc62fe18f0f6dd"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.10+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "e94c92c7bf4819685eb80186d51c43e71d4afa17"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.76.5+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "ac7b73d562b8f4287c3b67b4c66a5395a19c1ae8"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.2"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "8b72179abc660bfab5e28472e019392b97d0985c"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.4"

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.InplaceOps]]
deps = ["LinearAlgebra", "Test"]
git-tree-sha1 = "50b41d59e7164ab6fda65e71049fee9d890731ff"
uuid = "505f98c9-085e-5b2c-8e89-488be7bf1f34"
version = "0.3.0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5fdf2fe6724d8caabf43b557b84ce53f3b7e2f6b"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.0.2+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "88a101217d7cb38a7b481ccd50d21876e1d1b0e0"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.1"
weakdeps = ["Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsUnitfulExt = "Unitful"

[[deps.IntervalSets]]
git-tree-sha1 = "dba9ddf07f77f60450fe5d2e2beb9854d9a49bd0"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.10"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "68772f49f54b479fa88ace904f6127f0a3bb2e46"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.12"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "a53ebe394b71470c7f97c2e7e170d51df21b17af"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.7"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3336abae9a713d2210bb57ab484b1e065edd7d23"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.2+0"

[[deps.KernelAbstractions]]
deps = ["Adapt", "Atomix", "InteractiveUtils", "LinearAlgebra", "MacroTools", "PrecompileTools", "Requires", "SparseArrays", "StaticArrays", "UUIDs", "UnsafeAtomics", "UnsafeAtomicsLLVM"]
git-tree-sha1 = "c7753cc3febe006708ce6798482004241f7d890b"
uuid = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
version = "0.9.17"

    [deps.KernelAbstractions.extensions]
    EnzymeExt = "EnzymeCore"

    [deps.KernelAbstractions.weakdeps]
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "fee018a29b60733876eb557804b5b109dd3dd8a7"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.8"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVM]]
deps = ["CEnum", "LLVMExtra_jll", "Libdl", "Preferences", "Printf", "Requires", "Unicode"]
git-tree-sha1 = "ddab4d40513bce53c8e3157825e245224f74fae7"
uuid = "929cbde3-209d-540e-8aea-75f648917ca0"
version = "6.6.0"

    [deps.LLVM.extensions]
    BFloat16sExt = "BFloat16s"

    [deps.LLVM.weakdeps]
    BFloat16s = "ab4f0b2a-ad5b-11e8-123f-65d77653426b"

[[deps.LLVMExtra_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl", "TOML"]
git-tree-sha1 = "88b916503aac4fb7f701bb625cd84ca5dd1677bc"
uuid = "dad2f222-ce93-54a1-a47d-0025e8a3acab"
version = "0.0.29+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LRUCache]]
git-tree-sha1 = "b3cc6698599b10e652832c2f23db3cab99d51b59"
uuid = "8ac3fa9e-de4c-5943-b1dc-09c6b5f20637"
version = "1.6.1"
weakdeps = ["Serialization"]

    [deps.LRUCache.extensions]
    SerializationExt = ["Serialization"]

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "f428ae552340899a935973270b8d98e5a31c49fe"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.1"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LeftChildRightSiblingTrees]]
deps = ["AbstractTrees"]
git-tree-sha1 = "fb6803dafae4a5d62ea5cab204b1e657d9737e7f"
uuid = "1d6d02ad-be62-4b6b-8a6d-2f90e265016e"
version = "0.2.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtask]]
deps = ["FunctionWrappers", "LRUCache", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "345a40c746404dd9cb1bbc368715856838ab96f2"
uuid = "6f1fad26-d15e-5dc8-ae53-837a1d7b8c9f"
version = "0.8.6"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e5edc369a598dfde567269dc6add5812cfa10cd5"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.39.3+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogDensityProblems]]
deps = ["ArgCheck", "DocStringExtensions", "Random"]
git-tree-sha1 = "f9a11237204bc137617194d79d813069838fcf61"
uuid = "6fdf6af0-433a-55f7-b3ed-c6c6e0b8df7c"
version = "2.1.1"

[[deps.LogDensityProblemsAD]]
deps = ["DocStringExtensions", "LogDensityProblems", "Requires", "SimpleUnPack"]
git-tree-sha1 = "9c50732cd0f188766b6217ed6a2ebbdaf9890029"
uuid = "996a588d-648d-4e1f-a8f0-a84b347e47b1"
version = "1.7.0"

    [deps.LogDensityProblemsAD.extensions]
    LogDensityProblemsADADTypesExt = "ADTypes"
    LogDensityProblemsADEnzymeExt = "Enzyme"
    LogDensityProblemsADFiniteDifferencesExt = "FiniteDifferences"
    LogDensityProblemsADForwardDiffBenchmarkToolsExt = ["BenchmarkTools", "ForwardDiff"]
    LogDensityProblemsADForwardDiffExt = "ForwardDiff"
    LogDensityProblemsADReverseDiffExt = "ReverseDiff"
    LogDensityProblemsADTrackerExt = "Tracker"
    LogDensityProblemsADZygoteExt = "Zygote"

    [deps.LogDensityProblemsAD.weakdeps]
    ADTypes = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
    BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "18144f3e9cbe9b15b070288eef858f71b291ce37"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.27"
weakdeps = ["ChainRulesCore", "ChangesOfVariables", "InverseFunctions"]

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.MCMCChains]]
deps = ["AbstractMCMC", "AxisArrays", "Dates", "Distributions", "Formatting", "IteratorInterfaceExtensions", "KernelDensity", "LinearAlgebra", "MCMCDiagnosticTools", "MLJModelInterface", "NaturalSort", "OrderedCollections", "PrettyTables", "Random", "RecipesBase", "Statistics", "StatsBase", "StatsFuns", "TableTraits", "Tables"]
git-tree-sha1 = "d0ce57aa5ebbdb456bac3bc5a2ca15cd06ec5f1b"
uuid = "c7f686f2-ff18-58e9-bc7b-31028e88f75d"
version = "6.0.5"

[[deps.MCMCDiagnosticTools]]
deps = ["AbstractFFTs", "DataAPI", "DataStructures", "Distributions", "LinearAlgebra", "MLJModelInterface", "Random", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Tables"]
git-tree-sha1 = "8ba8b1840d3ab5b38e7c71c23c3193bb5cbc02b5"
uuid = "be115224-59cd-429b-ad48-344e309966f0"
version = "0.3.10"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "72dc3cf284559eb8f53aa593fe62cb33f83ed0c0"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.0.0+0"

[[deps.MLJModelInterface]]
deps = ["Random", "ScientificTypesBase", "StatisticalTraits"]
git-tree-sha1 = "14bd8088cf7cd1676aa83a57004f8d23d43cd81e"
uuid = "e80e1ace-859a-464e-9ed9-23947d8ae3ea"
version = "1.9.5"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.MicroCollections]]
deps = ["BangBang", "InitialValues", "Setfield"]
git-tree-sha1 = "629afd7d10dbc6935ec59b32daeb33bc4460a42e"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.4"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.MultivariateStats]]
deps = ["Arpack", "LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI", "StatsBase"]
git-tree-sha1 = "68bf5103e002c44adfd71fea6bd770b3f0586843"
uuid = "6f286f6a-111f-5878-ab1e-185364afe411"
version = "0.10.2"

[[deps.NNlib]]
deps = ["Adapt", "Atomix", "ChainRulesCore", "GPUArraysCore", "KernelAbstractions", "LinearAlgebra", "Pkg", "Random", "Requires", "Statistics"]
git-tree-sha1 = "877f15c331337d54cf24c797d5bcb2e48ce21221"
uuid = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
version = "0.9.12"

    [deps.NNlib.extensions]
    NNlibAMDGPUExt = "AMDGPU"
    NNlibCUDACUDNNExt = ["CUDA", "cuDNN"]
    NNlibCUDAExt = "CUDA"
    NNlibEnzymeCoreExt = "EnzymeCore"

    [deps.NNlib.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    cuDNN = "02a925ec-e4fe-4b08-9a7e-0d78e3d38ccd"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NamedArrays]]
deps = ["Combinatorics", "DataStructures", "DelimitedFiles", "InvertedIndices", "LinearAlgebra", "Random", "Requires", "SparseArrays", "Statistics"]
git-tree-sha1 = "6d42eca6c3a27dc79172d6d947ead136d88751bb"
uuid = "86f7a689-2022-50b4-a561-43c23ac3c673"
version = "0.10.0"

[[deps.NaturalSort]]
git-tree-sha1 = "eda490d06b9f7c00752ee81cfa451efe55521e21"
uuid = "c020b1a1-e9b0-503a-9c33-f039bfc54a85"
version = "1.0.0"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "ded64ff6d4fdd1cb68dfcbb818c69e144a5b2e4c"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.16"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "6a731f2b5c03157418a20c12195eb4b74c8f8621"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.13.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "51901a49222b09e3743c65b8847687ae5fc78eb2"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "60e3045590bd104a16fefb12836c00c0ef8c7f8c"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.13+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optimisers]]
deps = ["ChainRulesCore", "Functors", "LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "264b061c1903bc0fe9be77cb9050ebacff66bb63"
uuid = "3bd65402-5787-11e9-1adc-39752487f4e2"
version = "0.3.2"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "862942baf5663da528f66d24996eb6da85218e76"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.0"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "ccee59c6e48e6f2edf8a5b64dc817b6729f99eb5"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.39.0"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "71a22244e352aa8c5f0f2adde4150f62368a3f2e"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.58"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "88b895d13d53b5577fd53379d913b9ab9ac82660"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "80d919dee55b9c50e8d9e2da5eeafff3fe58b539"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.4"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "00099623ffee15972c16111bcf84c58a0051257c"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.9.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "37b7bb7aabf9a085e0044307e1717436117f2b3b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.5.3+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9b23c31e76e333e6fb4c1595ae6afa74966a729e"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.9.4"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "4743b43e5a9c4a2ede372de7061eed81795b12e7"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.7.0"

[[deps.RandomNumbers]]
deps = ["Random", "Requires"]
git-tree-sha1 = "043da614cc7e95c703498a491e2c21f58a2b8111"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.5.3"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "Requires", "SparseArrays", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "27ee1c03e732c488ecce1a25f0d7da9b5d936574"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.3.3"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.Roots]]
deps = ["Accessors", "ChainRulesCore", "CommonSolve", "Printf"]
git-tree-sha1 = "754acd3031a9f2eaf6632ba4850b1c01fe4460c1"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.1.2"

    [deps.Roots.extensions]
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Roots.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "6aacc5eefe8415f47b3e34214c1d79d2674a0ba2"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.12"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SciMLBase]]
deps = ["ADTypes", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FillArrays", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables", "TruncatedStacktraces"]
git-tree-sha1 = "09324a0ae70c52a45b91b236c62065f78b099c37"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.15.2"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseZygoteExt = "Zygote"

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLOperators]]
deps = ["ArrayInterface", "DocStringExtensions", "Lazy", "LinearAlgebra", "Setfield", "SparseArrays", "StaticArraysCore", "Tricks"]
git-tree-sha1 = "51ae235ff058a64815e0a2c34b1db7578a06813d"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.7"

[[deps.ScientificTypesBase]]
git-tree-sha1 = "a8e18eb383b5ecf1b5e6fc237eb39255044fd92b"
uuid = "30f210dd-8aff-4c5f-94ba-8e64358c1161"
version = "3.0.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "0e7508ff27ba32f26cd459474ca2ede1bc10991f"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SimpleUnPack]]
git-tree-sha1 = "58e6353e72cde29b90a69527e56df1b5c3d8c437"
uuid = "ce78b400-467f-4804-87d8-8f486da07d0a"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SparseInverseSubset]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "52962839426b75b3021296f7df242e40ecfc0852"
uuid = "dc90abb0-5640-4711-901d-7e5b23a2fada"
version = "0.1.2"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "e08a62abc517eb79667d0a29dc08a3b589516bb5"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.15"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "bf074c045d3d5ffd956fa0a461da38a44685d6b2"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.3"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.StatisticalTraits]]
deps = ["ScientificTypesBase"]
git-tree-sha1 = "30b9236691858e13f167ce829490a68e1a597782"
uuid = "64bff920-2084-43da-a3e6-9bb72801c0c9"
version = "3.2.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "1d77abd07f617c4868c33d4f5b9e1dbb2643c9cf"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.2"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "cef0472124fab0695b58ca35a77c6fb942fdab8a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.1"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StatsPlots]]
deps = ["AbstractFFTs", "Clustering", "DataStructures", "Distributions", "Interpolations", "KernelDensity", "LinearAlgebra", "MultivariateStats", "NaNMath", "Observables", "Plots", "RecipesBase", "RecipesPipeline", "Reexport", "StatsBase", "TableOperations", "Tables", "Widgets"]
git-tree-sha1 = "3b1dcbf62e469a67f6733ae493401e53d92ff543"
uuid = "f3b207a7-027a-5e70-b257-86293d7955fd"
version = "0.15.7"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a04cabe79c5f01f4d723cc6704070ada0b9d46d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.4"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "f4dc295e983502292c4c3f951dbb4e985e35b3be"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.18"
weakdeps = ["Adapt", "GPUArraysCore", "SparseArrays", "StaticArrays"]

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = "GPUArraysCore"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.SymbolicIndexingInterface]]
git-tree-sha1 = "be414bfd80c2c91197823890c66ef4b74f5bf5fe"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableOperations]]
deps = ["SentinelArrays", "Tables", "Test"]
git-tree-sha1 = "e383c87cf2a1dc41fa30c093b2a19877c83e1bc1"
uuid = "ab02a1b2-a7df-11e8-156e-fb1833f50b87"
version = "1.2.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "cb76cf677714c095e535e3501ac7954732aeea2d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.11.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TerminalLoggers]]
deps = ["LeftChildRightSiblingTrees", "Logging", "Markdown", "Printf", "ProgressLogging", "UUIDs"]
git-tree-sha1 = "f133fab380933d042f6796eda4e130272ba520ca"
uuid = "5d786b92-1e48-4d6f-9151-6b4477ca9bed"
version = "0.1.7"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tracker]]
deps = ["Adapt", "DiffRules", "ForwardDiff", "Functors", "LinearAlgebra", "LogExpFunctions", "MacroTools", "NNlib", "NaNMath", "Optimisers", "Printf", "Random", "Requires", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "5c942be30a85ac75d14e9e527b55504031e1bbd3"
uuid = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
version = "0.2.31"
weakdeps = ["PDMats"]

    [deps.Tracker.extensions]
    TrackerPDMatsExt = "PDMats"

[[deps.TranscodingStreams]]
git-tree-sha1 = "54194d92959d8ebaa8e26227dbe3cdefcdcd594f"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.3"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "ConstructionBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "3064e780dbb8a9296ebb3af8f440f787bb5332af"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.80"

    [deps.Transducers.extensions]
    TransducersBlockArraysExt = "BlockArrays"
    TransducersDataFramesExt = "DataFrames"
    TransducersLazyArraysExt = "LazyArrays"
    TransducersOnlineStatsBaseExt = "OnlineStatsBase"
    TransducersReferenceablesExt = "Referenceables"

    [deps.Transducers.weakdeps]
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    LazyArrays = "5078a376-72f3-5289-bfd5-ec5146d43c02"
    OnlineStatsBase = "925886fa-5bf2-5e8e-b522-a9147a512338"
    Referenceables = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.Turing]]
deps = ["ADTypes", "AbstractMCMC", "AdvancedHMC", "AdvancedMH", "AdvancedPS", "AdvancedVI", "BangBang", "Bijectors", "DataStructures", "Distributions", "DistributionsAD", "DocStringExtensions", "DynamicPPL", "EllipticalSliceSampling", "ForwardDiff", "Libtask", "LinearAlgebra", "LogDensityProblems", "LogDensityProblemsAD", "MCMCChains", "NamedArrays", "Printf", "Random", "Reexport", "Requires", "SciMLBase", "Setfield", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "afb5bb484e67bb4179507baff464da9c4d18b307"
uuid = "fce5fe82-541a-59a6-adf8-730c64b5f9a0"
version = "0.30.5"

    [deps.Turing.extensions]
    TuringDynamicHMCExt = "DynamicHMC"
    TuringOptimExt = "Optim"

    [deps.Turing.weakdeps]
    DynamicHMC = "bbc10e6e-7c05-544b-b16e-64fede858acb"
    Optim = "429524aa-4258-5aef-a3af-852621145aeb"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "3c793be6df9dd77a0cf49d80984ef9ff996948fa"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.19.0"
weakdeps = ["ConstructionBase", "InverseFunctions"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.UnsafeAtomics]]
git-tree-sha1 = "6331ac3440856ea1988316b46045303bef658278"
uuid = "013be700-e6cd-48c3-b4a1-df204f14c38f"
version = "0.2.1"

[[deps.UnsafeAtomicsLLVM]]
deps = ["LLVM", "UnsafeAtomics"]
git-tree-sha1 = "323e3d0acf5e78a56dfae7bd8928c989b4f3083e"
uuid = "d80eeb9a-aca5-4d75-85e5-170c8b632249"
version = "0.1.3"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.Widgets]]
deps = ["Colors", "Dates", "Observables", "OrderedCollections"]
git-tree-sha1 = "fcdae142c1cfc7d89de2d11e08721d0f2f86c98a"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.6"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "07e470dabc5a6a4254ffebc29a1b3fc01464e105"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.5+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "37195dcb94a5970397ad425b95a9a26d0befce3a"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.6.0+0"

[[deps.Xorg_libICE_jll]]
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "e5becd4411063bdcac16be8b66fc2f9f6f1e8fe5"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.0.10+1"

[[deps.Xorg_libSM_jll]]
deps = ["Libdl", "Pkg", "Xorg_libICE_jll"]
git-tree-sha1 = "4a9d9e4c180e1e8119b5ffc224a7b59d3a7f7e18"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.3+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.ZygoteRules]]
deps = ["ChainRulesCore", "MacroTools"]
git-tree-sha1 = "27798139afc0a2afa7b1824c206d5e87ea587a00"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.5"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a68c9655fbe6dfcab3d972808f1aafec151ce3f8"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.43.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.7.0+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "1ea2ebe8ffa31f9c324e8c1d6e86b4165b84a024"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╠═103e5ba0-cfdc-11ee-13b1-cf53dfdd9a3b
# ╠═2c4d5454-b734-4ecb-8f71-d38496951307
# ╟─e87dd9a4-f4ed-46d7-9f9d-dae0cadb7d05
# ╟─0f6cf892-b218-4740-a298-43f47e51acae
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
# ╠═ff813a28-5ab7-4b49-8b05-b4b6a73b94a5
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
# ╠═35bc2b8f-716d-4a85-8e60-7242709448c3
# ╟─0c652efb-0edb-4a4e-a14c-7d5ac6eb5ffa
# ╠═ab09925c-4689-4f87-be44-f8eca12a5c9c
# ╟─69f2c623-632f-4df4-8566-a48574c0f9fb
# ╟─1a24999d-1ab7-4220-8e22-b1ac33ad03ac
# ╟─5e6c6135-f4f9-447c-8f7a-a49ec25e414a
# ╠═3fca01f3-ac6e-4fdc-b8ca-48040ca562f3
# ╟─138bf97c-9295-440d-99f4-24f7fe58b774
# ╠═87b984ee-2daa-4cca-b5eb-60c4aa3b2fb5
# ╠═76b0886a-fc91-40ef-89d5-dcf5d261907f
# ╠═06e83a23-fa3a-49f7-8f43-d46f3285e943
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
# ╠═575ad310-8cbd-4b05-aee0-d84fb55ad6ed
# ╟─4a86bf84-434a-4ef3-a9c6-4c72df5068e0
# ╠═e6201f93-d5fc-4f33-9e3c-f7d8eb0146b5
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
# ╠═34cf42a4-b38d-4a94-89ee-6fb85c44acbd
# ╠═efefaff1-e6c7-46fa-bfa7-b4e6dec67641
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
