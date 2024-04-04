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

# ╔═╡ 2d346892-cb15-11ee-2d81-73e08fcc3288
# ╠═╡ skip_as_script = true
#=╠═╡
begin
    using Pkg
	Pkg.activate("..")
end
  ╠═╡ =#

# ╔═╡ cf8dbe5e-d340-4ebb-8ad7-f483be07fffd
	using PlutoUI, Plots, LinearAlgebra, Markdown, Random, LaTeXStrings

# ╔═╡ 4d96d7ea-1216-46d6-80fe-5859640eb9c9
using Distributions

# ╔═╡ 2b79770c-00a8-4c3b-b39e-b31400593ea4
using Sobol

# ╔═╡ ba4460b1-2855-4695-a3d8-4976666ef164
using Turing

# ╔═╡ cfc49c56-f223-4b78-b6da-4e0c83a16bbf
md"""
# Monte Carlo methods and modeling with probability distributions

In this chapter, we study random sampling to obtain certain numerical results that might be challenging to obtain in other ways. Within the context of modelling and simulation, these methods are often referred to as *Monte Carlo methods*. 

For example, suppose you want to estimate $\pi$ using a Monte Carlo method. First, you draw the unit circle ($x^2+y^2=1$). Then, you randomly pick a point in the unit square, i.e., in $[0,1]\times[0,1]$ (like throwing a dart!). The unit square has an area of $1\times 1=1$ and contains a quarter of the circle with an area of $(1/2)^2\pi=\pi/4$. Hence, the probability of a random throw falling in this circle quarter is $\pi/4$. In a simulation using 1000 throws, 751 ended up in the quarter circle. So, our estimate would be 751/1000 = 0.751, which is not too far from $\pi/4\approx 0.7853$. We can improve our estimation accuracy as much as we like by just performing more samples.

The following recipe gives the general Monte Carlo method:
1. randomly sample a parameter value $\theta\sim P(\theta)$;  
2. perform a simulation (called a *throw* or *simulation*) to generate an using the model (i.e. $P(X\mid \theta)$); 
3. repeat steps 1 & 2 times to obtain a statistically representative sample;  
4. combine all simulation results for statistical analysis.

"""

# ╔═╡ c469b254-ffc4-458d-8676-db135e11b306
md"n throws : $(@bind n_pi Slider(1000:1000:100_000, default=5000, show_value=true))"

# ╔═╡ 5e81b02d-a879-4023-a7fe-e6f187a558f1
n_pi  # number of throws

# ╔═╡ 9b3f076f-e50b-4f23-a282-91fac47de21b
md"Generate points in $[0,1]\times[0,1]$."

# ╔═╡ 5db3ed2d-ee29-4cd5-b438-b60d550a6ff8
x_unif = rand(n_pi)

# ╔═╡ 7d1e1a93-510c-4f81-bb09-d87e3b40e136
y_unif = rand(n_pi)

# ╔═╡ 5821aaf1-744d-452a-96d2-6516c7423fff
in_circle = x_unif.^2 .+ y_unif.^2 .≤ 1

# ╔═╡ cb36a752-da05-4da8-8dc8-1c2d34a11ad5
sum(in_circle)

# ╔═╡ 7c244dc0-6c7b-451e-b7c0-da25a77da657
pi_est = 4sum(in_circle) / length(in_circle)

# ╔═╡ fe471a6a-5c29-49f6-afba-83be7ab5d770
π  # true value

# ╔═╡ e0f96ed2-bc51-467e-be3e-8617293e33e2
md"""
## Revision of basic probability theory

### Random variables and probability distributions
Let us start with revising the basics of probability theory. We will work with random or stochastic variables, usually denoted with a capital letter such as $X$ or $Y$. These can either be *continuous*, i.e., defined over the reals or some continuous subset of them, or can be *discrete*. The outcome of a random variable, for example when one performs an experiment to gain a measurement, will be denoted using the corresponding lowercase letter, here $x$ and $y$. 

Outcomes of random variables are assigned probabilities via probability functions. For discrete random variables, these are called *probability mass functions* (PMF):

$$p_X(x) = P(X=x)\,,$$

which should satisfy that for every $x$ it holds that  $p_X(x)\ge 0$ (probabilities are non-zero) and $\sum_x p_X(x) = 1$ (probabilities should be normalized). For example, for a simple fair die, the outcomes $x$ are $\{1,2,3,4,5,6\}$ with respective probabilities $\{1/6,1/6,1/6,1/6,1/6,1/6\}$.

Much of our work will however deal with continuous variables, which are characterized using *probability density functions* (PDF) $f_X(x)$, which likewise are non-negative and normalized:

$$\int_{-\infty}^\infty f_X(x) \mathrm{d}x = 1\,.$$

For example, we take a Weibull distribution:
"""

# ╔═╡ c0677edb-29a4-4c4d-8ca7-9d81eac05c85
dist = Weibull(2, 7)

# ╔═╡ d786da44-5cfc-4c36-ae41-af6f623d5a4e
pdf(dist, 5.0)  # density in x=5.0

# ╔═╡ bce61d18-767d-4ae4-a211-759ef3e4da44
f_X(x) = pdf(dist, x)  # handy function for the pdf

# ╔═╡ 1f9b17a8-d8f1-4256-bebc-d5f92e085fe8
f_X(5.0)

# ╔═╡ fce1c003-6483-44f1-bdbe-c6363d09807b
md"""
A probability density function is a bit harder to relate to probabilities. For example, $f_X(x)$ generally does not give the probability that $X=x$, though $f_X(x)/f_X(x')$ can be interpreted as the ratio probabilities of both both events. In general, probabilities are obtained by integrating to compute the probability of $X$ giving a value in some region $R$:

$$\int_R f_X(x) \mathrm{d}x = P(X\in R)\,.$$

Those who followed a statistics course might have ingrained that for a standard normal distribution, the probability of seeing a value in $[-2,2]$ is 0.954 (hence an approximate 95% confidence interval is the estimate $\pm$ two times the standard error).

Often, it is easier to work with *cumulative density functions* (CDF), which directly yield probabilities:

$$F_X(x) = \int_{-\infty}^x f_X(s) \mathrm{d}s = P(X \le x)\,,$$

from which one can directly compute useful quantities, such as the probability of observing a value in some interval $[a,b]$:

$$P(a < X\le b) = \int_a^b f_X(x) \mathrm{d}x= F_X(b)-F_X(a)\,.$$
"""

# ╔═╡ dd8c1370-5fb9-4ab7-8729-ea756762788d
cdf(dist, 5.0)  # P(X ≤ 5)

# ╔═╡ 1882bc65-abca-45b7-8f53-e76360b8dacd
F_X(x) = cdf(dist, x)

# ╔═╡ d04ea704-c9c5-4d57-8908-a7b19b5567a0
F_X(5) - F_X(1)  # P(1 < X ≤ 5)

# ╔═╡ 2538dd4a-96a3-4d4f-b03e-82d2e8b9ea1f
md"x : $(@bind xpdf Slider(0:0.2:20, default=5, show_value=true))"

# ╔═╡ fff75ceb-b21c-49fe-b1c0-0d2f3213b37d
md"""
Most probability distributions have one or more parameters, which are generally lumped as $\theta$. We can write this explicitly as $f_X(x;\theta)$ if we want to stress the parameterization. Familiar examples include the mean and standard deviation of a normal distribution ($\mu$ and $\sigma$) or the rate parameter $\lambda$ of a Poisson distribution. These parameters are either known, or we have to estimate them from data using, for example, the maximum likelihood principle. Or, as we will see later, we can assign probability distributions to the parameters themselves to encode our beliefs what these parameters would be. We will usually omit the dependency on the parameter and just write $f_X(x)$ if the parameters are known.

Distributions can also be defined over two or more variables. For example, the *joint distribution* of the ensemble $X,Y$ would be written as $f_{XY}(x,y)$. We can obtain the *marginal distribution* of $X$ by integrating over all values of $Y$:

$$F_X(x) = \int_{-\infty}^\infty f_{XY}(x,y)\mathrm{d}y\,.$$
"""

# ╔═╡ 6ae427c7-a65e-4f9b-b5a7-8cb523fcd382
md"""
### Expected value and other moments
A key concept is the *expected value* of a distribution:

$$E[X] = \int_{-\infty}^\infty xf_X(x) \mathrm{d}x\,.$$
"""

# ╔═╡ 090d0c6f-de9f-4ce8-befa-4ae9dbe1b38e
mean(dist)

# ╔═╡ 70cfc225-dd83-472c-b745-5266026c3bfc
md"""
Think of this as the center of mass or the first moment of the distribution. One often uses the symbol $\mu$ to denote the expected value. For discrete variables, replace the integral sign with a sum. The expected value is a linear operator (because integrating is a linear operator), so it holds that

$$E[aX+bY] = aE[X]+bE[Y]\,,$$

for any constants $a$ and $b$.

We are often interested in the expected value of functions of our random variable:

$$E[g(x)] = \int_{-\infty}^\infty g(x)f_X(x) \mathrm{d}x\,.$$

The above identity is sometimes called the *Law of the Unconscious Statistician*. For example, if $X$ represents the distribution of the radii of particles in an emulsion, we might compute the average area or volume. However, $g(\cdot)$ might be much more complex! It could represent a simulator described by a large system of ordinary differential equations, for example, a complex, multi-step industrial process in food industry to process particles into delicious foods. The reader might already feel that in this case it becomes rather hard to compute this integral analytically.
"""

# ╔═╡ 85e08b91-10f9-4808-ac53-195f90387688
md"""
In addition to the expected value, which measures the location of the distribution, we often want to know the "spread" of the distribution, this is given by the *variance* (often denoted by $\sigma^2$), which is defined as the expected squared difference with the mean:

$$\text{Var}[X] = E[(X-E[X])^2]\,.$$

A little algebra yields the useful formula 

$$\text{Var}[X] = E[X^2] - E[X]^2\,.$$
"""

# ╔═╡ 36b475e9-a16b-464f-a423-a4fbfbf21ef0
var(dist)

# ╔═╡ f1ee69eb-d8e5-4a6b-b408-a52d33a09e52
md"Variance is not a linear operator. If $X$ and $Y$ are *independent*, it holds that

$$\text{Var}[aX+bY] = a^2\text{Var}[X]+b^2\text{Var}[Y]\,.$$

Note that the units of variance also differ. If $X$ is given in $\mu$m, then $\text{Var}[X]$ is given in $\mu$m$^2$. For this reason, one often reports the *standard deviation* $\sigma$, which is the square root of the variance.  The latter is easier to interpret, though the variance is easier to compute with."

# ╔═╡ ab4afb22-28e1-4d9e-8d89-a1a9f0741195
std(dist)

# ╔═╡ 4d6e67d1-42b8-481f-b9fe-3df524d4b07b
sqrt(var(dist))  # same!

# ╔═╡ 2a30c7eb-de92-46a1-ac79-b15529e8f463
md"The package `Distributions` also have several other convienient functions."

# ╔═╡ 83db8998-125e-4cae-b6fe-1cdddfe3f868
median(dist)

# ╔═╡ 9f525fef-80b3-4417-85ea-e318ad387985
mode(dist)

# ╔═╡ d249188a-a298-44b8-8dcb-691753331439
md"""
### Conditional probabilities
Conditional probabilities allow to update probabilities of events in light of new information. For example:
- Given your probability of having COVID, how does this change when you test positive?
- What is the chance or rolling a six using a fair die, given that you have thrown an even number?
- What is the distribution of length for people who weigh 75 kg?
- What is the distribution of length, given that everyone is greater than 1.8 meter?

Given two events $A$ (for example, having COVID) and $B$ (for example, testing positive for COVID), the *conditional probability* of $A$ given $B$ is defined as

$$P(A\mid B) = \frac{P(A\cap B)}{P(B)}\,,$$

with $P(A\cap B)$ the chance of $A$ and $B$ occurring simultaneously. In the space of possible events, conditional probabilities "zoom in" to worlds where $B$ happens, hence $P(B)$ is used for renormalizing to probability of $A$ and $B$ happening.

When modelling, we often know conditional probabilities better than joint probabilities, so this version variant is often used:

$$P(A\cap B)= P(B\mid A)P(A)\,.$$

Keeping our example, we often know from health statistics the $P(A)$, the probability of a random person having COVID, and likewise the reliability a test, the probability of having a positive test given that you have the disease is often also known. Their product is the probability of both having the disease and testing positive. 

We know the the initial probability of a disease ($P(A)$) and the chance a diseased person will test positive  ($P(B\mid A)$) is also known. The thing that we are interested in though, is the chance of having the disease, given that we test positive ($P(A\mid B)$)! This probability can be computed by applying the rule for conditional probability twice:

$$P(A\mid B) = \frac{P(B\mid A)P(A)}{P(B)}\,.$$

The above equation is the famous *Bayes' theorem*, one of the most powerful equations in applied mathematics. It allows to flip conditional probabilities and make inferences based on data. Astute readers will notice however that applying this rule is not immediately obvious, we don't know $P(B)$ (chance of testing positive), which will depend on the probability that one has the disease or not and the probabilities and how reliable the test is when the patient has COVID or not! In your course of probability theory, you have seen that you can compute $P(A)$ using the law of total probability for simple cases as this example. For this course, we would like to stress the following:

> We have $P(A\mid B)\propto P(B\mid A)P(A)$, though the normalization constant may be hard or even impossible to compute! We will need clever algorithms to condition.

Conditional distributions can also be derived from the joint distribution:

$$f_{X\mid Y=y}(x,y) = \frac{f_{XY}(x,y)}{f_Y(y)}=\frac{f_{Y\mid X=x}(x,y)f_X(x)}{f_Y(y)}\,.$$

Note that to compute the denominator, we need to evaluate
$$\int_{-\infty}^\infty f_{XY}(x,y)\mathrm{d}x\,.$$

We can also define conditional expectations (or other quantities):

$$E[X\mid Y=y]=\int_{-\infty}^\infty xf_{X\mid Y=y}(x) \mathrm{d}x\,.$$
Often, we shall use $X$ as a variable of interest and $Y$ the outcome of some noisy, indirect measurement of relating to $X$. Conditional probabilities will allow us to learn about $X$ given $Y$.
"""

# ╔═╡ 3b9e8f64-d307-48bc-9a71-1ab9accd134e
md"""
## Sampling from a distribution

### The Law of Large Numbers

The process of generating a value of a distribution, either computationally or by performing an experiment, is called *sampling*. A sample $x$ is thus a realization of a random variable $X$ and we will denote this as either
$$x\sim X$$
or 
$$x\sim f_X(x)\,.$$
depending on the context. In this course, we will be somewhat flexible with our notation.

If you take a very large number of samples, you expect their histogram (in case of a continuous distribution) of their frequencies closely match the PDF or PMF. One of the most important concepts in probability theory (or science in general) is the *Law of Large Numbers* (LLN):

> The average of a **large** number of samples from a probability distribution closely matches the expected value of that distribution.

For an average, one often uses the following notation:

$$\bar{X}_n =\frac{X_1+ X_2 +\ldots+X_n}{n} = \frac{1}{n}\sum_{i=1}^n X_i\,,$$

where $X_1, X_2 ,\ldots,X_n$ represent independent and identically distributed (i.i.d.) draws from $X$. The LLN then states

$$\bar{X}_n \rightarrow E[X]\,.$$

There are two version of this law, the weak and the strong law of large numbers:
- **Weak law** means convergence in probability: $\lim_{n\rightarrow\infty}P(|\bar{X}_n-\mu| < \varepsilon)=1$, for any $\varepsilon > 0$. So no matter how small you take a margin, the average will certainly get closer to the expected value.
- **Strong law** means that the sample average will almost surely converge to the true average: $P(\lim_{n\rightarrow\infty}\bar{X}_n=\mu)=1$.

Note that the LLN only holds if the mean and variance are finite. For example, what would happen if you wanted to estimate the average mass of objects in our solar system by collecting a large number of them (cars, asteroids, particles, etc.)? Might some objects skew the mean?
"""

# ╔═╡ 1a991ac1-d4ed-4774-a9e2-a5f27f99e1ff
rand(dist)  # one sample from a distribution

# ╔═╡ b5b55083-2fe3-4620-9759-8a2329154143
rand(dist, 10)  # ten samples from a distribution

# ╔═╡ f625d656-35de-47cf-985b-3ae9fa54c5d4
md"""
So, if we have an easy way of sampling from a distribution, we can closely estimate expected values of our distribution. Even better, the Law of the Unconscious Statistician gives us a way of computing $E[g(X)]$ by merely applying our function $g(\cdot)$ to the samples before averaging. The implications are numerous:
- if you want to estimate the probability of $X\in R$, just count the number of occurrences where the samples are in this region;
- the sample variance will also approximate the true variance;
- a histogram will converge to the PDF, as the bin heights correspond the to expected fraction of observations in this interval.

The difference between this course and your statistics course, is that while statisticians are very worried about how close their sample average is to the true population average (because experiments are expensive, so there usually are a limited number of observations) we will have a much more causal life view. Since our simulations are done _in silico_, they are (relatively) cheap, so we assume that we have enough observations that are close to the expected value we are interested in."""

# ╔═╡ 482d1da0-3415-425c-84de-69f7fb821934
g(r) = 4π * r^2  # area sphere

# ╔═╡ 7008afc0-785e-45a8-9b6d-5b8f26bae09d
dr = 1e-4

# ╔═╡ 64e0265f-a621-495f-b73a-0af23b0ef1ac
Ā_integrate = sum(r->g(r) * pdf(dist, r) * dr, 0:dr:50)

# ╔═╡ 18a1b129-eac5-4080-8af9-7d9fbf800107
g(mean(dist))  # g(E[X])

# ╔═╡ 5a1522f7-411d-40f5-92e8-b5a154b8d735
Ā_sampling = mean(g, rand(dist, 100_000))  # E[g(x)]

# ╔═╡ 8a493f5e-1f1b-4cf5-91e1-5504593cee9d
md"""
### Generating random numbers

Given a probability distribution, we often want to *sample* from it, meaning that we would like draw random numbers that follow the given probability density or mass function. This practically means that if you make a histogram of a large number of samples from the distribution, you expect this to follow the PDF closely. Usually, one departs from random numbers from the uniform distribution in $[0, 1]$. For example, to generate a Bernoulli variable with $p=0.1$, one just takes a threshold on the uniformly randomly distributed number.

For simple distributions, one can use *inverse transform sampling*  where one uses the inverse of the CDF $F_X$ (which can either be computed analytically or represented numerically) to transform a number $U\sim \text{Unif}(0,1)$ into $F^{-1}_X(U)$. The resulting sample follows the distribution of $X$, as can easily be shown:

$$P(F^{-1}_X(U)\le x) = P(U\le F_X(x))=F_X(x)\,.$$
Libraries to work with probability distributions usually implement optimized methods for sampling from the standard distributions. 

In addition to inverse transform sampling, there exist a plethora of specialized sampling methods. Consider the normal distribution as an important special case. The *Box–Muller transform* generates two standard normally distributed values from two independent samples $U_1$ and $U_2$ from $\text{Unif}(0,1)$:

$$Z_1=\sqrt{-2\log U_1}\cos(2\pi U_2)$$
$$Z_2=\sqrt{-2\log U_2}\sin(2\pi U_1)$$
"""

# ╔═╡ dc7eb2b3-b160-4943-84e9-8cf20d687d93
md"Standard normal numbers are available in Base julia:"

# ╔═╡ 2d3e09ff-3766-46d5-abdd-4132c4fa06b9
randn()

# ╔═╡ 3ecfb175-0c26-419e-9770-d406cd1dd2dc
rand(10)

# ╔═╡ 9ba45ba1-f2d3-4d08-a3e7-2ed8730e3d90
md"Of course, you can also access them via Distributions:"

# ╔═╡ e67e0b43-1346-4d7c-8bcc-ca53949412b8
normal = Normal(0, 1)  # standard normal

# ╔═╡ ac78fddb-0b24-4956-8137-44615721bc06
rand(normal)

# ╔═╡ 44ca88b0-4256-4123-93f8-82d1ba3b3b77
rand(normal, 10)

# ╔═╡ 896da840-ad94-40d7-8a65-be92b174e88b
md"""

In the next chapter, we will consider how to sample from much more complex, composite distributions.

One might wonder how to obtain samples of $U$ in the first place! This is especially challenging on a computer, which is a deterministic machine. In practice, the random number one obtains are *pseudorandom numbers*, which result from an algorithm that generates a sequence of numbers that look random (i.e., have no statistically detectable patterns) but can be fully predicted if one has the algorithm and the starting value. This starting value is called the *seed*. One of the most frequently used psuedorandom generators is the *Mersenne Twister*, which generates a sequence of numbers based on the eponymous Mersenne primes (primes that can be written as $2^n-1$). For numerical simulations or statistical analyses using pseudorandom numbers, one often fixes the seed at the beginning of the script such that the results can be reproduced exactly. (Some students forget to set different seeds when running simulations of different cores of a cluster computer and learn to their horror that their overnight Monte Carlos simulations have all returned exactly the same result). 

For some specialized statistical studies or when safety is of concern, such as cryptography, pseudorandom numbers might not suffice, they require *true random numbers*. Such numbers can only be obtained by having a physical proces that is stochastic, such as throwing a die. One can buy *hardware random number generators* that convert a noise source, such as electrical noise or quantum processes, into random numbers. Some websites, such as https://www.random.org/ offer true random numbers based on atmospheric fluctuations for lotteries, passwords or other applications. True random number are, just like pseudorandom numbers, expensive to generate, so their free access is limited.

For some applications, (pseudo-)random numbers are too random! This is because there are likely to cover the space uniformly in expected values (and hence only for large samples). For some Monte-Carlo application, such as estimating  a volume in a space, one needs numbers that cover the space a bit more homogeneously. These are called low-discrepancy or *quasi-random numbers*. These would fail tests to check if they are random because they often look too 'regular'. One method to generate such numbers is the *Sobol method*.  If we compare pseudo- and quasi-random numbers for estimating $\pi$ using simple rejection sampling, we see that the latter often converges faster. 
"""

# ╔═╡ a8f762eb-e3b6-434a-9526-8f58ae227ae5
md"""
### Convergence of the mean

Finally, let us consider how many samples we would need so that the average is "close" to the expected value. Remember, the spread of a random variable can be quantified using the variance:

$$
\text{Var}[\bar{X}_n] = \frac{1}{n^2}\sum_{i=1}^n \text{Var}[X_i] = \frac{\text{Var}[X]}{n}\,.
$$
Again, tactically assuming finite variances, we see that the variance of the mean scales inversely with the sample size or the standard error on the mean scales with $1/\sqrt{n}$. Note, however, that it does not scale with the dimension of the problem. In practice, this means that, to estimate an expected value to one additional significant digit (meaning reducing the standard error with a factor 10), we must generate a 100 times greater sample size! Looking at the graph of the standard error with sample size, we can draw the following conclusion about sampling methods:

> Sampling methods can easily estimate an expected value up to an order of magnitude below the standard error (tens or hundreds of samples). When large precision is needed, you must generate extreme sample sizes!

What does this mean, practically speaking?
- Suppose you want to numerically compute an definite integral (i.e. an area), it is very difficult to estimate this accurately using stochastic sampling methods so you better use deterministic numerical integration, such as Gaussian quadrature methods. Such methods can get an extremely impressive exponential error decay in the number of samples for smooth functions: $\mathcal{O}(e^{-n})$.
- However, when simulating, we might already be happy with a modest number of observations. For example, suppose $X$ represents the distribution of rainfall duration and intensity and $g(\cdot)$ is a complicated simulator of the hydraulics of the region, including ecohydrology, sewer systems etc. Obviously, there is a lot of variation in $X$ and $g(X)$. If we take a very modest number of independent 12 samples, we have a standard error of 1/4 of the standard error of $g(X)$, meaning we have a 95% confidence interval with a width of a single standard deviation. For these complex problems (and keeping the inherent error in the models into account) this might already suffice to guide management decisions!

More quantitatively, one can bound probabilities using concentration bounds, such as Chebyshev's inequality:

$$P(|X-\mu|\ge k\sigma) \le \frac{1}{k^2}\,,$$

or, specifically for using averages:

$$P(|\bar{X}_n-\mu|\ge k\sigma) \le \frac{1}{nk^2}\,,$$

These bounds are extremely general, though in practice weak to the point of being nearly useless. For example, if we want to compute the probability that throwing a fair coin a 100 times gives exactly 50 times head ($p=0.5\times 100$ in a binomial distribution with $\sigma=\sqrt{100\times(0.5\times 0.5)}=5$) we would find that the probability is less than $1/0.1^2=100$, which is not very helpful to say the least... 

The *central limit theorem*, which states that sample averages of i.i.d. distributed values converge to a normal distribution with corresponding mean and standard deviation, is much more accurate. So, here $X$ is approximately distributed as Normal(50, 5). Here, we can easily estimate the probability $P(49.5<X\le 55.5) \approx 0.04839$. By the way, the true probability of $P(X=50)\approx0.07959$. Statisticians seems to use 30 as the magical sample size where the central limit theory kicks in (a poor man's infinity, if you will). In practice, the convergence is greatly determined by how "bell-shaped" the initial distribution is. """

# ╔═╡ fe24583f-afc7-43ee-901f-c81585677ecd
md"""## Elementary probability distributions

With the above refresher on probability, let us look at some of the most important distributions to describe physical and biological processes. These distributions will serve as building blocks to construct flexible, hierarchical, multivariate distributions. For every distribution, we give the PMF/PDF, the first moment and the parameters with its typical use. There is no need to memorize the exact distribution formullas. Those you frequently use will stick! Rather, see this as a handy reference sheet for modeling. The Julia library `Distributions` contains all of these and more, together with the functions to compute the PDF, CDF, mean, variance and to sample from them. 
"""

# ╔═╡ 2db26aa9-7516-400c-905f-b6b710820667
md"""
### Distributions over integers

| Distribution | PMF                                          | Support                           | Parameters                       | $E[X]$          | Meaning                                                                      | Example                                        |
| ------------ | -------------------------------------------- | --------------------------------- | -------------------------------- | --------------- | ---------------------------------------------------------------------------- | ---------------------------------------------- |
| Bernoulli    | $P(X=k)=p^k(1-p)^{1-k}$                      | $k\in\{0,1\}$                     | $p\in[0,1]$                      | $p$             | Outcome of an experiment that asks a yes-no question.                        | Coin flip.                                     |
| Binomial     | $P(X=k) = \binom{n}{k}\, p^k (1-p)^{n-k}$    | $k \in \{0, 1, ..., n\}$          | $n \in \mathbb{N}, p \in [0, 1]$ | $np$            | Number of successes in a fixed number of independent trials.                 | Number of heads in $n$ coin flips.             |
| Poisson      | $P(X=k) = \frac{e^{-\lambda} \lambda^k}{k!}$ | $k \in \{0, 1, 2, ...\}$          | $\lambda > 0$                    | $\lambda$       | Number of events occurring in a fixed interval of time or space.             | Number of mutations in a long sequence of DNA. |
| Geometric    | $P(X=k) = (1-p)^{k-1}p$                      | $k \in \{1, 2, ...\}$             | $p \in (0,1]$                    | $\frac{1}{p}$   | Number of trials needed until the first success in Bernoulli trials.         | Number of coin flips until the first head.     |
| Uniform      | $P(X=k) = \frac{1}{b-a+1}$                   | $k \in \{a, a+1,\ldots, b-1, b\}$ | $a, b$                           | $\frac{a+b}{2}$ | All outcomes in the integers $\{a, a+1,\ldots, b-1, b\}$ are equally likely. | Rolling a fair six-sided die.                  |

**Example:** In a given period of time, a biologist is interested in counting the number of fish caught in a specific location. The average rate of fish caught per hour is $\lambda = 3$. The biologist can model the number of fish caught in a certain time interval using the Poisson distribution. The parameter $\lambda$ represents the average rate of events (in this case, catching fish) in a fixed time interval. The PMF of the Poisson distribution allows the biologist to calculate the probability of observing a specific number of fish caught within that time interval. This distribution is useful for studying rare events where the average rate is known and where each event is independent of others, such as counting occurrences of diseases in a population or arrivals at a service counter.
"""

# ╔═╡ 8ae351eb-c1b8-4934-aa48-023f8feba75c
md"""
### Distributions over unbounded real numbers

| Distribution | PDF                                                                                                                                                     | Support            | Parameters                       | $E[X]$                                   | Meaning                                                                                   | Example                                                     |
| ------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------ | -------------------------------- | ---------------------------------------- | ----------------------------------------------------------------------------------------- | ----------------------------------------------------------- |
| Normal       | $f_X(x) = \frac{1}{\sqrt{2\pi}\sigma} e^{-\frac{(x-\mu)^2}{2\sigma^2}}$                                                                                 | $x \in \mathbb{R}$ | $\mu \in \mathbb{R}, \sigma > 0$ | $\mu$                                    | Models the distribution of a continuous random variable with symmetric bell-shaped curve. | Yearly rainfall                                             |
| Student's t  | $f_X(x) = \frac{\Gamma\left(\frac{\nu+1}{2}\right)}{\sqrt{\nu\pi}\Gamma\left(\frac{\nu}{2}\right)} \left(1 + \frac{x^2}{\nu}\right)^{-\frac{\nu+1}{2}}$ | $x \in \mathbb{R}$ | $\nu > 0$                        | $0$, for $\nu > 1$, otherwise undefined. | Bell-shaped curve with potential unbounded variance.                                      | Gene expression differences in small groups                 |
| Cauchy       | $f_X(x) = \frac{1}{\pi\gamma(1 + (\frac{x-x_0}{\gamma})^2)}$                                                                                            | $x \in \mathbb{R}$ | $x_0$, $\gamma$                  | Undefined                                | Heavy-tailed distribution, ratio of two normal distributions.                             | Exteme events, such as the annual maximum one-day rainfall  |
| Laplace      | $f_X(x) = \frac{1}{2b} e^{-\frac{x-\mu}{b}}$                                                                                                            | $x \in \mathbb{R}$ | $\mu \in \mathbb{R}, b > 0$      | $\mu$                                    | Double exponential distibution.                                                           | Modelling noise in electronic devices or signal processing. |
The normal or Gaussian distribution is by far the most used probability distribution for unbounded real numbers. The central limit theorem and other justifications are used why many naturally occuring phenomena follow a normal distribution. Indeed, the main bulk of a PDF is often bell-shaped. However, a striking feature of the normal distribution are its very light tail: the log-probability-density decreases quadratically. Encountering events that occurs more than $5\sigma$ from the expected value should only occur with a frequency of $6\times 10^{-7}$, which (as who has invested in the stock market can attest) not very realistic. Extreme events can occur much more frequently than a normal distribution would suggest. Many of these alternative distributions, such as the Laplace distribution have much ticker tails and are more suitable to describe processes with extreme events, such as in climate change.
"""

# ╔═╡ af7119a0-2470-4bf0-8b30-59b38ddd8d5d
md"""
### Distributions over bounded real numbers

| Distribution  | PDF                                                                                                                                                           | Support        | Parameters                                | $E[X]$                         | Meaning                                                                                                  | Example                                                  |
| ------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------- | ----------------------------------------- | ------------------------------ | -------------------------------------------------------------------------------------------------------- | -------------------------------------------------------- |
| Exponential   | $f_X(x) = \lambda e^{-\lambda x}$                                                                                                                             | $x \geq 0$     | $\lambda > 0$                             | $\frac{1}{\lambda}$            | Models the time between independent Poisson events.                                                      | Time between two cell divisions.                         |
| Gamma         | $f_X(x) = \frac{\lambda^k}{\Gamma(k)} x^{k-1} e^{-\lambda x}$                                                                                                 | $x \geq 0$     | $k > 0, \lambda > 0$                      | $\frac{k}{\lambda}$            | Sum of $k$ independent exponentially distributed random variables.                                       | Copy number of constitutively expressed protein.         |
| Inverse Gamma | $f_X(x) = \frac{\lambda^{\alpha}}{\Gamma(\alpha)} \frac{1}{x^{\alpha+1}} e^{-\frac{\lambda}{x}}$                                                              | $x > 0$        | $\alpha > 0, \lambda > 0$                 | $\frac{\lambda}{\alpha-1}$     | Reciprocal of a gamma-distributed random variable.                                                       | Modeling the precision of a normal distribution.         |
| Log-normal    | $f_X(x) = \frac{1}{x \sigma \sqrt{2\pi}} e^{-\frac{(\ln x - \mu)^2}{2\sigma^2}}$                                                                              | $x > 0$        | $\mu \in \mathbb{R}, \sigma > 0$          | $e^{\mu + \frac{\sigma^2}{2}}$ | Distribution of a random variable whose logarithm is normally distributed.                               | Distribution of particle sizes in a powder or a colloid. |
| Uniform       | $f_X(x) = \frac{1}{b-a}$                                                                                                                                      | $x \in [a, b]$ | $a, b \in \mathbb{R}, a < b$              | $\frac{a+b}{2}$                | All outcomes in the interval $[a, b]$ are equally likely.                                                | Generation of random numbers within a range.             |
| Triangular    | $f_X(x) = \begin{cases} \frac{2(x-a)}{(b-a)(c-a)} & \text{for } a \leq x < c \\ \frac{2}{b-a} & \text{for } c \leq x < b \\ 0 & \text{otherwise} \end{cases}$ | $x \in [a, b]$ | $a, b, c \in \mathbb{R}, a \leq c \leq b$ | $\frac{a+b+c}{3}$              | Continuous probability distribution with lower and upper limits.                                         | Estimation of time to complete a task.                   |
| Beta          | $f_X(x) = \frac{x^{\alpha-1}(1-x^{\beta-1}}{B(\alpha,\beta)}$                                                                                                 | $x\in [0,1]$   | $\alpha>0, \beta > 0$                     | $\frac{\alpha}{\alpha+\beta}$  | Distribution over a success probability of a Bernoulli, having seen $\alpha$ success and $\beta$ misses. | Belief of a fraction of mutants.                         |

"""

# ╔═╡ 4290b2c6-1f5f-4780-9f32-08ae949b666c
md"""
### MaxEnt: building distributions

You might wonder how one can come up with sensible probability distributions. One way is using the *maximum entropy principe* (MaxEnt). For any probability function, one can compute the information entropy as:

$$H(X)=-\int_{-\infty}^\infty f_X(x) \log(f_X(x))\,,$$

or,  for discrete probability distributions

$$H(X) = -\sum_ip_i\log(p_i)\,,$$

which measures the average degree of "uncertainty", "information" or "surprise" in the distribution. The MaxEnt method finds distributions that have a maximal entropy, given some fixed properties, such as the support, the mean, etc. For example, the normal distribution is the unique distribution with a mean $\mu$ and a variance $\sigma^2$ with the largest entropy. 

Maximum entropy can be motivated by:
1. It yields the least informative distributions with the largest uncertainty given the data.
2. Nature tends to generate empirical distributions with high entropy.
3. It just works.

Most of the distributions you know and love can be obtained this way. The MaxEnt principle is used a lot in the life sciences. For example, in ecology in species distribution modelling, one often takes the distribution with the largest entropy that fits the data.
"""

# ╔═╡ 07949347-3e7b-499f-9ca6-c62c131e2a65
md"## Combining simple distributions into complex ones"

# ╔═╡ 8e7fa947-480d-4037-88a9-652ed7cc2cbd
md"""
### Joint distributions of independent variables

The easiest way of combining univariate distributions of random variables into a multivariate distribution is by taking an ensemble of independent variables. The joint PDF can then be constructed from the product of the individual variables. So, combining two random variables $X$ and $Y$ the way results in the joint distribution of $X\times Y$:

$$f_{X\times Y}(x,y) = f_X(x)f_Y(y)\,.$$

By definition, conditioning on one variable always leads to the distribution of the other variable because the variables are independent. This follows from the definition:

$$f_{Y\mid X=x}(y) = \frac{f_{X\times Y}(x,y)}{f_X(x)}= \frac{ f_X(x)f_Y(y)}{f_X(x)} = f_Y(y)\,.$$

Using the `Distributions` library, one can create a product distribution of two or more distributions as `product_distribution([dist1, dist2])` (mind the brackets `[]`, the input is a vector of distributions).
"""

# ╔═╡ 1d883714-026b-4206-8ee3-a58e89a70c52
distX = Laplace(8, 6)

# ╔═╡ 52fd8a26-93a1-4c83-a5fd-4247cdff3629
distY = TriangularDist(-2, 2, 1)

# ╔═╡ 432c5cd3-a359-4f80-9a78-d8198b100524
dist_prod = product_distribution([distX, distY])

# ╔═╡ e495c24a-76d5-4cba-b1f9-76d3e25912d8
md" x conditioning : $(@bind xslice Slider(-10:0.05:25, default=16, show_value=true))"

# ╔═╡ 46b4ef15-bd34-4b59-aef5-5dc4e1193ec3
md"y conditioning : $(@bind yslice Slider(-2:0.01:2, default=0, show_value=true))"

# ╔═╡ ffa84522-dd02-47f8-ac3b-b1db08806b2d
md"""
### Mixtures of probability distribitions

A mixture distribution is an elegant way of combining several simple distributions for a variable (or set of variables) into a more complex distribution. Given $k$ PDFs with as as input $x$, $f_{X_1}(x), \ldots f_{X_k}(x)$, we have the mixture density as:

$$f(x) = \sum_{i=1}^kw_if_{X_i}(x)\,,$$

where the $w_i$'s with $w_i\ge 0$ and $\sum_iw_i=1$ are the mixing coefficients. Think of this as picking one of the PDFs randomly with a probability given by the weight according to categorical distribution  (hence, why it was needed to that this could act as a probability vector) and use this one to generate generate a observation according to that distribution. So $w_i$ can be seen as a prior distribution for the components. Mixtures are flexible powerful models! For example, while the humble normal distribution is only a simple unimodal bump, one can combine several normal distributions into a multimodal distribution. For example:
- if you model properties of two populations, such as body weight in a flock of birds containing males and females or canopy size of a group of different species of trees;
- your model might contain a single, sharp peak given by a Gaussian, which you mix with a second distribution with a much larger variance to account for fat tails.

In `Distributions`, one can easily compose a mixture using `MixtureModel([d1, d2, ...], [w1, w2, ...])`.

"""

# ╔═╡ 0252887f-b133-426e-82fc-43ac091a7de6
md"w₁ : $(@bind w1 Slider(0:0.05:1, default=0.2, show_value=true))"

# ╔═╡ 8bfd6e24-27ce-4caf-b68c-b13d89f08337
w = [w1, 1-w1]

# ╔═╡ e79c0f96-02ca-4ebf-bcc5-0f6e0dc3d7d6
d1 = Normal(-2, 1)

# ╔═╡ 5a4a0cb6-93b2-4194-8a65-eb84e396b999
d2 = Normal(3, 2)

# ╔═╡ c580ca85-78ff-4cbc-b033-f0c5acde3193
dist_mixture = MixtureModel([d1, d2], w)

# ╔═╡ f1906adc-79f7-44da-a085-486e8f46e027
md"""
Using Bayes' rule, one can estimate the probability that an observation originates from one of the $k$ distributions:

$$P(\text{sample $x$ originates from component $i$}) = \frac{w_i \times f_{X_i}(x)}{\sum_{i=1}^kw_if_{X_i}(x)}\propto w_i \times f_{X_i}(x)\,.$$
"""

# ╔═╡ d6cd5202-0171-40f9-bc17-d001dadee567
md" x : $(@bind xm Slider(-5:0.1:10, default=0, show_value=true))"

# ╔═╡ 89178521-511e-4c53-a2e1-c1d4dad03a23
md"Mixtures allow for creating complex distributions of more simple variables. We have seen that product distributions are an easy way of combining different random variables $X$, $Y$ into a joint distribution $X\times Y$. They are not very exciting in themselves as there is no dependence between the two variables. However, combining multiple product distributions $X_1\times Y_1$, $X_2\times Y_2$, ... into a mixture, we can introduce some sophisticated dependencies between the variables! This way of combining distributions are sometimes called sum-product networks, for the operations that are used to combine them."

# ╔═╡ 67b154c9-3c14-4eb2-98f2-19f9e3f74d16
begin
	dmv1 = product_distribution([Normal(-2), Normal(-3)])
	dmv2 = product_distribution([Normal(2, 0.7), Normal(3, 2)])
	mvmixture = MixtureModel([dmv1, dmv2], [0.6, 0.4])
end

# ╔═╡ 9381abbc-3715-4f1a-97cf-2487706e13f6
md"Or, we can go very crazy!"

# ╔═╡ 36093484-bdd1-4b1b-bd68-39c7b18a9f6b
spn = MixtureModel(
	[MixtureModel(
		[product_distribution([Normal(10, 1), Normal(20, 1.5)]),
				product_distribution([Normal(20, 1.5), Normal(21, 0.8)]),
				product_distribution([Uniform(5, 25), Uniform(7, 9)]),
				product_distribution([TriangularDist(12, 17), TriangularDist(14, 16, 15.5)]),
				product_distribution([Uniform(2, 28), LogNormal(log(27), 0.1)]),
			],[0.3, 0.25, 0.25, 0.15, 0.05]),
		product_distribution([Uniform(0, 30), Uniform(0, 30)])], [0.99, 0.01]);

# ╔═╡ 1efe9a9e-cd1b-4fe2-8392-2ce42c871a5f
md"""### The multivariate normal distribution

The normal distribution is the most used (and, arguably, misused) probability univariate distribution. Similarly, the multivariate normal distribution (MVN) is undoubtedly the most important multivariate probability distribution. A key property of the normal distribution is that any linear combination of two normally distributed variables follows a normal distribution. The MVN defines random vectors for which arbitrary linear projections follow a normal distribution. The pdf in $k$ dimensions is given by

$$f_{\mathbf {X} }(x_{1},\ldots ,x_{k})={\frac {\exp \left(-{\frac {1}{2}}\left({\mathbf {x} }-{\boldsymbol {\mu }}\right)^{\intercal }{\boldsymbol {\Sigma }}^{-1}\left({\mathbf {x} }-{\boldsymbol {\mu }}\right)\right)}{\sqrt {(2\pi )^{k}|{\boldsymbol {\Sigma }}|}}}$$

with $\boldsymbol{\mu}$ is the vector with the mean and $\boldsymbol{\Sigma}$ the covariance matrix, which has to be symmetric and positive-definite (all eigenvalues must be positive) to be valid.

The MVN can again be motivated from a higher-dimensional case of the Central Limit Theorem or be derived using the MaxEnt principle. The MVN is powerful, as it encodes dependencies between the variables using the covariance matrix. It looks like a well-defined peak that extends in certain directions. Its contours are ellipsis. Intriguingly, the MVN is often both too flexible and not flexible enough:
- For a high number of dimensions $k$, one needs to define $k$ parameters to pinpoint the mean and $k(k-1)/2$ parameters for $\boldsymbol{\Sigma}$. For this reason, the covariance is sometimes considered diagonal (only $k$ parameters) or a scaled identity matrix (one parameter).
- Despite the large number of parameters, the MVN can only capture quite simple patterns: single peaks with only linear dependency between the variables. 

The MVN is key for many, many applied mathematics methods such as statistics, machine learning and control theory, for example, as a building block for the Kalman filter and Gaussian processes. The main reason is that linear transformations of an MVN random vector remain an MVN distributed. Likewise, conditioning is also particularly easy. As this course will predominantly use the sampling approach, we refer to standard works on probability, statistics or machine learning for these formulas. 
"""

# ╔═╡ a1aed96d-8839-4c40-9f18-1ca453be740b
μ = [1, 1]

# ╔═╡ fb6698d6-c565-42f7-ab34-8a4c860f7ef4
md"""
σ₁: $(@bind σ₁ Slider(0.1:0.2:2, show_value=true, default=1))

σ₂: $(@bind σ₂ Slider(0.1:0.2:2, show_value=true, default=2))

ρ: $(@bind ρ Slider(-0.95:0.05:0.95, show_value=true, default=0.9))
"""

# ╔═╡ e0c947b7-306b-42d3-b004-5f37702be73a
Σ = [σ₁^2 σ₁*σ₂*ρ;
	 σ₁*σ₂*ρ σ₂^2]

# ╔═╡ af983208-c64e-4efe-9322-000b2ce09072
md"""
### Bayesian hierarchical modeling

Bayesian hierarchical modeling* is a probabilistic approach approach where the parameters of one probability distribution are themselves treated as random variables and modeled using another distribution. In our notation, this would mean that:

$$\theta \sim f_\theta(\theta)$$

and

$$X\sim f_X(x;\theta)$$

The final joint PDF is given by (see rule for conditional probability):

$$f_{X,\theta}(x, \theta) = f_X(x;\theta)\,f_\theta(\theta)\,.$$

This allows for the incorporation of hierarchical structures, capturing dependencies and variability at multiple levels. This way we create a joint distribution where $X$ depends on $\theta$. It enables the modeling of complex relationships and uncertainties in data by nesting probability distributions within each other. This approach is particularly useful when dealing with hierarchical or structured data, such as grouped data, repeated measurements, or data with varying levels of aggregation. By representing parameters as random variables, Bayesian hierarchical modeling provides a flexible framework for inference and prediction, allowing one to account for uncertainty at different levels of analysis and incorporate prior knowledge effectively.

For example, suppose one orders packages of chili pepper seeds, usually sold in packages of 10 seeds. Assuming these seeds are similar, every seed has a germination probability $p$, and the number of seeds that germinate can be modeled using a Bernoulli distribution $X\sim$Binom(10, $p$). Now, the germination probability $p$ itself depends on the cultivar (hotter peppers typically have a lower chance of germinating). Suppose that over all cultivars, the germination probability $p\sim$Beta(8, 3). We will typically refer to this as a *prior distribution*. You would be free to invent a more complex model for $p$, maybe taking the specific cultivar into account or even considering the conditions in which the seeds are planted. To take a sample of the ensemble $p, X$ we just follow the flow:
1. sample a $p\sim$Beta(8, 3)
2. sample $X\sim$Binom(10, $p$).
You can see the resulting (marginal distribution) for the number of seeds germinated.

Though simple, this model above already illustrates something pretty powerful. We have used our prior knowledge about some process (germination success) and have linked this to a concrete experiment that generates data. Using sampling, we can get an accurate picture of the germination frequency over all possible values of $p$ (without integrating!). It already hints to the inverse problem: given that we sow the ten seeds and observe, for example, eight germinating **can we obtain a better distribution for $p$, given this new data**? We will spend much effort to answer these kinds of questions.

"""

# ╔═╡ e295e21f-40a7-47ae-b680-4831b3d0cdbb
dist_p = Beta(8, 3)

# ╔═╡ cc058abb-e7b6-433e-aacf-13295a401c5f
md"We use a `Turing` model to encode our distribution."

# ╔═╡ f24fffb4-deb7-428a-ba21-0e4999962624
md"How many packages in this large sample have more than 8 seeds germinating?"

# ╔═╡ 95502761-d344-448a-9259-a1a6079a4a8b
md"""
Let us consider a second example. Precipitation intensity and rainfall duration are linked: more intense rainfalls are shorter. Suppose we model the intensity $X$ (in mm/h) and the rainfall duration $Y$ in h as follows using exponential distributions:

$$X\sim Exp(10)$$

$$Y\sim Exp(2/(2+X))$$

We can again sample for the ensemble $X,Y$. By sampling, we can get insight in the total amount of water that pours from the sky during a rainfall (i.e., product of intensity and duration $XY$). (Question: is most of the rain generated by longer, low-intensity rain or short, high-intensity rain)?
"""

# ╔═╡ f3f5595a-e82d-4a9c-b6ab-8dadb553056f
@model function rainfall()
	# we can just use names instead of symbols
	intensity ~ Exponential(10)  
	duration ~ Exponential(2 / (2 + intensity))
end

# ╔═╡ d999582d-5c0e-42df-850f-f7bbe4ccb636
rainfall_samples = [rand(rainfall()) for i in 1:10_000]

# ╔═╡ 4b96292b-f0ad-48ac-99cd-360c0708ef72
rainfall_amount = prod.(rainfall_samples)

# ╔═╡ cb903ece-6f36-4840-ad53-1b095f189070
md"""
Many textbooks describe such hierarchical models a *graphical models*. These are networks where the nodes represent distributions and the links the connections between them. There are many variants, such as Markov random fields (undirected) or (directed) *Bayesian networks*. We focus only on the latter, where distributions are placed in a directed acyclic graph as seen below. 

![](https://upload.wikimedia.org/wikipedia/commons/thumb/e/e2/Example_of_a_Directed_Graph.svg/440px-Example_of_a_Directed_Graph.svg.png)

You will find it likely not too tricky to sample from Bayesian networks by sampling from the 'root' distributions and working yourself through the graph. We will however use a systematic approach were we will encode our hierarchical distribution in a called "probabilistic program". Probabilistic programming allows one to implement complex hierarchical probability distributions and has convenient methods for sampling from these (the topic of this chapter) and perform inference (the topic of next chapter). There are many probabilistic programming languages. We will use Turing, which is a domain-specific language in Julia.
"""

# ╔═╡ 0681ddf3-0aba-42e2-831d-c5058c417009
md"Here is a Turing version to estimate $\pi$:"

# ╔═╡ 875cb743-7c6e-4aed-893d-03b969dc4c35
@model function circle_throw()
	# generate a random point in the [-1, 1] × [-1, 1] square
	x ~ Uniform(-1, 1)
	y ~ Uniform(-1, 1)
	# check if (x, y) is in the circle
	in_circle ~ Dirac(x^2 + y^2 ≤ 1)  # Dirac stores a value
end

# ╔═╡ c12e165f-14b1-4df3-bc88-28ada3bf5add
rand(circle_throw())

# ╔═╡ b76245bf-a5a4-41f4-a20b-0327280e68e6
4mean([rand(circle_throw())[:in_circle] for i in 1:10_000])

# ╔═╡ 9197faa7-603f-4626-b886-a6f71e4111f0
md"Note that we have used the special `Dirac` distribution, which has a probability mass of 1 when its input is true and zero elsewhere. It is a convenient way of sotoring a value."

# ╔═╡ 43305b47-9e34-40dd-aa5f-fb54410fc45e
md"""
## Putting it all together

With what we have seen up to now, we approximately answer many interesting questions of probability. To summarize what we have seen:
1. Using the ideas of Bayesian hierarchical modeling, mixture modeling and the like, we can easily create a probabilistic program of complex, joint distributions.
2. Simply by using the `rand` function, we can generate a large number of samples from this distribution.
3. Most questions that we can think of, can be phrased as expected values, which the LLN allows us to estimate by averaging from a large sample.

Let us revise the pepper seed germination example. What if we want to know the probability that more than eight seed germinate? Using conventional probability, this would boils down to integrating over all values of $p$

$$P(X\ge 9) = \sum_{k=9}^{10}\int_0^1f_X(k;p)f_p(p)\mathrm{d}p\,,$$
which does not look very fun to compute! However, suppose we generate 10,000 samples from this model (so over all possible values of $p$), we can estimate this probability as the fraction in which nine or more seed germinated! We can use the the function `count` to estimate this probability: `count(pX->pX[:X]>8, seed_sample)`, which evaluates (in this run) to 2912, so $P(X>8)\approx 0.2912$.  An alternative way would be using the `mean` function, where `mean(g, X)` is equivalent to `mean(g.(X))`. So the same result is obtained using `mean(pX->pX[:X]>8, seed_sample)`. A third way is to use the function `filter`, which works the same way, but just keeps the subset of samples that satisfy the indicator function. So the length of `filter(pX->pX[:X]>8, seed_sample)` retains only samples where $X>8$.

The `filter` function can be used to create conditional distributions. For example, can we obtain the distribution of $p$ given that $X=8$? This subsample can be obtained using `filter(pX->pX[:X]==8, seed_sample)`. We can plot the histogram of `p` in this subsample to view this distribution! You might find it easier to use this list comprehension to extract these conditional values:
"""

# ╔═╡ b1e85d42-199e-4874-843c-f5aa81728c38
md"$(@bind ngerm Slider(0:10, default=8))"

# ╔═╡ 6d8c2a59-baf3-45e4-8811-adc7857e1732
ngerm

# ╔═╡ a0e20b67-3d84-44cf-b12e-8c0b2edba982
md"""
In this syntax, you can iterate over all values in `seed_sample` and chose which value to retain under what condition! We can compute the conditional mean and variance of `p_cond`, compute the posterior likelihood that $p$ is greater than some value and much more!

Let us systematically go over the main probabilistic quantities you can compute using a sample $\{x_1,\ldots, x_n\}$ drawn i.i.d. from $X$:
- **Expected values**: Let $g$ be any function, we can estimate $E[g(X)]$ as $\frac{1}{n}\sum_{i=1}^ng(x_i)$, or, in Julia either as `mean(g, Xs)` or, equivalently `mean(g.(Xs))`. For example, the average number of seeds that  grow from a package of 10.
- **Probabilities**: In the most general form, we can an event using a function $b(x)$, which either returns 1 if $x$ satisfies the condition of the event or 0 otherwise. This can be computed using expected vales using `mean(b, Xs)`, which computes $P(b(X)=1)=E[b(X)]$. For example, to compute the probability $P(5 < X < 8)$ one would use `mean(x->5 < x < 8, Xs)` or `mean(5 .< Xs .< 8)`.
- **Conditional expectations** of the form $E[g(x)\mid b(x)=1]$, which can be attained using `filter` and `mean`: `mean(g, filter(b, Xs))`. 
"""

# ╔═╡ 1ad364c5-b65d-48d7-a829-e284c5e8eebc
md"""
## Appendix 🐉
"""

# ╔═╡ 5614c4b1-2c78-4fee-a480-3998a8e1ce6e
TableOfContents()

# ╔═╡ efd7bcba-f393-4659-ac5b-006740085144
cummean(x) = cumsum(x) ./ (1:length(x))

# ╔═╡ 0b9060d0-9c5a-4806-8d68-20e4ffd8f864
dist2pdf(distribution) = x -> pdf(distribution, x)

# ╔═╡ 7d71303b-defc-4e5f-beee-2bf684674ecf
dist2cdf(distribution) = x -> cdf(distribution, x)

# ╔═╡ 57052ba3-151c-4dda-904a-dac9c421241c
plots = Dict()  # storing all the figures

# ╔═╡ 41967e5b-9068-4318-acb6-d958b608c8f7
let
	p = plot(x->sqrt(1-x^2), 0, 1, aspect_ratio=:equal, xlim=[0,1], ylim=[0, 1], lw=2, fillrange=zero, fillalpha=0.3, xlab=L"x", ylab=L"y", label="circle", title="Monte Carlo for estimating π/4")
	scatter!(x_unif[in_circle], y_unif[in_circle], ms=0.8, markercolor="orange", markerstrokewidth=0, label="in circle")
	scatter!(x_unif[.!in_circle], y_unif[.!in_circle], ms=0.8, markercolor="green", markerstrokewidth=0, label="out of circle")
	plots["pi_sample"] = p
	p
end

# ╔═╡ b9be819e-2db9-4952-9453-98bf815b00e8
let
	pfx = plot(f_X, 0, 20, label="pdf", xlab=L"x", lw=2, ylab=L"f_X(x)")
	plot!(pfx, zero, 0, xpdf, fillrange=f_X, fillalpha=0.5, label="P(X≤$xpdf)", title="Probability density function of Weibull distributon")
	pFx = plot(F_X, 0, 20, label="cdf", lw=2)
	scatter!(pFx, [xpdf], [F_X(xpdf)], label="P(X≤$xpdf)", xlab=L"x", ylab=L"F_X(x)",title="Cumulative density function of Weibull distributon")
	p = plot(pfx, pFx, layout=(2,1))
	plots["pdf_cdf"] = p
	plots["pdf"] = pfx
	plots["cdf"] = pFx
	p
end

# ╔═╡ 7779254a-83b9-4849-957c-af4b62c3916a
let
	ps = []
	for n in [10, 100, 1000, 10_000]
		p = histogram(rand(dist, n), label="", title="n=$n", normalize=true)
		plot!(x -> pdf(dist, x), 0, 20, label="", lw=2)
		push!(ps, p)
		plots["lln_$(n)"] = p
	end
	plots["lln_vert"] = plot(ps..., layout=(4,:))
	plot(ps...)
end

# ╔═╡ 73898b6c-d8bd-49a9-ba36-a068e62a7ec3
let
	d = Binomial(30, 0.3)
	#d = dist
	p = hline([mean(d)], linestyle=:dash, lw=2, color=:red, label=L"E[X]", xlab=L"n", ylabel=L"\bar{X}_n", title="The Law of Large Numbers")
	for i in 1:5
		plot!(cummean(rand(d, 250)), lw=2, label="", alpha=0.6)
	end
	#ylims!(8, 10)
	plots["llm_conv"] = p
	p
end

# ╔═╡ d7e88a1a-73e9-48ff-bb51-2292346f88e0
let
	p = scatter(rand(200), rand(200), label="", title="200 pseudo-random numbers", aspect_ratio=:equal, xlims=[0,1])
	plots["PRN"] = p
end

# ╔═╡ 74e50a3d-7623-45c2-8d3a-2c34f363177e
begin
	s = SobolSeq(2)
	p = reduce(hcat, next!(s) for i = 1:200)'
	pl = scatter(p[:,1], p[:,2], label="", title="200 Sobol quasi-random numbers", aspect_ratio=:equal, xlims=[0,1])
	plots["QRN"] = pl
end

# ╔═╡ a3c7cd83-4488-40a1-8f0f-b79be73afb44
@model function seed_germ(n=10)
	p ~ Beta(8, 3)
	X ~ Binomial(n, p)
end

# ╔═╡ 8f13d578-f400-4885-a071-a54259b93f12
Xp = rand(seed_germ())  # sample from the seed distribution

# ╔═╡ 788193d3-742b-4c27-aad4-41a9024ff0d0
Xp[:p]  # value of p

# ╔═╡ dc668eb0-899a-42dd-901b-bbbde5d3d0dd
Xp[:X]  # value of X

# ╔═╡ 03cedb9e-f74e-4662-abc8-e3cc87cad581
seed_sample = [rand(seed_germ()) for i in 1:10_000]  # many samples

# ╔═╡ c6e08f0d-e30e-4d30-989e-2463289ee038
germination_counts = [Xp[:X] for Xp in seed_sample]  # extract X

# ╔═╡ 42867e88-cd02-4634-a225-a9bf0ae46b70
count(>(8), germination_counts)

# ╔═╡ 13cf931e-5333-46ff-a00a-ea94c0b2e509
mean(>(8), germination_counts)  # P(X > 8)

# ╔═╡ 6d565e90-7aa6-480f-96fb-25c745ec6974
count(((p, X),)->X>8, seed_sample)

# ╔═╡ 79be0044-9856-4b3b-80e3-6c362e23b9f4
p_cond = [p for (p, X) in seed_sample if X==ngerm]

# ╔═╡ 18e788ca-5b74-4a0b-81f0-5f307857e7e2
mean(>(0.9), p_cond)  # probability of p > 0.9 given x

# ╔═╡ 037a2175-22b2-48c9-9664-3743b11e8f7c
mean(p_cond), var(p_cond)  # mean and variance of the conditional distr. of p

# ╔═╡ 6fe05275-f1f7-4b3f-bbc2-8077bf981391
length(p_cond)

# ╔═╡ 67594833-14c8-423c-94ce-91ce2a44b6c9
let
	Random.seed!(12)
	psob = reduce(hcat, next!(s) for i = 1:1024)'
	ppseu = rand(1024, 2)
	
	pi_sob = pi .- 4cummean(norm.(eachrow(psob)) .≤ 1) .|> abs
	pi_pseu = pi .- 4cummean(norm.(eachrow(ppseu)) .≤ 1) .|> abs
	
	p = plot(pi_pseu, yscale=:log10, label="pseudo-random", xlab=L"n", ylab=L"|E[X]-\bar{X}_n|", lw=2)
	plot!(pi_sob, label="quasi-random", lw=2)
	title!("Error estimating π using sampling")
	plots["pi_sampling_conv"] = p
end

# ╔═╡ e3d9a5b7-73ec-4c16-8799-2e2a4f99c018
let
	p = plot(n->1/√(n), 1, 100, lw=2, label=L"1/\sqrt{n}", ylab=L"\sigma_{\bar{X}_n}", xlabel=L"n", title="Standard error of the mean")
	plots["ste_mean"] = p
end

# ╔═╡ e94fa51a-e1a7-4a1b-88ef-76b38d3a6d74
let
	p = plot(n->1/√(n), 1, 1000_000, lw=2, yscale=:log10, xscale=:log10, label=L"1/\sqrt{n}", ylab=L"\sigma_{\bar{X}_n}", xlabel=L"n", title="Standard error of the mean (log-scale)")
	plots["ste_mean_log"] = p
end

# ╔═╡ 222cdcae-5b51-4fdc-a794-11471449ebcd
let
# Binomial
	pbin = scatter(k->pdf(Binomial(20, .5), k), 0:25, label="n=25, p=0.5", xlab=L"k", ylab=L"P(X=k)")
	scatter!(pbin, k->pdf(Binomial(20, .25), k), 0:25, label="n=25, p=0.25", marker=:^)
	title!("Binomial distribution")
	plots["binomial_distr"] = pbin
	
	# Poisson
	ppois = scatter(k->pdf(Poisson(1), k), 0:20, label="λ=1", xlab=L"k", ylab=L"P(X=k)")
	scatter!(ppois, k->pdf(Poisson(5), k), 0:20, label="λ=5", marker=:^)
	scatter!(ppois, k->pdf(Poisson(10), k), 0:20, label="λ=10", marker=:v)
	title!("Poisson distribution")
	plots["poisson_distr"] = ppois
	
	# Geometric
	pgeo = scatter(k->pdf(Geometric(0.1), k), 1:25, label="p=0.1", xlab=L"k", ylab=L"P(X=k)")
	scatter!(pgeo, k->pdf(Geometric(0.2), k), 0:25, label="p=0.2", marker=:^)
	scatter!(pgeo, k->pdf(Geometric(0.05), k), 0:25, label="p=0.05", marker=:v)
	title!("Geometric distribution")
	plots["geom_distr"] = pgeo
	
	# Geometric
	punif = scatter(k->pdf(DiscreteUniform(5, 10), k), 0:30, label="a=5, b=10", xlab=L"k", ylab=L"P(X=k)")
	scatter!(punif, k->pdf(DiscreteUniform(10, 25), k), 0:30, label="a=10, b=25", marker=:^)
	title!("Uniform distribution")
	plots["unif_discr_distr"] = punif
	
	plot(pbin, ppois, pgeo, punif)
end

# ╔═╡ ebe190e6-da7a-4c88-bd31-b5dc53f6a3e1
let
	# normal
	pnorm = plot(xlabel=L"x", ylabel=L"f_X(x)", title="Normal distribution")
	for (μ, σ) in [(0, 1), (2, 1), (0, 2), (-1, 2)]
		plot!(x->pdf(Normal(μ, σ),x),-5, 5, lw=2, label="N($μ, $σ)", ls=:auto)
	end
	plots["normal_distr"] = pnorm




end

# ╔═╡ 3de53132-5b89-483a-872e-4cf1142fc90b
let
		# Laplace
	plaplace = plot(xlabel=L"x", ylabel=L"f_X(x)", title="Laplace distribution")
	for (μ, θ) in [(0, 1), (2, 1), (0, 2), (-1, 2)]
		plot!(x->pdf(Laplace(μ, θ),x),-8, 8, lw=2, label="Laplace($μ, $θ)", ls=:auto)
	end
	plaplace

	plots["laplace_distr"] = plaplace
end

# ╔═╡ 4360ca84-c0bd-4c7a-afbc-edfef5a5f976
let
	dist_trian = TriangularDist(-2, 2) 
	dist_trian2 = TriangularDist(-2, 2, 1)  # non-symmetric

	p = plot(dist2pdf(dist_trian), -3, 3, xlab=L"x", label="Triang(-2, 2)", lw=2,ylab=L"f_X(x)", title="Triangular distribution")
	plot!(dist2pdf(dist_trian2), -3, 3, label="Triang(-2, 2, 1)", lw=2, ls=:auto)
	plot!(dist2pdf(TriangularDist(-3, 1, -2)), -3, 3, label="Triang(-3, 1, -2)", lw=2, ls=:auto)
	plots["triangular_distr"] = p
end

# ╔═╡ dbdf1da5-67fd-47ea-b363-f03e1e0cdeb0
let
	p = plot(dist2pdf(Exponential(1)), 0, 8, xlab=L"x", label="Exp(1)", lw=2,ylab=L"f_X(x)", title="Exponential distribution")
	plot!(dist2pdf(Exponential(2)), 0, 8, label="Exp(2)", lw=2, ls=:dash)
	plot!(dist2pdf(Exponential(1/2)), 0, 8, label="Exp(1/2)", lw=2, ls=:dot)
	plots["exp_distr"] = p
end

# ╔═╡ 0b3e16bd-032b-4630-8cb3-57c63bc44572
let
	p = plot(dist2cdf(Exponential(1)), 0, 8, xlab=L"x", label="Exp(1)", lw=2,ylab=L"F_X(x)", title="Exponential distribution (CDF)")
	plot!(dist2cdf(Exponential(2)), 0, 8, label="Exp(2)", lw=2, ls=:dash)
	plot!(dist2cdf(Exponential(1/2)), 0, 8, label="Exp(1/2)", lw=2, ls=:dot)
	plots["exp_distr_CDF"] = p
end

# ╔═╡ 05152721-dbb5-4e26-86ba-2fc16c03dac8
let
	xm = 4
	p = plot(dist2pdf(LogNormal()), 0, xm, xlab=L"x", label="Log-normal(0, 1)", lw=2,ylab=L"f_X(x)", title="Log-normal distribution")
	plot!(dist2pdf(LogNormal(1, 1)), 0, xm, label="log-normal(1, 1)", lw=2, ls=:auto)
	plot!(dist2pdf(LogNormal(0, 2)), 0, xm, label="log-normal(0, 2)", lw=2, ls=:auto)
	plot!(dist2pdf(LogNormal(1, 2)), 0, xm, label="log-normal(1, 2)", lw=2, ls=:auto)
	plots["lognorm_distr"] = p
end

# ╔═╡ adbb17b4-56d2-49ca-8bdb-ecd2d3f418cd
let
	p = plot(dist2pdf(Beta(1, 1)), 0, 1, xlab=L"x", label="Beta(1, 1)", lw=2,ylab=L"f_X(x)", title="Beta distribution")
	plot!(dist2pdf(Beta(4, 4)), 0, 1, label="Beta(4, 4)", lw=2, ls=:auto)
	plot!(dist2pdf(Beta(8, 2)), 0, 1, label="Beta(8, 2)", lw=2, ls=:auto)
	plot!(dist2pdf(Beta(1, 3)), 0, 1, label="Beta(1, 3)", lw=2, ls=:auto)
	plots["beta_distr"] = p
end

# ╔═╡ f2faec3d-3d44-43e9-bf03-37ba68e37300
let
	p = contourf(-10:0.05:25, -2:0.01:2, (x,y)->pdf(dist_prod, [x,y]), color=:speed, xlab=L"x", ylab=L"y")
	vline!([xslice], label="Y | X=$xslice", color="orange", lw=2)
	hline!([yslice], label="X | Y=$yslice", color="blue", ls=:dash, lw=2)
	title!("Joint PDF of a Laplace and Triangular distribution")
	plots["prod_distr"] = p
end

# ╔═╡ 4bd5e387-ac06-4f4b-822f-cbc89e7509f9
plots["prod_cond_y"] = plot(y->pdf(distY, y), -2:0.01:2, label="Y | X=$xslice", color="orange", lw=2, xlabel="y", )

# ╔═╡ c7b2fc1e-9dd3-42e9-aa55-f9aea4781784
plots["prod_cond_x"] = plot(x->pdf(distX, x), -10:0.05:25, label="f_{X | Y=$yslice}", color="blue", ls=:dash, lw=2, xlabel="x", )

# ╔═╡ 94dd1214-e63c-44d5-b5c3-9ca834a89ebf
let
	p = plot(x->pdf(dist_mixture, x), -5, 10, label="mixture", lw=2, xlab=L"x")
	plot!(x->pdf(d1, x), -5, 10, label="component 1", ls=:dash)
	plot!(x->pdf(d2, x), -5, 10, label="component 2", ls=:dash)
	plots["mixture"] = p
end

# ╔═╡ fc3a7a70-89e5-43e4-8dae-b7d2d8b54355
let
	p = plot(x->pdf(dist_mixture, x), -5, 10, label="mixture", lw=2, xlab=L"x", ylab=L"f_X(x)")
	plot!(x->pdf(d1, x), -5, 10, label="component 1", ls=:dash)
	plot!(x->pdf(d2, x), -5, 10, label="component 2", ls=:dash)
	vline!([xm], color=:gold, lw=2, label="x=$xm")
	plots["cond_mixture"] = p
end

# ╔═╡ e61e45a0-d9ec-4ab6-8c1b-fc8579a9a9b3
let
	likelihood = [pdf(d1, xm), pdf(d2, xm)]
	posterior = likelihood .* w ./ pdf(dist_mixture, xm)
	plots["mixture_likelihood"] = bar(posterior, xticks=1:length(w), xlabel="component", ylabel="posterior probability", label="X=$xm", color=:gold)
end

# ╔═╡ ebf13b85-7fa0-4859-88d1-a5189c5bbd40
let

	plots["gaussian_mixture"] = contourf(-5:0.1:5, -5:0.1:8, (x,y)->pdf(mvmixture, [x,y]), color=:speed, xlab=L"x", ylab=L"y", title="mixture of two MVN")
end

# ╔═╡ 40fbbc43-4a73-4d25-a37c-b39798d60761
plots["spn"] = contourf(0:0.1:30, 0:0.1:30, (x,y)->logpdf(spn, [x,y]), color=:speed, xlab=L"x", ylab=L"y")

# ╔═╡ ccbcce90-c46f-4f61-8907-2ccc3e692184
let
	mvn = MultivariateNormal(μ, Σ)
	p = contourf(-5:0.05:5, -5:0.05:5, (x,y)->pdf(mvn, [x,y]), color=:speed, xlab=L"x", ylab=L"y")
	title!("PDF of a MVN")
	plots["mvn"] = p
end

# ╔═╡ 6602c372-158d-4b40-b7fe-aea7e534bd39
plots["beta"] = plot(p->pdf(dist_p, p), 0, 1, lw=2, label="Beta(8,3)", xlab=L"x", ylab=L"f_p(p)")

# ╔═╡ 77664b93-4b90-43b3-a0cc-0b09c7bd9854
plots["seeds"] = histogram(germination_counts, xticks=0:10,
		xlabel=L"k", ylabel="frequency in sample",
		label="germination succes X", title="Numer of pepper seeds that germinated")

# ╔═╡ cda5faf1-8d74-4679-a335-7a5974010f5a
plots["rainfall"] = scatter(first.(rainfall_samples), last.(rainfall_samples), alpha=0.2, xlab="intensity (mm/h)", ylab="duration (h)", label="")

# ╔═╡ f7ca0af5-dcdf-4e75-9063-2d194738be64
plots["rain_cond"] = histogram(rainfall_amount, xlab="rainfall amount (mm)", ylab="frequency", label="intensity * duration")

# ╔═╡ cbb71827-67cd-43aa-8a36-757bdf98e16a
plots["seeds_cond"] = histogram(p_cond, ylabel="frequency", xlab=L"p", label="p given n=$ngerm", xlims=[0,1])

# ╔═╡ a9a8f40f-da04-475b-aeaa-07add585ad9c
plots # dictionary of all the figures, for saving

# ╔═╡ Cell order:
# ╠═2d346892-cb15-11ee-2d81-73e08fcc3288
# ╠═cf8dbe5e-d340-4ebb-8ad7-f483be07fffd
# ╠═4d96d7ea-1216-46d6-80fe-5859640eb9c9
# ╟─cfc49c56-f223-4b78-b6da-4e0c83a16bbf
# ╟─c469b254-ffc4-458d-8676-db135e11b306
# ╠═5e81b02d-a879-4023-a7fe-e6f187a558f1
# ╟─9b3f076f-e50b-4f23-a282-91fac47de21b
# ╠═5db3ed2d-ee29-4cd5-b438-b60d550a6ff8
# ╠═7d1e1a93-510c-4f81-bb09-d87e3b40e136
# ╟─41967e5b-9068-4318-acb6-d958b608c8f7
# ╠═5821aaf1-744d-452a-96d2-6516c7423fff
# ╠═cb36a752-da05-4da8-8dc8-1c2d34a11ad5
# ╠═7c244dc0-6c7b-451e-b7c0-da25a77da657
# ╠═fe471a6a-5c29-49f6-afba-83be7ab5d770
# ╟─e0f96ed2-bc51-467e-be3e-8617293e33e2
# ╠═c0677edb-29a4-4c4d-8ca7-9d81eac05c85
# ╠═d786da44-5cfc-4c36-ae41-af6f623d5a4e
# ╠═bce61d18-767d-4ae4-a211-759ef3e4da44
# ╠═1f9b17a8-d8f1-4256-bebc-d5f92e085fe8
# ╟─fce1c003-6483-44f1-bdbe-c6363d09807b
# ╠═dd8c1370-5fb9-4ab7-8729-ea756762788d
# ╠═1882bc65-abca-45b7-8f53-e76360b8dacd
# ╠═d04ea704-c9c5-4d57-8908-a7b19b5567a0
# ╟─b9be819e-2db9-4952-9453-98bf815b00e8
# ╟─2538dd4a-96a3-4d4f-b03e-82d2e8b9ea1f
# ╟─fff75ceb-b21c-49fe-b1c0-0d2f3213b37d
# ╟─6ae427c7-a65e-4f9b-b5a7-8cb523fcd382
# ╠═090d0c6f-de9f-4ce8-befa-4ae9dbe1b38e
# ╟─70cfc225-dd83-472c-b745-5266026c3bfc
# ╟─85e08b91-10f9-4808-ac53-195f90387688
# ╠═36b475e9-a16b-464f-a423-a4fbfbf21ef0
# ╟─f1ee69eb-d8e5-4a6b-b408-a52d33a09e52
# ╠═ab4afb22-28e1-4d9e-8d89-a1a9f0741195
# ╠═4d6e67d1-42b8-481f-b9fe-3df524d4b07b
# ╟─2a30c7eb-de92-46a1-ac79-b15529e8f463
# ╠═83db8998-125e-4cae-b6fe-1cdddfe3f868
# ╠═9f525fef-80b3-4417-85ea-e318ad387985
# ╟─d249188a-a298-44b8-8dcb-691753331439
# ╟─3b9e8f64-d307-48bc-9a71-1ab9accd134e
# ╠═1a991ac1-d4ed-4774-a9e2-a5f27f99e1ff
# ╠═b5b55083-2fe3-4620-9759-8a2329154143
# ╟─f625d656-35de-47cf-985b-3ae9fa54c5d4
# ╠═482d1da0-3415-425c-84de-69f7fb821934
# ╠═7008afc0-785e-45a8-9b6d-5b8f26bae09d
# ╠═64e0265f-a621-495f-b73a-0af23b0ef1ac
# ╠═18a1b129-eac5-4080-8af9-7d9fbf800107
# ╠═5a1522f7-411d-40f5-92e8-b5a154b8d735
# ╟─7779254a-83b9-4849-957c-af4b62c3916a
# ╟─73898b6c-d8bd-49a9-ba36-a068e62a7ec3
# ╟─8a493f5e-1f1b-4cf5-91e1-5504593cee9d
# ╟─dc7eb2b3-b160-4943-84e9-8cf20d687d93
# ╠═2d3e09ff-3766-46d5-abdd-4132c4fa06b9
# ╠═3ecfb175-0c26-419e-9770-d406cd1dd2dc
# ╟─9ba45ba1-f2d3-4d08-a3e7-2ed8730e3d90
# ╠═e67e0b43-1346-4d7c-8bcc-ca53949412b8
# ╠═ac78fddb-0b24-4956-8137-44615721bc06
# ╠═44ca88b0-4256-4123-93f8-82d1ba3b3b77
# ╟─896da840-ad94-40d7-8a65-be92b174e88b
# ╟─d7e88a1a-73e9-48ff-bb51-2292346f88e0
# ╠═2b79770c-00a8-4c3b-b39e-b31400593ea4
# ╟─74e50a3d-7623-45c2-8d3a-2c34f363177e
# ╟─67594833-14c8-423c-94ce-91ce2a44b6c9
# ╟─a8f762eb-e3b6-434a-9526-8f58ae227ae5
# ╟─e3d9a5b7-73ec-4c16-8799-2e2a4f99c018
# ╟─e94fa51a-e1a7-4a1b-88ef-76b38d3a6d74
# ╟─fe24583f-afc7-43ee-901f-c81585677ecd
# ╟─2db26aa9-7516-400c-905f-b6b710820667
# ╟─222cdcae-5b51-4fdc-a794-11471449ebcd
# ╟─8ae351eb-c1b8-4934-aa48-023f8feba75c
# ╟─ebe190e6-da7a-4c88-bd31-b5dc53f6a3e1
# ╟─3de53132-5b89-483a-872e-4cf1142fc90b
# ╟─4360ca84-c0bd-4c7a-afbc-edfef5a5f976
# ╟─af7119a0-2470-4bf0-8b30-59b38ddd8d5d
# ╟─dbdf1da5-67fd-47ea-b363-f03e1e0cdeb0
# ╟─0b3e16bd-032b-4630-8cb3-57c63bc44572
# ╟─05152721-dbb5-4e26-86ba-2fc16c03dac8
# ╟─adbb17b4-56d2-49ca-8bdb-ecd2d3f418cd
# ╟─4290b2c6-1f5f-4780-9f32-08ae949b666c
# ╟─07949347-3e7b-499f-9ca6-c62c131e2a65
# ╟─8e7fa947-480d-4037-88a9-652ed7cc2cbd
# ╠═1d883714-026b-4206-8ee3-a58e89a70c52
# ╠═52fd8a26-93a1-4c83-a5fd-4247cdff3629
# ╠═432c5cd3-a359-4f80-9a78-d8198b100524
# ╟─f2faec3d-3d44-43e9-bf03-37ba68e37300
# ╟─e495c24a-76d5-4cba-b1f9-76d3e25912d8
# ╟─4bd5e387-ac06-4f4b-822f-cbc89e7509f9
# ╟─46b4ef15-bd34-4b59-aef5-5dc4e1193ec3
# ╟─c7b2fc1e-9dd3-42e9-aa55-f9aea4781784
# ╟─ffa84522-dd02-47f8-ac3b-b1db08806b2d
# ╟─0252887f-b133-426e-82fc-43ac091a7de6
# ╟─8bfd6e24-27ce-4caf-b68c-b13d89f08337
# ╠═e79c0f96-02ca-4ebf-bcc5-0f6e0dc3d7d6
# ╠═5a4a0cb6-93b2-4194-8a65-eb84e396b999
# ╠═c580ca85-78ff-4cbc-b033-f0c5acde3193
# ╟─94dd1214-e63c-44d5-b5c3-9ca834a89ebf
# ╟─f1906adc-79f7-44da-a085-486e8f46e027
# ╟─d6cd5202-0171-40f9-bc17-d001dadee567
# ╟─fc3a7a70-89e5-43e4-8dae-b7d2d8b54355
# ╟─e61e45a0-d9ec-4ab6-8c1b-fc8579a9a9b3
# ╟─89178521-511e-4c53-a2e1-c1d4dad03a23
# ╠═67b154c9-3c14-4eb2-98f2-19f9e3f74d16
# ╟─ebf13b85-7fa0-4859-88d1-a5189c5bbd40
# ╟─9381abbc-3715-4f1a-97cf-2487706e13f6
# ╠═36093484-bdd1-4b1b-bd68-39c7b18a9f6b
# ╟─40fbbc43-4a73-4d25-a37c-b39798d60761
# ╟─1efe9a9e-cd1b-4fe2-8392-2ce42c871a5f
# ╠═a1aed96d-8839-4c40-9f18-1ca453be740b
# ╟─fb6698d6-c565-42f7-ab34-8a4c860f7ef4
# ╟─e0c947b7-306b-42d3-b004-5f37702be73a
# ╟─ccbcce90-c46f-4f61-8907-2ccc3e692184
# ╟─af983208-c64e-4efe-9322-000b2ce09072
# ╠═ba4460b1-2855-4695-a3d8-4976666ef164
# ╠═e295e21f-40a7-47ae-b680-4831b3d0cdbb
# ╟─6602c372-158d-4b40-b7fe-aea7e534bd39
# ╟─cc058abb-e7b6-433e-aacf-13295a401c5f
# ╠═a3c7cd83-4488-40a1-8f0f-b79be73afb44
# ╠═8f13d578-f400-4885-a071-a54259b93f12
# ╠═788193d3-742b-4c27-aad4-41a9024ff0d0
# ╠═dc668eb0-899a-42dd-901b-bbbde5d3d0dd
# ╠═03cedb9e-f74e-4662-abc8-e3cc87cad581
# ╠═c6e08f0d-e30e-4d30-989e-2463289ee038
# ╟─77664b93-4b90-43b3-a0cc-0b09c7bd9854
# ╟─f24fffb4-deb7-428a-ba21-0e4999962624
# ╠═6d565e90-7aa6-480f-96fb-25c745ec6974
# ╠═42867e88-cd02-4634-a225-a9bf0ae46b70
# ╠═13cf931e-5333-46ff-a00a-ea94c0b2e509
# ╠═18e788ca-5b74-4a0b-81f0-5f307857e7e2
# ╟─95502761-d344-448a-9259-a1a6079a4a8b
# ╠═f3f5595a-e82d-4a9c-b6ab-8dadb553056f
# ╠═d999582d-5c0e-42df-850f-f7bbe4ccb636
# ╟─cda5faf1-8d74-4679-a335-7a5974010f5a
# ╠═4b96292b-f0ad-48ac-99cd-360c0708ef72
# ╟─f7ca0af5-dcdf-4e75-9063-2d194738be64
# ╟─cb903ece-6f36-4840-ad53-1b095f189070
# ╟─0681ddf3-0aba-42e2-831d-c5058c417009
# ╠═875cb743-7c6e-4aed-893d-03b969dc4c35
# ╠═c12e165f-14b1-4df3-bc88-28ada3bf5add
# ╠═b76245bf-a5a4-41f4-a20b-0327280e68e6
# ╟─9197faa7-603f-4626-b886-a6f71e4111f0
# ╟─43305b47-9e34-40dd-aa5f-fb54410fc45e
# ╟─b1e85d42-199e-4874-843c-f5aa81728c38
# ╠═6d8c2a59-baf3-45e4-8811-adc7857e1732
# ╠═79be0044-9856-4b3b-80e3-6c362e23b9f4
# ╠═037a2175-22b2-48c9-9664-3743b11e8f7c
# ╠═6fe05275-f1f7-4b3f-bbc2-8077bf981391
# ╟─cbb71827-67cd-43aa-8a36-757bdf98e16a
# ╟─a0e20b67-3d84-44cf-b12e-8c0b2edba982
# ╟─1ad364c5-b65d-48d7-a829-e284c5e8eebc
# ╠═5614c4b1-2c78-4fee-a480-3998a8e1ce6e
# ╠═efd7bcba-f393-4659-ac5b-006740085144
# ╠═0b9060d0-9c5a-4806-8d68-20e4ffd8f864
# ╠═7d71303b-defc-4e5f-beee-2bf684674ecf
# ╠═57052ba3-151c-4dda-904a-dac9c421241c
# ╠═a9a8f40f-da04-475b-aeaa-07add585ad9c
