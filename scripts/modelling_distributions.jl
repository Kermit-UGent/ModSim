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
using Distributions

# ╔═╡ c175d1ab-28be-4766-bc0c-20cbc72f191d
using PlutoUI, Plots, LinearAlgebra, Markdown, Random, LaTeXStrings

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
md"n throws : $(@bind n_pi Slider(1000:1000:100_000, default=1000, show_value=true))"

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

# ╔═╡ 41967e5b-9068-4318-acb6-d958b608c8f7
begin
	plot(x->sqrt(1-x^2), 0, 1, aspect_ratio=:equal, xlim=[0,1], ylim=[0, 1], lw=2, fillrange=zero, fillalpha=0.3, xlab=L"x", ylab=L"y", label="circle", title="Monte Carlo for estimating π/4")
	scatter!(x_unif[in_circle], y_unif[in_circle], ms=0.8, markercolor="orange", markerstrokewidth=0, label="in circle")
	scatter!(x_unif[.!in_circle], y_unif[.!in_circle], ms=0.8, markercolor="green", markerstrokewidth=0, label="out of circle")
end

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

# ╔═╡ b9be819e-2db9-4952-9453-98bf815b00e8
let
	pfx = plot(f_X, 0, 20, label="pdf", xlab=L"x", lw=2, ylab=L"f_X(x)")
	plot!(pfx, zero, 0, xpdf, fillrange=f_X, fillalpha=0.5, label="P(X≤$xpdf)")
	pFx = plot(F_X, 0, 20, label="cdf", lw=2)
	scatter!(pFx, [xpdf], [F_X(xpdf)], label="P(X≤$xpdf)", xlab=L"x", ylab=L"F_X(x)")
	plot(pfx, pFx, layout=(2,1))
end

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

# ╔═╡ 7779254a-83b9-4849-957c-af4b62c3916a
let
	ps = []
	for n in [10, 100, 1000, 10_000]
		p = histogram(rand(dist, n), label="", title="n=$n", normalize=true)
		plot!(x -> pdf(dist, x), 0, 20, label="", lw=2)
		push!(ps, p)
	end
	plot(ps...)
end

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

# ╔═╡ d7e88a1a-73e9-48ff-bb51-2292346f88e0
begin
	scatter(rand(200), rand(200), label="", title="200 pseudo-random numbers", aspect_ratio=:equal, xlims=[0,1])
end

# ╔═╡ 74e50a3d-7623-45c2-8d3a-2c34f363177e
begin
	s = SobolSeq(2)
	p = reduce(hcat, next!(s) for i = 1:200)'
	scatter(p[:,1], p[:,2], label="", title="200 Sobol quasi-random numbers", aspect_ratio=:equal, xlims=[0,1])
end

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

# ╔═╡ e3d9a5b7-73ec-4c16-8799-2e2a4f99c018
plot(n->1/√(n), 1, 100, lw=2, label=L"1/\sqrt{n}", ylab=L"\sigma_{\bar{X}_n}", xlabel=L"n", title="Standard error of the mean")

# ╔═╡ e94fa51a-e1a7-4a1b-88ef-76b38d3a6d74
plot(n->1/√(n), 1, 1000_000, lw=2, yscale=:log10, xscale=:log10, label=L"1/\sqrt{n}", ylab=L"\sigma_{\bar{X}_n}", xlabel=L"n", title="Standard error of the mean (log-scale)")

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

# ╔═╡ 222cdcae-5b51-4fdc-a794-11471449ebcd
let
# Binomial
	pbin = scatter(k->pdf(Binomial(20, .5), k), 0:25, label="n=25, p=0.5", xlab=L"k", ylab=L"P(X=k)")
	scatter!(pbin, k->pdf(Binomial(20, .25), k), 0:25, label="n=25, p=0.25", marker=:^)
	title!("Binomial distribution")
	
	# Poisson
	ppois = scatter(k->pdf(Poisson(1), k), 0:20, label="λ=1", xlab=L"k", ylab=L"P(X=k)")
	scatter!(ppois, k->pdf(Poisson(5), k), 0:20, label="λ=5", marker=:^)
	scatter!(ppois, k->pdf(Poisson(10), k), 0:20, label="λ=10", marker=:v)
	title!("Poisson distribution")
	
	# Geometric
	pgeo = scatter(k->pdf(Geometric(0.1), k), 1:25, label="p=0.1", xlab=L"k", ylab=L"P(X=k)")
	scatter!(pgeo, k->pdf(Geometric(0.2), k), 0:25, label="p=0.2", marker=:^)
	scatter!(pgeo, k->pdf(Geometric(0.05), k), 0:25, label="p=0.05", marker=:v)
	title!("Geometric distribution")
	
	# Geometric
	punif = scatter(k->pdf(DiscreteUniform(5, 10), k), 0:30, label="a=5, b=10", xlab=L"k", ylab=L"P(X=k)")
	scatter!(punif, k->pdf(DiscreteUniform(10, 25), k), 0:30, label="a=10, b=25", marker=:^)
	title!("Uniform distribution")
	
	plot(pbin, ppois, pgeo, punif)
end

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

# ╔═╡ ebe190e6-da7a-4c88-bd31-b5dc53f6a3e1
let
	# normal
	pnorm = plot(xlabel=L"x", ylabel=L"f_X(x)", title="Normal")
	for (μ, σ) in [(0., 1.), (2., 1.), (0.0, 2.0), (-1, 2)]
		plot!(x->pdf(Normal(μ, σ),x),-5, 5, lw=2, label="μ=$μ, σ=$σ")
	end
	pnorm

	# normal
	plaplace = plot(xlabel=L"x", ylabel=L"f_X(x)", title="Laplace")
	for (μ, θ) in [(0., 1.), (2., 1.), (0.0, 2.0), (-1, 2)]
		plot!(x->pdf(Laplace(μ, θ),x),-8, 8, lw=2, label="μ=$μ, θ=$θ")
	end
	plaplace

	plot(pnorm, plaplace)

end

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
### Jont distributions of independent variables

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

# ╔═╡ 4bd5e387-ac06-4f4b-822f-cbc89e7509f9
plot(y->pdf(distY, y), -2:0.01:2, label="Y | X=$xslice", color="orange", lw=2, xlabel="y", )

# ╔═╡ 46b4ef15-bd34-4b59-aef5-5dc4e1193ec3
md"y conditioning : $(@bind yslice Slider(-2:0.01:2, default=0, show_value=true))"

# ╔═╡ f2faec3d-3d44-43e9-bf03-37ba68e37300
let
	contourf(-10:0.05:25, -2:0.01:2, (x,y)->pdf(dist_prod, [x,y]), color=:speed, xlab=L"x", ylab=L"y")
	vline!([xslice], label="Y | X=$xslice", color="orange", lw=2)
	hline!([yslice], label="X | Y=$yslice", color="blue", ls=:dash, lw=2)
	title!("Joint PDF of a Laplace and Triangular distribution")
end

# ╔═╡ c7b2fc1e-9dd3-42e9-aa55-f9aea4781784
plot(x->pdf(distX, x), -10:0.05:25, label="f_{X | Y=$yslice}", color="blue", ls=:dash, lw=2, xlabel="x", )

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

# ╔═╡ 94dd1214-e63c-44d5-b5c3-9ca834a89ebf
let
	plot(x->pdf(dist_mixture, x), -5, 10, label="mixture", lw=2, xlab=L"x")
	plot!(x->pdf(d1, x), -5, 10, label="component 1", ls=:dash)
	plot!(x->pdf(d2, x), -5, 10, label="component 2", ls=:dash)
end

# ╔═╡ f1906adc-79f7-44da-a085-486e8f46e027
md"""
Using Bayes' rule, one can estimate the probability that an observation originates from one of the $k$ distributions:

$$P(\text{sample $x$ originates from component $i$}) = \frac{w_i \times f_{X_i}(x)}{\sum_{i=1}^kw_if_{X_i}(x)}\propto w_i \times f_{X_i}(x)\,.$$
"""

# ╔═╡ d6cd5202-0171-40f9-bc17-d001dadee567
md" x : $(@bind xm Slider(-5:0.1:10, default=0, show_value=true))"

# ╔═╡ fc3a7a70-89e5-43e4-8dae-b7d2d8b54355
let
	plot(x->pdf(dist_mixture, x), -5, 10, label="mixture", lw=2, xlab=L"x", ylab=L"f_X(x)")
	plot!(x->pdf(d1, x), -5, 10, label="component 1", ls=:dash)
	plot!(x->pdf(d2, x), -5, 10, label="component 2", ls=:dash)
	vline!([xm], color=:gold, lw=2, label="x=$xm")
end

# ╔═╡ e61e45a0-d9ec-4ab6-8c1b-fc8579a9a9b3
let
	likelihood = [pdf(d1, xm), pdf(d2, xm)]
	posterior = likelihood .* w ./ pdf(dist_mixture, xm)
	bar(posterior, xticks=1:length(w), xlabel="component", ylabel="posterior probability", label="X=$xm", color=:gold)
end

# ╔═╡ 89178521-511e-4c53-a2e1-c1d4dad03a23
md"Mixtures allow for creating complex distributions of more simple variables. We have seen that product distributions are an easy way of combining different random variables $X$, $Y$ into a joint distribution $X\times Y$. They are not very exciting in themselves as there is no dependence between the two variables. However, combining multiple product distributions $X_1\times Y_1$, $X_2\times Y_2$, ... into a mixture, we can introduce some sophisticated dependencies between the variables! This way of combining distributions are sometimes called sum-product networks, for the operations that are used to combine them."

# ╔═╡ 67b154c9-3c14-4eb2-98f2-19f9e3f74d16
begin
	dmv1 = product_distribution([Normal(-2), Normal(-3)])
	dmv2 = product_distribution([Normal(2, 0.7), Normal(3, 2)])
	mvmixture = MixtureModel([dmv1, dmv2], [0.6, 0.4])
end

# ╔═╡ ebf13b85-7fa0-4859-88d1-a5189c5bbd40
let

	contourf(-5:0.1:5, -5:0.1:8, (x,y)->pdf(mvmixture, [x,y]), color=:speed, xlab=L"x", ylab=L"y", title="mixture of two MVN")
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

# ╔═╡ 40fbbc43-4a73-4d25-a37c-b39798d60761
contourf(0:0.1:30, 0:0.1:30, (x,y)->logpdf(spn, [x,y]), color=:speed, xlab=L"x", ylab=L"y")

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

# ╔═╡ ccbcce90-c46f-4f61-8907-2ccc3e692184
let
	mvn = MultivariateNormal(μ, Σ)
	contourf(-5:0.05:5, -5:0.05:5, (x,y)->pdf(mvn, [x,y]), color=:speed, xlab=L"x", ylab=L"y")
	title!("PDF of a MVN")
end

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

# ╔═╡ 6602c372-158d-4b40-b7fe-aea7e534bd39
plot(p->pdf(dist_p, p), 0, 1, lw=2, label="Beta(8,3)", xlab=L"x", ylab=L"f_p(p)")

# ╔═╡ cc058abb-e7b6-433e-aacf-13295a401c5f
md"We use a `Turing` model to encode our distribution."

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

# ╔═╡ 77664b93-4b90-43b3-a0cc-0b09c7bd9854
histogram(germination_counts, xticks=0:10,
		xlabel=L"k", ylabel="frequency in sample",
		label="germination succes X", title="Numer of pepper seeds that germinated")

# ╔═╡ f24fffb4-deb7-428a-ba21-0e4999962624
md"How many packages in this large sample have more than 8 seeds germinating?"

# ╔═╡ 6d565e90-7aa6-480f-96fb-25c745ec6974
count(((p, X),)->X>8, seed_sample)

# ╔═╡ 42867e88-cd02-4634-a225-a9bf0ae46b70
count(>(8), germination_counts)

# ╔═╡ 13cf931e-5333-46ff-a00a-ea94c0b2e509
mean(>(8), germination_counts)  # P(X > 8)

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

# ╔═╡ cda5faf1-8d74-4679-a335-7a5974010f5a
scatter(first.(rainfall_samples), last.(rainfall_samples), alpha=0.2, xlab="intensity (mm/h)", ylab="duration (h)", label="")

# ╔═╡ 4b96292b-f0ad-48ac-99cd-360c0708ef72
rainfall_amount = prod.(rainfall_samples)

# ╔═╡ f7ca0af5-dcdf-4e75-9063-2d194738be64
histogram(rainfall_amount, xlab="rainfall amount (mm)", ylab="frequency", label="intensity * duration")

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

# ╔═╡ 79be0044-9856-4b3b-80e3-6c362e23b9f4
p_cond = [p for (p, X) in seed_sample if X==ngerm]

# ╔═╡ 18e788ca-5b74-4a0b-81f0-5f307857e7e2
mean(>(0.9), p_cond)  # probability of p > 0.9 given x

# ╔═╡ 037a2175-22b2-48c9-9664-3743b11e8f7c
mean(p_cond), var(p_cond)  # mean and variance of the conditional distr. of p

# ╔═╡ 6fe05275-f1f7-4b3f-bbc2-8077bf981391
length(p_cond)

# ╔═╡ cbb71827-67cd-43aa-8a36-757bdf98e16a
histogram(p_cond, ylabel="frequency", xlab=L"p", label="p given n=$ngerm", xlims=[0,1])

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

# ╔═╡ 73898b6c-d8bd-49a9-ba36-a068e62a7ec3
let
	d = Binomial(30, 0.3)
	#d = dist
	p = hline([mean(d)], linestyle=:dash, lw=2, color=:red, label=L"E[X]", xlab=L"n", ylabel=L"\bar{X}_n", title="The Law of Large Numbers")
	for i in 1:5
		plot!(cummean(rand(d, 250)), lw=2, label="", alpha=0.6)
	end
	#ylims!(8, 10)
	p
end

# ╔═╡ 67594833-14c8-423c-94ce-91ce2a44b6c9
let
	Random.seed!(12)
	psob = reduce(hcat, next!(s) for i = 1:1024)'
	ppseu = rand(1024, 2)
	
	pi_sob = pi .- 4cummean(norm.(eachrow(psob)) .≤ 1) .|> abs
	pi_pseu = pi .- 4cummean(norm.(eachrow(ppseu)) .≤ 1) .|> abs
	
	plot(pi_pseu, yscale=:log10, label="pseudo-random", xlab=L"n", ylab=L"|E[X]-\bar{X}_n|", lw=2)
	plot!(pi_sob, label="quasi-random", lw=2)
	title!("Error estimating π using sampling")
end

# ╔═╡ 57052ba3-151c-4dda-904a-dac9c421241c


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Markdown = "d6f4376e-aef5-505a-96c1-9c027394607a"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Sobol = "ed01d8cd-4d21-5b2a-85b4-cc3bdc58bad4"
Turing = "fce5fe82-541a-59a6-adf8-730c64b5f9a0"

[compat]
Distributions = "~0.25.107"
LaTeXStrings = "~1.3.1"
Plots = "~1.39.0"
PlutoUI = "~0.7.55"
Sobol = "~1.5.0"
Turing = "~0.30.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.0"
manifest_format = "2.0"
project_hash = "546a4285cc5ca9ea0d4c38b968585b216c5aeec6"

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
git-tree-sha1 = "0f748c81756f2e5e6854298f11ad8b2dfae6911a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.0"

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
deps = ["AbstractMCMC", "AxisArrays", "Dates", "Distributions", "IteratorInterfaceExtensions", "KernelDensity", "LinearAlgebra", "MCMCDiagnosticTools", "MLJModelInterface", "NaturalSort", "OrderedCollections", "PrettyTables", "Random", "RecipesBase", "Statistics", "StatsBase", "StatsFuns", "TableTraits", "Tables"]
git-tree-sha1 = "d28056379864318172ff4b7958710cfddd709339"
uuid = "c7f686f2-ff18-58e9-bc7b-31028e88f75d"
version = "6.0.6"

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

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

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
git-tree-sha1 = "9e8fed0505b0c15b4c1295fd59ea47b411c019cf"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.2"

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

[[deps.Sobol]]
deps = ["DelimitedFiles", "Random"]
git-tree-sha1 = "5a74ac22a9daef23705f010f72c81d6925b19df8"
uuid = "ed01d8cd-4d21-5b2a-85b4-cc3bdc58bad4"
version = "1.5.0"

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
# ╠═2d346892-cb15-11ee-2d81-73e08fcc3288
# ╠═c175d1ab-28be-4766-bc0c-20cbc72f191d
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
# ╟─2b79770c-00a8-4c3b-b39e-b31400593ea4
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
# ╟─af7119a0-2470-4bf0-8b30-59b38ddd8d5d
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
# ╠═57052ba3-151c-4dda-904a-dac9c421241c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
