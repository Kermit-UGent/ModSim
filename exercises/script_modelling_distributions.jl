using Distributions
using Plots, LinearAlgebra, Random, LaTeXStrings


# Random variables and probability distributions

dist = Weibull(2, 7)

# probability density function
pdf(dist, 0.5)
pdf(dist, 5.0)
pdf(dist, 17.0)

f_X(x) = pdf(dist, x)

f_X(5.0)

plot(f_X, xlim=(0, 20), ylabel=L"f_X(x)", xlabel=L"x", label="pdf")

# cumulative density function
cdf(dist, 5.0)          # P(X <= 5)

F_X(x) = cdf(dist, x)

F_X(5) - F_X(1)         # P(1 < X <= 5)

plot(F_X, xlim=(0, 20), ylabel=L"F_X(x)", xlabel=L"x", label="cdf")

# Expected value and other moments

mean(dist)

var(dist)

std(dist)

sqrt(var(dist))

median(dist)

mode(dist)

# Conditional probabilities

# ...


# Sampling from a distribution

# The Law of Large Numbers

rand(dist)

rand(dist, 10)

g(r) = 4π * r^2

dr = 1e-4

# Here we 'integrate' g(r)*f_X(r) from 0 to 50 in steps of dr
# Hier heeft r de distributie dist
A_integrate = sum(r->g(r) * pdf(dist, r) * dr, 0:dr:50)
# r-> betekent: neem de r uit 0:dr:50
# dus g(r)*f_X(r)*dr wordt gesommeerd voor alle r's uit 0:dr:50

g(mean(dist))

# Hier nemen we 100000 willekeurige samples uit dist
# en steken die in de functie g om dan van al die waarden
# de gemiddelde van de bepalen.
A_sampling = mean(g, rand(dist, 100_000))   # E[g(x)]

# Generating random numbers

randn()

rand(10)    # 10 uniformly distributed numbers

# Make a normal distribution with μ = 0.0 and σ = 1.0
normal = Normal(0, 1)

# Take a random sample from this distribution
rand(normal)

# Take 10 random sample from this distribution
rand(normal, 10)

# Combining simple distributions into complex ones

# Joint distributions of independent variables

distX = Laplace(8, 6)        # midden op 8, maat voor breedte is 6
g_X(x) = pdf(distX, x)
plot(g_X, xlim=(-4, 20))

distY = TriangularDist(-2, 2, 1)    # driehoek met punten (-2, 0), (2, 0) en (1, 0.5)
g_Y(y) = pdf(distY, y)
plot(g_Y, xlim=(-4, 20))
