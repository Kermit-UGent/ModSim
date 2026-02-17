### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 7e26d33b-2401-4e11-a457-bde8316df2f0
using Pkg; Pkg.activate("..")

# ╔═╡ 058d1250-30a2-11f0-05e1-4b134c2349ec
using Markdown, InteractiveUtils

# ╔═╡ c070007f-b4b8-4edc-8d04-700bf7e323de
using StatsPlots, StatsBase

# ╔═╡ 7e450e64-a8d6-47e9-9c2a-8bf8a4e2e71e
using ModelingToolkit, OrdinaryDiffEq, StochasticDiffEq

# ╔═╡ d1185b3f-5fa3-4d84-9a7f-c105b7c456ac
using ModelingToolkit: t_nounits as t, D_nounits as D

# ╔═╡ 6846bf57-0737-4ffd-87cf-955ecb4bd4de
using PlutoUI; TableOfContents()

# ╔═╡ fd3b7b7e-b573-4c14-afe1-a9a9083f39fe
md"""
# Introduction: Solving SDE problems with ModelingToolkit
"""

# ╔═╡ ba1450a1-5061-4352-994f-a308c053c3ce
# https://www.goldavenue.com/_next/image?url=https%3A%2F%2Fcdn.sanity.io%2Fimages%2Fyrlg7o79%2Fproduction%2F2974fbf9471a13b1a69294cfd81839266c779303-2560x1638.jpg&w=750&q=75
md"""
![Stock and volatility](https://www.goldavenue.com/_next/image?url=https%3A%2F%2Fcdn.sanity.io%2Fimages%2Fyrlg7o79%2Fproduction%2F2974fbf9471a13b1a69294cfd81839266c779303-2560x1638.jpg&w=750&q=75)
"""

# ╔═╡ e2ba0bc5-7147-462a-8eb9-76bf7145a51e
md"""
**Stochastic Differential Equations** (SDEs) are mathematical equations used to model
systems influenced by random noise. They extend **Ordinary Differential Equations** (ODEs) by incorporating terms that represent **stochastic processes**, typically in the form of a Wiener process or Brownian motion. SDEs are widely used in various fields, such as physics, biology, finance, and engineering, to describe the evolution of systems under uncertainty or with inherent randomness.
"""

# ╔═╡ 00b7ef23-1cfd-4923-918a-4d62de6ef17d
md"""
We will illustrate the concepts of SDE problems using a simplified **Heston model**.\
"""

# ╔═╡ 9352137e-5210-41fe-b6fe-f5e3c51d7eeb
md"""
In finance, the Heston model, named after Steven L. Heston, is a mathematical model that describes the evolution of the volatility of an underlying stock. It is a stochastic volatility model: such a model assumes that the volatility of the stock is not constant, nor even deterministic, but follows a random process.

- *Stock*: the goods or merchandise or asset (something having value)
- *Volatility*: a tendency to change quickly and unpredictably
"""

# ╔═╡ 843c9131-4eb3-4cdd-9136-a99171173b36
md"""
Our simplified Heston model assumes that $S$, the price of the stock, is determined by a stochastic process:

$$\cfrac{dS}{dt} = \mu\,S\,dt + \sqrt{V}\,S\,\cfrac{dW}{dt}$$

where the volatility $\sqrt{V}$ is given by a Cox-Ingersoll-Ross (CIR) process:

$$\cfrac{dV}{dt} = \kappa\,(\theta - V)\,dt - \sigma\,\sqrt{V}\,\cfrac{dW}{dt}$$

and $W$ is a Wiener process (i.e., continuous random walk). The value $V$, being the square of the volatility, is called the instantaneous *variance*.
"""

# ╔═╡ 93584181-33dc-4733-b548-945f5b80aaab
md"""
There are five parameters:
-  $\mu$: the expected return (drift of the stock)
-  $\kappa$: the rate at which the variance of the price reverts to long-term mean $\theta$
-  $\theta$: long-term mean of the variance of the price.
-  $\sigma$: the volatility of the volatility
"""

# ╔═╡ 60c709a9-4943-4eee-991b-95b4209500de
md"""
The values of the parameters are summerized in the following table:

| Parameter | Value |
|:---------- |:---------- |
| `μ`    | 0.03 |
| `κ`    | 0.90 |
| `θ`    | 0.04 |
| `σ`    | 0.08 |

We will use `100.0` and `0.04` as initial values for `S` and `V` respectively.
"""

# ╔═╡ 877ceee2-9e53-4ca9-9c75-8cb6252e6fef
md"""
## Implementation of the system
"""

# ╔═╡ 8d6c3bf9-35cf-48b9-8ce6-7ef711b6fe03
md"""
We will create a MTK model for the aforementioned problem in order to simulate the evolution of the price of the stock $S$ with time as a Stochastic Differential Equation (SDE) problem with the noisy factors in the equations above.
"""

# ╔═╡ 3ff63798-3337-4764-87a9-fe2dd037232c
md"""
We define the variables `S` and `V` with `@variables` and we can immediately add default initial conditions.
"""

# ╔═╡ cdaf2c50-1a58-4f04-ad7d-cf7b0bf39f96
@variables S(t)=100.0 V(t)=0.04

# ╔═╡ 762847df-5bba-4f95-96aa-d30fc8852daa
md"""
We define the parameters `μ`, `κ`, `θ` and `σ` with `@parameters` and we can immediately add their default values. In addition we add the scaling parameter `n` and we set it to `0.5`. We this parameter you can literally scale the noise.
"""

# ╔═╡ 944f859c-24fb-4f1e-a5a1-2ee248e99965
@parameters μ=0.03 κ=0.90 θ=0.04 σ=0.08 n=0.5

# ╔═╡ ff673d4f-11b1-4e52-a63e-19e007e03b5f
md"""
Next, we instantiate a source `B` of Brownian noise (a Wiener process by default) with `@brownian`.
"""

# ╔═╡ 48b0b8d1-a569-48c6-b2a8-10de2273752f
@brownian B

# ╔═╡ f890b80d-61eb-4eb7-89f2-2f733fd00c33
md"""
Now we set up the model equations including the noise terms `sqrt(V)*S*n*B` and `-σ*sqrt(V)*n*B`, and put them in a vector/array. We have named the array `eqs_spsv`.
"""

# ╔═╡ 1bc1cdba-43c6-4455-8bc7-2751f5538021
eqs_spsv = [
	D(S) ~ μ*S + sqrt(V)*S*n*B,
	D(V) ~ κ*(θ - V) - σ*sqrt(V)*n*B
]

# ╔═╡ f5b1f97d-9409-4d92-8a6c-2ecfc3400340
md"""
Now we are ready to build a system of equations with `@mtkbuild` and the function `System`. We provide the following arguments to `System`: the model equations and `t`. We have named our system `sys_spsv`.
"""

# ╔═╡ 68ca9da8-a309-4fe0-a754-61305f8e70a9
@mtkbuild sys_spsv = System(eqs_spsv, t)

# ╔═╡ 55bbe791-089f-4b01-acfa-c2fd6a18edf0
md"""
## Simulating the system as an SDE problem
"""

# ╔═╡ 4baab57b-ca86-46e4-b70a-ff0c4ff91500
md"""
If you want to modify the default initial conditions or parameter values, you can define new vectors for them. For the sake of completeness, we will define new vectors `u0` and `parms` but using the same default values as were set before.
"""

# ╔═╡ 915df4bd-6ab3-4267-80e2-d4613d7d46f7
md"""
### Setting initial conditions
"""

# ╔═╡ 4dd3594d-033a-479f-b3b6-48458d7f7a21
u0 = [S=>100.0, V=>0.04]

# ╔═╡ 76080f37-ddb2-404d-ac40-246e4eff94f1
md"""
### Setting the timespan
"""

# ╔═╡ 1f02b734-2519-44bc-b304-924fce19c735
md"""
We will simulate the price of the stock in the time interval `(0.0, 5.0)`.
"""

# ╔═╡ fa33322d-e9f3-4065-99c0-c356ac871d19
tspan = (0.0, 5.0)

# ╔═╡ 4315f33a-27f3-41fc-9ef1-8645d69c0452
md"""
### Setting parameter values
"""

# ╔═╡ c97ab6e7-bc7b-41d4-b8ad-94c1fb9ec64c
parms = [μ=>0.03, κ=>0.90, θ=>0.04, σ=>0.08, n=>0.5]

# ╔═╡ 44aa6e33-9517-442b-97a9-db2b15f8f32a
md"""
### Creating an SDEProblem
"""

# ╔═╡ 301599b6-bb60-49a2-90ec-1bd33adc3748
md"""
We will create the SDE problem using the function `SDEProblem`. As with the function `ODEProblem` you need to provide the same kind of arguments: the name of the MTK system, the vector of initial values, the time span and the vector of parameters. Since the initial conditions and the parameter values were given as default (cf. when we defined them with `@variables` and `@parameters`), we could just put empty vectors (cf. `[]`) if we like. We have named our SDE problem `sprob_spsv`.
"""

# ╔═╡ fd3faa3e-f530-46a8-901d-3d32e2aad5e4
sprob_spsv = SDEProblem(sys_spsv, u0, tspan, parms) # using values in u0 and parms
# sprob_spsv = SDEProblem(sys_spsv, [], tspan, [])  # using the default values

# ╔═╡ 128438ca-57e4-4424-9cb7-44b282f67ac7
md"""
### Solving the SDEProblem

There are many solving methods available for solving SDE problems. You can find a [list of methods here](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/#Full-List-of-Methods). We will simply use the first one in this list, namely `EM()`, with the time step option `dt=0.01` that will introduce some randomness at every time step. We have stored the solution in `ssol_spsv`.
"""

# ╔═╡ db240bb7-5011-4112-b55e-e05b7ccb769d
ssol_spsv = solve(sprob_spsv, EM(), dt=0.01)

# ╔═╡ ec980078-64d8-490b-9606-64f38fd3e4ef
md"""
We will only plot the variable `S` using the function `plot` with the option `idxs=[S]`. In our simplified model, the graph of `V` will be exactly the negative of `S` but on a much smaller scale.
"""

# ╔═╡ c85137a4-edc1-40b0-ac90-02125dc5106a
plot(ssol_spsv; idxs=[S])

# ╔═╡ 9b32d5cb-f05b-4775-b3f7-6d286e13fcb5
md"""
## Simulating the system as an EnsembleProblem.

In order to see to have an idea of the extend of the stochastic effect on the solutions, we can create a so-called *EnsembleProblem*. This allows us to plot many possible solutions in one plot.

In order to create an *EnsembleProblem*, you need to create an *SDEProblem* first. Since we already have our *SDEProblem* called `sprob_spsv`, we can readily create an *EnsembleProblem* from this. All you need to do is call the function `EnsembleProblem` with the name of the SDE problem object (cf. `sprob_spsv`) as argument.
"""

# ╔═╡ 53629d40-01df-4bfd-bb29-05506cc23a39
esprob_spsv = EnsembleProblem(sprob_spsv)

# ╔═╡ decfa049-d033-4049-868e-fa61669e2024
md"""
### Solving the EnsembleProblem

Solving the ensemble problem can be done with our, yet familiar, function `solve` as we did when solving the SDE problem, but now we need to provide how many trajectories (simulations) you want to make. If you want an ensemble of 100 simulations, you can put `trajectories=100`.
"""

# ╔═╡ e4fdecc6-39a8-41c2-ac5c-2c9401393a4a
essol_spsv = solve(esprob_spsv, EM(), dt=0.01, trajectories=100)

# ╔═╡ da73ed53-cd93-4b03-a14d-ba06a52275e5
md"""
We can now plot all simulations of `S` simply as before with the function `plot`.
"""

# ╔═╡ d419be32-54e4-4262-9656-c6690ed1cef1
plot(essol_spsv; idxs=[S])

# ╔═╡ d85efd34-e78e-4e0e-9989-dd0b849bbe5e
md"""
If you want to inspect the end value of $S$ of the 32th simulation, you can do as follow:
"""

# ╔═╡ 2d3d3fcd-fa00-4cf2-a0ba-14a9b412195a
essol_spsv[32][S][end]

# ╔═╡ 528e0012-c275-4d88-9fce-8624af804750
md"""
In order to have an idea of the spread of the end value of S after 5 units of time for all simulations, we will store all the end values in a vector `S_end_vals` using a `for`-loop as illustrated in the following code-cell:
"""

# ╔═╡ a487e6a0-bb20-4ff6-9e38-242c49078f13
begin
	S_end_vals = []                # define an empty vector 
	for i = 1:length(essol_spsv)   # loop index-wise over the end values 
		# append the end value of S of the i-th simulation to S_end_vals
		append!(S_end_vals, essol_spsv[i][S][end])
	end
end

# ╔═╡ 4e6dab85-b377-4989-8ad7-80fc45ee78aa
md"""
Next, we will plot a histogram using the end values. Hereby, we will group our data in `bins` of size 5 in the range 0 to 200 (cf. the option `bins=0:5:200`).
"""

# ╔═╡ b1eaa791-0213-44ab-9df9-e6cee2f9ef2d
histogram(S_end_vals; bins=0:5:200)

# ╔═╡ b54452ac-3792-425b-886a-0ce218498d9b
md"""
You can determine the mean value of all end values of S in the following way:
"""

# ╔═╡ 95a28316-028e-431a-920c-e5ffe28bc8a9
mean(S_end_vals)

# ╔═╡ 078cd51d-bea4-4408-a35b-255928d39e32
md"""
You can determine the standard deviation of all end values of S in the following way:
"""

# ╔═╡ c5bf45b5-9c5f-4f19-96b7-afde835dde19
std(S_end_vals)

# ╔═╡ Cell order:
# ╠═058d1250-30a2-11f0-05e1-4b134c2349ec
# ╠═7e26d33b-2401-4e11-a457-bde8316df2f0
# ╠═c070007f-b4b8-4edc-8d04-700bf7e323de
# ╠═7e450e64-a8d6-47e9-9c2a-8bf8a4e2e71e
# ╠═d1185b3f-5fa3-4d84-9a7f-c105b7c456ac
# ╠═6846bf57-0737-4ffd-87cf-955ecb4bd4de
# ╟─fd3b7b7e-b573-4c14-afe1-a9a9083f39fe
# ╟─ba1450a1-5061-4352-994f-a308c053c3ce
# ╟─e2ba0bc5-7147-462a-8eb9-76bf7145a51e
# ╟─00b7ef23-1cfd-4923-918a-4d62de6ef17d
# ╟─9352137e-5210-41fe-b6fe-f5e3c51d7eeb
# ╟─843c9131-4eb3-4cdd-9136-a99171173b36
# ╟─93584181-33dc-4733-b548-945f5b80aaab
# ╟─60c709a9-4943-4eee-991b-95b4209500de
# ╟─877ceee2-9e53-4ca9-9c75-8cb6252e6fef
# ╟─8d6c3bf9-35cf-48b9-8ce6-7ef711b6fe03
# ╟─3ff63798-3337-4764-87a9-fe2dd037232c
# ╠═cdaf2c50-1a58-4f04-ad7d-cf7b0bf39f96
# ╟─762847df-5bba-4f95-96aa-d30fc8852daa
# ╠═944f859c-24fb-4f1e-a5a1-2ee248e99965
# ╟─ff673d4f-11b1-4e52-a63e-19e007e03b5f
# ╠═48b0b8d1-a569-48c6-b2a8-10de2273752f
# ╟─f890b80d-61eb-4eb7-89f2-2f733fd00c33
# ╠═1bc1cdba-43c6-4455-8bc7-2751f5538021
# ╟─f5b1f97d-9409-4d92-8a6c-2ecfc3400340
# ╠═68ca9da8-a309-4fe0-a754-61305f8e70a9
# ╟─55bbe791-089f-4b01-acfa-c2fd6a18edf0
# ╟─4baab57b-ca86-46e4-b70a-ff0c4ff91500
# ╟─915df4bd-6ab3-4267-80e2-d4613d7d46f7
# ╠═4dd3594d-033a-479f-b3b6-48458d7f7a21
# ╟─76080f37-ddb2-404d-ac40-246e4eff94f1
# ╟─1f02b734-2519-44bc-b304-924fce19c735
# ╠═fa33322d-e9f3-4065-99c0-c356ac871d19
# ╟─4315f33a-27f3-41fc-9ef1-8645d69c0452
# ╠═c97ab6e7-bc7b-41d4-b8ad-94c1fb9ec64c
# ╟─44aa6e33-9517-442b-97a9-db2b15f8f32a
# ╟─301599b6-bb60-49a2-90ec-1bd33adc3748
# ╠═fd3faa3e-f530-46a8-901d-3d32e2aad5e4
# ╟─128438ca-57e4-4424-9cb7-44b282f67ac7
# ╠═db240bb7-5011-4112-b55e-e05b7ccb769d
# ╟─ec980078-64d8-490b-9606-64f38fd3e4ef
# ╠═c85137a4-edc1-40b0-ac90-02125dc5106a
# ╟─9b32d5cb-f05b-4775-b3f7-6d286e13fcb5
# ╠═53629d40-01df-4bfd-bb29-05506cc23a39
# ╟─decfa049-d033-4049-868e-fa61669e2024
# ╠═e4fdecc6-39a8-41c2-ac5c-2c9401393a4a
# ╟─da73ed53-cd93-4b03-a14d-ba06a52275e5
# ╠═d419be32-54e4-4262-9656-c6690ed1cef1
# ╟─d85efd34-e78e-4e0e-9989-dd0b849bbe5e
# ╠═2d3d3fcd-fa00-4cf2-a0ba-14a9b412195a
# ╟─528e0012-c275-4d88-9fce-8624af804750
# ╠═a487e6a0-bb20-4ff6-9e38-242c49078f13
# ╟─4e6dab85-b377-4989-8ad7-80fc45ee78aa
# ╠═b1eaa791-0213-44ab-9df9-e6cee2f9ef2d
# ╟─b54452ac-3792-425b-886a-0ce218498d9b
# ╠═95a28316-028e-431a-920c-e5ffe28bc8a9
# ╟─078cd51d-bea4-4408-a35b-255928d39e32
# ╠═c5bf45b5-9c5f-4f19-96b7-afde835dde19
