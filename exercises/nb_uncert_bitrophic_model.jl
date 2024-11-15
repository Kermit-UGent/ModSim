### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 7d2ce89c-de3e-4fff-8f23-8756351b377e
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 079bc7cc-0c53-11ef-2323-8b475737e481
using Markdown

# ╔═╡ 3d189b53-8aff-4ca0-a784-11b0bca02374
using InteractiveUtils

# ╔═╡ d89c7b51-a94c-42e5-94f1-2dc1e641dfb3
using Catalyst

# ╔═╡ 49c93319-f76d-485e-a2e2-e06210e25df3
using DifferentialEquations, Plots

# ╔═╡ b8c575fd-a94b-4221-8751-74c51f10ec7d
using Measurements

# ╔═╡ bedf7419-590c-47fd-8089-b890dded0468
md"
### Exercise: Bitrophic model - Uncertainty analysis
"

# ╔═╡ 8f795d92-b842-4392-8bc3-956080fcc6a8
md"
In one of the previous practica we were introduced to a bitrophic model in which the dynamic relationship between a field crop $C$ and a voracious insect population $A$ within an ecosystem was modelled.

$$\begin{eqnarray*}
\frac{dC}{dt} &= \theta C \left(1-\frac{C}{k}\right)-fCA \\
\frac{dA}{dt} &= \phi f CA -(1 + p)\, \mu A
\end{eqnarray*}$$
"

# ╔═╡ bc2e3724-de5f-4455-a4be-019c6a8accae
md"
The *reaction network object* for this model could be set-up as:
"

# ╔═╡ 90e9da81-7576-43bb-b4ce-1293a48bf3ec
bitrophic_model = @reaction_network begin
	θ*(1-C/k), C --> 2C
	f, C+A --> (1+ϕ)*A
    (1+p)*μ, A --> 0
end

# ╔═╡ 465e4ea9-7066-4929-98c4-e5c378ef34d8
md"

Assume the uncertainties in the following parameter values:

-  $\theta=0.20 \pm 0.02\;d^{-1}$
-  $\phi=0.20 \pm 0.02$
-  $p=3.0 \pm 0.2$ 

and that the uncertainty in the other parameters values: $k=4000\;kg/ha$, $f=0.001\;ha/(kg\,d)$ and $\mu=0.1\;d^{-1}$ are negligible. Suppose that at the beginning of a season, $100\;kg$ of the crop and $0.5\;kg$ of insects per $ha$ are present. Perform an uncertainty analysis by plotting the uncertainty bands on the simulation results of $C$ and $A$ in a timespan of $[0, 200]\,days$.

Interpret your results.
"

# ╔═╡ 1129972e-3e3d-4b3b-8c2c-3c57eaa3728a
md"
Initialize a vector `u₀` with the initial conditions, and set the timespan:
"

# ╔═╡ 79802166-efbd-4944-ae58-fb26a5d9ad4b
# u₀ = missing               # Uncomment and complete the instruction
u0 = [:C => 100, :A => 0.5]

# ╔═╡ 146855a0-5cda-4783-89af-7918c067e952
# tspan = missing
tspan = (0.0, 200.0)         # Uncomment and complete the instruction

# ╔═╡ 4e8fea9e-e16c-4b17-8a6c-52bd8721cfd8
md"
We initialize a vector `params_uncert` with the parameter values and their corresponding uncertainty:
"

# ╔═╡ 0588cc16-a4b1-4b3f-8ba9-742680996e11
# params_uncert = missing     # Uncomment and complete the instruction
params_uncert = [:θ => 0.20±0.02, :k => 4000.0, :f => 0.001, :ϕ => 0.2±0.02, :p => 3.0±0.2, :μ => 0.1]

# ╔═╡ 6893ccdb-7f6c-4ec9-8d6b-520e784acd03
md"
We create the corresponding ODE problem and store it in `oprob_uncert`:
"

# ╔═╡ 60a6fb20-93d2-47d3-8ee1-3087d48c569d
# oprob_uncert = missing       # Uncomment and complete the instruction
oprob_uncert = ODEProblem(bitrophic_model, u0, tspan, params_uncert, combinatoric_ratelaws=false)

# ╔═╡ 0739d816-0456-4cc2-beb7-2e926cf2e229
md"
We solve the ODE problem. Use `Tsit5()` and `saveat=2.0`. Store the solution in `osol_uncert`:
"

# ╔═╡ f829fc66-3df8-4ac1-a1bb-8eb2ed722423
# osol_uncert = missing          # Uncomment and complete the instruction
osol_uncert = solve(oprob_uncert, Tsit5(), saveat=2.0)

# ╔═╡ d4d6f69c-844b-488a-9170-b09e7c22dd25
md"
Plot the results (simulation of the output variables $C$ and $A$ together with their uncertainty band):
"

# ╔═╡ 2152a72b-31cd-44b9-90d8-615610e8d261
# missing                        # Uncomment and complete the instruction
plot(osol_uncert)

# ╔═╡ 600dd2e9-aed9-41da-a510-9ba78788f207
md"
Try to relate the local sensitivity analysis to the uncertainty analysis. Hence, study the effect of the individual parameter uncertainties on the output variables $C$ and $A$ and compare with your local sensitivity results of the corresponding parameter.

In order to do that, analyse the effect on the uncertainty bands for $C$ and $A$ by taking one uncertainty on a parameter at a time. In other words, analyse the subsequent cases seperately:
- Assume uncertainty only in $\theta$
- Assume uncertainty only in $\phi$
- Assume uncertainty only in $p$

Draw your conclusions:
- Answer: missing
"

# ╔═╡ Cell order:
# ╠═079bc7cc-0c53-11ef-2323-8b475737e481
# ╠═3d189b53-8aff-4ca0-a784-11b0bca02374
# ╠═7d2ce89c-de3e-4fff-8f23-8756351b377e
# ╠═d89c7b51-a94c-42e5-94f1-2dc1e641dfb3
# ╠═49c93319-f76d-485e-a2e2-e06210e25df3
# ╠═b8c575fd-a94b-4221-8751-74c51f10ec7d
# ╠═bedf7419-590c-47fd-8089-b890dded0468
# ╠═8f795d92-b842-4392-8bc3-956080fcc6a8
# ╠═bc2e3724-de5f-4455-a4be-019c6a8accae
# ╠═90e9da81-7576-43bb-b4ce-1293a48bf3ec
# ╠═465e4ea9-7066-4929-98c4-e5c378ef34d8
# ╠═1129972e-3e3d-4b3b-8c2c-3c57eaa3728a
# ╠═79802166-efbd-4944-ae58-fb26a5d9ad4b
# ╠═146855a0-5cda-4783-89af-7918c067e952
# ╠═4e8fea9e-e16c-4b17-8a6c-52bd8721cfd8
# ╠═0588cc16-a4b1-4b3f-8ba9-742680996e11
# ╠═6893ccdb-7f6c-4ec9-8d6b-520e784acd03
# ╠═60a6fb20-93d2-47d3-8ee1-3087d48c569d
# ╠═0739d816-0456-4cc2-beb7-2e926cf2e229
# ╠═f829fc66-3df8-4ac1-a1bb-8eb2ed722423
# ╠═d4d6f69c-844b-488a-9170-b09e7c22dd25
# ╠═2152a72b-31cd-44b9-90d8-615610e8d261
# ╠═600dd2e9-aed9-41da-a510-9ba78788f207
