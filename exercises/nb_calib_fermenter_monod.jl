### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ c54dae10-60af-4141-b56d-ed61cb0ced8a
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 245ca9d0-10f9-11ef-0ef6-a73594e96db9
using Markdown

# ╔═╡ 78a25bef-31e5-45ef-b0ba-b9a8c8a9edeb
using InteractiveUtils

# ╔═╡ 16438e07-1b2b-467e-822a-081d19cae92b
using Catalyst, DifferentialEquations, Plots

# ╔═╡ 295caa68-db27-4c9b-bc34-86ab088fec24
using Turing

# ╔═╡ a0269a25-269c-49c1-aa05-047fa3a959d4
using StatsPlots, StatsBase

# ╔═╡ dc6e9bdc-dae0-43aa-b624-f9314d1d9884
using LinearAlgebra

# ╔═╡ 31243ea7-1f0f-490f-8886-b1b7ab7ae5b4
using Optim

# ╔═╡ 2f0a4c62-3441-4c63-9bb9-383e7f554eb5
md"""
### Exercise: Fermenter - Monod kinetics - Calibration
"""

# ╔═╡ 595ea8ee-bc67-4696-9232-982612fb554d
md"""
In one of the previous practica we were introduced to a fermenter in which biomass $X$ [$g/L$] grows by breaking down substrate $S$ [$g/L$]. The reactor is fed with a inlet flow rate $Q_{in}$ [$L/h$], which consist of a (manipulable) input concentration of substrate $S_{in}$ [$g/L$]. This process was modelled using Monod kinetics, resulting in the model below:

$$\begin{eqnarray*}
%S \xrightarrow[\quad\quad]{\beta} Y \, X
S \xrightarrow[\quad\quad]{r} Y \, X \quad\quad\quad\quad r = \mu \, X \quad \textrm{with} \quad \mu = \mu_{max} \, \cfrac{S}{S + K_s}
\end{eqnarray*}$$
"""

# ╔═╡ 824db995-7a66-4719-a534-7e0f6dec90b5
md"""
The *reaction network object* for this model could be set-up as:
"""

# ╔═╡ 245c2636-95da-4c76-8b03-c4d20bbabb48
fermenter_monod = @reaction_network begin
	X * mm(S, μmax, Ks), S => Y*X
	# Alternative:
	# X * mm(S, μmax, Ks), S --> Y*X
    Q/V, (S, X) --> 0
    Q/V*Sin, 0 --> S
end

# ╔═╡ 68f9ecb3-15b0-4a53-8864-5dac13a89e95
parameters(fermenter_monod)

# ╔═╡ de8ddc14-8f82-403d-8f42-29673ef2a722
md"""
which resulted in the following differential equations:

$$\begin{eqnarray*}
\cfrac{dS}{dt} &=& \cfrac{Q}{V} \left(S_{in} - S \right) - \mu_{max}\cfrac{S}{S + K_s} X\\
\cfrac{dX}{dt} &=& -\cfrac{Q}{V} X + Y \mu_{max}\cfrac{S}{S + K_s} X
\end{eqnarray*}$$
"""

# ╔═╡ b7b7d58f-d406-4596-b834-ced6d8fada83
md"""
Suppose that during an experiment measurement data has been collected of the substrate $S$ and biomass $X$ concentration at an interval of $5\;h$ within $100\;h$:
"""

# ╔═╡ 99c6f31a-0968-4804-9980-71fcc1af1f49
S_meas = [1.0e-5, 0.0047, 0.00796, 0.01056, 0.01214, 0.01325, 0.01344, 0.01338, 0.0115, 0.00917, 0.00604, 0.00458, 0.00438, 0.00342, 0.00323, 0.00329, 0.00312, 0.00314, 0.00319, 0.00299, 0.00311]

# ╔═╡ bf4ad873-e0fe-415c-9e78-fe0b5ac1414e
X_meas = [0.00052, 0.00042, 0.00074, 0.00078, 0.00122, 0.00159, 0.00242, 0.00372, 0.00534, 0.0077, 0.00935, 0.00997, 0.01114, 0.01144, 0.01264, 0.01276, 0.01183, 0.01319, 0.01256, 0.01277, 0.01377]

# ╔═╡ 1dae5875-f405-4ecb-8b7b-3c3f22b549bb
t_meas = 0.0:5.0:100.0

# ╔═╡ 6c481447-28c6-4530-bf2c-64762121bc71
md"""
Make a scatter plot of the measured data for both $S$ and $X$. Use the following options:
- `label=\"S meas\", color=:blue` for $S$, and
- `label=\"X meas\", color=:red` for $X$.
"""

# ╔═╡ 918fd524-81fa-4aff-a403-37402e47235b
# Uncomment and complete the instruction
# begin
# 	missing
#   missing
# end
begin
	scatter(t_meas, S_meas, label="S meas", color=:blue)
    scatter!(t_meas, X_meas, label="X meas", color=:red) # notice exclamation mark !
end

# ╔═╡ ef977370-06ee-4a73-85e2-609a744167d3
md"
We have previously used the following parameter values:

-  $\mu_{max} = 0.40\;h^{-1}$, $K_s = 0.015\;g/L$, $S_{in} = 0.22\;g/L$
-  $Y = 0.67$, $Q = 2.0\;L/h$, $V = 40.0\;L$

Furthermore, suppose that at $t = 0\;h$ no substrate $S$ is present in the reactor but that there is initially some biomass with a concetration of $0.0005\;g/L$.

Calibrate the parameter values for $\mu_{max}$ and $K_s$ using the aforementioned measurement data for $S$ and $X$ in a timespan of $[0, 100]\,h$. Take the values above as initial values for $\mu_{max}$ and $K_s$.
"

# ╔═╡ 7a227eaf-18d0-44f4-ac4b-f529e81c7471
md"""
Create an `ODEProblem`. Use the aforementioned values as initial values for the problem.
"""

# ╔═╡ 6375478f-1af9-4fd2-b6f3-101a6f796f2d
# u0 = missing      # Uncomment and complete the instruction
u0 = [:S => 0.0, :X => 0.0005]

# ╔═╡ 38fe8304-af61-40a7-ac86-480dfb892185
# tspan = missing   # Uncomment and complete the instruction
tspan = (0.0, 100)

# ╔═╡ 87482f88-8413-4820-9613-7941f3d61bd7
# params = missing  # Uncomment and complete the instruction
params = [:μmax => 0.40, :Ks => 0.015, :Y => 0.67, :Q => 2, :V => 40, :Sin => 0.022]

# ╔═╡ 94f3bd7b-5c2c-4661-a0ab-2cdaf2cd6743
# oprob = missing   # Uncomment and complete the instruction
oprob = ODEProblem(fermenter_monod, u0, tspan, params, combinatoric_ratelaws=false)

# ╔═╡ f6a8f134-6db0-4d74-8af5-82826347d8f0
md"""
Declare the Turing model.
"""

# ╔═╡ 4c28a66a-ee2c-42a2-95c7-ea4ddb6a232d
# Uncomment and complete the instruction
# @model function fermenter_inference(t_meas, S, X)
    # σ_S ~ missing
    # σ_X ~ missing
    # μmax ~ missing
    # Ks ~ missing
	# params = missing
	# oprob = missing
    # osol = missing
    # S ~ missing
    # X ~ missing
# end
@model function fermenter_inference(t_meas, S, X)
	σ_S ~ InverseGamma()
	σ_X ~ InverseGamma()
    μmax ~ LogNormal()
	Ks ~ LogNormal()
    # μmax ~ Normal(0.40, 0.1)
	# Ks ~ Normal(0.015, 0.01)
    params = [:μmax => μmax, :Ks => Ks, :Y => 0.67, :Q => 2, :V => 40, :Sin => 0.022]
	oprob = ODEProblem(fermenter_monod, u0, tspan, params)
	# , combinatoric_ratelaws=false
    osol = solve(oprob, Tsit5(), saveat=t_meas)
	S ~ MvNormal(osol[:S], σ_S^2 * I)
	X ~ MvNormal(osol[:X], σ_X^2 * I)
end

# ╔═╡ 3136b15d-5078-4bcd-954b-e89bcb8aed1b
md"""
Provide the measurements to the Turing model.
"""

# ╔═╡ 6a508a62-61b9-4273-8e45-b26f594e8da9
# fermenter_inf = missing       # Uncomment and complete the instruction
fermenter_inf = fermenter_inference(t_meas, S_meas, X_meas)

# ╔═╡ 63420055-55f8-4def-8b0e-11ea61483010
md"""
Optimize the priors ($\sigma_S$, $\sigma_X$, $\mu_{max}$ and $K_s$). Do this with `MLE` method and Nelder-Mead. Store the optimization results in `results_mle`.
"""

# ╔═╡ d52c9da8-d8a4-4db0-ac6d-6d16ccf4775c
# results_map = missing           # Uncomment and complete the instruction
results_mle = optimize(fermenter_inf, MLE(), NelderMead())

# ╔═╡ e1b0ee01-f16c-40e9-a0f9-80072d690936
md"""
Visualize a summary of the optimized parameters.
"""

# ╔═╡ f2d7daf8-8218-446d-b1d2-e9e05aeadfd9
# missing        # Uncomment and complete the instruction
results_mle |> coeftable

# ╔═╡ 23d58bb1-d077-402e-8bee-3866c68e069a
md"""
Get the optimized values and assign them to `μmax_opt`, `Ks_opt` and `Sin_opt`.
"""

# ╔═╡ 7b3a3677-b251-43c1-b125-6d6ff1a11ea3
# μmax_opt = missing              # Uncomment and complete the instruction
μmax_opt = coef(results_mle)[:μmax]

# ╔═╡ fa77bcbe-2ddc-4113-8f6a-4a18d219da9e
# Ks_opt = missing                # Uncomment and complete the instruction
Ks_opt = coef(results_mle)[:Ks]

# ╔═╡ 05d13a48-adc8-4e24-a6e4-be24af2c7a59
md"""
Make a plot of $S$ and $X$ simulated with the optimized parameter values.
"""

# ╔═╡ 57ef3824-7c20-4876-8d83-665cb4f97a58
md"""
Set up parameter values with optimized parameter values:
"""

# ╔═╡ 75cf59ed-af8e-4e8a-8ed2-1f3bf4d386d0
# params_opt = missing         # Uncomment and complete the instruction
params_opt = [:μmax => μmax_opt, :Ks => Ks_opt, :Y => 0.67, :Q => 2, :V => 40, :Sin => 0.022]

# ╔═╡ 4e8870dc-2da6-4b80-82d6-26c7ceedad7d
md"""
Create an ODEProblem and solve it:
"""

# ╔═╡ 853c1a92-d50f-4b05-9ed3-d3ee1656665a
# oprob_opt = missing         # Uncomment and complete the instruction
oprob_opt = ODEProblem(fermenter_monod, u0, tspan, params_opt)
# , combinatoric_ratelaws=false

# ╔═╡ f45e8124-e942-438e-99c5-3032ccc01454
# osol_opt = missing          # Uncomment and complete the instruction
osol_opt = solve(oprob_opt, Tsit5(), saveat=0.5)

# ╔═╡ 5a39b0e0-1ea1-4854-8e68-66d0d4bbf25c
md"""
Plot $S$ and $X$ simulated with the optimized parameter values together with the measured data.
"""

# ╔═╡ d0156099-ad03-4711-ac0f-94882fb78266
# Uncomment and complete the instruction
# begin
#   missing
#   missing
#   missing
# end
begin
  plot(osol_opt, labels=["S sim" "X sim"], xlabel="t")
  scatter!(t_meas, S_meas, label="S meas", color=:blue)
  scatter!(t_meas, X_meas, label="X meas", color=:red)
end

# ╔═╡ Cell order:
# ╠═245ca9d0-10f9-11ef-0ef6-a73594e96db9
# ╠═78a25bef-31e5-45ef-b0ba-b9a8c8a9edeb
# ╠═c54dae10-60af-4141-b56d-ed61cb0ced8a
# ╠═16438e07-1b2b-467e-822a-081d19cae92b
# ╠═295caa68-db27-4c9b-bc34-86ab088fec24
# ╠═a0269a25-269c-49c1-aa05-047fa3a959d4
# ╠═dc6e9bdc-dae0-43aa-b624-f9314d1d9884
# ╠═31243ea7-1f0f-490f-8886-b1b7ab7ae5b4
# ╟─2f0a4c62-3441-4c63-9bb9-383e7f554eb5
# ╟─595ea8ee-bc67-4696-9232-982612fb554d
# ╟─824db995-7a66-4719-a534-7e0f6dec90b5
# ╠═245c2636-95da-4c76-8b03-c4d20bbabb48
# ╠═68f9ecb3-15b0-4a53-8864-5dac13a89e95
# ╟─de8ddc14-8f82-403d-8f42-29673ef2a722
# ╠═b7b7d58f-d406-4596-b834-ced6d8fada83
# ╠═99c6f31a-0968-4804-9980-71fcc1af1f49
# ╠═bf4ad873-e0fe-415c-9e78-fe0b5ac1414e
# ╠═1dae5875-f405-4ecb-8b7b-3c3f22b549bb
# ╟─6c481447-28c6-4530-bf2c-64762121bc71
# ╠═918fd524-81fa-4aff-a403-37402e47235b
# ╟─ef977370-06ee-4a73-85e2-609a744167d3
# ╟─7a227eaf-18d0-44f4-ac4b-f529e81c7471
# ╠═6375478f-1af9-4fd2-b6f3-101a6f796f2d
# ╠═38fe8304-af61-40a7-ac86-480dfb892185
# ╠═87482f88-8413-4820-9613-7941f3d61bd7
# ╠═94f3bd7b-5c2c-4661-a0ab-2cdaf2cd6743
# ╟─f6a8f134-6db0-4d74-8af5-82826347d8f0
# ╠═4c28a66a-ee2c-42a2-95c7-ea4ddb6a232d
# ╟─3136b15d-5078-4bcd-954b-e89bcb8aed1b
# ╠═6a508a62-61b9-4273-8e45-b26f594e8da9
# ╟─63420055-55f8-4def-8b0e-11ea61483010
# ╠═d52c9da8-d8a4-4db0-ac6d-6d16ccf4775c
# ╟─e1b0ee01-f16c-40e9-a0f9-80072d690936
# ╠═f2d7daf8-8218-446d-b1d2-e9e05aeadfd9
# ╟─23d58bb1-d077-402e-8bee-3866c68e069a
# ╠═7b3a3677-b251-43c1-b125-6d6ff1a11ea3
# ╠═fa77bcbe-2ddc-4113-8f6a-4a18d219da9e
# ╟─05d13a48-adc8-4e24-a6e4-be24af2c7a59
# ╟─57ef3824-7c20-4876-8d83-665cb4f97a58
# ╠═75cf59ed-af8e-4e8a-8ed2-1f3bf4d386d0
# ╟─4e8870dc-2da6-4b80-82d6-26c7ceedad7d
# ╠═853c1a92-d50f-4b05-9ed3-d3ee1656665a
# ╠═f45e8124-e942-438e-99c5-3032ccc01454
# ╟─5a39b0e0-1ea1-4854-8e68-66d0d4bbf25c
# ╠═d0156099-ad03-4711-ac0f-94882fb78266
