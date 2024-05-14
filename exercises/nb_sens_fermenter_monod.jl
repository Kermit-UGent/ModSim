### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 3ef93246-657d-4e77-9bf0-8380c64bfcfd
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 55cdebd2-0881-11ef-2722-91de1447877a
using Markdown

# ╔═╡ 03ae0690-06a0-4276-9f00-d07b206fe124
using InteractiveUtils

# ╔═╡ a355b0ba-baaf-49f4-a5dc-965364a884f0
using Catalyst

# ╔═╡ 00fd6d49-f561-42e9-9413-d33af92f83dc
using DifferentialEquations, Plots

# ╔═╡ 7ae714c4-d25d-4f9f-ab3d-cc067db9c156
using SciMLSensitivity

# ╔═╡ 31d294d1-3a1f-41db-abff-54f2a67c7ed9
md"
### Exercise: Fermenter - Monod kinetics - Sensitivity analysis
"

# ╔═╡ 5ffe7dcb-620d-4f22-95fe-2f77cda6fbe7
md"
In one of the previous practica we were introduced to a fermenter in which biomass $X$ [$g/L$] grows by breaking down substrate $S$ [$g/L$]. The reactor is fed with a inlet flow rate $Q_{in}$ [$L/h$], which consist of a (manipulable) input concentration of substrate $S_{in}$ [$g/L$]. This process was modelled using Monod kinetics, resulting in the model below:

$$\begin{eqnarray*}
%S \xrightarrow[\quad\quad]{\beta} Y \, X
S \xrightarrow[\quad\quad]{r} Y \, X \quad\quad\quad\quad r = \mu \, X
\end{eqnarray*}$$

where

$$\mu = \mu_{max} \, \cfrac{S}{S + K_s}$$
"

# ╔═╡ 6ec6da23-853b-4129-94cf-67b5cadb1f95
md"
The *reaction network object* model for this problem could be defined as:
"

# ╔═╡ 935ca610-7a7a-4692-8908-fc26abb880b4
fermenter_monod = @reaction_network begin
    X * mm(S, μmax, Ks), S --> Y*X
    Q/V, (S, X) --> ∅
    Q/V*Sin, ∅ --> S
end

# ╔═╡ 79e6056a-881c-442f-8989-5bc284d3d777
md"
which resulted in the following differential equations:
"

# ╔═╡ fa93e2c3-8b43-418e-ba24-406645b2e397
md"
$$\begin{eqnarray*}
\cfrac{dS}{dt} &=& \cfrac{Q}{V} \left(S_{in} - S \right) - \mu_{max}\cfrac{S}{S + K_s} X\\
\cfrac{dX}{dt} &=& -\cfrac{Q}{V} X + Y \mu_{max}\cfrac{S}{S + K_s} X
\end{eqnarray*}$$
"

# ╔═╡ 911b22bf-4cec-455d-a7f7-967bb55afea9
md"
Check out the order of the parameters:
"

# ╔═╡ be15ae00-163c-44e5-bc33-e939ec63ed05
# missing                           # Uncomment and complete the instruction
parameters(fermenter_monod)

# ╔═╡ 55f1d688-0c53-481b-9965-5e92ca87ad83
md"
The parameter values are $\mu_{max} = 0.30$, $K_s = 0.15$, $Y = 0.80$, $Q = 2.0$, $V = 40.0$ and $S_{in} = 2.2\;g/L$. Suppose that at $t$ no substrate $S$ is present in the reactor but that there is initially some biomass with a concetration of $0.01\;g/L$. Compute and plot the following in a timespan of $[0, 100]\,h$:

- The sensitivity of $S$ and $X$ to $\mu_{max}$.
- The sensitivity of $S$ and $X$ to $K_s$.
- The sensitivity of $S$ and $X$ to $S_{in}$.

Interpret your results.
"

# ╔═╡ 6c4e3c09-4b84-4f5c-8739-2ac18e6f2af6
md"
Initialize a vector `u₀` with the initial conditions, set the timespan and initialize a vector `param` with the parameter values:
"

# ╔═╡ 2ee277e5-ce4a-4ade-be0e-9bba7a4dc08c
# u₀ = missing                 # Uncomment and complete the instruction
u₀ = [:S => 0.0, :X => 0.01]

# ╔═╡ 3fdc6b17-cdeb-4dc5-8886-9d3a62caac8d
# tspan = missing              # Uncomment and complete the instruction
tspan = (0.0, 100.0)

# ╔═╡ 79b0eb65-5a0f-40b3-aa97-4088421c562e
# params = missing             # Uncomment and complete the instruction
params = [:μmax => 0.30, :Ks => 0.15, :Y => 0.80, :Q => 2, :V => 40, :Sin => 2.2]

# ╔═╡ f0f4fa14-6f99-4f21-a743-be61e08444a7
md"
Create the ODE problem and store it in `oprob`. Next, solve the ODE problem using `Tsit5()` and `saveat=0.5`, and store the solution in `osol`. Finally plot the results.
"

# ╔═╡ 8b2f23f6-80b2-4e63-942e-e5cd17d8ba72
# oprob = missing              # Uncomment and complete the instruction
oprob = ODEProblem(fermenter_monod, u₀, tspan, params)

# ╔═╡ 89a31c32-88a4-479f-a688-ffcb75ee8e91
# osol = missing                # Uncomment and complete the instruction
osol = solve(oprob, Tsit5(), saveat=0.5)

# ╔═╡ 51a9b7e6-8ad9-477d-9596-ffd614df2c79
# missing                        # Uncomment and complete the instruction
plot(osol)

# ╔═╡ cb484c24-501d-4e8c-81e5-6deb259ee948
md"
Create the ODE forward sensitivity problem and store it in `oprob_sens`:
"

# ╔═╡ 6567e6de-0067-4b7a-96f4-8ff773163b9c
# oprob_sens = missing           # Uncomment and complete the instruction
oprob_sens = ODEForwardSensitivityProblem(oprob.f, [0.0, 0.01], tspan, [0.3, 0.15, 0.80, 2.0, 40.0, 2.2])

# ╔═╡ 57052e7c-a93a-4cb3-a1d7-c8edcb2d438a
md"
Solve the *ODE forward sensitivity problem* using `Tsit5()` and `saveat=0.5`, and store the solution in `osol_sens`:
"

# ╔═╡ b4fa0025-b8a0-427d-ab96-39bb7833e67b
# osol_sens = missing            # Uncomment and complete the instruction
osol_sens = solve(oprob_sens, Tsit5(), saveat=0.5)

# ╔═╡ cf48664e-8e99-4c6d-8a92-8ac1f6e0b499
md"
Extract the sensitivity functions. Store the simulation results in the variable `u` and the sensitivities in the variable `dp`:
"

# ╔═╡ 0e66bf05-5924-4cff-bf11-ffc137ef4592
# u, dp = missing                # Uncomment and complete the instruction
u, dp = extract_local_sensitivities(osol_sens)

# ╔═╡ 97c3b8e7-31d5-4bad-8fc8-9c29641c449b
md"
Select by indexing and assigment the sensitivities of $S$ and $X$ to $\mu_{max}$, $K_s$ and $S_{in}$, to the variables `sens_μmax`, `sens_Ks` and `sens_Sin` respectively. Don't forget to transpose where necessary.

**Remark:**
- each of the elements `pd[i]` will contain two sensitivity functions (one of $S$ and one of $X$) to the `i`-th parameter. You can find the parameter indices of $\mu_{max}$, $K_s$ and $S_{in}$ by calling the function `parameters` in conjunction with the *reaction network* model name.
"

# ╔═╡ 2579de54-85d0-450e-8bbf-d116ea0b4ea0
# Absolute sensitivity of S and X to μmax
# sens_μmax = missing              # Uncomment and complete the instruction
sens_μmax = dp[1]'

# ╔═╡ ea539e22-2c83-4032-805f-d75e04939ce9
# Absolute sensitivity of S and X to Ks
# sens_Ks   = missing              # Uncomment and complete the instruction
sens_Ks   = dp[2]'

# ╔═╡ 53bae966-90db-47d5-86ba-6664b4dfbc5f
# Absolute sensitivity of S and X to Sin
# sens_Sin  = missing              # Uncomment and complete the instruction
sens_Sin  = dp[6]'

# ╔═╡ 2fdbbd60-c2ca-4d29-935f-252dee50ca4f
md"
Plot both sensitivity functions in separate notebook cells (with appropriate title and label):
"

# ╔═╡ 0144d423-cd18-4604-bf61-3bbf479717b4
# Absolute sensitivity of S and X to μmax
#  missing                         # Uncomment and complete the instruction
plot(osol_sens.t, sens_μmax, title="Absolute sensitivity of S and X to μmax", label=["dS/dμmax" "dX/dμmax"])

# ╔═╡ 98100622-9d84-413a-8e4b-49e69a07f5c0
# Absolute sensitivity of S and X to Ks
# missing                          # Uncomment and complete the instruction
plot(osol_sens.t, sens_Ks, title="Absolute sensitivity of S and X to Ks", label=["dS/dKs" "dX/dKs"])


# ╔═╡ 85519f77-b106-4e22-bbd8-438fb2b8eb7c
# Absolute sensitivity of S and X to Sin
# missing                          # Uncomment and complete the instruction
plot(osol_sens.t, sens_Sin, title="Absolute sensitivity for S and X on Ks", label=["dS/dSin" "dX/dSin"])


# ╔═╡ 5aaecc5b-6303-4050-813e-d2185b6e686b
md"
Draw your conclusions.
"

# ╔═╡ Cell order:
# ╠═55cdebd2-0881-11ef-2722-91de1447877a
# ╠═03ae0690-06a0-4276-9f00-d07b206fe124
# ╠═3ef93246-657d-4e77-9bf0-8380c64bfcfd
# ╠═a355b0ba-baaf-49f4-a5dc-965364a884f0
# ╠═00fd6d49-f561-42e9-9413-d33af92f83dc
# ╠═7ae714c4-d25d-4f9f-ab3d-cc067db9c156
# ╠═31d294d1-3a1f-41db-abff-54f2a67c7ed9
# ╠═5ffe7dcb-620d-4f22-95fe-2f77cda6fbe7
# ╠═6ec6da23-853b-4129-94cf-67b5cadb1f95
# ╠═935ca610-7a7a-4692-8908-fc26abb880b4
# ╠═79e6056a-881c-442f-8989-5bc284d3d777
# ╠═fa93e2c3-8b43-418e-ba24-406645b2e397
# ╠═911b22bf-4cec-455d-a7f7-967bb55afea9
# ╠═be15ae00-163c-44e5-bc33-e939ec63ed05
# ╠═55f1d688-0c53-481b-9965-5e92ca87ad83
# ╠═6c4e3c09-4b84-4f5c-8739-2ac18e6f2af6
# ╠═2ee277e5-ce4a-4ade-be0e-9bba7a4dc08c
# ╠═3fdc6b17-cdeb-4dc5-8886-9d3a62caac8d
# ╠═79b0eb65-5a0f-40b3-aa97-4088421c562e
# ╠═f0f4fa14-6f99-4f21-a743-be61e08444a7
# ╠═8b2f23f6-80b2-4e63-942e-e5cd17d8ba72
# ╠═89a31c32-88a4-479f-a688-ffcb75ee8e91
# ╠═51a9b7e6-8ad9-477d-9596-ffd614df2c79
# ╠═cb484c24-501d-4e8c-81e5-6deb259ee948
# ╠═6567e6de-0067-4b7a-96f4-8ff773163b9c
# ╠═57052e7c-a93a-4cb3-a1d7-c8edcb2d438a
# ╠═b4fa0025-b8a0-427d-ab96-39bb7833e67b
# ╠═cf48664e-8e99-4c6d-8a92-8ac1f6e0b499
# ╠═0e66bf05-5924-4cff-bf11-ffc137ef4592
# ╠═97c3b8e7-31d5-4bad-8fc8-9c29641c449b
# ╠═2579de54-85d0-450e-8bbf-d116ea0b4ea0
# ╠═ea539e22-2c83-4032-805f-d75e04939ce9
# ╠═53bae966-90db-47d5-86ba-6664b4dfbc5f
# ╠═2fdbbd60-c2ca-4d29-935f-252dee50ca4f
# ╠═0144d423-cd18-4604-bf61-3bbf479717b4
# ╠═98100622-9d84-413a-8e4b-49e69a07f5c0
# ╠═85519f77-b106-4e22-bbd8-438fb2b8eb7c
# ╠═5aaecc5b-6303-4050-813e-d2185b6e686b
