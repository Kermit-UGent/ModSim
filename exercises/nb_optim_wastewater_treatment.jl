### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 2b26c3b2-df08-4b24-a08d-23717248c10d
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ f08fa69c-a744-11ef-0e79-3daf5bf297ea
using Markdown

# ╔═╡ e6c97e37-d062-4b65-96a4-bac0dab220d8
using InteractiveUtils

# ╔═╡ 0407d891-a46d-4deb-a21a-23833acbcb87
using Catalyst, DifferentialEquations, Plots

# ╔═╡ e48dc930-be03-47b2-b9e3-16e854782aec
using Turing

# ╔═╡ 39cc0864-7399-431b-bfc2-70196866c41b
using StatsBase

# ╔═╡ 9d23f4bf-c3bd-4e79-923c-4b3b77b806a3
using Optim

# ╔═╡ c1bc698d-41ee-45e6-b17d-29f0d53557a1
md"""
#### Exercise: Wastewater treatment - Optimisation
"""

# ╔═╡ dfe77a8c-a8db-46d7-9a4d-3b00413d383b
md"""
Consider a wastewater treatment plant where wastewater circulates in cylindrical tanks so that microorganisms can break down the organic material present. At the top of such a tank with volume $V\;[m^3]$, wastewater enters at a flow rate $q\;[m^3/h]$. The concentration of organic material in the inflow is known and equal to $C_{in}\;[kg/m^3]$. At the bottom of the tank, wastewater and microorganisms leave the tank at the same flow rate $q\;[m^3/h]$ so that the volume of wastewater in the tank remains constant.

The concentration of organic material in the tank is denoted as $C\;[kg/m^3]$ and the concentration of microorganisms is denoted as $X\;[kg/m^3]$. The microorganisms in the tank break down the organic material at a rate proportional to $r\cfrac{K_s}{K_s+C}\;[m^3\,h^{-1}\,kg^{-1}]$ with yield coefficient $Y$. $K_s\;[kg/m^3]$ is the concentration of $C$ where the rate is half its maximum rate and $r\;[m^3\,h^{-1}\,kg^{-1}]$ is the maximum growth rate coefficient. Furthermore, the microorganisms degrade with rate coefficient $k_d$. In the middle of the tank there is a mixing system that ensures that the wastewater and microorganisms are well mixed. This means that the concentration in the outflow is equal to the concentration in the tank: $C_{out} = C$ and $X_{out} = X$. The system of differential equations describing the change in the concentrations $C(t)$ and $X(t)$ is given by:

$$\cfrac{dC}{dt} = \cfrac{q}{V}\left(C_{in} - C\right) - r\cfrac{K_s}{K_s+C}\,C\,X$$
$$\cfrac{dX}{dt} = -\cfrac{q}{V}X -k_d\,X + Y\,r\cfrac{K_s}{K_s+C}\,C\,X$$

The initial concentrations and the parameter values are summarised in the following tables:

|   $C_0$   |   $X_0$   |
|:---------:|:---------:|
|   $4.0$  |   $0.01$  |

|   $q$    |   $V$    |   $r$    | $C_{in}$  |   $K_s$   |   $k_d$   |   $Y$    |
|:--------:|:--------:|:--------:|:---------:|:---------:|:---------:|:--------:|
|   $5.0$  |   $50$   |  $0.4$   |   $3.0$   |   $5.2$   |  $0.10$   |   $1.2$  |

The amount of organic waste that is being breakdown by the microorganisms depends on the flow rate $q$. First (Part 1), we will simulate the system with the parameters given above. Second (Part 2), we will optimize the value of the flow rate $q$ so that the concentration of organic waste in the tank is at most $0.28\; kg\,m^{-3}$.
"""

# ╔═╡ c0c83df7-a9cc-4bde-b6ec-038a423b0d90
md"""
#### Part 1
"""

# ╔═╡ bbd50cc5-032a-4219-bcee-91145935a7c4
md"""
##### Implementation of the system
"""

# ╔═╡ a5f79c20-f62e-4df2-be79-b4f2141ced5e
md"""
Create a *reaction network object* model for the aforementioned problem. Name it `wastewater_treatment`.

Tips:
- You can use the repressive Michaelis-Menten function `mmr(C, r, Ks)` for $r\cfrac{K_s}{K_s+C}$.

"""

# ╔═╡ 8708de16-3532-4352-b211-c092f95c82d3
wastewater_treatment = @reaction_network begin
    @parameters q=5.0 V=50 r=0.4 Cin=3.0 Ks=5.2 kd=0.10 Y=1.2
    @species C(t)=4 X(t)=1.0e-2
    # mmr(C, r, Ks), C + X --> (1+Y)*X
    # r*Ks/(C+Ks), C + X --> (1+Y)*X
	mmr(C, r, Ks)*X, C --> Y*X
    q/V*Cin, 0 --> C
    q/V, C --> 0
    q/V, X --> 0
    kd, X --> 0
end

# ╔═╡ 10c73294-a32b-4aa3-80a8-10785e5eab8f
md"""
Convert the system to a symbolic differential equation model and verify your system of differential equations.
"""

# ╔═╡ fee917dd-7ab5-4fda-b1b7-87ee61e21f19
osys  = convert(ODESystem, wastewater_treatment)

# ╔═╡ 08ebcb95-8603-4579-879e-810b1494b013
md"""
##### Setting up initial conditions, timespan and parameter values
"""

# ╔═╡ d248f64e-ebba-4443-9c49-ff0290aa7810
md"""
Initialize a vector `u0` with the initial conditions:
"""

# ╔═╡ fe02a755-5b00-4d80-a511-fec115b42964
u0 = [:C=>3.0, :X=>0.5]

# ╔═╡ 67481927-0d03-4da9-af6c-9afa409fc006
md"""
Set the timespan:
"""

# ╔═╡ fadd372a-a665-4b16-9b6d-e32cb7f25d7f
tspan = (0.0, 72.0)

# ╔═╡ 734e4d51-95a7-464e-9a23-5ad6c8715d65
md"""
Initialize a vector `params` with the parameter values:
"""

# ╔═╡ 15ce9889-a437-46c8-9062-74b8d234a8bd
params = [:q=>5, :V=>50, :r=>0.4, :Cin=>3.0, :Ks=>5.2, :kd=>0.10, :Y=>1.2]

# ╔═╡ b8a48461-3882-45f6-980c-38d650ac52c7
md"""
##### Creating an ODE problem, solve the problem and plot results
"""

# ╔═╡ a4c57b64-6a7d-4bd4-8bb2-578923e184d2
md"""
Create the ODE problem and store it in `oprob`:
"""

# ╔═╡ b1e18139-5277-4f06-b1f8-b0f5f11c41d8
oprob = ODEProblem(wastewater_treatment, u0, tspan, params)

# ╔═╡ ad6d8fe6-e62f-4c67-8d63-4ee13b928ad0
md"""
Solve the ODE problem. Use `Tsit5()` and `saveat=0.1`. Store the solution in `osol`:
"""

# ╔═╡ e2ffba9e-aaf2-4540-84cf-8b7297ae9285
osol = solve(oprob, Tsit5(), saveat=0.1)

# ╔═╡ 70871ee8-b0a4-4a9a-af39-5a63459b55f7
md"""
Plot the results.
"""

# ╔═╡ 34309734-3751-47e0-a602-d113ffaae510
begin
    plot(osol, ylim=(0, 4), linewidth=2)
    plot!([tspan[1], tspan[2]], [0.28, 0.28], linestyle=:dash, linewidth=2, linecolor=:green, label="")
end

# ╔═╡ a0e735ad-09c2-4aa8-bc41-b294a9d56ea8
md"""
Check out the end value of the organic waste.
"""

# ╔═╡ f62898d5-1b8d-4350-8655-78aa3decb2a2
osol[:C][end]

# ╔═╡ 7eb5df9c-a475-4812-81c3-e43484c82242
md"""
#### Part 2
"""

# ╔═╡ b6bac48a-4a3d-47e4-90ea-788ca20dadff
@model function wastewater_treatment_inference(C_val)
	q ~ Uniform(0, 5.0)
    u0 = [:C=>3.0, :X=>0.5]
    tspan = (0.0, 72.0)  # the time interval to solve on
    params = [:q=>q, :V=>50, :r=>0.4, :Cin=>3.0, :Ks=>5.2, :kd=>0.10, :Y=>1.2]
    oprob = ODEProblem(wastewater_treatment, u0, tspan, params)
    osol = solve(oprob, Tsit5(), saveat=0.1)
    C_val ~ Normal(osol[:C][end], 1e-3)
end

# ╔═╡ 2df409ef-bd95-4ac3-a2b8-c5e17c490eba
C_val = 0.28

# ╔═╡ afc035be-075b-464b-8ba2-20235082f005
results_mle = optimize(wastewater_treatment_inference(C_val), MLE(), NelderMead())

# ╔═╡ 98a157a1-8c20-474d-acb8-00373ee6d224
q_opt = coef(results_mle)[:q]

# ╔═╡ e275df05-5c77-4c17-ad2e-503574596c31
params_opt = [:q=>q_opt, :V=>50, :r=>0.4, :Cin=>3.0, :Ks=>5.2, :kd=>0.10, :Y=>1.2]

# ╔═╡ a4388f06-1223-4815-a557-9b9c3ec232bb
oprob_opt = ODEProblem(wastewater_treatment, u0, tspan, params_opt)

# ╔═╡ ec7bf654-b275-4cfd-a819-d82bdc1be93b
osol_opt = solve(oprob_opt, Tsit5(), saveat=0.1)

# ╔═╡ 81429279-4190-41d1-a72a-20da0ce90528
begin
    plot(osol_opt, ylim=(0, 4), linewidth=2)
    plot!([tspan[1], tspan[2]], [C_val, C_val], linestyle=:dash, linewidth=2, linecolor=:green)
end

# ╔═╡ Cell order:
# ╠═f08fa69c-a744-11ef-0e79-3daf5bf297ea
# ╠═e6c97e37-d062-4b65-96a4-bac0dab220d8
# ╠═2b26c3b2-df08-4b24-a08d-23717248c10d
# ╠═0407d891-a46d-4deb-a21a-23833acbcb87
# ╠═e48dc930-be03-47b2-b9e3-16e854782aec
# ╠═39cc0864-7399-431b-bfc2-70196866c41b
# ╠═9d23f4bf-c3bd-4e79-923c-4b3b77b806a3
# ╠═c1bc698d-41ee-45e6-b17d-29f0d53557a1
# ╠═dfe77a8c-a8db-46d7-9a4d-3b00413d383b
# ╠═c0c83df7-a9cc-4bde-b6ec-038a423b0d90
# ╠═bbd50cc5-032a-4219-bcee-91145935a7c4
# ╠═a5f79c20-f62e-4df2-be79-b4f2141ced5e
# ╠═8708de16-3532-4352-b211-c092f95c82d3
# ╠═10c73294-a32b-4aa3-80a8-10785e5eab8f
# ╠═fee917dd-7ab5-4fda-b1b7-87ee61e21f19
# ╠═08ebcb95-8603-4579-879e-810b1494b013
# ╠═d248f64e-ebba-4443-9c49-ff0290aa7810
# ╠═fe02a755-5b00-4d80-a511-fec115b42964
# ╠═67481927-0d03-4da9-af6c-9afa409fc006
# ╠═fadd372a-a665-4b16-9b6d-e32cb7f25d7f
# ╠═734e4d51-95a7-464e-9a23-5ad6c8715d65
# ╠═15ce9889-a437-46c8-9062-74b8d234a8bd
# ╠═b8a48461-3882-45f6-980c-38d650ac52c7
# ╠═a4c57b64-6a7d-4bd4-8bb2-578923e184d2
# ╠═b1e18139-5277-4f06-b1f8-b0f5f11c41d8
# ╠═ad6d8fe6-e62f-4c67-8d63-4ee13b928ad0
# ╠═e2ffba9e-aaf2-4540-84cf-8b7297ae9285
# ╠═70871ee8-b0a4-4a9a-af39-5a63459b55f7
# ╠═34309734-3751-47e0-a602-d113ffaae510
# ╠═a0e735ad-09c2-4aa8-bc41-b294a9d56ea8
# ╠═f62898d5-1b8d-4350-8655-78aa3decb2a2
# ╠═7eb5df9c-a475-4812-81c3-e43484c82242
# ╠═b6bac48a-4a3d-47e4-90ea-788ca20dadff
# ╠═2df409ef-bd95-4ac3-a2b8-c5e17c490eba
# ╠═afc035be-075b-464b-8ba2-20235082f005
# ╠═98a157a1-8c20-474d-acb8-00373ee6d224
# ╠═e275df05-5c77-4c17-ad2e-503574596c31
# ╠═a4388f06-1223-4815-a557-9b9c3ec232bb
# ╠═ec7bf654-b275-4cfd-a819-d82bdc1be93b
# ╠═81429279-4190-41d1-a72a-20da0ce90528
