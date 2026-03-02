### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 2b26c3b2-df08-4b24-a08d-23717248c10d
using Pkg; Pkg.activate("..")

# ╔═╡ f08fa69c-a744-11ef-0e79-3daf5bf297ea
using Markdown, InteractiveUtils

# ╔═╡ 0407d891-a46d-4deb-a21a-23833acbcb87
using ModelingToolkit, OrdinaryDiffEq

# ╔═╡ 120914e8-afd7-410d-9b57-01a8f24bdf59
using ModelingToolkit: t_nounits as t, D_nounits as D

# ╔═╡ e48dc930-be03-47b2-b9e3-16e854782aec
using Turing, StatsPlots, StatsBase, Optim

# ╔═╡ 6458329f-73dd-4cb0-8da4-90678875a1f0
using PlutoUI; TableOfContents()

# ╔═╡ c1bc698d-41ee-45e6-b17d-29f0d53557a1
md"""
# Exercise: Wastewater treatment - Optimisation
"""

# ╔═╡ dfe77a8c-a8db-46d7-9a4d-3b00413d383b
md"""
Consider a wastewater treatment plant where wastewater circulates through cylindrical tanks, allowing microorganisms to break down the organic material present. At the top of such a tank with volume $V\;[\mathrm{m^3}]$, wastewater enters at a flow rate $q\;[\mathrm{m^3/h}]$. The concentration of organic material in the inflow is known and equal to $C_{in}\;[\mathrm{kg/m^3}]$. At the bottom of the tank, wastewater and microorganisms leave the tank at the same flow rate $q\;[\mathrm{m^3/h}]$ so that the volume of wastewater in the tank remains constant.

The concentration of organic material in the tank is denoted as $C\;[\mathrm{kg/m^3}]$ and the concentration of microorganisms is denoted as $X\;[\mathrm{kg/m^3}]$. The microorganisms in the tank break down the organic material at a rate proportional to $r\cfrac{K_s}{K_s+C}\;[\mathrm{m^3\,h^{-1}\,kg^{-1}}]$ with yield coefficient $Y$. The factor $K_s\;[\mathrm{kg/m^3}]$ is the concentration of $C$ where the rate is half its maximum rate and $r\;[\mathrm{m^3\,h^{-1}\,kg^{-1}}]$ is the maximum growth rate coefficient. Furthermore, the microorganisms degrade with a rate coefficient $k_d$. In the middle of the tank, a mixing system ensures that wastewater and microorganisms are thoroughly mixed. This means that the concentration in the outflow is equal to the concentration in the tank: $C_{out} = C$ and $X_{out} = X$. The system of differential equations describing the change in the concentrations $C(t)$ and $X(t)$ is given by:

$$\cfrac{dC}{dt} = \cfrac{q}{V}\left(C_{in} - C\right) - r\cfrac{K_s}{K_s+C}\,C\,X$$
$$\cfrac{dX}{dt} = -\cfrac{q}{V}X -k_d\,X + Y\,r\cfrac{K_s}{K_s+C}\,C\,X$$

The initial concentrations and the parameter values are summarised in the following tables:

|   $C_0$   |   $X_0$   |
|:---------:|:---------:|
|   $3.0$  |   $0.5$  |

|   $q$    |   $V$    |   $r$    | $C_{in}$  |   $K_s$   |   $k_d$   |   $Y$    |
|:--------:|:--------:|:--------:|:---------:|:---------:|:---------:|:--------:|
|   $5.0$  |   $50$   |  $0.4$   |   $3.0$   |   $5.2$   |  $0.10$   |   $1.2$  |

The amount of organic waste being broken down by microorganisms depends on the flow rate $q$. First (Part 1), we will simulate the system with the parameters given above. Second (Part 2), we will optimize the value of the flow rate $q$ so that the concentration of organic waste in the tank is at most $0.28\;\mathrm{kg\,m^{-3}}$.
"""

# ╔═╡ c0c83df7-a9cc-4bde-b6ec-038a423b0d90
md"""
## Part 1

In this part, we will simulate the system with the parameters given above.
"""

# ╔═╡ bbd50cc5-032a-4219-bcee-91145935a7c4
md"""
### Implementation of the system
"""

# ╔═╡ 93d49247-8e37-4abc-b792-f380c8d59883
md"""
Model the system by means of ModelingToolkit.
"""

# ╔═╡ d48dae1d-a8a2-4b23-b4b8-c840be19ff66
md"""
Define the variables and assign them to their default values.
"""

# ╔═╡ 08d82f11-068a-4ab3-a664-46af49f1ff11
# @variables missing
@variables C(t)=3.0 X(t)=0.5

# ╔═╡ dc2189d6-ca6f-4edd-be41-49ea5a6d556e
md"""
Define the parameters and assign them to their default values.
"""

# ╔═╡ 8dc5c4fd-687f-4dd9-ab80-bf80e40af0df
# @parameters missing
@parameters q=5.0 V=50 r=0.4 Cin=3.0 Ks=5.2 kd=0.10 Y=1.2

# ╔═╡ 1b3c329c-4cd6-40cc-97bc-e8f829eb04d6
md"""
Set up the equations for the change in $C$ and $X$.
"""

# ╔═╡ 799f699c-1c94-406e-bd84-a345b96144c5
# change_C = missing
change_C = D(C) ~ (q/V)*(Cin - C) - r*Ks/(Ks + C)*C*X

# ╔═╡ 7a60bdbb-a5ac-40c2-9229-f8daa4cdcd92
# change_X = missing
change_X = D(X) ~ -(q/V)*X - kd*X + Y*r*Ks/(Ks + C)*C*X

# ╔═╡ 452734e6-5e51-4488-a0c2-88b86f310f24
md"""
Build the ODE system and name it `sys_ww_treat`.
"""

# ╔═╡ abd9fcba-568f-4617-a138-a74a06850e9a
# @mtkbuild missing
@mtkbuild sys_ww_treat = ODESystem([change_C, change_X], t)

# ╔═╡ 08ebcb95-8603-4579-879e-810b1494b013
md"""
##### Setting up initial conditions, timespan and parameter values
"""

# ╔═╡ d248f64e-ebba-4443-9c49-ff0290aa7810
md"""
Initialize a vector `u0` with the initial conditions:
"""

# ╔═╡ fe02a755-5b00-4d80-a511-fec115b42964
# u0 = missing
u0 = [:C=>3.0, :X=>0.5]

# ╔═╡ 67481927-0d03-4da9-af6c-9afa409fc006
md"""
Set the timespan:
"""

# ╔═╡ fadd372a-a665-4b16-9b6d-e32cb7f25d7f
# tspan = missing
tspan = (0.0, 72.0)

# ╔═╡ 734e4d51-95a7-464e-9a23-5ad6c8715d65
md"""
Initialize a vector `parms` with the parameter values:
"""

# ╔═╡ 15ce9889-a437-46c8-9062-74b8d234a8bd
# parms = missing
parms = [:q=>5, :V=>50, :r=>0.4, :Cin=>3.0, :Ks=>5.2, :kd=>0.10, :Y=>1.2]

# ╔═╡ b8a48461-3882-45f6-980c-38d650ac52c7
md"""
### Creating an ODE problem, solve the problem and plot results
"""

# ╔═╡ a4c57b64-6a7d-4bd4-8bb2-578923e184d2
md"""
Create the ODE problem and store it in `oprob`:
"""

# ╔═╡ b1e18139-5277-4f06-b1f8-b0f5f11c41d8
# oprob = missing
oprob = ODEProblem(sys_ww_treat, u0, tspan, parms)
# or, since we will use the default values here for now:
# oprob = ODEProblem(sys_ww_treat, [], tspan, [])

# ╔═╡ ad6d8fe6-e62f-4c67-8d63-4ee13b928ad0
md"""
Solve the ODE problem. Use `Tsit5()` and `saveat=0.1`. Store the solution in `osol`:
"""

# ╔═╡ e2ffba9e-aaf2-4540-84cf-8b7297ae9285
# osol = missing
osol = solve(oprob, Tsit5(), saveat=0.1)

# ╔═╡ 70871ee8-b0a4-4a9a-af39-5a63459b55f7
md"""
Plot the results. Use `ylim=(0, 4)` and `linewidth=2` as options.
"""

# ╔═╡ 34309734-3751-47e0-a602-d113ffaae510
# begin
#     missing
#     plot!([tspan[1], tspan[2]], [0.28, 0.28],
# 		linestyle=:dash, linewidth=2, linecolor=:green, label="")
# end
begin
    plot(osol, ylim=(0, 4), linewidth=2)
    plot!([tspan[1], tspan[2]], [0.28, 0.28],
		linestyle=:dash, linewidth=2, linecolor=:green, label="C=0.28")
end

# ╔═╡ a0e735ad-09c2-4aa8-bc41-b294a9d56ea8
md"""
Check out the end value of the organic waste.
"""

# ╔═╡ f62898d5-1b8d-4350-8655-78aa3decb2a2
# missing
osol[:C][end]

# ╔═╡ 7eb5df9c-a475-4812-81c3-e43484c82242
md"""
## Part 2

In this part, we will optimize the value of the flow rate $q$ so that the concentration of organic waste in the tank is at most $0.28\;\mathrm{kg\,m^{-3}}$.
"""

# ╔═╡ 5a695734-677f-4bf6-a703-8e22382b7529
md"""
Declare the Turing model function. Sample the flow rate $q$ prior from an uniform distribution in the range $[0, 5]\;\mathrm{kg\,m^{-3}}$. Suppose therein that the desired end value of the organic waste (i.e. $0.28\;\mathrm{kg\,m^{-3}}$) is normally distributed with  mean the end value obtained from the solution and standard deviation $10^{-3}\;\mathrm{kg\,m^{-3}}$.
"""

# ╔═╡ b6bac48a-4a3d-47e4-90ea-788ca20dadff
# @model function wastewater_treatment_inference(C_val)
#     q ~ missing
#     u0 = missing
#     tspan = missing
#     params = missing
#     oprob = missing
#     osol = missing
#     C_val ~ missing
# end
@model function ww_treat_inf(C_val)
	q ~ Uniform(0, 5.0)
    u0 = [:C=>3.0, :X=>0.5]
    tspan = (0.0, 72.0)  # the time interval to solve on
    parms = [:q=>q]      # only parameter q has a new value
    oprob = ODEProblem(sys_ww_treat, u0, tspan, parms)
	# or, since initial conditions don't change:
	# oprob = ODEProblem(sys_ww_treat, [], tspan, parms)
    osol = solve(oprob, Tsit5(), saveat=0.1)
    C_val ~ Normal(osol[:C][end], 1e-3)
end

# ╔═╡ b3a40556-0c00-4f6d-8cd9-c5fca79d8bbf
md"""
Define the desired value for the organic waste with the variable name `C_val`.
"""

# ╔═╡ 2df409ef-bd95-4ac3-a2b8-c5e17c490eba
# missing
C_val = 0.28

# ╔═╡ ee1ffc12-55a1-47ef-ac5b-33148706a09b
md"""
Provide the measurement to the Turing model and optimize the prior. Do this with `MLE` method and Nelder-Mead. Store the optimization results in `results_mle`.
"""

# ╔═╡ afc035be-075b-464b-8ba2-20235082f005
# results_mle = missing
results_mle = optimize(ww_treat_inf(C_val), MLE(), NelderMead())

# ╔═╡ 3ee8121e-3e78-4901-a32d-f04d0c6a0996
md"""
Get the optimized value for $q$ and assign it to `q_opt`.
"""

# ╔═╡ 98a157a1-8c20-474d-acb8-00373ee6d224
# q_opt = missing
q_opt = coef(results_mle)[:q]

# ╔═╡ ceb146c9-a09a-458b-b7d8-3bb7d3de38e0
md"""
Set up parameter values with the optimized parameter value.
"""

# ╔═╡ e275df05-5c77-4c17-ad2e-503574596c31
# parms_opt = missing
parms_opt = [:q=>q_opt]     # only parameter q has a new value

# ╔═╡ cd515dad-44fb-4af2-b933-805ef76be9b3
md"""
Create an ODEProblem and solve it. Use `Tsit5()` and `saveat=0.1`.
"""

# ╔═╡ a4388f06-1223-4815-a557-9b9c3ec232bb
# oprob_opt = missing
oprob_opt = ODEProblem(sys_ww_treat, u0, tspan, parms_opt)

# ╔═╡ ec7bf654-b275-4cfd-a819-d82bdc1be93b
# osol_opt = missing
osol_opt = solve(oprob_opt, Tsit5(), saveat=0.1)

# ╔═╡ 82809c26-4cab-405e-8107-a8a43e81f699
md"""
Plot $C$ and $X$ simulated with the optimized parameter value. Use `ylim=(0, 4)` and `linewidth=2` as options. The dashed line indicates $C = 0.28\;\mathrm{kg\,m^{-3}}$.
"""

# ╔═╡ 81429279-4190-41d1-a72a-20da0ce90528
# begin
#     missing
#     plot!([tspan[1], tspan[2]], [C_val, C_val],
# 		linestyle=:dash, linewidth=2, linecolor=:green, label="C=0.28")
# end
begin
    plot(osol_opt, ylim=(0, 4), linewidth=2)
    plot!([tspan[1], tspan[2]], [C_val, C_val],
		linestyle=:dash, linewidth=2, linecolor=:green, label="C=0.28")
end

# ╔═╡ 6589acfd-1d81-4c10-adea-34ca7fa1ab5d
md"""
Draw your conclusion.

- Conclusion: missing
"""

# ╔═╡ Cell order:
# ╟─c1bc698d-41ee-45e6-b17d-29f0d53557a1
# ╠═f08fa69c-a744-11ef-0e79-3daf5bf297ea
# ╠═2b26c3b2-df08-4b24-a08d-23717248c10d
# ╠═0407d891-a46d-4deb-a21a-23833acbcb87
# ╠═120914e8-afd7-410d-9b57-01a8f24bdf59
# ╠═e48dc930-be03-47b2-b9e3-16e854782aec
# ╠═6458329f-73dd-4cb0-8da4-90678875a1f0
# ╟─dfe77a8c-a8db-46d7-9a4d-3b00413d383b
# ╟─c0c83df7-a9cc-4bde-b6ec-038a423b0d90
# ╟─bbd50cc5-032a-4219-bcee-91145935a7c4
# ╟─93d49247-8e37-4abc-b792-f380c8d59883
# ╟─d48dae1d-a8a2-4b23-b4b8-c840be19ff66
# ╠═08d82f11-068a-4ab3-a664-46af49f1ff11
# ╟─dc2189d6-ca6f-4edd-be41-49ea5a6d556e
# ╠═8dc5c4fd-687f-4dd9-ab80-bf80e40af0df
# ╟─1b3c329c-4cd6-40cc-97bc-e8f829eb04d6
# ╠═799f699c-1c94-406e-bd84-a345b96144c5
# ╠═7a60bdbb-a5ac-40c2-9229-f8daa4cdcd92
# ╟─452734e6-5e51-4488-a0c2-88b86f310f24
# ╠═abd9fcba-568f-4617-a138-a74a06850e9a
# ╟─08ebcb95-8603-4579-879e-810b1494b013
# ╟─d248f64e-ebba-4443-9c49-ff0290aa7810
# ╠═fe02a755-5b00-4d80-a511-fec115b42964
# ╟─67481927-0d03-4da9-af6c-9afa409fc006
# ╠═fadd372a-a665-4b16-9b6d-e32cb7f25d7f
# ╟─734e4d51-95a7-464e-9a23-5ad6c8715d65
# ╠═15ce9889-a437-46c8-9062-74b8d234a8bd
# ╟─b8a48461-3882-45f6-980c-38d650ac52c7
# ╟─a4c57b64-6a7d-4bd4-8bb2-578923e184d2
# ╠═b1e18139-5277-4f06-b1f8-b0f5f11c41d8
# ╟─ad6d8fe6-e62f-4c67-8d63-4ee13b928ad0
# ╠═e2ffba9e-aaf2-4540-84cf-8b7297ae9285
# ╟─70871ee8-b0a4-4a9a-af39-5a63459b55f7
# ╠═34309734-3751-47e0-a602-d113ffaae510
# ╟─a0e735ad-09c2-4aa8-bc41-b294a9d56ea8
# ╠═f62898d5-1b8d-4350-8655-78aa3decb2a2
# ╟─7eb5df9c-a475-4812-81c3-e43484c82242
# ╟─5a695734-677f-4bf6-a703-8e22382b7529
# ╠═b6bac48a-4a3d-47e4-90ea-788ca20dadff
# ╟─b3a40556-0c00-4f6d-8cd9-c5fca79d8bbf
# ╠═2df409ef-bd95-4ac3-a2b8-c5e17c490eba
# ╟─ee1ffc12-55a1-47ef-ac5b-33148706a09b
# ╠═afc035be-075b-464b-8ba2-20235082f005
# ╟─3ee8121e-3e78-4901-a32d-f04d0c6a0996
# ╠═98a157a1-8c20-474d-acb8-00373ee6d224
# ╟─ceb146c9-a09a-458b-b7d8-3bb7d3de38e0
# ╠═e275df05-5c77-4c17-ad2e-503574596c31
# ╟─cd515dad-44fb-4af2-b933-805ef76be9b3
# ╠═a4388f06-1223-4815-a557-9b9c3ec232bb
# ╠═ec7bf654-b275-4cfd-a819-d82bdc1be93b
# ╟─82809c26-4cab-405e-8107-a8a43e81f699
# ╠═81429279-4190-41d1-a72a-20da0ce90528
# ╠═6589acfd-1d81-4c10-adea-34ca7fa1ab5d
