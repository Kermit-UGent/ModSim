### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 55675f3d-2fae-4a97-a0a0-ead29a6352e6
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 2b010e5c-1121-11ef-16fe-a5e3317122e4
using Markdown

# ╔═╡ 6750d246-8e7a-4cfb-810e-d1100aa4fdef
using InteractiveUtils

# ╔═╡ 4947b0fd-13be-4f6a-b605-ed35b509d7ff
using Catalyst, DifferentialEquations, Plots

# ╔═╡ ea02aff2-7fb8-4b8f-8d0f-0bb3c6150708
using Optim

# ╔═╡ 55d5400d-1777-4918-a030-b94cb9a59f63
md"
### Exercise: Irrigation experiment - Calibration
"

# ╔═╡ 8f1afdec-b78d-4aba-a74f-cd3e4b35fab1
md"
In one of the previous practica we were introduced to an irrigation experiment carried out on a soil column consisting of two layers of soil, each with specific soil characteristics. However, here the volume of water per unit of time, $R$, irrigated evenly over the soil column, will be kept constant at $5\;mm\,h^{-1}$ in these new experiments.

The water falls on the upper layer and percolates to the lower layer. The relative moisture content in both layers (i.e., relative to their residual moisture contents) is denoted by $S_1$ and $S_2$. 

A model description of the relative moisture content in both soil layers is given by:

$$\begin{align}
\frac{dS_1}{dt} &= R\left(1-\cfrac{S_{1,res}}{S_{max}}\right) - \cfrac{R}{S_{max}}S_1 - \cfrac{k}{S_{max}}S_1 \\
\frac{dS_2}{dt} &= \cfrac{k}{S_{max}}S_1 - v \,S_2^2
\end{align}$$

where $v = 10^{-3}\;h^{-1}\,mm^{-1}$ and $S_{1,res}=10 \;mm$. Previously, we also assumed $k = 3\;mm\,h^{-1}$ and $S_{max} = 150\;mm$.
"

# ╔═╡ ad42d3a7-6d83-4362-aa3f-31628a1db9b2
md"
The *reaction network object* for this model could be set up as:
"

# ╔═╡ dc26abff-f8ab-4881-9acf-7b325b386a16
irrigation_mod = @reaction_network begin
    k/Smax, S₁ --> S₂
    v, 2S₂ --> ∅
    R * (1 - S₁res / Smax), ∅ --> S₁
    R/Smax, S₁ --> ∅
end

# ╔═╡ e5d7520d-fd8c-48c0-bd36-826766212217
md"
In order to have better estimates the parameters $k$ and $S_{max}$, two experiments were conducted, each with a different initial condition:

1. Starting from zero relative moisture content in both soil layers.
2. Starting from a relative moisture content of $140\;mm$ in the top layer, and $135\;mm$ in the bottom layer.

The measurement data consist of measurements of the relative moisture contents $S_1$ and $S_2$ measured at intervals of $10\;h$ within a timespan of $150\;h$.
"

# ╔═╡ 73c9b5fb-4f56-4bde-beb4-387651409c1b
md"
The measurement data for the 1st experiment are:
"

# ╔═╡ 9f94c63e-628f-4ff3-ad29-0f90d32dfcb1
S1_meas1 = [0.2, 35.94, 52.49, 66.86, 60.66, 67.81, 73.22, 71.31, 72.94, 64.08, 70.11, 68.53, 70.54, 63.63, 67.39, 62.84]

# ╔═╡ 3d406f41-62fa-4a31-8e6b-a621a06a118c
S2_meas1 = [0.63, 6.2, 17.67, 22.96, 35.41, 44.08, 43.5, 53.34, 47.57, 47.77, 43.96, 52.22, 46.67, 46.74, 46.46, 39.92]

# ╔═╡ 68b6158a-a918-4809-bf81-b554bc70c6d0
md"
The measurement data for the 2nd experiment are:
"

# ╔═╡ 620729c8-62c2-4f0c-8684-713033a208bd
S1_meas2 = [137.96, 106.15, 90.15, 84.64, 76.15, 75.73, 73.32, 68.48, 70.06, 69.36, 70.91, 72.13, 76.25, 74.34, 74.93, 71.58]

# ╔═╡ 5b83a98c-9c2b-459c-9d7f-7cc75d3bf70e
S2_meas2 = [124.08, 80.14, 60.15, 50.12, 49.66, 47.78, 46.56, 48.41, 42.7, 43.72, 49.03, 51.91, 48.24, 46.14, 51.22, 43.78]

# ╔═╡ 9c568f28-8985-4a4a-a7aa-0010bbe37dc8
md"
For both experiments:
"

# ╔═╡ 0dc6fa2c-1eb2-4877-9155-dc7ea6cf6f18
t_meas = 0:10:150

# ╔═╡ 9115e49a-edaa-4b8c-b88e-05095f5c142f
md"
In addition, assume the follwoing measurement errors:

$$\sigma_{S_1} = 4\;mm \quad\quad\quad \sigma_{S_2} = 6\;mm$$

We will call them `S1_sigma` and `S2_sigma`:
"

# ╔═╡ 4502a26e-b8b5-4f50-b49c-6b845536a0ea
S1_sigma = 4.0

# ╔═╡ 1eee4185-6346-4748-92f8-a1aab11faeed
S2_sigma = 6.0

# ╔═╡ a4160054-3dfe-4595-81cb-94db4dd2fe20
md"
We can make a scatter plot of the measured data for both $S_1$ and $S_2$ for the 1st and 2nd experiments in the following way:
"

# ╔═╡ cef4b9a8-b5bf-4a2d-8a8f-5d8f85534859
begin
  scatter(t_meas, S1_meas1, label="S1 meas", color=:blue, title="Experment 1")
  scatter!(t_meas, S2_meas1, label="S2 meas", color=:red, ylims=(0, 150))
end

# ╔═╡ fc2cabd7-e778-4211-bf87-b5c11ca054c9
begin
  scatter(t_meas, S1_meas2, label="S1 meas", color=:blue, title="Experment 2")
  scatter!(t_meas, S2_meas2, label="S2 meas", color=:red, ylims=(0, 150))
end

# ╔═╡ c0b2db7b-0632-4008-9cff-d5fbf3e59807
md"
Calibrate the parameter values for $k$ and $S_{max}$ using the aforementioned measurement data for $S_1$ and $S_2$ in a timespan of $[0, 150]\,h$. Take the values from above as initial values for $k$ and $S_{max}$.
"

# ╔═╡ f37538b7-2627-4a0a-b682-1eed0134fd97
md"
Set the timespan:
"

# ╔═╡ f4d49b1d-9105-4050-a9d4-196fa00a0591
# tspan = missing          # Uncomment and complete the instruction
tspan = (0.0, 150.0)

# ╔═╡ 923d04ce-b4d2-44b0-afff-7062c4628ad0
md"
Set up the objective function and name it `Jtheta_irrigation`. Your objective function should take both experiments into account.
"

# ╔═╡ 017ae458-1047-462c-860f-8014e56e54d9
# Uncomment and complete the instruction
# function Jtheta_irrigation(thetas)
#   params = missing
#   # Experiment 1
#   u₀ = missing
#   oprob = missing
#   osol = missing
#   S1_sol = missing
#   S2_sol = missing
#   J1 = missing
#   # Experiment 2
#   u₀ = missing
#   oprob = missing
#   osol = missing
#   S1_sol = missing
#   S2_sol = missing
#   J2 = missing
#   return J1+J2
# end
function Jtheta_irrigation(thetas)
  params = [:k => thetas[1], :Smax => thetas[2], :v => 1e-3, :R => 5, :S₁res => 10]
  # Experiment 1
  u₀ = [0.0, 0.0]
  oprob = ODEProblem(irrigation_mod, u₀, tspan, params)
  osol = solve(oprob, Tsit5(), saveat=t_meas)
  S1_sol = osol[:S₁]
  S2_sol = osol[:S₂]
  J1 = (1/S1_sigma^2)*sum(abs2, S1_sol - S1_meas1) +
       (1/S2_sigma^2)*sum(abs2, S2_sol - S2_meas1)
  # Experiment 2
  u₀ = [140.0, 135.0]
  oprob = ODEProblem(irrigation_mod, u₀, tspan, params)
  osol = solve(oprob, Tsit5(), saveat=t_meas)
  S1_sol = osol[:S₁]
  S2_sol = osol[:S₂]
  J2 = (1/S1_sigma^2)*sum(abs2, S1_sol - S1_meas2) +
       (1/S2_sigma^2)*sum(abs2, S2_sol - S2_meas2)
  return J1+J2
end

# ╔═╡ 7e0c672e-c78c-42a6-915e-de2984188cfc
md"
Create a vector `thetas_init` with the initial values for $k$ and $S_{max}$:
"

# ╔═╡ 88e91c8e-8a44-4e99-a791-9c6aac78a677
# thetas_init = missing       # Uncomment and complete the instruction
thetas_init = [3.0, 150.0]

# ╔═╡ df933ae8-1f51-4467-93a7-33f153e5e4f8
md"
Minimize the value of the objective function by optimizing the parameter values. Store the optimization results in `results`.
"

# ╔═╡ 6a2f746e-23c8-4315-b780-ac5caa38e098
# results = missing          # Uncomment and complete the instruction
results = optimize(Jtheta_irrigation, thetas_init, NelderMead())

# ╔═╡ f7f47956-7c3b-44cc-bff7-fb7d32af874a
md"
Store the optimized values in the variable `thetas_opt`:
"

# ╔═╡ 8f24c1a2-abeb-43c4-a926-f0c8e10eb8ba
# thetas_opt = missing       # Uncomment and complete the instruction
thetas_opt = Optim.minimizer(results)

# ╔═╡ 0702dc6c-7d56-4fa4-8e46-24650dd31bbe
md"
Determine the minimum value of the objective function:
"

# ╔═╡ c6e981c9-f3bb-4538-8e4f-36f5af879fee
# missing                    # Uncomment and complete the instruction
Optim.minimum(results)

# ╔═╡ f15a1df5-047a-4f46-9419-8492ac1248e0
md"
Set up parameter values with optimized parameter values:
"

# ╔═╡ 650d2773-255f-455e-823d-8534481cfebf
# params = missing           # Uncomment and complete the instruction
params = [:k => thetas_opt[1], :Smax => thetas_opt[2], :v => 1e-3, :R => 5, :S₁res => 10]

# ╔═╡ 4026773f-ac5b-433e-bd9d-2122242861fd
md"
Simulate $S_1$ and $S_2$ for the 1st experiment using the optimized parameter values. Therefore, initialize a vector `u₀₁` with the initial conditions for the 1st experiment, create an ODEProblem and solve it.
"

# ╔═╡ 8c4e4106-f5a6-4532-b6b5-0832cf9c091b
# u₀₁ = missing               # Uncomment and complete the instruction
u₀₁ = [0.0, 0.0]

# ╔═╡ 7ed06684-8802-4f39-8091-36c615a2d35c
# oprob1 = missing            # Uncomment and complete the instruction
oprob1 = ODEProblem(irrigation_mod, u₀₁, tspan, params)

# ╔═╡ c399cca9-ee91-4f4b-bd25-b8185174eed9
# osol1 = missing             # Uncomment and complete the instruction
osol1 = solve(oprob1, Tsit5(), saveat=1.0)

# ╔═╡ dfd2ac98-5cdc-4627-b6cf-71b33c0ff0d4
md"
Plot the simulation results $S_1$ and $S_2$ for the 1st experiment together with the corresponding measured data.
"

# ╔═╡ aab21043-a0fa-410b-82c3-5bb703065448
# Uncomment and complete the instruction
# begin
#   missing
#   missing
#   missing
# end
begin
  plot(osol1, labels=["S1 sim" "S2 sim"], xlabel="t")
  scatter!(t_meas, S1_meas1, label="S1 meas", color=:blue, title="Experment 1")
  scatter!(t_meas, S2_meas1, label="S2 meas", color=:red, ylims=(0, 150))
end

# ╔═╡ 8c7e0c75-01d2-4e21-a6cf-7ad70e0c6aae
md"
Simulate $S_1$ and $S_2$ for the 2nd experiment using the optimized parameter values. Therefore, initialize a vector `u₀₂` with the initial conditions for the 2nd experiment, create an ODEProblem and solve it.
"

# ╔═╡ 497741a2-acc9-4d48-bbf8-eb0d711831d2
# u₀₂ = missing          # Uncomment and complete the instruction
u₀₂ = [140.0, 135.0]

# ╔═╡ ea09016c-6234-4e32-b91c-b4375fd46eac
# oprob2 = missing       # Uncomment and complete the instruction
oprob2 = ODEProblem(irrigation_mod, u₀₂, tspan, params)

# ╔═╡ 4584c051-0bb9-4813-b7fd-08ca83904144
# osol2 = missing        # Uncomment and complete the instruction
osol2 = solve(oprob2, Tsit5(), saveat=1.0)

# ╔═╡ 42e5d3f9-314b-4f7a-a72b-1ee2acfec890
md"
Plot the simulation results $S_1$ and $S_2$ for the 2nd experiment together with the corresponding measured data.
"

# ╔═╡ 190fa7b1-7376-4911-8f07-ad323193a867
# Uncomment and complete the instruction
# begin
#   missing
#   missing
#   missing
# end
begin
  plot(osol2, labels=["S1 sim" "S2 sim"], xlabel="t")
  scatter!(t_meas, S1_meas2, label="S1 meas", color=:blue, title="Experment 2")
  scatter!(t_meas, S2_meas2, label="S2 meas", color=:red, ylims=(0, 150))
end

# ╔═╡ Cell order:
# ╠═2b010e5c-1121-11ef-16fe-a5e3317122e4
# ╠═6750d246-8e7a-4cfb-810e-d1100aa4fdef
# ╠═55675f3d-2fae-4a97-a0a0-ead29a6352e6
# ╠═4947b0fd-13be-4f6a-b605-ed35b509d7ff
# ╠═ea02aff2-7fb8-4b8f-8d0f-0bb3c6150708
# ╠═55d5400d-1777-4918-a030-b94cb9a59f63
# ╠═8f1afdec-b78d-4aba-a74f-cd3e4b35fab1
# ╠═ad42d3a7-6d83-4362-aa3f-31628a1db9b2
# ╠═dc26abff-f8ab-4881-9acf-7b325b386a16
# ╠═e5d7520d-fd8c-48c0-bd36-826766212217
# ╠═73c9b5fb-4f56-4bde-beb4-387651409c1b
# ╠═9f94c63e-628f-4ff3-ad29-0f90d32dfcb1
# ╠═3d406f41-62fa-4a31-8e6b-a621a06a118c
# ╠═68b6158a-a918-4809-bf81-b554bc70c6d0
# ╠═620729c8-62c2-4f0c-8684-713033a208bd
# ╠═5b83a98c-9c2b-459c-9d7f-7cc75d3bf70e
# ╠═9c568f28-8985-4a4a-a7aa-0010bbe37dc8
# ╠═0dc6fa2c-1eb2-4877-9155-dc7ea6cf6f18
# ╠═9115e49a-edaa-4b8c-b88e-05095f5c142f
# ╠═4502a26e-b8b5-4f50-b49c-6b845536a0ea
# ╠═1eee4185-6346-4748-92f8-a1aab11faeed
# ╠═a4160054-3dfe-4595-81cb-94db4dd2fe20
# ╠═cef4b9a8-b5bf-4a2d-8a8f-5d8f85534859
# ╠═fc2cabd7-e778-4211-bf87-b5c11ca054c9
# ╠═c0b2db7b-0632-4008-9cff-d5fbf3e59807
# ╠═f37538b7-2627-4a0a-b682-1eed0134fd97
# ╠═f4d49b1d-9105-4050-a9d4-196fa00a0591
# ╠═923d04ce-b4d2-44b0-afff-7062c4628ad0
# ╠═017ae458-1047-462c-860f-8014e56e54d9
# ╠═7e0c672e-c78c-42a6-915e-de2984188cfc
# ╠═88e91c8e-8a44-4e99-a791-9c6aac78a677
# ╠═df933ae8-1f51-4467-93a7-33f153e5e4f8
# ╠═6a2f746e-23c8-4315-b780-ac5caa38e098
# ╠═f7f47956-7c3b-44cc-bff7-fb7d32af874a
# ╠═8f24c1a2-abeb-43c4-a926-f0c8e10eb8ba
# ╠═0702dc6c-7d56-4fa4-8e46-24650dd31bbe
# ╠═c6e981c9-f3bb-4538-8e4f-36f5af879fee
# ╠═f15a1df5-047a-4f46-9419-8492ac1248e0
# ╠═650d2773-255f-455e-823d-8534481cfebf
# ╠═4026773f-ac5b-433e-bd9d-2122242861fd
# ╠═8c4e4106-f5a6-4532-b6b5-0832cf9c091b
# ╠═7ed06684-8802-4f39-8091-36c615a2d35c
# ╠═c399cca9-ee91-4f4b-bd25-b8185174eed9
# ╠═dfd2ac98-5cdc-4627-b6cf-71b33c0ff0d4
# ╠═aab21043-a0fa-410b-82c3-5bb703065448
# ╠═8c7e0c75-01d2-4e21-a6cf-7ad70e0c6aae
# ╠═497741a2-acc9-4d48-bbf8-eb0d711831d2
# ╠═ea09016c-6234-4e32-b91c-b4375fd46eac
# ╠═4584c051-0bb9-4813-b7fd-08ca83904144
# ╠═42e5d3f9-314b-4f7a-a72b-1ee2acfec890
# ╠═190fa7b1-7376-4911-8f07-ad323193a867
