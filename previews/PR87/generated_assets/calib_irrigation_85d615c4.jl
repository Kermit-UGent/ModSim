### A Pluto.jl notebook ###
# v0.20.21

#> [frontmatter]
#> order = "33"
#> title = "6. Calibration irrigation"
#> tags = ["exercises"]
#> layout = "layout.jlhtml"
#> description = "Calibration irrigation"
#> 
#>     [[frontmatter.author]]
#>     name = "Gauthier Vanhaelewyn"

using Markdown
using InteractiveUtils

# ╔═╡ 55675f3d-2fae-4a97-a0a0-ead29a6352e6
# Running this yourself? Point this at your own environment —
# we advise one shared project in the parent folder: Pkg.activate("..")
using Pkg; Pkg.activate("../../pluto-deployment-environment")

# ╔═╡ 2b010e5c-1121-11ef-16fe-a5e3317122e4
using Markdown, InteractiveUtils

# ╔═╡ 4947b0fd-13be-4f6a-b605-ed35b509d7ff
using ModelingToolkit, OrdinaryDiffEq

# ╔═╡ 61d14819-ba44-40fe-95a9-9d7b0bf3dc33
using ModelingToolkit: t_nounits as t, D_nounits as D

# ╔═╡ f6e77c8d-de11-4b9d-93c6-45bdcfbbbf9b
using StatsPlots, StatsBase, Turing

# ╔═╡ 9345dd8f-0a60-4aaf-a27f-ef8bf860f495
using LinearAlgebra, Optim

# ╔═╡ 55d5400d-1777-4918-a030-b94cb9a59f63
md"
# Exercise: Irrigation experiment - Calibration
"

# ╔═╡ 8f1afdec-b78d-4aba-a74f-cd3e4b35fab1
md"""
In one of the previous practica we were introduced to an irrigation experiment carried out on a soil column consisting of two layers of soil, each with specific soil characteristics. However, here the volume of water per unit of time, $r$, irrigated evenly over the soil column, will be kept constant at $5\;mm\,h^{-1}$ in these new experiments.

The water falls on the upper layer and percolates to the lower layer. The relative moisture content in both layers (i.e., relative to their residual moisture contents) is denoted by $S_1$ and $S_2$. 

A model description of the relative moisture content in both soil layers is given by:

$$\begin{align}
\frac{dS_1}{dt} &= r\left(1-\cfrac{S_{1,res}}{S_{max}}\right) - \cfrac{r}{S_{max}}S_1 - \cfrac{k}{S_{max}}S_1 \\
\frac{dS_2}{dt} &= \cfrac{k}{S_{max}}S_1 - v \,S_2^2
\end{align}$$

where $v = 10^{-3}\;h^{-1}\,mm^{-1}$ and $S_{1,res}=10 \;mm$. Previously, we also assumed $k = 3\;mm\,h^{-1}$ and $S_{max} = 150\;mm$.
"""

# ╔═╡ 9f6ad49c-cfbe-46e2-a2dd-c01bb67eeb61
# @variables missing

# ╔═╡ fb2c7db4-e52f-45df-afa0-aeb1db78c849
# @parameters missing

# ╔═╡ c1773f17-825f-414f-9761-5744b34f1b71
# change_S1 = missing

# ╔═╡ 6e2165d4-ef72-4c37-bac0-2dfac417a192
# change_S2 = missing

# ╔═╡ 59654ad7-f42c-443f-bb34-6674db72606f
# @mtkbuild sys_irrigation = missing

# ╔═╡ e5d7520d-fd8c-48c0-bd36-826766212217
md"""
In order to have better estimates the parameters $k$ and $S_{max}$, two experiments were conducted, each with a different initial condition:

1. Starting from zero relative moisture content in both soil layers.
2. Starting from a relative moisture content of $140\;mm$ in the top layer, and $135\;mm$ in the bottom layer.

The measurement data consist of measurements of the relative moisture contents $S_1$ and $S_2$ measured at intervals of $10\;h$ within a timespan of $150\;h$.
"""

# ╔═╡ 73c9b5fb-4f56-4bde-beb4-387651409c1b
md"""
The measurement data for the 1st experiment are:
"""

# ╔═╡ 9f94c63e-628f-4ff3-ad29-0f90d32dfcb1
S1_meas1 = [0.2, 35.94, 52.49, 66.86, 60.66, 67.81, 73.22, 71.31, 72.94, 64.08, 70.11, 68.53, 70.54, 63.63, 67.39, 62.84]

# ╔═╡ 3d406f41-62fa-4a31-8e6b-a621a06a118c
S2_meas1 = [0.63, 6.2, 17.67, 22.96, 35.41, 44.08, 43.5, 53.34, 47.57, 47.77, 43.96, 52.22, 46.67, 46.74, 46.46, 39.92]

# ╔═╡ 68b6158a-a918-4809-bf81-b554bc70c6d0
md"""
The measurement data for the 2nd experiment are:
"""

# ╔═╡ 620729c8-62c2-4f0c-8684-713033a208bd
S1_meas2 = [137.96, 106.15, 90.15, 84.64, 76.15, 75.73, 73.32, 68.48, 70.06, 69.36, 70.91, 72.13, 76.25, 74.34, 74.93, 71.58]

# ╔═╡ 5b83a98c-9c2b-459c-9d7f-7cc75d3bf70e
S2_meas2 = [124.08, 80.14, 60.15, 50.12, 49.66, 47.78, 46.56, 48.41, 42.7, 43.72, 49.03, 51.91, 48.24, 46.14, 51.22, 43.78]

# ╔═╡ 9c568f28-8985-4a4a-a7aa-0010bbe37dc8
md"""
For both experiments:
"""

# ╔═╡ 0dc6fa2c-1eb2-4877-9155-dc7ea6cf6f18
t_meas = 0:10:150

# ╔═╡ a4160054-3dfe-4595-81cb-94db4dd2fe20
md"""
We can make a scatter plot of the measured data for both $S_1$ and $S_2$ for the 1st and 2nd experiments in the following way:
"""

# ╔═╡ cef4b9a8-b5bf-4a2d-8a8f-5d8f85534859
begin
  scatter(t_meas, S1_meas1, label="S1 meas", color=:blue, title="Experiment 1")
  scatter!(t_meas, S2_meas1, label="S2 meas", color=:red, ylims=(0, 150))
end

# ╔═╡ fc2cabd7-e778-4211-bf87-b5c11ca054c9
begin
  scatter(t_meas, S1_meas2, label="S1 meas", color=:blue, title="Experiment 2")
  scatter!(t_meas, S2_meas2, label="S2 meas", color=:red, ylims=(0, 150))
end

# ╔═╡ c0b2db7b-0632-4008-9cff-d5fbf3e59807
md"""
Calibrate the parameter values for $k$ and $S_{max}$ using the aforementioned measurement data for $S_1$ and $S_2$ in a timespan of $[0, 150]\,h$. Take the values from above as initial values for $k$ and $S_{max}$.
"""

# ╔═╡ a65a0997-9945-436e-925b-8fc02055f61a
# tspan = missing

# ╔═╡ 923d04ce-b4d2-44b0-afff-7062c4628ad0
md"""
Declare the Turing model. Make sure you take both experiments into account for optimizing $k$ and $S_{max}$. Based on literature, you can assume that the value of $S_{max}$ lies somewhere between 100 and 200 mm.
"""

# ╔═╡ 481eb8b9-5de2-4f68-b06a-ec18e054c9f5
# @model function irrigation_inference(t_meas)
# 	σ_S1 ~ missing
# 	σ_S2 ~ missing
# 	k ~ missing
# 	Smax ~ missing
# 	parms = missing
#   # For experiment 1:
# 	u01 = missing
# 	oprob1 = missing
# 	osol1 = missing
# 	S1_ex1 ~ missing
# 	S2_ex1 ~ missing
#   # For experiment 2:
# 	u02 = missing
# 	oprob2 = missing
# 	osol2 = missing
# 	S1_ex2 ~ missing
# 	S2_ex2 ~ missing
# end

# ╔═╡ df933ae8-1f51-4467-93a7-33f153e5e4f8
md"""
Instantiate the model and condition it with the measurements of $S_1$ and $S_2$ from both experiments:
"""

# ╔═╡ 0e2aa675-9e09-4e06-b5f8-118707ee652a
# irrigation_inf = missing   

# ╔═╡ f7f47956-7c3b-44cc-bff7-fb7d32af874a
md"""
Optimize the priors ($\sigma_{S1}$, $\sigma_{S2}$, $k$ and $S_{max}$). Do this with `MLE` method and Nelder-Mead. Store the optimization results in `results_mle`.
"""

# ╔═╡ 8c254d5a-225b-4772-9fdd-e9f700495fbd
# results_mle = missing      

# ╔═╡ f15a1df5-047a-4f46-9419-8492ac1248e0
md"""
Visualize a summary of the optimized parameters. Beware that this may take a lot of time...
"""

# ╔═╡ 00d944e4-2c88-4a5d-b809-69f435df4684
# missing           

# ╔═╡ 89eb31ef-b24f-44c8-bbe5-19101d859937
md"""
Get the optimized values and assign them to `k_opt` and `Smax_opt`.
"""

# ╔═╡ 92daa779-3373-40c0-b308-23e75e6674b6
# k_opt = missing               

# ╔═╡ 35ab6ee5-fcd7-4dcc-9909-cc918fb1fe80
# Smax_opt = missing            

# ╔═╡ 4026773f-ac5b-433e-bd9d-2122242861fd
md"""
Make plots of $S_1$ and $S_2$ for both experiments simulated with the optimized parameter values.
"""

# ╔═╡ 8aa60652-eb9f-4dd3-ab06-0ce3dd261fe6
md"""
Set up parameter values with optimized parameter values:
"""

# ╔═╡ 97d53e48-590a-485b-bcf3-edc6a6124faf
# params_opt = missing         

# ╔═╡ dfd2ac98-5cdc-4627-b6cf-71b33c0ff0d4
md"""
Plot the simulation results $S_1$ and $S_2$ for the 1st experiment together with the corresponding measured data. Therefore initialize a vector `u01` with initial conditions for the 1st experiment.
"""

# ╔═╡ 95ace332-52c0-46c3-ae28-d038320ed2c8
# u01 = missing                

# ╔═╡ 6ae63a13-d5ae-4dfb-b88d-be295b11a472
# oprob1_opt = missing         

# ╔═╡ bc6505ca-a61d-467f-afe6-47792a510ad5
# osol1_opt = missing          

# ╔═╡ 67e423ea-e941-45bf-af4f-3fdecb648fbc
# begin
#   missing
#   missing
#   missing
# end

# ╔═╡ 8c7e0c75-01d2-4e21-a6cf-7ad70e0c6aae
md"""
Plot the simulation results $S_1$ and $S_2$ for the 1st experiment together with the corresponding measured data. Therefore initialize a vector `u02` with initial conditions for the 2nd experiment.
"""

# ╔═╡ a7040b8e-c240-415b-8a9a-4a1a137398d4
# u02 = missing                  

# ╔═╡ fe8f4961-68bd-42dc-a3f5-6692e918e241
# oprob2_opt = missing           

# ╔═╡ 7f280230-7846-4529-a2ff-a81a2b9480bf
# osol2_opt = missing            

# ╔═╡ ad9818a9-ccbe-4645-8b91-0c3fa773632a
# begin
#   missing
#   missing
#   missing
# end

# ╔═╡ 3c243670-2ba7-4396-8c6d-084726636741
md"""
!!! question
	Do your simulations fit well the measurements?
"""

# ╔═╡ 1aeb0d86-b276-4bd6-9811-faf4a297ae6f
md"- Answer: missing"

# ╔═╡ Cell order:
# ╠═2b010e5c-1121-11ef-16fe-a5e3317122e4
# ╠═55675f3d-2fae-4a97-a0a0-ead29a6352e6
# ╠═4947b0fd-13be-4f6a-b605-ed35b509d7ff
# ╠═61d14819-ba44-40fe-95a9-9d7b0bf3dc33
# ╠═f6e77c8d-de11-4b9d-93c6-45bdcfbbbf9b
# ╠═9345dd8f-0a60-4aaf-a27f-ef8bf860f495
# ╟─55d5400d-1777-4918-a030-b94cb9a59f63
# ╟─8f1afdec-b78d-4aba-a74f-cd3e4b35fab1
# ╠═9f6ad49c-cfbe-46e2-a2dd-c01bb67eeb61
# ╠═fb2c7db4-e52f-45df-afa0-aeb1db78c849
# ╠═c1773f17-825f-414f-9761-5744b34f1b71
# ╠═6e2165d4-ef72-4c37-bac0-2dfac417a192
# ╠═59654ad7-f42c-443f-bb34-6674db72606f
# ╟─e5d7520d-fd8c-48c0-bd36-826766212217
# ╟─73c9b5fb-4f56-4bde-beb4-387651409c1b
# ╠═9f94c63e-628f-4ff3-ad29-0f90d32dfcb1
# ╠═3d406f41-62fa-4a31-8e6b-a621a06a118c
# ╟─68b6158a-a918-4809-bf81-b554bc70c6d0
# ╠═620729c8-62c2-4f0c-8684-713033a208bd
# ╠═5b83a98c-9c2b-459c-9d7f-7cc75d3bf70e
# ╟─9c568f28-8985-4a4a-a7aa-0010bbe37dc8
# ╠═0dc6fa2c-1eb2-4877-9155-dc7ea6cf6f18
# ╟─a4160054-3dfe-4595-81cb-94db4dd2fe20
# ╠═cef4b9a8-b5bf-4a2d-8a8f-5d8f85534859
# ╠═fc2cabd7-e778-4211-bf87-b5c11ca054c9
# ╟─c0b2db7b-0632-4008-9cff-d5fbf3e59807
# ╠═a65a0997-9945-436e-925b-8fc02055f61a
# ╟─923d04ce-b4d2-44b0-afff-7062c4628ad0
# ╠═481eb8b9-5de2-4f68-b06a-ec18e054c9f5
# ╟─df933ae8-1f51-4467-93a7-33f153e5e4f8
# ╠═0e2aa675-9e09-4e06-b5f8-118707ee652a
# ╟─f7f47956-7c3b-44cc-bff7-fb7d32af874a
# ╠═8c254d5a-225b-4772-9fdd-e9f700495fbd
# ╟─f15a1df5-047a-4f46-9419-8492ac1248e0
# ╠═00d944e4-2c88-4a5d-b809-69f435df4684
# ╟─89eb31ef-b24f-44c8-bbe5-19101d859937
# ╠═92daa779-3373-40c0-b308-23e75e6674b6
# ╠═35ab6ee5-fcd7-4dcc-9909-cc918fb1fe80
# ╟─4026773f-ac5b-433e-bd9d-2122242861fd
# ╟─8aa60652-eb9f-4dd3-ab06-0ce3dd261fe6
# ╠═97d53e48-590a-485b-bcf3-edc6a6124faf
# ╟─dfd2ac98-5cdc-4627-b6cf-71b33c0ff0d4
# ╠═95ace332-52c0-46c3-ae28-d038320ed2c8
# ╠═6ae63a13-d5ae-4dfb-b88d-be295b11a472
# ╠═bc6505ca-a61d-467f-afe6-47792a510ad5
# ╠═67e423ea-e941-45bf-af4f-3fdecb648fbc
# ╟─8c7e0c75-01d2-4e21-a6cf-7ad70e0c6aae
# ╠═a7040b8e-c240-415b-8a9a-4a1a137398d4
# ╠═fe8f4961-68bd-42dc-a3f5-6692e918e241
# ╠═7f280230-7846-4529-a2ff-a81a2b9480bf
# ╠═ad9818a9-ccbe-4645-8b91-0c3fa773632a
# ╟─3c243670-2ba7-4396-8c6d-084726636741
# ╠═1aeb0d86-b276-4bd6-9811-faf4a297ae6f
