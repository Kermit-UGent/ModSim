### A Pluto.jl notebook ###
# v0.20.4

#> [frontmatter]
#> order = "7"
#> title = "EXTRA. ODE fermentor monod"
#> date = "2025-02-07"
#> tags = ["exercises"]
#> description = "ODE model of a fermentor monod"
#> layout = "layout.jlhtml"
#> 
#>     [[frontmatter.author]]
#>     name = "Gauthier Vanhaelewyn"

using Markdown
using InteractiveUtils

# ╔═╡ e99680dc-73af-40aa-bf57-a06d3a7372be
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate("../../pluto-deployment-environment")
end

# ╔═╡ 2e58f4ae-f711-11ee-2598-7f3a6f2e2013
using Markdown

# ╔═╡ 539c1823-16d0-4355-97b6-fe9f0b106864
using InteractiveUtils

# ╔═╡ e66518ee-b6f6-4cca-a224-30e01cffddbe
using Catalyst

# ╔═╡ bd648109-f042-42de-9e0e-017b502fab95
using DifferentialEquations, Plots

# ╔═╡ 7856d878-8586-4cfd-9cf6-d61234450e41
md"
### Exercise: Fermenter - Monod kinetics
"

# ╔═╡ 8500c35e-6bc3-4900-81bf-7705ddd61532
md"""
In a fermenter reactor biomass grows on substrate. The reactor is fed with a inlet flow rate $Q_{in}$ [$L/h$], which consist of a (manipulable) input concentration of substrate $S_{in}$ [$g/L$]. Inside the reactor, biomass, with a concentration of $X$ [$g/L$], is produced through **Monod** kinetics:

$$\begin{eqnarray*}
%S  \xrightarrow[\quad\quad]{\beta} Y \, X
S \xrightarrow[\quad\quad]{r} Y \, X \quad\quad\quad\quad r = \mu \, X
\end{eqnarray*}$$

where

$$\mu = \mu_{max} \, \cfrac{S}{S + K_s}$$

is called the specific growth rate [$h^{-1}$]. Therein, $\mu_{max}$ is the maximum speficic growth rate, and $K_s$ [$g/L$] is the so-called *half-velocity constant* (i.e. the value of $S$ when $\mu/\mu_{max} = 0.5$). Futhermore, $Y$ [$gX/gS$] is the yield coefficient which is defined here by the amount of produced biomass by consumption of one unit of substrate. The reactor is drained with an outlet flow $Q$ [$L/h$], which consist of the current concentrations of substrate $S$ [$g/L$] and biomass $X$ [$g/L$] inside the reactor. The volume $V$ [$L$] of the reactor content is kept constant by setting $Q_{in} = Q$.
"""

# ╔═╡ f1350528-07a5-4860-ad2d-627588186abc
md"""
Create a *reaction network object* model for the aforementioned problem in order to simulate the evolution of substrate $S$ and biomass $X$ with time. Name it `fermenter_monod`.

Tip: The specific growth rate $\mu = \mu_{max} \, \cfrac{S}{S + K_s}$ can be implemented with `mm(S, μmax, Ks)`. The function `mm` stands for the Michaelis-Menten kinetics, whcih is equivalent to Monod kinetics.
"""

# ╔═╡ 331a34f4-89d4-4193-896c-c14ab0bf04e7
# fermenter_monod = @reaction_network begin
#     ...        # Y*X is created from one S at a rate mm(S, μmax, Ks)*X
#     ...        # S is created at a rate Q/V*Sin
#     ...        # S and X are degraded at a rate Q/V*S
# end

# ╔═╡ 55746566-2d46-4475-851a-02b7fad87a1a
md"""
Convert the system to a symbolic differential equation model and verify, by analyzing the differential equation, that your model is correctly implemented.

Keep in mind that `mm(S, μmax, Ks)` stands for $\mu_{max} \, \cfrac{S}{S + K_s}$.
"""

# ╔═╡ ec9cb3bd-f5ed-4ab0-9b3d-b875692227ac
# osys = missing

# ╔═╡ 67117a27-dcea-4b43-b962-9ad9fd07f4f4
md"""
The parameter values are $\mu_{max} = 0.40$, $K_s = 0.015$, $Y = 0.67$, $Q = 2.0$, $V = 40.0$ and $S_{in} = 0.02\;g/L$. Suppose that at $t=0\;h$ no substrate $S$ is present in the reactor but that there is initially some biomass with a concetration of $0.0005\;g/L$. Simulate the evolution of $S$ and $X$ during $200$ hours.
"""

# ╔═╡ d13e6e38-037e-4812-85e9-2c18bed360f6
md"""
Initialize a vector `u0` with the initial conditions:
"""

# ╔═╡ 4b556cf0-8fad-434d-be56-dc1848d898ae
# u0 = missing            # Uncomment and complete the instruction

# ╔═╡ ea55d648-7575-43c3-a385-5f4979996ef2
md"""
Set the timespan for the simulation:
"""

# ╔═╡ 1365c12e-e662-4858-983b-02ba94cd9f0d
# tspan = missing         # Uncomment and complete the instruction

# ╔═╡ 3941bd60-a83c-4f72-84b3-28e28cb845d0
md"""
Initialize a vector `param` with the parameter values:
"""

# ╔═╡ d6c1316a-cf96-43d1-854a-f25925cf4a55
# params = missing         # Uncomment and complete the instruction

# ╔═╡ 4926b941-c3b6-4804-b4a4-11e13e5186f2
md"""
Create the ODE problem and store it in `oprob`:
"""

# ╔═╡ ab2a9842-6a9c-46bd-812b-db01629d6a1c
# oprob = missing           # Uncomment and complete the instruction

# ╔═╡ b6a526bd-6ee5-442b-9fb8-3fbe1e280dd4
md"""
#### Part 1

Solve the ODE problem. Use `Tsit5()` and `saveat=0.5`. Store the solution in `osol1`:
"""

# ╔═╡ 1f62e66d-571f-41ca-9f02-f36a8ca10ab9
# osol1 = missing           # Uncomment and complete the instruction

# ╔═╡ 37ced6e3-b435-4546-a720-a1ec1af23a65
md"""
Plot the results:
"""

# ╔═╡ d609bed4-94cf-4167-80fe-924501a5835c
# missing             # Uncomment and complete the instruction

# ╔═╡ 3c638fc4-8ebc-4e26-984e-b4513035287e
md"""
Inspect the final values in both the $S$ and $X$ vector.\
Tip: use something like: `(osol1[...][...], osol1[...][...])`
"""

# ╔═╡ 25ef069e-5a08-4e85-acea-a5b64e0890f6
# (osol1[...][...], osol1[...][...])   # Uncomment and complete the instruction

# ╔═╡ d8234337-7516-4817-8ce5-194af694e3e3
md"""
We will now show you how to determine the steady state values for $S$ and $X$ under the current conditions (cf. current initial values and current parameter values).
"""

# ╔═╡ b63fa707-6023-44fc-aac9-25c76e075f90
md"""
First, we initialize a vector `u_guess1` with the final values for $S$ and $X$:
"""

# ╔═╡ 8d96d79f-6ddb-4bba-a6ea-821588e13107
# u_guess1 = missing

# ╔═╡ e8969045-27ac-460e-86a2-7494903534e8
md"""
Then we make a so-called SteadyStateProblem based on the ODEProblem but now with `u_guess1` as initial conditions! Finally we use `solve` to solve the steady state problem. The outputs are the steady state values for $S$ and $X$ which we have denoted as `Seq1` and `Xeq1`.
"""

# ╔═╡ 1777503e-b793-4be2-b80b-b4edcd7041b5
# Seq1, Xeq1 = missing

# ╔═╡ 1d5a2118-96e9-49b2-9931-5d4b201cb8f5
md"""
Next, we can just inspect these values:
"""

# ╔═╡ cb57997e-c3ec-47e0-b9a0-b10aa9f5608d
# missing

# ╔═╡ bb2d06f8-1940-4585-961a-54068da50e91
md"""
Interpret the results. Ask yourself the following questions:

1. Explain why $S$ first increases and then decreases while $X$ only increases during the first 50 hours.
"""

# ╔═╡ 7c29c97d-dea5-4ea3-b5aa-bdccfe93939c
md"- Answer: missing"

# ╔═╡ 34018fc4-af3e-4f5a-9f24-a73293af2e85
md"""
2. What are the steady state values of $S$ and $X$.
"""

# ╔═╡ eebef095-ead8-4193-889b-53cbcab84514
md"- Answer: missing"

# ╔═╡ 9bb450c6-5499-42f6-8356-bdc4985b74e7
md"""
#### Part 2

Suppose that the substrate inlet concentration $S_{in}$ suddenly increases to $0.022\;g/L$ at $t = 100\;h$. Simulate the evolution of $S$ and $X$.
"""

# ╔═╡ 0298953a-b90c-41cd-8613-cb47ce752e43
md"""
Create the *condition* that contains the timepoint for the sudden change in $S_{in}$. Store it in `condition2`:
"""

# ╔═╡ c85e505d-99c6-4616-b7c1-42c05b4894fc
# condition2 = missing               # Uncomment and complete the instruction

# ╔═╡ 7c96bead-9f7b-4e84-abb9-9b6651208667
md"""
Make a new *reaction system* where the discrete event is included. Name it `fermenter_monod2`.
"""

# ╔═╡ 439dbeef-55b9-4fa4-aef5-fa0bb5a2ccf1
# @named fermenter_monod2 = missing     # Uncomment and complete the instruction

# ╔═╡ 43614c69-3fb5-4bef-b2a2-9805a5545fb8
md"""
Complete the new *reaction system*. Name it `fermenter_monod2_com`.
"""

# ╔═╡ 532306b5-2a71-4fcf-93ac-9dc52457a3f9
# fermenter_monod2_com = missing        # Uncomment and complete the instruction

# ╔═╡ 9e8f4500-0b6a-47f0-a1f5-a74daea9d117
md"""
Create the ODE problem and store it in `oprob2`:
"""

# ╔═╡ 7af72709-2f82-4971-8342-f02943f947c8
# oprob2 = missing                      # Uncomment and complete the instruction

# ╔═╡ e019f797-a6ad-4f8f-8f9e-69db00ed3c39
md"""
Solve the ODE problem. Make a deepcopy and use `Tsit5()` and `saveat=0.5`. Store the solution in `osol2`:
"""

# ╔═╡ 5f77450b-aa96-41b0-8017-a3d29fd7023a
# osol2 = missing               # Uncomment and complete the instruction

# ╔═╡ 310a78a5-94ce-4a29-b7a1-37831ce5c64e
md"""
Plot the results:
"""

# ╔═╡ 85742dc1-24cb-42d1-a70b-70f1be6b6c1e
# missing

# ╔═╡ 7f358845-c9fc-4e77-9882-94506f9338d6
md"""
Calculate the state state values for $S$ and $X$.
"""

# ╔═╡ 1102146d-abbe-43ce-9602-c863d1a91071
md"""
Inspect the final values in both the $S$ and $X$ vector.\
Tip: use something like: `(osol2[...][...], osol2[...][...])`
"""

# ╔═╡ 196edf3a-b220-4a54-8137-b136b509617e
# (osol2[...][...], osol2[...][...])   # Uncomment and complete the instruction

# ╔═╡ 3c7bf6f2-4ffd-4678-a930-94cf1322ba9f
md"""
Initialize a vector `u_guess2` with the final values for $S$ and $X$:
"""

# ╔═╡ 56eaa343-03b6-4cae-868b-e5c36ac66546
# u_guess2 = missing                    # Uncomment and complete the instruction

# ╔═╡ 4d8c05d8-ee29-4cb8-84cd-4342cb1db289
md"""
Initialize a vector `param_mod` with the parameter values. Notice that all parameter values will be the same, **except** the one of $S_{in}$.
"""

# ╔═╡ f121efc5-4e64-4e82-8672-2765ad85443e
# params_mod = missing                   # Uncomment and complete the instruction

# ╔═╡ 9fe054e0-cf21-49ba-a777-a8200b34b7dd
md"""
Make and solve the steady state problem. Call the output values `Seq2` and `Xeq2`.
"""

# ╔═╡ 0b992750-a446-447b-b2a1-26658c11c0bf
# Seq2, Xeq2 = missing                    # Uncomment and complete the instruction

# ╔═╡ 4ed59602-ad9d-4aae-8a21-83dabbfd3846
md"""
Inspect those values.
"""

# ╔═╡ 47a63fd8-f805-4c7d-8695-c9ee6550f24f
# missing                                 # Uncomment and complete the instruction

# ╔═╡ 4d962d0f-da41-405e-9438-733d5668cde3
md"""
Interpret the results. Ask yourself the following questions:

1. Can you clearly see the effect of the increase in $S_{in}$?
"""

# ╔═╡ 0bdf383c-e742-46ef-baba-11059e64f9c9
md"- Answer: missing"

# ╔═╡ d5e8e268-2dc1-4979-a800-9b733a7f7818
md"""
2. Find the steady state values of $S$ and $X$. Is the steady state value of $S$ influenced by the increase of $S_{in}$? Show how you can deduce that from the differential equations.
"""

# ╔═╡ 4c141769-65af-4821-a8d3-0e263eccaf8f
md"- Answer: missing"

# ╔═╡ 716362c9-54c4-49ef-bccc-d69425792c63
md"""
3. Can you explain why $X$ increased permanently?
"""

# ╔═╡ 8a0e6b54-bf50-4fc8-a7b9-fbfd2deb8d06
md"- Answer: missing"

# ╔═╡ 53980767-a84f-44f2-a878-2a7d57e0e2ae
md"""
#### Part 3

Suppose that the inlet/outlet flow $Q$ is suddenly doubled at $t = 100\;h$. Simulate the evolution of $S$ and $X$.
"""

# ╔═╡ 31b64d91-f8ff-413a-9c0c-402fe2215a81
md"""
Create the *condition* that contains the timepoint for the sudden change in $Q$. Store it in `condition3`:
"""

# ╔═╡ 6e771231-abf1-41c0-9aa8-ac7220f2a9cd
# condition3 = missing             # Uncomment and complete the instruction

# ╔═╡ 97b0de4f-c300-44f1-97fc-804d3263d8b5
md"""
Make a new *reaction system* where the discrete event is included. Name it `fermenter_monod3`.
"""

# ╔═╡ 837e27f1-a2d8-4d2c-aaa5-e82b6761e4fd
# @named fermenter_monod3 = missing    # Uncomment and complete the instruction

# ╔═╡ d8002843-03c6-4fa4-b8e5-b42eac27588c
md"""
Complete the new *reaction system*. Name it `fermenter_monod3_com`.
"""

# ╔═╡ 508f1dfe-3a92-4d80-b48d-e84a8738f97f
# fermenter_monod3_com =missing        # Uncomment and complete the instruction

# ╔═╡ 5d4c2573-4e57-455b-bdcc-1cee79b08ce2
md"""
Create the ODE problem and store it in `oprob3`:
"""

# ╔═╡ dd388e88-53af-48d3-800e-09b5c182a83b
# oprob3 = missing                      # Uncomment and complete the instruction

# ╔═╡ 5133d846-e6a6-4b50-9ce1-cb91cf04cbd1
md"""
Solve the ODE problem. Make a deepcopy and use `Tsit5()` and `saveat=0.5`. Store the solution in `osol3`:
"""

# ╔═╡ b3bc3348-a524-41e6-9fdb-865055246cd9
# osol3 = missing                  # Uncomment and complete the instruction

# ╔═╡ 39ef225d-3222-4998-bfd4-5ff88f74a0f9
md"""
Plot the results:
"""

# ╔═╡ 76849cf5-170b-452b-acdd-c4017feaad18
# missing                          # Uncomment and complete the instruction

# ╔═╡ 040f040f-aa28-442e-9f44-2a897e22ed4f
md"""
Interpret the results. Ask yourself the following questions:

1. Can you clearly see the effect of doubling of $Q$?
"""

# ╔═╡ 1b792c99-0165-49a3-8466-91082aa514bc
md"- Answer: missing"

# ╔═╡ f9e232c0-a74b-48b9-854d-f507717f8cdd
md"""
2. Can you argue, by means of reasoning, why $S$ increases and $X$ decreases?
"""

# ╔═╡ 7f773577-99f7-4aee-b53c-08cd9a25a236
md"- Answer: missing"

# ╔═╡ Cell order:
# ╠═2e58f4ae-f711-11ee-2598-7f3a6f2e2013
# ╠═539c1823-16d0-4355-97b6-fe9f0b106864
# ╠═e99680dc-73af-40aa-bf57-a06d3a7372be
# ╟─7856d878-8586-4cfd-9cf6-d61234450e41
# ╟─8500c35e-6bc3-4900-81bf-7705ddd61532
# ╠═e66518ee-b6f6-4cca-a224-30e01cffddbe
# ╟─f1350528-07a5-4860-ad2d-627588186abc
# ╠═331a34f4-89d4-4193-896c-c14ab0bf04e7
# ╟─55746566-2d46-4475-851a-02b7fad87a1a
# ╠═ec9cb3bd-f5ed-4ab0-9b3d-b875692227ac
# ╠═bd648109-f042-42de-9e0e-017b502fab95
# ╟─67117a27-dcea-4b43-b962-9ad9fd07f4f4
# ╟─d13e6e38-037e-4812-85e9-2c18bed360f6
# ╠═4b556cf0-8fad-434d-be56-dc1848d898ae
# ╟─ea55d648-7575-43c3-a385-5f4979996ef2
# ╠═1365c12e-e662-4858-983b-02ba94cd9f0d
# ╟─3941bd60-a83c-4f72-84b3-28e28cb845d0
# ╠═d6c1316a-cf96-43d1-854a-f25925cf4a55
# ╟─4926b941-c3b6-4804-b4a4-11e13e5186f2
# ╠═ab2a9842-6a9c-46bd-812b-db01629d6a1c
# ╟─b6a526bd-6ee5-442b-9fb8-3fbe1e280dd4
# ╠═1f62e66d-571f-41ca-9f02-f36a8ca10ab9
# ╟─37ced6e3-b435-4546-a720-a1ec1af23a65
# ╠═d609bed4-94cf-4167-80fe-924501a5835c
# ╟─3c638fc4-8ebc-4e26-984e-b4513035287e
# ╠═25ef069e-5a08-4e85-acea-a5b64e0890f6
# ╟─d8234337-7516-4817-8ce5-194af694e3e3
# ╟─b63fa707-6023-44fc-aac9-25c76e075f90
# ╠═8d96d79f-6ddb-4bba-a6ea-821588e13107
# ╟─e8969045-27ac-460e-86a2-7494903534e8
# ╠═1777503e-b793-4be2-b80b-b4edcd7041b5
# ╟─1d5a2118-96e9-49b2-9931-5d4b201cb8f5
# ╠═cb57997e-c3ec-47e0-b9a0-b10aa9f5608d
# ╟─bb2d06f8-1940-4585-961a-54068da50e91
# ╟─7c29c97d-dea5-4ea3-b5aa-bdccfe93939c
# ╟─34018fc4-af3e-4f5a-9f24-a73293af2e85
# ╟─eebef095-ead8-4193-889b-53cbcab84514
# ╟─9bb450c6-5499-42f6-8356-bdc4985b74e7
# ╟─0298953a-b90c-41cd-8613-cb47ce752e43
# ╠═c85e505d-99c6-4616-b7c1-42c05b4894fc
# ╟─7c96bead-9f7b-4e84-abb9-9b6651208667
# ╠═439dbeef-55b9-4fa4-aef5-fa0bb5a2ccf1
# ╟─43614c69-3fb5-4bef-b2a2-9805a5545fb8
# ╠═532306b5-2a71-4fcf-93ac-9dc52457a3f9
# ╟─9e8f4500-0b6a-47f0-a1f5-a74daea9d117
# ╠═7af72709-2f82-4971-8342-f02943f947c8
# ╟─e019f797-a6ad-4f8f-8f9e-69db00ed3c39
# ╠═5f77450b-aa96-41b0-8017-a3d29fd7023a
# ╟─310a78a5-94ce-4a29-b7a1-37831ce5c64e
# ╠═85742dc1-24cb-42d1-a70b-70f1be6b6c1e
# ╟─7f358845-c9fc-4e77-9882-94506f9338d6
# ╟─1102146d-abbe-43ce-9602-c863d1a91071
# ╠═196edf3a-b220-4a54-8137-b136b509617e
# ╟─3c7bf6f2-4ffd-4678-a930-94cf1322ba9f
# ╠═56eaa343-03b6-4cae-868b-e5c36ac66546
# ╟─4d8c05d8-ee29-4cb8-84cd-4342cb1db289
# ╠═f121efc5-4e64-4e82-8672-2765ad85443e
# ╟─9fe054e0-cf21-49ba-a777-a8200b34b7dd
# ╠═0b992750-a446-447b-b2a1-26658c11c0bf
# ╟─4ed59602-ad9d-4aae-8a21-83dabbfd3846
# ╠═47a63fd8-f805-4c7d-8695-c9ee6550f24f
# ╟─4d962d0f-da41-405e-9438-733d5668cde3
# ╟─0bdf383c-e742-46ef-baba-11059e64f9c9
# ╟─d5e8e268-2dc1-4979-a800-9b733a7f7818
# ╟─4c141769-65af-4821-a8d3-0e263eccaf8f
# ╟─716362c9-54c4-49ef-bccc-d69425792c63
# ╟─8a0e6b54-bf50-4fc8-a7b9-fbfd2deb8d06
# ╟─53980767-a84f-44f2-a878-2a7d57e0e2ae
# ╟─31b64d91-f8ff-413a-9c0c-402fe2215a81
# ╠═6e771231-abf1-41c0-9aa8-ac7220f2a9cd
# ╟─97b0de4f-c300-44f1-97fc-804d3263d8b5
# ╠═837e27f1-a2d8-4d2c-aaa5-e82b6761e4fd
# ╟─d8002843-03c6-4fa4-b8e5-b42eac27588c
# ╠═508f1dfe-3a92-4d80-b48d-e84a8738f97f
# ╟─5d4c2573-4e57-455b-bdcc-1cee79b08ce2
# ╠═dd388e88-53af-48d3-800e-09b5c182a83b
# ╟─5133d846-e6a6-4b50-9ce1-cb91cf04cbd1
# ╠═b3bc3348-a524-41e6-9fdb-865055246cd9
# ╟─39ef225d-3222-4998-bfd4-5ff88f74a0f9
# ╠═76849cf5-170b-452b-acdd-c4017feaad18
# ╟─040f040f-aa28-442e-9f44-2a897e22ed4f
# ╟─1b792c99-0165-49a3-8466-91082aa514bc
# ╟─f9e232c0-a74b-48b9-854d-f507717f8cdd
# ╟─7f773577-99f7-4aee-b53c-08cd9a25a236
