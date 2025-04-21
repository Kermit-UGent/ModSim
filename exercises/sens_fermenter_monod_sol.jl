### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 3ef93246-657d-4e77-9bf0-8380c64bfcfd
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate("..")
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
using ForwardDiff

# ╔═╡ 31d294d1-3a1f-41db-abff-54f2a67c7ed9
md"""
### Exercise: Fermenter - Monod kinetics - Sensitivity analysis
"""

# ╔═╡ 5ffe7dcb-620d-4f22-95fe-2f77cda6fbe7
md"""
In one of the previous practicals we were introduced to a fermenter in which biomass $X$ [$g/L$] grows by breaking down substrate $S$ [$g/L$]. The reactor is fed with an inlet flow rate $Q_{in}$ [$L/h$], which consists of a (manipulable) input concentration of substrate $S_{in}$ [$g/L$]. This process was modelled using Monod kinetics, resulting in the model below:

$$\begin{eqnarray*}
S + X \xrightarrow[\quad\quad]{k} (1 + Y) \, X \quad\quad\quad\quad \textrm{with} \quad k = \cfrac{\mu_{max}}{S + K_s}
\end{eqnarray*}$$
"""

# ╔═╡ 6ec6da23-853b-4129-94cf-67b5cadb1f95
md"""
The *reaction network object* model for this problem could be defined as:
"""

# ╔═╡ 935ca610-7a7a-4692-8908-fc26abb880b4
fermenter_monod = @reaction_network begin
    μmax/(S+Ks), S + X --> (1 + Y)*X
    # Alternatives:
    # mm(S, μmax, Ks)*X, S => Y*X
	# mm(S, μmax, Ks)*X, S + X => (1 + Y)*X
    Q/V, (S, X) --> 0
    Q/V*Sin, 0 --> S
end

# ╔═╡ 79e6056a-881c-442f-8989-5bc284d3d777
md"""
which resulted in the following differential equations:
"""

# ╔═╡ fa93e2c3-8b43-418e-ba24-406645b2e397
md"""
$$\begin{eqnarray*}
\cfrac{dS}{dt} &=& \cfrac{Q}{V} \left(S_{in} - S \right) - \mu_{max}\cfrac{S}{S + K_s} X\\
\cfrac{dX}{dt} &=& -\cfrac{Q}{V} X + Y \mu_{max}\cfrac{S}{S + K_s} X
\end{eqnarray*}$$
"""

# ╔═╡ 06730f54-7293-43f2-b772-84eec3e5528a
md"""
Convert the system to a symbolic differential equation model and verify, by analyzing the differential equation, that your model is correctly implemented.

In case you want to use the `mm` function, keep in mind that `mm(S, μmax, Ks)` stands for $\mu_{max} \, \cfrac{S}{S + K_s}$.
"""

# ╔═╡ 7f8b7a2e-bc65-4b51-ad59-bd7ac98604dd
# osys = missing                    # Uncomment and complete the instruction
osys = convert(ODESystem, fermenter_monod)

# ╔═╡ 55f1d688-0c53-481b-9965-5e92ca87ad83
md"""
The parameter values are $\mu_{max} = 0.40$, $K_s = 0.015$, $Y = 0.67$, $Q = 2.0$, $V = 40.0$ and $S_{in} = 0.022\;g/L$. Suppose that at $t=0$ no substrate $S$ is present in the reactor but that there is initially some biomass with a concetration of $0.0005\;g/L$.\
Compute the following in a timespan of $[0, 100]\,h$:

- The sensitivities of $S$ and $X$ wrt. $\mu_{max}$, $K_s$ and $S_{in}$.

Plot the following:
- A figure with the sensitivity functions of $S$ and $X$ wrt. $S_{in}$.
- A figure with the sensitivity functions of $S$ wrt. $\mu_{max}$, $K_s$ and $S_{in}$.
- A figure with the sensitivity functions of $X$ wrt. $\mu_{max}$, $K_s$ and $S_{in}$.

Interpret your results.
"""

# ╔═╡ 6c4e3c09-4b84-4f5c-8739-2ac18e6f2af6
md"""
Initialize a vector `u0` with the initial conditions, set the timespan and initialize a vector `param` with the parameter values:
"""

# ╔═╡ 2ee277e5-ce4a-4ade-be0e-9bba7a4dc08c
# u0 = missing                 # Uncomment and complete the instruction
u0 = [:S => 0.0, :X => 0.0005]

# ╔═╡ 3fdc6b17-cdeb-4dc5-8886-9d3a62caac8d
# tspan = missing              # Uncomment and complete the instruction
tspan = (0.0, 100.0)

# ╔═╡ 0139da85-02e3-4021-9b39-84af7e68d428
md"""
For the sake of clarity, we will use the variables `μmax`, `Ks` and `Sin` to store the parameter values that are used for the calculation of the sensitivity functions.
"""

# ╔═╡ 0f995929-4d2b-4a7a-8da1-04e4d501385f
# μmax = missing               # Uncomment and complete the instruction
μmax = 0.40

# ╔═╡ 262e8346-df6d-49bf-9186-92f5afb421e0
# Ks = missing                 # Uncomment and complete the instruction
Ks = 0.015

# ╔═╡ baa777d2-abb8-45f8-87aa-b3b17c8dc07c
# Sin = missing                # Uncomment and complete the instruction
Sin = 0.022

# ╔═╡ 79b0eb65-5a0f-40b3-aa97-4088421c562e
# params = missing             # Uncomment and complete the instruction
params = [:μmax => μmax, :Ks => Ks, :Y => 0.67, :Q => 2, :V => 40, :Sin => Sin]

# ╔═╡ f0f4fa14-6f99-4f21-a743-be61e08444a7
md"""
Create the ODE problem and store it in `oprob`. Next, solve the ODE problem using `Tsit5()` and `saveat=0.5`, and store the solution in `osol`. Finally plot the results.
"""

# ╔═╡ 8b2f23f6-80b2-4e63-942e-e5cd17d8ba72
# oprob = missing;              # Uncomment and complete the instruction
oprob = ODEProblem(fermenter_monod, u0, tspan, params);

# ╔═╡ 89a31c32-88a4-479f-a688-ffcb75ee8e91
# osol = missing                # Uncomment and complete the instruction
osol = solve(oprob, Tsit5(), saveat=0.5)

# ╔═╡ 51a9b7e6-8ad9-477d-9596-ffd614df2c79
# missing                        # Uncomment and complete the instruction
plot(osol)

# ╔═╡ 693844d0-3858-4861-bae0-b47e78809f17
md"""
Write a solution function with as argument a vector of the parameters (with the values for which we want to calculate the sensitivity), and that returns the outputs.
"""

# ╔═╡ 9622f7ca-f71a-4ad9-a309-d7d10a1c3e3b
# Uncomment and complete the instruction
# function fermenter_monod_sim(params)
	# μmax, Ks, Sin = missing
	# u0 = missing
	# tspan = missing
	# params = missing
	# oprob = missing;
	# osol = missing
	# return missing
# end
function fermenter_monod_sim(params)
	μmax, Ks, Sin = params
	u0 = [:S => 0.0, :X => 0.0005]
	tspan = (0.0, 100.0)
	params = [:μmax => μmax, :Ks => Ks, :Y => 0.67, :Q => 2, :V => 40, :Sin => Sin]
	oprob = ODEProblem(fermenter_monod, u0, tspan, params, combinatoric_ratelaws=false)
	osol = solve(oprob, Tsit5(), saveat=0.5)
	return osol
end

# ╔═╡ 4a5971b1-f4d0-43b6-805f-e17f5052ae92
md"""
Make two functions based on the solution function, where each returns a single output; hence, one function that returns the output $S$, and another function that returns the output $X$.
"""

# ╔═╡ f40c6402-3c28-4a7d-b629-83507a9f29bd
# fermenter_monod_sim_S(params) = missing # Uncomment and complete the instruction
fermenter_monod_sim_S(params) = fermenter_monod_sim(params)[:S]

# ╔═╡ 3ae5bd00-2e06-4789-aab3-d897824d5e29
# fermenter_monod_sim_X(params) = missing # Uncomment and complete the instruction
fermenter_monod_sim_X(params) = fermenter_monod_sim(params)[:X]

# ╔═╡ 4bd2bcca-9c42-4333-b062-2aaa9f7be3fe
md"""
Make the time vector.
"""

# ╔═╡ fbd98975-aa32-46ae-8db0-0e65cdf48309
# t_vals = missing              # Uncomment and complete the instruction
t_vals = 0:0.5:100.0

# ╔═╡ 93791eb3-1eaa-4146-90b5-c4811fb3485b
md"""
Compute the two outputs $S$ and $X$ for the given parameter values.
"""

# ╔═╡ dc0557d6-81b9-4759-8ed7-3129f60c6dc3
# S_sim = missing               # Uncomment and complete the instruction
S_sim = fermenter_monod_sim_S([μmax, Ks, Sin])

# ╔═╡ 95bc683c-f6e6-4b42-b90b-b5a863edd4d5
# X_sim = missing               # Uncomment and complete the instruction
X_sim = fermenter_monod_sim_X([μmax, Ks, Sin])

# ╔═╡ fa970c0e-fb3b-486f-bbc1-345d44f8f0da
md"""
Using `ForwardDiff.jacobian` to compute the sensitivities for the single ouputs $S$ and $X$. Hence, you need to call `ForwardDiff.jacobian` twice.
"""

# ╔═╡ 64354302-f4cc-4592-9302-5db0f5bccb2e
# sens_S = missing                # Uncomment and complete the instruction
sens_S = ForwardDiff.jacobian(fermenter_monod_sim_S, [μmax, Ks, Sin])

# ╔═╡ 49a94b9a-a543-495e-b4f1-c8579e59304d
# sens_X = missing                # Uncomment and complete the instruction
sens_X = ForwardDiff.jacobian(fermenter_monod_sim_X, [μmax, Ks, Sin])

# ╔═╡ 9cace6c1-e678-4dd7-8705-92a55eb32fa9
md"""
Extract the (absolute) sensitivities of the outputs on the different parameters.
"""

# ╔═╡ a6dc2b60-6a0a-4140-892e-02cde8dc79d3
# Uncomment and complete the instruction
# begin
# 	sens_S_on_μmax = missing
# 	sens_S_on_Ks   = missing
# 	sens_S_on_Sin  = missing
# end;
begin
	sens_S_on_μmax = sens_S[:, 1]
	sens_S_on_Ks   = sens_S[:, 2]
	sens_S_on_Sin  = sens_S[:, 3]
end;

# ╔═╡ f806c243-9032-46b7-add3-4714344691c7
# Uncomment and complete the instruction
# begin
# 	sens_X_on_μmax = missing
# 	sens_X_on_Ks   = missing
# 	sens_X_on_Sin  = missing
# end;
begin
	sens_X_on_μmax = sens_X[:, 1]
	sens_X_on_Ks   = sens_X[:, 2]
	sens_X_on_Sin  = sens_X[:, 3]
end;

# ╔═╡ 5bf3a62d-d2aa-4653-8ee8-e90caa9504e8
md"""
Compute the normalized sensitivities.
"""

# ╔═╡ 76846731-929c-408f-a3de-970581c497e9
# Uncomment and complete the instruction
# begin
# 	sens_S_on_μmax_rel = missing
# 	sens_S_on_Ks_rel   = missing
# 	sens_S_on_Sin_rel  = missing
# end
begin
	sens_S_on_μmax_rel = sens_S_on_μmax .* μmax ./ S_sim
	sens_S_on_Ks_rel   = sens_S_on_Ks .* Ks ./ S_sim
	sens_S_on_Sin_rel  = sens_S_on_Sin .* Sin ./ S_sim
end;

# ╔═╡ b6c57444-547c-4e82-8526-6a30566e07c5
# Uncomment and complete the instruction
# begin
# 	sens_X_on_μmax_rel = missing
# 	sens_X_on_Ks_rel   = missing
# 	sens_X_on_Sin_rel  = missing
# end;
begin
	sens_X_on_μmax_rel = sens_X_on_μmax .* μmax ./ X_sim
	sens_X_on_Ks_rel   = sens_X_on_Ks .* Ks ./ X_sim
	sens_X_on_Sin_rel  = sens_X_on_Sin .* Sin ./ X_sim
end;

# ╔═╡ 5388c2a7-5a11-4da8-be09-46045cde8a4e
md"""
Plot the sensitivity functions of $S$ and $X$ wrt. $S_{in}$. Provide a suitable title (`title="..."`), labels (`label=["..." "..."]`) and an x-label (`xlabel="..."`), and set the line width to 2 (`linewidth=...`).
"""

# ╔═╡ db840c76-a6c6-49fb-a0bb-d9149f947bc0
# missing               # Uncomment and complete the instruction
plot(t_vals, [sens_S_on_Sin_rel, sens_X_on_Sin_rel], title="Normalized sensitivities", label=["S on Sin" "X on Sin"], xlabel="Time (hours)", linewidth=2)

# ╔═╡ 02314e9f-6186-4c97-9b28-03a4a1a15668
maximum(sens_X_on_Sin_rel)

# ╔═╡ d41375ef-6958-4705-a417-4c6a491232ee
md"""
Interpret your results. Try to answer the following question(s):
"""

# ╔═╡ f6be14f3-e0d0-41dc-bc7b-54174b1ea5ea
md"""
!!! question
	Which output variable, $S$ or $X$, is most sensitive to $S_{in}$ in steady state?
"""

# ╔═╡ b2583ee2-412b-4a0a-a4d0-fe476e0510e5
md"- Answer: missing"
#=
- X is most sensitive wrt. Sin because the higher Sin, the more X can be produced.
=#

# ╔═╡ dbbf694b-942a-4757-a255-f67b208ba03b
md"""
!!! question
	Why is the sensitivity function of $S$ wrt. $S_{in}$ at first positive but then becomes zero?
"""

# ╔═╡ f10ac0f2-6b14-4805-b649-56c63f5b523a
md"- Answer: missing"
#=
- In the beginning the substrate S is positively affected by Sin because S enters the tank through Sin and only little biomass X is present, so the biomass cannot consume the substrate very fast.
=#

# ╔═╡ be89600a-4927-4afc-9813-d8a70adb2852
md"""
Plot the sensitivity functions of $S$ wrt. $\mu_{max}$, $K_s$ and $S_{in}$. Provide a suitable title (`title="..."`), labels (`label=["..." "..." "..."]`) and an x-label (`xlabel="..."`), and set the line width to 2 (`linewidth=...`).
"""

# ╔═╡ c0223da4-9959-48d0-b607-633b2e82986c
# missing                  # Uncomment and complete the instruction
plot(t_vals, [sens_S_on_μmax_rel, sens_S_on_Ks_rel, sens_S_on_Sin_rel], title="Normalized sensitivities", label=["S on μmax" "S on Ks" "S on Sin"], xlabel="Time (hours)", linewidth=2)

# ╔═╡ 260ede0e-2584-4fe9-ae92-69990ca8f854
md"""
Interpret your results. Try to answer the following question(s):
"""

# ╔═╡ ff86a29f-9308-473b-aa1c-dfd4af8179c7
md"""
!!! question
	Which parameter, $\mu_{max}$, $K_s$ or $S_{in}$, affects the output $S$ the most in steady state? Why is this?
"""

# ╔═╡ 555228b2-8075-49ac-a7ef-f51aa51dd95d
md"- Answer: missing"
#=
- It seems like μmax is affecting S the most in steady state and its influence is negative, hence, the larger μmax, the smaller S.
=#

# ╔═╡ 6704e43b-30bc-4168-93db-fb853c9761fa
md"""
!!! question
	Why is the sensitivity function of $S$ wrt. $\mu_{max}$ negative?
"""

# ╔═╡ 4d2094d8-98b3-4ac9-9612-942b5cd18ce0
md"- Answer: missing"
#=
- The sensitivity function of S on μmax is negative, because the larger μmax, the larger the reaction rate r (S => Y*X). Hence, more X will be produced, so more S will be consumed. Remember that r = μmax*S*X/(S + Ks) and μmax is in the numerator.
=#

# ╔═╡ dc8ef7cb-ae97-401c-9a25-12439e4363c6
md"""
!!! question
	Why is the sensitivity function of $S$ wrt. $K_s$ positive? Does this correspond to the meaning of the half-saturation constant $K_s$ or substrate concentration that allows to achieve half the maximum growth rate? *In other words: do higher values of $K_s$ support growth at low substrate concentrations?*
"""

# ╔═╡ 66712fb2-50ea-41be-bcca-305319d158e9
md"- Answer: missing"
#=
- The sensitivity function of S on Ks is positive, because the larger Ks, the smaller the reaction rate r (S => Y*X). Hence, less X will be produced, so less S will be consumed. Remember that r = μmax*S*X/(S + Ks) and Ks is in the denominator.
=#

# ╔═╡ 16a84fdb-8ce2-45b9-bfb7-7f4e1284a1d7
md"""
Plot the sensitivity functions of $X$ on $\mu_{max}$, $K_s$ and $S_{in}$. Provide a suitable title (`title="..."`), labels (`label=["..." "..." "..."]`) and an x-label (`xlabel="..."`), and set the line width to 2 (`linewidth=...`).
"""

# ╔═╡ 53134149-0bf7-41c1-9b35-e5037744211f
# missing
plot(t_vals, [sens_X_on_μmax_rel, sens_X_on_Ks_rel, sens_X_on_Sin_rel], title="Normalized sensitivities", label=["X on μmax" "X on Ks" "X on Sin"], xlabel="Time (hours)", linewidth=2)

# ╔═╡ 0355bbf6-853d-4bd8-b32f-98fdbd2b761c
md"""
Interpret your results. Try to answer the following question(s):
"""

# ╔═╡ 355ca6a7-466b-4969-ab48-28e2257f9810
md"""
!!! question
	Which parameter, $\mu_{max}$, $K_s$ or $S_{in}$, affects the output $X$ the most in steady state?
"""

# ╔═╡ 4e08baf1-eeb6-49be-9751-204dd9ae1103
md"- Answer: missing"
#=
- It seems like Sin is affecting X the most in steady state and its influence is positive, hence, the larger Sin, the larger X.
=#

# ╔═╡ 182e0c2b-9aa4-4514-96e8-b6a9f53b371b
md"""
!!! question
	Why is the sensitivity function of $X$ wrt. $\mu_{max}$ positive? How does this compare to substrate $S$?
"""

# ╔═╡ 05ec320e-ac87-45af-904e-5d69639e1f27
md"- Answer: missing"
#=
- The sensitivity function of X on μmax is positive, because the larger μmax, the larger the reaction rate r (S => Y*X). Hence, more X will be produced. Remember that r = μmax*S*X/(S + Ks) and μmax is in the numerator.
=#

# ╔═╡ 22e98c77-0e9f-444a-a2f0-b940940996b7
md"""
!!! question
	Why is the sensitivity function of $X$ wrt. $K_s$ negative? Compare it to that of substrate $S$.
"""

# ╔═╡ 91795413-c9a2-43cb-a7b0-d9e055e1cf94
md"- Answer: missing"
#=
- The sensitivity function of X on Ks is negative, because the larger Ks, the smaller the reaction rate r (S => Y*X). Hence, less X will be produced. Remember that r = μmax*S*X/(S + Ks) and Ks is in the denominator.
=#

# ╔═╡ Cell order:
# ╠═55cdebd2-0881-11ef-2722-91de1447877a
# ╠═03ae0690-06a0-4276-9f00-d07b206fe124
# ╠═3ef93246-657d-4e77-9bf0-8380c64bfcfd
# ╠═a355b0ba-baaf-49f4-a5dc-965364a884f0
# ╠═00fd6d49-f561-42e9-9413-d33af92f83dc
# ╠═7ae714c4-d25d-4f9f-ab3d-cc067db9c156
# ╟─31d294d1-3a1f-41db-abff-54f2a67c7ed9
# ╟─5ffe7dcb-620d-4f22-95fe-2f77cda6fbe7
# ╟─6ec6da23-853b-4129-94cf-67b5cadb1f95
# ╠═935ca610-7a7a-4692-8908-fc26abb880b4
# ╟─79e6056a-881c-442f-8989-5bc284d3d777
# ╟─fa93e2c3-8b43-418e-ba24-406645b2e397
# ╟─06730f54-7293-43f2-b772-84eec3e5528a
# ╠═7f8b7a2e-bc65-4b51-ad59-bd7ac98604dd
# ╟─55f1d688-0c53-481b-9965-5e92ca87ad83
# ╟─6c4e3c09-4b84-4f5c-8739-2ac18e6f2af6
# ╠═2ee277e5-ce4a-4ade-be0e-9bba7a4dc08c
# ╠═3fdc6b17-cdeb-4dc5-8886-9d3a62caac8d
# ╟─0139da85-02e3-4021-9b39-84af7e68d428
# ╠═0f995929-4d2b-4a7a-8da1-04e4d501385f
# ╠═262e8346-df6d-49bf-9186-92f5afb421e0
# ╠═baa777d2-abb8-45f8-87aa-b3b17c8dc07c
# ╠═79b0eb65-5a0f-40b3-aa97-4088421c562e
# ╟─f0f4fa14-6f99-4f21-a743-be61e08444a7
# ╠═8b2f23f6-80b2-4e63-942e-e5cd17d8ba72
# ╠═89a31c32-88a4-479f-a688-ffcb75ee8e91
# ╠═51a9b7e6-8ad9-477d-9596-ffd614df2c79
# ╟─693844d0-3858-4861-bae0-b47e78809f17
# ╠═9622f7ca-f71a-4ad9-a309-d7d10a1c3e3b
# ╟─4a5971b1-f4d0-43b6-805f-e17f5052ae92
# ╠═f40c6402-3c28-4a7d-b629-83507a9f29bd
# ╠═3ae5bd00-2e06-4789-aab3-d897824d5e29
# ╟─4bd2bcca-9c42-4333-b062-2aaa9f7be3fe
# ╠═fbd98975-aa32-46ae-8db0-0e65cdf48309
# ╟─93791eb3-1eaa-4146-90b5-c4811fb3485b
# ╠═dc0557d6-81b9-4759-8ed7-3129f60c6dc3
# ╠═95bc683c-f6e6-4b42-b90b-b5a863edd4d5
# ╟─fa970c0e-fb3b-486f-bbc1-345d44f8f0da
# ╠═64354302-f4cc-4592-9302-5db0f5bccb2e
# ╠═49a94b9a-a543-495e-b4f1-c8579e59304d
# ╟─9cace6c1-e678-4dd7-8705-92a55eb32fa9
# ╠═a6dc2b60-6a0a-4140-892e-02cde8dc79d3
# ╠═f806c243-9032-46b7-add3-4714344691c7
# ╟─5bf3a62d-d2aa-4653-8ee8-e90caa9504e8
# ╠═76846731-929c-408f-a3de-970581c497e9
# ╠═b6c57444-547c-4e82-8526-6a30566e07c5
# ╟─5388c2a7-5a11-4da8-be09-46045cde8a4e
# ╠═db840c76-a6c6-49fb-a0bb-d9149f947bc0
# ╠═02314e9f-6186-4c97-9b28-03a4a1a15668
# ╟─d41375ef-6958-4705-a417-4c6a491232ee
# ╟─f6be14f3-e0d0-41dc-bc7b-54174b1ea5ea
# ╠═b2583ee2-412b-4a0a-a4d0-fe476e0510e5
# ╟─dbbf694b-942a-4757-a255-f67b208ba03b
# ╠═f10ac0f2-6b14-4805-b649-56c63f5b523a
# ╟─be89600a-4927-4afc-9813-d8a70adb2852
# ╠═c0223da4-9959-48d0-b607-633b2e82986c
# ╟─260ede0e-2584-4fe9-ae92-69990ca8f854
# ╟─ff86a29f-9308-473b-aa1c-dfd4af8179c7
# ╠═555228b2-8075-49ac-a7ef-f51aa51dd95d
# ╟─6704e43b-30bc-4168-93db-fb853c9761fa
# ╠═4d2094d8-98b3-4ac9-9612-942b5cd18ce0
# ╟─dc8ef7cb-ae97-401c-9a25-12439e4363c6
# ╠═66712fb2-50ea-41be-bcca-305319d158e9
# ╟─16a84fdb-8ce2-45b9-bfb7-7f4e1284a1d7
# ╠═53134149-0bf7-41c1-9b35-e5037744211f
# ╟─0355bbf6-853d-4bd8-b32f-98fdbd2b761c
# ╟─355ca6a7-466b-4969-ab48-28e2257f9810
# ╠═4e08baf1-eeb6-49be-9751-204dd9ae1103
# ╟─182e0c2b-9aa4-4514-96e8-b6a9f53b371b
# ╠═05ec320e-ac87-45af-904e-5d69639e1f27
# ╟─22e98c77-0e9f-444a-a2f0-b940940996b7
# ╠═91795413-c9a2-43cb-a7b0-d9e055e1cf94
