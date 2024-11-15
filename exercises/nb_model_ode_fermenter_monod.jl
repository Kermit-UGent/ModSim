### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ e99680dc-73af-40aa-bf57-a06d3a7372be
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate(".")
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
md"
In a fermenter reactor biomass grows on substrate. The reactor is fed with a inlet flow rate $Q_{in}$ [$L/h$], which consist of a (manipulable) input concentration of substrate $S_{in}$ [$g/L$]. Inside the reactor, biomass, with a concentration of $X$ [$g/L$], is produced through **Monod** kinetics:

$$\begin{eqnarray*}
%S  \xrightarrow[\quad\quad]{\beta} Y \, X
S \xrightarrow[\quad\quad]{r} Y \, X \quad\quad\quad\quad r = \mu \, X
\end{eqnarray*}$$

where

$$\mu = \mu_{max} \, \cfrac{S}{S + K_s}$$

is called the specific growth rate [$h^{-1}$]. Therein, $\mu_{max}$ is the maximum speficic growth rate, and $K_s$ [$g/L$] is the so-called *half-velocity constant* (i.e. the value of $S$ when $\mu/\mu_{max} = 0.5$). Futhermore, $Y$ [$gX/gS$] is the yield coefficient which is defined here by the amount of produced biomass by consumption of one unit of substrate. The reactor is drained with an outlet flow $Q$ [$L/h$], which consist of the current concentrations of substrate $S$ [$g/L$] and biomass $X$ [$g/L$] inside the reactor. The volume $V$ [$L$] of the reactor content is kept constant by setting $Q_{in} = Q$.
"

# ╔═╡ f1350528-07a5-4860-ad2d-627588186abc
md"
Create a *reaction network object* model for the aforementioned problem in order to simulate the evolution of substrate $S$ and biomass $X$ with time. Name it `fermenter_firstorder`.

Tip: The specific growth rate $\mu = \mu_{max} \, \cfrac{S}{S + K_s}$ can be implemented with `mm(S, μmax, Ks)`. The function `mm` stands for the Michaelis-Menten kinetics, whcih is equivalent to Monod kinetics.
"


# ╔═╡ 331a34f4-89d4-4193-896c-c14ab0bf04e7
# fermenter_monod = @reaction_network begin
#     ...        # Y*X is created from one S at a rate mm(S, μmax, Ks)
#     ...        # S is created at a rate Q/V*Sin
#     ...        # S and X are degraded at a rate Q/V*S
# end
fermenter_monod = @reaction_network begin
	# Y*X is created from one S at a rate X * mm(S, μmax, Ks)
	# or, when S and X meet, this results in Y*X and X
    # X * mm(S, μmax, Ks), S --> Y*X
	mm(S, μmax, Ks), S + X --> (1 + Y)*X
    Q/V, (S, X) --> 0               # S and X are degraded at a rate Q/V*S
    Q/V*Sin, 0 --> S                # S is created at a rate Q/V*Sin 
end

# ╔═╡ 42190228-40d3-48e8-b52f-156a0c7cbddc
parameters(fermenter_monod)

# ╔═╡ 55746566-2d46-4475-851a-02b7fad87a1a
md"
Convert the system to a symbolic differential equation model and verify, by analyzing the differential equation, that your model is correctly implemented.

Keep in mind that `mm(S, μmax, Ks)` stands for $\mu_{max} \, \cfrac{S}{S + K_s}$.
"

# ╔═╡ ec9cb3bd-f5ed-4ab0-9b3d-b875692227ac
# osys = ...
osys = convert(ODESystem, fermenter_monod)

# ╔═╡ 67117a27-dcea-4b43-b962-9ad9fd07f4f4
md"
The parameter values are $\mu_{max} = 0.30$, $K_s = 0.15$, $Y = 0.80$, $Q = 2.0$, $V = 40.0$ and $S_{in} = 2.2\;g/L$. Suppose that at $t=0\;h$ no substrate $S$ is present in the reactor but that there is initially some biomass with a concetration of $0.01\;g/L$. Simulate the evolution of $S$ and $X$ during 200 hours.
"

# ╔═╡ d13e6e38-037e-4812-85e9-2c18bed360f6
md"
Initialize a vector `u₀` with the initial conditions:
"

# ╔═╡ 4b556cf0-8fad-434d-be56-dc1848d898ae
# u0 = ...            # Uncomment and complete the instruction
u0 = [:S => 0.0, :X => 0.01]

# ╔═╡ ea55d648-7575-43c3-a385-5f4979996ef2
md"
Set the timespan for the simulation:
"

# ╔═╡ 1365c12e-e662-4858-983b-02ba94cd9f0d
# tspan = ...         # Uncomment and complete the instruction
tspan = (0.0, 200.0)

# ╔═╡ 3941bd60-a83c-4f72-84b3-28e28cb845d0
md"
Initialize a vector `param` with the parameter values:
"

# ╔═╡ d6c1316a-cf96-43d1-854a-f25925cf4a55
# params = ...         # Uncomment and complete the instruction
params = [:μmax => 0.30, :Ks => 0.15, :Y => 0.80, :Q => 2, :V => 40, :Sin => 2.2]

# ╔═╡ 4926b941-c3b6-4804-b4a4-11e13e5186f2
md"
Create the ODE problem and store it in `oprob`:
"

# ╔═╡ ab2a9842-6a9c-46bd-812b-db01629d6a1c
# oprob = ...           # Uncomment and complete the instruction
oprob = ODEProblem(fermenter_monod, u0, tspan, params, combinatoric_ratelaws=false)

# ╔═╡ b6a526bd-6ee5-442b-9fb8-3fbe1e280dd4
md"
#### Part 1

Solve the ODE problem. Use `Tsit5()` and `saveat=0.5`. Store the solution in `osol1`:
"

# ╔═╡ 1f62e66d-571f-41ca-9f02-f36a8ca10ab9
# osol1 = ...           # Uncomment and complete the instruction
osol1 = solve(oprob, Tsit5(), saveat=0.5)

# ╔═╡ 37ced6e3-b435-4546-a720-a1ec1af23a65
md"
Plot the results:
"

# ╔═╡ d609bed4-94cf-4167-80fe-924501a5835c
# ...
plot(osol1)

# ╔═╡ bb2d06f8-1940-4585-961a-54068da50e91
md"
Interpret the results. Ask yourself the following questions:

1. Explain why $S$ first increases and then decreases while $X$ only increases during the first 50 hours.
- Answer: ...
2. Find the steady state values of $S$ and $X$.
- Answer: ...
"
#=
1. At the start there is no substrate present in the reactor vessel. So, the substrate concentration basically only increases. At the start, the amount of biomass X is very small, so the convertion of S to X will be little and furthermore X is continuously being drained due to the outlet flow. Once S and X become substantially larger, the the draining of S due to the outlet flow will be larger and more S will be converted to X. Hence, S will first increase and then decrease, while X will continue to increase before reaching equilibrium.
=#

# ╔═╡ 9bb450c6-5499-42f6-8356-bdc4985b74e7
md"
#### Part 2

Suppose that the substrate inlet concentration $S_{in}$ suddenly increases to $2.8\;g/L$ at $t = 100\;h$. Simulate the evolution of $S$ and $X$.
"

# ╔═╡ 0298953a-b90c-41cd-8613-cb47ce752e43
md"
Create the *condition* that contains the timepoint for the sudden change in $S_{in}$. Store it in `condition2`:
"

# ╔═╡ c44e0ae9-addb-4518-a3ee-559a24852c3a
# condition = ...           # Uncomment and complete the instruction
condition2 = [100.0]


# ╔═╡ 36a3b978-57ff-4d03-ae31-0b07a23fbfc1
md"
Determine the index number of the relevant parameter that needs to be modified in the model:
"

# ╔═╡ 73554f2f-f40f-4e23-b4bd-cef0c538f3d1
# ...
parameters(fermenter_monod)

# ╔═╡ 96a811c3-5b7e-4e26-9bbc-1470f71d755c
md"
Create a function called `affect2!`, that will be called by the solver at the timepoint(s) stored in `condition2` in order to alter the relevant parameter value:
"

# ╔═╡ a90cc966-cc6b-422c-9557-d2c8bfe9c37c
# function affect2!(integrator)
#     ...
# end
function affect2!(integrator)
    integrator.ps[:Sin]=2.8
end

# ╔═╡ 407eadc6-410a-48f2-ac61-c5e4d8ad7789
md"
Create the callback function using `condition2` and `affect2!`. Store it in `cb2`:
"

# ╔═╡ cef1051e-e699-43aa-85eb-1cbbd324f4dd
cb2 = PresetTimeCallback(condition2, affect2!)

# ╔═╡ 23d353d3-90ca-4148-a398-22dc20a52136
md"
Solve the ODE problem. Use `Tsit5()` and `saveat=0.5`. Store the solution in `osol2`:
"

# ╔═╡ 2cf06170-902f-44af-a1db-ea6567118c73
osol2 = solve(deepcopy(oprob), Tsit5(), saveat=0.5, callback=cb2)

# ╔═╡ 310a78a5-94ce-4a29-b7a1-37831ce5c64e
md"
Plot the results:
"

# ╔═╡ 85742dc1-24cb-42d1-a70b-70f1be6b6c1e
plot(osol2)

# ╔═╡ 196edf3a-b220-4a54-8137-b136b509617e
osol2[:S][end], osol2[:X][end]

# ╔═╡ 56eaa343-03b6-4cae-868b-e5c36ac66546
u_guess = [:S => osol2[:S][end], :X => osol2[:X][end]]

# ╔═╡ f121efc5-4e64-4e82-8672-2765ad85443e
params_mod = [:Y => 0.80, :Q => 2.0, :V => 40, :Sin => 2.8, :μmax => 0.30, :Ks => 0.15]

# ╔═╡ 0b992750-a446-447b-b2a1-26658c11c0bf
Seq, Xeq = solve(SteadyStateProblem(ODEProblem(fermenter_monod, u_guess, tspan, params_mod, combinatoric_ratelaws=false)))

# ╔═╡ 47a63fd8-f805-4c7d-8695-c9ee6550f24f
Seq, Xeq

# ╔═╡ 4d962d0f-da41-405e-9438-733d5668cde3
md"
Interpret the results. Ask yourself the following questions:

1. Can you clearly see the effect of the increase in $S_{in}$?
- Answer: ...
2. Find the steady state values of $S$ and $X$. Is the steady state value of $S$ influenced by the increase of $S_{in}$? Show how can you deduce that from the differential equations.
- Answer: ...
3. Can you explain why $X$ increased permanently?
- Answer: ...
"
#=
1. Yes, there is a temporary slight increase in substrate concentration $S$, while there is a permanent change in biomass concentration $X$.
2. The steady state value of S is seemingly not influenced by Sin. The steady state value of S (Seq) can be deduced from the second differential equation:
Seq = (Q/V)*Ks/(Y*mumax - Q/V)
and is independend of Sin.
3. The increase of X is due to the higher substrate concentration in the inlet flow (Sin). If there is more substrate S available, more biomass X will be produced.
=#

# ╔═╡ 53980767-a84f-44f2-a878-2a7d57e0e2ae
md"
#### Part 3

Suppose that the inlet/outlet flow $Q$ is suddenly doubled at $t = 100\;h$. Simulate the evolution of $S$ and $X$.
"

# ╔═╡ 31b64d91-f8ff-413a-9c0c-402fe2215a81
md"
Create the *condition* that contains the timepoint for the sudden change in $Q$. Store it in `condition3`:
"

# ╔═╡ 1e6a4e2f-044a-4797-b8ab-6a3756a0fea0
# condition3 = ...           # Uncomment and complete the instruction
condition3 = [100.0]

# ╔═╡ 772bcf35-0308-4106-9268-d844850db3a3
md"
Determine the order of the relevant parameter that needs to be modified in the model:
"

# ╔═╡ c5e0c018-c586-445c-bb15-083e2e059cb3
# ...                       # Uncomment and complete the instruction
parameters(fermenter_monod)

# ╔═╡ 0e8fb79d-4655-4a57-9b76-f151f4098fdb
md"
Create a function called `affect3!`, that will be called by the solver at the timepoint(s) stored in `condition3` in order to alter the relevant parameter value:
"

# ╔═╡ 589bdea5-673a-46fd-861b-f4fef40aea62
# Uncomment and complete the instruction
# function affect3!(integrator)
#     ...
# end
function affect3!(integrator)
    integrator.ps[:Q] *= 2.0
end

# ╔═╡ a6024052-df1b-4bff-abb4-b9cda5a4c08e
md"
Create the callback function using `condition3` and `affect3!`. Store it in `cb3`:
"

# ╔═╡ f688a7a0-1b2b-4e60-9d1a-444f09576176
# cb3 = ...           # Uncomment and complete the instruction
cb3 = PresetTimeCallback(condition3, affect3!)

# ╔═╡ 1fcb6283-eefa-4aa4-a7a5-1e90e8d6d814
md"
Solve the ODE problem. Use `Tsit5()` and `saveat=0.5`. Store the solution in `osol3`:
"

# ╔═╡ 03eadcc9-29f4-4141-bf20-94a5d62f8655
# osol3 = ...         # Uncomment and complete the instruction
osol3 = solve(deepcopy(oprob), Tsit5(), saveat=0.5, callback=cb3)

# ╔═╡ 39ef225d-3222-4998-bfd4-5ff88f74a0f9
md"
Plot the results:
"

# ╔═╡ 76849cf5-170b-452b-acdd-c4017feaad18
# ...
plot(osol3)

# ╔═╡ 040f040f-aa28-442e-9f44-2a897e22ed4f
md"
Interpret the results. Ask yourself the following questions:

1. Can you clearly see the effect of doubling of $Q$?
- Answer: ...
2. Can you argue, by means of reasoning, why $S$ increases and $X$ decreases?
- Answer: ...
"
#=
1. Yes, at t=100, S slightly increases and X slightly decreases
2. If the flow is increased, biomass X will leave the reactor vessel faster. Therefore the biomass concentration will decrease. Furthermore, this means that there will be less biomass X available to consume substrate S. Hence, the substrate concentration will increase.
=#

# ╔═╡ Cell order:
# ╠═2e58f4ae-f711-11ee-2598-7f3a6f2e2013
# ╠═539c1823-16d0-4355-97b6-fe9f0b106864
# ╠═e99680dc-73af-40aa-bf57-a06d3a7372be
# ╠═7856d878-8586-4cfd-9cf6-d61234450e41
# ╠═8500c35e-6bc3-4900-81bf-7705ddd61532
# ╠═e66518ee-b6f6-4cca-a224-30e01cffddbe
# ╠═f1350528-07a5-4860-ad2d-627588186abc
# ╠═331a34f4-89d4-4193-896c-c14ab0bf04e7
# ╠═42190228-40d3-48e8-b52f-156a0c7cbddc
# ╠═55746566-2d46-4475-851a-02b7fad87a1a
# ╠═ec9cb3bd-f5ed-4ab0-9b3d-b875692227ac
# ╠═bd648109-f042-42de-9e0e-017b502fab95
# ╠═67117a27-dcea-4b43-b962-9ad9fd07f4f4
# ╠═d13e6e38-037e-4812-85e9-2c18bed360f6
# ╠═4b556cf0-8fad-434d-be56-dc1848d898ae
# ╠═ea55d648-7575-43c3-a385-5f4979996ef2
# ╠═1365c12e-e662-4858-983b-02ba94cd9f0d
# ╠═3941bd60-a83c-4f72-84b3-28e28cb845d0
# ╠═d6c1316a-cf96-43d1-854a-f25925cf4a55
# ╠═4926b941-c3b6-4804-b4a4-11e13e5186f2
# ╠═ab2a9842-6a9c-46bd-812b-db01629d6a1c
# ╠═b6a526bd-6ee5-442b-9fb8-3fbe1e280dd4
# ╠═1f62e66d-571f-41ca-9f02-f36a8ca10ab9
# ╠═37ced6e3-b435-4546-a720-a1ec1af23a65
# ╠═d609bed4-94cf-4167-80fe-924501a5835c
# ╠═bb2d06f8-1940-4585-961a-54068da50e91
# ╠═9bb450c6-5499-42f6-8356-bdc4985b74e7
# ╠═0298953a-b90c-41cd-8613-cb47ce752e43
# ╠═c44e0ae9-addb-4518-a3ee-559a24852c3a
# ╠═36a3b978-57ff-4d03-ae31-0b07a23fbfc1
# ╠═73554f2f-f40f-4e23-b4bd-cef0c538f3d1
# ╠═96a811c3-5b7e-4e26-9bbc-1470f71d755c
# ╠═a90cc966-cc6b-422c-9557-d2c8bfe9c37c
# ╠═407eadc6-410a-48f2-ac61-c5e4d8ad7789
# ╠═cef1051e-e699-43aa-85eb-1cbbd324f4dd
# ╠═23d353d3-90ca-4148-a398-22dc20a52136
# ╠═2cf06170-902f-44af-a1db-ea6567118c73
# ╠═310a78a5-94ce-4a29-b7a1-37831ce5c64e
# ╠═85742dc1-24cb-42d1-a70b-70f1be6b6c1e
# ╠═196edf3a-b220-4a54-8137-b136b509617e
# ╠═56eaa343-03b6-4cae-868b-e5c36ac66546
# ╠═f121efc5-4e64-4e82-8672-2765ad85443e
# ╠═0b992750-a446-447b-b2a1-26658c11c0bf
# ╠═47a63fd8-f805-4c7d-8695-c9ee6550f24f
# ╠═4d962d0f-da41-405e-9438-733d5668cde3
# ╠═53980767-a84f-44f2-a878-2a7d57e0e2ae
# ╠═31b64d91-f8ff-413a-9c0c-402fe2215a81
# ╠═1e6a4e2f-044a-4797-b8ab-6a3756a0fea0
# ╠═772bcf35-0308-4106-9268-d844850db3a3
# ╠═c5e0c018-c586-445c-bb15-083e2e059cb3
# ╠═0e8fb79d-4655-4a57-9b76-f151f4098fdb
# ╠═589bdea5-673a-46fd-861b-f4fef40aea62
# ╠═a6024052-df1b-4bff-abb4-b9cda5a4c08e
# ╠═f688a7a0-1b2b-4e60-9d1a-444f09576176
# ╠═1fcb6283-eefa-4aa4-a7a5-1e90e8d6d814
# ╠═03eadcc9-29f4-4141-bf20-94a5d62f8655
# ╠═39ef225d-3222-4998-bfd4-5ff88f74a0f9
# ╠═76849cf5-170b-452b-acdd-c4017feaad18
# ╠═040f040f-aa28-442e-9f44-2a897e22ed4f
