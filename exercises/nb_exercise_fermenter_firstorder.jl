### A Pluto.jl notebook ###
# v0.19.38

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
### Exercise: Fermenter - First order kinetics
"

# ╔═╡ 8500c35e-6bc3-4900-81bf-7705ddd61532
md"
In a fermenter reactor biomass grows on substrate. The reactor is fed with a inlet flow rate $Q_{in}$ [$L/h$], which consist of a (manipulable) input concentration of substrate $S_{in}$ [$g/L$]. Inside the reactor, biomass, with a concentration of $X$ [$g/L$], is produced through first-order kinetics:

$$\begin{eqnarray*}
S \xrightarrow[\quad\quad]{\beta} Y \, X
%S \xrightarrow[\quad\quad]{\text{r}} Y \, X \quad\quad\quad\quad r = \beta \, S
\end{eqnarray*}$$

with $\beta$ [$h^{-1}$] the reaction rate constant, and $Y$ [$gX/gS$] the yield coefficient which is defined here by the amount of produced biomass by consumption of one unit of substrate. Futhermore, the reactor is drained with an outlet flow $Q$ [$L/h$], which consist of the current concentrations of substrate $S$ [$g/L$] and biomass $X$ [$g/L$] inside the reactor. The volume $V$ [$L$] of the reactor content is kept constant by setting $Q_{in} = Q$.
"

# ╔═╡ f1350528-07a5-4860-ad2d-627588186abc
md"
Create a *reaction network object* model for the aforementioned problem in order to simulate the evolution of substrate $S$ and biomass $X$ with time. Name it `fermenter_firstorder`.
"


# ╔═╡ 331a34f4-89d4-4193-896c-c14ab0bf04e7
# fermenter_firstorder = @reaction_network begin
#     ...        # Y*X is created from one S at a rate β
#     ...        # S is created at a rate Q/V*Sin
#     ...        # S and X are degraded at a rate Q/V*S
# end
fermenter_firstorder = @reaction_network begin
    β, S --> Y*X        # Y*X is created from one S at a rate β
    Q/V*Sin, 0 --> S    # S is created at a rate Q/V*Sin(t)
    Q/V, (S, X) --> 0   # S and X are degraded at a rate Q/V*S
end

# ╔═╡ 55746566-2d46-4475-851a-02b7fad87a1a
md"
Convert the system to a symbolic differential equation model and verify, by analyzing the differential equation, that your model is correctly implemented.
"

# ╔═╡ ec9cb3bd-f5ed-4ab0-9b3d-b875692227ac
# osys = ...
osys = convert(ODESystem, fermenter_firstorder)

# ╔═╡ 67117a27-dcea-4b43-b962-9ad9fd07f4f4
md"
The parameter values are $\beta = 0.98$, $Y = 0.80$, $Q = 2.0$, $V = 40.0$ and $S_{in} = 2.2\;g/L$. With the latter values, the fermenter reactor is in steady state operation with concentrations for substrate $S = 0.1068\;g/L$ and biomass $X = 1.6746\;g/L$. Suppose that at timepoint $t = 20\;h$, the concentration of substrate in the inlet flow (cf. $S_{in}$) is suddently increased to $3.4\;g/L$. Simulate the evolution of $S$ and $X$ during 120 hours.
"


# ╔═╡ d13e6e38-037e-4812-85e9-2c18bed360f6
md"
Initialize a vector `u₀` with the initial conditions:
"

# ╔═╡ 4b556cf0-8fad-434d-be56-dc1848d898ae
# u0 = ...            # Uncomment and complete the instruction
u0 = [:S => 0.1068, :X => 1.6746]

# ╔═╡ ea55d648-7575-43c3-a385-5f4979996ef2
md"
Set the timespan for the simulation:
"

# ╔═╡ 1365c12e-e662-4858-983b-02ba94cd9f0d
# tspan = ...         # Uncomment and complete the instruction
tspan = (0.0, 120.0)

# ╔═╡ 3941bd60-a83c-4f72-84b3-28e28cb845d0
md"
Initialize a vector `param` with the parameter values:
"

# ╔═╡ d6c1316a-cf96-43d1-854a-f25925cf4a55
# params = ...         # Uncomment and complete the instruction
params = [:β => 0.98, :Y => 0.80, :Q => 2, :V => 40, :Sin => 2.2]

# ╔═╡ 4926b941-c3b6-4804-b4a4-11e13e5186f2
md"
Create the ODE problem and store it in `oprob`:
"

# ╔═╡ ab2a9842-6a9c-46bd-812b-db01629d6a1c
# oprob = ...           # Uncomment and complete the instruction
oprob = ODEProblem(fermenter_firstorder, u0, tspan, params, combinatoric_ratelaws=false)

# ╔═╡ 31b64d91-f8ff-413a-9c0c-402fe2215a81
md"
Create the *condition* that contains the timepoint for the sudden change in $S_{in}$. Store it in `condition`:
"

# ╔═╡ 1e6a4e2f-044a-4797-b8ab-6a3756a0fea0
# condition = ...           # Uncomment and complete the instruction
condition = [20.0]

# ╔═╡ 772bcf35-0308-4106-9268-d844850db3a3
md"
Determine the index number of the relevant parameter that needs to be modified in the model:
"

# ╔═╡ c5e0c018-c586-445c-bb15-083e2e059cb3
# ...                       # Uncomment and complete the instruction
parameters(fermenter_firstorder)

# ╔═╡ 0e8fb79d-4655-4a57-9b76-f151f4098fdb
md"
Create a function called `affect!`, that will be called by the solver at the timepoint(s) stored in `condition` in order to alter the relevant parameter value:
"

# ╔═╡ 589bdea5-673a-46fd-861b-f4fef40aea62
# Uncomment and complete the instruction
# function affect!(integrator)
#     ...
# end
function affect!(integrator)
    integrator.p[5] = 3.4
end

# ╔═╡ a6024052-df1b-4bff-abb4-b9cda5a4c08e
md"
Create the callback function using `condition` and `affect!`. Store it in `cb`:
"

# ╔═╡ f688a7a0-1b2b-4e60-9d1a-444f09576176
# cb = ...           # Uncomment and complete the instruction
cb = PresetTimeCallback(condition, affect!)

# ╔═╡ 1fcb6283-eefa-4aa4-a7a5-1e90e8d6d814
md"
Solve the ODE problem. Use `Tsit5()` and `saveat=0.5`. Store the solution in `osol`:
"

# ╔═╡ 03eadcc9-29f4-4141-bf20-94a5d62f8655
# osol = ...         # Uncomment and complete the instruction
osol = solve(deepcopy(oprob), Tsit5(), saveat=0.5, callback=cb)

# ╔═╡ 39ef225d-3222-4998-bfd4-5ff88f74a0f9
md"
Plot the results:
"

# ╔═╡ 76849cf5-170b-452b-acdd-c4017feaad18
# ...
plot(osol)

# ╔═╡ 040f040f-aa28-442e-9f44-2a897e22ed4f
md"
Interpret the results. Ask yourself the following questions:

1. Can you clearly see the effect of the increase in $S_{in}$?
- Answer: ...
2. Can you argue, by means of reasoning or by determining and analyzing the operating point, why the increase of $X$ is larger than the increase of $S$?
- Answer: ...
"


# ╔═╡ Cell order:
# ╠═2e58f4ae-f711-11ee-2598-7f3a6f2e2013
# ╠═539c1823-16d0-4355-97b6-fe9f0b106864
# ╠═e99680dc-73af-40aa-bf57-a06d3a7372be
# ╠═7856d878-8586-4cfd-9cf6-d61234450e41
# ╠═8500c35e-6bc3-4900-81bf-7705ddd61532
# ╠═e66518ee-b6f6-4cca-a224-30e01cffddbe
# ╠═f1350528-07a5-4860-ad2d-627588186abc
# ╠═331a34f4-89d4-4193-896c-c14ab0bf04e7
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
