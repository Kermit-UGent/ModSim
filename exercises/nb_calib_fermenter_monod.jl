### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ c54dae10-60af-4141-b56d-ed61cb0ced8a
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate("..")
end
  ╠═╡ =#

# ╔═╡ 245ca9d0-10f9-11ef-0ef6-a73594e96db9
using Markdown

# ╔═╡ 78a25bef-31e5-45ef-b0ba-b9a8c8a9edeb
using InteractiveUtils

# ╔═╡ 16438e07-1b2b-467e-822a-081d19cae92b
using Catalyst, DifferentialEquations, Plots

# ╔═╡ 31243ea7-1f0f-490f-8886-b1b7ab7ae5b4
using Optim

# ╔═╡ 2f0a4c62-3441-4c63-9bb9-383e7f554eb5
md"
### Exercise: Fermenter - Monod kinetics - Calibration
"

# ╔═╡ 595ea8ee-bc67-4696-9232-982612fb554d
md"
In one of the previous practica we were introduced to a fermenter in which biomass $X$ [$g/L$] grows by breaking down substrate $S$ [$g/L$]. The reactor is fed with a inlet flow rate $Q_{in}$ [$L/h$], which consist of a (manipulable) input concentration of substrate $S_{in}$ [$g/L$]. This process was modelled using Monod kinetics, resulting in the model below:

$$\begin{eqnarray*}
%S \xrightarrow[\quad\quad]{\beta} Y \, X
S \xrightarrow[\quad\quad]{r} Y \, X \quad\quad\quad\quad r = \mu \, X \quad \textrm{with} \quad \mu = \mu_{max} \, \cfrac{S}{S + K_s}
\end{eqnarray*}$$

Note that $Y$ is the yield and not a species.
"

# ╔═╡ 824db995-7a66-4719-a534-7e0f6dec90b5
md"
The *reaction network object* for this model could be set-up as:
"

# ╔═╡ 245c2636-95da-4c76-8b03-c4d20bbabb48
fermenter_monod = @reaction_network begin
    X * mm(S, μmax, Ks), S => Y*X
    Q/V, (S, X) --> ∅
    Q/V*Sin, ∅ --> S
end

# ╔═╡ de8ddc14-8f82-403d-8f42-29673ef2a722
md"
which resulted in the following differential equations:

$$\begin{eqnarray*}
\cfrac{dS}{dt} &=& \cfrac{Q}{V} \left(S_{in} - S \right) - \mu_{max}\cfrac{S}{S + K_s} X\\
\cfrac{dX}{dt} &=& -\cfrac{Q}{V} X + Y \mu_{max}\cfrac{S}{S + K_s} X
\end{eqnarray*}$$
"

# ╔═╡ 94a8834d-9528-41a2-b9d1-54f87833adca
convert(ODESystem, fermenter_monod, combinatoric_ratelaws=false)

# ╔═╡ b7b7d58f-d406-4596-b834-ced6d8fada83
md"
Suppose that during an experiment measurement data has been collected of the substrate $S$ and biomass $X$ concentration at an interval of $2\;h$ within $100\;h$:
"

# ╔═╡ 99c6f31a-0968-4804-9980-71fcc1af1f49
S_meas = [0.002, 0.204, 0.423, 0.581, 0.676, 0.873, 0.95, 1.131, 1.228, 1.243, 1.341, 1.386, 1.439, 1.443, 1.451, 1.481, 1.227, 1.119, 1.007, 0.823, 0.801, 0.64, 0.576, 0.514, 0.485, 0.467, 0.458, 0.424, 0.415, 0.438, 0.463, 0.375, 0.451, 0.413, 0.42, 0.408, 0.382, 0.463, 0.4, 0.367, 0.375, 0.416, 0.381, 0.408, 0.382, 0.366, 0.378, 0.416, 0.375, 0.386, 0.393]

# ╔═╡ bf4ad873-e0fe-415c-9e78-fe0b5ac1414e
X_meas = [0.035, 0.014, 0.007, 0.01, 0.004, 0.039, 0.074, 0.004, 0.0, 0.013, 0.043, 0.06, 0.078, 0.1, 0.188, 0.184, 0.343, 0.506, 0.551, 0.622, 0.813, 0.967, 0.832, 1.051, 1.089, 1.177, 1.228, 1.166, 1.227, 1.111, 1.044, 1.146, 1.251, 1.27, 1.29, 1.192, 1.359, 1.283, 1.097, 1.216, 1.23, 1.25, 1.227, 1.07, 1.206, 1.404, 1.311, 1.378, 1.263, 1.331, 1.241]

# ╔═╡ 1dae5875-f405-4ecb-8b7b-3c3f22b549bb
t_meas = 0.0:2.0:100.0

# ╔═╡ f6f576fa-c8d7-4a12-af09-e1bc15c24de4
md"
In addition, assume the follwoing measurement errors:

$$\sigma_S = 0.05\;g/L \quad\quad\quad \sigma_X = 0.10\;g/L$$

We will call them `S_sigma` and `X_sigma`:
"

# ╔═╡ 397588a2-905b-4d43-bb42-24590dd285b8
S_sigma = 0.05

# ╔═╡ f23a7cee-b40b-4863-8792-a77fa5976730
X_sigma = 0.10

# ╔═╡ 6c481447-28c6-4530-bf2c-64762121bc71
md"
We can make a scatter plot of the measured data for both $S$ and $X$ in the following way:
"

# ╔═╡ 918fd524-81fa-4aff-a403-37402e47235b
begin
	scatter(t_meas, S_meas, label="S meas", color=:blue)
    scatter!(t_meas, X_meas, label="X meas", color=:red) # notice exclamation mark !
end

# ╔═╡ 27811996-7689-421d-adb4-ff379893c834
md"
Remark the following:

- In order to execute multiple command in one notebook cell, we need to place these commands in a `begin`-`end`-block.
- In order to plot the data for $X$ in the same plot as the data for $S$, we need to place an *exclamation* mark `!` after the second `scatter`-command.
"

# ╔═╡ ef977370-06ee-4a73-85e2-609a744167d3
md"
We have previously used the following parameter values:

-  $\mu_{max} = 0.30\;h^{-1}$, $K_s = 0.15\;g/L$, $S_{in} = 2.2\;g/L$
-  $Y = 0.80$, $Q = 2.0\;L/h$, $V = 40.0\;L$

Furthermore, suppose that at $t = 0\;h$ no substrate $S$ is present in the reactor but that there is initially some biomass with a concetration of $0.01\;g/L$.

Calibrate the parameter values for $\mu_{max}$, $K_s$ and $S_{in}$ using the aforementioned measurement data for $S$ and $X$ in a timespan of $[0, 100]\,h$. Take the values above as initial values for $\mu_{max}$, $K_s$ and $S_{in}$.
"

# ╔═╡ 7a227eaf-18d0-44f4-ac4b-f529e81c7471
md"
Initialize a vector `u₀` with the initial conditions, and set the timespan:
"

# ╔═╡ 6375478f-1af9-4fd2-b6f3-101a6f796f2d
# u₀ = missing      # Uncomment and complete the instruction
u₀ = [:S => 0.0, :X => 0.01]

# ╔═╡ 38fe8304-af61-40a7-ac86-480dfb892185
# tspan = missing
tspan = (0.0, 100)  # Uncomment and complete the instruction

# ╔═╡ f6a8f134-6db0-4d74-8af5-82826347d8f0
md"
Set up the objective function and name it `Jtheta_fermenter_monod`:
"

# ╔═╡ 4c28a66a-ee2c-42a2-95c7-ea4ddb6a232d
# Uncomment and complete the instruction
# function Jtheta_fermenter_monod(thetas)
#     params = missing
#     oprob = missing
#     sol = missing
#     S_sol = missing
#     X_sol = missing
#     J = missing
#     return J
# end
function Jtheta_fermenter_monod(thetas)           # function return to be minimized
    params = [:μmax => thetas[1], :Ks => thetas[2], :Y => thetas[3],
		      :Q => 2, :V => 40, :Sin => 2.2]
    oprob = ODEProblem(fermenter_monod, u₀, tspan, params,
		               combinatoric_ratelaws=false)
    sol = solve(oprob, Tsit5(), saveat=t_meas)
    S_sol = sol[:S]
    X_sol = sol[:X]
    J = (1/S_sigma^2)*sum(abs2, S_sol - S_meas) + 
	    (1/X_sigma^2)*sum(abs2, X_sol - X_meas)
    return J
end

# ╔═╡ 3136b15d-5078-4bcd-954b-e89bcb8aed1b
md"
Create a vector `thetas_init` with the initial values for $\mu_{max}$, $K_s$ and $S_{in}$:
"

# ╔═╡ 6a508a62-61b9-4273-8e45-b26f594e8da9
# thetas_init = missing       # Uncomment and complete the instruction
thetas_init = [0.30, 0.15, 0.80]

# ╔═╡ 63420055-55f8-4def-8b0e-11ea61483010
md"
Minimize the value of the objective function by optimizing the parameter values. Store the optimization results in `results`.
"

# ╔═╡ d52c9da8-d8a4-4db0-ac6d-6d16ccf4775c
# results = missing           # Uncomment and complete the instruction
results = optimize(Jtheta_fermenter_monod, thetas_init, NelderMead())

# ╔═╡ e1b0ee01-f16c-40e9-a0f9-80072d690936
md"
Store the optimized values in the variable `thetas_opt`:
"

# ╔═╡ f2d7daf8-8218-446d-b1d2-e9e05aeadfd9
# thetas_opt = missing        # Uncomment and complete the instruction
thetas_opt = Optim.minimizer(results)

# ╔═╡ 23d58bb1-d077-402e-8bee-3866c68e069a
md"
Determine the minimum value of the objective function:
"

# ╔═╡ 7b3a3677-b251-43c1-b125-6d6ff1a11ea3
# missing                     # Uncomment and complete the instruction
Optim.minimum(results)

# ╔═╡ 05d13a48-adc8-4e24-a6e4-be24af2c7a59
md"
Set up parameter values with optimized parameter values:
"

# ╔═╡ a9be9ad2-8933-4e8d-a110-e93564330578
# params = missing            # Uncomment and complete the instruction
params = [:μmax => thetas_opt[1], :Ks => thetas_opt[2], :Y => thetas_opt[3], :Q => 2, :V => 40, :Sin => 2.2]

# ╔═╡ 4e8870dc-2da6-4b80-82d6-26c7ceedad7d
md"
Create an ODEProblem and solve it:
"

# ╔═╡ fef343a8-da6f-47cc-bb07-a629b5e9566e
# oprob = missing             # Uncomment and complete the instruction
oprob = ODEProblem(fermenter_monod, u₀, tspan, params, combinatoric_ratelaws=false)

# ╔═╡ 428fdb55-4d18-4937-aa0e-90e0294cd6d8
# osol = missing               # Uncomment and complete the instruction
osol = solve(oprob, Tsit5(), saveat=0.5)

# ╔═╡ 5a39b0e0-1ea1-4854-8e68-66d0d4bbf25c
md"
Plot $S$ and $X$ simulated with the optimized parameter values together with the measured data that was used to find the optimized values.
"

# ╔═╡ 28b89b26-e5ee-40b0-8b45-a925c104c580
# Uncomment and complete the instruction
# begin
#     missing
#     missing
#     missing
# end
begin
    plot(osol, labels=["S sim" "X sim"], xlabel="t (h)")
    scatter!(t_meas, S_meas, label="S meas", color=:blue)
    scatter!(t_meas, X_meas, label="X meas", color=:red)
end

# ╔═╡ Cell order:
# ╠═245ca9d0-10f9-11ef-0ef6-a73594e96db9
# ╠═78a25bef-31e5-45ef-b0ba-b9a8c8a9edeb
# ╠═c54dae10-60af-4141-b56d-ed61cb0ced8a
# ╠═16438e07-1b2b-467e-822a-081d19cae92b
# ╠═31243ea7-1f0f-490f-8886-b1b7ab7ae5b4
# ╠═2f0a4c62-3441-4c63-9bb9-383e7f554eb5
# ╠═595ea8ee-bc67-4696-9232-982612fb554d
# ╠═824db995-7a66-4719-a534-7e0f6dec90b5
# ╠═245c2636-95da-4c76-8b03-c4d20bbabb48
# ╠═de8ddc14-8f82-403d-8f42-29673ef2a722
# ╠═94a8834d-9528-41a2-b9d1-54f87833adca
# ╠═b7b7d58f-d406-4596-b834-ced6d8fada83
# ╠═99c6f31a-0968-4804-9980-71fcc1af1f49
# ╠═bf4ad873-e0fe-415c-9e78-fe0b5ac1414e
# ╠═1dae5875-f405-4ecb-8b7b-3c3f22b549bb
# ╠═f6f576fa-c8d7-4a12-af09-e1bc15c24de4
# ╠═397588a2-905b-4d43-bb42-24590dd285b8
# ╠═f23a7cee-b40b-4863-8792-a77fa5976730
# ╠═6c481447-28c6-4530-bf2c-64762121bc71
# ╠═918fd524-81fa-4aff-a403-37402e47235b
# ╠═27811996-7689-421d-adb4-ff379893c834
# ╠═ef977370-06ee-4a73-85e2-609a744167d3
# ╠═7a227eaf-18d0-44f4-ac4b-f529e81c7471
# ╠═6375478f-1af9-4fd2-b6f3-101a6f796f2d
# ╠═38fe8304-af61-40a7-ac86-480dfb892185
# ╠═f6a8f134-6db0-4d74-8af5-82826347d8f0
# ╠═4c28a66a-ee2c-42a2-95c7-ea4ddb6a232d
# ╠═3136b15d-5078-4bcd-954b-e89bcb8aed1b
# ╠═6a508a62-61b9-4273-8e45-b26f594e8da9
# ╠═63420055-55f8-4def-8b0e-11ea61483010
# ╠═d52c9da8-d8a4-4db0-ac6d-6d16ccf4775c
# ╠═e1b0ee01-f16c-40e9-a0f9-80072d690936
# ╠═f2d7daf8-8218-446d-b1d2-e9e05aeadfd9
# ╠═23d58bb1-d077-402e-8bee-3866c68e069a
# ╠═7b3a3677-b251-43c1-b125-6d6ff1a11ea3
# ╠═05d13a48-adc8-4e24-a6e4-be24af2c7a59
# ╠═a9be9ad2-8933-4e8d-a110-e93564330578
# ╠═4e8870dc-2da6-4b80-82d6-26c7ceedad7d
# ╠═fef343a8-da6f-47cc-bb07-a629b5e9566e
# ╠═428fdb55-4d18-4937-aa0e-90e0294cd6d8
# ╠═5a39b0e0-1ea1-4854-8e68-66d0d4bbf25c
# ╠═28b89b26-e5ee-40b0-8b45-a925c104c580
