### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ e99680dc-73af-40aa-bf57-a06d3a7372be
using Pkg; Pkg.activate("..")

# ╔═╡ 2e58f4ae-f711-11ee-2598-7f3a6f2e2013
using Markdown, InteractiveUtils

# ╔═╡ 2eddb05c-680b-4b4c-b969-df2688591762
using StatsPlots, PlutoUI; TableOfContents()

# ╔═╡ e66518ee-b6f6-4cca-a224-30e01cffddbe
using OrdinaryDiffEq, Catalyst

# ╔═╡ 7856d878-8586-4cfd-9cf6-d61234450e41
md"""
# Exercise: Fermenter - First order kinetics
"""

# ╔═╡ 8500c35e-6bc3-4900-81bf-7705ddd61532
md"""
In a fermenter reactor biomass $X$ grows on substrate $S$. The reactor is fed with a inlet flow rate $Q_{in}$ [$L/h$], which consist of a (manipulable) input concentration of substrate $S_{in}$ [$g/L$]. Inside the reactor, biomass, with a concentration of $X$ [$g/L$], is produced through first-order kinetics (first-order in $S$):

$$\begin{eqnarray*}
S \xrightarrow[\quad\quad]{\beta} Y \, X
\end{eqnarray*}$$

with $\beta$ [$h^{-1}$] the reaction rate constant, and $Y$ [$gX/gS$] the yield coefficient which is defined here by the amount of produced biomass by consumption of one unit of substrate. The reaction rate is actually:

$$\begin{eqnarray*}
r = \beta \, S
\end{eqnarray*}$$

Futhermore, the reactor is drained with an outlet flow $Q$ [$L/h$], which consist of the current concentrations of substrate $S$ [$g/L$] and biomass $X$ [$g/L$] inside the reactor. The volume $V$ [$L$] of the reactor content is kept constant by setting $Q_{in} = Q$.
"""

# ╔═╡ f1350528-07a5-4860-ad2d-627588186abc
md"""
Create a *reaction network object* model for the aforementioned problem in order to simulate the evolution of substrate $S$ and biomass $X$ with time. Name it `fermenter_firstorder`.
"""


# ╔═╡ 331a34f4-89d4-4193-896c-c14ab0bf04e7
fermenter_firstorder = @reaction_network begin
    β, S --> Y*X        # Y*X is created from one S at a rate β
	# β*S, S => Y*X        # Y*X is created from one S at a rate β
    Q/V*Sin, 0 --> S    # S is created at a rate Q/V*Sin(t)
    Q/V, (S, X) --> 0   # S and X are degraded at a rate Q/V*S
end

# ╔═╡ 55746566-2d46-4475-851a-02b7fad87a1a
md"""
Convert the system to a symbolic differential equation model and verify, by analyzing the differential equation, that your model is correctly implemented.
"""

# ╔═╡ ec9cb3bd-f5ed-4ab0-9b3d-b875692227ac
osys = convert(ODESystem, fermenter_firstorder, combinatoric_ratelaws=false)

# ╔═╡ 67117a27-dcea-4b43-b962-9ad9fd07f4f4
md"""
The parameter values are $\beta = 0.98$, $Y = 0.80$, $Q = 2.0$, $V = 40.0$ and $S_{in} = 2.2\;g/L$. With the latter values, the fermenter reactor is in steady state operation with concentrations for substrate $S = 0.1068\;g/L$ and biomass $X = 1.6746\;g/L$. Suppose that at timepoint $t = 20\;h$, the concentration of substrate in the inlet flow (cf. $S_{in}$) is suddently increased to $3.4\;g/L$. Simulate the evolution of $S$ and $X$ during 120 hours.
"""


# ╔═╡ d13e6e38-037e-4812-85e9-2c18bed360f6
md"""
Initialize a vector `u0` with the initial conditions:
"""

# ╔═╡ 4b556cf0-8fad-434d-be56-dc1848d898ae
u0 = [:S => 0.1068, :X => 1.6746]

# ╔═╡ ea55d648-7575-43c3-a385-5f4979996ef2
md"""
Set the timespan for the simulation:
"""

# ╔═╡ 1365c12e-e662-4858-983b-02ba94cd9f0d
tspan = (0.0, 120.0)

# ╔═╡ 3941bd60-a83c-4f72-84b3-28e28cb845d0
md"""
Initialize a vector `parms` with the parameter values:
"""

# ╔═╡ d6c1316a-cf96-43d1-854a-f25925cf4a55
parms = [:β => 0.98, :Y => 0.80, :Q => 2, :V => 40, :Sin => 2.2]

# ╔═╡ eeb8ec6e-154e-4fe9-8b5b-edbe71914985
md"""
Create the *condition* that contains the timepoint for the sudden change in $S_{in}$. Store it in `condition`:
"""

# ╔═╡ 7ca8efaa-97b6-46f2-b4d3-6ca8aa97dda7
condition = [20] => [fermenter_firstorder.Sin ~ 3.4]

# ╔═╡ d864bfc3-05b2-483b-9a55-da026931703f
md"""
Make a new *reaction system* where the discrete event is included. Name it `fermenter_firstorder_c`.
"""

# ╔═╡ d758b918-d184-473f-9b9e-ed6da7b0f088
@named fermenter_firstorder_c = ReactionSystem(equations(fermenter_firstorder), discrete_events=condition)

# ╔═╡ ef596d8c-efdc-4b5c-9584-335d799acfe8
md"""
Complete the new *reaction system*. Name it `fermenter_firstorder_c_com`.
"""

# ╔═╡ c07b3121-4c27-454c-b9e4-7dfa27371ebd
fermenter_firstorder_c_com = complete(fermenter_firstorder_c)

# ╔═╡ c6e81f41-a244-48c9-9d18-b3b9e0984fbb
md"""
Create the ODE problem and store it in `oprob`:
"""

# ╔═╡ ed56f8d6-2260-4829-9190-69b60b7d7599
oprob = ODEProblem(fermenter_firstorder_c_com, u0, tspan, parms, combinatoric_ratelaws=false);

# ╔═╡ 8b73b16b-7f7d-4d2e-a1c2-7e1adf2336e9
md"""
Solve the ODE problem. Make a deepcopy and use `Tsit5()` and `saveat=0.5`. Store the solution in `osol`:
"""

# ╔═╡ 433db8d1-f038-4e7b-9133-90bfeccabd07
osol = solve(deepcopy(oprob), Tsit5(), saveat=0.5)

# ╔═╡ fc824241-f718-46a8-b50b-a680470c062b
md"""
Plot the results:
"""

# ╔═╡ e2609226-5f2f-4802-89e3-0efeac740081
plot(osol)

# ╔═╡ c400711b-e531-4bfb-ba11-49184b868349
md"""
!!! questions
	Interpret the results. Ask yourself the following questions:
	1. Can you clearly see the effect of the increase in $S_{in}$?
	2. Can you argue, by means of reasoning or by determining and analyzing the operating point, why the increase of $X$ is larger than the increase of $S$?
"""

# ╔═╡ 975728cf-e064-4e59-92b5-1e83a14c99a1
md"""
Answers:
1. missing
2. missing
"""
#=
1. Yes, both S and X increase at t=20. One can see that for this set of parameters, the increase of X is larger than the increase of S.
2. Denote the operating point for S and X as Sw and Xw (w from working). These would be: (Sw, Xw) = ( Q/V*Sin/(Q/V+β), Y*β*Sin/(Q/V+β) ). In our case Y*β >> Q/V, so that any change in Sin, inflict a larger change in X than in S. On the otherhand, since if S is increased in the inlet flow, then S is increased also in the outlet flow, so S would not significantly change. If S is increased, then more X will be produced and the smaller the ratio Q/V, the more X will be produced.
=#

# ╔═╡ Cell order:
# ╠═2e58f4ae-f711-11ee-2598-7f3a6f2e2013
# ╠═e99680dc-73af-40aa-bf57-a06d3a7372be
# ╠═2eddb05c-680b-4b4c-b969-df2688591762
# ╠═e66518ee-b6f6-4cca-a224-30e01cffddbe
# ╟─7856d878-8586-4cfd-9cf6-d61234450e41
# ╟─8500c35e-6bc3-4900-81bf-7705ddd61532
# ╟─f1350528-07a5-4860-ad2d-627588186abc
# ╠═331a34f4-89d4-4193-896c-c14ab0bf04e7
# ╟─55746566-2d46-4475-851a-02b7fad87a1a
# ╠═ec9cb3bd-f5ed-4ab0-9b3d-b875692227ac
# ╟─67117a27-dcea-4b43-b962-9ad9fd07f4f4
# ╟─d13e6e38-037e-4812-85e9-2c18bed360f6
# ╠═4b556cf0-8fad-434d-be56-dc1848d898ae
# ╟─ea55d648-7575-43c3-a385-5f4979996ef2
# ╠═1365c12e-e662-4858-983b-02ba94cd9f0d
# ╟─3941bd60-a83c-4f72-84b3-28e28cb845d0
# ╠═d6c1316a-cf96-43d1-854a-f25925cf4a55
# ╟─eeb8ec6e-154e-4fe9-8b5b-edbe71914985
# ╠═7ca8efaa-97b6-46f2-b4d3-6ca8aa97dda7
# ╟─d864bfc3-05b2-483b-9a55-da026931703f
# ╠═d758b918-d184-473f-9b9e-ed6da7b0f088
# ╟─ef596d8c-efdc-4b5c-9584-335d799acfe8
# ╠═c07b3121-4c27-454c-b9e4-7dfa27371ebd
# ╟─c6e81f41-a244-48c9-9d18-b3b9e0984fbb
# ╠═ed56f8d6-2260-4829-9190-69b60b7d7599
# ╟─8b73b16b-7f7d-4d2e-a1c2-7e1adf2336e9
# ╠═433db8d1-f038-4e7b-9133-90bfeccabd07
# ╟─fc824241-f718-46a8-b50b-a680470c062b
# ╠═e2609226-5f2f-4802-89e3-0efeac740081
# ╠═c400711b-e531-4bfb-ba11-49184b868349
# ╠═975728cf-e064-4e59-92b5-1e83a14c99a1
