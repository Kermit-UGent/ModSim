### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 715c1fae-2a4d-11f0-0569-bd3495497826
using Pkg; Pkg.activate("..")

# ╔═╡ 34856dc0-e961-4e3e-a996-1966e4191e35
using Plots, PlutoUI; TableOfContents()

# ╔═╡ bc7f57f2-bb4d-44d4-b693-f4cf412b470a
using DifferentialEquations, ModelingToolkit

# ╔═╡ 45c84728-ff11-4345-9c71-fab084e3159f
using ModelingToolkit: t_nounits as t, D_nounits as D

# ╔═╡ 4252c39a-ab0f-4392-8e77-1c2e2c3e36e6
md"""
# Exercise: Irrigation experiment
"""

# ╔═╡ 22a32fa7-d9e3-4ed0-bc69-0f82b229a524
md"""
An irrigation experiment is carried out on a soil column consisting of two layers of soil, each with specific soil characteristics. An adjustable volume of water per unit of time, $r$, is irrigated evenly over the soil column, starting with $5\;mm\,h^{-1}$ ($mm$ indicates a volume of water: $1\;mm = 10^{-3}\,m^3$). After $60\;h$ the added flow rate is increased to $10\;mm\,h^{-1}$. The water falls on the upper layer and percolates to the lower layer. The **relative moisture content** in both layers (i.e., relative to their residual moisture contents) is denoted by $S_1$ and $S_2$. Initially a moisture content of $30\;mm$ is present in the upper layer (cf. $S_1$) and of $25\;mm$ in the lower layer (cf. $S_2$). The residual moisture content in the upper layer is $S_{1,res}=10 \;mm$.

A model description of the relative moisture content in both soil layers is given by:

$$\begin{align}
\frac{dS_1}{dt} &= r\left(1-\cfrac{S_{1,res}}{S_{max}}\right) - \cfrac{r}{S_{max}}S_1 - \cfrac{k}{S_{max}}S_1 \\
\frac{dS_2}{dt} &= \cfrac{k}{S_{max}}S_1 - v \,S_2^2
\end{align}$$

Here $S_{max}$ ($150\;mm$) is the saturated water quantity for the top soil layer, $k$ is the percolation ratio ($3\;mm\,h^{-1}$) and $v$ is the flow factor into the groundwater ($10^{-3}\;h^{-1}\,mm^{-1}$).

Three measurements are made over the duration ($150\;h$) of the experiment:
- The excess running water (runoff): $R = r \cfrac{S_1 + S_{1,res}}{S_{max}}$,
- The underground outflow into groundwater: $O = v\,S_2^2$,
- The amount of percolation to deeper soil layers: $P = \cfrac{k}{S_{max}} S_1$.

The latter three are called *observables*.
"""

# ╔═╡ 1bebe2ff-603a-42ba-8526-ab100fa3880d
md"""
Model the aforementioned system of differential equations using ModelingToolkit.
"""

# ╔═╡ f5389efd-941f-488b-8d39-407f223277ec
md"""
## Setting up the equations
"""

# ╔═╡ 832e65cb-5029-48f5-9ed2-41c46e759a6c
md"""
Define the variables $S_1$ and $S_2$ with their corresponding initial values. Consider the observables $R$, $O$ and $P$ also as variables, as well as the 'parameter' $r$. Mention the dependency on the time $t$.
"""

# ╔═╡ 51fbf7e9-00f6-4864-9844-eabe8fdc6f1f
# @variables missing
@variables S₁(t)=30 S₂(t)=25 R(t) O(t) P(t) r(t)

# ╔═╡ 3c396ae5-bbbf-40f4-9203-f7fdce4cfca8
md"""
Define the parameters for this model and assign their corresponding values. Mind that $r$ is a variable and, hence, should not be listed as a parameter.

| Parameter | Value | Unit |
|:---------- |:---------- |:------------|
| `k`    | 3 | $mm\,h^{-1}$  |
| `v`    | 1.0e-3 | $h^{-1}\,mm^{-1}$  |
| `S₁res`| 10 | $mm$  |
| `Smax` | 150 | $mm$  |
"""

# ╔═╡ 69b88c66-1b0c-4881-83cc-44b66b6b3935
# @parameters missing
@parameters k=3 v=1.0e-3 S₁res=10 Smax=150

# ╔═╡ ffffbbcb-29d2-4c70-8c68-97be874608e8
md"""
Set up de model equations.
"""

# ╔═╡ d7e2b0d6-cdbc-49cd-a2fa-d055c6476b5b
# change_S1 = missing
change_S1 = D(S₁) ~ r*(1 - S₁res/Smax) - k/Smax*S₁ - k/Smax*S₁

# ╔═╡ 64516f9d-3ce2-4fee-a5b7-3f89e5df992f
# change_S2 = missing
change_S2 = D(S₂) ~ k/Smax*S₁ - v*S₂^2

# ╔═╡ e326c749-2b10-4f15-b10b-49b4b7253c92
md"""
Set up an equation for the variable $r$.

!!! hint
	The equation contains two terms: a constant term $5$ and a another term that will add $5$ when $t > 60$.
"""

# ╔═╡ caa7a0bb-4d0e-450a-8057-8efb4fbfdf58
# eq_flow_rate = missing
eq_flow_rate = r ~ 5.0 + 5.0*(t > 60.0)

# ╔═╡ 3905ea2d-bd3c-4043-a802-0d68386011c4
md"""
Set up the equations for the runoff $R$, the outflow $O$ and the percolation $P$.
"""

# ╔═╡ 289135c1-b97b-46b3-87b2-a0a02538fa95
# eq_runoff = missing
eq_runoff = R ~ r*(S₁+S₁res)/Smax

# ╔═╡ 1523b4c5-0c02-4848-bfe3-68695f5a30b7
# eq_outflow = missing
eq_outflow = O ~ v*S₂^2

# ╔═╡ 817dad8d-99fe-415c-87dd-a52ad4b2bf9f
# eq_percolation = missing
eq_percolation = P ~ k*S₁/Smax

# ╔═╡ a97c2de8-1873-45f5-a708-3b34f9a216ae
md"""
Bundle all equations.
"""

# ╔═╡ 8236ec1e-5236-4c3f-a5ed-63e7aef63d01
# eqs_irrigation = missing
eqs_irrigation = [change_S1, change_S2, eq_flow_rate, eq_runoff, eq_outflow, eq_percolation]

# ╔═╡ e0456a9a-dbb3-4f9c-889c-7e8dfc5e082b
md"""
## Building the ODE system
"""

# ╔═╡ 872ca902-45b6-496e-ae96-fa1b1c6e16bd
md"""
Build a system of equations with `@mtkbuild`. Name it `sys_irrigation`.
"""

# ╔═╡ 9c70c339-c2ac-4983-8564-4a02fdb73e8a
# @mtkbuild missing
@mtkbuild sys_irrigation = ODESystem(eqs_irrigation, t)

# ╔═╡ 5840e484-447d-4ed7-bdf1-98c71141ba55
md"""
## Create and solve the ODE problem
"""

# ╔═╡ 72f30ddc-d442-49e9-8b18-ac94cb561a7f
md"""
Create the ODE problem and store it in `oprob_irrigation`:
"""

# ╔═╡ fa1c55a1-3020-44f0-88b4-3e089fd691ec
# oprob_irrigation = missing
oprob_irrigation = ODEProblem(sys_irrigation, [], (0.0, 150.0), []);

# ╔═╡ a8c3ae37-a2d0-438e-8acf-36ba19f981ac
md"""
Solve the ODE problem. Use `Tsit5()`, `saveat=1.0` and `reltol=1e-9`. Store the solution in `osol_irrigation`:
"""

# ╔═╡ 0e382f7c-cf9c-4ed2-b51a-32539a9eb8dc
# osol_irrigation = missing
osol_irrigation = solve(oprob_irrigation, Tsit5(), saveat=1.0, reltol=1e-9)

# ╔═╡ 9ec65f94-4577-4f54-9565-e7fd2061dfd0
md"""
## Plotting results
"""

# ╔═╡ cebac4b0-1ca1-4e5c-adb6-73fbe2b3ad1f
md"""
Plot $R$, $O$ and $P$. Use the option `idxs=[R, O, P]` for that.
"""

# ╔═╡ c9aa4a15-5fd4-4598-ad46-6680585965f7
# missing
plot(osol_irrigation; idxs=[R, O, P])

# ╔═╡ 8954cc45-2feb-40fe-a622-7ccd13661ad9
md"""
Check out the solution for $S_1$ and $S_2$.
"""

# ╔═╡ ff59a8f8-dc05-4b2e-be7f-3cef50ea5605
# missing
osol_irrigation[S₁]

# ╔═╡ 86c80ab7-a7d9-4224-817b-8b84de3885a9
# missing
osol_irrigation[S₂]

# ╔═╡ fd18fb8a-f0d5-417a-9854-8d8f0f0bb3ce
md"""
What is the value of $S_1$ for $t=100\;h$?
"""

# ╔═╡ 90586677-2b9a-4fef-b9ef-098eb288ee64
# missing
osol_irrigation[S₁][100]    # juist

# ╔═╡ 19d5382d-a657-4993-bc87-75e9edfd1dbd
md"""
Calculate the value of $R$ for $t=100\;h$ using the value for $S_1$ at $t=100\;h$ and the parameters.
"""

# ╔═╡ ccf44020-d61a-4bb0-addd-4b27a9dbc0ad
# missing
(osol_irrigation[S₁][100]+10.0)*10.0/150.0

# ╔═╡ 7cac3f5c-684f-4edc-b55f-6cf343556291
md"""
!!! question
	Does this value correspond to the value you can determine from the plot of $R$?
"""

# ╔═╡ 9f3cce18-4af3-4270-803b-b6b1c773bb4c
md"""
Answers: missing
"""

# ╔═╡ 0f774f27-6251-46d3-a9e1-64e3d20ea9be
md"""
Make a plot of $S_1$ and $S_2$.
"""

# ╔═╡ 2ac192f4-fa55-44bd-9b26-854a0f3266ca
# missing
plot(osol_irrigation)

# ╔═╡ Cell order:
# ╟─4252c39a-ab0f-4392-8e77-1c2e2c3e36e6
# ╠═715c1fae-2a4d-11f0-0569-bd3495497826
# ╠═34856dc0-e961-4e3e-a996-1966e4191e35
# ╠═bc7f57f2-bb4d-44d4-b693-f4cf412b470a
# ╠═45c84728-ff11-4345-9c71-fab084e3159f
# ╟─22a32fa7-d9e3-4ed0-bc69-0f82b229a524
# ╟─1bebe2ff-603a-42ba-8526-ab100fa3880d
# ╟─f5389efd-941f-488b-8d39-407f223277ec
# ╟─832e65cb-5029-48f5-9ed2-41c46e759a6c
# ╠═51fbf7e9-00f6-4864-9844-eabe8fdc6f1f
# ╟─3c396ae5-bbbf-40f4-9203-f7fdce4cfca8
# ╠═69b88c66-1b0c-4881-83cc-44b66b6b3935
# ╟─ffffbbcb-29d2-4c70-8c68-97be874608e8
# ╠═d7e2b0d6-cdbc-49cd-a2fa-d055c6476b5b
# ╠═64516f9d-3ce2-4fee-a5b7-3f89e5df992f
# ╟─e326c749-2b10-4f15-b10b-49b4b7253c92
# ╠═caa7a0bb-4d0e-450a-8057-8efb4fbfdf58
# ╟─3905ea2d-bd3c-4043-a802-0d68386011c4
# ╠═289135c1-b97b-46b3-87b2-a0a02538fa95
# ╠═1523b4c5-0c02-4848-bfe3-68695f5a30b7
# ╠═817dad8d-99fe-415c-87dd-a52ad4b2bf9f
# ╟─a97c2de8-1873-45f5-a708-3b34f9a216ae
# ╠═8236ec1e-5236-4c3f-a5ed-63e7aef63d01
# ╟─e0456a9a-dbb3-4f9c-889c-7e8dfc5e082b
# ╟─872ca902-45b6-496e-ae96-fa1b1c6e16bd
# ╠═9c70c339-c2ac-4983-8564-4a02fdb73e8a
# ╟─5840e484-447d-4ed7-bdf1-98c71141ba55
# ╟─72f30ddc-d442-49e9-8b18-ac94cb561a7f
# ╠═fa1c55a1-3020-44f0-88b4-3e089fd691ec
# ╟─a8c3ae37-a2d0-438e-8acf-36ba19f981ac
# ╠═0e382f7c-cf9c-4ed2-b51a-32539a9eb8dc
# ╟─9ec65f94-4577-4f54-9565-e7fd2061dfd0
# ╟─cebac4b0-1ca1-4e5c-adb6-73fbe2b3ad1f
# ╠═c9aa4a15-5fd4-4598-ad46-6680585965f7
# ╟─8954cc45-2feb-40fe-a622-7ccd13661ad9
# ╠═ff59a8f8-dc05-4b2e-be7f-3cef50ea5605
# ╠═86c80ab7-a7d9-4224-817b-8b84de3885a9
# ╟─fd18fb8a-f0d5-417a-9854-8d8f0f0bb3ce
# ╠═90586677-2b9a-4fef-b9ef-098eb288ee64
# ╟─19d5382d-a657-4993-bc87-75e9edfd1dbd
# ╠═ccf44020-d61a-4bb0-addd-4b27a9dbc0ad
# ╟─7cac3f5c-684f-4edc-b55f-6cf343556291
# ╠═9f3cce18-4af3-4270-803b-b6b1c773bb4c
# ╟─0f774f27-6251-46d3-a9e1-64e3d20ea9be
# ╠═2ac192f4-fa55-44bd-9b26-854a0f3266ca
