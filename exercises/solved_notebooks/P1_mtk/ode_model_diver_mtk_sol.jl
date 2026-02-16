### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 758484ec-9945-11f0-2834-e37380616e53
using Pkg; Pkg.activate("..")

# ╔═╡ 0ee79db9-1934-4c18-b0ed-28057a852614
using Markdown, InteractiveUtils

# ╔═╡ e430b466-ca2f-42f9-b9cf-f73c429506f6
using StatsPlots, PlutoUI; TableOfContents()

# ╔═╡ 8d55ad5c-cdfa-45dc-9abc-81cc7f1899be
using OrdinaryDiffEq, ModelingToolkit

# ╔═╡ 7d8cd22d-fc4e-4f4d-9640-fbe4f0fc075b
using ModelingToolkit: t_nounits as t, D_nounits as D

# ╔═╡ de2a2b81-1896-4faa-8936-37e0a176978e
md"""
# Exercise: Diver
"""

# ╔═╡ 047f4c07-527f-416d-b5c6-f9f7a26122e5
solution(text) = Markdown.MD(Markdown.Admonition("hint", "Solution", [text]));

# ╔═╡ ece4a09b-781a-48d1-994c-856674b6ef71
md"""
A diver is pulled up from a depth of $30\;m$ to a ship with a constant velocity $v$. The pressure change that takes place in the body is proportional to the difference between the ambient pressure (at a certain depth) and the internal pressure in the diver's body through a coefficient gradient of $c_2\;[s^{−1}]$. We are interested in the rate of pressure change in the body, as it has to stay below a certain critical value (estimated to be $0.02\;bar/s$) if we want to avoid the Caisson disease.
"""

# ╔═╡ 92e3c710-9cd2-4441-a395-f4303d3c122e
# PlutoUI.LocalResource("fig/diver_position.png")
# https://users.ugent.be/~gvhaelew/fig/diver_position.png
md"""
![Diver](https://users.ugent.be/~gvhaelew/fig/diver_position.png)
"""

# ╔═╡ 25ab9b2b-c594-4a48-ab77-910897597782
md"""
Consider a $Z$-axis pointing upwards with the water surface at position $z=0\;m$. Hence, the position of the diver (beneath the water surface) is negative ($z < 0$). The pressure that the water exerts on the diver's body is given by Pascal's law: $p = p_a - \rho g z$, with $p_a$ the air pressure. It is easy to verify that $p = p_a$ at $z = 0$.

**Important remaks:**
- When the diver is residing for a relatively long period of time at a certain position $z$, its internal body pressure $p_b$ will equal the ambient pressure $p$ ($= p_a - \rho \, g \, z$).
- When the diver is (suddenly) pulled up, it takes time for the internal body pressure of the diver to adapt to the (new) ambient pressure. Assume a linear relationship between the rate of change in internal body pressure and the difference between the ambient and the internal body pressure: $$\cfrac{dp_b}{dt} = c_2 \left(p - p_b\right) \, ,$$ with $c_2$ the adaption coefficient.
- Since the adaption of $p_b$ is not instant (cf. adaptation coefficient $c_2$), the factor $p - p_b$ will be negative, and hence, also $\cfrac{dp_b}{dt}$. It is the absolute value of $\cfrac{dp_b}{dt}$ that must stay below a certain critical value in order to avoid the Caisson disease.
"""

# ╔═╡ 8bba353f-557c-4b42-a12c-666fa35fbd47
md"""
## Part 1: constant speed
"""

# ╔═╡ 79bdf2b6-1b34-4543-902d-58b14adc9c92
md"""
#### Deriving the model equations
"""

# ╔═╡ 28a230f6-78be-4d0c-9258-bf7b8290e33a
md"""
!!! task
	Derive the model equations the position $z$ and the body pressure $p_b$.
"""

# ╔═╡ bccfb4eb-5b3d-4027-b367-d6f035f5db29
md"""
!!! hint
	- Set up:
	``` math
	\begin{align}
	\cfrac{d z}{dt} &= \cdots \\
	\cfrac{d p_b}{dt} &= \cdots
	\end{align}
	```
"""

# ╔═╡ 12744c81-58da-43b4-88f3-33190a57d79c
solution(md"""
``` math
\begin{align}
	\cfrac{d z}{dt} &= v \\
	\cfrac{d p_b}{dt} &= c_2 \left( p_a - \rho g z - p_b\right)
\end{align}
```
""")

# ╔═╡ 3c1ce356-19b0-4ac2-ad01-24be3f2baaa3
md"""
### Setting up the equations
"""

# ╔═╡ 67b148b0-777f-4d34-b5af-c8f92e869ebe
md"""
Consider for this part the velocity $v$ being constant with a value of $0.15\;m/s$. Set up a model for the position $z\;[m]$ and the body pressure $p_b\;[bar]$. Include in your model equations also a variable keeping track of the rate of change in body pressure $dp_b/dt\;[bar/s]$.
"""

# ╔═╡ c58a62a3-6519-4a81-a851-04ce6ec8a9db
md"""
Define variables for the position $z$, the body pressure $p_b$ and the rate of change in body pressure $dp_b/dt$.
"""

# ╔═╡ 7be849c3-d841-4967-936b-2e42dc332de5
# @variables missing
@variables z(t) pb(t) dpbdt(t)

# ╔═╡ d548f583-5f9b-4921-8e89-85995b892e51
md"""
Define the parameters for this model and assign their corresponding values.

| Parameter | Value | Unit | Meaning |
|:---------- |:---------- |:------------|:------------|
| $p_a$    | 101325*1e-5 | $bar$  |    Air pressure   |
| $\rho$    | 1000.0 | $kg/m^3$  |  water density  |
| $g$    | 9.81 | $m/s^2$  |   gravitational constant |
| $c_2$    | 0.05 | $1/s$  |   coefficient gradient |
| $v$    | 0.15 | $m/s$  |   pull-up velocity |
"""

# ╔═╡ fb444a48-63ec-4bf4-aa55-8e8ae04640cd
# @parameters missing
@parameters pa=101325*1e-5 ρ=1000.0 g=9.81 c₂=0.05 v_const=0.15

# ╔═╡ 1916a1ca-565b-4692-82c2-5484a486e435
md"""
Define the initial diver position.
"""

# ╔═╡ 9a289ab4-87b5-452e-890e-2cad5c16ab4a
# z₀ = missing
z₀ = -30.0

# ╔═╡ 3a4501ed-79dc-4b6e-ba09-df08d78f8924
md"""
Define the initial body pressure. Don't forget to include `*1e-5` in the hydrostatic term in order to have everything in $bar$.
"""

# ╔═╡ 07b8e2e3-3fad-4602-a3f0-55461ed1b5fd
# pb₀ = missing
pb₀ = pa - ρ*g*z₀*1e-5

# ╔═╡ 962f3a24-cecf-4c99-9850-4706d41ed70f
md"""
In order to have an idea of this initial body pressure, it is preferable to calculate it in a `let` ... `end` block. Everything in this block will be local scope and won't interfere with variable names defined elsewhere.
"""

# ╔═╡ 3b460c68-cd31-4adf-985a-057508b40105
# let
# 	pa=101325*1e-5; ρ=1000.0; g=9.81;
# 	missing
# end
let
	pa=101325*1e-5; ρ=1000.0; g=9.81;
	pb₀ = pa - ρ*g*z₀*1e-5
end

# ╔═╡ fd38ceaf-ca26-43a7-8e8b-52aa63f389a9
md"""
Set up the equation for the rate of change in position.
"""

# ╔═╡ 3c1b7030-9e95-421f-b818-ede261ff8ff9
# eq1_position = missing
eq1_position = D(z) ~ v_const

# ╔═╡ 57aa99f6-bfa0-46a5-9565-ecc4ea04dffb
md"""
Set up the equation for the rate of change in body pressure. Don't forget to include `*1e-5` in one of the terms in order to have the pressure in $bar$.
"""

# ╔═╡ ba1b330b-d57d-4f40-841d-a05840a4c059
# eq1_body_pressure = missing
eq1_body_pressure = D(pb) ~ c₂*(pa - ρ*g*z*1e-5 - pb)

# ╔═╡ 3a43debf-ed0b-4fcf-8aaf-496de76c3914
md"""
Set up the equation that will keep track of the rate of change in body pressure.
"""

# ╔═╡ e59c36db-83c3-4ac8-81a1-2a5c175c75ff
# eq1_dpbdt = missing
eq1_dpbdt = dpbdt ~ c₂*(pa - ρ*g*z*1e-5 - pb)

# ╔═╡ f42f55a5-4063-4f14-ae30-0d0468e32f61
md"""
Bundle all equations.
"""

# ╔═╡ 32f4e0d7-b512-409a-9831-fb3073fe3c3d
# eqns1_diver = missing
eqns1_diver = [eq1_position, eq1_body_pressure, eq1_dpbdt]

# ╔═╡ ef793c81-0e60-46af-b681-376265f333e8
md"""
### Building the ODE system
"""

# ╔═╡ 4ca01607-0a0a-4ada-96de-44286565f6be
md"""
Build the model and include the (continuous) event that when $z$ hits zero, the pull-up velocity must become zero as well.
"""

# ╔═╡ fafc4370-569b-4e05-9683-8e54340713d2
# @mtkbuild missing
@mtkbuild sys1_diver = ODESystem(eqns1_diver, t; continuous_events = [[z ~ 0]=>[v_const ~ 0]])

# ╔═╡ e8e2b1a6-a436-45c0-a4d3-1a7baf0ec99d
md"""
### Create and solve the ODE problem
"""

# ╔═╡ d42afda6-edd4-4fb6-93b9-8663adc0c96d
md"""
Create the ODE problem, specifying the initial values for the diver's position and internal body pressure. Use a simulation time span of `500.0` seconds. Use an empty vector `[]` for the parameters argument to use their default values.
"""

# ╔═╡ df8fda53-7226-409f-89a4-7375a26468e6
# oprob1_diver = missing
oprob1_diver = ODEProblem(sys1_diver, [z=>z₀, pb=>pb₀], (0.0, 500.0), [])

# ╔═╡ 84f365dd-04f8-4763-b5a3-95226eba9a98
md"""
Solve the ODE problem. Make a deepcopy of the ODE problem, use `Tsit5()` and `saveat=1`.
"""

# ╔═╡ 700103f1-e0a7-4d09-8b1e-c152e9485c65
# sol1_diver = missing
sol1_diver = solve(deepcopy(oprob1_diver), Tsit5(), saveat=1)  # , reltol=1e-9

# ╔═╡ b45e8caf-e0e6-4484-b858-1129dd5b7834
md"""
### Plotting results
"""

# ╔═╡ df0f1383-425e-4f4b-8e30-54e072fe481e
md"""
Plot the position of the diver over time. Include a y-label (cf. `ylabel="..."`) and a title (cf. `title="..."`) in all plots.
"""

# ╔═╡ abfe5216-48b7-4c57-8c51-79c49493e419
# missing
plot(sol1_diver; idxs=[z], ylabel="z [m]",title="Position")

# ╔═╡ 53dd78a4-40d3-48e7-94fb-731a53f6eb79
md"""
Plot the body pressure of the diver over time.
"""

# ╔═╡ dd525893-627f-4e40-8674-e80cb51da255
# missing
plot(sol1_diver; idxs=[pb], ylabel="pb [bar]", title="Body pressure")

# ╔═╡ 6f7d3199-d223-4069-a41c-aa774566d990
md"""
Plot the rate of change in body pressure over time.
"""

# ╔═╡ 97eff5f5-5c6f-4cf6-9d3e-8ceb216eab63
# missing
plot(sol1_diver; idxs=[dpbdt], ylabel="dpb/dt [bar/s]", title="Rate of change in body pressure")

# ╔═╡ f3d946f8-54e0-4e48-b6a9-4b5eb89ca412
md"""
Calculate the maximum rate of change in the body pressure in absolute value.
"""

# ╔═╡ c8481d75-945a-4d89-b3d2-f32f5e65bef2
# missing
maximum(abs.(sol1_diver[dpbdt]))

# ╔═╡ 087a0928-e8f9-48a4-b946-c517a404d3dc
md"""
!!! question
	Is this maximum rate of change in body pressure for a diver that is pulled up at a speed of $0.15\;m/s$ a safe value?
"""

# ╔═╡ a3ec7991-b130-4917-b5f7-2b77a45fffb2
md"""
Answer: missing
"""

# ╔═╡ a41f2eec-5a73-4510-b189-f8d11d2aa771
md"""
## Part 2: external forces
"""

# ╔═╡ 1d07eb14-dc16-42b0-9707-2499cdcd106e
md"""
Suppose that the diver is now attached with a rope to a motor on a boat crane to be pulled up.

Four forces now act on the diver: gravity, a frictional force, the force exerted by acceleration and the Archimedean force. The frictional force on the diver is proportional to the speed of the diver through a coefficient $c_1,[kg/s]$, while the upward Archimedean force is given by $F_{arch}=\rho V g$, with $V$ the volume of the diver.
"""

# ╔═╡ 91734c50-f848-448b-ad13-5b26a41fa0eb
md"""
### Deriving the model equations
"""

# ╔═╡ 9bb2ec9d-9661-4808-85b9-2b46066a5659
md"""
In this case the velocity $v$ won't be constant in the beginning but will be governed by the external force.

**You need to take into account the following forces:**
- The external force: $\vec{F}_{ext} = F_{ext}\, \vec{e}_z$. This force is directed upward.
- The buoyuancy force: $\vec{F}_b = \rho \, V \, g \, \vec{e}_z$. This force is always directed upward.
- The weight of the diver: $\vec{F}_g = -m\,g\,\vec{e}_z$. This force is always directed downward.
- The friction force between the diver and the water: $\vec{F}_w = -c_1\,v \,\vec{e}_z$. This force is always opposite the movement direction of the diver. In our case the diver is pulled up, hence, this force is directed downward (remember that $v=\cfrac{dz}{dt} > 0$).
"""

# ╔═╡ 07e2d7c1-05b6-4309-a524-782dfe6bc14c
md"""
!!! task
	Derive the model equations for the position $z$, the velocity $v$ and the body pressure $p_b$.
"""

# ╔═╡ 7f32442f-3028-457c-9013-39248fb26cb9
md"""
!!! hint
	- Start setting up:
	``` math
	\begin{align}
	m\cfrac{d^2 z}{dt^2} &= \cdots \\
	\cfrac{d p_b}{dt} &= \cdots
	\end{align}
	```
	and work toward:
	``` math
	\begin{align}
	\cfrac{d z}{dt} &= \cdots \\
	m \cfrac{d v}{dt} &= \cdots \\
	\cfrac{d p_b}{dt} &= \cdots
	\end{align}
	```
"""

# ╔═╡ 0dbe6003-b50b-41af-bee7-201bf67bf7ca
solution(md"""
``` math
\begin{align}
	\cfrac{d z}{dt} &= v \\
	m \cfrac{d v}{dt} &= F_{ext} + \left(\rho V - m \right)g - c_1 v\\
	\cfrac{d p_b}{dt} &= c_2 \left( p_a - \rho g z - p_b\right)
\end{align}
```
""")

# ╔═╡ 04d9a2cb-5f73-4ac0-8f93-1733f7f71840
md"""
### Setting up the equations
"""

# ╔═╡ 180a67b1-a03e-4dbc-9331-d82e688e49d0
md"""
Define an additional variable for the velocity $v$.
"""

# ╔═╡ c78d0438-ba53-4f37-ac5e-9a32172c501b
# @variables missing
@variables v(t)

# ╔═╡ 588e1ad2-181d-4a1f-b974-a7f5f0f70c28
md"""
Define the additional parameters for this model and assign their corresponding values.

| Parameter | Value | Unit | Meaning |
|:---------- |:---------- |:------------|:------------|
| $m$    | 100 | $kg$  |    mass of the diver   |
| $V$    | 0.082 | $m^3$  |  volume of the diver  |
| $c_1$    | 20.0 | $kg/s$  |   friction coefficient |
"""

# ╔═╡ 1d3a4fba-9593-4503-8656-ffbceb8046ec
# @parameters missing
@parameters m=100.0 V=0.082 c₁=20.0

# ╔═╡ e88c9c4e-d50b-436a-a810-76fe40ddc24d
md"""
In order to have an idea of what external force you need to pull up the diver at a constant velocity of 0.15 m/s, you can calculate it by setting $\cfrac{dv}{dt}$ to zero in the equation for the rate of change in the velocity, and solving for $F_{ext}$.
"""

# ╔═╡ 2d5c2f05-c902-4b6d-8b73-38d092a11fd8
# let
# 	ρ=1000.0; g=9.81; v_const=0.15;
# 	m=100.0; V=0.082; c₁=20.0;
# 	missing
# end
let
	ρ=1000.0; g=9.81; v_const=0.15;
	m=100.0; V=0.082; c₁=20.0;
	Fext = c₁*v_const - ρ*g*V + m*g
end

# ╔═╡ 9db90eab-683d-4c7c-a0b9-5c7904dd0aa2
md"""
Set up the additional parameter for the external force with a default value of 180.
"""

# ╔═╡ c8b34161-ecb6-4be7-a4c4-d0381198846c
# @parameters missing
@parameters Fext=180.0

# ╔═╡ acf469f8-9765-47e9-b4f5-55715a4075c6
md"""
Set up the equation for the rate of change in position. Multiply the right handside of the equation with $(z < 0)$ to make sure that after $z$ hits $0$, $z$ remains constant and, hence, also $v$ remains $0$.
"""

# ╔═╡ 81bba0e7-fc3a-4304-987a-099a43fb38eb
# eq2_position = missing
eq2_position = D(z) ~ v*(z<0)

# ╔═╡ 03711ec2-bab0-4917-a40a-01ac2904612b
md"""
Set up the equation for the rate of change in velocity. This includes the external force, the buoyuancy force, the weight of the diver and the frictional force.
"""

# ╔═╡ 5a002b19-5986-4adb-84d9-90e95c655a32
# eq2_velocity = missing
eq2_velocity = m*D(v) ~ Fext + (ρ*V - m)*g - c₁*v
# eq2_velocity = D(v) ~ Fext/m + (ρ*V/m - 1)*g - c₁/m*v   # alternative

# ╔═╡ 4b3de655-674e-4a13-83f3-f09b6443b26a
md"""
Set up the equation for the rate of change in body pressure. Don't forget to include `*1e-5` in one of the terms in order to have everything in $bar$.
"""

# ╔═╡ c50c8b9c-0e1c-4f7a-800a-33220e83dfed
# eq2_body_pressure = missing
eq2_body_pressure = D(pb) ~ c₂*(pa - ρ*g*z*1e-5 - pb)

# ╔═╡ b6c8b269-b215-4213-8248-4a0dee79a43a
md"""
Set up the equation that will keep track of the rate of change in body pressure.
"""

# ╔═╡ bbf20a87-77f2-4818-b0bb-ad67b30a5bbc
# eq2_dpbdt = missing
eq2_dpbdt = dpbdt ~ c₂*(pa - ρ*g*z*1e-5 - pb)

# ╔═╡ 37f10ef2-d46f-499e-bbad-f0df56b70540
md"""
Bundle all equations.
"""

# ╔═╡ 9e29b21e-20c9-4b98-84d8-f76b42c627c8
# eqns2_diver = missing
eqns2_diver = [eq2_position, eq2_velocity, eq2_body_pressure, eq2_dpbdt]

# ╔═╡ e3483144-c9c5-4fdb-8053-fe0c08ae8765
md"""
### Building the ODE system
"""

# ╔═╡ 9827f9e4-d514-454e-893f-42f9c97b5a43
md"""
Build the model and include the (continuous) event that when $z$ hits zero, $v$ must become zero as well.
"""

# ╔═╡ 70a6a572-c25a-40bd-824a-80b21957efbb
# @mtkbuild missing
@mtkbuild sys2_diver = ODESystem(eqns2_diver, t; continuous_events = [[z ~ 0]=>[v ~ 0]])

# ╔═╡ fd3ae56c-34fe-4615-b73e-088e490fbb59
md"""
### Create and solve the ODE problem
"""

# ╔═╡ 4bb1e0f4-2544-437b-8523-4e5547caee25
md"""
Create the ODE problem. The initial conditions for the position and body pressure are the same as before. For $v$, assume that the diver is initially at rest. Use a simulation time span of `500.0` seconds. Use the default values of the parameters again.
"""

# ╔═╡ 26eb5e36-4026-41b8-b2fc-4be950727e29
# oprob2_diver = missing
oprob2_diver = ODEProblem(sys2_diver, [z=>z₀, v=>0, pb=>pb₀], (0.0, 500.0), [])

# ╔═╡ da2cd981-9b9e-42dc-b1b3-9ac09b53c5e1
md"""
Solve the ODE problem. Make a deepcopy of the ODE problem, use `Tsit5()` and `saveat=1`.
"""

# ╔═╡ 61fdda96-1c27-42aa-b6d3-dd0b2f9dae6c
# sol2_diver = missing
sol2_diver = solve(oprob2_diver, Tsit5(), saveat=1, reltol=1e-9)

# ╔═╡ 790ce55b-fffb-49d7-8013-755cb647cf7d
md"""
### Plotting results
"""

# ╔═╡ d2e5c67a-fb7a-4cf1-a849-ff11734703fe
md"""
Plot the position of the diver over time. Include a y-label (cf. `ylabel="..."`) and a title (cf. `title="..."`) in all plots.
"""

# ╔═╡ 5e1faec8-d39b-4d41-be7e-7de077b807fe
# missing
plot(sol2_diver; idxs=[z], ylabel="z [m]",title="Position")

# ╔═╡ dd8ba91e-7e9c-4857-8682-b01dbc7a595e
md"""
Plot the velocity of the diver over time.
"""

# ╔═╡ 4f9c301a-a30e-4373-a409-70af890d52c3
# missing
plot(sol2_diver; idxs=[v], ylabel="v [m/s]", title="Velocity")

# ╔═╡ 60845e85-fb5e-4e88-b880-c51c49166684
md"""
Plot the body pressure of the diver over time.
"""

# ╔═╡ 8e0fd1bd-4141-49ad-b7bc-13612378e47e
# missing
plot(sol2_diver; idxs=[pb], ylabel="pb [bar]", title="Body pressure")

# ╔═╡ 58ba6450-56e4-48a8-8629-a84b56a53269
md"""
Plot the rate of change in body pressure over time.
"""

# ╔═╡ a6fc94ca-76c7-478a-82ef-256ad748dac0
# missing
plot(sol2_diver; idxs=[dpbdt], ylabel="dpb/dt [bar/s]", title="Rate of change in body pressure")

# ╔═╡ 5cf2ebb6-0a84-40aa-9593-43370e665dc1
md"""
Calculate the maximum rate of change in the body pressure in absolute value.
"""

# ╔═╡ 852e421f-0a7d-42aa-b8d4-d903407e170b
# missing
maximum(abs.(sol2_diver[dpbdt]))

# ╔═╡ c10eb86b-9637-458d-aa60-40e66924fbac
md"""
!!! question
	Is this maximum rate of change in body pressure for a diver that is pulled up with $F_{ext} = 180.0\;N$ a safe value?
"""

# ╔═╡ 66178148-0e65-497d-b001-a6421faee5f1
md"""
Answer: missing
"""

# ╔═╡ 0f8ba58d-2632-4892-8d9c-b95a905aad6a
md"""
## Part 3: plot of maximum dpb/dt vs Fext
"""

# ╔═╡ 942faab9-9d07-4322-ba1a-d1b614f8abd7
md"""
Make a plot of the maximum rate of change in body pressure versus the external force in the range [177, 187] $N$ with a step size of 0.5 $N$. Use the same initial conditions and time span as in Part 2 but iteratively modify the parameter value of $F_{ext}$. Instead of `Tsit5()`, use now the `Rosenbrock32()` solver.
"""

# ╔═╡ dd93fc6c-d0f2-431d-a0fc-1c744cf85850
md"""
!!! hints
	- Append during each iteration the maximum rate of change in body pressure to the vector `dpbdt_max_vals`.
	- Define a range object for the external forces as `Fext_vals = 178:0.5:184`.
	- Use the model `sys2_diver` from Part 2, but name the ODE problem as `oprob3_diver`.
	- Change the parameter value of `Fext` at each iteration step.
"""

# ╔═╡ 693c4b30-8bcc-4ee2-9557-6fcd638a7888
# begin
# 	dpbdt_max_vals = [];
# 	Fext_vals = missing
# 	for Fext_val in Fext_vals
# 		oprob3_diver = missing
# 		sol3_diver = missing
# 		append!(dpbdt_max_vals, missing)
# 	end
# end
begin
	dpbdt_max_vals = [];
	Fext_vals = 178:0.5:184
	for Fext_val in Fext_vals
		oprob3_diver = ODEProblem(sys2_diver, [z=>z₀, v=>0, pb=>pb₀], (0.0, 500.0), [Fext=>Fext_val])
		sol3_diver = solve(oprob3_diver, Rosenbrock32(), saveat=1)
		append!(dpbdt_max_vals, maximum(abs.(sol3_diver[dpbdt])))
	end
end

# ╔═╡ ad818822-45c3-40eb-ba6d-f0279aa0b457
md"""
Plot the maximum rate of change in body pressure vs the external force.
"""

# ╔═╡ 3ed42f54-55c7-48b8-83c0-dec8fd089457
# missing
plot(Fext_vals, dpbdt_max_vals)

# ╔═╡ c38e90f9-6e27-463e-a759-0f1153a6f6be
# dpbdt_max_vals .>= 0.02

# ╔═╡ 1828807f-ae09-425a-a310-98a22b7aa3b2
# i002 = findfirst(==(true), dpbdt_max_vals .>= 0.02)

# ╔═╡ 68d90aac-153b-4ed3-a7fc-76edd3ea96b6
# Fext_vals[i002]

# ╔═╡ 7a0ae4a1-85cc-4000-92ac-8cce972ac999
md"""
!!! task
	Incept graphically from the plot the value of the external force so that the rate of change in body pressure is 0.02.
"""

# ╔═╡ 40fbc323-049f-40b9-b44d-26fce4135c19
md"""
Response: missing
"""

# ╔═╡ 9c611e9f-d355-48f3-b5d6-c8acd186151a
md"""
## Part 4: more accurate value of Fext
"""

# ╔═╡ c06096c4-3f44-4ddc-9f25-f916d27e2a5b
md"""
Find a more accurate value for the external force so that the maximum rate of change in body pressure is 0.02.
"""

# ╔═╡ 2e521e5f-ef58-4d5b-9fd6-2a87c77cb276
md"""
!!! hints
	- Use a `while` loop.
	- The start value of the external force is defined as `Fext_val = 180.0`, its step size as `DFext_val = 0.05`.
	- Assign within the `while` loop the value of the maximum rate of change in body pressure to `dpbdt_max_val`.
	- Use the model `sys2_diver` from Part 2, but name the ODE problem as `oprob4_diver`.
	- Within the `while` loop you need to put `global` before `dpbdt_max_val` and `Fext_val` because they were define outside the `while` loop.
"""

# ╔═╡ 52d10672-20a4-4516-a632-d49e54d8fc4b
# begin
# 	Fext_val = missing
# 	DFext_val = missing
# 	dpbdt_max_val = 0.00
# 	while missing
# 		oprob4_diver = missing
# 		sol4_diver = missing
# 		global dpbdt_max_val = missing
# 		global Fext_val += missing
# 	end
# end
begin
	Fext_val = 180.0
	DFext_val = 0.05
	dpbdt_max_val = 0.00
	while dpbdt_max_val < 0.02
		oprob4_diver = ODEProblem(sys2_diver, [z=>z₀, v=>0, pb=>pb₀], (0.0, 500.0), [Fext=>Fext_val])   # consider remake
		sol4_diver = solve(oprob4_diver, Rosenbrock32(), saveat=1, reltol=1e-9)
		global dpbdt_max_val = maximum(abs.(sol4_diver[dpbdt]))
		global Fext_val += DFext_val
	end
end

# ╔═╡ 54c340f7-5f7e-4bd1-bcfd-9b10c0fb47ff
md"""
Show the value of the retrieved external force.
"""

# ╔═╡ b87163b4-9efa-46ba-9420-94979c49abc6
# missing
Fext_val

# ╔═╡ 1092d818-d725-4a34-96f9-adc1f686c039
md"""
## Part 5: above the water surface
"""

# ╔═╡ dfecaf29-1be5-479b-bcd2-9d21643e682c
md"""
Suppose that the diver is now pulled up to a height of $5\;m$ above the water surface. The velocity of the diver should be constant and the same as when the diver hits the surface.
"""

# ╔═╡ 93f63315-7f6d-43cb-b95d-d12ce103f9ff
md"""
### Setting up the equations
"""

# ╔═╡ 499331aa-5062-4f9a-a4ff-9f7597be6268
md"""
Set up the equation for the rate of change in position. Multiply the right handside of the equation with `(z<5)` to make sure that after $z$ hits $5$, $z$ remains constant and, hence, also $v$ remains $0$.
"""

# ╔═╡ 1e34dfa5-b9ca-479c-8a86-fdd064dfe408
# eq5_position = missing
eq5_position = D(z) ~ v*(z<5)

# ╔═╡ 8947f891-b3c5-4472-b312-dba04face3d5
md"""
Set up the equation for the rate of change in velocity. This equation only holds when $z < 0$, hence multiply it with `(z<0)`. When $0 \leq z < 5$, the velocity should remain constant, hence, $\cfrac{dv}{dt} = 0$.
"""

# ╔═╡ 0fa60678-66c6-442a-8ab0-63c3060dbd6e
# eq5_velocity = missing
eq5_velocity = m*D(v) ~ (Fext + (ρ*V - m)*g - c₁*v)*(z<0)
# eq5_velocity = D(v) ~ (Fext/m + (ρ*V/m - 1)*g - c₁/m*v)*(z<0)   # alternative

# ╔═╡ 1a9e0ba4-52b6-4dc0-844b-88949e69befc
md"""
Set up the equation for the rate of change in body pressure. The term with $\rho\,g\,z$ has to do with the hydrostatic pressure and should vanish when the diver is above the water surface. Hence, include the factor `(z<0)` to this term.
"""

# ╔═╡ 61959657-80cd-4b73-b57b-4dc11a9fd6d2
# eq5_body_pressure = missing
eq5_body_pressure = D(pb) ~ c₂*(pa - ρ*g*z*1e-5*(z<0) - pb)

# ╔═╡ 70b8a56e-6377-44bc-b5a4-af2a006d3dd3
md"""
Set up the equation that will keep track of the rate of change in body pressure.
"""

# ╔═╡ 5cb47905-7379-4453-900b-92bf98a94ad1
# eq5_dpbdt = missing
eq5_dpbdt = dpbdt ~ c₂*(pa - ρ*g*z*1e-5*(z<0) - pb)

# ╔═╡ c807bf81-afc7-46b0-b073-26a27da9a55c
md"""
Bundle all equations.
"""

# ╔═╡ 3b11cd98-0317-4a9f-87af-7195c71e4fbc
# eqns5_diver = missing
eqns5_diver = [eq5_position, eq5_velocity, eq5_body_pressure, eq5_dpbdt]

# ╔═╡ 9fa0030c-6099-4fad-9d5b-800d5e197cf5
md"""
### Building the ODE system
"""

# ╔═╡ 6d13fdf7-7214-483b-9de3-9fd7accab4a1
md"""
Build the model and include the (continuous) event that when $z$ hits $5$, $v$ must become zero.
"""

# ╔═╡ 6a87b974-4ce6-45b0-b7e8-ce007c45659c
# @mtkbuild missing
@mtkbuild sys5_diver = ODESystem(eqns5_diver, t; continuous_events = [[z ~ 5.0]=>[v ~ 0]])

# ╔═╡ abad5c96-060d-4a37-8da0-6d6b97dd1926
md"""
### Create and solve the ODE problem
"""

# ╔═╡ 68a2f5cf-3e0b-4dbc-ad66-ab728cbe7984
md"""
Create the ODE problem. All initial conditions are the same as before. Use a simulation time span of `500.0` seconds. Use the default parameter values.
"""

# ╔═╡ 0a447fbc-c31d-4eac-b272-65818efb62d7
# oprob5_diver = missing
oprob5_diver = ODEProblem(sys5_diver, [z=>z₀, v=>0, pb=>pb₀], (0.0, 500.0), [])

# ╔═╡ 56ba3f8d-2c63-4b04-a389-37ff62e125d5
md"""
Solve the ODE problem. Make a deepcopy of the ODE problem, use `Tsit5()`, `saveat=1` and `reltol=1e-9`.
"""

# ╔═╡ 1fd0aa0c-94d0-465d-a213-f788fa9572b7
# sol5_diver = missing
sol5_diver = solve(oprob5_diver, Tsit5(), saveat=1, reltol=1e-9)

# ╔═╡ 269b83aa-e7d2-4b04-b141-40997b6ba2a9
md"""
### Plotting results
"""

# ╔═╡ 251ac174-9ed8-41dc-be7c-db919e29515e
md"""
Plot the position of the diver over time. Include y-labels (cf. `ylabel="..."`) and titles (cf. `title="..."`).
"""

# ╔═╡ 893e77ec-6a36-47f2-9c84-6e723fe1a542
# missing
plot(sol5_diver; idxs=[z], ylabel="z [m]", title="Position")

# ╔═╡ 905f0f37-5443-434d-899d-a78e46bb1344
# missing
plot(sol5_diver; idxs=[v], ylabel="v [m/s]", title="Velocity")

# ╔═╡ 0f9b59cc-a8c6-4443-9080-aae4d8cd967b
# missing
plot(sol5_diver; idxs=[pb], ylabel="pb [bar]", title="Body pressure")

# ╔═╡ ff722e23-a476-4149-832f-ca875112a302
# missing
plot(sol5_diver; idxs=[dpbdt], ylabel="dpb/dt [bar/s]", title="Rate of change in body pressure")

# ╔═╡ 2a5e302e-2c44-4fdc-bd9f-cbfef0df5763
md"""
!!! questions
	1. Which plots do you expect to be slightly different compared to those in Part 2.
	2. Are the plots that are different according to your expectations?
"""

# ╔═╡ c0335f76-d167-4444-9bcc-27f0dc825db2
md"""
Answers:
1. missing
2. missing
"""

# ╔═╡ Cell order:
# ╟─de2a2b81-1896-4faa-8936-37e0a176978e
# ╠═0ee79db9-1934-4c18-b0ed-28057a852614
# ╠═758484ec-9945-11f0-2834-e37380616e53
# ╠═e430b466-ca2f-42f9-b9cf-f73c429506f6
# ╠═8d55ad5c-cdfa-45dc-9abc-81cc7f1899be
# ╠═7d8cd22d-fc4e-4f4d-9640-fbe4f0fc075b
# ╟─047f4c07-527f-416d-b5c6-f9f7a26122e5
# ╟─ece4a09b-781a-48d1-994c-856674b6ef71
# ╟─92e3c710-9cd2-4441-a395-f4303d3c122e
# ╟─25ab9b2b-c594-4a48-ab77-910897597782
# ╟─8bba353f-557c-4b42-a12c-666fa35fbd47
# ╟─79bdf2b6-1b34-4543-902d-58b14adc9c92
# ╟─28a230f6-78be-4d0c-9258-bf7b8290e33a
# ╟─bccfb4eb-5b3d-4027-b367-d6f035f5db29
# ╟─12744c81-58da-43b4-88f3-33190a57d79c
# ╟─3c1ce356-19b0-4ac2-ad01-24be3f2baaa3
# ╟─67b148b0-777f-4d34-b5af-c8f92e869ebe
# ╟─c58a62a3-6519-4a81-a851-04ce6ec8a9db
# ╠═7be849c3-d841-4967-936b-2e42dc332de5
# ╟─d548f583-5f9b-4921-8e89-85995b892e51
# ╠═fb444a48-63ec-4bf4-aa55-8e8ae04640cd
# ╟─1916a1ca-565b-4692-82c2-5484a486e435
# ╠═9a289ab4-87b5-452e-890e-2cad5c16ab4a
# ╟─3a4501ed-79dc-4b6e-ba09-df08d78f8924
# ╠═07b8e2e3-3fad-4602-a3f0-55461ed1b5fd
# ╟─962f3a24-cecf-4c99-9850-4706d41ed70f
# ╠═3b460c68-cd31-4adf-985a-057508b40105
# ╟─fd38ceaf-ca26-43a7-8e8b-52aa63f389a9
# ╠═3c1b7030-9e95-421f-b818-ede261ff8ff9
# ╟─57aa99f6-bfa0-46a5-9565-ecc4ea04dffb
# ╠═ba1b330b-d57d-4f40-841d-a05840a4c059
# ╟─3a43debf-ed0b-4fcf-8aaf-496de76c3914
# ╠═e59c36db-83c3-4ac8-81a1-2a5c175c75ff
# ╟─f42f55a5-4063-4f14-ae30-0d0468e32f61
# ╠═32f4e0d7-b512-409a-9831-fb3073fe3c3d
# ╟─ef793c81-0e60-46af-b681-376265f333e8
# ╟─4ca01607-0a0a-4ada-96de-44286565f6be
# ╠═fafc4370-569b-4e05-9683-8e54340713d2
# ╟─e8e2b1a6-a436-45c0-a4d3-1a7baf0ec99d
# ╟─d42afda6-edd4-4fb6-93b9-8663adc0c96d
# ╠═df8fda53-7226-409f-89a4-7375a26468e6
# ╟─84f365dd-04f8-4763-b5a3-95226eba9a98
# ╠═700103f1-e0a7-4d09-8b1e-c152e9485c65
# ╟─b45e8caf-e0e6-4484-b858-1129dd5b7834
# ╟─df0f1383-425e-4f4b-8e30-54e072fe481e
# ╠═abfe5216-48b7-4c57-8c51-79c49493e419
# ╟─53dd78a4-40d3-48e7-94fb-731a53f6eb79
# ╠═dd525893-627f-4e40-8674-e80cb51da255
# ╟─6f7d3199-d223-4069-a41c-aa774566d990
# ╠═97eff5f5-5c6f-4cf6-9d3e-8ceb216eab63
# ╟─f3d946f8-54e0-4e48-b6a9-4b5eb89ca412
# ╠═c8481d75-945a-4d89-b3d2-f32f5e65bef2
# ╟─087a0928-e8f9-48a4-b946-c517a404d3dc
# ╠═a3ec7991-b130-4917-b5f7-2b77a45fffb2
# ╟─a41f2eec-5a73-4510-b189-f8d11d2aa771
# ╟─1d07eb14-dc16-42b0-9707-2499cdcd106e
# ╟─91734c50-f848-448b-ad13-5b26a41fa0eb
# ╟─9bb2ec9d-9661-4808-85b9-2b46066a5659
# ╟─07e2d7c1-05b6-4309-a524-782dfe6bc14c
# ╟─7f32442f-3028-457c-9013-39248fb26cb9
# ╟─0dbe6003-b50b-41af-bee7-201bf67bf7ca
# ╟─04d9a2cb-5f73-4ac0-8f93-1733f7f71840
# ╟─180a67b1-a03e-4dbc-9331-d82e688e49d0
# ╠═c78d0438-ba53-4f37-ac5e-9a32172c501b
# ╟─588e1ad2-181d-4a1f-b974-a7f5f0f70c28
# ╠═1d3a4fba-9593-4503-8656-ffbceb8046ec
# ╟─e88c9c4e-d50b-436a-a810-76fe40ddc24d
# ╠═2d5c2f05-c902-4b6d-8b73-38d092a11fd8
# ╟─9db90eab-683d-4c7c-a0b9-5c7904dd0aa2
# ╠═c8b34161-ecb6-4be7-a4c4-d0381198846c
# ╟─acf469f8-9765-47e9-b4f5-55715a4075c6
# ╠═81bba0e7-fc3a-4304-987a-099a43fb38eb
# ╟─03711ec2-bab0-4917-a40a-01ac2904612b
# ╠═5a002b19-5986-4adb-84d9-90e95c655a32
# ╟─4b3de655-674e-4a13-83f3-f09b6443b26a
# ╠═c50c8b9c-0e1c-4f7a-800a-33220e83dfed
# ╟─b6c8b269-b215-4213-8248-4a0dee79a43a
# ╠═bbf20a87-77f2-4818-b0bb-ad67b30a5bbc
# ╟─37f10ef2-d46f-499e-bbad-f0df56b70540
# ╠═9e29b21e-20c9-4b98-84d8-f76b42c627c8
# ╟─e3483144-c9c5-4fdb-8053-fe0c08ae8765
# ╟─9827f9e4-d514-454e-893f-42f9c97b5a43
# ╠═70a6a572-c25a-40bd-824a-80b21957efbb
# ╟─fd3ae56c-34fe-4615-b73e-088e490fbb59
# ╟─4bb1e0f4-2544-437b-8523-4e5547caee25
# ╠═26eb5e36-4026-41b8-b2fc-4be950727e29
# ╟─da2cd981-9b9e-42dc-b1b3-9ac09b53c5e1
# ╠═61fdda96-1c27-42aa-b6d3-dd0b2f9dae6c
# ╟─790ce55b-fffb-49d7-8013-755cb647cf7d
# ╟─d2e5c67a-fb7a-4cf1-a849-ff11734703fe
# ╠═5e1faec8-d39b-4d41-be7e-7de077b807fe
# ╟─dd8ba91e-7e9c-4857-8682-b01dbc7a595e
# ╠═4f9c301a-a30e-4373-a409-70af890d52c3
# ╟─60845e85-fb5e-4e88-b880-c51c49166684
# ╠═8e0fd1bd-4141-49ad-b7bc-13612378e47e
# ╟─58ba6450-56e4-48a8-8629-a84b56a53269
# ╠═a6fc94ca-76c7-478a-82ef-256ad748dac0
# ╟─5cf2ebb6-0a84-40aa-9593-43370e665dc1
# ╠═852e421f-0a7d-42aa-b8d4-d903407e170b
# ╟─c10eb86b-9637-458d-aa60-40e66924fbac
# ╠═66178148-0e65-497d-b001-a6421faee5f1
# ╟─0f8ba58d-2632-4892-8d9c-b95a905aad6a
# ╟─942faab9-9d07-4322-ba1a-d1b614f8abd7
# ╟─dd93fc6c-d0f2-431d-a0fc-1c744cf85850
# ╠═693c4b30-8bcc-4ee2-9557-6fcd638a7888
# ╟─ad818822-45c3-40eb-ba6d-f0279aa0b457
# ╠═3ed42f54-55c7-48b8-83c0-dec8fd089457
# ╠═c38e90f9-6e27-463e-a759-0f1153a6f6be
# ╠═1828807f-ae09-425a-a310-98a22b7aa3b2
# ╠═68d90aac-153b-4ed3-a7fc-76edd3ea96b6
# ╟─7a0ae4a1-85cc-4000-92ac-8cce972ac999
# ╠═40fbc323-049f-40b9-b44d-26fce4135c19
# ╟─9c611e9f-d355-48f3-b5d6-c8acd186151a
# ╟─c06096c4-3f44-4ddc-9f25-f916d27e2a5b
# ╟─2e521e5f-ef58-4d5b-9fd6-2a87c77cb276
# ╠═52d10672-20a4-4516-a632-d49e54d8fc4b
# ╟─54c340f7-5f7e-4bd1-bcfd-9b10c0fb47ff
# ╠═b87163b4-9efa-46ba-9420-94979c49abc6
# ╟─1092d818-d725-4a34-96f9-adc1f686c039
# ╟─dfecaf29-1be5-479b-bcd2-9d21643e682c
# ╟─93f63315-7f6d-43cb-b95d-d12ce103f9ff
# ╟─499331aa-5062-4f9a-a4ff-9f7597be6268
# ╠═1e34dfa5-b9ca-479c-8a86-fdd064dfe408
# ╟─8947f891-b3c5-4472-b312-dba04face3d5
# ╠═0fa60678-66c6-442a-8ab0-63c3060dbd6e
# ╟─1a9e0ba4-52b6-4dc0-844b-88949e69befc
# ╠═61959657-80cd-4b73-b57b-4dc11a9fd6d2
# ╟─70b8a56e-6377-44bc-b5a4-af2a006d3dd3
# ╠═5cb47905-7379-4453-900b-92bf98a94ad1
# ╟─c807bf81-afc7-46b0-b073-26a27da9a55c
# ╠═3b11cd98-0317-4a9f-87af-7195c71e4fbc
# ╟─9fa0030c-6099-4fad-9d5b-800d5e197cf5
# ╟─6d13fdf7-7214-483b-9de3-9fd7accab4a1
# ╠═6a87b974-4ce6-45b0-b7e8-ce007c45659c
# ╟─abad5c96-060d-4a37-8da0-6d6b97dd1926
# ╟─68a2f5cf-3e0b-4dbc-ad66-ab728cbe7984
# ╠═0a447fbc-c31d-4eac-b272-65818efb62d7
# ╟─56ba3f8d-2c63-4b04-a389-37ff62e125d5
# ╠═1fd0aa0c-94d0-465d-a213-f788fa9572b7
# ╟─269b83aa-e7d2-4b04-b141-40997b6ba2a9
# ╟─251ac174-9ed8-41dc-be7c-db919e29515e
# ╠═893e77ec-6a36-47f2-9c84-6e723fe1a542
# ╠═905f0f37-5443-434d-899d-a78e46bb1344
# ╠═0f9b59cc-a8c6-4443-9080-aae4d8cd967b
# ╠═ff722e23-a476-4149-832f-ca875112a302
# ╟─2a5e302e-2c44-4fdc-bd9f-cbfef0df5763
# ╠═c0335f76-d167-4444-9bcc-27f0dc825db2
