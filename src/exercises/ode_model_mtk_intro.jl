### A Pluto.jl notebook ###
# v0.20.13

#> [frontmatter]
#> order = "1"
#> title = "1. ODE_model_MTK_intro"
#> tags = ["exercises"]
#> layout = "layout.jlhtml"
#> description = "Introduction to ModelingToolkit"

using Markdown
using InteractiveUtils

# ╔═╡ 1c46508e-2354-11f0-08db-d373e43929bd
using Pkg; Pkg.activate("..")

# ╔═╡ 2d125aad-334c-4210-a921-1acacefc5cfa
using StatsPlots, PlutoUI; TableOfContents()
# packages needed for example to make plots

# ╔═╡ d0f32979-01f6-4f07-acb3-627077c7c029
using OrdinaryDiffEq, ModelingToolkit
# packages needed to solve differential equations with ModelingToolkit

# ╔═╡ 43c4f402-45a9-48d5-b1c6-8a6ea8107ce2
using ModelingToolkit: t_nounits as t, D_nounits as D
# t will be the symbol for the time
# D will be the differentiation operator

# ╔═╡ f998d6c0-0791-4484-b088-a810d6e89a40
md"""
# Introduction to ModelingToolkit (ODE)
"""

# ╔═╡ 8f48f7a2-c9d7-4293-9966-4d52cc4e6b1e
TableOfContents() # displays that table of contents on the right

# ╔═╡ ee985cd8-86d1-462a-b4f1-3c25be215b57
md"""
**ModelingToolkit.jl** is a Julia package for **symbolic and numeric modeling of complex systems**, such as those described by for example **differential equations** and **algebraic equations**. It’s part of the **SciML ecosystem** (Scientific Machine Learning).

At its core, it lets you:

1. **Define models symbolically** – You describe variables, parameters, and equations symbolically, rather than writing them out procedurally.
2. **Automatically generate efficient code** – The toolkit simplifies, differentiates, and compiles your model into optimized Julia code for simulation or analysis.
3. **Compose models** – You can build complex systems by connecting smaller subsystems (useful in fields like biology, mechanics, or electrical circuits).
4. **Perform advanced analyses** – It supports symbolic simplification, Jacobian and sensitivity computation, structural analysis, and automatic differentiation.

In short, **ModelingToolkit** bridges the gap between symbolic math and high-performance simulation — enabling you to describe models declaratively, and then solve them efficiently using Julia’s differential equation solvers.

This is a barebones tutorial to ModelingToolkit, helping you to build simple models based on ODEs. See the [documentations](https://docs.sciml.ai/ModelingToolkit/stable/) for more details.

Below we will illustrate the use of ModelingToolkit for three different cases.
"""

# ╔═╡ f4db6481-80c9-4d56-b5a4-50e620d95143
md"""
## Case 1: Conical tank

![](https://i.ibb.co/7xYxhD7f/conical-tank.png)
"""

# ╔═╡ 5c1b77bb-5ef3-4434-a66d-89c86ca30216
md"""
As the figure here above illustrates, a conical tank with height $H$ and base radius $R$ has a constant inlet flow rate $q_{in}=q$. At the bottom, the outlet flow rate is proportional (proportionality factor: $k$) to the height $h$ of the liquid in the tank. At the liquid's surface level in the tank its radius is denoted as $r$. Hence, the liquid's volume is $V = \frac{1}{3}r^2 h \pi$. It also hold that $r=\frac{Rh}{H}$. The rate of change in $V$ is given by the following differential equation:
``` math
\cfrac{d V(t)}{dt} = q_{in} - q_{out} = q - k h(t)
```
We are not going to substitute $r$ in the equation for the volume and/or simplify $\frac{dV(t)}{dt}$. Instead, we are going to define the relevant variables that change over time and define their relation to eachother.
"""

# ╔═╡ e12df3dd-f896-4713-9328-551c949c34b1
md"""
### Defining variables and parameters
"""

# ╔═╡ 2c3e60ef-f738-43b4-bf14-0d723a389f09
md"""
We define variables for the height $h$, the radius $r$ and the volume $V$ as follow. Mind that you need to explicitely write their dependency on the time $t$. Here we have also provided a default initial value (initial condition) for $V$. *This is optional as later when creating the ODE problem you can overwrite these values*. The initial values for $h$ and $r$ will be calculated when solving the problem.
"""

# ╔═╡ 41af0644-4611-4880-950d-0b6ff2c58a15
@variables h(t) r(t) V(t)=0.1

# ╔═╡ ea09ad52-79f6-43c1-9389-e8efe2c73e21
md"""
We define the parameters. Optionally, as with the variables, you can also provide default values for the parameters and/or overwrite them later when creating the ODE problem. Here we have provided default parameter values for $H$ and $R$.
"""

# ╔═╡ 38602497-c591-4e22-af60-f9622096067f
@parameters k q H=5.0 R=2.0

# ╔═╡ c6f8bea8-948e-42f8-9730-d36d0fc3adb0
md"""
### Defining the equations
"""

# ╔═╡ 33d63418-855c-4ce9-98d9-df8e0184e752
md"""
Now we will define all necessary equations.
"""

# ╔═╡ 638e2962-115d-4b95-947e-ace3f2495bc7
md"""
We define the equation that relates $r$ to $h$. The equation is named `eq_radius`. You can choose an arbitrary name but use one that makes sense.  Mind that between the lefthandside and the righthandside a `~` (tilde) is used and not a `=` (equal sign). The `=` is only used to *name* the equation. Furthermore, write the name of the variables without explicitely stating their dependency on the time.
"""

# ╔═╡ f98ecb14-9df2-4ba3-a1c7-3ab75e3e8055
eq_radius = r ~ h * R / H

# ╔═╡ 3c4ac6b0-08d2-4f93-bd98-b52a03e9af7d
md"""
Next, we define the equation that relates $V$ to $r$ and $h$.
"""

# ╔═╡ c18b1610-eadf-474e-8da9-9c1f21845bc1
eq_volume = V ~ π * r^2 * h / 3

# ╔═╡ f6cb6648-95bd-4a8e-85ad-626b97bb6cb4
md"""
Finally, we will define the differential equation that relates the rate of change in $V$ to the inlet flow rate and the outlet flow rate. The symbolic operator `D` takes the symbolic differentiation of `V` in `D(V)`.
"""

# ╔═╡ e87ea53c-8629-4f49-af3b-1c753940e27d
change_vol = D(V) ~ q - h * k   # change in volume = inflow - outflow

# ╔═╡ 74fc9ede-c032-4d5a-b508-c6de4021f460
md"""
After all equations have been defined, we will bundle them in a vector using square brackets as below.
"""

# ╔═╡ 8e6efab4-ac91-40bd-9805-88e0dabb2cd1
eqns_tank = [eq_radius, eq_volume, change_vol]

# ╔═╡ 168f79ba-ca24-45af-9ac4-dfa10d48a3d7
md"""
### Building the MTK system and creating the ODE problem
"""

# ╔═╡ 9d302a6c-3b1b-47b8-85df-a658baaca050
md"""
We will build an ODE system for MTK (ModelingToolkit). You need to provide the list of equations (cf. here: `eqns_tank`) and the symbol `t` for the time to the function `ODESystem`. Optionally you can also provide so-called *continuous* and *discrete* events, but this will be illustrated later. The name of the MTK model (or system) given is here `tank`.
"""

# ╔═╡ c6bfc46c-7ca5-4d9e-88fb-3e7065ec6bb0
@mtkbuild tank = ODESystem(eqns_tank, t)

# ╔═╡ c3fce45c-f4d6-49e7-82cf-79a58412bc57
md"""
After having built the MTK model we need to create the ODE problem. You need to provide the MTK model (cf. here: `tanks`), a vector with the initial conditions, a time span and a vector of parameter values to the function `ODEProblem`. Because $V$, $r$ and $h$ are related to each other, we provide the initial condition for $V$ only and provide guessed values for $r$ and $h$.

Some remarks:
- Since the initial value for $V$ was set when defining the variables, we could as well have written `[]` instead of `[V=>0.1]` if you want to use the default value. You can overwrite the initial value here if you wish.
- Beware that $V$, $r$ and $h$ are related with eachother and that only one of them must be given an initial value. The other two values need to be 'guessed' (cf. the keyword argument `guesses=...`).
- We have intentionally not specified the parameters `R` and `H` in the vector of parameters, because we want to use their default values. You can overwrite the default parameter values here if you wish.
"""

# ╔═╡ d6a153d4-d41f-47d4-88b0-b10ac2e22b04
tank_prob = ODEProblem(tank, [V=>0.1],
	(0.0, 100.0), [q=>0.1, k=>0.05],
	guesses=[r=>0.1, h=>0.1])

# ╔═╡ 816fda45-635b-4354-9437-c2efd7374782
md"""
### Solving and plotting
"""

# ╔═╡ f6e0db08-723d-4b55-ac40-6828105ed6a7
md"""
We will now solve the ODE problem using the function `solve`. The function expects two inputs: an `ODEProblem` and an ODE solver. It is also possible to call `solve` with just the `ODEProblem`, in which case it will try to automatically select a suitable solver for the ODE problem. This function also has many keyword arguments, which you can find in the docs page of the function (type `?solve` in a cell or write in directly in the `🔍 Live docs` in the bottom right corner). One we will use often is `saveat`, which decides at what timesteps the solution will be returned. For example, if you provide the option `saveat=1`, then the solution will be approximated every time unit.
"""

# ╔═╡ 9492a458-2e7c-4828-85a3-20ed266e794c
tank_sol = solve(tank_prob)
# tank_sol = solve(tank_prob, saveat=1)

# ╔═╡ 61e22259-4a15-44cc-ab2b-e21030502135
md"""
If you want to check the number of time instances at which the solution was approximated, you can use the function `length`:
"""

# ╔═╡ 66947904-6b80-4279-ae2b-a3c20886062d
length(tank_sol)

# ╔═╡ baf2cc98-0a10-4c30-8d32-66b93b4c37aa
md"""
If you plan on further analyzing the solution for one of the variables, you can access their values by indexing the solution object with the symbolic variable. For example:
"""

# ╔═╡ 258ad571-6e2d-470f-ae5f-5c5ac34afbab
tank_sol[V]

# ╔═╡ b5ef6331-0239-4c4c-a7e2-8df1d1077b07
md"You can also use the variable name as a `Symbol` (created by adding a leading `:`):"

# ╔═╡ 3240e248-06ad-4720-b0d0-14665a0ee128
tank_sol[V] == tank_sol[:V]

# ╔═╡ 9bbaf346-0195-405d-bf3a-8360d49f03ab
md"""
You can plot the evolution of $V$, $h$ and $r$ using the function `plot` and providing the name of the solution object. The keyword argument `idxs` allows you to choose what variables to plot.
"""

# ╔═╡ 99158cd8-753c-4d6b-9288-a78dfc1fa86a
plot(tank_sol, idxs=[V, r, h])
# plot(tank_sol)

# ╔═╡ ce573c77-ff84-446e-9460-1aff19b324f2
md"""
You can calculate the steady state (equilibrium) values by solving a steady state problem in the way below. The first argument to `SteadyStateProblem` is the name of the MTK system/model, the second argument consists of inital guesses of the steady state values and the third argument consists of the parameter values.
"""

# ╔═╡ 08a72b39-81c8-4b3a-90a7-761bf6398e80
equil_val_tank = solve(SteadyStateProblem(tank, [V=>1.0, h=>2.0, r=>0.5], [q=>0.1, k=>0.05]))

# ╔═╡ 24087790-80c0-41fb-af35-1b49078705e7
md"""
Then you can retrieve the steady state value of each variable in the following way:
"""

# ╔═╡ f3afa935-f194-42f0-989c-48593f8ef8eb
equil_val_tank[V]

# ╔═╡ fb554cd3-8974-4cb2-9ed1-71fdfb52911a
equil_val_tank[h]

# ╔═╡ 12fc3604-1871-4197-bf95-5674dc9d9afb
md"""
You can also provide an expression to `idxs` consisting of variables and/or parameters. Below, we show the ratio $V/h$.
"""

# ╔═╡ 8c9e806b-ee8b-4c84-b9bb-113e1d8ab4e4
plot(tank_sol, idxs=V/h)

# ╔═╡ b108c4a6-ee94-4773-8b88-fe02b70addfb
md"""
## Case 2: Mass-and-spring system

![](https://i.ibb.co/MDTGngfq/damped-spring.png)
"""

# ╔═╡ abda4552-c4d2-4b30-92bf-c17c5bfd1535
md"""
In this case we will simulate an under-damped harmonic oscillation. In an under-damped harmonic oscillator, the system experiences oscillatory motion while the amplitude gradually decreases over time due to damping. The deviation $y(t)$ describes the displacement of the oscillator from its equilibrium position as a function of time, typically following an exponentially decaying sinusoidal form. The total mechanical energy $E(t)$, composed of both kinetic and potential contributions, also decays exponentially as energy is dissipated by the damping mechanism. Studying $y(t)$ and $E(t)$ provides insight into how the oscillation amplitude and system energy evolve under the influence of damping.

The motion of an under-damped harmonic oscillator is governed by Newton’s second law, leading to the differential equation
``` math
m y''(t) = -k y(t) - \mu y'(t)
```
where $m$ is the mass, $\mu$ the damping coefficient, and $k$ the spring constant. The under-damped condition corresponds to $\mu^2 < 4mk$, resulting in oscillatory motion with an exponentially decaying amplitude. We will make sure that the parameter values in this example meets this condition.

The total mechanical energy is given by
``` math
E(t) = \cfrac{1}{2} m y'(t)^2 + \cfrac{1}{2} k y(t)^2
```
and, due to damping, it decreases over time as energy is continuously dissipated by the resistive force.
"""

# ╔═╡ ab2c91b9-56e8-4935-91a2-06519c2dcf0a
md"""
### Defining variables and parameters
"""

# ╔═╡ 497fe740-385e-4bd1-b5c6-1203b00eb1ab
md"""
We define variables for the position $y$ and the total energy $E$ as follow. Don't forget to mention the dependency on the time $t$. Here, you could provide a default initial condition for $y$ but this is optional. $E$ doesn't need an initial condition because $E$ will be calculated at each iteration step when $y$ and $y'$ were computed.
"""

# ╔═╡ a0f5eeee-3b5f-4588-bd31-1cff6c993f6b
@variables y(t) E(t)  # position and energy

# ╔═╡ b1277bf9-9958-4c07-aba1-77a9e5252b20
md"""
We define the parameters. Optionally, as with the variables, you can also provide default values for the parameters and/or overwrite them later when creating the ODE problem. An other option is to provide an informative *description* string with the meaning of each of the parameters.
"""

# ╔═╡ 26921be4-4438-4be7-b6c0-d9ea28af0865
md"""
!!! note
	We can't define a parameter `k` here because we already defined one in the previous section, and Pluto doesn't allow variables to be defined multiple times in different places (this guarantees the code returns the same output no matter what order it is run in). Because of this, we call it `k_s` instead.
"""

# ╔═╡ d3954593-63ae-4dfe-b752-01339b44adcd
@parameters μ [description="friction"] m [description="mass"] k_s [description="spring constant"]

# ╔═╡ f67ffde7-622a-45bf-b987-5802a7413610
md"""
### Defining the equations
"""

# ╔═╡ b195f1cb-c87b-407b-9840-a2e3c3da4fba
md"""
In the differential equation (Newton's second law) we use the second derivative of $y$: $y''=\cfrac{d^2y}{dt^2}$. You can use the following notations to write the symbolic second derivative of $y$:
"""

# ╔═╡ b4c29986-4d70-4551-bef3-4b79f032c3ce
D(D(y))  # second-order derivative wrt time

# ╔═╡ 268ab13e-2f15-4b1e-90b7-cbbcf5496163
(D^2)(y)  # second-order derivative wrt time

# ╔═╡ b9a3e9c0-577b-45ae-82bd-79d33879343f
md"""
We will define the second-order differential equation for the position (i.e. the equation based on Newton’s second law).
"""

# ╔═╡ ed48e70a-a608-46a1-9339-0bb3a93a4ea1
# spring_eq = m * D(D(y)) ~ - k_s * y - μ * D(y)
eq_spring = m * (D^2)(y) ~ - k_s * y - μ * D(y)

# ╔═╡ caada693-9aca-41ac-bd37-bc3fee641125
md"""
Next, we will define the expression for the total energy of the system.
"""

# ╔═╡ dc7e26c6-2c88-475d-a4ab-64dbb4f9a956
exp_energy = E ~ m*D(y)^2/2 + y^2*k_s/2  # we can add quantities to keep track off. 

# ╔═╡ ecc6d921-ae18-4a25-b2f8-3636ae8fc85c
md"""
### Building the MTK system and creating the ODE problem
"""

# ╔═╡ 3e7eeaec-ee16-4078-829e-61eb9acba0e8
md"""
We will now build an ODE system for MTK (ModelingToolkit). You will need to provide the vector of equations and the symbol `t` for the time to the function `ODESystem`. In addition we will introduce a discrete event: at time t=`40`, the mass is reduced by 75%. You could image that part of the mass fell from the spring during the movement. The name of the MTK model (or system) given is here `spring_ode`.
If you have multiple discrete events you can include them in the following way:
`[[...]=>[...~...], [...]=>[...~...], ...]`.
"""

# ╔═╡ 75a3db3f-818a-4cb8-948a-dfbcd0582a45
@mtkbuild spring_ode = ODESystem([eq_spring, exp_energy], t; discrete_events=[[40]=>[m~(1-0.75)*m]])

# ╔═╡ 2bcd3e8f-1915-4247-81a4-f7803c55e252
md"""
As you can notice, the second order ODE has been converted into a system of two first order ODEs. A new variable `y_t` was introduced such that `y_t` $= \cfrac{dy}{dt}$. This new variable is nothing else than $y'$.
"""

# ╔═╡ 9746751b-144f-4216-8916-d64dce6ce747
md"""
After having built the MTK model we need to create the ODE problem. You need to provide the name of the MTK model (cf. here: `spring_ode`), a vector with the initial conditions, a time span and a vector of parameter values to the function `ODEProblem`. Remark that since the second order ODE was converted into a system of two ODEs, you need to provide two initial conditions: one for $y$ and one for $y'$. The symbolic notation for $y'$ is `D(y)` and the latter should be used while creating the ODE problem.
"""

# ╔═╡ 30b3ed77-0709-4410-8d56-d576d69b5e0b
spring_prob = ODEProblem(spring_ode, [y=>2.0, D(y)=>-1.0], (0.0, 100.), [m=>3, k_s=>0.6, μ=>1e-1])

# ╔═╡ 53913925-13d2-4a70-b310-8c19ae628581
md"""
### Solving and plotting
"""

# ╔═╡ 85d111eb-5452-43de-a674-a6242e7dc5ca
md"""
We will now solve the ODE problem using the function `solve`. In this case we have provided a solver (cf. `Tsit5()`) and a relative tolerance to be met by the solver (cf. `reltol=1e-9`). The solver `Tsit5()` is a recommended solver for non-stiff problems. See the [documentation](https://docs.sciml.ai/DiffEqDocs/dev/solvers/ode_solve/) for more details.
"""

# ╔═╡ 15ba5472-72da-4ae4-9438-fe8a1d88d9c0
sol_spring = solve(deepcopy(spring_prob), Tsit5(), reltol=1e-9)

# ╔═╡ cf253b87-8a15-4e4e-b372-f427793143bb
md"""
!!! note
	When working with events, it is good practice to take a `deepcopy` of the ODE problem before solving it. This is because events can change parameter values of the problem while solving, which do not reset when solving is done. Therefore solving the same problem a second time would use a changed set of initial parameters, and give different results. Copying the problem before solving prevents this issue.
"""

# ╔═╡ 4b04ccc3-94c4-4331-a781-d02be313ef73
md"""
You can plot the evolution of $y$ and $y'$ using the function `plot` by just providing the name of the solution object.
"""

# ╔═╡ d29f022c-924b-4d26-ab8c-65f11eb85b6a
plot(sol_spring)

# ╔═╡ ab4aca36-5d8a-43b5-9e7a-c9c6f8288fef
md"""
You can clearly notice some change in the oscillatory at t=40.
"""

# ╔═╡ f954d1ff-d9dd-4ba2-8dc2-a5ccd87e1ce3
md"""
If you want to see the evolution of $E$, you can provide it to the `idxs` keyword argument.
"""

# ╔═╡ ff3645ab-12d1-46ac-b4d8-ce431ba27491
plot(sol_spring, idxs=E)

# ╔═╡ 27670f6a-9f2c-4768-a1ea-986f9ba5301c
md"""

## Case 3: Lotka-Volterra

Classical model for prey-predator relations
"""

# ╔═╡ 718d6e46-7e5e-4f08-a3c6-f14c9bd64cb2
md"""
The **Lotka–Volterra model**, also known as the **predator–prey model**, is a pair of first-order, nonlinear differential equations that describe the dynamic interaction between two biological species: one as a prey population and the other as its predator. Developed independently by Alfred J. Lotka and Vito Volterra in the 1920s, the model captures the cyclical nature of population sizes — with predator numbers rising and falling in response to changes in prey abundance, and vice versa. Despite its simplicity, the Lotka–Volterra framework remains a cornerstone of theoretical ecology, providing insight into population oscillations, stability, and the balance of ecosystems.

In this case we will consider a population of rabbits ($R$) as prey and foxes ($F$) as predators. Suppose that their evolution over time is governed by the following equations:
``` math
\begin{align}
\cfrac{dR}{dt} &= \alpha \left(1 - \cfrac{R}{K}\right) R - \beta R F \\
\cfrac{dF}{dt} &= \gamma R F - \delta F
\end{align}
```
The rabbit population grows at rate with $\alpha$ but is limited by a carrying capacity $K$, while predation by foxes reduces it at a rate $\beta$. The fox population increases proportionally to the number of hunts ($\gamma R F$) and declines naturally at rate $\delta$.
"""

# ╔═╡ 8644c5c2-efe3-42b1-86f9-ef3c3e37a4cd
md"""
### Defining variables and parameters
"""

# ╔═╡ 4537da11-be37-42d2-b4d8-5e08272e9f65
md"""
We will define the variables for this model. For the rabbits we will use the symbol 🐰 and for the foxes the symbol 🦊. In order to get the first symbol, type a back slash `\` and then type `:rabbit:` followed by the TAB-key. The second symbol you can get analogously, type `\` and then type `:fox_face:` followed by the TAB-key. If you hit the TAB-key before finishing the name of the symbol, you can see different options. Don't forget to specify the dependence on the time `t`.
"""

# ╔═╡ 3100925e-d5dc-4d56-af84-ffab2059a4c5
@variables 🐰(t) 🦊(t)

# ╔═╡ 454c13d8-6938-4a92-b553-db83d47793c8
md"""
Next, we will define the parameters. To get the greek letters is similar to getting the rabbit and foxes symbols. For example, to get α, type `\` and then `alpha` followed by the TAB-key. The other ones you can get with `beta`, `gamma` and `delta`.
"""

# ╔═╡ 298c0f1b-4237-4c17-993b-a2969257790f
@parameters α β γ δ K

# ╔═╡ 62266eb6-a6b1-404a-926f-870925e9f429
md"""
### Defining the equations
"""

# ╔═╡ 087eb724-7f8a-4c99-be7a-5924a449cd29
md"""
Now we will define the equations. We will readily define and bundle both equations simultanuously and name the vector of equations `LV_eqs`.
"""

# ╔═╡ 130f2711-9ac8-491f-b580-28fe0889257a
LV_eqs = [
	D(🐰) ~ α * 🐰 * (1 - 🐰 / K) - β * 🐰 * 🦊,
	D(🦊) ~ γ *  🐰 * 🦊 - δ * 🦊 
]

# ╔═╡ 8181f5ca-3774-4cf2-bcf3-24d86774b88e
md"""
### Building the MTK system and creating the ODE problem
"""

# ╔═╡ 6d10c03d-0795-4458-88a6-f7484efe9478
md"""
In this example we will introduce a so-called continuous event. The continuous event is formulated here below. It states that when the population of rabbits hits 300, then its population is brought back to 100. You can imagine that if you have enough rabbits, say 300, hunters will hunt the rabbits and bring its population back to 100. In our case there is only one continuous event, but it is possible to include multiple continuous events in the following way:
`[[...~...]=>[...~...], [...~...]=>[...~...], ...]`
"""

# ╔═╡ 1a920b39-00f6-4455-b226-b6961af5c6d1
rabbitmanagment = [[🐰 ~ 300.0] => [🐰 ~ 100.0]]

# ╔═╡ 3f2170d3-be36-415e-8c65-4f492baccac0
md"""
We will now build an ODE system for MTK. You will need to provide the list of equations, the symbol `t` for the time and the continuous event to the function `ODESystem`.
"""

# ╔═╡ 1cdf417d-daea-4684-89bd-7deb951c6b3a
@mtkbuild lv_sys = ODESystem(LV_eqs, t; continuous_events=rabbitmanagment)

# ╔═╡ b7299cae-799b-4840-8501-e935dca06183
md"""
Next, we will create the ODE problem by providing the name of the MTK system, the initial conditions for 🐰 and 🦊, the time span and the parameter values.
"""

# ╔═╡ 08bfcddc-703f-4ef4-a44e-6386eb91de76
# lv_prob = ODEProblem(lv_sys, [🐰=>1.0, 🦊=>1e-2], (0, 100), [α=>0.6, β=>0.6, γ=>0.04, δ=>0.5, K=>1000])
lv_prob = ODEProblem(lv_sys, [🐰=>100, 🦊=>4], (0.0, 100.0), [α=>0.6, β=>0.016, γ=>0.004, δ=>0.6, K=>1000])

# ╔═╡ f5d834ca-e79a-4f96-bebf-049171c625a4
md"""
### Solving and plotting
"""

# ╔═╡ 1eebace8-9332-4555-b6a9-1f386059816b
md"""
We will now solve the ODE problem using the function `solve`. In this case we have provided a solver (cf. `Tsit5()`) and a relative tolerance to be met by the solver (cf. `reltol=1e-9`). Also notice the `deepcopy` of the ODE problem when dealing with events.
"""

# ╔═╡ 5c927b82-1133-4508-996a-7ed7111cc4d5
sol_LV = solve(deepcopy(lv_prob), Tsit5(), reltol=1e-9)

# ╔═╡ d57502bc-46e2-4268-a8cc-09ccb8f7c21d
md"""
You can plot the evolution of 🐰 and 🦊 using the function `plot` by just providing the name of the solution object. You can clearly see that when 🐰 hits 300, its population is brought back to 100. This happens three times. In this period of time the population of 🦊 has grown in a way that the 🐰 cannot reach a population of 300 anymore. Instead, both populations go to steady state values.

We have provided written labels this time because the fancy symbols don't come through in the legend of the plot.
"""

# ╔═╡ 9db55fe4-9f67-4d23-85ec-fbafe21f09a7
plot(sol_LV, label=["rabbits" "foxes"])

# ╔═╡ a218fb57-a64b-47ca-b24a-df1a549aaa9c
md"""
Here below we calculate the steady state values for the rabbits and the foxes.
"""

# ╔═╡ a6de7f90-5349-44db-b681-dc35433344e2
equil_val_LV = solve(SteadyStateProblem(lv_sys, [🐰=>100, 🦊=>40], [α=>0.6, β=>0.016, γ=>0.004, δ=>0.6, K=>1000]))

# ╔═╡ 055e5609-ab93-483a-a494-d3587802c92b
md"""
The steady state value of the rabbits:
"""

# ╔═╡ 3b3d890f-7648-4d8e-abfb-bb1ef386bbea
equil_val_LV[🐰]

# ╔═╡ 68b6c3ea-9093-49bb-8e4f-b4f9544b1ad6
md"""
The steady state value of the foxes:
"""

# ╔═╡ 43a0eb8c-af3d-4cc3-bc5a-18a4ce121a40
equil_val_LV[🦊]

# ╔═╡ 8d7f21e3-6512-4fc2-a290-82b7bbeefa7b
md"""
You can also plot the ratio of the rabbits over the foxes in the following way.
"""

# ╔═╡ cc65899f-f16c-4194-b2d8-bcfe3dcbd8c8
plot(sol_LV, idxs=🐰/🦊, label="rabbits/foxes")

# ╔═╡ Cell order:
# ╟─f998d6c0-0791-4484-b088-a810d6e89a40
# ╠═1c46508e-2354-11f0-08db-d373e43929bd
# ╠═2d125aad-334c-4210-a921-1acacefc5cfa
# ╠═d0f32979-01f6-4f07-acb3-627077c7c029
# ╠═43c4f402-45a9-48d5-b1c6-8a6ea8107ce2
# ╠═8f48f7a2-c9d7-4293-9966-4d52cc4e6b1e
# ╟─ee985cd8-86d1-462a-b4f1-3c25be215b57
# ╟─f4db6481-80c9-4d56-b5a4-50e620d95143
# ╟─5c1b77bb-5ef3-4434-a66d-89c86ca30216
# ╟─e12df3dd-f896-4713-9328-551c949c34b1
# ╟─2c3e60ef-f738-43b4-bf14-0d723a389f09
# ╠═41af0644-4611-4880-950d-0b6ff2c58a15
# ╟─ea09ad52-79f6-43c1-9389-e8efe2c73e21
# ╠═38602497-c591-4e22-af60-f9622096067f
# ╟─c6f8bea8-948e-42f8-9730-d36d0fc3adb0
# ╟─33d63418-855c-4ce9-98d9-df8e0184e752
# ╟─638e2962-115d-4b95-947e-ace3f2495bc7
# ╠═f98ecb14-9df2-4ba3-a1c7-3ab75e3e8055
# ╟─3c4ac6b0-08d2-4f93-bd98-b52a03e9af7d
# ╠═c18b1610-eadf-474e-8da9-9c1f21845bc1
# ╟─f6cb6648-95bd-4a8e-85ad-626b97bb6cb4
# ╠═e87ea53c-8629-4f49-af3b-1c753940e27d
# ╟─74fc9ede-c032-4d5a-b508-c6de4021f460
# ╠═8e6efab4-ac91-40bd-9805-88e0dabb2cd1
# ╟─168f79ba-ca24-45af-9ac4-dfa10d48a3d7
# ╟─9d302a6c-3b1b-47b8-85df-a658baaca050
# ╠═c6bfc46c-7ca5-4d9e-88fb-3e7065ec6bb0
# ╟─c3fce45c-f4d6-49e7-82cf-79a58412bc57
# ╠═d6a153d4-d41f-47d4-88b0-b10ac2e22b04
# ╟─816fda45-635b-4354-9437-c2efd7374782
# ╟─f6e0db08-723d-4b55-ac40-6828105ed6a7
# ╠═9492a458-2e7c-4828-85a3-20ed266e794c
# ╟─61e22259-4a15-44cc-ab2b-e21030502135
# ╠═66947904-6b80-4279-ae2b-a3c20886062d
# ╟─baf2cc98-0a10-4c30-8d32-66b93b4c37aa
# ╠═258ad571-6e2d-470f-ae5f-5c5ac34afbab
# ╟─b5ef6331-0239-4c4c-a7e2-8df1d1077b07
# ╠═3240e248-06ad-4720-b0d0-14665a0ee128
# ╟─9bbaf346-0195-405d-bf3a-8360d49f03ab
# ╠═99158cd8-753c-4d6b-9288-a78dfc1fa86a
# ╟─ce573c77-ff84-446e-9460-1aff19b324f2
# ╠═08a72b39-81c8-4b3a-90a7-761bf6398e80
# ╟─24087790-80c0-41fb-af35-1b49078705e7
# ╠═f3afa935-f194-42f0-989c-48593f8ef8eb
# ╠═fb554cd3-8974-4cb2-9ed1-71fdfb52911a
# ╟─12fc3604-1871-4197-bf95-5674dc9d9afb
# ╠═8c9e806b-ee8b-4c84-b9bb-113e1d8ab4e4
# ╟─b108c4a6-ee94-4773-8b88-fe02b70addfb
# ╟─abda4552-c4d2-4b30-92bf-c17c5bfd1535
# ╟─ab2c91b9-56e8-4935-91a2-06519c2dcf0a
# ╟─497fe740-385e-4bd1-b5c6-1203b00eb1ab
# ╠═a0f5eeee-3b5f-4588-bd31-1cff6c993f6b
# ╟─b1277bf9-9958-4c07-aba1-77a9e5252b20
# ╟─26921be4-4438-4be7-b6c0-d9ea28af0865
# ╠═d3954593-63ae-4dfe-b752-01339b44adcd
# ╟─f67ffde7-622a-45bf-b987-5802a7413610
# ╟─b195f1cb-c87b-407b-9840-a2e3c3da4fba
# ╠═b4c29986-4d70-4551-bef3-4b79f032c3ce
# ╠═268ab13e-2f15-4b1e-90b7-cbbcf5496163
# ╟─b9a3e9c0-577b-45ae-82bd-79d33879343f
# ╠═ed48e70a-a608-46a1-9339-0bb3a93a4ea1
# ╟─caada693-9aca-41ac-bd37-bc3fee641125
# ╠═dc7e26c6-2c88-475d-a4ab-64dbb4f9a956
# ╟─ecc6d921-ae18-4a25-b2f8-3636ae8fc85c
# ╟─3e7eeaec-ee16-4078-829e-61eb9acba0e8
# ╠═75a3db3f-818a-4cb8-948a-dfbcd0582a45
# ╟─2bcd3e8f-1915-4247-81a4-f7803c55e252
# ╟─9746751b-144f-4216-8916-d64dce6ce747
# ╠═30b3ed77-0709-4410-8d56-d576d69b5e0b
# ╟─53913925-13d2-4a70-b310-8c19ae628581
# ╟─85d111eb-5452-43de-a674-a6242e7dc5ca
# ╠═15ba5472-72da-4ae4-9438-fe8a1d88d9c0
# ╟─cf253b87-8a15-4e4e-b372-f427793143bb
# ╟─4b04ccc3-94c4-4331-a781-d02be313ef73
# ╠═d29f022c-924b-4d26-ab8c-65f11eb85b6a
# ╟─ab4aca36-5d8a-43b5-9e7a-c9c6f8288fef
# ╟─f954d1ff-d9dd-4ba2-8dc2-a5ccd87e1ce3
# ╠═ff3645ab-12d1-46ac-b4d8-ce431ba27491
# ╟─27670f6a-9f2c-4768-a1ea-986f9ba5301c
# ╟─718d6e46-7e5e-4f08-a3c6-f14c9bd64cb2
# ╟─8644c5c2-efe3-42b1-86f9-ef3c3e37a4cd
# ╟─4537da11-be37-42d2-b4d8-5e08272e9f65
# ╠═3100925e-d5dc-4d56-af84-ffab2059a4c5
# ╟─454c13d8-6938-4a92-b553-db83d47793c8
# ╠═298c0f1b-4237-4c17-993b-a2969257790f
# ╟─62266eb6-a6b1-404a-926f-870925e9f429
# ╟─087eb724-7f8a-4c99-be7a-5924a449cd29
# ╠═130f2711-9ac8-491f-b580-28fe0889257a
# ╟─8181f5ca-3774-4cf2-bcf3-24d86774b88e
# ╟─6d10c03d-0795-4458-88a6-f7484efe9478
# ╠═1a920b39-00f6-4455-b226-b6961af5c6d1
# ╟─3f2170d3-be36-415e-8c65-4f492baccac0
# ╠═1cdf417d-daea-4684-89bd-7deb951c6b3a
# ╟─b7299cae-799b-4840-8501-e935dca06183
# ╠═08bfcddc-703f-4ef4-a44e-6386eb91de76
# ╟─f5d834ca-e79a-4f96-bebf-049171c625a4
# ╟─1eebace8-9332-4555-b6a9-1f386059816b
# ╠═5c927b82-1133-4508-996a-7ed7111cc4d5
# ╟─d57502bc-46e2-4268-a8cc-09ccb8f7c21d
# ╠═9db55fe4-9f67-4d23-85ec-fbafe21f09a7
# ╟─a218fb57-a64b-47ca-b24a-df1a549aaa9c
# ╠═a6de7f90-5349-44db-b681-dc35433344e2
# ╟─055e5609-ab93-483a-a494-d3587802c92b
# ╠═3b3d890f-7648-4d8e-abfb-bb1ef386bbea
# ╟─68b6c3ea-9093-49bb-8e4f-b4f9544b1ad6
# ╠═43a0eb8c-af3d-4cc3-bc5a-18a4ce121a40
# ╟─8d7f21e3-6512-4fc2-a290-82b7bbeefa7b
# ╠═cc65899f-f16c-4194-b2d8-bcfe3dcbd8c8
