### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 768f7b5e-f6df-11f0-a8ab-2920af629627
# ╠═╡ show_logs = false
using Pkg; Pkg.activate("..")

# ╔═╡ 0ede7402-6b92-456c-aaf9-89e1b272d63d
using Turing, StatsPlots

# ╔═╡ 59ce4733-40b6-4102-892c-13b7f46d9d11
using PlutoUI; TableOfContents()

# ╔═╡ 866cdda4-1863-4a8b-95fc-2e65849f6346
md"# Sampling notebook #1: Intro"

# ╔═╡ e4b9d0c5-776a-4fa0-ad95-3ab632e22609
md"""
This notebook will guide you through the basics of sampling in Julia.

To start off, let's load the required packages.
"""

# ╔═╡ 139b8552-bb44-4547-9c05-f5d435683e31
md"## Problem description"

# ╔═╡ 1785baee-b26b-4b2a-a3de-38efea3eb7b5
md"""
Most students of our beautiful campus Coupure had to endure the bitter task of cycling through the rain. What's worse, even when the rain has stopped, one is not safe from the wetness: get too close behind a fellow cyclist and their back wheel will pelt your face with dirty water droplets. In order to prevent this tragic fate, we must find how far we must stay away from the cyclist in front of us to keep our faces safe. We will do this by simulating the trajectories of the water droplets and **estimating the probability that they will coincide with our face**.
"""

# ╔═╡ 2adece60-5f0d-408f-82eb-51bf047a6d0f
md"""
![Problem illustration](https://i.imgur.com/7TDkD08.png)
> **Figure 1**: An illustration of the problem.
"""

# ╔═╡ e1ad094e-d10d-4cca-9b83-f653e02f4d58
md"## The model"

# ╔═╡ 034660af-152f-4850-a1d8-7e624ba3b8da
md"""
We consider a simple model for a droplet's motion: it is launched into the air by the back wheel giving it an initial velocity, and then falls down due to gravity. We neglect other factors such as wind resistance.
"""

# ╔═╡ 09875352-5f8c-4c89-a978-b704bf1ad6c4
md"""
After simulating the droplet's motion, we need to check whether it hits our face. We'll assume it's $1.5$m to $1.7$m above the ground and at some horizontal distance $x_f$ away from the cyclist in front of us. We then get hit if the droplet's altitude is in the range [1.5, 1.7] when it has travelled a distance $x_f$.
"""

# ╔═╡ 03f4935b-ac2b-4f73-9cee-ea04f700df83
md"### Mathematical description"

# ╔═╡ 543cea22-f592-41ef-a6d2-7d35063c387f
md"""
We can describe the horizontal position $x$ and vertical position $y$ through time $t$ using simple equations from [the physics of projectile motion](https://en.wikipedia.org/wiki/Projectile_motion#Trajectory_in_vacuum) (see also **Figure 2** for an illustration):

```math
\begin{align}
x(t) &= \mathrm{cos}(α) \, v_0 \, t \, ,
\\ y(t) &= \mathrm{sin}(α) \, v_0 \, t - \frac{g}{2} \, t^2 \, ,
\end{align}
```

where the parameters are:
-  $v_0$ is the initial velocity.
-  $\alpha$ is the launching angle.
-  $g = 9.8$ is the gravitational acceleration.
"""

# ╔═╡ ee9cd05c-7fa5-4412-b3b4-32b0a1774c87
md"""
Checking for droplet collision then comes down to finding the time $t_f$ the droplet has travelled a distance $x_f$
```math
t_f = \frac{x_f}{\mathrm{cos}(\alpha \, v_0)} \, ,
```
evaluating it's y-coordinate there
```math
y_d = \mathrm{sin}(α) \, v_0 \, t_f - \frac{g}{2} \, t_f^2 \, ,
```
and checking if it's in the range of our face
```math
y_d ∈ [1.5, 1.7] \, .
```
"""

# ╔═╡ af5a6f7d-b017-4cb6-b54a-a0defbfe7cdf
md"""
![Model illustration](https://i.imgur.com/fMyFGn4.png)
> **Figure 2**: An illustration of the droplet launching model. Note: one vector represents a velocity while the other represents an acceleration!
"""

# ╔═╡ e8821ac2-baa3-4674-868c-3a5b8593cc95
md"### The variables"

# ╔═╡ 4811e1af-9ba4-49b3-a3a1-5d8ff8d5af92
md"""
To reiterate, we have three variables in our model:
1. The initial droplet velocity $v_0$.
1. The droplet launching angle $\alpha$.
1. The horizontal position of our face $x_f$.

The first two of these are **stochastic variables**: they don't have a single, exact value, but can rather take a range of different values following some underlying distribution. This can be either because we simply don't know the real value, or there really are random fluctuations on the value of the variable. For our model, we need to specify what these distributions are. We will assume for now we magically know the distributions of these variables, though they are usually inferred from data as we will see next practical.
"""

# ╔═╡ 01156ef7-b194-4372-8bc8-de9fad5e12f1
md"For the initial droplet velocity $v_0$, we will assume a normal distribution with mean 5 [m/s] and a standard deviation of 1:"

# ╔═╡ f30290bd-9d31-4939-b5d6-f89f066775d2
plot(
	Normal(5, 1),
	xlabel = "v0 [m/s]", # set x label
	ylabel = "Probability density", # set y label
	legend = false, # remove legend from figure
	title = "Initial velocity prior" # set title
)

# ╔═╡ dab23eba-3eb1-434b-aef5-e08730e35114
md"For the launching angle $\alpha$, we will assume a symmetric **triangular distribution** going from $0$ radians to $\pi/2$ radians:"

# ╔═╡ 299e402c-634a-4b18-b867-489bafaf4564
plot(
	TriangularDist(0, pi/2), 
	xlabel = "α [rad]", # set x label
	ylabel = "Probability density", # set y label
	legend = false, # remove legend from figure
	title = "Launching angle prior" # set title
)

# ╔═╡ 314ba0ec-015f-488d-ac63-780b1dbf5bc2
md"The latter is an **input variable**: when cycling, we can choose how far away we stay from the cyclist in front of us. The x-position of our face $x_f$ is therefore an input to the model."

# ╔═╡ c070f3c9-dd80-49f9-b279-c77f46316116
md"## Defining the model"

# ╔═╡ 84a08354-7bf3-4783-9285-a7e8850e798b
md"""
Turing models are defined as Julia functions preceded by the `@model` macro.
Inside of them, you can define random variables with the `var ~ Distribution(params)` syntax, aside from doing the usual programming stuff. 
"""

# ╔═╡ 8f94567c-b87b-4eb4-8441-c6434da393a8
md"First, we define the constant parameter `g` outside of our model, so we can use it later on for plotting:" 

# ╔═╡ 07c7645a-2b56-4089-90c9-736d7794a109
const g = 9.8

# ╔═╡ 43582e68-1289-4d5a-8895-dc06ec2b9b5b
md"""
!!! note
	Using a global variable (a variable not defined inside a function) inside of a function, such as our Turing model, will slow down your code. We can fix this by declaring that the value of this global variable is constant with `const`. But don't worry: this is no programming course, so it's no problem if you forget this on your exam.
"""

# ╔═╡ 40f4edce-13c4-4063-98f4-23fd948c58dd
md"""
Our model can be defined as follows:
"""

# ╔═╡ f29fd729-0f05-49b5-9280-cf88edf02faf
# note: x_f is an input variable, so we add it as an input to the model
@model function splash(x_f) 
	# define the initial velocity
	v0 ~ Normal(5, 1)
	# define the launching angle
	α ~ TriangularDist(0, pi/2)
		
	# get the time it takes the droplet to travel the distance x_f
	t_f = x_f / (cos(α) * v0)
	# get the droplet's altitude at this time
	y_s = sin(α)*v0*t_f - g/2 * t_f^2
	# to respect the ground, set altitude to 0 if it would be negative
	y_s = max(0, y_s)
	
	return y_s 
end

# ╔═╡ 22144610-0c37-4501-b14c-83895d8b0eae
md"This function does not yet know the value of the input variable $x_f$, and is  therefore an incomplete model. Calling it with a value for all its inputs will return a complete Turing model. Let's set it to 1 for now."

# ╔═╡ 92d8ccc8-cd5e-4d0b-a8bf-4be301784a18
splash_model = splash(1)

# ╔═╡ 6e1d6cfb-47ed-4c0f-bc5f-c327c58ae128
md"""
!!! note
	Some models don't have any input variables. In this case, you still need to call it to create a "completed" model, though you will simply call it with no input arguments.
"""

# ╔═╡ 4e2bf07a-5f7c-4f33-b7e7-b481a6f8cef7
md"## Sampling the model"

# ╔═╡ cda45ee9-612e-4b2e-aca7-5da34a23dd02
md"""
Now that we have defined our model, we can use it to answer our question. We will estimate the probability a droplet hits our face from a distance $x_f$ by sampling a large number of different trajectories and checking in which fraction the droplet hits us.
"""

# ╔═╡ cbad0f13-aa5e-4e74-bf51-482437ff7ad3
md"### Generating samples"

# ╔═╡ c00cfd87-e68b-4c13-a4a4-2b4fdded49f8
md"We sample different sets of values for the model's parameters using the `sample` function:"

# ╔═╡ aea98d63-9277-4c07-9db5-9c94d940ef84
chain = sample(
	splash_model, # the Turing model we sample
	Prior(), # the sampling algorithm - we will always use `Prior()` this practical
	500 # the amount of samples
);

# ╔═╡ 7e67c4bc-e146-484a-8284-50127c3fbea1
chain

# ╔═╡ f648fd3a-057d-4e73-9331-e87c60c014a7
md"""
!!! note
	The `logprior` column gives the log prior probability of a sample. For example, for a sample with $v_0 = 6.1$ and $\alpha = 0.9$, this column will equal $\mathrm{log}(P(v_0 = 6.1) + P(\alpha = 0.9))$, where $P(v_0 = 6.1)$ equals the probability density of its prior distribution, a Normal(5, 1), at $x = 6.1$ (and analogous for $\alpha$).
"""

# ╔═╡ 1efb24f2-b988-437b-a7b7-2b3b7746dc07
md"### Getting the model output"

# ╔═╡ 46f17c07-9aec-4fee-8590-4eeb5a13b93e
md"To answer our question, we are interested in the `y_s` variable returned by the model. We can extract the output variable from our chain using the `generated_quantities` function:"

# ╔═╡ 071fc405-7a33-4dab-a69a-1e68a0fb35ad
y_s_samples = generated_quantities(splash_model, chain)

# ╔═╡ 20e99145-e80b-4259-9bca-f7e14c2d355e
md"It's often useful to visualize the distribution of our variable of interest. We can use the `histogram` function for this:"

# ╔═╡ 45d44cf5-2289-4947-8d4e-1dcddb5447de
histogram(y_s_samples, bins = 15)

# ╔═╡ 8d179748-40aa-4e88-b277-9234fe938251
md"Checking whether the droplet altitude matches our face's can be done with a simple inequality operation:"

# ╔═╡ 633a2b4c-379e-4e18-bbcd-b450930196a5
face_hit_samples = [1.5 <= y_s <= 1.7 for y_s in y_s_samples]

# ╔═╡ ad8a6ea2-30ba-4104-824f-12a7ea8e9718
md"Finally, we can estimate the probability of our face being hit by taking the mean amount of times we got hit in our sample."

# ╔═╡ e1f9d1e9-8a0c-4b28-ab67-a89e3287b50f
mean(face_hit_samples)

# ╔═╡ 1c9e2de1-21a0-4d4d-918d-9e2f0fc4d912
md"""And that's the answer to our question! If we're 1m behind the cyclist in front of us, we can expect to get hit by about ~1% of the droplets they launch (exact number depends on the random samples). You can manually change the value of $x_f$ to investigate for other distances, or scroll down to the [Essentials](#Essentials) section for a handy slider."""

# ╔═╡ bf5642f4-09d3-491d-8eb9-47fcdb67c977
md"""
!!! note
	Do you like your code sleek? You can also calculate this probability from your samples with just one line: `mean(1.5 .<= y_s_samples .<= 1.7)`
"""

# ╔═╡ 3067af6e-df69-48a0-a024-a6e219dcfc7c
md"### Getting the stochastic variables"

# ╔═╡ cbc04c5f-aa93-4340-a61e-1583d20d86b1
md"Often, it's enough to get a sample of the variable of interest we `return` in our Turing model. Sometimes, however, it can be useful to also get sampled values for the stochastic variables used to calculate the variable of interest. In this case, we could make a nice plot of the droplet's trajectories if we have all the sampled values of $v_0$ and $\alpha$."

# ╔═╡ d36784b8-9f3b-4434-88fe-03b75b94754e
md"Getting the values of a model's stochastic variable can be done simply by indexing the chain with the name of that variable as a `String` or `Symbol` (whichever you prefer):"

# ╔═╡ dad51c0a-b93f-4637-b095-244fd5513838
begin 
	v0_samples = chain["v0"]
	α_samples = chain[:α]
end

# ╔═╡ 9093b274-d329-4e4c-a374-d3bc625fa2a6
md"""
To plot the trajectories, we want an equation for $y$ in function of $x$. We can do this by substituting $x(t)$ into the equation for $y(t)$ to acquire the following:
"""

# ╔═╡ 0f5fbf8e-3e65-4212-9068-1d9feb09a484
md"""
```math
y(x) = \frac{\mathrm{sin}(α)}{\mathrm{cos}(α)} \, x - \frac{g}{2} \, \left( \frac{x}{\mathrm{cos}(α) \, v_0} \right)^2
```
"""

# ╔═╡ c9df5202-37f1-4e93-ab94-546ac7bb9444
md"""
!!! note
	You can also use this equation to find $y(x_f)$ in the Turing model.
"""

# ╔═╡ 4faf9589-a815-4223-92f9-34264fb22401
md"First we define our trajectory function in function of $x$ and our stochastic variables."

# ╔═╡ 577df19a-7da2-40b1-b1d4-38c3aec16ae7
y(x, α, v0) = sin(α)/cos(α)*x - g/2*(x / (cos(α) * v0))^2

# ╔═╡ 18f6c707-4369-4360-94a7-ee4a7f4e57d2
md"Then we get the specific trajectory for every sample of our chain by filling in the stochastic variables with the values we sampled."

# ╔═╡ 5d0bee7a-8fd0-44ff-adb5-6068477aaf21
trajectories = [
	# note: every trajectory in our vector is a function of x!
	# so we will write an anonymous function x -> function(x,...,...) 
	# to plot the trajectory in function of x
	x -> y(x, α_samples[i], v0_samples[i])
	for i in 1:length(v0_samples)
];

# ╔═╡ 00572e1a-d1e6-4c84-b1b2-b53404087dd8
md"Finally, we plot all sampled trajectories."

# ╔═╡ f550ed7a-1bba-40b5-9338-cdeb0520aae9
begin
	plot(trajectories, xlims = (0, 3), ylims = (0, 2), label = false,
		 color = :blue, alpha = 0.3
	)
	plot!([1, 1], [1.5, 1.7], linewidth = 5, color = :orange, label = "face")
end

# ╔═╡ 3b465fc0-3a9d-4e67-8d4a-a14b3cd61f4e
md"""
!!! note
	We set a lot of keyword arguments of our plot to make the figure pretty here, but don't panic: this is still no programming course, so we will always clearly provide which one(s) you need if you need to plot something. If you ever do want to look for one yourself, you can find them on the function's help page (`?plot` or the `🔍 Live docs` in the bottom right corner).
"""

# ╔═╡ befb8f05-6ab7-4594-9cf1-71df47c1df03
md"## Essentials"

# ╔═╡ 52af5de0-b516-4947-84f7-dde63e4c513e
md"""
The most essential code for the first practical is reiterated here without long explanations to provide an easy reference for making the practical exercises. Additionally, $x_f$ has been assigned to a slider so you can more easily explore at what distance your face is safe.

Side note: the code is wrapped in a let block so Pluto won't complain about the same variable names being used again.
"""

# ╔═╡ 2ceede7e-8f0c-420c-bde2-f16b726cd204
md"""
!!! note
	If you enjoy playing with the $x_f$ slider, the code will run much faster if you comment out the problem definition in the code below (select the code and press `CTRL+/`). This will make it reuse the (same) model defined previously and saves the time to run the model generation.
"""

# ╔═╡ fbabbc25-6924-4584-9236-3707ab03e859
@bind x_f Slider(0:0.1:5, default = 1.0, show_value = true)

# ╔═╡ 9a7926ac-615a-422f-b306-e21c56b8ae1c
# ╠═╡ show_logs = false
let
	# model definition (comment out for speed)
	g = 9.8
	@model function splash(x_f)
		v0 ~ Normal(5, 1) 
		α ~ TriangularDist(0, pi/2)
	
		t_f = x_f / (cos(α) * v0)
		y_s = sin(α)*v0*t_f - g/2 * t_f^2
		y_s = max(0, y_s)
		return y_s 
	end

	# get samples
	splash_model = splash(x_f)
	chain = sample(splash_model, Prior(), 500);

	y_s_samples = generated_quantities(splash_model, chain)
	v0_samples = chain["v0"]
	α_samples = chain[:α]

	# calculate interesting things
	prob_face_hit = mean([1.5 <= y_s <= 1.7 for y_s in y_s_samples])
	trajectories = [
		x -> sin(α_samples[i]) / cos(α_samples[i]) * x -
		g / 2 * ( x / (cos(α_samples[i])*v0_samples[i]) )^2
		for i in 1:length(v0_samples)
	]

	plot(trajectories, xlims = (0, 5), ylims = (0, 2), 
		 legend = false, color = :blue, alpha = 0.3,
		 title = "Probability of getting hit ≈ $(round(100*prob_face_hit, digits = 3))%"
	)
	annotate!(x_f, 1.6, text("☺", 30))
end

# ╔═╡ Cell order:
# ╟─866cdda4-1863-4a8b-95fc-2e65849f6346
# ╟─e4b9d0c5-776a-4fa0-ad95-3ab632e22609
# ╠═768f7b5e-f6df-11f0-a8ab-2920af629627
# ╠═0ede7402-6b92-456c-aaf9-89e1b272d63d
# ╠═59ce4733-40b6-4102-892c-13b7f46d9d11
# ╟─139b8552-bb44-4547-9c05-f5d435683e31
# ╟─1785baee-b26b-4b2a-a3de-38efea3eb7b5
# ╟─2adece60-5f0d-408f-82eb-51bf047a6d0f
# ╟─e1ad094e-d10d-4cca-9b83-f653e02f4d58
# ╟─034660af-152f-4850-a1d8-7e624ba3b8da
# ╟─09875352-5f8c-4c89-a978-b704bf1ad6c4
# ╟─03f4935b-ac2b-4f73-9cee-ea04f700df83
# ╟─543cea22-f592-41ef-a6d2-7d35063c387f
# ╟─ee9cd05c-7fa5-4412-b3b4-32b0a1774c87
# ╟─af5a6f7d-b017-4cb6-b54a-a0defbfe7cdf
# ╟─e8821ac2-baa3-4674-868c-3a5b8593cc95
# ╟─4811e1af-9ba4-49b3-a3a1-5d8ff8d5af92
# ╟─01156ef7-b194-4372-8bc8-de9fad5e12f1
# ╟─f30290bd-9d31-4939-b5d6-f89f066775d2
# ╟─dab23eba-3eb1-434b-aef5-e08730e35114
# ╟─299e402c-634a-4b18-b867-489bafaf4564
# ╟─314ba0ec-015f-488d-ac63-780b1dbf5bc2
# ╟─c070f3c9-dd80-49f9-b279-c77f46316116
# ╟─84a08354-7bf3-4783-9285-a7e8850e798b
# ╟─8f94567c-b87b-4eb4-8441-c6434da393a8
# ╠═07c7645a-2b56-4089-90c9-736d7794a109
# ╟─43582e68-1289-4d5a-8895-dc06ec2b9b5b
# ╟─40f4edce-13c4-4063-98f4-23fd948c58dd
# ╠═f29fd729-0f05-49b5-9280-cf88edf02faf
# ╟─22144610-0c37-4501-b14c-83895d8b0eae
# ╠═92d8ccc8-cd5e-4d0b-a8bf-4be301784a18
# ╟─6e1d6cfb-47ed-4c0f-bc5f-c327c58ae128
# ╟─4e2bf07a-5f7c-4f33-b7e7-b481a6f8cef7
# ╟─cda45ee9-612e-4b2e-aca7-5da34a23dd02
# ╟─cbad0f13-aa5e-4e74-bf51-482437ff7ad3
# ╟─c00cfd87-e68b-4c13-a4a4-2b4fdded49f8
# ╠═aea98d63-9277-4c07-9db5-9c94d940ef84
# ╠═7e67c4bc-e146-484a-8284-50127c3fbea1
# ╟─f648fd3a-057d-4e73-9331-e87c60c014a7
# ╟─1efb24f2-b988-437b-a7b7-2b3b7746dc07
# ╟─46f17c07-9aec-4fee-8590-4eeb5a13b93e
# ╠═071fc405-7a33-4dab-a69a-1e68a0fb35ad
# ╟─20e99145-e80b-4259-9bca-f7e14c2d355e
# ╠═45d44cf5-2289-4947-8d4e-1dcddb5447de
# ╟─8d179748-40aa-4e88-b277-9234fe938251
# ╠═633a2b4c-379e-4e18-bbcd-b450930196a5
# ╟─ad8a6ea2-30ba-4104-824f-12a7ea8e9718
# ╠═e1f9d1e9-8a0c-4b28-ab67-a89e3287b50f
# ╟─1c9e2de1-21a0-4d4d-918d-9e2f0fc4d912
# ╟─bf5642f4-09d3-491d-8eb9-47fcdb67c977
# ╟─3067af6e-df69-48a0-a024-a6e219dcfc7c
# ╟─cbc04c5f-aa93-4340-a61e-1583d20d86b1
# ╟─d36784b8-9f3b-4434-88fe-03b75b94754e
# ╠═dad51c0a-b93f-4637-b095-244fd5513838
# ╟─9093b274-d329-4e4c-a374-d3bc625fa2a6
# ╟─0f5fbf8e-3e65-4212-9068-1d9feb09a484
# ╟─c9df5202-37f1-4e93-ab94-546ac7bb9444
# ╟─4faf9589-a815-4223-92f9-34264fb22401
# ╠═577df19a-7da2-40b1-b1d4-38c3aec16ae7
# ╟─18f6c707-4369-4360-94a7-ee4a7f4e57d2
# ╠═5d0bee7a-8fd0-44ff-adb5-6068477aaf21
# ╟─00572e1a-d1e6-4c84-b1b2-b53404087dd8
# ╠═f550ed7a-1bba-40b5-9338-cdeb0520aae9
# ╟─3b465fc0-3a9d-4e67-8d4a-a14b3cd61f4e
# ╟─befb8f05-6ab7-4594-9cf1-71df47c1df03
# ╟─52af5de0-b516-4947-84f7-dde63e4c513e
# ╟─2ceede7e-8f0c-420c-bde2-f16b726cd204
# ╠═fbabbc25-6924-4584-9236-3707ab03e859
# ╠═9a7926ac-615a-422f-b306-e21c56b8ae1c
