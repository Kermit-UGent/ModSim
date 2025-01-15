### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 886c7932-da4b-45cc-ba73-8047389e4895
using Pkg; Pkg.activate("..");

# ╔═╡ 80bc0e86-5ad3-4d61-9600-8dc05b86599d
using Turing, StatsPlots

# ╔═╡ f2327bc5-35ff-4d2d-b2dd-5a3964a417be
n_samples = 1000

# ╔═╡ 1ff166c0-b639-11ef-2143-ad11b6834d4c
md"# 1: Buffon's needles"

# ╔═╡ 5ec0a0ee-0c5f-437e-be85-162805af4731
md"""
A wise man once said: ["there is no greater joy than estimating π"](https://en.wikipedia.org/wiki/Approximations_of_%CF%80). Next to throwing darts at the unit square, another method to accomplish this is using [Buffon's needle problem](https://en.wikipedia.org/wiki/Buffon%27s_needle_problem).

The experiment is as follows: consider a floor with parallel lines all a distance of 1 away from eachother. Now drop a needle of length 1 (and width ~0) on the floor with a **random position and angle**. What is the probability $P_{cross}$ that the needle will cross one of the lines?
"""

# ╔═╡ 1d74c5dc-8cae-4c2d-9a17-d9a6c52a9169
md"The following image illustrates the problem (imagine $l$ = $t$ = 1) for two needles, where `a` crosses a line and `b` does not."

# ╔═╡ 4a67aba8-6015-4c95-9b5b-ca4fd2bed1bd
md"""
![Buffon's needles](https://upload.wikimedia.org/wikipedia/commons/thumb/5/58/Buffon_needle.svg/1920px-Buffon_needle.svg.png)
"""

# ╔═╡ 2382f49a-2e76-4e30-b1d0-eca298779d73
md"""
Using sampling magic, it's not difficult to make an estimate of this probability, $\hat{P}_{cross}$. Solving the problem analytically shows that the exact value is:
```math
P_{cross} = \frac{2}{\pi}
```

Therefore, our estimator for π is:
```math
π = \frac{2}{\hat{P}_{cross}}
```
"""

# ╔═╡ 2f277afb-eb92-44bf-b4fb-b994fefb7266
md"""
!!! question
	Estimate π using the Buffon's needle approximation.
"""

# ╔═╡ 870b3d09-de00-4b9c-a67f-a38a4cc289aa
@model function needles()
    x0 ~ Uniform(0, 1)
    theta ~ Uniform(0, pi)
    x_end = x0 + cos(theta)
    crossing = x_end < 0 || x_end > 1
    return crossing
end

# ╔═╡ af516e05-5f26-49f9-a1ca-d1f3d4f8ea99
needle_model = needles()

# ╔═╡ a50f848d-4d27-47e0-bfb6-76ef620660c0
sp_n = [needle_model() for i in 1:n_samples]

# ╔═╡ 4fb459cc-8449-4a3c-9164-7b89c3348c4e
mean(sp_n)

# ╔═╡ f2061e56-1de3-4be1-8c3b-5bcae649bf4d
pi_est = 2 / mean(sp_n)

# ╔═╡ c8941726-9e81-47fc-9b7e-cb3b5c0c61ca
md"# 2: Attraction"

# ╔═╡ 431023df-3724-4325-b0ac-96dbf5e4fd20
md"""
Following a course on electromagnetism will teach one that computing the net force between 2 arbitrary shapes can be a terrifying task. Tragedy has it then, that this is a very general problem with application from making fusion reactors to space travel. We can ease the pain by turning it into a sampling problem.

We'll start in a humble manner and simulate **the gravitational force between 2 cubes**. You can do this by randomly sampling a point from both cubes and using the [formula for gravitational force](https://en.wikipedia.org/wiki/Gravity?variant=zh-tw#Newton's_theory_of_gravitation). You can simplify the problem with $m_1 = m_2 = 1$ and [$G = 1$](https://www.reddit.com/media?url=https%3A%2F%2Fpreview.redd.it%2Fqnn7guyqtug41.jpg%3Fauto%3Dwebp%26s%3Da1cd045d4f4704b8439fb70585fad1b1d39c79f6).
"""

# ╔═╡ 058d2a08-7d40-46c5-8c5b-2b5ce9d2eb18
md"""
!!! questions
	- Given 2 cubes of size 1 with their centres a distance of 1.1 away from eachother, what is the estimated net force between the two? Is this the same as if you had treated the cubes as point masses?
	- How many samples do you need to estimate this force reliably? Define a reliable estimator as one having a standard deviation of 10% of the mean value. Visualise the distribution of the estimator.
"""

# ╔═╡ 06d0e92f-1f07-41fd-b6ee-e94eb539627d
@model function cubeforce(m1 = [0.0, 0.0, 0.0], m2 = [1.1, 0.0, 0.0])
	# defining m1 and m2 is not required but it looks nice
    x1 ~ Uniform(m1[1] - 1/2, m1[1] + 1/2)
    y1 ~ Uniform(m1[2] - 1/2, m1[2] + 1/2)
    z1 ~ Uniform(m1[3] - 1/2, m1[3] + 1/2)

	x2 ~ Uniform(m2[1] - 1/2, m2[1] + 1/2)
    y2 ~ Uniform(m2[2] - 1/2, m2[2] + 1/2)
    z2 ~ Uniform(m2[3] - 1/2, m2[3] + 1/2)

	r_squared = (x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2
	G = 1 # for simplicity
	F = G / r_squared
	
    return F
end

# ╔═╡ 43119060-bd35-4c4c-8831-d5c5df1d8dd5
cubemodel = cubeforce();

# ╔═╡ d5229ed1-8473-42e5-b1f9-98fb547d3e85
force_sp = [cubemodel() for _ in 1:10_000];

# ╔═╡ 420120c3-1c8c-46e4-a7c0-0f2405dd159a
pointmass_force = 1 / 1.1^2

# ╔═╡ 1ce8f808-b917-4831-89b4-4c318f05724d
mean(force_sp)

# ╔═╡ a324a615-7066-4dbb-9410-95e1e1c46ac4
histogram(force_sp)

# ╔═╡ e38e8f51-90db-4136-84e2-f06cd03d502a
required_samples = 150 # manually change until standard deviation is about 0.1

# ╔═╡ a7e975a9-beb5-4169-9d22-3e970be2838b
averaged_force_sp = [mean([cubemodel() for _ in 1:required_samples]) for _ in 1:1000];

# ╔═╡ 7ff5cb25-47bd-432f-a404-359abaf44e3a
std(averaged_force_sp)

# ╔═╡ 99c7f6d6-eada-491b-b966-fdb1195fc111
histogram(averaged_force_sp)

# ╔═╡ 616ca5c1-22ec-4ff9-97fb-2e49157eefd2
md"# 2⭐: Donut hole (bonus exercise)"

# ╔═╡ 957e1821-12c0-4d1a-86ac-7f24873a5860
md"""
The previous exercise promised the computation of forces between 2 arbitrary shapes. While cubes are a nice shape, they're not exactly exciting. This time, let's instead estimate the net force between a cube and a [torus](https://en.wikipedia.org/wiki/Torus), and use the estimation to animate the resulting movement of the cube.

The function `movecube` written below will make the required animation. The only thing it still lacks is the correct implementation of the function `donutforce`. This is a **Turing model** which does the following:
- It takes the centre of a cube as input, e.g. `[3.0, 0.0, 1.4]`
- It returns the force between a random point of this cube and that of a torus with centre [0.0, 0.0, 0.0], $R$ = 1, $r$ = 0.1, and rotated so a point moving along the x-axis does not touch it (see figure below)
"""

# ╔═╡ 74076669-ecd0-4b70-a001-cea5de0c8e91
@model function donutforce(m)
	r_big = 1.0
	r_small_max = 0.1

	theta_big ~ Uniform(0, 2*pi)
    theta_small ~ Uniform(0, 2*pi)
    r_small ~ Uniform(0, r_small_max/2)

    x1 = r_small * cos(theta_small)
    y1 = r_big * sin(theta_big) + r_small * sin(theta_small)
    z1 = r_big * cos(theta_big) + r_small * cos(theta_small)

	size = 1.0
	
	x2 ~ Uniform(m[1] - size/2, m[1] + size/2)
    y2 ~ Uniform(m[2] - size/2, m[2] + size/2)
    z2 ~ Uniform(m[3] - size/2, m[3] + size/2)

	r_squared = (x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2
	G = 1.0
	F = G / r_squared
	
    return F
end

# ╔═╡ 71915c43-fa93-47c8-8807-11a9ee5e9daa
"""
	movecube(num_samples, steps)

Animates the movement of a cube under influence of a torus centered around the origin.

# Inputs
- `num_points`: Number of samples that is used to estimate the net force.
- `steps`: Number of steps to simulate over

# Optinal inputs
- `m0`: The initial position of the cube
"""
function movecube(num_samples, steps; m0 = [1.0, 0.0, 0.0])
    mid = m0 # middle of cube
	vel = zeros(3) # velocity

	tm = zeros(3) # torus middle
	thetas = pi:0.1:3*pi
	ntp = length(thetas) # number of torus points
	
	torus_xs = zeros(ntp)
	torus_ys = sin.(thetas)
	torus_zs = cos.(thetas)

    anim = @animate for step in 1:steps

		# estimate force
        force = [donutforce(mid)() for _ in 1:num_samples] |>
			mean
		force_direction = - mid
		force_vec = force * force_direction

		# move cube according to force
        vel += force_vec # F ~ a, assuming equal mass
		mid += vel
		
		# plot

		scatter(torus_xs[(ntp÷2):end], torus_ys[(ntp÷2):end],
			torus_zs[(ntp÷2):end], color = :purple) # torus (back)
		scatter!([mid[1]], [mid[2]], [mid[3]], 
			xlim = [-5, 5], ylim = [-2, 2], zlim = [-2, 2],
			legend = false, color = :pink, shape = :square) # cube
		scatter!(torus_xs[1:(ntp÷2)], torus_ys[1:(ntp÷2)], 
			torus_zs[1:(ntp÷2)], color = :purple) # torus (front)
    end

	return gif(anim, fps = 12)
end

# ╔═╡ 70583a24-9b58-4262-8f8a-ee2e1af9560e
movecube(100, 50, m0 = [3.0, 0.0, 0.0])

# ╔═╡ Cell order:
# ╠═886c7932-da4b-45cc-ba73-8047389e4895
# ╠═80bc0e86-5ad3-4d61-9600-8dc05b86599d
# ╠═f2327bc5-35ff-4d2d-b2dd-5a3964a417be
# ╟─1ff166c0-b639-11ef-2143-ad11b6834d4c
# ╟─5ec0a0ee-0c5f-437e-be85-162805af4731
# ╟─1d74c5dc-8cae-4c2d-9a17-d9a6c52a9169
# ╟─4a67aba8-6015-4c95-9b5b-ca4fd2bed1bd
# ╟─2382f49a-2e76-4e30-b1d0-eca298779d73
# ╟─2f277afb-eb92-44bf-b4fb-b994fefb7266
# ╠═870b3d09-de00-4b9c-a67f-a38a4cc289aa
# ╠═af516e05-5f26-49f9-a1ca-d1f3d4f8ea99
# ╠═a50f848d-4d27-47e0-bfb6-76ef620660c0
# ╠═4fb459cc-8449-4a3c-9164-7b89c3348c4e
# ╠═f2061e56-1de3-4be1-8c3b-5bcae649bf4d
# ╟─c8941726-9e81-47fc-9b7e-cb3b5c0c61ca
# ╟─431023df-3724-4325-b0ac-96dbf5e4fd20
# ╟─058d2a08-7d40-46c5-8c5b-2b5ce9d2eb18
# ╠═06d0e92f-1f07-41fd-b6ee-e94eb539627d
# ╠═43119060-bd35-4c4c-8831-d5c5df1d8dd5
# ╠═d5229ed1-8473-42e5-b1f9-98fb547d3e85
# ╠═420120c3-1c8c-46e4-a7c0-0f2405dd159a
# ╠═1ce8f808-b917-4831-89b4-4c318f05724d
# ╠═a324a615-7066-4dbb-9410-95e1e1c46ac4
# ╠═e38e8f51-90db-4136-84e2-f06cd03d502a
# ╠═a7e975a9-beb5-4169-9d22-3e970be2838b
# ╠═7ff5cb25-47bd-432f-a404-359abaf44e3a
# ╠═99c7f6d6-eada-491b-b966-fdb1195fc111
# ╟─616ca5c1-22ec-4ff9-97fb-2e49157eefd2
# ╟─957e1821-12c0-4d1a-86ac-7f24873a5860
# ╟─71915c43-fa93-47c8-8807-11a9ee5e9daa
# ╠═74076669-ecd0-4b70-a001-cea5de0c8e91
# ╠═70583a24-9b58-4262-8f8a-ee2e1af9560e
