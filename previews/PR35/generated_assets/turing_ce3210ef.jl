### A Pluto.jl notebook ###
# v0.20.4

#> [frontmatter]
#> order = "4"
#> title = "Turing Cheat Sheet"
#> date = "2025-01-29"
#> tags = ["cheat sheets"]
#> description = "Turing Cheat Sheet"
#> layout = "layout.jlhtml"
#> 
#>     [[frontmatter.author]]
#>     name = "Bram Spanoghe"
#>     [[frontmatter.author]]
#>     name = "Michiel Stock"

using Markdown
using InteractiveUtils

# ╔═╡ 1cd55133-8bb2-4883-953d-77824c3df54f
using Pkg; Pkg.activate("../../pluto-deployment-environment")

# ╔═╡ 13ff3e4f-dae2-4e08-946a-857f672e1049
using Distributions

# ╔═╡ c5eb9fd3-e518-4f96-9581-c0032864c7ec
using Turing

# ╔═╡ f20a633b-9e9f-40f8-8a1a-eec1890a5652
using StatsPlots

# ╔═╡ c11a0059-c0a5-4789-90b1-41b7d7c46d1c
md"""
!!! important
	When running this notebook locally, deactivate or delete the above cell."""

# ╔═╡ 387c26b2-ff21-4d85-b6eb-0cc0c8ee7347
md"# `Turing` cheatsheet"

# ╔═╡ 157d654f-1a17-4eeb-ac0d-cdb6ceac7729
md"## Distributions.jl"

# ╔═╡ 7833de91-b850-41bd-96b4-4bd706473343
distr = LogNormal(2.0, 1.0) # Define a LogNormal distribution with mean 2.0 and standard deviation 1.0

# ╔═╡ e6394f7a-cefc-4734-854f-72440380973e
md"### Basic statistics"

# ╔═╡ e6e377a8-edee-4ae6-8f1b-d30243ef139c
mean(distr) # Calculate the mean of the distribution

# ╔═╡ 914f4cbb-3594-4c1c-955c-f50fa2116476
var(distr) # Calculate the variance

# ╔═╡ f426160a-4872-41b1-9676-9c20f8e56528
std(distr) # Calculate the standard deviation

# ╔═╡ b70257c0-bc4c-4090-b5c5-bb0e97b42479
quantile(distr, [0.25, 0.5, 0.75]) # Calculate the quartiles

# ╔═╡ a2f3f4b5-ccfa-4ea6-ac98-67473cef00b0
md"### Evaluate probability density and cumulative probability"

# ╔═╡ 99efa065-6d95-487c-a17b-8cbefd2ba34e
pdf(distr, 2.0) # Probability density at x = 2.0

# ╔═╡ 79409039-2e23-4278-8c42-49f019a5eadd
cdf(distr, 5.0) # Probability that a random variable is less than 5.0

# ╔═╡ a46067bd-c80e-46b5-b691-4fdb5ccb73f0
md"### Sampling random values"

# ╔═╡ 073779a6-9850-4a33-bc0a-961e701567ef
rand(distr) # Draw a single random sample

# ╔═╡ 187d3922-7f1a-46dc-a3c3-a60294484238
mysample = rand(distr, 1000) # Generate 1000 random samples

# ╔═╡ 460920bf-6bbd-40a3-817b-ea15c812305d
md"### Calculate statistics from samples"

# ╔═╡ e0c5c5f8-5468-4dd5-a0e7-845822657c2a
mean(mysample) # Approximate the mean using the sample

# ╔═╡ f567cde8-7c1d-480b-8d0b-0cf8570ac079
std(mysample) # Approximate the standard deviation

# ╔═╡ cf744a4c-bee2-4b0b-b731-073f632a402e
md"### Calculate probabilities using samples"

# ╔═╡ 2f130b63-0ea8-4acc-80e3-4c6c218fb7ff
mean(x -> x^2 > 5, mysample) # P(X^2 > 5): Method 1 - Anonymous function and mean

# ╔═╡ 4d86b8d1-1c18-4da0-a6ce-97bb2169b51a
mean(mysample.^2 .> 5) # P(X^2 > 5): Method 2 - Boolean operations

# ╔═╡ 6c8c2963-be4f-4fdc-bc27-3185fc84c68b
filtered_sample = filter(x -> x^2 > 5, mysample) # P(X^2 > 5): Method 3 - Filtering

# ╔═╡ 31d88376-619b-45e5-b216-9b78d2310e9d
length(filtered_sample) / length(mysample)

# ╔═╡ 4379606f-0e88-47e2-ae4c-de76d29b5e00
md"### Other calculations with samples"

# ╔═╡ 228dc51e-bd89-49d0-86bb-1b02c464aa49
mean(sin, mysample) # Approximate E[sin(X)] using the sample (more efficient)

# ╔═╡ 884435f8-d679-4337-998d-37ae1ef41ee0
mean(sin.(mysample)) # Same

# ╔═╡ d6fe2fbe-93e3-460c-9528-6fe01cbcc69a
md"## Turing.jl"

# ╔═╡ ba244c3f-f960-4af8-b755-11c8a636d053
@model function mymodel()
	x ~ Exponential(2.0) # Exponential prior for x
	y ~ Truncated(Normal(1., x), 0.0, 10.0) # Truncated Normal for y, dependent on x
	z ~ Poisson(y) # Poisson distribution for z, dependent on y
	return z^2 / y # computed result (optional)
end

# ╔═╡ 697b913c-8521-42df-81e5-1394c5fa8512
md"### Sampling"

# ╔═╡ 50e0fb1c-7c1e-441a-89c6-7d6311064487
xyzmodel = mymodel() # Build the sampling model

# ╔═╡ 9778ba2d-316b-4af9-b2d8-639f8d0f9b23
xyz = rand(xyzmodel) # Generate a single sample (x, y, z)

# ╔═╡ 4782d1ec-fec4-4988-9a85-5ae0776bf986
xyz[:x], xyz[:y], xyz[:z] # Extract the individual variables

# ╔═╡ f085a715-4271-410b-aa2b-668200cfa5bf
mysamples = [rand(xyzmodel) for i in 1:1000] # Generate 1000 samples

# ╔═╡ f82ec171-dd3e-4221-8898-5aca1a5ad673
xyzmodel() # random sample of the result (z^2 / y)

# ╔═╡ 108bb5fd-e2d7-4c18-b1b5-ba9dc1782819
md"### Calculate probabilities using samples:"

# ╔═╡ 2cf531dc-4eb3-420e-9bc0-c49b6444362c
mean(xyz -> xyz[:z] == 0, mysamples) 
	# P(Z=0): Method 1 - Anonymous function and mean

# ╔═╡ 068d7aab-0119-486c-9d0c-ca20716edfa2
length(filter(xyz -> xyz[:z] == 0, mysamples)) / length(mysamples) 
	# P(Z=0): Method 2 - Filtering

# ╔═╡ 8a560071-3fc4-483e-a793-59119fa519de
mean([rand(xyzmodel)[:z] == 0 for i in 1:1000]) 
	# P(Z=0): Method 3 - Boolean Operations on samples

# ╔═╡ f89e246e-4d6b-4973-981b-825126a5d0cc
mean(xyz -> xyz[:z] == 0, filter(xyz -> xyz[:x] > 1, mysamples)) 
	# P(Z=0 | x > 1): Method 1 - Filtering and mean

# ╔═╡ c0f7a9fa-923a-4c7a-84fc-07466ca25b4a
mean([xyz[:z] == 0 for xyz in mysamples if xyz[:x] > 1]) 
	# P(Z=0 | x > 1): Method 2 - Boolean operations on samples

# ╔═╡ 89b1794e-fff4-4961-9743-a15dda6e694d
md"### Inference"

# ╔═╡ 837b652b-6f3a-4cc0-92c4-add2b543f2b0
xyzmodel_cond = xyzmodel | (z=3.0,) # Condition the model on Z = 3

# ╔═╡ eb97d036-7b8b-471b-8b6f-27bfeda6ec57
logprior(xyzmodel_cond, (x=1.3, y=0.3))

# ╔═╡ 6179e4b6-d227-4d77-b210-1953db31dd9f
loglikelihood(xyzmodel_cond, (x=1.3, y=0.3))

# ╔═╡ a139ae59-a5ad-4c67-9dda-ca8922edd5d9
logjoint(xyzmodel_cond, (x=1.3, y=0.3)) # log-prior + log-likelihood

# ╔═╡ 1e76681a-0405-4f1d-a984-76a5339f5bed
chain = sample(xyzmodel_cond, NUTS(), 10_000) # Obtain samples from posterior

# ╔═╡ 8193a275-00ac-402a-875e-ba6e393cff74
summarize(chain) # Summarize the chain (means, quantiles, etc.)

# ╔═╡ a275fb02-94b7-46cd-acf2-2b56a4a5310e
quantile(chain) # Quantiles, default 2.5%, 25.0%, 50.0%, 75.0%, 97.5%

# ╔═╡ bef94063-6fe5-4cf9-91db-0bc2530fb2eb
generated_quantities(xyzmodel, chain) # generates the result (z^2 / y) based on the

# ╔═╡ daa390d6-2a1a-494e-80cc-fdf41234dbd4
md"### Plotting"

# ╔═╡ 82661a71-76bb-4f06-91a1-1d6f926556b2
plot(chain) # Create diagnostic plots of the chain (traceplot, etc.)

# ╔═╡ 69457dcb-cac1-4d03-b817-e13294dda8e2
chain_x = chain[:x] # Extract samples for ’x’

# ╔═╡ bb27747c-f3e2-48e5-9d56-00ab229adb6d
chain_y = chain[:y]

# ╔═╡ 1d6817da-d23e-40d6-8192-fc838ced9541
histogram(chain_x, title="Histogram of x | z=3") # Plot posterior of ’x’

# ╔═╡ b2999345-b85c-4edf-a156-8469b5ada1f3
md"### Calculations on posterior samples"

# ╔═╡ 11888ba7-c494-4d08-876e-898420d62a7c
mean(log, chain_x) # Approximate E[log(X) | Z=3]

# ╔═╡ 500dba58-49d9-4ab0-839b-95b7c41e3e2a
mean(chain_x .> chain_y) # Approximate P(X > Y | Z=3)

# ╔═╡ Cell order:
# ╠═1cd55133-8bb2-4883-953d-77824c3df54f
# ╟─c11a0059-c0a5-4789-90b1-41b7d7c46d1c
# ╟─387c26b2-ff21-4d85-b6eb-0cc0c8ee7347
# ╟─157d654f-1a17-4eeb-ac0d-cdb6ceac7729
# ╠═13ff3e4f-dae2-4e08-946a-857f672e1049
# ╠═7833de91-b850-41bd-96b4-4bd706473343
# ╟─e6394f7a-cefc-4734-854f-72440380973e
# ╠═e6e377a8-edee-4ae6-8f1b-d30243ef139c
# ╠═914f4cbb-3594-4c1c-955c-f50fa2116476
# ╠═f426160a-4872-41b1-9676-9c20f8e56528
# ╠═b70257c0-bc4c-4090-b5c5-bb0e97b42479
# ╟─a2f3f4b5-ccfa-4ea6-ac98-67473cef00b0
# ╠═99efa065-6d95-487c-a17b-8cbefd2ba34e
# ╠═79409039-2e23-4278-8c42-49f019a5eadd
# ╟─a46067bd-c80e-46b5-b691-4fdb5ccb73f0
# ╠═073779a6-9850-4a33-bc0a-961e701567ef
# ╠═187d3922-7f1a-46dc-a3c3-a60294484238
# ╟─460920bf-6bbd-40a3-817b-ea15c812305d
# ╠═e0c5c5f8-5468-4dd5-a0e7-845822657c2a
# ╠═f567cde8-7c1d-480b-8d0b-0cf8570ac079
# ╟─cf744a4c-bee2-4b0b-b731-073f632a402e
# ╠═2f130b63-0ea8-4acc-80e3-4c6c218fb7ff
# ╠═4d86b8d1-1c18-4da0-a6ce-97bb2169b51a
# ╠═6c8c2963-be4f-4fdc-bc27-3185fc84c68b
# ╠═31d88376-619b-45e5-b216-9b78d2310e9d
# ╟─4379606f-0e88-47e2-ae4c-de76d29b5e00
# ╠═228dc51e-bd89-49d0-86bb-1b02c464aa49
# ╠═884435f8-d679-4337-998d-37ae1ef41ee0
# ╟─d6fe2fbe-93e3-460c-9528-6fe01cbcc69a
# ╠═c5eb9fd3-e518-4f96-9581-c0032864c7ec
# ╠═ba244c3f-f960-4af8-b755-11c8a636d053
# ╟─697b913c-8521-42df-81e5-1394c5fa8512
# ╠═50e0fb1c-7c1e-441a-89c6-7d6311064487
# ╠═9778ba2d-316b-4af9-b2d8-639f8d0f9b23
# ╠═4782d1ec-fec4-4988-9a85-5ae0776bf986
# ╠═f085a715-4271-410b-aa2b-668200cfa5bf
# ╠═f82ec171-dd3e-4221-8898-5aca1a5ad673
# ╟─108bb5fd-e2d7-4c18-b1b5-ba9dc1782819
# ╠═2cf531dc-4eb3-420e-9bc0-c49b6444362c
# ╠═068d7aab-0119-486c-9d0c-ca20716edfa2
# ╠═8a560071-3fc4-483e-a793-59119fa519de
# ╠═f89e246e-4d6b-4973-981b-825126a5d0cc
# ╠═c0f7a9fa-923a-4c7a-84fc-07466ca25b4a
# ╟─89b1794e-fff4-4961-9743-a15dda6e694d
# ╠═837b652b-6f3a-4cc0-92c4-add2b543f2b0
# ╠═eb97d036-7b8b-471b-8b6f-27bfeda6ec57
# ╠═6179e4b6-d227-4d77-b210-1953db31dd9f
# ╠═a139ae59-a5ad-4c67-9dda-ca8922edd5d9
# ╠═1e76681a-0405-4f1d-a984-76a5339f5bed
# ╠═8193a275-00ac-402a-875e-ba6e393cff74
# ╠═a275fb02-94b7-46cd-acf2-2b56a4a5310e
# ╠═bef94063-6fe5-4cf9-91db-0bc2530fb2eb
# ╟─daa390d6-2a1a-494e-80cc-fdf41234dbd4
# ╠═f20a633b-9e9f-40f8-8a1a-eec1890a5652
# ╠═82661a71-76bb-4f06-91a1-1d6f926556b2
# ╠═69457dcb-cac1-4d03-b817-e13294dda8e2
# ╠═bb27747c-f3e2-48e5-9d56-00ab229adb6d
# ╠═1d6817da-d23e-40d6-8192-fc838ced9541
# ╟─b2999345-b85c-4edf-a156-8469b5ada1f3
# ╠═11888ba7-c494-4d08-876e-898420d62a7c
# ╠═500dba58-49d9-4ab0-839b-95b7c41e3e2a
