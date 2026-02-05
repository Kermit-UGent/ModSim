### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ a18638d0-138d-11f0-0571-33060f7da7ba
using Pkg; Pkg.activate("..")

# ╔═╡ fe265128-33a2-40d1-ab71-bc5c53979a54
using Turing, StatsPlots

# ╔═╡ c870d28d-7b2e-48f8-8990-b2cc0943cb09
using Optim, StatsBase

# ╔═╡ 2a73d23c-73fc-4845-b359-c6fe22077140
using PlutoUI

# ╔═╡ fb07436c-cc5d-4721-8a80-d7f7201721d7
md"# Model selection"

# ╔═╡ a397f48e-4228-435e-af13-c2bc71c8cb05
TableOfContents()

# ╔═╡ a6635715-52e3-44e7-9c00-0be751f830d6
md"## Who's that distribution?"

# ╔═╡ 473c003a-a34b-4e92-8a82-99450514d755
md"""
You decide to turn your life around and invest all your money into **clams**, or more specifically, **pearl farming**. Before setting up your full-scale farm, you decide to test the pearl-producing capabilities of different species of mollusk. You cultivate 10 different species, wait a year, and collect and measure the resulting pearls.

You want to compare the species by **fitting a distribution** to the pearl sizes. This way you can compare average size, expected deviation and the probability to get a really big pearl. **However, you don't know what distribution the pearl sizes follow.**
Since they're positive real numbers, 2 good candidates are the `Exponential` and `LogNormal` distributions.
"""

# ╔═╡ a2dec4e7-7c16-4ae5-974e-8c01ec31507a
md"""
!!! question
	For every molluks species, does the data follow an Exponential or a LogNormal distribution?
"""

# ╔═╡ c4e09491-fd0a-43e7-8e74-04832a045a48
md"""
![Picture of a black pearl in its shell](https://upload.wikimedia.org/wikipedia/commons/thumb/2/24/Black_pearl_and_his_shell.jpg/1280px-Black_pearl_and_his_shell.jpg)
> Source: Brocken Inaglory (Wikipedia)
"""

# ╔═╡ aba7efa9-a0c3-4da4-82a2-bb7efc6b9e22
md"### Data"

# ╔═╡ c8b56408-0710-4646-bb2b-e336e7f689c2
# ╠═╡ disabled = true
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	
function generate_point()
	firstdistr = rand() < 0.5
	if firstdistr
		medist = Exponential(rand(Uniform(0.1, 10)))
	else
		medist = LogNormal(rand(Uniform(0.1, log(10))), rand(Uniform(0.1, 1.0)))
	end
	n_samples = rand(Poisson(15))
	samples = rand(medist, n_samples) .|> x -> round(x, digits = 2)
	return samples
end

distr_data = [generate_point() for _ in 1:10];
	
end;
  ╠═╡ =#

# ╔═╡ 19dd461d-cc03-41c9-a3d1-fcd477d08eb0
distr_data = [[5.23, 2.79, 5.81, 4.36, 7.46, 4.46, 0.83, 6.45, 6.2, 6.53, 6.24, 8.72, 3.15], [1.12, 1.04, 0.09, 0.06, 0.67, 0.33, 0.41, 0.87, 1.23, 4.28, 7.46, 1.21, 0.19, 0.3, 0.59, 1.74, 0.66, 5.97, 0.3, 1.43, 1.11], [0.79, 3.37], [6.84, 11.28, 9.32, 6.27, 6.73, 10.28, 13.69, 8.32, 6.95], [0.48, 8.69, 3.92], [1.53, 1.83, 1.86, 0.87, 1.53, 2.51, 2.14, 1.82, 0.28, 3.57, 0.42, 1.67, 2.39, 4.18], [6.0, 2.37, 14.05, 4.01, 8.51, 5.29, 5.24, 18.01, 2.65, 8.91, 6.37, 2.54], [0.58, 2.41, 12.87, 14.67, 3.97, 13.8, 2.54, 4.7, 17.6, 18.3, 11.16, 0.81, 18.86, 2.3], [1.07, 0.6, 2.24, 0.02, 13.28, 4.88, 0.22, 18.54, 2.81, 2.97, 9.29, 2.98, 23.94, 0.39, 29.25, 1.05, 5.52, 0.39, 4.81, 3.73, 0.49], [8.39, 10.45, 1.93, 12.18, 3.26, 5.12, 8.3, 4.09, 20.41, 0.61, 18.31]];

# ╔═╡ 01a530db-5e18-4540-99a9-eecd9e61e1af
md"You can choose the mollusk species here and see the data for its pearl sizes."

# ╔═╡ 6d39fbde-81c5-4786-99e1-b587bd15f94c
md"Mollusk species"

# ╔═╡ fce554e2-a6da-4046-ab4e-9153e27aef7a
@bind distr_index Slider(1:10, show_value = true)

# ╔═╡ 7b1ce376-d796-4913-a488-ee921bd13855
pearlsizes = distr_data[distr_index]

# ╔═╡ 712756bc-0408-45a6-931b-eec07b6e052f
histogram(pearlsizes, bins = 0:ceil(maximum(pearlsizes)))

# ╔═╡ b4de93ea-a919-41c2-bb9c-bc1beccc9ea8
md"### Model definition"

# ╔═╡ 5a6f8d49-3841-46ff-a826-88ba0bed868b
md"""
We need to define a model for the two candidate distributions. The likelihood was already given above. For the priors, you can assume the following:
- Exponential model
  - μ ~ `Uniform(0, 10)`
- LogNormal model
  - μ ~ `Uniform(0, log(10))`
  - σ ~ `Uniform(0, 1)`
"""

# ╔═╡ 88edd0b1-738b-4a86-8995-abdfcc2bcc85
md"""
!!! note
	The `LogNormal` distribution is a bit weird: `LogNormal(μ, σ)` gives the distribution of **the exponential** of a normally distributed value with mean μ and standard deviation σ:
	```math
	\begin{gather}
	X \sim \text{Normal}(μ, σ) \, ,
	\\ \Rightarrow \text{exp}(X) \sim \text{LogNormal}(μ, σ) \, .
	\end{gather}
	```
	This means that **μ** is not actually the mean of a `LogNormal(μ, σ)`, but something closer to **log(μ)** (it's complicated). Hence the `log(10)` in the prior above. 
"""

# ╔═╡ 05622800-ccd9-4aa8-9a95-a1fbe34aa76c
@model function expon(num_pearls)
	μ_exp ~ missing
	pearls = zeros(num_pearls)
	for i in 1:num_pearls
		pearls[i] ~ missing
	end
end

# ╔═╡ c5b042bb-c1d9-4b0e-a2bd-a2c367019da3
@model function lognorm(num_pearls)
	μ_lognorm ~ missing
	σ_lognorm ~ missing
	pearls = zeros(num_pearls)
	for i in 1:num_pearls
		pearls[i] ~ missing
	end
end

# ╔═╡ c72cc55c-a5f7-47f5-b5c3-49af9a5db50d
md"Instantiate the models and condition them on the available data."

# ╔═╡ 88306e52-d975-4bd5-a8f0-7e9293c9fb82
expmodel = missing

# ╔═╡ a7cf9228-7515-4d53-b689-bc711283fecc
lognormmodel = missing

# ╔═╡ 811ee6ce-bb0f-4bcf-95be-6de95ab9166b
md"### Maximum likelihood"

# ╔═╡ 392ad62d-b7ef-4951-b729-01e41684cbcf
md"""
Determine the maximum likelihood estimation (MLE) of the parameter values given the data, using the `NelderMead()` algorithm. Plot the fitted parameters on the data for a visual comparison.
"""

# ╔═╡ 3b16c753-a2ba-4b59-a776-b35c4dd2e927
exp_res = missing

# ╔═╡ 1ce90419-083a-475a-8334-32d73dafb606
exp_mean = missing

# ╔═╡ 654d0ec1-92b0-4220-91a8-607978f1cf4f
lognorm_res = missing

# ╔═╡ 17e90886-1a41-4c19-b93a-96d8bd7a075f
lognorm_mean = missing

# ╔═╡ 27883caa-5b08-4f93-8d1a-ae3d85faa9d1
lognorm_spread = missing

# ╔═╡ 62f000f3-6f21-4c1e-ab37-29a122b84088
begin
	histogram(pearlsizes, normalize = :pdf)
	# add plot of best fit exponential distribution
end

# ╔═╡ 7dc215e7-83ef-430e-9f99-24cada961964
begin
	histogram(pearlsizes, normalize = :pdf)
	# add plot of best fit LogNormal distribution
end

# ╔═╡ abffd9fb-dd0a-4145-8b7d-8e233e498735
md"### Bayes factor"

# ╔═╡ 1adee574-150d-4df6-94e0-6aa20ad92389
md"""
Compare both models using the Bayes factor $K$. Start off by calculating the model evidence $P(D \mid M)$ of the data $D$ for each model $M$, approximating the integral with a [Riemann sum](https://en.wikipedia.org/wiki/Riemann_sum):
"""

# ╔═╡ 461c10b2-13f6-442c-ae29-a5f5116a0d01
md"""
```math
P(D \mid M) =\int_{\theta\in\Theta} P(D \mid M, \theta) \, P(\theta) \, d \theta \approx \sum_{i} P(D \mid M, \theta_i) \, P(\theta_i) \, \Delta \theta_i
```
"""

# ╔═╡ feee9ca3-05c5-480b-8cb9-19a3e0c536c1
md"""
The figure below illustrates the different probabilities involved. The red curve is the product of the two curves above, and the area underneath it is the model evidence we want to calculate.
"""

# ╔═╡ 259e3755-eac2-4bba-8c15-a54e326f365d
prior_exp(m) = exp(logprior(expmodel, (μ_exp=m,)));

# ╔═╡ 85a07f1b-df87-4531-ba2a-ccae822379c4
likelihood_exp(m) = exp(loglikelihood(expmodel, (μ_exp=m,)));

# ╔═╡ 8fce7eb3-1371-4d3b-b9cb-8730e18d2be4
posterior_exp(m) = exp(logjoint(expmodel, (μ_exp=m,))) ;
	# prior * likelihood: not yet normalized with evidence!

# ╔═╡ 304b0c2e-7169-403d-baf3-141b0b35f560
let
	xs = 0.1:0.1:15
	ys = [posterior_exp(x) for x in xs]

	p_likelihood = plot(x -> likelihood_exp(x), xlims = (0, 15),
		label = "Likelihood: P(D | M, μ)", color = :blue, width = 2
	)
	p_prior = plot(prior_exp, label = "Prior: P(μ)", color = :cyan, width = 2, xlims = (0, 15))
	p_post = plot(xs, ys, label = "Unnormalized posterior: P(D| M)",
		color = :red, width = 2, line = :dash, xlims = (0, 15), xlabel = "μ_exp",
		ribbon = (ys, zeros(length(xs))), 
		yticks = round.(0:maximum(ys)/10:maximum(ys), sigdigits = 1)
	)
	plot(p_likelihood, p_prior, p_post, ylabel = "density", plottitle = "Evidence", layout = (3, 1))

end

# ╔═╡ 6d39160d-5e5a-4702-bd3f-617d6e7c5513
Δm = 0.1

# ╔═╡ 89547bb1-9027-47c2-a8c0-0212a833879f
begin
	evidence_exp = 0.0
	for m in 0.1:Δm:10
		likelihood_per_point = [missing for pearlsize in pearlsizes]
		likelihood = missing
		prior = missing
		evidence_exp += likelihood * prior * Δm
	end
	println(evidence_exp)
end

# ╔═╡ 3a0d70d3-92d6-4254-b77b-58b75947548d
Δs = 0.01

# ╔═╡ c3fc3e53-3f60-4444-a99c-3a0ef7ce5a94
begin
	evidence_lognorm = 0.0
	for m in 0.1:Δm:log(10)
		for s in 0.1:Δs:1.0
			likelihood_per_point = [
				missing
				for pearlsize in pearlsizes
			]
			likelihood = missing
			prior = missing
			evidence_lognorm += likelihood * prior * Δm * Δs
		end
	end
	println(evidence_lognorm)
end

# ╔═╡ 89683590-7f48-4ee8-bc5c-46d37fc0067b
md"""
Now calculate the Bayes factor as follows:
```math
K = \frac{P(M_2 \mid D)}{P(M_1 \mid D)} = \frac{P(D \mid M_2) \, P(M_2)}{P(D \mid M_1) \, P(M_1)}
```
"""

# ╔═╡ 8c63219d-e520-40b7-a4b3-c6e2ffc7b17e
P_M_exp = 0.5

# ╔═╡ 84f99a01-f3e5-4197-998d-2bdec0d53fb1
P_M_lognorm = 1 - P_M_exp

# ╔═╡ 37f977de-f1a9-4a83-85f1-8dadf7910025
bayes_factor = missing

# ╔═╡ a7fd0d4d-e0de-4f4b-84a8-64c94a300fba
md"""
Another comparison we can make between the models is calculating whether the first model is the correct one:
```math
\begin{align}
P(M_1 \mid D) &= \frac{P(D \mid M_1) \, P(M_1)}{P(D)} \, ,
\\&= \frac{P(D \mid M_1) \, P(M_1)}{P(D \mid M_1) \, P(M_1) + P(D \mid M_2) \, P(M_2)} \, .
\end{align}
```
"""

# ╔═╡ ce13fb27-8b5a-4efd-b415-4d0bfd6d3a6a
P_M_exp_cond_D = missing

# ╔═╡ 5b34c046-bc87-474d-8c20-b5732c0de219
md"### AIC"

# ╔═╡ cc4e9632-53ab-48c9-b179-ca2b10cbfabc
md"""
Using the likelihoods calculated above, calculate the Akaike Information Criterion (AIC) for both models:
"""

# ╔═╡ 964706fe-e48c-4cf6-8e80-484867719bf6
md"""
```math
\text{AIC} = 2 k - 2 \, \text{log}(L)
```
"""

# ╔═╡ 03a7e04d-c16f-468d-ba84-314effa6d838
md"""
!!! tip
	To get your model's best possible **AIC** value, you need the highest possible loglikelihood. By definition, this corresponds with your **MLE**. If `opt_res` is the variable returned by the `optimize` function, you can get the correspondig maximal loglikelihood using `opt_res.lp`.
"""

# ╔═╡ 2e1fb72c-fea9-4de1-bed5-ee947306bc65
AIC(num_params, loglikelihood) = missing

# ╔═╡ 3553ed7d-07a0-418a-8d5f-6791e62d9e10
AIC_exp = missing

# ╔═╡ 257f0675-c5fe-45c5-aba0-ac808efb0ead
AIC_lognorm = missing

# ╔═╡ 1936f7bb-ad6b-4759-b37e-a9b73a8dabc2
md"### BIC"

# ╔═╡ 3867fe78-e317-4d89-9a0b-45068597c015
md"""
Do the same for the (dissapointingly non-Bayesian) Bayesian Information Criterion (BIC):
"""

# ╔═╡ cd826270-0b8e-4437-931a-5cb6f01a1e56
md"""
```math
\text{BIC} = k \, \text{log}(n) - 2 \, \text{log}(L)
```
"""

# ╔═╡ 01648f73-8b75-4158-a7ea-00e55ece7548
BIC(num_observations, num_params, loglikelihood) = missing

# ╔═╡ 1e685af1-6abe-4849-b07b-2ab1da058edc
BIC_exp = missing

# ╔═╡ 23a7d39c-0c0f-43fb-88be-ddfaba988a2c
BIC_lognorm = missing

# ╔═╡ 4cf6a2d0-e58c-4b80-abf9-525097d75f35
md"## Overlapping cells"

# ╔═╡ a1f2199f-29ee-4908-9ad6-bdebcf043de7
md"""
When counting cells, overlapping cells are a common cause of errors. Here we will tackle a simplified version of the problem where we try to distinguish whether a point cloud originates from one or two circles.
"""

# ╔═╡ 5407d952-c5a9-4145-8e4a-4b3bb65e6c19
md"""
![Overlapping cell picture](https://media.springernature.com/full/springer-static/image/art%3A10.1007%2Fs11334-022-00478-y/MediaObjects/11334_2022_478_Fig1_HTML.png?as=webp)
> **Source:** Efficient detection and partitioning of overlapped red blood cells using image processing approach (Dhar 2022)
"""

# ╔═╡ 54dc6f54-fb0e-4d80-b32c-c5ce8f8974c6
md"### Data"

# ╔═╡ bb19f97a-4778-4fc5-8739-1f03bc5416a8
cell_data = [[0.68 -1.34 -0.53 0.5 -1.85 0.68 0.57 -1.55 0.16 -0.04 1.06 1.34 1.41 -1.67 -1.56 -0.51; -0.36 -1.73 -0.4 -1.5 -0.97 1.62 -3.71 -1.98 0.9 1.55 -1.82 -4.56 2.46 2.18 -1.23 -1.06], [-0.3 -0.29 0.99 2.58 -0.38 -2.16 -1.51 -0.36 0.9 1.27 -0.3 0.77 -0.6 -0.94; 0.73 -0.63 -1.67 0.39 2.15 0.29 -0.91 -2.4 -0.18 2.23 2.05 1.49 -0.16 0.49], [-3.23 -1.51 -2.78 1.1 2.52 0.76 -1.34 -3.79 0.39 0.76; 0.08 0.63 -0.11 2.2 1.48 2.94 -0.82 -0.87 0.38 2.21], [-1.56 0.53 1.02 -0.53 -2.08 -1.22 -0.12 1.04 -0.95 0.74; -0.18 0.04 -1.19 -0.76 -0.58 -0.69 0.88 -1.1 -0.93 1.72], [0.11 0.57 -2.06 1.59 1.45 1.11 -2.2 1.24 0.89 0.67 -0.17 1.21 -0.89 1.01 -0.01 1.9 1.26 -1.48 0.6 -0.74 1.6; -0.45 -0.56 0.53 -0.45 -2.05 2.68 -1.75 0.35 -0.67 -0.44 -0.4 -0.79 -2.12 -2.59 -1.31 -1.66 0.54 -0.2 -3.03 -0.16 -0.56], [-0.21 0.36 -0.89 -0.83 -0.36 -1.75 -2.84 0.46 1.1 3.34 -1.61 0.08; -0.38 -2.23 0.27 -1.6 -2.72 -1.87 -1.48 -0.1 -0.83 0.26 0.46 0.57], [3.15 -0.1 0.77 1.62 -0.5 0.28 0.66 -0.01 1.93 -0.15 -0.94 -0.42 1.79 0.27 -0.01 1.7 0.96 2.35 1.61; -0.05 0.28 0.06 -1.26 1.64 -0.48 0.42 1.47 1.05 0.03 -0.65 -0.74 0.26 0.89 1.43 -0.83 -1.55 -0.48 1.72], [2.16 1.24 3.64 -1.18 1.11 2.4 1.19 1.14 1.26 1.11 0.95 2.14 1.88 1.5 2.43 0.64 1.84 0.05 -0.83 1.5 4.44; -1.13 -0.33 -0.98 0.34 -3.2 0.41 0.77 0.1 1.33 0.76 -0.73 -2.07 0.64 -1.96 -0.7 -1.34 0.84 -2.28 -0.95 -0.28 0.24], [0.38 2.4 2.14 -0.65 -0.23 1.37 0.7 0.74 -0.17 2.53 -1.42 -0.03; 1.25 -2.24 0.0 1.12 -2.23 0.93 -0.86 0.89 -1.61 0.93 -1.51 1.58], [1.53 0.05 -0.39 -1.14 0.04 0.36 0.78 -3.02 -0.28 -2.49 -0.3 -0.55 -1.58 -0.24; 2.5 1.84 -0.67 -1.69 1.57 0.57 1.96 -1.94 3.22 1.47 0.57 0.45 -0.23 0.93]];

# ╔═╡ eb552ae5-f32d-4383-8aae-99999ce42552
md"You can choose the cell picture and visualize the data here."

# ╔═╡ d3f98c65-de66-4a3c-b087-c2ef34340110
md"Picture idx"

# ╔═╡ 0533c731-82b2-4dab-8c5e-d5913ee0f4f3
@bind picture_idx Slider(1:length(cell_data), show_value = true)

# ╔═╡ a4c89c22-636c-4621-ac28-b285cf2ecbef
xs, ys = eachrow(cell_data[picture_idx]);

# ╔═╡ 272327fd-a587-4c14-80cc-d581ff2d7f27
scatter(xs, ys, xlims = (-5, 5), ylims = (-5, 5))

# ╔═╡ f11ca2fe-dc50-41b5-bf9c-299f1b18a9e2
md"### Model definition"

# ╔═╡ 65a23129-7fc0-48aa-a9f0-c9f885b7e4c2
md"""
The model for one cell is defined as follows:
- The points originate from one pointcloud with a centre (`xm`, `ym`).
- `xm` and `ym` both follow a standard Normal distribution.
- All x-values follow a Normal distribution around `xm` with $σ = 1$.
- All y-values follow a Normal distribution around `ym` with $σ = 1$.
"""

# ╔═╡ c8da1bc6-4455-4577-98bb-d7dd41ff4f06
@model function singlecell(n) # n is number of points
	xm ~ missing
	ym ~ missing

	missing
end

# ╔═╡ 523d6280-7147-4652-b593-cdd802d80b4e
md"""
The model for two cells is very similar:
- The points originate from one of two pointclouds, one with centre (`xm1`, `ym1`), the other with centre (`xm2`, `ym2`).
- `xm1`, `ym1`, `xm2` and `ym2` all follow standard Normal distributions.
- All x-values follow either a Normal distribution ($σ = 1$) around `xm1` or `xm2`, with equal chance for either.
- The same idea goes for the y-values.
"""

# ╔═╡ 740afb3b-4bd6-4516-9a2a-7bbe5f19ccf0
md"""
!!! hint
	To model the likelihood, consider the humble `MixtureModel`.
"""

# ╔═╡ 166f1c0c-1615-47e7-8538-f19cbdaa6923
@model function doublecell(n) # n is number of points
	missing
end

# ╔═╡ 36093548-5b29-4bd0-959a-befdd4da3de5
md"Instantiate and condition the models."

# ╔═╡ a0502199-3e7f-4b50-b383-2e1bb4ff9d41
singlemodel = missing

# ╔═╡ 153ccedb-b437-4906-a2ef-1745b0dbf53e
doublemodel = missing

# ╔═╡ 29e3dc3a-f2c4-45a6-80a0-01bde4d41d98
md"### Maximum likelihood"

# ╔═╡ 5f92804b-49b2-4fa3-a10a-38c66feb3ce8
function plotsinglecell(xm, ym; bounds = 5)
	mydist = MvNormal([xm, ym], [1.0 0.0; 0.0 1.0])
	
	xs = -bounds:0.1:bounds
	ys = -bounds:0.1:bounds
	
	f(x,y) = pdf(mydist, [x, y])
	contourf(xs, ys, f, xlims = (-bounds, bounds), ylims = (-bounds, bounds),
		color = :viridis, aspect_ratio = :equal, legend = false,
		title = "Single cell model"
	)
end

# ╔═╡ 0977b0e1-092a-4934-9df3-674cdf12b367
function plotdoublecell(xm1, xm2, ym1, ym2; bounds = 5)
	mydist = MixtureModel(
		[
			MvNormal([xm1, ym1], [1.0 0.0; 0.0 1.0]),
			MvNormal([xm2, ym2], [1.0 0.0; 0.0 1.0]),
		]
	)
	
	xs = -bounds:0.1:bounds
	ys = -bounds:0.1:bounds
	
	f(x,y) = pdf(mydist, [x, y])
	contourf(xs, ys, f, xlims = (-bounds, bounds), ylims = (-bounds, bounds),
		color = :viridis, aspect_ratio = :equal, legend = false,
		title = "Two cells model"
	)
end

# ╔═╡ 15b5d7d1-cf91-4024-a7fd-b8b8f0563dff
md"""
Determine the maximum likelihood estimation (MLE) of the parameter values given the data, using the `NelderMead()` algorithm.
"""


# ╔═╡ 8942d992-5084-4c2c-9cbc-c442b0381922
singleres = missing

# ╔═╡ 99a8a4b4-bbcd-4321-ba41-4f307c04af8e
begin
	single_xm = missing
	single_ym = missing
end

# ╔═╡ 51ce7374-c7d8-4a32-815d-1c2ff9ba9970
doubleres = missing

# ╔═╡ 5a0860a2-8a61-4abb-acaf-3ea3cad0c286
begin
	double_xm1 = missing
	double_xm2 = missing
	double_ym1 = missing
	double_ym2 = missing
end

# ╔═╡ 73d259c6-f776-4586-a24f-3ecc368e28ae
md"Visualise the results"

# ╔═╡ 1a00a07e-db38-4761-b881-7f475881ff4f
begin
	plotsinglecell(single_xm, single_ym)
	scatter!(xs, ys)
end

# ╔═╡ 629280c7-b6c8-423d-92e9-815eb76a78f0
begin
	plotdoublecell(double_xm1, double_xm2, double_ym1, double_ym2)
	scatter!(xs, ys)
end

# ╔═╡ 4e28f87c-0953-48c4-8e3a-d7e853cae816
md"### AIC"

# ╔═╡ ec4ff409-f28c-4640-8bdb-abe21252afaf
md"Using the MLE results from the previous section, determine the AIC of both models. You can use the implementation from previous exercise."

# ╔═╡ fef14ffe-3399-48d8-a0cb-998826cdfda4
AIC_single = missing

# ╔═╡ 40077482-7c98-4ca6-b441-ed3ebd893369
AIC_double = missing

# ╔═╡ Cell order:
# ╟─fb07436c-cc5d-4721-8a80-d7f7201721d7
# ╠═a18638d0-138d-11f0-0571-33060f7da7ba
# ╠═fe265128-33a2-40d1-ab71-bc5c53979a54
# ╠═c870d28d-7b2e-48f8-8990-b2cc0943cb09
# ╠═2a73d23c-73fc-4845-b359-c6fe22077140
# ╠═a397f48e-4228-435e-af13-c2bc71c8cb05
# ╟─a6635715-52e3-44e7-9c00-0be751f830d6
# ╟─473c003a-a34b-4e92-8a82-99450514d755
# ╟─a2dec4e7-7c16-4ae5-974e-8c01ec31507a
# ╟─c4e09491-fd0a-43e7-8e74-04832a045a48
# ╟─aba7efa9-a0c3-4da4-82a2-bb7efc6b9e22
# ╟─c8b56408-0710-4646-bb2b-e336e7f689c2
# ╟─19dd461d-cc03-41c9-a3d1-fcd477d08eb0
# ╟─01a530db-5e18-4540-99a9-eecd9e61e1af
# ╟─6d39fbde-81c5-4786-99e1-b587bd15f94c
# ╟─fce554e2-a6da-4046-ab4e-9153e27aef7a
# ╠═7b1ce376-d796-4913-a488-ee921bd13855
# ╟─712756bc-0408-45a6-931b-eec07b6e052f
# ╟─b4de93ea-a919-41c2-bb9c-bc1beccc9ea8
# ╟─5a6f8d49-3841-46ff-a826-88ba0bed868b
# ╟─88edd0b1-738b-4a86-8995-abdfcc2bcc85
# ╠═05622800-ccd9-4aa8-9a95-a1fbe34aa76c
# ╠═c5b042bb-c1d9-4b0e-a2bd-a2c367019da3
# ╟─c72cc55c-a5f7-47f5-b5c3-49af9a5db50d
# ╠═88306e52-d975-4bd5-a8f0-7e9293c9fb82
# ╠═a7cf9228-7515-4d53-b689-bc711283fecc
# ╟─811ee6ce-bb0f-4bcf-95be-6de95ab9166b
# ╟─392ad62d-b7ef-4951-b729-01e41684cbcf
# ╠═3b16c753-a2ba-4b59-a776-b35c4dd2e927
# ╠═1ce90419-083a-475a-8334-32d73dafb606
# ╠═654d0ec1-92b0-4220-91a8-607978f1cf4f
# ╠═17e90886-1a41-4c19-b93a-96d8bd7a075f
# ╠═27883caa-5b08-4f93-8d1a-ae3d85faa9d1
# ╠═62f000f3-6f21-4c1e-ab37-29a122b84088
# ╠═7dc215e7-83ef-430e-9f99-24cada961964
# ╟─abffd9fb-dd0a-4145-8b7d-8e233e498735
# ╟─1adee574-150d-4df6-94e0-6aa20ad92389
# ╟─461c10b2-13f6-442c-ae29-a5f5116a0d01
# ╟─feee9ca3-05c5-480b-8cb9-19a3e0c536c1
# ╟─304b0c2e-7169-403d-baf3-141b0b35f560
# ╠═259e3755-eac2-4bba-8c15-a54e326f365d
# ╠═85a07f1b-df87-4531-ba2a-ccae822379c4
# ╠═8fce7eb3-1371-4d3b-b9cb-8730e18d2be4
# ╠═6d39160d-5e5a-4702-bd3f-617d6e7c5513
# ╠═89547bb1-9027-47c2-a8c0-0212a833879f
# ╠═3a0d70d3-92d6-4254-b77b-58b75947548d
# ╠═c3fc3e53-3f60-4444-a99c-3a0ef7ce5a94
# ╟─89683590-7f48-4ee8-bc5c-46d37fc0067b
# ╠═8c63219d-e520-40b7-a4b3-c6e2ffc7b17e
# ╠═84f99a01-f3e5-4197-998d-2bdec0d53fb1
# ╠═37f977de-f1a9-4a83-85f1-8dadf7910025
# ╟─a7fd0d4d-e0de-4f4b-84a8-64c94a300fba
# ╠═ce13fb27-8b5a-4efd-b415-4d0bfd6d3a6a
# ╟─5b34c046-bc87-474d-8c20-b5732c0de219
# ╟─cc4e9632-53ab-48c9-b179-ca2b10cbfabc
# ╟─964706fe-e48c-4cf6-8e80-484867719bf6
# ╟─03a7e04d-c16f-468d-ba84-314effa6d838
# ╠═2e1fb72c-fea9-4de1-bed5-ee947306bc65
# ╠═3553ed7d-07a0-418a-8d5f-6791e62d9e10
# ╠═257f0675-c5fe-45c5-aba0-ac808efb0ead
# ╟─1936f7bb-ad6b-4759-b37e-a9b73a8dabc2
# ╟─3867fe78-e317-4d89-9a0b-45068597c015
# ╟─cd826270-0b8e-4437-931a-5cb6f01a1e56
# ╠═01648f73-8b75-4158-a7ea-00e55ece7548
# ╠═1e685af1-6abe-4849-b07b-2ab1da058edc
# ╠═23a7d39c-0c0f-43fb-88be-ddfaba988a2c
# ╟─4cf6a2d0-e58c-4b80-abf9-525097d75f35
# ╟─a1f2199f-29ee-4908-9ad6-bdebcf043de7
# ╟─5407d952-c5a9-4145-8e4a-4b3bb65e6c19
# ╟─54dc6f54-fb0e-4d80-b32c-c5ce8f8974c6
# ╟─bb19f97a-4778-4fc5-8739-1f03bc5416a8
# ╟─eb552ae5-f32d-4383-8aae-99999ce42552
# ╟─d3f98c65-de66-4a3c-b087-c2ef34340110
# ╟─0533c731-82b2-4dab-8c5e-d5913ee0f4f3
# ╠═a4c89c22-636c-4621-ac28-b285cf2ecbef
# ╟─272327fd-a587-4c14-80cc-d581ff2d7f27
# ╟─f11ca2fe-dc50-41b5-bf9c-299f1b18a9e2
# ╟─65a23129-7fc0-48aa-a9f0-c9f885b7e4c2
# ╠═c8da1bc6-4455-4577-98bb-d7dd41ff4f06
# ╟─523d6280-7147-4652-b593-cdd802d80b4e
# ╟─740afb3b-4bd6-4516-9a2a-7bbe5f19ccf0
# ╠═166f1c0c-1615-47e7-8538-f19cbdaa6923
# ╟─36093548-5b29-4bd0-959a-befdd4da3de5
# ╠═a0502199-3e7f-4b50-b383-2e1bb4ff9d41
# ╠═153ccedb-b437-4906-a2ef-1745b0dbf53e
# ╟─29e3dc3a-f2c4-45a6-80a0-01bde4d41d98
# ╟─5f92804b-49b2-4fa3-a10a-38c66feb3ce8
# ╟─0977b0e1-092a-4934-9df3-674cdf12b367
# ╟─15b5d7d1-cf91-4024-a7fd-b8b8f0563dff
# ╠═8942d992-5084-4c2c-9cbc-c442b0381922
# ╠═99a8a4b4-bbcd-4321-ba41-4f307c04af8e
# ╠═51ce7374-c7d8-4a32-815d-1c2ff9ba9970
# ╠═5a0860a2-8a61-4abb-acaf-3ea3cad0c286
# ╟─73d259c6-f776-4586-a24f-3ecc368e28ae
# ╠═1a00a07e-db38-4761-b881-7f475881ff4f
# ╠═629280c7-b6c8-423d-92e9-815eb76a78f0
# ╟─4e28f87c-0953-48c4-8e3a-d7e853cae816
# ╟─ec4ff409-f28c-4640-8bdb-abe21252afaf
# ╠═fef14ffe-3399-48d8-a0cb-998826cdfda4
# ╠═40077482-7c98-4ca6-b441-ed3ebd893369
