### A Pluto.jl notebook ###
# v0.20.4

#> [frontmatter]
#> order = "5"
#> title = "example estrogen"
#> date = "2025-02-07"
#> tags = ["project"]
#> description = "Project example estrogen"
#> layout = "layout.jlhtml"
#> 
#>     [[frontmatter.author]]
#>     name = "Michiel Stock"

using Markdown
using InteractiveUtils

# ╔═╡ 21357c48-f35d-11ee-23f8-2534bb1d82f4
begin
    
    using Pkg
    Pkg.activate("../../pluto-deployment-environment")
    
	# make this cell invisible when you are finished
	title = "Estrogen Estimation"
	names = ["Vo Orbeeld", "Pro Ject"]

	x = 4
	academic_year = "202$(x)_202$(x+1)"

	email_main_person = "mail@domain.be"

	using PlutoUI  # interactivity
	using Random # set seed
    using Turing # sampling
	using StatsPlots # plots
	TableOfContents()
end;

# ╔═╡ 14c7e803-c0ff-4211-8b60-2c9c246934dd
md"""
# $title

**$(join(names, ", ", " and "))**
"""

# ╔═╡ c33de800-7aa3-438a-85d0-22b8fba1c8f3
md"""
Note: This project is an adaption of [the following blog post](https://www.oxinabox.net/2022/11/11/Estimating-Estrogen.html) by Dr. Frames Catherine White. Students are (sadly) not permitted to copy existing blog posts for their own projects. 
"""

# ╔═╡ e7629b47-8514-4f82-b3d0-1c51b79dbadf
# ╠═╡ disabled = true
#=╠═╡
begin
	using Pkg
	Pkg.activate("..")
end
  ╠═╡ =#

# ╔═╡ 6948149f-854e-4b3c-b51c-099dd221ab83
md"""
## Abstract

"""

# ╔═╡ 4328d3b7-c175-4bb6-891e-b22af044baac
md"""
As a trans-femme on HRT, I would like to know the concentrations of estradiol in my blood at all hours of day. This is useful as the peak, the trough and average all have effects. However, I only get blood tests a finite number of times per day – usually once. I am not a medical doctor, but I am the kind of doctor who can apply scientific modelling to the task of estimating curves based on limited observations. I am honestly surprised no one has done this. The intersection of trans folk and scientific computing is non-trivial. After all, the hardest problem in computer science is gender dysphoria.
"""

# ╔═╡ 89551690-500d-4e37-ae20-5beb71cc87ac
md"""
## Model
"""

# ╔═╡ 8a6f6d76-0db0-4253-94c5-14b8184b75a8
md"""
In [this blog post](https://web.archive.org/web/20230128040153/http://transascity.org/sublingual-versus-oral-estrogen/) on Sublingual versus Oral Estrogen they approximated the estradiol function with a linear to the peak then an exponential decay.

$$c(t) = \begin{cases}\frac{c_\max t}{t_\max} & \text{if } t \le t_\max\\ c_\max 2^{-(t-t_\max)/t_{1/2}} & \text{if }t_\max < t\,. \end{cases}$$
"""

# ╔═╡ 5b58dc75-6d62-4ba8-b3c3-ddde107d8790
estrogen_conc(t, c_max = 100, t_max = 3, halflife = 3) = ifelse(t < t_max, c_max/t_max * t, c_max * 2^(-(t-t_max)/halflife))

# ╔═╡ bf645a16-e64e-4051-9b31-bc92148adda9
# ╠═╡ disabled = true
#=╠═╡
function estrogen_conc(t, c_max = 100, t_max = 3, halflife = 3)
    if t < t_max
        c = c_max/t_max * t
    else
        c = c_max * 2^(-(t-t_max)/halflife)
    end
    return c
end
  ╠═╡ =#

# ╔═╡ 39ba24b9-0ba1-49a9-b31f-d7cb98d45303
plot(estrogen_conc, xlims = (0, 24), title = "Estrogen concentration model", xlabel = "t (h)", ylabel = "Estrogen concentration c (pg/ml)", legend = nothing)

# ╔═╡ 0212cb1a-664c-43a5-9d03-932168e98a0a
md"""
The curve defines the current blood concentration of estradiol c at time t hours after application of the gel. It is described by 3 parameters:

- `c_max`: the peak concentration.
- `t_max`: the time it takes to reach peak concentration.
- `halflife`: the time it takes for the concentration to half after reaching peak.
"""

# ╔═╡ 2d306ecd-ccee-4a18-9023-0c2f93a959cf
md"""
It’s broadly biologically plausible. We expect a fast initial absorption, that should end at some point in few few hours. Since it is fast and short, it doesn’t really matter what we model it with, so linear is fine. Then we expect a tail off as it is consumed. It makes sense for the rate of absorption to be related to the quantity remaining – which suggests some exponential. We see this kind of thing very frequently in biological systems. This all might be nonsense, I am no systems biologist.
"""

# ╔═╡ 638d0a38-ffe8-499b-84a7-5b49f8c26e47
md"""
## Calibration
"""

# ╔═╡ e4cbfc65-db9a-47cd-9f99-b908045603eb
md"""
Järvinen et al. (1997) give 3 curves for single dose. I am going to plot the data from Järvinen et al against curves using my formula, best fit by my own inspection. I am downshifting all the data from Järvinen et al by 25 pg/mL, as that data was from post-menopausal cis women, who produce about 25 pg/mL of estradiol on their own before you take into account HRT. We only want to model the HRT component.
"""

# ╔═╡ 120058d5-18d9-4f0f-871b-8d5f57133f4a
md"### Visual inspection"

# ╔═╡ 7b17c6ce-614e-4831-9b20-3d6b250797bd
t_obs = [0, 1, 2, 3, 4, 6, 8, 10, 12, 16, 24];

# ╔═╡ 0fb4b2e2-463a-4621-af36-d6de04db6542
c_obs = [
	[0, 25, 100, 132, 90, 82, 60, 55, 32, 15, 4],
	[25, 35, 70, 75, 55, 45, 35, 32, 22, 15, 4],
	[5, 14, 17, 20, 12, 10, 5, 2, 5, 5, 4]
];

# ╔═╡ f3835020-f588-4031-bf0b-e2aba72a4e82
color_palette = [RGB(91/255, 206/255, 250/255) RGB(245/255, 169/255, 184/255) RGB(1, 1, 1)]; # custom colors for the occassion

# ╔═╡ 132d574c-76e8-48b8-b9be-f6733de4b037
p_data = scatter(t_obs, c_obs, color = color_palette, bg = :lightgray, label = ["A200 obs" "A400 obs" "Amax obs"], xlabel = "t (h)", ylabel = "Estrogen concentration (pg/ml)")

# ╔═╡ 446a89ee-a31e-404c-aca2-66d281c85912
estimated_funcs = [
	t -> estrogen_conc(t, 132, 3, 3.5), t -> estrogen_conc(t, 100, 2.5, 3.5), t -> estrogen_conc(t, 20, 3.5, 2.7)
]

# ╔═╡ a5c0f97e-d62a-4efd-a61c-04d92527f1b5
plot!(p_data, estimated_funcs, color = color_palette, label = ["A200 pred" "A400 pred" "Amax pred"])

# ╔═╡ a015d9cd-c089-404c-a007-a76ef13d31ef
md"""
By looking at these plots, it seems a pretty decent model. Of-course with enough degrees of freedom, you can fit an elephant. However, we have 10 points and only 3 degrees of freedom, of which we only varied 2 of them across the 3 datasets. So it seems like we are good.
"""

# ╔═╡ 4152a577-c91e-483b-abf6-54370ecec988
md"""
Now I just fit those curves by eye. We can find the the most likely parameters via least squares regression. But really we are not after a single curve at all. We are interested in distributions over possible curves, given the observations. These tell use the possible realities that would explain what we are seeing.
"""

# ╔═╡ f27f10a6-ec2b-40f2-8015-a060b41a1053
md"### Bayesian inference"

# ╔═╡ 65d69035-d508-4430-84d4-8573f58033b3
md"""
To begin with lets think about our priors. These are our beliefs about the values the parameters might take before we look at the data.

- `c_max` is somewhere between 0 and 500 pg/mL (ie. 0-1835 pmol/L). If your E2 is above that something is very wrong. For now let’s not assume anything more and just go with a Uniform distribution. Though perhaps we could do something smarter hand select something that tailed off nicely towards the ends.
- `t_max` is somewhere between 1 and 4 hours, we know this because the instruction say don’t let anyone touch you for the first hour (so its definitely still absorbing then), and common wisdom is to not wash the area for at least 4 hours – so it must be done but then. If we use a Triangular distribution it has some push towards the center.
- `halflife`, we know this has to be positive, since otherwise it would not decay. Being log-normal makes sense since it appears in an exponential. We would like it to have mode of 3.5 since that is what by eye we saw fit the curves all nicely (probably bad Bayesian cheating here) and because that means it is mostly all decayed by 24 hours – it can’t all that much higher usually since otherwise wouldn’t need daily doses, nor that much lower since in that case would need multiple doses per day. To set the mode to 3.5 we use `LogNormal(log(3.5)+1, 1)`
"""

# ╔═╡ 94775f83-0df3-4b1d-b1e2-8ff824177cd1
plot(
    [Uniform(0, 500), TriangularDist(1, 4), LogNormal(log(3.5)+1, 1)],
    legend = false, linewidth = 2, title = ["c_max" "t_max" "halflife"], 
	layout = (1, 3)
)

# ╔═╡ fc3ef0fe-a8a9-4acf-94d1-ac8548fbbf56
md"""
The other component we will want is an error term. We want to express our observations of the concentration as being noisy samples from a normal distribution centered on the actual curve we are estimating. So we need an error term which will allow some wiggle room about that curve, without throwing off the inference for the real parameters. 

We will define a variable called `err` which is the standard deviation of this error term. Our prior on this error term should be positive with a peak at 0 and rapidly tailing off. Gamma(1, 1) meets our requirement.
"""

# ╔═╡ b88b22e8-531e-41af-b74b-82b498cbfab0
plot(Gamma(1,1), title="err", legend=false)

# ╔═╡ 402d619a-0ab6-4ddc-9189-ebe35841b5cc
@model function single_dose(t_obs, c_obs)
    c_max ~ Uniform(0, 500)
    t_max ~ TriangularDist(1, 4)
    halflife ~ LogNormal(log(3.5) + 1, 1)
    
    err ~ Gamma(1, 1)
    
    for i in eachindex(t_obs)
        c_pred = estrogen_conc(t_obs[i], c_max, t_max, halflife)
        c_obs[i] ~ Normal(c_pred, err)
    end
end

# ╔═╡ 709392a8-78fa-48ea-8dc0-d92118651af6
chain = sample(single_dose(t_obs, c_obs[1]), NUTS(), 2000);

# ╔═╡ 07e1d932-f49f-4fcf-9d95-eaea3254b99a
plot(chain)

# ╔═╡ 7f19f097-c048-496d-ba9b-2b672b12b7cf
md"So let’s look at the distribution over curves (as represented by samples)."

# ╔═╡ cf42c732-67cb-4899-8ac0-897aaa39da2f
md"""
We see this nice kinda clear and fairly small range of values for the parameters: `c_max`, `t_max`, `halflife`. The error term, `err`, is quite large
"""

# ╔═╡ 70d2441c-d2c6-4f76-b553-c4a1ebfeb6dd
md"### Using less datapoints"

# ╔═╡ ba55f46a-6cb8-4e12-9ecc-46c965c705b5
md"""
Now that we have shown we can do inference to find distributions over parameters that fit the data let’s get on to a more realistic task. No one gets blood tests every few hours outside of an experimental data gathering exercise. The most frequent blood tests I have heard of is every 2 weeks, and most are more like every 3-6 months. So what we are really interested in is inferring what could be happening with blood levels from a single observation.
"""

# ╔═╡ 919e9ea3-f584-41ec-ba29-a9fe210aa7de
chain_1 = sample(single_dose([8], [60]), NUTS(), 2000);

# ╔═╡ f7f9d65a-2f09-486d-8afd-4e4a15294f1b
md"""
So that’s actually really informative. There are a range of possible explanations. From a very small t_max and a large c_max meaning it peaked early and has tailed off a lot, to the more likely ones which look more like the kind of curves we were seeing based on the experimental data with more frequent measurements.
"""

# ╔═╡ 928c6102-5bbd-4072-971e-f09537884b30
md"""
We can add more observation points and cut-down the number of realities we might be in. This realistically is actually a practical thing to do. You can see in the following plots that if we add a reading of 60 at 8 hours after application we break the possible universes into two possible sets of explanations. One set where the 3 hour reading is while it is still rising, and one set where it is falling.
"""

# ╔═╡ e520b02c-eebc-4c03-866b-b00a0f5ee1a3
chain_2 = sample(single_dose([3, 8], [100, 60]), NUTS(), 2000);

# ╔═╡ 756df9f3-3b20-4ff6-9990-a8c431f8849a
chain_3 = sample(single_dose([1, 3, 8], [50, 100, 60]), NUTS(), 2000);

# ╔═╡ 60c63f1b-5a27-4539-a486-78c083457b0e
md"""
## Conclusion
"""

# ╔═╡ 40c387d4-4570-400f-b8e6-a6f35664f3c2
md"""
This is just a first look at this topic. I imagine I might return to it again in the future. Here are some extra things we might like to look at:

- Determining optimal times to test: 3 readings will not always capture the curve, different times may be more informative than others, especially when we consider the error level vs the signal level.
- Average levels: the distribution of average level is likely fairly collapsed – multiple different sets of parameter values can lead to same average level.
- Multi-day: Since estrogen doesn’t hit zero at 24 hours can model across days. Can also include a term for variation in when it was applied in the day since people are not that consistent. Multiday is crucial for making the model realistic
- Dose changes: extending beyond multiday, people change there does, and we know higher dose leads to higher levels so we can insert that prior knowledge.


Probabilistic programming is a cool technique for working on pharmacodynamics. It lets us handle the fact that we have many unknowns about people’s individual biology, while still narrowing down a possible set of worlds they might live in.
"""

# ╔═╡ 890afefc-f42b-4d74-b775-6dee5e5f0c2b
md"## Appendix"

# ╔═╡ 8c0ed057-e5de-4f2e-ad65-f9235b9c1b7d
function plot_estrogen_estimations(chain, t_obs, c_obs)
    c_max, t_max, halflife = [
		chain[param] for param in [:c_max, :t_max, :halflife]
	]
    est_funcs = [
		t -> estrogen_conc(t, c_max[i], t_max[i], halflife[i])
		for i in eachindex(c_max)
	]
    
    plot(est_funcs, color = color_palette[2], xlims = (0, 25), ylims = (0, 200),
		label = nothing, opacity = 0.01);
    scatter!(t_obs, c_obs, color = color_palette[1], label = nothing)
end

# ╔═╡ 1b9fdff3-dc7a-4969-a027-9c61dc8fc2a4
plot_estrogen_estimations(chain, t_obs, c_obs[1])

# ╔═╡ 1c2bd15b-44f9-4cb7-9b81-435fcf43bb5e
plot_estrogen_estimations(chain_1, [8], [60])

# ╔═╡ 30355b66-de47-4712-b2a8-27302c5e2931
plot_estrogen_estimations(chain_2, [3, 8], [100, 60])

# ╔═╡ ad3760d6-1e0f-4a63-bb23-e04b4a4023d8
plot_estrogen_estimations(chain_3, [1, 3, 8], [50, 100, 60])

# ╔═╡ 946cf45a-329b-42ca-9f77-a3c7bdc5611a
md"""
### References
"""

# ╔═╡ e59da868-a7f5-48f0-af1b-8a990f146f20
md"""
Järvinen, A., Granander, M., Nykänen, S., Laine, T., Geurts, P., & Viitanen, A. (1997). Steady‐state pharmacokinetics of oestradiol gel in post‐menopausal women: effects of application area and washing. BJOG: An International Journal of Obstetrics & Gynaecology, 104, 14-18.
"""

# ╔═╡ Cell order:
# ╟─14c7e803-c0ff-4211-8b60-2c9c246934dd
# ╟─c33de800-7aa3-438a-85d0-22b8fba1c8f3
# ╟─e7629b47-8514-4f82-b3d0-1c51b79dbadf
# ╟─6948149f-854e-4b3c-b51c-099dd221ab83
# ╟─4328d3b7-c175-4bb6-891e-b22af044baac
# ╟─89551690-500d-4e37-ae20-5beb71cc87ac
# ╟─8a6f6d76-0db0-4253-94c5-14b8184b75a8
# ╠═5b58dc75-6d62-4ba8-b3c3-ddde107d8790
# ╟─bf645a16-e64e-4051-9b31-bc92148adda9
# ╟─39ba24b9-0ba1-49a9-b31f-d7cb98d45303
# ╟─0212cb1a-664c-43a5-9d03-932168e98a0a
# ╟─2d306ecd-ccee-4a18-9023-0c2f93a959cf
# ╟─638d0a38-ffe8-499b-84a7-5b49f8c26e47
# ╟─e4cbfc65-db9a-47cd-9f99-b908045603eb
# ╟─120058d5-18d9-4f0f-871b-8d5f57133f4a
# ╠═7b17c6ce-614e-4831-9b20-3d6b250797bd
# ╠═0fb4b2e2-463a-4621-af36-d6de04db6542
# ╟─f3835020-f588-4031-bf0b-e2aba72a4e82
# ╟─132d574c-76e8-48b8-b9be-f6733de4b037
# ╟─446a89ee-a31e-404c-aca2-66d281c85912
# ╟─a5c0f97e-d62a-4efd-a61c-04d92527f1b5
# ╟─a015d9cd-c089-404c-a007-a76ef13d31ef
# ╟─4152a577-c91e-483b-abf6-54370ecec988
# ╟─f27f10a6-ec2b-40f2-8015-a060b41a1053
# ╟─65d69035-d508-4430-84d4-8573f58033b3
# ╟─94775f83-0df3-4b1d-b1e2-8ff824177cd1
# ╟─fc3ef0fe-a8a9-4acf-94d1-ac8548fbbf56
# ╟─b88b22e8-531e-41af-b74b-82b498cbfab0
# ╠═402d619a-0ab6-4ddc-9189-ebe35841b5cc
# ╠═709392a8-78fa-48ea-8dc0-d92118651af6
# ╠═07e1d932-f49f-4fcf-9d95-eaea3254b99a
# ╟─7f19f097-c048-496d-ba9b-2b672b12b7cf
# ╟─1b9fdff3-dc7a-4969-a027-9c61dc8fc2a4
# ╟─cf42c732-67cb-4899-8ac0-897aaa39da2f
# ╟─70d2441c-d2c6-4f76-b553-c4a1ebfeb6dd
# ╟─ba55f46a-6cb8-4e12-9ecc-46c965c705b5
# ╠═919e9ea3-f584-41ec-ba29-a9fe210aa7de
# ╠═1c2bd15b-44f9-4cb7-9b81-435fcf43bb5e
# ╟─f7f9d65a-2f09-486d-8afd-4e4a15294f1b
# ╟─928c6102-5bbd-4072-971e-f09537884b30
# ╠═e520b02c-eebc-4c03-866b-b00a0f5ee1a3
# ╠═30355b66-de47-4712-b2a8-27302c5e2931
# ╠═756df9f3-3b20-4ff6-9990-a8c431f8849a
# ╠═ad3760d6-1e0f-4a63-bb23-e04b4a4023d8
# ╟─60c63f1b-5a27-4539-a486-78c083457b0e
# ╟─40c387d4-4570-400f-b8e6-a6f35664f3c2
# ╟─890afefc-f42b-4d74-b775-6dee5e5f0c2b
# ╟─8c0ed057-e5de-4f2e-ad65-f9235b9c1b7d
# ╟─21357c48-f35d-11ee-23f8-2534bb1d82f4
# ╟─946cf45a-329b-42ca-9f77-a3c7bdc5611a
# ╟─e59da868-a7f5-48f0-af1b-8a990f146f20
