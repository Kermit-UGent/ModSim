### A Pluto.jl notebook ###
# v0.20.4

#> [frontmatter]
#> order = "3"
#> title = "example antibiotics"
#> date = "2025-02-07"
#> tags = ["project"]
#> description = "Project example antibiotics"
#> layout = "layout.jlhtml"
#> 
#>     [[frontmatter.author]]
#>     name = "Michiel Stock"

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

# ╔═╡ 21357c48-f35d-11ee-23f8-2534bb1d82f4
begin

    using Pkg
    Pkg.activate("../../pluto-deployment-environment")
	
	# make this cell invisible when you are finished
	title = "Modelling the effect of antibiotics on microbial resistence"
	names = ["Michiel"]

	academic_year = "2023_2024"

	email_main_person = "mail@domain.be"

	using PlutoUI  # interactivity
	using Plots  # plotting
	TableOfContents()
end

# ╔═╡ 1e078c8b-2df9-4646-b629-bac7f151e937
using Catalyst, DifferentialEquations

# ╔═╡ 158c4ad8-a38b-4828-812f-1abce60695dc
using SciMLSensitivity

# ╔═╡ 6785c1af-28e5-46b4-ba34-a8ab65be7ff8
using ForwardDiff

# ╔═╡ 14c7e803-c0ff-4211-8b60-2c9c246934dd
md"""
# $title

**$(join(names, ", ", " and "))**
"""

# ╔═╡ 6948149f-854e-4b3c-b51c-099dd221ab83
md"""
## Abstract

The discovery of antibiotics is one of the greatest medical advancements of the 20th century. In this project, we use a simple ordinary differential equation (ODE) system to model the effect of antibiotic dosing on a system containing susceptible and resistant bacteria. Bacteria grow with a simple completion model (susceptible bacteria grow faster due to the fitness cost associated with being resistant). Antibiotics can be added to the system. Higher concentrations kill bacteria more effectively, though antibiotics are quickly removed from the system. When choosing doses and dosing time, can we maximally reduce bacterial infection while limiting our total use?
"""

# ╔═╡ 89551690-500d-4e37-ae20-5beb71cc87ac
md"""
## Model

In this system, we model three variables: the number of susceptible ($S(t)$) and resistant ($R(t)$) bacteria and the concentration of antibiotics ($C(t)$). 

The following processes take place:
- the growth rate of the bacteria depends on their total density according to a logistic growth: $\mu = r(1-\frac{S+R}{K})$.
- susceptible bacteria have a fitness cost of $a$, limiting their growth rate $\mu(1-a)$;
- both $S$ and $R$ are removed according to a first-order process with rate $\theta$;
- the antibiotic kills off bacteria according to a first-order kinetics, with the rate determined by a Hill function. Resistant bacteria have twice as high $K_s$.
- antibiotics leave the system (degradation and removal) with a rate of $g$.
"""

# ╔═╡ aed58771-86c1-4928-a1ba-b7d2f503b188
antibiotics = @reaction_network begin
	@species S(t)=500.0 R(t)=50.0 C(t)=0.0
	r * (1-(S+R)/K), S --> 2S  # growth susceptible bacteria
	r * (1-(S+R)/K) * (1-a), R --> 2R  # growth resistant bacteria (with fitness cost)
	θ, (S, R) --> ∅  # natural removal of the bacteria
	hill(C, v, Ks, 4), S --> ∅  # killing S bacteria via antibiotics 
	hill(C, v, 2Ks, 4), R --> ∅  # killing R bacteria via antibiotics
	g, C --> ∅  # removal of antibiotics
end

# ╔═╡ 33a6968f-d8df-4f56-baf3-6e0e88e54c23
md"These reactions form the following system of ordinary differential equations:"

# ╔═╡ 83b50ef8-9b68-428a-bf93-1583d032e9aa
convert(ODESystem, antibiotics)

# ╔═╡ ccb31a3a-f75e-479d-bec7-75cf590652c5
parameters(antibiotics)

# ╔═╡ f284fda5-aaf7-4121-a8f0-b996d501bec5
md"## Simulation and analysis"

# ╔═╡ aaf3d52e-7e6b-4742-b9df-5440b148e7ee
md"### One dose of antibiotics"

# ╔═╡ a2164307-6d18-4983-9aaa-25c8f5828655
md"We set some sensible parameter values:"

# ╔═╡ eab8adc0-e3c7-4693-990b-cb47e3d9c2cf
pars = [:r=>2.7, :K=>1e3, :θ=>0.2, :a=>0.2, :g=>log(2), :v=>5.3, :Ks=>4.] 

# ╔═╡ 556243a1-c387-4606-94d9-2dc06dfdc5f0
md"We can simulate the system. Let us assume an initial concentration of antibiotics $C_0$ at $t=0$ and see what the effect is."

# ╔═╡ 98451eb8-8755-47fe-b8c6-ce2959256dd8
tspan = (0.0, 50.0)

# ╔═╡ b62c37ac-8d33-49ac-8aa6-18d959bc6933
@bind C0 Slider(0:100, show_value=true, default=20)

# ╔═╡ 797ce7d7-c49d-4a00-a8e3-058b6de4d3a9
oprob = ODEProblem(antibiotics, [:C=>C0], tspan, pars);

# ╔═╡ aa32c5e2-0f67-4910-a612-d958fa13169c
plot(solve(oprob, RK4()))

# ╔═╡ 4679507b-e6e4-4d39-8610-77141f7b82e0
md"""
We note that:
- When no antibiotic is present, the susceptible bacteria quickly take over.
- At low concentrations (2-20), the resistant bacteria are killed off, while the resistant bacteria take (temporarily) over.
- High antibiotic concentrations (> 40) show nearly complete eradication of both types of bacteria.

We note that after a while (when all the antibiotics have left the system), the bacteria quickly return to total capacity.

What if we give a second dose of antibiotic at a different time?
"""

# ╔═╡ 49e97512-4150-49e2-bd9b-c3bfd7c375dc
md"### Second dose of antibiotics"

# ╔═╡ 50794b21-e90f-438d-ad9c-eaddfab14ce7
function dose!(integrator, D)
	integrator[:C] += D
end

# ╔═╡ 7f9c8c7d-90ef-496c-b907-300feba5edfc
obprob2 = ODEProblem(antibiotics, [], tspan, pars);  # no AB in the system

# ╔═╡ bf81ccd3-4cdf-4a9c-b334-f7339fed8dfe
@bind D Slider(1.0:100.0, show_value=true)

# ╔═╡ 438de4bb-6d15-4d76-b9a5-95b51924cbc2
@bind tdose Slider(1.0:35.0, show_value=true)

# ╔═╡ 37ecd681-bd01-4162-94cf-97f6d8d82a8e
ps_cb = PresetTimeCallback(tdose, I->dose!(I, D));

# ╔═╡ 0d464337-8a70-47eb-8c08-cc2a7d388cff
begin
	plot(solve(oprob, callback=ps_cb))
	vline!([tdose], label="dosing time", ls=:dash, title="Second dosing with concentration of $D at day $tdose")
end

# ╔═╡ 0bfc00dd-cb88-4166-8cb0-2583882975f4
md"A strong second dose after about eight days can keep the population in check for a while. However, if we give a lower dose, the resistant bacteria will strongly dominate!"

# ╔═╡ 89798b17-a1f0-43ff-9b3f-42cae0d6020a
md"### Sensitivity analysis"

# ╔═╡ a0bff583-53ee-4a76-8dea-93711b93ec4c
md"We can explore the system further by performing a local sensitivity analysis. Let us explore two parameters: the effectivity of the antibiotics (given by parameter $K_s$, the concentration where it at half its maximal effectivty) and $a$, the fitness cost of the resistant bacteria. We only consider a single initial dose of antibiotics."

# ╔═╡ 1dc98bcc-c991-43f2-884e-ec6c20ee4554
begin
	tsteps = 0:0.1:30

	S((Ks, a)) = solve(remake(oprob, p=[:Ks=>Ks, :a=>a]), RK4(), saveat=tsteps)[:S]
	R((Ks, a)) = solve(remake(oprob, p=[:Ks=>Ks, :a=>a]), RK4(), saveat=tsteps)[:R]

	sens_S((Ks, a)) = ForwardDiff.jacobian(S, [Ks, a])  # sensitivity for susceptible bacteria
	sens_R((Ks, a)) = ForwardDiff.jacobian(R, [Ks, a])  # sensitivity for resistant bacteria

end

# ╔═╡ c321d1b8-e22b-4e28-adfc-dc2b20d28083
@bind Ks Slider(1.0:1.0:20, show_value=true, default=4)

# ╔═╡ 9f856137-5a46-4a6d-be12-aa38f9e20642
@bind a Slider(0:0.1:0.9, show_value=true, default=0.1)

# ╔═╡ 5600d0fa-80ef-40c2-ae12-570c348d1f14
md"Below is a simulation with the given parameters:"

# ╔═╡ 7c18a27d-a738-4d04-b712-d87fd91828af
plot(tsteps, [S([Ks, a]) R([Ks, a])], label=["S" "R"], xlab="t", title="Bacterial load with Ks=$Ks and a=$a")

# ╔═╡ 3e1c9cd4-aa4e-4f29-8788-45708b94406f
md"Next, we perform a sensitivity analysis for Ks and a, respectively."

# ╔═╡ 752a90ee-3e92-4426-b33c-9d65b8c465dc
plot(tsteps, [sens_S((Ks, a))[:,1] sens_R((Ks, a))[:,1]], label=["S" "R"], xlab="t", title="Sensitivity for Ks")

# ╔═╡ d45219d2-183e-4ca1-89bc-34d0eae1419b
md"We see here that $K_s$ greatly positively impacts both types of bacteria shortly after the the antibiotics has been removed. The greater $K_s$, the higher the antibiotics concentration needs to be to substantionally effect the bacterial density. The effect is the greatest for resistant bacteria. At longer time intervals, a small increase in $K_s$ give a positive effect on the susceptible bacteria and a negative effect on the resistant ones. Due to competition, this is a zero-sum game and if the susceptible bacteria are less harmed, they can easier take back a large share of the system."

# ╔═╡ d5287aba-93ab-440a-9cb5-97d879f5af9b
plot(tsteps, [sens_S((Ks, a))[:,2] sens_R((Ks, a))[:,2]], label=["S" "R"], xlab="t", title="Sensitivity for a")

# ╔═╡ c93ed717-dcfd-46f5-83c1-91d8eb7fa943
md"Above, we see the effect of the fitness cost $a$, the decrease in growth rate the resistant bacteria show. Though the fitness cost only directly impacts the resistant bacteria, both are impacted due to competition. Just after the initial boom when the antibiotics have dissipated, there is a large negative local minimum for the resistant bacteria. After that, we see a strong negative effect for the resistant bacteria and, conversely, a positive effect for the susceptible ones."

# ╔═╡ 60c63f1b-5a27-4539-a486-78c083457b0e
md"""
## Conclusion

This toy model illustrates the effect of antibiotic treatment on a mixed population of susceptible and resistant bacteria. It highlights the importance of a sufficiently strong initial dose to prevent resistant bacteria from dominating the population. The model also demonstrates that the timing of a second dose can significantly influence the outcome, with a well-timed second dose potentially keeping the bacterial population in check.

However, this model is a simplification of real-world scenarios. It assumes a simple first-order removal of antibiotics, whereas actual pharmacokinetics involve absorption, distribution, metabolism, and excretion. Incorporating these factors would enhance the model's accuracy. Additionally, the model could be improved by including a more realistic immune response, capable of eliminating bacteria at low concentrations. Future work could also explore the effects of multiple doses or continuous antibiotic infusions, as well as the impact of stochastic fluctuations in bacterial growth and antibiotic effects.
"""

# ╔═╡ 890afefc-f42b-4d74-b775-6dee5e5f0c2b
md"## Appendix"

# ╔═╡ Cell order:
# ╟─14c7e803-c0ff-4211-8b60-2c9c246934dd
# ╟─6948149f-854e-4b3c-b51c-099dd221ab83
# ╟─89551690-500d-4e37-ae20-5beb71cc87ac
# ╠═1e078c8b-2df9-4646-b629-bac7f151e937
# ╠═aed58771-86c1-4928-a1ba-b7d2f503b188
# ╟─33a6968f-d8df-4f56-baf3-6e0e88e54c23
# ╠═83b50ef8-9b68-428a-bf93-1583d032e9aa
# ╠═ccb31a3a-f75e-479d-bec7-75cf590652c5
# ╟─f284fda5-aaf7-4121-a8f0-b996d501bec5
# ╟─aaf3d52e-7e6b-4742-b9df-5440b148e7ee
# ╟─a2164307-6d18-4983-9aaa-25c8f5828655
# ╠═eab8adc0-e3c7-4693-990b-cb47e3d9c2cf
# ╟─556243a1-c387-4606-94d9-2dc06dfdc5f0
# ╠═98451eb8-8755-47fe-b8c6-ce2959256dd8
# ╠═b62c37ac-8d33-49ac-8aa6-18d959bc6933
# ╠═797ce7d7-c49d-4a00-a8e3-058b6de4d3a9
# ╠═aa32c5e2-0f67-4910-a612-d958fa13169c
# ╟─4679507b-e6e4-4d39-8610-77141f7b82e0
# ╟─49e97512-4150-49e2-bd9b-c3bfd7c375dc
# ╠═50794b21-e90f-438d-ad9c-eaddfab14ce7
# ╠═7f9c8c7d-90ef-496c-b907-300feba5edfc
# ╠═bf81ccd3-4cdf-4a9c-b334-f7339fed8dfe
# ╠═37ecd681-bd01-4162-94cf-97f6d8d82a8e
# ╠═438de4bb-6d15-4d76-b9a5-95b51924cbc2
# ╟─0d464337-8a70-47eb-8c08-cc2a7d388cff
# ╟─0bfc00dd-cb88-4166-8cb0-2583882975f4
# ╟─89798b17-a1f0-43ff-9b3f-42cae0d6020a
# ╠═158c4ad8-a38b-4828-812f-1abce60695dc
# ╠═6785c1af-28e5-46b4-ba34-a8ab65be7ff8
# ╟─a0bff583-53ee-4a76-8dea-93711b93ec4c
# ╟─1dc98bcc-c991-43f2-884e-ec6c20ee4554
# ╠═c321d1b8-e22b-4e28-adfc-dc2b20d28083
# ╠═9f856137-5a46-4a6d-be12-aa38f9e20642
# ╟─5600d0fa-80ef-40c2-ae12-570c348d1f14
# ╟─7c18a27d-a738-4d04-b712-d87fd91828af
# ╟─3e1c9cd4-aa4e-4f29-8788-45708b94406f
# ╟─752a90ee-3e92-4426-b33c-9d65b8c465dc
# ╟─d45219d2-183e-4ca1-89bc-34d0eae1419b
# ╟─d5287aba-93ab-440a-9cb5-97d879f5af9b
# ╟─c93ed717-dcfd-46f5-83c1-91d8eb7fa943
# ╟─60c63f1b-5a27-4539-a486-78c083457b0e
# ╟─890afefc-f42b-4d74-b775-6dee5e5f0c2b
# ╟─21357c48-f35d-11ee-23f8-2534bb1d82f4
