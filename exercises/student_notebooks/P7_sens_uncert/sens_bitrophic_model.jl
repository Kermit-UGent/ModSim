### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ a3c1055e-b44f-4412-9ca0-f8ec5f972494
begin
	# add this cell if you want the notebook to use the environment from where the Pluto server is launched
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ fd357fe0-0920-11ef-02e4-25a84575c6a2
using Markdown

# ╔═╡ 0ddc1082-bf8d-4357-a056-2fa861824b47
using InteractiveUtils

# ╔═╡ cce8fcde-069f-49b9-aadd-6abccbcfbace
using ModelingToolkit

# ╔═╡ ab75267c-d28f-47ce-88bb-5aeb885398e4
using ModelingToolkit: t_nounits as t, D_nounits as D

# ╔═╡ 88aefc2c-5114-48ff-843b-e0d8f08ba3b5
using Turing

# ╔═╡ b48e4b46-3a07-4ff0-b617-6bc1125f6a1e
using OrdinaryDiffEq, StatsPlots

# ╔═╡ ddc035a2-3b03-4a54-807c-05a271384ca7
using ForwardDiff

# ╔═╡ aaa54f95-466f-463e-958a-fc02560fd7d1
using GlobalSensitivity

# ╔═╡ 2293569d-491b-4b7f-a065-0c17d43b56d5
md"""
# Hermione's PhD on plant insect interactions
"""


# ╔═╡ 9bc7b562-e1d2-4395-beb6-92f413589928
md"""
During her PhD Hermione has studied in her first year the plant Cucurbita pepo in a greenhouse. She determined the following parameters.
- the carrying capacity $k$ of her experimental setup is 4000 kg/ha
- the crop growth rate θ is 0.2 $d^{-1}$

In the second year she studied the leaf-cutter ant (Atta cephalotes) and the effect of a neurotoxin on the ant from which she determined the following parameters. 
- The mortality rate μ equals 0.1 $d^{-1}$
- the insecticide effect factor $p$ equals 3

In her last year, she is interested in the interactions between Cucurbita pepo and the leaf-cutter ant in her experimental setup. She knows that the relationship between a field crop (C) and a voracious insect (A) is described by the bitrophic model:

$$\begin{eqnarray*}
\frac{dC}{dt} &= \theta C \left(1-\frac{C}{k}\right)-fCA \\
\frac{dA}{dt} &= \phi f CA -(1 + p)\, \mu A
\end{eqnarray*}$$

| Parameter | Value | Description |
| :--- | :---: | :--- |
| $\theta$ | $0.2$ | Growth rate of $C$ |
| $k$ | $4000$ | Carrying capacity |
| $p$ | $3$ | Insecticide effect |
| $\mu$ | $0.1$ | Natural mortality rate |
| $f$ | **?** | feeding rate |
| $\phi$ | **?** | Conversion efficiency |

Based on a literature review she retrieves some initial ideas of the distributions of the feeding rate $f$ and the conversion $ϕ$.
- The feeding rate $f$ is LogNormal distributed around $f=0.001\;ha/(kg\,d)$ as `LogNormal(log(0.001), 0.5)`. In literature the values for $f$ are in the range of `[0.0005, 0.002]`.
- The conversion efficiency $ϕ$ is somewhere between 0.07 and 0.3 with a peak around 0.15.

**Her final goal is to be able to predict the pumpkin harvest in natural habitats with ants after 200 days.** She wants to know if in natural habitats the yield can reach over 2000 kg/ha, because that makes it profitable for farmers.
"""


# ╔═╡ 07495e3b-9b66-4f58-9b23-db461751e3ca
md"""
## Roadmap of Hermione's third year
"""


# ╔═╡ 39eeb5cf-2c10-40d3-a628-4778d593776a
md"""
1. Gain insights in the overall sensitivity of the model. What parameters have a major influence on the pumpkin harvest?
2. Determine the influence of the unknown parameters on the harvested pumpkin. What are the chances a farmer has less than 2000 kg/ha yield?
3. After the first year, the farmer reports a yield of 2400 kg/ha. Can you a priori know which parameter ($ϕ$ or $f$) will be better estimated using this extra measurement? Secondly use this measurement to update your beliefs about these parameters.
"""


# ╔═╡ 60a25144-e430-45cf-b250-3d196bdb4723
md"""
!!! tip "Step 0"
	Implement the model, use the default values `f=0.001`and `ϕ=0.15`. The initial crop density is `100 kg/ha` and the intial ant population is `0.5 kg/ha`. We want to simulate the harvest after `200 days`.
"""


# ╔═╡ fdd9a533-c441-4b77-bba8-687bddf3e4cb


# ╔═╡ f741a2eb-48dc-4ff4-817e-849a6dd839d3
md"""
## Part I: What determines the pumpkin harvest?
"""


# ╔═╡ ed3c3b7e-bb51-442b-a455-eebff351b227
md"""
!!! tip "Question 1"
    Gain insights in the overall sensitivity of the model. What is the potential influence of $f$ and $ϕ$ on the final harvest?
"""



# ╔═╡ 592a9822-f8e9-40ed-acdf-ee93975b6784


# ╔═╡ e296e434-017c-45ce-bfe0-76601e277612
md"""
## Part II: How reliable/variable are the outcomes?
"""


# ╔═╡ 8f8c1f5d-a5c5-4da0-8230-d81b9ce7dab5
md"""
!!! tip "Question 2"
	Determine the influence of the unknown parameters on the harvested pumpkin. What are the chances a farmer has less than 2000 kg/ha yield?
"""


# ╔═╡ e9f132c6-4828-43bc-9c7b-7b93515d4322


# ╔═╡ af93263f-fda7-4cad-8f3d-74ab4a5b345d
md"""
## Part III: Can we reduce uncertainty after year one?
"""


# ╔═╡ fe08b97c-08ad-4a19-8272-b047bc652ec6
md"""
!!! tip "Question 3"
	After the first year, the farmer reports a yield of 2400 kg/ha. Can you a priori know which parameter ($ϕ$ or $f$) will be better estimated using this extra measurement? Secondly use this measurement to update your beliefs about these parameters. What is the updated chance that the farmer will have less than 2000 kg/ha yield?

	Make use of the Metropolis-Hastings algorithm and use 5000 samples. You can assume a standard deviation of 50 kg/ha on the measured harvest yield.
"""


# ╔═╡ 81aa36db-b28f-4dff-9513-cfc3b80daaab


# ╔═╡ Cell order:
# ╠═fd357fe0-0920-11ef-02e4-25a84575c6a2
# ╠═0ddc1082-bf8d-4357-a056-2fa861824b47
# ╠═a3c1055e-b44f-4412-9ca0-f8ec5f972494
# ╠═cce8fcde-069f-49b9-aadd-6abccbcfbace
# ╠═ab75267c-d28f-47ce-88bb-5aeb885398e4
# ╠═88aefc2c-5114-48ff-843b-e0d8f08ba3b5
# ╠═b48e4b46-3a07-4ff0-b617-6bc1125f6a1e
# ╠═ddc035a2-3b03-4a54-807c-05a271384ca7
# ╠═aaa54f95-466f-463e-958a-fc02560fd7d1
# ╟─2293569d-491b-4b7f-a065-0c17d43b56d5
# ╟─9bc7b562-e1d2-4395-beb6-92f413589928
# ╟─07495e3b-9b66-4f58-9b23-db461751e3ca
# ╟─39eeb5cf-2c10-40d3-a628-4778d593776a
# ╟─60a25144-e430-45cf-b250-3d196bdb4723
# ╠═fdd9a533-c441-4b77-bba8-687bddf3e4cb
# ╟─f741a2eb-48dc-4ff4-817e-849a6dd839d3
# ╟─ed3c3b7e-bb51-442b-a455-eebff351b227
# ╠═592a9822-f8e9-40ed-acdf-ee93975b6784
# ╟─e296e434-017c-45ce-bfe0-76601e277612
# ╟─8f8c1f5d-a5c5-4da0-8230-d81b9ce7dab5
# ╠═e9f132c6-4828-43bc-9c7b-7b93515d4322
# ╟─af93263f-fda7-4cad-8f3d-74ab4a5b345d
# ╟─fe08b97c-08ad-4a19-8272-b047bc652ec6
# ╠═81aa36db-b28f-4dff-9513-cfc3b80daaab
