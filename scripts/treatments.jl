### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 6cf979c3-4cf8-4add-91f5-9f3999e4572f
using Pkg; Pkg.activate("..")

# ╔═╡ 4bb3199e-fcf8-11ef-2850-7dc8622013c6
using Turing, StatsPlots

# ╔═╡ e95bf2da-18b1-42d3-b527-d91c096301ff
using Optim

# ╔═╡ 2befdefd-5385-4617-a28f-826f8cc2b6e1
using StatsBase

# ╔═╡ 4320444b-feaf-4cd0-aac6-f26921117448
md"""

## Treatments

We are testing two potential treatments for COVID-19. Treatment A is effective for 28 in 30 people, treatment B is effective for 17 in 20 people. **Which treatment is more effective?**
"""

# ╔═╡ 6b50e9e2-f536-4a61-b1e0-7b1324c4f0e3
N_A₊ = 28  # recovered treatment A

# ╔═╡ e7cc1f89-b36f-462a-b708-01d24f2b2629
N_A₋ = 2  # not recovered treatment A

# ╔═╡ bc0c3e4e-4a88-4ebb-b98e-beaa7c4f0259
N_B₊ = 17  # recovered treatment B

# ╔═╡ 27c88c3a-c974-431d-98cf-fcab420d78be
N_B₋ = 3  # not recovered treatment B

# ╔═╡ a7118c6f-8af6-40ac-951d-f7d1547d6373
@model function treatment(N_A₊, N_A₋, N_B₊, N_B₋)
    # prior distributions of probability of success
    p_A ~ Uniform(0, 1)  # probability that A is effective
    p_B ~ Uniform(0, 1)  # probability that B is effective
    # number of cases for each treatment
    N_A = N_A₊ + N_A₋
    N_B = N_B₊ + N_B₋
    # sample successes
    N_A₊ ~ Binomial(N_A, p_A)
    N_B₊ ~ Binomial(N_B, p_B)
	return
end

# ╔═╡ ec62ca20-7335-469e-8acd-b3b27930f73a
chain = sample(treatment(N_A₊, N_A₋, N_B₊, N_B₋), NUTS(), 1_000);

# ╔═╡ 6c97eabc-8627-472c-b7e5-e7619e9a1bf1
summarize(chain)

# ╔═╡ 2736dbf9-79b5-422b-b33f-7b880acd9efc
plot(chain)

# ╔═╡ be6f849d-5edf-44e9-9c64-98f26b512a6a
quantile(chain)

# ╔═╡ abeaac18-31b2-4912-946f-bc7109127610
Δ = chain[:p_A] - chain[:p_B]

# ╔═╡ e4e2aa95-74d0-4f37-8a8c-dc1b81a4e4c2
histogram(Δ)

# ╔═╡ 162e05c6-23a7-4537-83a3-94e36a81faa2
p_A_groter_B = mean(Δ .> 0)

# ╔═╡ 10a56778-13d6-458a-b4ff-576c41ddac01
p_A_veel_groter_B = mean(Δ .> 0.10)

# ╔═╡ 86a0e46f-7376-44fd-ba18-2f6f989bcaa4
mle = optimize(treatment(N_A₊, N_A₋, N_B₊, N_B₋), MLE(), GradientDescent())

# ╔═╡ 02f13ba9-15b4-4d07-acc4-35288a99743f
coeftable(mle)

# ╔═╡ b0214284-e484-42a5-a82f-9c480bc96b1f


# ╔═╡ Cell order:
# ╠═6cf979c3-4cf8-4add-91f5-9f3999e4572f
# ╠═4bb3199e-fcf8-11ef-2850-7dc8622013c6
# ╟─4320444b-feaf-4cd0-aac6-f26921117448
# ╠═6b50e9e2-f536-4a61-b1e0-7b1324c4f0e3
# ╠═e7cc1f89-b36f-462a-b708-01d24f2b2629
# ╠═bc0c3e4e-4a88-4ebb-b98e-beaa7c4f0259
# ╠═27c88c3a-c974-431d-98cf-fcab420d78be
# ╠═a7118c6f-8af6-40ac-951d-f7d1547d6373
# ╠═ec62ca20-7335-469e-8acd-b3b27930f73a
# ╠═6c97eabc-8627-472c-b7e5-e7619e9a1bf1
# ╠═2736dbf9-79b5-422b-b33f-7b880acd9efc
# ╠═be6f849d-5edf-44e9-9c64-98f26b512a6a
# ╠═abeaac18-31b2-4912-946f-bc7109127610
# ╠═e4e2aa95-74d0-4f37-8a8c-dc1b81a4e4c2
# ╠═162e05c6-23a7-4537-83a3-94e36a81faa2
# ╠═10a56778-13d6-458a-b4ff-576c41ddac01
# ╠═e95bf2da-18b1-42d3-b527-d91c096301ff
# ╠═2befdefd-5385-4617-a28f-826f8cc2b6e1
# ╠═86a0e46f-7376-44fd-ba18-2f6f989bcaa4
# ╠═02f13ba9-15b4-4d07-acc4-35288a99743f
# ╠═b0214284-e484-42a5-a82f-9c480bc96b1f
