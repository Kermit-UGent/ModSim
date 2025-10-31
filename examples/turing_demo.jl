### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 60e456e0-f5d3-11ef-1905-3f48fb33b4b3
using Pkg; Pkg.activate("..")

# ╔═╡ e8fd6780-8011-490b-99b8-c84b30db053c
using Turing, StatsPlots

# ╔═╡ ad1ec40c-0491-4585-a402-ad3c4eca89ce
plot(Beta(25, 10))

# ╔═╡ 70ea121e-dbe1-457c-9870-3e076c066986
plot(Poisson(80))

# ╔═╡ cfa77a91-4577-4009-9567-2dff97b30004
md"""

# Gentle Turing intro

We want to generate a sample of a biomial distribution `X ~ Binom(N, P)`, though both the number of trials $N$ and the probability of success follow a distribution (a *prior*):

`P ~ Beta(25, 10)`

`N ~ Poisson(100)`

Sample 10,000 observations from this distribution and draw a histogram. Compute, based on the sample:
-  $E[X]$
-  $E[\log(X + 1)]$
-  $P(X \ge 60)$
-  $P( 60 \le X < 90)$
-  $E[X \mid P > 0.8]$
-  $P(X \ge 60 \mid N < 80)$ 
"""

# ╔═╡ cd7f7173-1003-4c1c-9f51-96e7c01dd446
@model function mijnmodel()
	P ~ Beta(25, 10)
	N ~ Poisson(100)
	X ~ Binomial(N, P)
end

# ╔═╡ 46681f14-831c-4cd5-b4a3-5a699529cda9
[rand(mijnmodel()) for i in 1:100]

# ╔═╡ 0cece9dc-ee58-423c-a7b6-43ad7fb52e18
steekproef = sample(mijnmodel(), Prior(), 10000)

# ╔═╡ 3f7dfd1d-ac20-4cf6-8b72-f8e87999c62f
summarize(steekproef) # summarize

# ╔═╡ 4b0f63d5-a626-4500-a983-4e95b69826a8
std(Beta(25, 10))

# ╔═╡ bd69dbd1-85b0-4775-b205-481c2d4b6bae
quantile(steekproef) # quantile

# ╔═╡ de06b284-4016-4e38-879f-0751bf2147b2
Ps = steekproef[:P]

# ╔═╡ e0b3ddb9-3b8d-4f88-b1f2-2b4c6766153d
Ns = steekproef[:N]

# ╔═╡ d211be44-34ba-489a-a793-93ade8f2c3d5
Xs = steekproef[:X]

# ╔═╡ 3c56da32-d278-4863-8081-0be4bac4d821
#histogram?
histogram(Xs)

# ╔═╡ 7ea78e6d-a98d-4fc7-81c4-35a211699b49
md"$E[X]$"

# ╔═╡ bcda9d4c-4b5f-4bcb-9404-c79765d0b9b3
mean(Xs)

# ╔═╡ ba489dd3-1156-4d7d-8694-9e6477bf7346
md"$E[\log(X + 1)]$"

# ╔═╡ ee773a01-2404-462b-82ad-48169c3ba2ed
mean(log.(Xs .+ 1))

# ╔═╡ fec6059f-a84a-4971-bdf4-a3fa73c1567c
md"$P(X \ge 60)$"

# ╔═╡ 9a4a1360-7b1b-4859-9adf-f38f59adde9d
Xs .>= 60

# ╔═╡ d74c65f9-e092-4199-9a35-ec23184f9592
mean(Xs .>= 60)

# ╔═╡ d450f129-e0b5-41e8-be41-70e28ced85f4
md"$P( 60 \le X < 90)$"

# ╔═╡ 53f4d609-5083-478f-8c0e-29000b928e0f
mean(60 .≤ Xs .< 90)

# ╔═╡ fd5e5de0-35c5-4b96-9d72-31d5134e5772
md"$E[X \mid P > 0.8]$"

# ╔═╡ 9995f935-22e3-492a-8747-9262ca3c17fb
mean(Xs[Ps .> 0.8])

# ╔═╡ 863fdcf9-cba6-4b51-b6cd-84e0206a8477
md"$P(X \ge 60 \mid N < 80)$"

# ╔═╡ dd71ebce-3687-4108-8373-b16617883879
mean(Xs[Ns .< 80] .≥ 60)

# ╔═╡ 2c688ce5-a657-4e81-b932-77cfaa44c9dd
mijnmodel() | (X=89,)

# ╔═╡ 76a3fb6a-2d68-4f9f-9d1f-c14ae82c1eee


# ╔═╡ d4177abc-549a-4760-a7db-6061e9fcb176


# ╔═╡ 5776a18d-0f76-4277-8204-6408bfe5c5cf


# ╔═╡ Cell order:
# ╠═60e456e0-f5d3-11ef-1905-3f48fb33b4b3
# ╠═e8fd6780-8011-490b-99b8-c84b30db053c
# ╠═ad1ec40c-0491-4585-a402-ad3c4eca89ce
# ╠═70ea121e-dbe1-457c-9870-3e076c066986
# ╟─cfa77a91-4577-4009-9567-2dff97b30004
# ╠═cd7f7173-1003-4c1c-9f51-96e7c01dd446
# ╠═46681f14-831c-4cd5-b4a3-5a699529cda9
# ╠═0cece9dc-ee58-423c-a7b6-43ad7fb52e18
# ╠═3f7dfd1d-ac20-4cf6-8b72-f8e87999c62f
# ╠═4b0f63d5-a626-4500-a983-4e95b69826a8
# ╠═bd69dbd1-85b0-4775-b205-481c2d4b6bae
# ╠═de06b284-4016-4e38-879f-0751bf2147b2
# ╠═e0b3ddb9-3b8d-4f88-b1f2-2b4c6766153d
# ╠═d211be44-34ba-489a-a793-93ade8f2c3d5
# ╠═3c56da32-d278-4863-8081-0be4bac4d821
# ╟─7ea78e6d-a98d-4fc7-81c4-35a211699b49
# ╠═bcda9d4c-4b5f-4bcb-9404-c79765d0b9b3
# ╟─ba489dd3-1156-4d7d-8694-9e6477bf7346
# ╠═ee773a01-2404-462b-82ad-48169c3ba2ed
# ╟─fec6059f-a84a-4971-bdf4-a3fa73c1567c
# ╠═9a4a1360-7b1b-4859-9adf-f38f59adde9d
# ╠═d74c65f9-e092-4199-9a35-ec23184f9592
# ╟─d450f129-e0b5-41e8-be41-70e28ced85f4
# ╠═53f4d609-5083-478f-8c0e-29000b928e0f
# ╟─fd5e5de0-35c5-4b96-9d72-31d5134e5772
# ╠═9995f935-22e3-492a-8747-9262ca3c17fb
# ╟─863fdcf9-cba6-4b51-b6cd-84e0206a8477
# ╠═dd71ebce-3687-4108-8373-b16617883879
# ╠═2c688ce5-a657-4e81-b932-77cfaa44c9dd
# ╠═76a3fb6a-2d68-4f9f-9d1f-c14ae82c1eee
# ╠═d4177abc-549a-4760-a7db-6061e9fcb176
# ╠═5776a18d-0f76-4277-8204-6408bfe5c5cf
