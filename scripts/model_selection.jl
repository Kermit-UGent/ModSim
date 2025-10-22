### A Pluto.jl notebook ###
# v0.20.3

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

# ╔═╡ 80bc6542-2f0a-11ef-3308-8b1269e840a1
begin
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ 04b15fa3-536f-4676-a41e-f1c8147a2c97
using PlutoUI, Plots, StatsPlots

# ╔═╡ 39866945-8759-4034-acdb-0a34a6a26b6d
using Turing, Optim, StatsBase, LaTeXStrings

# ╔═╡ b1653a62-2456-4ec5-bc4d-f0c6be776f1e
using Distributions

# ╔═╡ 37608790-16a3-471d-bf37-c10a9653ad50
using LinearAlgebra

# ╔═╡ 0376ff93-2c87-410a-8cf8-65d3c654b503
md"# Model Selection"

# ╔═╡ 2e4beb9f-068d-498d-99bc-eb7b348a544a
md"## Growth model"

# ╔═╡ 3fa5346f-9e04-4406-b3ed-9e2882a6b3ac
solar_power = [488, 560, 1770, 2115, 2352, 2883, 3045, 3085, 3149, 3563, 3528, 4259, 4678, 6413, 7193]  # in GWh, Belgian solar power production

# ╔═╡ 9212be24-0d15-4455-aa0d-eff226737cf9
years = 2009:2023

# ╔═╡ 26ddc2b4-591c-44a9-9373-e971cc2218a9
t, y = 0:length(solar_power)-1, solar_power

# ╔═╡ efc29b29-45e0-44a7-9ab3-168cc4a8a0d8
scatter(t, y)

# ╔═╡ 14d826a6-5506-4557-b8bc-85a70f961ea6
@model function linear(t, y)
	a ~ Turing.Flat()
	b ~  Turing.Flat()
	σsq ~ Turing.FlatPos(0)
	f = x -> a * x + b
	for i in 1:length(y)
		y[i] ~ Normal(f(t[i]), √(σsq))
	end
	return f
end

# ╔═╡ 6af89987-e8f6-4155-867a-38cafbdd113f
@model function quadratic(t, y)
	a ~ Turing.Flat()
	b ~ Turing.Flat()
	c ~ Turing.Flat()
	σsq ~ Turing.FlatPos(0)
	f = x -> a * x^2 + b * x + c
	for i in 1:length(y)
		y[i] ~ Normal(f(t[i]), √(σsq))
	end
	return f
end

# ╔═╡ 804e5e78-1831-4c04-a489-86eb7e8b23d1
@model function exponential(t, y)
	C ~ Turing.FlatPos(0)
	k ~  Turing.FlatPos(0)
	σsq ~ Turing.FlatPos(0)
	tmin = minimum(t)
	f = x -> C * exp(k * x)
	for i in 1:length(y)
		y[i] ~ Normal(f(t[i]), √(σsq))
	end
	return f
end

# ╔═╡ bb67020b-ce43-46d3-9dde-8fb35fcde7f6
@model function gompertz(t, y)
	a ~ FlatPos(0)
	b ~ FlatPos(0)
	c ~ FlatPos(0)
	σsq ~ Turing.FlatPos(0)
	f = x -> a * exp(-b * exp(-c * x))
	for i in 1:length(y)
		y[i] ~ Normal(f(t[i]), √(σsq))
	end
	return f
end

# ╔═╡ 3c67dadb-875b-4398-bcc4-db53b4eaf813
@model function verhulst(t, y)
	K ~ FlatPos(0)
	P0 ~ FlatPos(0)
	r ~ FlatPos(0)
	σsq ~ Turing.FlatPos(0)
	f = x -> K * P0 * exp(r*x) / (K + P0 * (exp(r*x) - 1))
	for i in 1:length(y)
		y[i] ~ Normal(f(t[i]), √(σsq))
	end
	return f
end

# ╔═╡ 4a0d5500-dec5-46c9-9009-44fc7e2f2366
@model function unstable_exponential(t, y)
	C ~ Turing.FlatPos(0)
	k ~  Turing.FlatPos(0)
	t0 ~ Turing.Flat()
	σsq ~ Turing.FlatPos(0)
	tmin = minimum(t)
	f = x -> C * (x-t0) * exp(-k * (x-t0))
	for i in 1:length(y)
		y[i] ~ Normal(f(t[i]), √(σsq))
	end
	return f
end

# ╔═╡ caaf6585-d370-4c81-a9fe-97eb0f3fdd66
model_zoo = Dict(
	"linear" => linear(t, y),
	"quadratic" => quadratic(t, y),
	"exponential" => exponential(t, y),
	"Gompertz" => gompertz(t, y),
	"Verhulst" => verhulst(t, y),
	"unstable exp." => unstable_exponential(t, y),
)

# ╔═╡ d3ca79db-e156-4eac-8175-aafa80c74098
logistic = t -> 13 * 2 * exp(1*t) / (13 + 2 * (exp(1*t) - 1))

# ╔═╡ e564f79f-bba1-499b-adcc-495d6dc9ece0
model_linear = linear(t, y)

# ╔═╡ a8c66ff0-c401-4308-9966-dccb570289dd
model_quadratic = quadratic(t, y)

# ╔═╡ ece8ac55-8c1e-4a9e-adfb-fd5ccb48bd99
model_exponential = exponential(t, y)

# ╔═╡ cd182100-51f9-4deb-8911-741975fde98b
model_gompertz = gompertz(t, y)

# ╔═╡ a1732042-4c15-49a5-9218-ab9924dbd5b5
model_verhulst = verhulst(t, y)

# ╔═╡ c5ce4945-b106-4b02-818d-4fd6993e7a1c
mle_verhulst = optimize(model_verhulst, MLE())

# ╔═╡ 083ec9d8-d5b8-4c12-8254-aab4c323af2e
coeftable(mle_verhulst)

# ╔═╡ 87218738-418d-4c70-bc54-b6796524affe
mle_exp = optimize(model_exponential, MLE(), NelderMead())

# ╔═╡ 87a4674c-09f5-47ca-a0b3-64a130595ba6
coeftable(mle_exp)

# ╔═╡ 65dfbd5c-af8a-40fd-b2e3-4841357efb58
model_unstab_exp = unstable_exponential(t, y)

# ╔═╡ e8054391-6151-4525-8bc9-2ce5ad4e9c28
optimize(model_verhulst, MLE())

# ╔═╡ 654f27e6-accd-4f70-b742-7b1ff80675be
optimize(model_quadratic, MLE())

# ╔═╡ 98af497f-33e6-480f-9914-97aff274a002
optimize(linear(t, y), MAP(), LBFGS())

# ╔═╡ 6beb4f72-3e3c-40c3-a047-55ee0cbb5881
mle_lin = optimize(linear(t, y), MLE(), BFGS())

# ╔═╡ 430afbe8-e576-4ed1-beae-cefe3f79d6a7
mle_lin.lp

# ╔═╡ 4b0eb819-4788-4ef0-a7f4-5fcb948b131a
params(mle_lin) .=> mle_lin.values

# ╔═╡ 96fdb096-e7fa-4721-8249-a63e87ab1244
NamedTuple(zip(params(mle_lin), mle_lin.values))

# ╔═╡ 7cba7c60-99ac-4b13-bdeb-e13487b0b6b6
f = generated_quantities(mle_lin.f.model, NamedTuple(zip(params(mle_lin), mle_lin.values)))

# ╔═╡ ee8980b5-528b-465d-855b-ad17c01ac6c5
pars = mle_lin.values

# ╔═╡ 8b89da5e-e842-4491-882a-d07cd04c7ecf
generated_quantities

# ╔═╡ ca8a9bad-26dd-4d08-a9f5-ca9ca683bc71
md" ## Coin example"

# ╔═╡ f61522ad-00a7-447a-a249-0d387050248c
n = 100

# ╔═╡ ddbf180d-f045-4ad3-88e1-c4ff9fd80ee0
pdf(Binomial(100, 0.5), 61)

# ╔═╡ 4a37e039-db45-4688-8d29-b1662218db60
x = 61

# ╔═╡ 0f9dde78-9921-4bae-a4a7-e784275b82d6
M1 = Binomial(n, 0.5)

# ╔═╡ 74ebbfe5-2e5f-4505-b5f5-8df50058da82
M2 = Binomial(n, 0.6)

# ╔═╡ f2eddd98-eaf1-4a7e-99a8-c05a34efe89d
pM1 = pdf(M1, x)

# ╔═╡ 52f00679-3529-410e-ad6c-dd9a36414d6e
pM2 = pdf(M2, x)

# ╔═╡ a5a02d4e-b2d9-4a2f-aa1b-5c6811d6fb6b
KA = pM2 / pM1

# ╔═╡ ba025f0f-e514-413e-8e7d-917bd64e47fe
prior_p = Beta(14, 10)

# ╔═╡ 8d39ea49-f38d-4be4-90e3-53078652ca45
dp = 0.01 # step size

# ╔═╡ cc250c89-5c4b-49be-8637-7763e4869605
# Riemann sum
pM3 = sum(p->pdf(Binomial(n, p), x) * pdf(Beta(14, 10), p) * dp, 0:dp:1)

# ╔═╡ 657dc886-07be-496c-80aa-848565fef4c3
KB = pM3 / pM1

# ╔═╡ bb3d17ad-139a-4906-8279-d12cac60a1da
md"## Fitting"

# ╔═╡ d29ee2bd-85c9-4e33-b434-fbe4346ec1bf
degree = 7

# ╔═╡ 2d0d98a3-cecc-46d0-ac39-db299cf58cfd
σ = 3

# ╔═╡ d4c49a72-e818-4836-8bca-d3cbc354baa5
N = 50

# ╔═╡ a8490458-7c12-489b-8d19-4ef3fd129e9f
xvals = 0:1:20

# ╔═╡ 3c9861b5-b73d-45ba-87ab-83d2d8b60a8d
xdense = 0:0.1:20

# ╔═╡ baec7e94-71eb-4e3e-ad39-6404316f6522
Y = xvals .* sin.(0.5xvals) .+ σ.* randn(length(xvals), N)

# ╔═╡ 6e4fb390-ca56-488e-8066-e22d058f89a7
B = [xi^d / factorial(d) for xi in xvals, d in 0:degree]

# ╔═╡ 3bfe5ba7-9044-4bbe-a2b8-2fbf7b9e8bd0
Blarge = [xi^d / factorial(d) for xi in xdense, d in 0:degree]

# ╔═╡ b2f98194-2255-4d90-9869-0aed381e5cac
Ŷ = Blarge * ((B'*B) \ B') * Y

# ╔═╡ 4d3d9d03-23bc-4b67-8fc6-b84ea1618360
Eŷ = mean(Ŷ, dims=2)[:,1]

# ╔═╡ 3be0872f-5e8d-4bbd-9161-1c3ec6014c93
@bind i Clock(max_value=N, repeat=true)

# ╔═╡ 8d9b8e11-91eb-47a4-beb9-d837e0f26161
let
	y = Y[:,i]
	ŷ = Ŷ[:,i]
	p = plot(xlab=L"x", ylab=L"y", title="Polynomial fit\n (degree=$degree)")
	scatter!(p, xvals, y, label="data", ylims=extrema(Y), color="yellow")
	for j in 1:N
		plot!(p, xdense, Ŷ[:,j], alpha=0.2, color="grey", label="", lw=0.5)
	end
	plot!(p, xdense, ŷ, alpha=0.7, color="blue", label="ŷ(x)", lw=2)
	plot!(p, xdense, Eŷ, alpha=0.7, color="red", label="E[ŷ(x)]", lw=2)
	
	p
end

# ╔═╡ 370ebfc9-d70d-4f5f-b55c-a9a1cebc7ae3
begin
	plot()
end

# ╔═╡ 45720d9f-3c25-45fb-b574-a4d4e0fe968a
poly_data = polymodel(xvals, yvals, degree)

# ╔═╡ 7dbb3666-fdfd-4f2e-96be-bce9922cf7b2
poly_map = optimize(poly_data, MAP())

# ╔═╡ 3ffd2e3c-e6ef-4674-b4f2-63b4336c59dd
chain_poly = sample(poly_data, NUTS(), 500);

# ╔═╡ 9c4fcf78-bec5-482d-9dfa-3891e22dfc2d
summarize(chain_poly)

# ╔═╡ 6c3699cc-2c8b-4bf9-a637-0843a3058800
poly_fits = generated_quantities(poly_data, chain_poly)

# ╔═╡ 3b69dd32-e549-4fa0-bae9-c8fe4ffd5716
plot(chain_poly)

# ╔═╡ 35b75da1-97c8-4df7-bdcd-0ee0e1a428ca
md"## Appendix"

# ╔═╡ 56f2b018-de9f-4ad2-abfd-920dc2169983
TableOfContents()

# ╔═╡ 31127a8d-365b-48c4-bfd8-2f968a72be09
plots = Dict()

# ╔═╡ 1fc1e71e-1a68-48ad-a8ed-f4d0966a686c
plots["solar_data"] = scatter(years, solar_power, label="", title="Annual photovoltaic power generation in Belgium", xlab="year", ylab="Production [GWh]", color="gold")

# ╔═╡ 51f46a70-99e7-4d6a-b095-8acfcefb962d
let
	p = plot(x->2x+1, 0, 5, lw=2, label="linear", xticks=[], yticks=[], xlab="t", title="Models for energy growth")
	plot!(x->0.7x^2-0.2x+1, 0, 5, lw=2, ls=:auto, label="quadratic")
	plot!(t->exp(t/2), 0, 5, lw=2, ls=:auto, label="exponential")
	plot!(t->10exp(-2exp(-2t)), 0, 5, lw=2, ls=:auto, label="Gompertz")
	plot!(logistic, 0, 5, lw=2, ls=:auto, label="Verhulst")
	plot!(t->10(t+1)*exp(-(t+1)/2), 0, 5, lw=2, ls=:auto, label="unstable exp.")
	plots["models"] = p
end

# ╔═╡ d4952d59-fcbd-4408-bbe8-ade98a2806a7
plots["prior_p"] = plot(p->pdf(prior_p, p), 0, 1, title="Prior on p", xlab=L"p", label="Beta(14, 10)", legend=:left, lw=2)

# ╔═╡ 6c1f29a8-936e-40d3-86fa-ad8424749aa4
# returns the model (generated quantities) of a MLE/MAP estimate of a model
get_model(mle) = generated_quantities(mle.f.model, NamedTuple(zip(params(mle), mle.values)))

# ╔═╡ bf79061b-5b34-4e0f-a452-a06de8ae5132
begin
	p = scatter(years, solar_power, label="", title="Annual photovoltaic power generation in Belgium", xlab="year", ylab="Production [GWh]", color="gold")

	aics = Dict{String,Float64}()
	bics  = Dict{String,Float64}()
	
	for (name, model) in model_zoo
		println(name * " :")
		mle = optimize(model, MLE())
		println("\tlog-likelihood = $(mle.lp)")
		map = optimize(model, MAP())
		println("\tlog-posterior = $(map.lp)")

		k = length(mle.values)

		aic = 2k - 2mle.lp
		bic = k * log(length(solar_power)) - 2mle.lp
		println("\tAIC : $aic")
		println("\tBIC : $bic")
		println("\tsigma : $(sqrt(last(mle.values)))")
		println("\tk : $k")

		aics[name] = aic
		bics[name] = bic

		f_mod = get_model(mle)
		tsteps = 0:0.1:20
		preds = f_mod.(tsteps)
		years_forecast = tsteps .+ minimum(years)
		plot!(p, years_forecast, preds, lw=2, alpha=0.6, ls=:auto, label=name)
	end
	plots["solar_forecast"] = p
		
end

# ╔═╡ bb021c5f-8157-4c94-a385-ee1c7a1dc5cb
begin
	aic_min = minimum(values(aics))
	p_model = Dict(name=>exp((aic_min-aic)/2) for (name, aic) in aics)
	ptot = sum(values(p_model))
	for (n, p) in p_model
		p_model[n] /= ptot
	end
	p_model
end

# ╔═╡ 59038a15-0ca8-4f77-a7d4-1840b80771ed
for (n, p) in p_model
	println("$n : $p")
end

# ╔═╡ 8b7b9fab-2725-4e0c-af68-0687ec5921c1
plots

# ╔═╡ 1eb4700c-20de-47a5-951d-3900d99461be
# ╠═╡ disabled = true
#=╠═╡
for (n, p) in plots
	savefig(p, "../figures/model_selection/$n.pdf")
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═04b15fa3-536f-4676-a41e-f1c8147a2c97
# ╠═39866945-8759-4034-acdb-0a34a6a26b6d
# ╠═80bc6542-2f0a-11ef-3308-8b1269e840a1
# ╠═0376ff93-2c87-410a-8cf8-65d3c654b503
# ╠═2e4beb9f-068d-498d-99bc-eb7b348a544a
# ╠═3fa5346f-9e04-4406-b3ed-9e2882a6b3ac
# ╠═9212be24-0d15-4455-aa0d-eff226737cf9
# ╠═1fc1e71e-1a68-48ad-a8ed-f4d0966a686c
# ╠═26ddc2b4-591c-44a9-9373-e971cc2218a9
# ╠═efc29b29-45e0-44a7-9ab3-168cc4a8a0d8
# ╠═14d826a6-5506-4557-b8bc-85a70f961ea6
# ╠═6af89987-e8f6-4155-867a-38cafbdd113f
# ╠═804e5e78-1831-4c04-a489-86eb7e8b23d1
# ╠═bb67020b-ce43-46d3-9dde-8fb35fcde7f6
# ╠═3c67dadb-875b-4398-bcc4-db53b4eaf813
# ╠═4a0d5500-dec5-46c9-9009-44fc7e2f2366
# ╠═caaf6585-d370-4c81-a9fe-97eb0f3fdd66
# ╠═51f46a70-99e7-4d6a-b095-8acfcefb962d
# ╠═d3ca79db-e156-4eac-8175-aafa80c74098
# ╠═bf79061b-5b34-4e0f-a452-a06de8ae5132
# ╠═bb021c5f-8157-4c94-a385-ee1c7a1dc5cb
# ╠═59038a15-0ca8-4f77-a7d4-1840b80771ed
# ╠═e564f79f-bba1-499b-adcc-495d6dc9ece0
# ╠═a8c66ff0-c401-4308-9966-dccb570289dd
# ╠═ece8ac55-8c1e-4a9e-adfb-fd5ccb48bd99
# ╠═cd182100-51f9-4deb-8911-741975fde98b
# ╠═a1732042-4c15-49a5-9218-ab9924dbd5b5
# ╠═c5ce4945-b106-4b02-818d-4fd6993e7a1c
# ╠═083ec9d8-d5b8-4c12-8254-aab4c323af2e
# ╠═87218738-418d-4c70-bc54-b6796524affe
# ╠═87a4674c-09f5-47ca-a0b3-64a130595ba6
# ╠═65dfbd5c-af8a-40fd-b2e3-4841357efb58
# ╠═e8054391-6151-4525-8bc9-2ce5ad4e9c28
# ╠═654f27e6-accd-4f70-b742-7b1ff80675be
# ╠═98af497f-33e6-480f-9914-97aff274a002
# ╠═6beb4f72-3e3c-40c3-a047-55ee0cbb5881
# ╠═430afbe8-e576-4ed1-beae-cefe3f79d6a7
# ╠═4b0eb819-4788-4ef0-a7f4-5fcb948b131a
# ╠═96fdb096-e7fa-4721-8249-a63e87ab1244
# ╠═7cba7c60-99ac-4b13-bdeb-e13487b0b6b6
# ╠═ee8980b5-528b-465d-855b-ad17c01ac6c5
# ╠═8b89da5e-e842-4491-882a-d07cd04c7ecf
# ╠═ca8a9bad-26dd-4d08-a9f5-ca9ca683bc71
# ╠═b1653a62-2456-4ec5-bc4d-f0c6be776f1e
# ╠═f61522ad-00a7-447a-a249-0d387050248c
# ╠═ddbf180d-f045-4ad3-88e1-c4ff9fd80ee0
# ╠═4a37e039-db45-4688-8d29-b1662218db60
# ╠═0f9dde78-9921-4bae-a4a7-e784275b82d6
# ╠═74ebbfe5-2e5f-4505-b5f5-8df50058da82
# ╠═f2eddd98-eaf1-4a7e-99a8-c05a34efe89d
# ╠═52f00679-3529-410e-ad6c-dd9a36414d6e
# ╠═a5a02d4e-b2d9-4a2f-aa1b-5c6811d6fb6b
# ╠═ba025f0f-e514-413e-8e7d-917bd64e47fe
# ╠═d4952d59-fcbd-4408-bbe8-ade98a2806a7
# ╠═8d39ea49-f38d-4be4-90e3-53078652ca45
# ╠═cc250c89-5c4b-49be-8637-7763e4869605
# ╠═657dc886-07be-496c-80aa-848565fef4c3
# ╠═bb3d17ad-139a-4906-8279-d12cac60a1da
# ╠═37608790-16a3-471d-bf37-c10a9653ad50
# ╠═d29ee2bd-85c9-4e33-b434-fbe4346ec1bf
# ╠═2d0d98a3-cecc-46d0-ac39-db299cf58cfd
# ╠═d4c49a72-e818-4836-8bca-d3cbc354baa5
# ╠═a8490458-7c12-489b-8d19-4ef3fd129e9f
# ╠═3c9861b5-b73d-45ba-87ab-83d2d8b60a8d
# ╠═baec7e94-71eb-4e3e-ad39-6404316f6522
# ╠═6e4fb390-ca56-488e-8066-e22d058f89a7
# ╠═3bfe5ba7-9044-4bbe-a2b8-2fbf7b9e8bd0
# ╠═b2f98194-2255-4d90-9869-0aed381e5cac
# ╠═4d3d9d03-23bc-4b67-8fc6-b84ea1618360
# ╠═3be0872f-5e8d-4bbd-9161-1c3ec6014c93
# ╠═8d9b8e11-91eb-47a4-beb9-d837e0f26161
# ╠═370ebfc9-d70d-4f5f-b55c-a9a1cebc7ae3
# ╠═45720d9f-3c25-45fb-b574-a4d4e0fe968a
# ╠═7dbb3666-fdfd-4f2e-96be-bce9922cf7b2
# ╠═3ffd2e3c-e6ef-4674-b4f2-63b4336c59dd
# ╠═9c4fcf78-bec5-482d-9dfa-3891e22dfc2d
# ╠═6c3699cc-2c8b-4bf9-a637-0843a3058800
# ╠═3b69dd32-e549-4fa0-bae9-c8fe4ffd5716
# ╠═35b75da1-97c8-4df7-bdcd-0ee0e1a428ca
# ╠═56f2b018-de9f-4ad2-abfd-920dc2169983
# ╠═31127a8d-365b-48c4-bfd8-2f968a72be09
# ╠═6c1f29a8-936e-40d3-86fa-ad8424749aa4
# ╠═8b7b9fab-2725-4e0c-af68-0687ec5921c1
# ╠═1eb4700c-20de-47a5-951d-3900d99461be
