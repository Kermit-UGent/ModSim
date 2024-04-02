### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 47f98a52-f2b2-11ea-0ac6-9b4b74b7b15d
using Plots, PlutoUI, LinearAlgebra

# ╔═╡ 5ac4989c-f2b3-11ea-3e8a-a79a4117e48c
using Symbolics

# ╔═╡ 6fe82e5a-f2b8-11ea-101d-d1ee5f9bd5b1
using ForwardDiff

# ╔═╡ 963c4910-f2b8-11ea-12bb-bb42539b1c93
using Zygote

# ╔═╡ 08cf2a80-f2b2-11ea-132d-2f26c934bbb2
md"""
# Numeric and automatic differentiation

STMO

![](https://github.com/MichielStock/STMO/blob/master/chapters/03.AutoDiff/Figures/logo.png?raw=true)

**Michiel Stock**

"""

# ╔═╡ b72a9871-3773-4a06-9a2f-09ce5a7d8496
md"## Motivation

Up to now, we confidently assumed that we would always be able to compute the derivative or gradient of any function. Despite differentiation being a relatively easy operation, it is frequenty not feasible (or desirable) to compute this by hand. *Numerical differentiation* can provide approximations of th derivate or gradient at a particular point. *Automatic differentiation* directly manipulates the computational graph to generate a function that computes the (exact) derivate. Such methods have advanced greatly in the last years and it is no exageration that their easy use in popular software libraries such as TenserFlow and PyTorch are a cornerstone of deep learning and other machine learning and scientific computing fields.
"

# ╔═╡ 5154a672-f2b2-11ea-1280-3ba5dadff7f9
md"""
## Definition of a derivative

$$\frac{\text{d}f(x)}{\text{d}x} = f'(x) = \lim _{h\to 0}{\frac {f(x+h)-f(x)}{h}}.$$

Derivation is in essence a mechanical process, following the rules below.

![](https://raw.githubusercontent.com/MichielStock/STMO/master/chapters/03.AutoDiff/Figures/derivatives.jpeg)

![](https://imgs.xkcd.com/comics/differentiation_and_integration.png)

When we work with function of several variables, we use *partial derivatives* (e.g. $\frac{\partial f(x, y)}{\partial x}$), indicating we keep all variables but $x$ fixed.
"""

# ╔═╡ f53fbf24-f2b2-11ea-2f59-d943884fee9d
md"""
Our running example:

$$f(x) = \log x + \frac{\sin x}{x}$$
"""

# ╔═╡ e6a55456-fe6f-11ea-328c-47fedb388194
md"**Assignment**: implement this function. (or another one)"

# ╔═╡ f7cd052e-f2bb-11ea-239a-47966e96e909
md"We will evaluate this function in `a`=$2."

# ╔═╡ 7244eba2-f2b3-11ea-045f-ffd208852eda
a = 2.0

# ╔═╡ 60058dbe-1514-4139-bb97-66b40fe47306
md"For the record:"

# ╔═╡ 401df7a4-f2b3-11ea-12a1-e9b46bcfcf1b
md"""
## Symbolic differentiation

Computing derivatives is relatively straightforward, as you have seen in basic calculus courses. Some software called **computer algebra systems** can manipulate symbolic expressions to generate derivatives.


For example, libraries that support symbolic differentiation:
- Maple
- Sympy (Python)
- Mathematica
- Maxima
- Symbolics (Julia)

Differentiation is *easy* compared to *integration* or *sampling*.

Advantages:
- exact derivatives!
- gives the formula for different evaluations.
    - insight in the system (i.e. variable are independent can be learned from second-order partial derivatives)
    - in some cases, closed-form solution extrema by solving $\frac{\text{d}f(x)}{\text{d}x}=0$
- no hyperparameters or tweaking: just works!

Disadvantages:
- some software not flexible enough (gradients, arrays, for-loops,...)
- sometimes explosion of terms: *expression swell*
- not always numerically optimal!
"""

# ╔═╡ 5db7350a-f2b3-11ea-293d-ffc572e0dbe8
@variables x  # define variable

# ╔═╡ 76ba97fe-b7c4-467b-85ae-9067ebcc8b14
Dx = Differential(x)  # differential operator

# ╔═╡ 92b65737-9ad7-422e-955c-e9adacd2032b
md"We can build this expression in a function:"

# ╔═╡ 67e7ea4d-0b76-4247-b22d-e20c8280f046
md"You can play with this dynamically:"

# ╔═╡ 45190ad2-d358-455a-89a1-4f9c3b3880fa
@bind n Slider(-10:10, show_value=true, default=2)

# ╔═╡ 732ad756-1707-48e7-9bfb-10707e64b118
n

# ╔═╡ c67f8b3f-0561-4f47-ade8-3c79d7920ae1
md"Derivative:"

# ╔═╡ 913234a2-f2b3-11ea-3501-e3bfc54f91fe
md"""
# Numerical differentiation

Finite difference approximation of the derivative/gradient based on a number of function evaluations.

Often based on the limit definition of a derivative. Theoretical analysis using Taylor approximation:

$$f(x + h) = f(x) + \frac{h}{1!}f'(x) + \frac{h^2}{2!}f''(x) + \frac{h^3}{3!}f^{(3)}(x)+\ldots$$

**Forward difference**

$$f'(x)\approx \frac{f(x+h) - f(x)}{h}$$

**Central difference**

$$f'(x)\approx \frac{f(x+h) - f(x-h)}{2h}$$

**Complex step method**

$$f'(x)\approx \frac{\text{Im}(f(x +ih))}{h}$$
"""

# ╔═╡ 0238579a-fe70-11ea-0968-0f11ae9557db
md"**Assignment**: Implement these functions. Hint `im` is the imaginary unit: i.e."

# ╔═╡ eef066e0-f2b3-11ea-3cfa-5f35ebf13f76
diff_fordiff(f, x; h=1e-10) = missing

# ╔═╡ f2158710-f2b3-11ea-16c9-dded5e1747d4
diff_centrdiff(f, x; h=1e-10) = missing

# ╔═╡ f54ee2b2-f2b3-11ea-24e8-f556ab9e3b85
diff_complstep(f, x; h=1e-10) = missing

# ╔═╡ 4f42234e-f2b4-11ea-2c86-93449b8dc60b
md"""
## Intermezzo: floats

Real numbers are always represented as floating point numbers in a computer.

![Encoding of a real number using a `Float32`.](https://raw.githubusercontent.com/MichielStock/STMO/master/chapters/03.AutoDiff/Figures/floats.png)

By default, Julia uses double precision floats (`Float64`). For brevity, let us take a look at the bit representation of a float. We use `Float32` for brevity's sake.
"""

# ╔═╡ 6adbc7fe-f2b4-11ea-1861-b398d517edfa
num = Float32(10.789)

# ╔═╡ 6cd1accc-f2b4-11ea-0077-e7685b77b1be
bitstring(num)

# ╔═╡ 77d4fd72-f2b4-11ea-0bd6-97fc2299ca2f
md"The first bit encodes the *sign*, here positive."

# ╔═╡ 773af248-f2b4-11ea-154a-cb67beebed3d
sign(num)

# ╔═╡ 8a5ee764-f2b4-11ea-23f3-27664d910c16
md"The next eight bits specify the *exponent*, the magnitude of the number."

# ╔═╡ 8ded9128-f2b4-11ea-0f44-9321077205de
exponent(num)

# ╔═╡ c3a1abd8-f2b4-11ea-0509-a5296ae8b361
md"While the final 23 bits specify the *mantissa*, a number between $[1,2]$ representing the precision."

# ╔═╡ c9f92e66-f2b4-11ea-322f-7778f149453d
significand(num)

# ╔═╡ cfb07cc6-f2b4-11ea-08a7-dd0b13d15b5e
md"These can be used to reconstrunct the number."

# ╔═╡ d9f2d114-f2b4-11ea-32b0-43fa1a872a97
md"The *machine precision* of a number can be retained using `eps`. This is the relative error."

# ╔═╡ de1fd700-f2b4-11ea-36dd-4bfba5ecc134
eps(1.2)

# ╔═╡ df14ce72-f2b4-11ea-307b-8fac2889681b
eps(1.2e10)

# ╔═╡ e584f976-f2b4-11ea-0245-67f6bcc4be0e
eps(1.2e-10)

# ╔═╡ eb68674c-f2b4-11ea-0cc2-f792d17a5005
md"""
This brings us with numerical issues we might encounter using numerical differentiation.

**First sin of numerical analysis**:

> *thou shalt not add small numbers to big numbers*

**second sin of numerical analysis**:

> *thou shalt not subtract numbers which are approximately equal*
"""

# ╔═╡ f4edad2c-f2b4-11ea-05fb-3d43c18e9d5c
md"### Back to numerical differentiation

Below is the absolute error using the three proposed methods."

# ╔═╡ 4e6b2b04-f2b5-11ea-0c52-5f32438e8c68
md"""
Advantages of numerical differentiation:
- easy to implement
- general, no assumptions needed

Disadvantages:
- not numerically stable (round-off errors)
- not efficient for gradients ($\mathcal{O}(n)$ evaluations for $n$-dimensional vectors)
"""

# ╔═╡ 55ff7230-f2b5-11ea-2898-db5e21c764fd
md"""## Approximations of multiplications with gradients

**Gradient-vector approximation** (the directional derivative)

$$\nabla f(\mathbf{x})^\intercal \mathbf{d} \approx \frac{f(\mathbf{x}+h\cdot\mathbf{d}) - f(\mathbf{x}-h\cdot\mathbf{d})}{2h}$$

**Hessian-vector approximation**

$$\nabla^2 f(\mathbf{x}) \mathbf{d} \approx \frac{\nabla f(\mathbf{x}+h\cdot\mathbf{d}) - \nabla f(\mathbf{x}-h\cdot\mathbf{d})}{2h}$$
"""

# ╔═╡ e5c12243-70ad-4e52-9bab-aae432203e21
md"Generate an example"

# ╔═╡ 67e9fa88-ff03-11ea-1177-c34331398d6c
md"$$g(\mathbf{x}) = \exp(-\mathbf{x}^\intercal A\mathbf{x})$$"

# ╔═╡ 8a738272-f2b5-11ea-2de6-f9ae6b30a466
md"Correct gradient and Hessian (by hand)"

# ╔═╡ c2924514-fe7f-11ea-2b01-81e9e0e1c491
h = 1e-10

# ╔═╡ d6c6d4ea-f2b7-11ea-0b41-d7f629cdcd2f
md"""
## Forward differentiation

Accumulation of the gradients along the *computational graph*.

![](https://raw.githubusercontent.com/MichielStock/STMO/master/chapters/03.AutoDiff/Figures/forwarddiff.png)

Forward differentiation computes the gradient from the inputs to the outputs.

### Differentiation rules

**Sum rule**:

$$\frac{\partial (f(x)+g(x))}{\partial x} =  \frac{\partial f(x)}{\partial x} + \frac{\partial f(x)}{\partial x}$$

**Product rule**:

$$\frac{\partial (f(x)g(x))}{\partial x} =  f(x)\frac{\partial g(x)}{\partial x} + g(x)\frac{\partial f(x)}{\partial x}$$

**Chain rule**:

$$\frac{\partial (g(f(x))}{\partial x} =  \frac{\partial g(u)}{\partial u}\mid_{u=f(x)} \frac{\partial f(x)}{\partial x}$$

## Dual numbers

Forward differentiation can be viewed as evaluating function using *dual numbers*, which can be viewed as truncated Taylor series:

$$v + \dot{v}\epsilon\,,$$

where $v,\dot{v}\in\mathbb{R}$ and $\epsilon$ a nilpotent number, i.e. $\epsilon^2=0$. For example, we have

$$(v + \dot{v}\epsilon) + (u + \dot{u}\epsilon) = (v+u) + (\dot{v} +\dot{u})\epsilon$$


$$(v + \dot{v}\epsilon)(u + \dot{u}\epsilon) = (vu) + (v\dot{u} +\dot{v}u)\epsilon\,.$$


These dual numbers can be used as

$$f(v+\dot{v}\epsilon) = f(v) + f'(v)\dot{v}\epsilon\,.$$
"""

# ╔═╡ fe18341e-f2b7-11ea-0236-03b44eebd328
struct Dual{T}
    v::T
    vdot::T
end

# ╔═╡ 3f2a4cc2-26fa-4c3e-bfb7-f3ef019b6d78
Base.show(io::IO, a::Dual) = print(io, "$(a.v) + $(a.vdot)ϵ")  # nice printing

# ╔═╡ 9e22156e-0880-11eb-24a8-3d0a41262c36
ϵ = Dual(0.0, 1.0)

# ╔═╡ ab39b554-0880-11eb-05de-71949f71b440
Dual(2.0, 3.0)

# ╔═╡ 0aba93a6-f2b8-11ea-1461-d96878b7a5c1
md"Let's implement some basic rules showing linearity."

# ╔═╡ 3f5a58da-f2b8-11ea-3ab0-5db5e9199bb1
md"And some more advanced ones, based on differentiation."

# ╔═╡ eaacb5a5-c748-4223-98e4-be73d3f2086f


# ╔═╡ 263e18dc-fe70-11ea-03a0-8fd5e031b06f
md"**Assignment**: complete this code."

# ╔═╡ 2f18c362-f2b8-11ea-1beb-719ff5cd5f8d
Base.:sin(a::Dual) = missing

# ╔═╡ 393f92a6-fe70-11ea-024e-ff30dcfc7103
Base.:cos(a::Dual) = missing

# ╔═╡ 329cb67e-f2b8-11ea-07b1-e13fd8ddbe42
Base.:exp(a::Dual) = missing

# ╔═╡ 36a2ca10-f2b8-11ea-0fae-53b552c76115
Base.:log(a::Dual) = missing

# ╔═╡ 555c8d88-f2b8-11ea-34d8-a17ce87ccdc1
md"This directly works for vectors!"

# ╔═╡ 6971e7e6-f2b8-11ea-2a76-9dda6b1455b6
md"In practice, we prefer to use a package to do this."

# ╔═╡ 8a254852-f2b8-11ea-142d-71bb270a5c01
md"""
Forward differentiation:

- exact gradients!
- computational complexity scales with **number of inputs**
- used when you have more outputs than inputs
"""

# ╔═╡ 91d0e03e-f2b8-11ea-07e5-db056636ac21
md"""
# Reverse differentiation

Compute the gradient from the output toward the inputs using the chain rule.

![](https://raw.githubusercontent.com/MichielStock/STMO/master/chapters/03.AutoDiff/Figures/reversediff.png)

Reverse differentiation:

- also exact!
- main workhorse for training artificial neural networks.
- efficient when more inputs than outputs (machine learning: thousands of parameters vs. one loss)
"""

# ╔═╡ a041fbf8-f2b8-11ea-28ba-0772497c649f
md"Constructing the derivative or gradient can be done by appending `'` to the function."

# ╔═╡ bdce5964-f2b8-11ea-03d1-57705d295040
md"More verbose, but exactly the same:"

# ╔═╡ e1057afe-fe6e-11ea-27fe-6def4d8e1f2d
md"This works for function with multiple inputs."

# ╔═╡ d6d310ee-f2b8-11ea-24e1-3b405f071545
md"Fuctions with vectorial input."

# ╔═╡ dbb1e68c-f2b8-11ea-1220-6f5715822e3d
md"Finding the Hessian:"

# ╔═╡ ee52d2b8-f2b8-11ea-15c4-19e9f2ad10e2
md"""
## Artificial neural networks

Multi-layer perceptron.

![](https://raw.githubusercontent.com/MichielStock/STMO/master/chapters/03.AutoDiff/Figures/ANN_example.png)

Forward differentiation.

![](https://raw.githubusercontent.com/MichielStock/STMO/master/chapters/03.AutoDiff/Figures/Forwardprop.png)

Reverse differentation or backpropagation.

![](https://raw.githubusercontent.com/MichielStock/STMO/master/chapters/03.AutoDiff/Figures/Backprop.png)

Returns effect of changing layer output on the loss. Can be related directly to the parameters!
"""

# ╔═╡ 213fc38a-f380-11ea-135e-6f72924c19aa
md"""
## Differentiating complex objects

Automatic differentiation can be used beyond machine learning and optimization:

- [physical engines](https://arxiv.org/abs/1611.01652) to learn robot control
- differentiating [protein](https://github.com/lupoglaz/TorchProteinLibrary) [structures](https://www.cell.com/cell-systems/fulltext/S2405-4712(19)30076-6)
- Sinkhorn algorithm
- [dynamic programming](https://arxiv.org/abs/1802.03676)
- [differential equations](https://julialang.org/blog/2019/01/fluxdiffeq)
- [molecules](https://www.nature.com/articles/s41592-021-01283-4)

Everything is computed by some straightforward and differentiable functions!
"""

# ╔═╡ 059679c0-d0ff-453b-b897-b889bbea326d
md"""
### Example: differentiable sequence alignment

One can differentiate the Needleman-Wunsch algorithm for sequence alignment. This yields insight in how the alignment score would change if one changes the alignment parameters.
"""

# ╔═╡ a20f0158-f32b-45e8-8bd5-6f693eca1860
abstract type Regularizer end

# ╔═╡ 23b0cc81-b551-446b-991f-54bd19ac63cd
struct NegEntropy{T<:Number} <: Regularizer
    γ::T
end

# ╔═╡ a9961630-08be-48f8-b1cc-4a9d67dfc519
function max_argmax(Ω::NegEntropy, x)
    γ = Ω.γ
    m = maximum(x)
    # stable log-sum-exp computation
    lse = m + γ * log(sum(exp, (x .- m) ./ γ))
    # stable softmax
    q = exp.((x .- lse) ./ γ)
    return lse, q  # smooth max and grad
end

# ╔═╡ e4d2350e-147b-4c3c-a9ff-6de66a54c48a
function ∇needleman_wunsch(Ω::Regularizer, θ, (cˢ, cᵗ))
    n, m = size(θ)
    D = zeros(n+1, m+1)        # initialize DP matrix
    D[2:n+1,1] .= -cumsum(cˢ)  # cost of starting with gaps in s
    D[1,2:m+1] .= -cumsum(cᵗ)  # cost of starting with gaps in t
    E = zeros(n+2, m+2)        # matrix for the gradient
    E[n+2,m+2] = 1.0
    Q = zeros(n+2, m+2, 3)     # matrix for backtracking
    Q[n+2,m+2,2] = 1.0  
    # forward pass, performing dynamic programming
    for i in 1:n, j in 1:m
        v, q = max_argmax(Ω, (D[i+1,j] - cˢ[i],   # gap in s
                              D[i,j] + θ[i,j],    # match
                              D[i,j+1] - cᵗ[j]))  # gap in t
        D[i+1,j+1] = v      # store smooth max
        Q[i+1,j+1,:] .= q   # store directions
    end
    v = D[n+1,m+1]        # get alignment score
    # backtracking through the directions to compute the gradient
    for i in n:-1:1, j in m:-1:1
        E[i+1,j+1] = Q[i+1,j+2,1] * E[i+1,j+2] +  
                     Q[i+2,j+2,2] * E[i+2,j+2] +
                     Q[i+2,j+1,3] * E[i+2,j+1]
    end
    return v, E[2:n+1,2:m+1]  # value and gradient
end

# ╔═╡ 501a0db3-538a-43f3-a769-158c824bca31
md"For example, consider two sequences and the Hamming distance."

# ╔═╡ ca64df96-082a-40dd-b82c-64016257fcab
s, t = "aattcaa", "atctaca"

# ╔═╡ ba43bc8b-27d6-490d-971e-7e440b991cfb
θ = [sᵢ==tᵢ for sᵢ in s, tᵢ in t]

# ╔═╡ 2d7dbec3-82eb-41f8-a959-b458e136d18c
heatmap(θ, flipy=true, yticks=(1:length(s), s), xticks=(1:length(t), t), title="theta")

# ╔═╡ 378500c1-9f19-4111-82cf-5a57c30d66ef
cˢ, cᵗ = 0.1ones(length(s)), 0.1ones(length(t))

# ╔═╡ 8586855b-2ad7-4fe4-af1a-2a6f1dd052b3
@bind logγ Slider(-8:0.2:4, default=-4)

# ╔═╡ 8877c563-d10a-41e4-be8b-546d27623b67
γ = exp(logγ)

# ╔═╡ a73e4d59-34f4-412f-afa4-f3b66244b35e
γ

# ╔═╡ 6707fd27-78f6-4663-b95d-a016d2c44c3b
# alignment value and score
v, E = ∇needleman_wunsch(NegEntropy(γ), θ, (cˢ, cᵗ))

# ╔═╡ 014ec92c-f2b8-11ea-3986-8d74a8c2458d
begin
	Base.:+(a::Dual, b::Dual) = Dual(a.v + b.v, a.vdot + b.vdot)
	Base.:-(a::Dual, b::Dual) = Dual(a.v - b.v, a.vdot - b.vdot)
	Base.:-(a::Dual, b::Real) = Dual(a.v - b, a.vdot)
	Base.:-(a::Dual) = Dual(-a.v, -a.vdot)
	Base.:*(a::Dual, b::Dual) = Dual(a.v * b.v, a.v * b.vdot + b.v * a.vdot)
	Base.:+(c::Real, b::Dual) = Dual(c + b.v, b.vdot)
	Base.:*(v::Real, b::Dual) = Dual(v, 0.0) * b
	Base.:*(a::Dual, b::Real) = v * a	
end

# ╔═╡ 9ecd1cd6-8412-4fe1-9d79-676f9830ec2c
imag_number = 2 + 3im  # im is a special type of number

# ╔═╡ af02692c-c1be-4f7c-a09f-fbd88b0c4bd8
real(imag_number), imag(imag_number) # extract real and imaginary parts

# ╔═╡ d56b74f2-f2b4-11ea-12b4-5f54d4999b38
sign(num) * significand(num) * 2^exponent(num)

# ╔═╡ 6fd08e42-f2b5-11ea-36e6-d7ff10412e2b
xvect = 2rand(10)

# ╔═╡ 4a6270a3-a271-40da-84d0-28247f287c22
2.0 + 3.0ϵ  # now this works!

# ╔═╡ 48e5da96-f2b8-11ea-0108-2104f2c8ac24
# you can derive this directly from the definition!
Base.:/(a::Dual, b::Dual) = Dual(a.v / b.v, (a.vdot * b.v - a.v * b.vdot) / b.v^2)

# ╔═╡ 2be8656c-f2b3-11ea-2244-5740c807deff
f(x) = log(x) + sin(x) / x

# ╔═╡ 21882f0e-8d54-41ae-84c6-7637ccc52c46
f(a)

# ╔═╡ d5edd644-fe7c-11ea-2aa8-c772ea2f1bde
f(x)

# ╔═╡ 64b28e4a-f2b3-11ea-0fc0-7d37b4794686
Dx(f(x))

# ╔═╡ 8a4a8499-af01-46ea-8196-fb2d718d0d17
df_sym = expand_derivatives(Dx(f(x)))  # this expands the derviatve operator

# ╔═╡ a38e2cce-851c-4065-8947-1926fa9ca5b4
df = build_function(df_sym, x) |> eval  #builds an expression and turns it into a function

# ╔═╡ 684ff7f4-f2b3-11ea-14f7-7d3f329ecd79
df(a)

# ╔═╡ 2ca60bca-f2b4-11ea-0fec-c50f4e306536
true_diff = df(a);  # save for checking

# ╔═╡ 7a6a18de-f2b3-11ea-1b3d-ef7a2b1b0434
if !ismissing(f(a))
plot(f, 1, 10, label="\$f(x)\$", xlabel="\$x\$", lw=2, color=:green)
plot!(df, 1, 10, label="\$f'(x)\$", lw=2, color=:orange)
end

# ╔═╡ 08a565e0-f2b4-11ea-2ab6-f5f3a1b3c22d
diff_fordiff(f, a)

# ╔═╡ 1055a996-f2b4-11ea-14d2-6116cbbaa6ba
diff_centrdiff(f, a)

# ╔═╡ 13db037a-f2b4-11ea-2a38-bd7c1c60f7d9
diff_complstep(f, a)

# ╔═╡ 4ee45daa-f2b8-11ea-10ea-f929853879b1
f(Dual(a, 1.0))  # works when rules are provided

# ╔═╡ 723fe4ca-f2b8-11ea-1440-a9e1cd213b66
ForwardDiff.derivative(f, a)

# ╔═╡ 9b5db398-f2b8-11ea-00bc-77de99792925
f'(a)  # that's it

# ╔═╡ c6728784-f2b8-11ea-12b5-d7eb039d77ce
Zygote.gradient(f, a)  # returns a tuple, since you can differentiate w.r.t. multiple arguments

# ╔═╡ c8a5c6f8-4ca1-4c7e-9cc1-d2e5c9e9da59
symexpr = (x + 2)^n / (x+2)

# ╔═╡ 997ec642-a865-4c97-9aa1-d4d9df9e1dd0
expand_derivatives(Dx(symexpr)) |> simplify

# ╔═╡ 688a1f40-f2b5-11ea-3fa8-ebb735323ddd
grad_vect(f, x, d; h=1e-10) = (f(x .+ h * d) - f(x .- h * d)) / (2h)

# ╔═╡ 6d6da8ba-f2b5-11ea-1064-df27ecea7d13
dvect = randn(10) / 10

# ╔═╡ 9a87095e-f2b5-11ea-334e-57a637e20c43
A = randn(10, 10) |> A -> (A * A' + I) / 100  # make sym and PD

# ╔═╡ 75ae48fe-f2b5-11ea-0614-912091cbfb6a
g(x) = exp(- sum(x .* (A * x)))

# ╔═╡ 8fe85430-f2b5-11ea-3fe7-5d9a0db03932
g(xvect)

# ╔═╡ b455e8ac-fe7f-11ea-1a8a-db85d11ccba7
grad_vect(g, xvect, dvect)

# ╔═╡ 76284520-f2b8-11ea-077e-43167cc7484c
ForwardDiff.gradient(g, xvect)

# ╔═╡ d3ae2c14-f2b8-11ea-1ac2-fb5e0b4b8500
g'(xvect)

# ╔═╡ df557ed2-f2b8-11ea-22ee-e1c90a2655bf
Zygote.hessian(g, xvect)

# ╔═╡ 7b95d73c-f2b5-11ea-1b29-858f7b6c99bf
∇g(x) = -2g(x) * A * x

# ╔═╡ b372c25c-f2b5-11ea-1c55-e5045e0dc5b7
∇g(xvect)

# ╔═╡ b6e48b12-f2b5-11ea-3bf2-81e5b2a3dc98
∇g(xvect)' * dvect

# ╔═╡ 836ca94a-f2b5-11ea-07e0-bd830f810fba
∇²g(x) = -2g(x) * A - 2A * x * ∇g(x)'

# ╔═╡ cd2821c4-fe7f-11ea-0c41-8feab49be07a
∇²g(xvect) * dvect

# ╔═╡ bebbd674-f2b5-11ea-19c6-4fe1257b0277
(∇g(xvect + h * dvect) - ∇g(xvect - h * dvect)) / 2h

# ╔═╡ 272bdcca-aeb1-405f-8251-33f13cddbdc5
(2.0+3.0ϵ) / (5.0+8.0ϵ)

# ╔═╡ 59eee760-f2b8-11ea-29f6-732292e2f95d
q(x) = 10.0 * x[1] * x[2] + x[1] * x[1] + sin(x[1]) / x[2]

# ╔═╡ 5ccd625e-f2b8-11ea-38c1-098e2f31ce53
q([1, 2])

# ╔═╡ 60b52f6e-f2b8-11ea-2afc-4be2d4ae3e19
q(Dual.([1, 2], [1, 0]))  # partial wrt x1

# ╔═╡ 63d5e0be-f2b8-11ea-16e8-d9528b4e8945
q(Dual.([1, 2], [0, 1]))  # partial wrt x2

# ╔═╡ 85a5c518-f2b8-11ea-2e12-25d16a12e7c8
ForwardDiff.gradient(q, [1, 2])

# ╔═╡ eb76e252-fe6e-11ea-0939-c51ff3c73827
Zygote.gradient((x1, x2) -> (x1^2 + x2^2 - 0.1x1*x2) / (x1 + 1.0),
					0.2, 0.3) 

# ╔═╡ 18e151c4-0b0e-485b-81d2-f18a2dd4592d
heatmap(E, flipy=true, yticks=(1:length(s), s), xticks=(1:length(t), t), title="gradient of NW", color=:speed)

# ╔═╡ 79221eee-f2be-11ea-3e43-4596f5ce5694
md"""
## Application: inverse kinematics of a robot arm

As an illustration of using gradients, let us study a [simple robot arm](https://appliedgo.net/roboticarm/). This arm consists of three joints that can be moved in the 2D plane. These lengths of the three segments are, respectivly 2, 1 and 1.5 meter.

The forward kinetics (i.e., going from the angles of the joints to the position of the arms) can be computed by the following functions.
"""

# ╔═╡ fb91a9c4-f2c0-11ea-38c4-5b554446e997
pos_second_joint((θ₁, θ₂, θ₃)) = [2cos(θ₁), 2sin(θ₁)]

# ╔═╡ 07c45d78-f2c3-11ea-131e-b9f0b67ca3cd
pos_third_joint((θ₁, θ₂, θ₃)) = pos_second_joint((θ₁, θ₂, θ₃)) + [cos( θ₂), sin( θ₂)]

# ╔═╡ 33023d10-f2c1-11ea-3375-e9ac23172ab6
pos_hand((θ₁, θ₂, θ₃)) = pos_third_joint((θ₁, θ₂, θ₃)) .+ [1.5cos(θ₃), 1.5sin(θ₃)]

# ╔═╡ 31cdefd2-f381-11ea-3037-e9c023ebc1e9
md"Suppose that there is a teapot at a certain location we want to grasp with the arm."

# ╔═╡ 4def0626-f2c1-11ea-2f01-31a8c42b121c
position_teapot = [-1, 2.3]

# ╔═╡ 46c54108-f381-11ea-141f-31a696a7a8ef
md"Below is a plot of the arm and the teapot. Can you turn the angles to bring the hand to the position of the teapot?"

# ╔═╡ 8049560e-f383-11ea-1c8e-39c629d20982
@bind θ₁ Slider(0:0.01:2π, show_value=true)

# ╔═╡ c7b5a1aa-f383-11ea-2a84-8da9c9e60dbc
@bind θ₂ Slider(0:0.01:2π, show_value=true)

# ╔═╡ cfad1a8c-f383-11ea-3dad-214f3da9485d
@bind θ₃ Slider(0:0.01:2π, show_value=true)

# ╔═╡ 83ddbe6a-f381-11ea-2415-255cbfb2590b
md"First complete a function to determine the distance between the arm and the teapot."

# ╔═╡ 30e9891c-f2c5-11ea-022d-5dbf9c064e43
distance_to_teapot(θ) = √(sum(abs2, position_teapot .- pos_hand(θ)))

# ╔═╡ 5812f84c-f2c1-11ea-0f6c-31b9fedc3704
begin
	psj = pos_second_joint((θ₁, θ₂, θ₃))
	ptj = pos_third_joint((θ₁, θ₂, θ₃))
	phand = pos_hand((θ₁, θ₂, θ₃))
	plot([0, psj[1], ptj[1], phand[1]], [0, psj[2],ptj[2], phand[2]], label="Arm")
	scatter!([0], [0], label="First joint (fixed)")
	scatter!([psj[1]], [psj[2]], label="Second joint")
	scatter!([ptj[1]], [ptj[2]], label="Second joint")
	scatter!([phand[1]], [phand[2]], label="Hand")
	scatter!([position_teapot[1]], [position_teapot[2]], m=:star, label="teapot")
	title!("Distance to teapot = $(distance_to_teapot((θ₁, θ₂, θ₃)))")
	xlims!(-5, 5)
	ylims!(-5, 5)
end

# ╔═╡ b07cb958-f381-11ea-3378-5517440464ca
md"Can you use gradients to position the arm?"

# ╔═╡ eaac2514-f8fe-11ea-2229-81f5b59f58fe
distance_to_teapot([θ₁, θ₂, θ₃])

# ╔═╡ fb044216-f381-11ea-2faf-fd53dad74950
distance_to_teapot'([θ₁, θ₂, θ₃])

# ╔═╡ 675b222c-f382-11ea-257a-afcb5a5afdd8
begin
	# initial value of θ
	θs = [0.0, 0.0, 0.0]
	for i in 1:500  # 500 steps, decaying
		θs .-= 1/i * distance_to_teapot'(θs)
	end
end

# ╔═╡ 818c1f9e-ff1b-11ea-097f-cba2f30af994
θs  # optimal values

# ╔═╡ 7aff4178-bcac-4cfa-ab22-de3506cf73a8
distance_to_teapot(θs)

# ╔═╡ 881dedd4-f385-11ea-342e-a3700e90ac66
md"""
## Exercise

Consider the *Wheeler's Ridge* function (also called "the banana function"):

$$f(\mathbf{x}) = -\exp(-(x_1 x_2 - a)^2 -(x_2 -a)^2)\,,$$

at the point $\mathbf{x}_0=[1.5, 1]^T$. We set $a=1.5$.

Implement this function.
"""


# ╔═╡ c1e19e26-f385-11ea-17b6-6d2d2e5cb7c0
fwr(x; a=1.5) = missing

# ╔═╡ fe18a0a0-f8fe-11ea-1595-4d95256dabca
fwr([1,2]) isa Missing || contour(-0:0.1:3, -0:0.1:3, (x1, x2) -> fwr([x1, x2]), color=:speed)

# ╔═╡ df5f20d8-f385-11ea-104d-9b9e7b094eb7
md"""

**Assignments**

1. Compute the gradient by hand.
2. Compute the gradient using symbolic differentiation.
2. Find the gradient and Hessian at $\mathbf{x}_0$ by numerical differentiation.
3. Compute the gradient and Hessian at $\mathbf{x}_0$ using automatic differentiation.
4. (optional) Use the function `quiver!` to draw the gradient as a vector field on the contour plot.

"""

# ╔═╡ e4a088f0-f385-11ea-00da-b12e40e4357f
x₀ = [1.5, 1.0]  # btw written as x\_0<TAB>

# ╔═╡ ec87028b-d9d4-4dc7-9dd5-36eaa1dd49d4
fwr(x₀)

# ╔═╡ ddd214d8-608a-4123-973f-b331c6ab331a
md"Symbolic"

# ╔═╡ 8f3fecac-f8fd-11ea-136d-1513f60ee3ed


# ╔═╡ f8863a18-f385-11ea-0965-0f074e58d9b2


# ╔═╡ 4c3817f6-396e-4485-a7dd-0f507cec124a
md"Numberic"

# ╔═╡ 78b72cb1-52f3-4dc3-a721-b8a89a24312e


# ╔═╡ 7180e310-a6f7-4feb-aff3-28e1c772b441


# ╔═╡ 02b42a67-7ebe-4ca2-ac55-aa481b30e997
md"Forward (dual numbers)"

# ╔═╡ 043e3e2f-3876-41c0-bf70-5ef53b5db986


# ╔═╡ 8ba368c9-3b04-469e-9556-1573cbc52185


# ╔═╡ 4da4a1b6-b7c8-432d-86ed-c45a499a1880


# ╔═╡ c573ba19-4be5-4afd-b9f5-8b75e9b18ad6
md"Reverse"

# ╔═╡ 639528b1-ecff-4a05-ae56-9a08b779edc1


# ╔═╡ 838292fc-2c0e-4bb9-9df9-6cd197c4a4fd


# ╔═╡ ca20cc75-3033-4bd1-bb5c-cae61c855de7


# ╔═╡ 34d48ec0-f2b9-11ea-147e-d3554a9fcc04
md"""
# References

- Gunes et. al. (2015) *Automatic differentiation in machine learning: a survey*
- Kochenderfer, M. J. and Wheeler, T., '*Algorithms for Optimization*'. MIT Press (2019)
"""

# ╔═╡ 0a576c65-3685-4987-93f1-980ace78d6f4
module Solution

	export diff_fordiff, diff_centrdiff, diff_complstep
	
	diff_fordiff(f, x; h=1e-10) = (f(x + h) - f(x)) / h
	diff_centrdiff(f, x; h=1e-10) = (f(x + h) - f(x - h)) / 2h
	diff_complstep(f, x; h=1e-10) = imag(f(x + im * h)) / h

	
end

# ╔═╡ 91aff94e-e77e-4810-82ad-fe5ab3eafdad
md"Show solution dual numbers: $(@bind show_dual_sol CheckBox())"

# ╔═╡ 8eb9ed90-981b-4094-baec-a521614fb3a1
if show_dual_sol 
	md"solution dual numbers:
```julia
	Base.:sin(a::Dual) = Dual(sin(a.v), cos(a.v) * a.vdot)
	Base.:cos(a::Dual) = Dual(cos(a.v), -sin(a.v) * a.vdot)
	Base.:exp(a::Dual) = Dual(exp(a.v), exp(a.v) * a.vdot)
	Base.:log(a::Dual) = Dual(log(a.v), 1.0 / a.v * a.vdot)
```"
end

# ╔═╡ d682439b-bc6b-456d-8b5d-2489d095e55f
md"Show solution quiver plot: $(@bind show_quiver_sol CheckBox())"

# ╔═╡ 2d4c0bf6-3936-449d-8341-56f9c5659a6e
if show_quiver_sol 
	md"solution quiver plot:
```julia
	x1vals, x2vals = -0:0.5:4, -0:0.5:3 
	
	contour(-0.1:0.1:4, -0.1:0.1:3, (x1, x2) -> fwr([x1, x2]),
			color=:speed)
	
	quiver!([x1 for x1 in x1vals for x2 in x2vals],
			[x2 for x1 in x1vals for x2 in x2vals],
			quiver=(x1,x2)->0.25.*fwr'((x1,x2)))
	xlabel!(\"x_1\")
	ylabel!(\"x_2\")
```"
end

# ╔═╡ fe4dfa5c-ea26-4f7c-aa68-7b3ec5f45bdc
begin
	myblue = "#304da5"
	mygreen = "#2a9d8f"
	myyellow = "#e9c46a"
	myorange = "#f4a261"
	myred = "#e76f51"
	myblack = "#50514F"

	mycolors = [myblue, myred, mygreen, myorange, myyellow]
end;

# ╔═╡ c25c9a89-69bd-4dd8-b710-5e66f68346b4
plots = Dict()

# ╔═╡ 2e00ef04-f2b3-11ea-1960-8136321c71db
plots["example_diff"] = plot(f, 1, 5, label="f(x)", color=myblue, lw=2)

# ╔═╡ fd9ed266-f2b4-11ea-3c53-3bd2fe30a7f1
let
	fexamp(x) = 64x*(1-x)*(1-2x)^2*(1-8x+8x^2)^2
	#dfexamp = diff(fexamp(x), x)
	dfexamp = build_function(expand_derivatives(Dx(fexamp(x))), x) |> eval
	error(diff, h; x=1.0) = max(abs(Float64(dfexamp(x)) - diff(fexamp, x, h=h)), 1e-50)
	stepsizes = map(t->10.0^t, -20:0.1:-1);
	p = plot(stepsizes, error.(Solution.diff_fordiff, stepsizes), label="forward difference",
    xscale=:log10, yscale=:log10, lw=2, legend=:bottomright, color=myblue)
	plot!(stepsizes, error.(Solution.diff_centrdiff, stepsizes), label="central difference", 		lw=2, color=myred)
	plot!(stepsizes, error.(Solution.diff_complstep, stepsizes), label="complex step", lw=2,
            color=myyellow)
	xlabel!("\$h\$")
	ylabel!("absolute error")
	plots["numdiff_error"] = p
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[compat]
ForwardDiff = "~0.10.36"
Plots = "~1.39.0"
PlutoUI = "~0.7.58"
Symbolics = "~5.14.1"
Zygote = "~0.6.69"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.0"
manifest_format = "2.0"
project_hash = "e209e8e05c069e5b2b28053138a2dc0c41198d32"

[[deps.ADTypes]]
git-tree-sha1 = "016833eb52ba2d6bea9fcb50ca295980e728ee24"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "0.2.7"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Preferences", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "d7832de8cf7af26abac741f10372080ac6cb73df"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.34.7"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "0f748c81756f2e5e6854298f11ad8b2dfae6911a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.0"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "cde29ddf7e5726c9fb511f340244ea3481267608"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.7.2"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "c5aeb516a84459e0318a02507d2261edad97eb75"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.7.1"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bijections]]
git-tree-sha1 = "c9b163bd832e023571e86d0b90d9de92a9879088"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.6"

[[deps.BitFlags]]
git-tree-sha1 = "2dc09997850d68179b69dafb58ae806167a32b1b"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.8"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a4c43f59baa34011e303e76f5c8c91bf58415aaf"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRules]]
deps = ["Adapt", "ChainRulesCore", "Compat", "Distributed", "GPUArraysCore", "IrrationalConstants", "LinearAlgebra", "Random", "RealDot", "SparseArrays", "SparseInverseSubset", "Statistics", "StructArrays", "SuiteSparse"]
git-tree-sha1 = "4e42872be98fa3343c4f8458cbda8c5c6a6fa97c"
uuid = "082447d4-558c-5d27-93f4-14fc19e9eca2"
version = "1.63.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "575cd02e080939a33b6df6c5853d14924c08e35b"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.23.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "59939d8a997469ee05c4b4944560a820f9ba0d73"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.4"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "c955881e3c981181362ae4088b35995446298b80"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.14.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.2+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "02d2316b7ffceff992f3096ae48c7829a8aa0638"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.3"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "6cbbd4d241d7e6579ab354737f4dd95ca43946e1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.1"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "260fd2400ed2dab602a7c15cf10c1933c59930a2"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.5"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "0f4b5d62a88d8f59003e43c25a8a90de9eb76317"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.18"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "7c302d7a5fec5214eb8a5a4c466dcf7a51fcf169"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.107"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "51b4b84d33ec5e0955b55ff4b748b99ce2c3faa9"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.6.7"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "0bb0a6f812213ecc8fbbcf472f4a993036858971"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.5.5"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random"]
git-tree-sha1 = "5b93957f6dcd33fc343044af3d48c215be2562f1"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.9.3"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "ff38ba61beff76b8f4acad8ab0c97ef73bb670cb"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.9+0"

[[deps.GPUArrays]]
deps = ["Adapt", "GPUArraysCore", "LLVM", "LinearAlgebra", "Printf", "Random", "Reexport", "Serialization", "Statistics"]
git-tree-sha1 = "85d7fb51afb3def5dcb85ad31c3707795c8bccc1"
uuid = "0c68f7d7-f131-5f86-a1c3-88cf8149b2d7"
version = "9.1.0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "2d6ca471a6c7b536127afccfa7564b5b39227fe0"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.5"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "27442171f28c952804dede8ff72828a96f2bfc1f"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.10"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "025d171a2847f616becc0f84c8dc62fe18f0f6dd"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.10+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "359a1ba2e320790ddbe4ee8b4d54a305c0ea2aff"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.0+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.Groebner]]
deps = ["AbstractAlgebra", "Combinatorics", "ExprTools", "Logging", "MultivariatePolynomials", "PrecompileTools", "PrettyTables", "Primes", "Printf", "Random", "SIMD", "TimerOutputs"]
git-tree-sha1 = "6b505ef15e55bdc5bb3ddbcfebdff1c9e67081e8"
uuid = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
version = "0.5.1"

[[deps.GroupsCore]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "6df9cd6ee79fc59feab33f63a1b3c9e95e2461d5"
uuid = "d5909c97-4eac-4ecc-a3dc-fdd0858a4120"
version = "0.4.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "8e59b47b9dc525b70550ca082ce85bcd7f5477cd"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.5"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "8b72179abc660bfab5e28472e019392b97d0985c"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.4"

[[deps.IRTools]]
deps = ["InteractiveUtils", "MacroTools", "Test"]
git-tree-sha1 = "5d8c5713f38f7bc029e26627b687710ba406d0dd"
uuid = "7869d1d1-7146-5819-86e3-90919afe41df"
version = "0.4.12"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "b8ffb903da9f7b8cf695a8bead8e01814aa24b30"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IntervalSets]]
git-tree-sha1 = "dba9ddf07f77f60450fe5d2e2beb9854d9a49bd0"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.10"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "a53ebe394b71470c7f97c2e7e170d51df21b17af"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.7"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3336abae9a713d2210bb57ab484b1e065edd7d23"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.2+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVM]]
deps = ["CEnum", "LLVMExtra_jll", "Libdl", "Preferences", "Printf", "Requires", "Unicode"]
git-tree-sha1 = "ab01dde107f21aa76144d0771dccc08f152ccac7"
uuid = "929cbde3-209d-540e-8aea-75f648917ca0"
version = "6.6.2"

    [deps.LLVM.extensions]
    BFloat16sExt = "BFloat16s"

    [deps.LLVM.weakdeps]
    BFloat16s = "ab4f0b2a-ad5b-11e8-123f-65d77653426b"

[[deps.LLVMExtra_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl", "TOML"]
git-tree-sha1 = "88b916503aac4fb7f701bb625cd84ca5dd1677bc"
uuid = "dad2f222-ce93-54a1-a47d-0025e8a3acab"
version = "0.0.29+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "ForwardDiff", "LinearAlgebra", "MacroTools", "PreallocationTools", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "d1f981fba6eb3ec393eede4821bca3f2b7592cd4"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.15.1"

[[deps.LambertW]]
git-tree-sha1 = "c5ffc834de5d61d00d2b0e18c96267cffc21f648"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.6"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "cad560042a7cc108f5a4c24ea1431a9221f22c1b"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.2"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "dae976433497a2f841baadea93d27e68f1a12a97"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.39.3+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0a04a1318df1bf510beb2562cf90fb0c386f58c4"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.39.3+1"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "18144f3e9cbe9b15b070288eef858f71b291ce37"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.27"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "769c9175942d91ed9b83fa929eee4fe6a1d128ad"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.4"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "302fd161eb1c439e4115b51ae456da4e9984f130"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.4.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "af81a32750ebc831ee28bdaaba6e1067decef51e"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.2"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "60e3045590bd104a16fefb12836c00c0ef8c7f8c"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.13+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "7b1a9df27f072ac4c9c7cbe5efb198489258d1f5"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "ccee59c6e48e6f2edf8a5b64dc817b6729f99eb5"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.39.0"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "71a22244e352aa8c5f0f2adde4150f62368a3f2e"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.58"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff", "Requires"]
git-tree-sha1 = "01ac95fca7daabe77a9cb705862bd87016af9ddb"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.13"

    [deps.PreallocationTools.extensions]
    PreallocationToolsReverseDiffExt = "ReverseDiff"

    [deps.PreallocationTools.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "88b895d13d53b5577fd53379d913b9ab9ac82660"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.3.1"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "cb420f77dc474d23ee47ca8d14c90810cafe69e7"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.6"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "37b7bb7aabf9a085e0044307e1717436117f2b3b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.5.3+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9b23c31e76e333e6fb4c1595ae6afa74966a729e"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.9.4"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RandomExtensions]]
deps = ["Random", "SparseArrays"]
git-tree-sha1 = "b8a399e95663485820000f26b6a43c794e166a49"
uuid = "fb686558-2515-59ef-acaa-46db3789a887"
version = "0.4.4"

[[deps.RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "Requires", "SparseArrays", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "27ee1c03e732c488ecce1a25f0d7da9b5d936574"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.3.3"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "6aacc5eefe8415f47b3e34214c1d79d2674a0ba2"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.12"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "d8911cc125da009051fb35322415641d02d9e37f"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.4.6"

[[deps.SciMLBase]]
deps = ["ADTypes", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FillArrays", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables", "TruncatedStacktraces"]
git-tree-sha1 = "09324a0ae70c52a45b91b236c62065f78b099c37"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.15.2"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseZygoteExt = "Zygote"

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLOperators]]
deps = ["ArrayInterface", "DocStringExtensions", "Lazy", "LinearAlgebra", "Setfield", "SparseArrays", "StaticArraysCore", "Tricks"]
git-tree-sha1 = "51ae235ff058a64815e0a2c34b1db7578a06813d"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.7"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SparseInverseSubset]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "52962839426b75b3021296f7df242e40ecfc0852"
uuid = "dc90abb0-5640-4711-901d-7e5b23a2fada"
version = "0.1.2"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "bf074c045d3d5ffd956fa0a461da38a44685d6b2"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.3"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "1d77abd07f617c4868c33d4f5b9e1dbb2643c9cf"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.2"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "cef0472124fab0695b58ca35a77c6fb942fdab8a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.1"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a04cabe79c5f01f4d723cc6704070ada0b9d46d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.4"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "f4dc295e983502292c4c3f951dbb4e985e35b3be"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.18"
weakdeps = ["Adapt", "GPUArraysCore", "SparseArrays", "StaticArrays"]

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = "GPUArraysCore"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.SymbolicIndexingInterface]]
git-tree-sha1 = "be414bfd80c2c91197823890c66ef4b74f5bf5fe"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.1"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "TimerOutputs", "Unityper"]
git-tree-sha1 = "669e43e90df46fcee4aa859b587da7a7948272ac"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "1.5.1"

[[deps.Symbolics]]
deps = ["ArrayInterface", "Bijections", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "Groebner", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "NaNMath", "PrecompileTools", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "SymbolicUtils"]
git-tree-sha1 = "8d28ebc206dec9e250e21b9502a2662265897650"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "5.14.1"

    [deps.Symbolics.extensions]
    SymbolicsSymPyExt = "SymPy"

    [deps.Symbolics.weakdeps]
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "cb76cf677714c095e535e3501ac7954732aeea2d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.11.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "f548a9e9c490030e545f72074a41edfd0e5bcdd7"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.23"

[[deps.TranscodingStreams]]
git-tree-sha1 = "14389d51751169994b2e1317d5c72f7dc4f21045"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.6"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "3c793be6df9dd77a0cf49d80984ef9ff996948fa"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.19.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "25008b734a03736c41e2a7dc314ecb95bd6bbdb0"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.6"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "07e470dabc5a6a4254ffebc29a1b3fc01464e105"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.5+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "31c421e5516a6248dfb22c194519e37effbf1f30"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.6.1+0"

[[deps.Xorg_libICE_jll]]
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "e5becd4411063bdcac16be8b66fc2f9f6f1e8fe5"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.0.10+1"

[[deps.Xorg_libSM_jll]]
deps = ["Libdl", "Pkg", "Xorg_libICE_jll"]
git-tree-sha1 = "4a9d9e4c180e1e8119b5ffc224a7b59d3a7f7e18"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.3+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.Zygote]]
deps = ["AbstractFFTs", "ChainRules", "ChainRulesCore", "DiffRules", "Distributed", "FillArrays", "ForwardDiff", "GPUArrays", "GPUArraysCore", "IRTools", "InteractiveUtils", "LinearAlgebra", "LogExpFunctions", "MacroTools", "NaNMath", "PrecompileTools", "Random", "Requires", "SparseArrays", "SpecialFunctions", "Statistics", "ZygoteRules"]
git-tree-sha1 = "4ddb4470e47b0094c93055a3bcae799165cc68f1"
uuid = "e88e6eb3-aa80-5325-afca-941959d7151f"
version = "0.6.69"

    [deps.Zygote.extensions]
    ZygoteColorsExt = "Colors"
    ZygoteDistancesExt = "Distances"
    ZygoteTrackerExt = "Tracker"

    [deps.Zygote.weakdeps]
    Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
    Distances = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ZygoteRules]]
deps = ["ChainRulesCore", "MacroTools"]
git-tree-sha1 = "27798139afc0a2afa7b1824c206d5e87ea587a00"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.5"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a68c9655fbe6dfcab3d972808f1aafec151ce3f8"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.43.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.7.0+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╟─08cf2a80-f2b2-11ea-132d-2f26c934bbb2
# ╟─b72a9871-3773-4a06-9a2f-09ce5a7d8496
# ╠═47f98a52-f2b2-11ea-0ac6-9b4b74b7b15d
# ╟─5154a672-f2b2-11ea-1280-3ba5dadff7f9
# ╟─f53fbf24-f2b2-11ea-2f59-d943884fee9d
# ╟─e6a55456-fe6f-11ea-328c-47fedb388194
# ╠═2be8656c-f2b3-11ea-2244-5740c807deff
# ╠═2e00ef04-f2b3-11ea-1960-8136321c71db
# ╟─f7cd052e-f2bb-11ea-239a-47966e96e909
# ╠═7244eba2-f2b3-11ea-045f-ffd208852eda
# ╟─60058dbe-1514-4139-bb97-66b40fe47306
# ╠═21882f0e-8d54-41ae-84c6-7637ccc52c46
# ╟─401df7a4-f2b3-11ea-12a1-e9b46bcfcf1b
# ╠═5ac4989c-f2b3-11ea-3e8a-a79a4117e48c
# ╠═5db7350a-f2b3-11ea-293d-ffc572e0dbe8
# ╠═d5edd644-fe7c-11ea-2aa8-c772ea2f1bde
# ╠═76ba97fe-b7c4-467b-85ae-9067ebcc8b14
# ╠═64b28e4a-f2b3-11ea-0fc0-7d37b4794686
# ╠═8a4a8499-af01-46ea-8196-fb2d718d0d17
# ╟─92b65737-9ad7-422e-955c-e9adacd2032b
# ╠═a38e2cce-851c-4065-8947-1926fa9ca5b4
# ╠═684ff7f4-f2b3-11ea-14f7-7d3f329ecd79
# ╟─2ca60bca-f2b4-11ea-0fec-c50f4e306536
# ╠═7a6a18de-f2b3-11ea-1b3d-ef7a2b1b0434
# ╟─67e7ea4d-0b76-4247-b22d-e20c8280f046
# ╟─45190ad2-d358-455a-89a1-4f9c3b3880fa
# ╠═732ad756-1707-48e7-9bfb-10707e64b118
# ╠═c8a5c6f8-4ca1-4c7e-9cc1-d2e5c9e9da59
# ╟─c67f8b3f-0561-4f47-ade8-3c79d7920ae1
# ╟─997ec642-a865-4c97-9aa1-d4d9df9e1dd0
# ╟─913234a2-f2b3-11ea-3501-e3bfc54f91fe
# ╟─0238579a-fe70-11ea-0968-0f11ae9557db
# ╠═9ecd1cd6-8412-4fe1-9d79-676f9830ec2c
# ╠═af02692c-c1be-4f7c-a09f-fbd88b0c4bd8
# ╠═eef066e0-f2b3-11ea-3cfa-5f35ebf13f76
# ╠═f2158710-f2b3-11ea-16c9-dded5e1747d4
# ╠═f54ee2b2-f2b3-11ea-24e8-f556ab9e3b85
# ╠═08a565e0-f2b4-11ea-2ab6-f5f3a1b3c22d
# ╠═1055a996-f2b4-11ea-14d2-6116cbbaa6ba
# ╠═13db037a-f2b4-11ea-2a38-bd7c1c60f7d9
# ╟─4f42234e-f2b4-11ea-2c86-93449b8dc60b
# ╠═6adbc7fe-f2b4-11ea-1861-b398d517edfa
# ╠═6cd1accc-f2b4-11ea-0077-e7685b77b1be
# ╟─77d4fd72-f2b4-11ea-0bd6-97fc2299ca2f
# ╠═773af248-f2b4-11ea-154a-cb67beebed3d
# ╟─8a5ee764-f2b4-11ea-23f3-27664d910c16
# ╠═8ded9128-f2b4-11ea-0f44-9321077205de
# ╟─c3a1abd8-f2b4-11ea-0509-a5296ae8b361
# ╠═c9f92e66-f2b4-11ea-322f-7778f149453d
# ╟─cfb07cc6-f2b4-11ea-08a7-dd0b13d15b5e
# ╠═d56b74f2-f2b4-11ea-12b4-5f54d4999b38
# ╟─d9f2d114-f2b4-11ea-32b0-43fa1a872a97
# ╠═de1fd700-f2b4-11ea-36dd-4bfba5ecc134
# ╠═df14ce72-f2b4-11ea-307b-8fac2889681b
# ╠═e584f976-f2b4-11ea-0245-67f6bcc4be0e
# ╟─eb68674c-f2b4-11ea-0cc2-f792d17a5005
# ╟─f4edad2c-f2b4-11ea-05fb-3d43c18e9d5c
# ╟─fd9ed266-f2b4-11ea-3c53-3bd2fe30a7f1
# ╟─4e6b2b04-f2b5-11ea-0c52-5f32438e8c68
# ╟─55ff7230-f2b5-11ea-2898-db5e21c764fd
# ╠═688a1f40-f2b5-11ea-3fa8-ebb735323ddd
# ╟─e5c12243-70ad-4e52-9bab-aae432203e21
# ╟─6d6da8ba-f2b5-11ea-1064-df27ecea7d13
# ╟─6fd08e42-f2b5-11ea-36e6-d7ff10412e2b
# ╟─9a87095e-f2b5-11ea-334e-57a637e20c43
# ╟─67e9fa88-ff03-11ea-1177-c34331398d6c
# ╠═75ae48fe-f2b5-11ea-0614-912091cbfb6a
# ╟─8a738272-f2b5-11ea-2de6-f9ae6b30a466
# ╠═7b95d73c-f2b5-11ea-1b29-858f7b6c99bf
# ╠═836ca94a-f2b5-11ea-07e0-bd830f810fba
# ╠═8fe85430-f2b5-11ea-3fe7-5d9a0db03932
# ╠═b372c25c-f2b5-11ea-1c55-e5045e0dc5b7
# ╠═b6e48b12-f2b5-11ea-3bf2-81e5b2a3dc98
# ╠═b455e8ac-fe7f-11ea-1a8a-db85d11ccba7
# ╠═c2924514-fe7f-11ea-2b01-81e9e0e1c491
# ╠═cd2821c4-fe7f-11ea-0c41-8feab49be07a
# ╠═bebbd674-f2b5-11ea-19c6-4fe1257b0277
# ╟─d6c6d4ea-f2b7-11ea-0b41-d7f629cdcd2f
# ╠═fe18341e-f2b7-11ea-0236-03b44eebd328
# ╠═3f2a4cc2-26fa-4c3e-bfb7-f3ef019b6d78
# ╠═9e22156e-0880-11eb-24a8-3d0a41262c36
# ╠═ab39b554-0880-11eb-05de-71949f71b440
# ╟─0aba93a6-f2b8-11ea-1461-d96878b7a5c1
# ╠═014ec92c-f2b8-11ea-3986-8d74a8c2458d
# ╠═4a6270a3-a271-40da-84d0-28247f287c22
# ╠═272bdcca-aeb1-405f-8251-33f13cddbdc5
# ╟─3f5a58da-f2b8-11ea-3ab0-5db5e9199bb1
# ╠═48e5da96-f2b8-11ea-0108-2104f2c8ac24
# ╠═eaacb5a5-c748-4223-98e4-be73d3f2086f
# ╟─263e18dc-fe70-11ea-03a0-8fd5e031b06f
# ╠═2f18c362-f2b8-11ea-1beb-719ff5cd5f8d
# ╠═393f92a6-fe70-11ea-024e-ff30dcfc7103
# ╠═329cb67e-f2b8-11ea-07b1-e13fd8ddbe42
# ╠═36a2ca10-f2b8-11ea-0fae-53b552c76115
# ╠═4ee45daa-f2b8-11ea-10ea-f929853879b1
# ╟─555c8d88-f2b8-11ea-34d8-a17ce87ccdc1
# ╠═59eee760-f2b8-11ea-29f6-732292e2f95d
# ╠═5ccd625e-f2b8-11ea-38c1-098e2f31ce53
# ╠═60b52f6e-f2b8-11ea-2afc-4be2d4ae3e19
# ╠═63d5e0be-f2b8-11ea-16e8-d9528b4e8945
# ╟─6971e7e6-f2b8-11ea-2a76-9dda6b1455b6
# ╠═6fe82e5a-f2b8-11ea-101d-d1ee5f9bd5b1
# ╠═723fe4ca-f2b8-11ea-1440-a9e1cd213b66
# ╠═76284520-f2b8-11ea-077e-43167cc7484c
# ╠═85a5c518-f2b8-11ea-2e12-25d16a12e7c8
# ╟─8a254852-f2b8-11ea-142d-71bb270a5c01
# ╟─91d0e03e-f2b8-11ea-07e5-db056636ac21
# ╠═963c4910-f2b8-11ea-12bb-bb42539b1c93
# ╟─a041fbf8-f2b8-11ea-28ba-0772497c649f
# ╠═9b5db398-f2b8-11ea-00bc-77de99792925
# ╟─bdce5964-f2b8-11ea-03d1-57705d295040
# ╠═c6728784-f2b8-11ea-12b5-d7eb039d77ce
# ╟─e1057afe-fe6e-11ea-27fe-6def4d8e1f2d
# ╠═eb76e252-fe6e-11ea-0939-c51ff3c73827
# ╟─d6d310ee-f2b8-11ea-24e1-3b405f071545
# ╠═d3ae2c14-f2b8-11ea-1ac2-fb5e0b4b8500
# ╟─dbb1e68c-f2b8-11ea-1220-6f5715822e3d
# ╠═df557ed2-f2b8-11ea-22ee-e1c90a2655bf
# ╟─ee52d2b8-f2b8-11ea-15c4-19e9f2ad10e2
# ╟─213fc38a-f380-11ea-135e-6f72924c19aa
# ╟─059679c0-d0ff-453b-b897-b889bbea326d
# ╠═a20f0158-f32b-45e8-8bd5-6f693eca1860
# ╠═23b0cc81-b551-446b-991f-54bd19ac63cd
# ╠═a9961630-08be-48f8-b1cc-4a9d67dfc519
# ╠═e4d2350e-147b-4c3c-a9ff-6de66a54c48a
# ╟─501a0db3-538a-43f3-a769-158c824bca31
# ╠═ca64df96-082a-40dd-b82c-64016257fcab
# ╠═ba43bc8b-27d6-490d-971e-7e440b991cfb
# ╟─2d7dbec3-82eb-41f8-a959-b458e136d18c
# ╠═378500c1-9f19-4111-82cf-5a57c30d66ef
# ╟─8586855b-2ad7-4fe4-af1a-2a6f1dd052b3
# ╟─8877c563-d10a-41e4-be8b-546d27623b67
# ╠═a73e4d59-34f4-412f-afa4-f3b66244b35e
# ╠═6707fd27-78f6-4663-b95d-a016d2c44c3b
# ╟─18e151c4-0b0e-485b-81d2-f18a2dd4592d
# ╟─79221eee-f2be-11ea-3e43-4596f5ce5694
# ╠═fb91a9c4-f2c0-11ea-38c4-5b554446e997
# ╠═07c45d78-f2c3-11ea-131e-b9f0b67ca3cd
# ╠═33023d10-f2c1-11ea-3375-e9ac23172ab6
# ╟─31cdefd2-f381-11ea-3037-e9c023ebc1e9
# ╠═4def0626-f2c1-11ea-2f01-31a8c42b121c
# ╟─46c54108-f381-11ea-141f-31a696a7a8ef
# ╟─5812f84c-f2c1-11ea-0f6c-31b9fedc3704
# ╠═8049560e-f383-11ea-1c8e-39c629d20982
# ╠═c7b5a1aa-f383-11ea-2a84-8da9c9e60dbc
# ╠═cfad1a8c-f383-11ea-3dad-214f3da9485d
# ╟─83ddbe6a-f381-11ea-2415-255cbfb2590b
# ╠═30e9891c-f2c5-11ea-022d-5dbf9c064e43
# ╟─b07cb958-f381-11ea-3378-5517440464ca
# ╠═eaac2514-f8fe-11ea-2229-81f5b59f58fe
# ╠═fb044216-f381-11ea-2faf-fd53dad74950
# ╠═675b222c-f382-11ea-257a-afcb5a5afdd8
# ╠═818c1f9e-ff1b-11ea-097f-cba2f30af994
# ╠═7aff4178-bcac-4cfa-ab22-de3506cf73a8
# ╟─881dedd4-f385-11ea-342e-a3700e90ac66
# ╠═c1e19e26-f385-11ea-17b6-6d2d2e5cb7c0
# ╠═fe18a0a0-f8fe-11ea-1595-4d95256dabca
# ╠═ec87028b-d9d4-4dc7-9dd5-36eaa1dd49d4
# ╟─df5f20d8-f385-11ea-104d-9b9e7b094eb7
# ╠═e4a088f0-f385-11ea-00da-b12e40e4357f
# ╟─ddd214d8-608a-4123-973f-b331c6ab331a
# ╠═8f3fecac-f8fd-11ea-136d-1513f60ee3ed
# ╠═f8863a18-f385-11ea-0965-0f074e58d9b2
# ╟─4c3817f6-396e-4485-a7dd-0f507cec124a
# ╠═78b72cb1-52f3-4dc3-a721-b8a89a24312e
# ╠═7180e310-a6f7-4feb-aff3-28e1c772b441
# ╟─02b42a67-7ebe-4ca2-ac55-aa481b30e997
# ╠═043e3e2f-3876-41c0-bf70-5ef53b5db986
# ╠═8ba368c9-3b04-469e-9556-1573cbc52185
# ╠═4da4a1b6-b7c8-432d-86ed-c45a499a1880
# ╟─c573ba19-4be5-4afd-b9f5-8b75e9b18ad6
# ╠═639528b1-ecff-4a05-ae56-9a08b779edc1
# ╠═838292fc-2c0e-4bb9-9df9-6cd197c4a4fd
# ╠═ca20cc75-3033-4bd1-bb5c-cae61c855de7
# ╟─34d48ec0-f2b9-11ea-147e-d3554a9fcc04
# ╟─0a576c65-3685-4987-93f1-980ace78d6f4
# ╟─91aff94e-e77e-4810-82ad-fe5ab3eafdad
# ╟─8eb9ed90-981b-4094-baec-a521614fb3a1
# ╟─d682439b-bc6b-456d-8b5d-2489d095e55f
# ╟─2d4c0bf6-3936-449d-8341-56f9c5659a6e
# ╟─fe4dfa5c-ea26-4f7c-aa68-7b3ec5f45bdc
# ╠═c25c9a89-69bd-4dd8-b710-5e66f68346b4
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
