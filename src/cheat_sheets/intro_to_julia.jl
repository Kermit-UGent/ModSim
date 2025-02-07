### A Pluto.jl notebook ###
# v0.20.4

#> [frontmatter]
#> order = "2"
#> title = "Intro to Julia"
#> date = "2025-01-28"
#> tags = ["cheat sheets"]
#> description = "General introduction to Julia"
#> layout = "layout.jlhtml"
#> 
#>     [[frontmatter.author]]
#>     name = "Daan Van Hauwermeiren"
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

# ╔═╡ 9a1cca86-0fff-4f8a-a033-8cf6db337503
begin
  using Pkg
  Pkg.activate("../../pluto-deployment-environment")
  using PlutoUI, Markdown; TableOfContents()
end

# ╔═╡ a15c46f7-b561-41c7-8da2-12f04465fb19
using LinearAlgebra  

# ╔═╡ 8f785f5f-163c-436e-8781-0b68c3aa422a
using Statistics

# ╔═╡ bf1385da-4ac2-11eb-3992-41abac921370
using Plots

# ╔═╡ fd21a9fa-4ab9-11eb-05e9-0d0963826b9f
md"""
# Notebook 1: Getting up and running

First of all, welcome to the course! We hope you enjoy the ride.
""" 

# ╔═╡ 2f210d77-c1fe-4854-b8cb-2c33dcf64292
md"""
## 0. Welcome to Pluto

We will do our exercises in the Pluto notebook environment. The Pluto notebooks are pure Julia alternatives to the Jupyter notebooks you might have worked with. They are fast and reactive and come equipped with their own package manager, making it easy to distribute them.

Cells are immediately executed in order of their *dependencies*, so not in the order that they appear. This can be confusing at first.

"""

# ╔═╡ fe9353dd-f2a5-4597-962d-df267c3c70b8
one = 2  # change me and everything is updated!

# ╔═╡ e4f67601-03c5-4920-a1fa-de4ec5947868
two = one / 2  # I depend on one, so I am executed second

# ╔═╡ 8dc1e292-43a1-459e-b69d-9ef620e121a1
three = two + one   # I am executed last

# ╔═╡ f7c7d372-53eb-402e-a0a9-f7a9a5bdb4de
md"To run a cell either press on the ▶ symbol or press `shift + ENTER`."

# ╔═╡ fdc8ef8f-6a52-4374-a343-4f44a419639f
# Hi there, you are not supposed to see me!

md"Press on the little 👁️ left of the cell to toggle a hidden cell."

# ╔═╡ 27e432a0-cf43-4ab4-b51d-bed405e9f791
md"You can only have one statement per line and all your variables need to have unique names. You can split statements in two lines or wrap them in a `begin ... end` block."

# ╔═╡ 264b75e7-2394-413f-b49a-333cb7bb20af
a = 2
a + 7

# ╔═╡ b1149a67-47cb-4596-9fd4-be3a8b136755
md"Pluto might seem strange at first, though its restrictions make it very flexible and allows to easily create interactivity!"

# ╔═╡ 9ef414d6-0b18-44f2-aec8-aa6d87ce9cde
md"period $(@bind period Slider(0.0:0.2:5, default=1, show_value=true))"

# ╔═╡ 4caf69d1-0788-4a55-8822-a5c2c54cca0b
mycos = t->cos(period * t);

# ╔═╡ 4e5a1441-5adf-43c5-b644-19125c1295ac
plot(mycos, 0, 8, xlab="t", title="Cos($period * t)")

# ╔═╡ 23d3c9cc-4abd-11eb-0cb0-21673effee6c
md"""## 1. The basics
*From zero to newbie.*
"""

# ╔═╡ 62c3b076-4ab7-11eb-0cf2-25cdf7d2540d
md"""
Let's get started with the basics. Some mathematical operations, """

# ╔═╡ 7bf5bdbe-4ab7-11eb-0d4b-c116e02cb9d9
1 + 2       # adding integers

# ╔═╡ 83306610-4ab7-11eb-3eb5-55a465e0abb9
1.0 + 2.0   # adding floats

# ╔═╡ 3fa0a15c-5008-11eb-13b5-a91b02c1eb2d
1 + 2.0     # adding a float to an integer...

# ╔═╡ 83311b8a-4ab7-11eb-0067-e57ceabdfe9d
2 / 4       # standard division

# ╔═╡ 833dbc66-4ab7-11eb-216d-f9900f95deb8
div(2, 4)   # Computes 2/4 truncated to an integer

# ╔═╡ 8342c042-4ab7-11eb-2136-497fc9e1b9c4
2 ÷ 4       # looks nicer but does exactly the same!

# ╔═╡ 619430de-2a22-4311-88a7-1d9a75b4e2b8
2^8         # raising to a power

# ╔═╡ 834d4cbc-4ab7-11eb-1f1a-df05b0c00d66
7 % 3       # get the remainder of the integer division

# ╔═╡ 8360ffac-4ab7-11eb-1162-f7a536eb0765
35 \ 7      # inverse division

# ╔═╡ 8365cb3e-4ab7-11eb-05c0-85f51cc9b018
1 // 3      # fractions, gives the result as a rational

# ╔═╡ 8370eaf0-4ab7-11eb-1cd3-dfeec9341c4b
1//2 + 1//4

# ╔═╡ 50bb93e6-5a6c-11eb-0a6c-d5d749857771
2.0 + 3.0im  # complex numbers

# ╔═╡ 8383f104-4ab7-11eb-38a5-33e59b1591f6
'c'        # characters (unicode)

# ╔═╡ 8387934a-4ab7-11eb-11b2-471b08d87b31
:symbol    # symbols, we will use this for parameters

# ╔═╡ 8bab2e50-5a6c-11eb-3c5f-a9f811483814
:ζ         # any LaTeX symbol

# ╔═╡ 9d2708ca-5a6c-11eb-1c0f-473f0e2b5363
:🎉        # or Unicode emoji

# ╔═╡ 8c14cb9a-4ab7-11eb-0666-b1d4aca00f97
md"variable assignment"

# ╔═╡ 93b5a126-4ab7-11eb-2f67-290ed869d44a
x = 2

# ╔═╡ 353efeea-6492-11eb-3d09-353d4dae491a
md"In the Pluto notebook environment you are currently working in, it is not possible to define the same variable in two cells. However, this is not standard Julia behaviour. You can see that redefining a variable is possible,"

# ╔═╡ 92f57780-6492-11eb-1264-bbef04a8ae99
begin
	variable1 = 2.0
	variable1 = 4.0
end;

# ╔═╡ a17ebeba-6492-11eb-07ba-f516a990affc
variable1

# ╔═╡ 3055711a-6493-11eb-252b-7f5d99115551
md"""
```julia
begin
 statement1
 statement2
end
```

Enable to wrap multiple statements, since only single-line statements are allowed in this notebook environment.

"""

# ╔═╡ f474d864-28b9-4299-b207-dca426554c46
md"Similarly, `let ... end` blocks allow you to define a separate envirionment. Everything you define in such a block is only available there."

# ╔═╡ c6ae05d2-7e8e-4683-b3e2-fe79d5e24e2f
let
	myprivatevar = 3.0
end

# ╔═╡ 9ce5dfc4-5715-4037-bf36-9e5f68c3273e
myprivatevar  # only available in the block

# ╔═╡ 962ae6d2-4ab7-11eb-14a2-c76a2221f544
τ = 1 / 37  # unicode variable names are allowed

# ╔═╡ 98d48302-4ab7-11eb-2397-710d0ae425f7
md"""

unicode! In most Julia editing environments, unicode math symbols can be typed when starting with a '\' and hitting '[TAB]'.

"""

# ╔═╡ acb3b57a-661c-11eb-2c6a-99793a47ff29
md"""> Unsure what the LaTeX name for a symbol is or how to type an emoiji? Just copy-paste it in the REPL with a `?` at the beginning, e.g., `?ζ` and it will tell you how to type it."""

# ╔═╡ cee8a766-4ab7-11eb-2bc7-898df2c9b1ff
# type \alpha  and <TAB>

# ╔═╡ e2c5b558-4ab7-11eb-09be-b354fc56cc6e
md"Operators are not needed for multiplication."

# ╔═╡ ec754104-4ab7-11eb-2a44-557e4304dd43
5x         # This works

# ╔═╡ f23a2d2a-4ab7-11eb-1e26-bb2d1d19829f
md"But strings are quite essential,"

# ╔═╡ fa836e88-4ab7-11eb-0ba6-5fc7372f32ab
mystery = "life, the universe and everything"

# ╔═╡ 0138ef46-4ab8-11eb-1813-55594927d661
md"and string interpolation is performed with `$`."

# ╔═╡ 0b73d66a-4ab8-11eb-06e9-bbe95285a69f
"The answer to $mystery is $(3*2*7)"

# ╔═╡ 6b6eb954-4ab8-11eb-17f9-ef3445d359a3
md"""
Printing can be done with `println()`. These Pluto notebooks distinguish the value of the evaluation or computation from what is printed. The latter is shown in a terminal.
"""

# ╔═╡ 94e3eb74-4ab8-11eb-1b27-573dd2f02b1d
println("The answer to $mystery is $(3*2*7)")

# ╔═╡ abf00a78-4ab8-11eb-1063-1bf4905ca250
md"""
repetitions of strings can be done using the operators `*` and `^`.
This use of `*` and `^` makes sense by analogy with multiplication and exponentiation. Just as `4^3` is equivalent to `4*4*4`, we expect `"Spam"^3` to be the same as `"Spam"*"Spam"*"Spam"`, and it is.
"""

# ╔═╡ be220a48-4ab8-11eb-1cd4-db99cd9db066
breakfast = "eggs"

# ╔═╡ cadaf948-4ab8-11eb-3110-259768055e85
abetterbreakfast = "SPAM"

# ╔═╡ cadb506e-4ab8-11eb-23ed-2d5f88fd30b0
breakfast * abetterbreakfast

# ╔═╡ caf56346-4ab8-11eb-38f5-41336c5b45a7
breakfast * abetterbreakfast^3 * breakfast

# ╔═╡ 046133a8-4ab9-11eb-0591-9de27d85bbca
md"""
Lots of handy `String`-operations are available in the standard library of Julia:
"""

# ╔═╡ 1f255304-4ab9-11eb-34f1-270fd5a95256
md"Unlike `Strings`, a `Char` value represents a single character and is surrounded by single quotes."

# ╔═╡ 34a18900-4ab9-11eb-17a0-1168dd9d06f9
'x'

# ╔═╡ 50d3f9ec-7d5a-4abd-94ba-5c9b0850fdb0
md"Similarly to Matlab, when using the REPL, Julia will print the result of every statement by default. To suppress this behaviour, just end the statement with a semicolon."

# ╔═╡ 15f8b7fe-4abd-11eb-2777-8fc8bf9d342e
var1 = 10;  # not printed...

# ╔═╡ efae58fa-5008-11eb-32fe-c3ae588d14f2
var1  # ...but still defined

# ╔═╡ 18f99e46-4abd-11eb-20a8-859cb1b12fe3
var2 = 20

# ╔═╡ b0893d91-a6b1-4742-8011-5df89fcc558e


# ╔═╡ 3a7954da-4abd-11eb-3c5b-858054b4d06b
md"""## 2. Logical statements

*From zero to one.*
"""


# ╔═╡ 8b17d538-4abd-11eb-0543-ab95c9548d6f
md"""**Boolean operators**

Julia uses `true` and `false` for Boolean variables.
"""

# ╔═╡ 29d34e64-5009-11eb-3301-f729150e17b2
I💖Julia = true 

# ╔═╡ 91a9d1a0-4abd-11eb-3337-71983f32b6ae
!true

# ╔═╡ 942d4202-4abd-11eb-1f01-dfe3df40a5b7
!false

# ╔═╡ 942dae0e-4abd-11eb-20a2-37d9c9882ba8
1 == 1

# ╔═╡ 943d9850-4abd-11eb-1cbc-a1bef988c910
2 == 1

# ╔═╡ 943de2ce-4abd-11eb-2410-31382ae9c74f
1 != 1

# ╔═╡ 9460c03c-4abd-11eb-0d60-4d8aeb5b0c1d
2 != 1

# ╔═╡ 946161f4-4abd-11eb-0ec5-df225dc140d0
1 < 10

# ╔═╡ 947d143a-4abd-11eb-067d-dff955c90407
1 > 10

# ╔═╡ 947fea8e-4abd-11eb-1d6a-2bc540f7a50e
2 <= 2  # or 2 ≤ 2  (\le<TAB>)

# ╔═╡ 948eff10-4abd-11eb-36d0-5183e882a9e2
2 >= 2  # or 2 ≥ 2  (\ge<TAB>)

# ╔═╡ 948f5032-4abd-11eb-3d1c-7da4cb64521c
# Comparisons can be chained
1 < 2 < 3

# ╔═╡ 94b520e6-4abd-11eb-3161-addf3b0e4f24
2 < 3 < 2

# ╔═╡ 94b78322-4abd-11eb-3006-454548efd164
# Logical operators
true && true

# ╔═╡ 94d28c80-4abd-11eb-08c0-717207e4c682
true || false

# ╔═╡ 9fe6e1a2-4abd-11eb-0c39-458ce94265c0
md"Likewise, we have the Boolean logic operators `&&` (AND), `||` (OR) and `⊻` (XOR, exclusive or)."

# ╔═╡ ae26ab9e-4abd-11eb-3270-33558dbdf663
true && true

# ╔═╡ b08dc886-4abd-11eb-1807-096a7e6fd6f9
true && false

# ╔═╡ b08e3a28-4abd-11eb-258a-a5a93b4b882c
true || false

# ╔═╡ b0a8dfe0-4abd-11eb-167d-2fc3974c7c92
false || false

# ╔═╡ b0a97e00-4abd-11eb-371c-e138aea17bb6
true ⊻ false

# ╔═╡ b0ccc252-4abd-11eb-048b-4bec3750bbf1
true ⊻ true

# ╔═╡ 60b066d8-5009-11eb-3b4c-8b8fa2f4831d
md"""
Chaining logic operators is frequently done in Julia as a short alternative for an `if` statement. The idea is if you use an `&&` statement, the second part is only evaluated if the first part is true! The inverse is true for `||`, where the second part is only evaluated if the first part is false."""

# ╔═╡ ec8744ba-000f-4435-a5f6-40ac83b4baa0
md"""
## 3. Vectors and matrices

Julia has powerful, flexible interfaces for vectors, matrices, and higher-order tensors.

A vector is defined with square brackets with the elements separated with a ",":
"""

# ╔═╡ 9aa0e36f-a096-461c-a766-acc88c6c7a92
v = [1, 3, 2]  # a vector of integers

# ╔═╡ b22afc5a-6fb3-4d60-b5ed-daef03c75ac0
v2 = [1.0, 2.0, 3.0]  # a vector of floats

# ╔═╡ 4b81b713-995c-416f-b80f-4f741574d906
v3 = [1.0, 2, 3]   # promotion occurs automatically to the most general type

# ╔═╡ c1eafd30-9368-46b0-b5ed-5a1d36fa2b3f
md"Matrices can also be defined with spaces to separate elements in rows and semicolumns to separate rows."

# ╔═╡ 3d58cc69-6b31-4ce0-8917-24104b72fe9a
A = [5 4 9; 1 2 7; 8 6 3]

# ╔═╡ 04c9a014-5d19-4f9b-adea-353637aca211
md"You can use spaces and brackets to combine matrices and vectors:"

# ╔═╡ 3abcee32-7df3-4420-b96e-b7fc0fbdac6e
[A v; 0 1 2 3]

# ╔═╡ 060fa76e-a62f-46ce-a1b9-a44455f1e5bc
md"Indexing is via **square brackets** (like python!) and the **index starts from 1** (like in Matlab and R)."

# ╔═╡ 53c40e65-3ab9-4e2b-bd20-0f2366476d4a
v[1]  # first element

# ╔═╡ aa438710-28db-4ee2-8638-0682e7adbbd3
v[0]  # does not exist... 

# ╔═╡ 393deb83-133d-4a51-8d61-9f5ac5bc6a26
v[end]  # last element

# ╔═╡ 4be43883-7447-4702-943c-5c0e96d17db9
A[2,3]  # two indices for matrices

# ╔═╡ e1704cec-2189-4f1a-913d-ffba9d2d6b16
A[1,:]  # first row

# ╔═╡ 4665e009-0eda-4fd0-b70b-29fa79a096fd
A[:,2]  # second column

# ╔═╡ ec6e9403-5579-41f1-bd48-8e6aab88e413
A[A.<5]  # conditional indexing, notice the "."

# ╔═╡ a43dff92-6bce-4653-99c4-e8fe4c00cab0
md"Many functions exist that process collections."

# ╔═╡ 55bdff5d-5779-426a-a2f8-7377e0ba9301
sum(A)

# ╔═╡ 2f355db5-b34f-4d93-80dd-2a8bbe4faad4
size(A)

# ╔═╡ 0adb146b-99f7-47f5-a5e8-791ceba25c15
sum(A, dims=1)  # sum over the rows

# ╔═╡ 78139b6f-76c4-42f2-8011-0d56122951bc
sort(v)

# ╔═╡ 90776288-e74e-4a62-9ed4-7e69e966e433
size(v)

# ╔═╡ 6ac8c154-5fb4-4aff-a835-512e8b74a1bb
length(v)

# ╔═╡ 3ba209f8-73d3-4752-9830-496fb301480f
count(isodd, A)  # count the number of odd elements in A

# ╔═╡ a7bb865d-9936-45da-8870-4ca89989f634
sum(sqrt, A)  # sum_ii √A_ij

# ╔═╡ 71a50984-f291-479d-aebf-f7e9921c050d
md"Many more advanced functions are available, for example linear algebra:"

# ╔═╡ 9b09eacc-aa4c-416b-90ed-8b201dd2a593
det(A)

# ╔═╡ f07f7808-8e22-4758-b8bb-7c8bc048d3c7
norm(v)

# ╔═╡ 1e83964c-a287-40fd-a634-ef4614e2def3
eigen(A)

# ╔═╡ 52095242-e845-4d44-b7ab-e1ee125bd7f7
mean(A)

# ╔═╡ 2839397e-06ca-47a6-bfb5-90bcd36b7ece
std(v)

# ╔═╡ 607a593b-d944-4894-927a-f36b9443ca38
md"Finally, there will be useful range objects, that define a linear range between begin and end values."

# ╔═╡ abc7f86d-8350-4531-8c4a-53517cba7ca9
1:100  # from 1 to 100

# ╔═╡ 5e70caa7-5718-4ddb-ac40-1e5795e516ac
0:0.1:1  # from 0 to 1 in steps of 0.1

# ╔═╡ 5f246dd0-7b9e-4c1c-9089-bb61f1100ffe
md"These work just like vectors."

# ╔═╡ 3f49f4fa-eeda-44bf-b1b7-d4ab7a3dc7c0
myrange = 10:0.2:809

# ╔═╡ 33b1455d-da89-40b6-8fba-c9cbb903c62a
myrange[87]

# ╔═╡ 64a05bbe-7325-4fe8-927c-08963977aa03
sum(myrange)

# ╔═╡ 54f38a8c-ed5b-4a1c-bdfb-54b99da357de
length(myrange)

# ╔═╡ fdb67aba-4ac0-11eb-1d4a-c354de54baa9
md"""## 4. Functions
Julia puts the fun in functions. User-defined functions can be declared as follows,
"""

# ╔═╡ 28f47a24-4ac1-11eb-271f-6b4de7311db3
function square(x)
  result = x * x
  return result
end

# ╔═╡ ae45c122-eff9-4cdc-aa50-9de1f95463f0
square(8)

# ╔═╡ 47338c78-4ac1-11eb-04d6-35c2361eaea6
md"Many of the functions we will need will be fairly simple equations. We can just define them in one line. A more condensed version of `square(x)`."

# ╔═╡ 463689b0-4ac1-11eb-1b0f-b7a239011c5c
s(x) = x * x

# ╔═╡ 53dd5e65-d4f4-4f56-9ebd-7bad62ebd716
s(8)

# ╔═╡ 0dbe0c34-500e-11eb-2633-67d8dc6b24c8
md"""
Functions are first-class and work just like any other variable! For example, you can give a function as an input in another function. 

In some cases, you might want to define an **anonymous function**, without giving them a name:"""

# ╔═╡ b531a7fa-0395-488a-8859-bd1d4bc3a14e
anfun = x -> x^2 - 2x - 8

# ╔═╡ 1380fc18-889c-4c71-b087-71e3c430a81d
md"This looks like a variable but can be used as a function:"

# ╔═╡ 9e89852e-8cf4-47f9-9712-5f207a2444d9
anfun(1.5)  # works just like any function

# ╔═╡ 42e76b94-d906-4706-a048-89937294b757
md"Why do we need this? Because we might want to define small functions on the fly."

# ╔═╡ cf58e15e-025a-4160-b04b-575e821efdee
count(x-> 4 < x^2 < 80, -100:100)  # count the numbers between -100 and 100, for which their square is between 4 and 80

# ╔═╡ bf8f49a0-c0bf-40cd-8bd0-f9ab8802ba09
md"""
> Complete the function `clip(x)`, which returns `x` if $0\le x \le 1$, `0` if $x<0$ and `1` if $x>1$.

"""

# ╔═╡ 0c693c24-4ac0-11eb-2329-c743dcc5039d
clip(x) = missing

# ╔═╡ 119d17c2-cd47-4827-9907-fe3d9f4e5a65
md"By default, a function is over the whole object. Using a `.`, you can use the function element-wise."

# ╔═╡ 2dc95df5-e4f7-4053-ae79-4125002e925b
square(A)  # A * A

# ╔═╡ 65ee39ab-3c2b-41fc-bc78-a082eb61089a
square.(A)  # each element squared

# ╔═╡ 49d3878a-3d70-48ff-94a5-ee2bac3844b1
square(v)  # square does not work for vectors

# ╔═╡ f038960c-6074-4216-8186-89275c207e3e
square.(v)  # element-wise works

# ╔═╡ 31b7328c-126b-45ab-a865-3125009ba9b2
exp(A)  # matrix exponential

# ╔═╡ 25dcf9e5-010e-40c9-a0d8-7488e9e94203
exp.(A)  # element-wise exponential

# ╔═╡ e4564f72-7bf6-47f2-a8c6-ca41c81146ab
A + 1  # won' t work

# ╔═╡ e8b365f5-0893-47e6-a213-d14eedba9070
A .+ 1  # add 1 to each element of A

# ╔═╡ 1c22b880-4abf-11eb-3f18-756c1198ccad
md"## 5. Control flow"

# ╔═╡ 37086212-4abf-11eb-3ec9-7f8dae57121e
md"The `if`, `else`, `elseif`-statement is instrumental to any programming language. Note that control flow is ended with an `end` statement. In constrast to Python, tabs are only for clarity but do not impact functionality."

# ╔═╡ 489421d8-4abf-11eb-0d5e-fd779cc918a1
if 4 > 3
  'A'
elseif 3 > 4
  'B'
else
  'C'
end

# ╔═╡ 2a5fca7c-4ac0-11eb-33a3-23d972ca27b8
md"## 6. Looping

Looping using a `for` loop can be done by iterating over a list or range. Don't forget to end with an `end` at the end.
"

# ╔═╡ 1ab81428-8df4-405d-ad60-d1c6e3f6e834
for i in 1:10
	println("$i squared = $(s(i))")
end

# ╔═╡ 3896642a-4ac0-11eb-2c7c-4f376ab82217
characters = ["Harry", "Ron", "Hermione"]

# ╔═╡ 3ef3faf8-4ac0-11eb-1965-fd23413e29f3
begin
	for char in characters
	  println("Character $char")
	end
end

# ╔═╡ 3916f50e-661d-11eb-0829-cb3821836fdf
md"We can use `enumerate` to generate an iterator of tuples containing the index and the values of an iterator."

# ╔═╡ 4118016e-4ac0-11eb-18bf-5de326782c87
begin
	for (i, char) in enumerate(characters)
	  println("$i. $char")
	end
end

# ╔═╡ 4119fbca-4ac0-11eb-1ea9-0bdd324214c5
pets = ["Hedwig", "Pig", "Crookshanks"]

# ╔═╡ 5ebe25c0-661d-11eb-389b-3d81570f7cf0
md"`zip` binds two or more iterators and yields tuples of the pairs."

# ╔═╡ 4139bf3c-4ac0-11eb-2b63-77a513149351
begin
	for (char, pet) in zip(characters, pets)
	  println("$char has $pet as a pet")
	end
end

# ╔═╡ de48a3f6-4f2f-11eb-314b-493546c37a21
 md"## 7. Macros
Macros provide a method to include generated code in the final body of a program. It is a way of generating a new output expression, given an unevaluated input expression. When your Julia program runs, it first parses and evaluates the macro, and the processed code produced by the macro is eventually evaluated like an ordinary expression.

Some nifty basic macros are `@time` and `@show`. `@time` prints the cpu time and memory allocations of an expression."

# ╔═╡ 85b96ff0-4ac2-11eb-077f-cf4aad8a3c24
@time square(10)

# ╔═╡ a11c2898-4ac2-11eb-24d3-6f8060b5fd65
md"""The `@show` macro is often useful for debugging purposes. It displays both the expression to be evaluated and its result, finally returning the value of the result."""

# ╔═╡ a686e67e-4ac2-11eb-228e-23524a3ddc59
@show 1 + 1

# ╔═╡ d50cced2-500d-11eb-2dcc-21fc50825f43
md"Macro's will be vital in the domain specific languages we use in this course. Remember, when you see an `@`, some code is changed into other code."

# ╔═╡ ad156892-4ac2-11eb-3634-a3783231e5a1
md"""## 8. Plotting

Quite essential for scientific programming is the visualisation of the results. `Plots` is the Julia package that handles a lot of the visualisation. `rand(10)` returns an array of 10 random floats between 0 and 1.
"""

# ╔═╡ e5eeb54f-49f4-4a24-a3f4-6d0b0da99c76
plot(rand(10))

# ╔═╡ d779956a-4ac2-11eb-39de-4b3cecace452
md"""When loading in a package for the first time Julia will have to precompile this package, hence this step can take some time."""

# ╔═╡ c7d2a048-4ac2-11eb-3902-b7c8505096ae
begin 
	plot(1:10, rand(10), label="first")
	plot!(1:10, rand(10), label="second")  # adding to current figure using plot!

	scatter!([1:10], randn(10), label="scatter")

	xlabel!("x")
	ylabel!("f(x)")
	title!("My pretty Julia plot")
end

# ╔═╡ cf35b2b2-4ac2-11eb-1ae6-5d3c108210df
plot(0:0.1:10, x -> sin(x) / x, xlabel="x", ylabel="sin(x)/x", color=:red, marker=:square, legend=:none) 
# notice the use of a symbol as an argument !

# ╔═╡ d1010f88-4ac2-11eb-0fa9-0902fef0cf9f
contour(-5:0.1:5, -10:0.1:10, (x, y) -> 3x^2-4y^2 + x*y/6)

# ╔═╡ d43d2244-7b1b-4360-9c5e-729f427fe9a7
md"You can also directly plot functions:"

# ╔═╡ 430b3fd0-58b8-443a-9093-24c085af422d
plot(sin, 0, 2pi)

# ╔═╡ 29b3b33c-92e9-4579-8f59-36e81887e283
md"Don't worry about making a tidy plot. For many objects (solutions of differential equations), the function `plot()` is overloaded, so we only have to `plot(sol)` for a pretty plot. More to follow!
"

# ╔═╡ 5ab03afc-a5fe-4d89-9327-81cc24afd4e0
md"""
> **Exercise: Stirling's approximation for factorials**

The factorial function,

$${\displaystyle n!=1\cdot 2\cdot 3\cdots (n-2)\cdot (n-1)\cdot n,}$$

is often used in combinatorics but also other mathematical areas. Especially for large numbers it can get quite inefficient to compute. Stirling's approximation is an approximation for factorials,

$${\displaystyle n!\sim {\sqrt {2\pi n}}\left({\frac {n}{e}}\right)^{n},}$$
	
Complete the function `stirling()` by implementing Stirling's approximation. 
"""

# ╔═╡ a3969292-57ff-11eb-059b-e9e931a30dc1
stirling(n) = missing

# ╔═╡ 17ef818e-52e6-400e-b839-fe8dd008aef7
md"You can add your approximation to the plot below."

# ╔═╡ dbd89309-f8a6-402f-8703-6f6867a4a19a
scatter(1:10, factorial.(1:10), xlab="n", label="n!", yscale=:log10)

# ╔═╡ cf65bee7-9f09-4b23-9898-07da6fdca98c
begin
	# Do NOT delete this cell!
	
	hint(text) = Markdown.MD(Markdown.Admonition("hint", "Hint", [text]));
	
	almost(text) = Markdown.MD(Markdown.Admonition("warning", "Almost there!", [text]));
	
	keep_working(text=md"The answer is not quite right.") = Markdown.MD(Markdown.Admonition("danger", "Keep working on it!", [text]));
	
	correct(text=md"Great! You got the right answer! Let's move on to the next section.") = Markdown.MD(Markdown.Admonition("correct", "Got it!", [text]));
	
	sol_stirling(n) = √(2π * n) * (n/exp(1))^n;
	
	md"" # Only the last evaluation is shown.
end

# ╔═╡ 060f0de1-1ced-4f10-9137-0267c3572ebc
hint(md"Check out `min`and `max`.")

# ╔═╡ 243d30cb-9072-4eb9-adc0-ff7dcecc2bb5
if !ismissing(clip(0.1)) 
	if clip(-1) == 0 && clip(0.25) ≈ 0.25 && clip(3.6) ≈ 1
		correct()
	else
		keep_working()
	end
end
    

# ╔═╡ 9e211fbf-4c43-439a-95ba-86b48767c053
if !ismissing(stirling(5)) 
	if sol_stirling(20) ≈ stirling(20)
		correct()
	else
		keep_working()
	end
end
    

# ╔═╡ Cell order:
# ╠═9a1cca86-0fff-4f8a-a033-8cf6db337503
# ╟─fd21a9fa-4ab9-11eb-05e9-0d0963826b9f
# ╟─2f210d77-c1fe-4854-b8cb-2c33dcf64292
# ╠═8dc1e292-43a1-459e-b69d-9ef620e121a1
# ╠═e4f67601-03c5-4920-a1fa-de4ec5947868
# ╠═fe9353dd-f2a5-4597-962d-df267c3c70b8
# ╟─f7c7d372-53eb-402e-a0a9-f7a9a5bdb4de
# ╟─fdc8ef8f-6a52-4374-a343-4f44a419639f
# ╠═27e432a0-cf43-4ab4-b51d-bed405e9f791
# ╠═264b75e7-2394-413f-b49a-333cb7bb20af
# ╟─b1149a67-47cb-4596-9fd4-be3a8b136755
# ╟─9ef414d6-0b18-44f2-aec8-aa6d87ce9cde
# ╟─4e5a1441-5adf-43c5-b644-19125c1295ac
# ╟─4caf69d1-0788-4a55-8822-a5c2c54cca0b
# ╟─23d3c9cc-4abd-11eb-0cb0-21673effee6c
# ╟─62c3b076-4ab7-11eb-0cf2-25cdf7d2540d
# ╠═7bf5bdbe-4ab7-11eb-0d4b-c116e02cb9d9
# ╠═83306610-4ab7-11eb-3eb5-55a465e0abb9
# ╠═3fa0a15c-5008-11eb-13b5-a91b02c1eb2d
# ╠═83311b8a-4ab7-11eb-0067-e57ceabdfe9d
# ╠═833dbc66-4ab7-11eb-216d-f9900f95deb8
# ╠═8342c042-4ab7-11eb-2136-497fc9e1b9c4
# ╠═619430de-2a22-4311-88a7-1d9a75b4e2b8
# ╠═834d4cbc-4ab7-11eb-1f1a-df05b0c00d66
# ╠═8360ffac-4ab7-11eb-1162-f7a536eb0765
# ╠═8365cb3e-4ab7-11eb-05c0-85f51cc9b018
# ╠═8370eaf0-4ab7-11eb-1cd3-dfeec9341c4b
# ╠═50bb93e6-5a6c-11eb-0a6c-d5d749857771
# ╠═8383f104-4ab7-11eb-38a5-33e59b1591f6
# ╠═8387934a-4ab7-11eb-11b2-471b08d87b31
# ╠═8bab2e50-5a6c-11eb-3c5f-a9f811483814
# ╠═9d2708ca-5a6c-11eb-1c0f-473f0e2b5363
# ╟─8c14cb9a-4ab7-11eb-0666-b1d4aca00f97
# ╠═93b5a126-4ab7-11eb-2f67-290ed869d44a
# ╟─353efeea-6492-11eb-3d09-353d4dae491a
# ╠═92f57780-6492-11eb-1264-bbef04a8ae99
# ╠═a17ebeba-6492-11eb-07ba-f516a990affc
# ╟─3055711a-6493-11eb-252b-7f5d99115551
# ╟─f474d864-28b9-4299-b207-dca426554c46
# ╠═c6ae05d2-7e8e-4683-b3e2-fe79d5e24e2f
# ╠═9ce5dfc4-5715-4037-bf36-9e5f68c3273e
# ╠═962ae6d2-4ab7-11eb-14a2-c76a2221f544
# ╟─98d48302-4ab7-11eb-2397-710d0ae425f7
# ╟─acb3b57a-661c-11eb-2c6a-99793a47ff29
# ╠═cee8a766-4ab7-11eb-2bc7-898df2c9b1ff
# ╟─e2c5b558-4ab7-11eb-09be-b354fc56cc6e
# ╠═ec754104-4ab7-11eb-2a44-557e4304dd43
# ╟─f23a2d2a-4ab7-11eb-1e26-bb2d1d19829f
# ╠═fa836e88-4ab7-11eb-0ba6-5fc7372f32ab
# ╟─0138ef46-4ab8-11eb-1813-55594927d661
# ╠═0b73d66a-4ab8-11eb-06e9-bbe95285a69f
# ╟─6b6eb954-4ab8-11eb-17f9-ef3445d359a3
# ╠═94e3eb74-4ab8-11eb-1b27-573dd2f02b1d
# ╟─abf00a78-4ab8-11eb-1063-1bf4905ca250
# ╠═be220a48-4ab8-11eb-1cd4-db99cd9db066
# ╠═cadaf948-4ab8-11eb-3110-259768055e85
# ╠═cadb506e-4ab8-11eb-23ed-2d5f88fd30b0
# ╠═caf56346-4ab8-11eb-38f5-41336c5b45a7
# ╟─046133a8-4ab9-11eb-0591-9de27d85bbca
# ╟─1f255304-4ab9-11eb-34f1-270fd5a95256
# ╠═34a18900-4ab9-11eb-17a0-1168dd9d06f9
# ╟─50d3f9ec-7d5a-4abd-94ba-5c9b0850fdb0
# ╠═15f8b7fe-4abd-11eb-2777-8fc8bf9d342e
# ╠═efae58fa-5008-11eb-32fe-c3ae588d14f2
# ╠═18f99e46-4abd-11eb-20a8-859cb1b12fe3
# ╠═b0893d91-a6b1-4742-8011-5df89fcc558e
# ╟─3a7954da-4abd-11eb-3c5b-858054b4d06b
# ╟─8b17d538-4abd-11eb-0543-ab95c9548d6f
# ╠═29d34e64-5009-11eb-3301-f729150e17b2
# ╠═91a9d1a0-4abd-11eb-3337-71983f32b6ae
# ╠═942d4202-4abd-11eb-1f01-dfe3df40a5b7
# ╠═942dae0e-4abd-11eb-20a2-37d9c9882ba8
# ╠═943d9850-4abd-11eb-1cbc-a1bef988c910
# ╠═943de2ce-4abd-11eb-2410-31382ae9c74f
# ╠═9460c03c-4abd-11eb-0d60-4d8aeb5b0c1d
# ╠═946161f4-4abd-11eb-0ec5-df225dc140d0
# ╠═947d143a-4abd-11eb-067d-dff955c90407
# ╠═947fea8e-4abd-11eb-1d6a-2bc540f7a50e
# ╠═948eff10-4abd-11eb-36d0-5183e882a9e2
# ╠═948f5032-4abd-11eb-3d1c-7da4cb64521c
# ╠═94b520e6-4abd-11eb-3161-addf3b0e4f24
# ╠═94b78322-4abd-11eb-3006-454548efd164
# ╠═94d28c80-4abd-11eb-08c0-717207e4c682
# ╟─9fe6e1a2-4abd-11eb-0c39-458ce94265c0
# ╠═ae26ab9e-4abd-11eb-3270-33558dbdf663
# ╠═b08dc886-4abd-11eb-1807-096a7e6fd6f9
# ╠═b08e3a28-4abd-11eb-258a-a5a93b4b882c
# ╠═b0a8dfe0-4abd-11eb-167d-2fc3974c7c92
# ╠═b0a97e00-4abd-11eb-371c-e138aea17bb6
# ╠═b0ccc252-4abd-11eb-048b-4bec3750bbf1
# ╟─60b066d8-5009-11eb-3b4c-8b8fa2f4831d
# ╟─ec8744ba-000f-4435-a5f6-40ac83b4baa0
# ╠═9aa0e36f-a096-461c-a766-acc88c6c7a92
# ╠═b22afc5a-6fb3-4d60-b5ed-daef03c75ac0
# ╠═4b81b713-995c-416f-b80f-4f741574d906
# ╟─c1eafd30-9368-46b0-b5ed-5a1d36fa2b3f
# ╠═3d58cc69-6b31-4ce0-8917-24104b72fe9a
# ╟─04c9a014-5d19-4f9b-adea-353637aca211
# ╠═3abcee32-7df3-4420-b96e-b7fc0fbdac6e
# ╟─060fa76e-a62f-46ce-a1b9-a44455f1e5bc
# ╠═53c40e65-3ab9-4e2b-bd20-0f2366476d4a
# ╠═aa438710-28db-4ee2-8638-0682e7adbbd3
# ╠═393deb83-133d-4a51-8d61-9f5ac5bc6a26
# ╠═4be43883-7447-4702-943c-5c0e96d17db9
# ╠═e1704cec-2189-4f1a-913d-ffba9d2d6b16
# ╠═4665e009-0eda-4fd0-b70b-29fa79a096fd
# ╠═ec6e9403-5579-41f1-bd48-8e6aab88e413
# ╟─a43dff92-6bce-4653-99c4-e8fe4c00cab0
# ╠═55bdff5d-5779-426a-a2f8-7377e0ba9301
# ╠═2f355db5-b34f-4d93-80dd-2a8bbe4faad4
# ╠═0adb146b-99f7-47f5-a5e8-791ceba25c15
# ╠═78139b6f-76c4-42f2-8011-0d56122951bc
# ╠═90776288-e74e-4a62-9ed4-7e69e966e433
# ╠═6ac8c154-5fb4-4aff-a835-512e8b74a1bb
# ╠═3ba209f8-73d3-4752-9830-496fb301480f
# ╠═a7bb865d-9936-45da-8870-4ca89989f634
# ╟─71a50984-f291-479d-aebf-f7e9921c050d
# ╠═a15c46f7-b561-41c7-8da2-12f04465fb19
# ╠═9b09eacc-aa4c-416b-90ed-8b201dd2a593
# ╠═f07f7808-8e22-4758-b8bb-7c8bc048d3c7
# ╠═1e83964c-a287-40fd-a634-ef4614e2def3
# ╠═8f785f5f-163c-436e-8781-0b68c3aa422a
# ╠═52095242-e845-4d44-b7ab-e1ee125bd7f7
# ╠═2839397e-06ca-47a6-bfb5-90bcd36b7ece
# ╟─607a593b-d944-4894-927a-f36b9443ca38
# ╠═abc7f86d-8350-4531-8c4a-53517cba7ca9
# ╠═5e70caa7-5718-4ddb-ac40-1e5795e516ac
# ╟─5f246dd0-7b9e-4c1c-9089-bb61f1100ffe
# ╠═3f49f4fa-eeda-44bf-b1b7-d4ab7a3dc7c0
# ╠═33b1455d-da89-40b6-8fba-c9cbb903c62a
# ╠═64a05bbe-7325-4fe8-927c-08963977aa03
# ╠═54f38a8c-ed5b-4a1c-bdfb-54b99da357de
# ╟─fdb67aba-4ac0-11eb-1d4a-c354de54baa9
# ╠═28f47a24-4ac1-11eb-271f-6b4de7311db3
# ╠═ae45c122-eff9-4cdc-aa50-9de1f95463f0
# ╟─47338c78-4ac1-11eb-04d6-35c2361eaea6
# ╠═463689b0-4ac1-11eb-1b0f-b7a239011c5c
# ╠═53dd5e65-d4f4-4f56-9ebd-7bad62ebd716
# ╟─0dbe0c34-500e-11eb-2633-67d8dc6b24c8
# ╠═b531a7fa-0395-488a-8859-bd1d4bc3a14e
# ╟─1380fc18-889c-4c71-b087-71e3c430a81d
# ╠═9e89852e-8cf4-47f9-9712-5f207a2444d9
# ╟─42e76b94-d906-4706-a048-89937294b757
# ╠═cf58e15e-025a-4160-b04b-575e821efdee
# ╟─bf8f49a0-c0bf-40cd-8bd0-f9ab8802ba09
# ╠═0c693c24-4ac0-11eb-2329-c743dcc5039d
# ╟─060f0de1-1ced-4f10-9137-0267c3572ebc
# ╟─243d30cb-9072-4eb9-adc0-ff7dcecc2bb5
# ╟─119d17c2-cd47-4827-9907-fe3d9f4e5a65
# ╠═2dc95df5-e4f7-4053-ae79-4125002e925b
# ╠═65ee39ab-3c2b-41fc-bc78-a082eb61089a
# ╠═49d3878a-3d70-48ff-94a5-ee2bac3844b1
# ╠═f038960c-6074-4216-8186-89275c207e3e
# ╠═31b7328c-126b-45ab-a865-3125009ba9b2
# ╠═25dcf9e5-010e-40c9-a0d8-7488e9e94203
# ╠═e4564f72-7bf6-47f2-a8c6-ca41c81146ab
# ╠═e8b365f5-0893-47e6-a213-d14eedba9070
# ╟─1c22b880-4abf-11eb-3f18-756c1198ccad
# ╟─37086212-4abf-11eb-3ec9-7f8dae57121e
# ╠═489421d8-4abf-11eb-0d5e-fd779cc918a1
# ╟─2a5fca7c-4ac0-11eb-33a3-23d972ca27b8
# ╠═1ab81428-8df4-405d-ad60-d1c6e3f6e834
# ╠═3896642a-4ac0-11eb-2c7c-4f376ab82217
# ╠═3ef3faf8-4ac0-11eb-1965-fd23413e29f3
# ╟─3916f50e-661d-11eb-0829-cb3821836fdf
# ╠═4118016e-4ac0-11eb-18bf-5de326782c87
# ╠═4119fbca-4ac0-11eb-1ea9-0bdd324214c5
# ╟─5ebe25c0-661d-11eb-389b-3d81570f7cf0
# ╠═4139bf3c-4ac0-11eb-2b63-77a513149351
# ╟─de48a3f6-4f2f-11eb-314b-493546c37a21
# ╠═85b96ff0-4ac2-11eb-077f-cf4aad8a3c24
# ╟─a11c2898-4ac2-11eb-24d3-6f8060b5fd65
# ╠═a686e67e-4ac2-11eb-228e-23524a3ddc59
# ╟─d50cced2-500d-11eb-2dcc-21fc50825f43
# ╟─ad156892-4ac2-11eb-3634-a3783231e5a1
# ╠═bf1385da-4ac2-11eb-3992-41abac921370
# ╠═e5eeb54f-49f4-4a24-a3f4-6d0b0da99c76
# ╟─d779956a-4ac2-11eb-39de-4b3cecace452
# ╠═c7d2a048-4ac2-11eb-3902-b7c8505096ae
# ╠═cf35b2b2-4ac2-11eb-1ae6-5d3c108210df
# ╠═d1010f88-4ac2-11eb-0fa9-0902fef0cf9f
# ╟─d43d2244-7b1b-4360-9c5e-729f427fe9a7
# ╠═430b3fd0-58b8-443a-9093-24c085af422d
# ╟─29b3b33c-92e9-4579-8f59-36e81887e283
# ╟─5ab03afc-a5fe-4d89-9327-81cc24afd4e0
# ╠═a3969292-57ff-11eb-059b-e9e931a30dc1
# ╟─9e211fbf-4c43-439a-95ba-86b48767c053
# ╟─17ef818e-52e6-400e-b839-fe8dd008aef7
# ╠═dbd89309-f8a6-402f-8703-6f6867a4a19a
# ╟─cf65bee7-9f09-4b23-9898-07da6fdca98c
