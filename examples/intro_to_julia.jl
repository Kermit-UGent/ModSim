### A Pluto.jl notebook ###
# v0.19.43

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

# ╔═╡ 9a1cca86-0fff-4f8a-a033-8cf6db337503
using PlutoUI, Markdown; TableOfContents()

# ╔═╡ a15c46f7-b561-41c7-8da2-12f04465fb19
using LinearAlgebra  

# ╔═╡ 8f785f5f-163c-436e-8781-0b68c3aa422a
using Statistics

# ╔═╡ bf1385da-4ac2-11eb-3992-41abac921370
using Plots

# ╔═╡ e97e5984-4ab9-11eb-3efb-9f54c6c307dd
# edit the code below to set your name and UGent email

student = (name = "Jenke Janssen", email = "Jenke.Janssen@UGent.be");

# press the ▶ button in the bottom right of this cell to run your edits
# or use Shift+Enter

# you might need to wait until all other cells in this notebook have completed running. 
# scroll down the page to see what's up

# ╔═╡ f089cbaa-4ab9-11eb-09d1-05f49911487f
md"""
Submission by: **_$(student.name)_**
"""

# ╔═╡ fd21a9fa-4ab9-11eb-05e9-0d0963826b9f
md"""
# Notebook 1: Getting up and running

First of all, welcome to the course, **$(student[:name])**! We hope you enjoy the ride.
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
md"You can only have one statement per line and all your variables need to have unaiue names. You can split statements in two lines or wrap them in a `begin ... end` block."

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
md"Similarly to Matlab, when using the REPL, Julia will print the result of every statement by default. To suppress this behavious, just end the statement with a semicolon."

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
md"Indexing is via **square brackets** (like python!) and the **index stqrts from 1** (like in Matlab and R)."

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
md"Many of the functions we will need will be fairly simple equations. We can just define them in one line. A more condensed version of of `square(x)`."

# ╔═╡ 463689b0-4ac1-11eb-1b0f-b7a239011c5c
s(x) = x * x

# ╔═╡ 53dd5e65-d4f4-4f56-9ebd-7bad62ebd716
s(8)

# ╔═╡ 0dbe0c34-500e-11eb-2633-67d8dc6b24c8
md"""
Functions are first-class and work just like any other variable! For example, you can give a function as an input in another function. 

In some cases, you might want to define an **anonymous function**, without given them a name:"""

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
md"The `if`, `else`, `elseif`-statement is instrumental to any programming language./ Note that control flow is ended with an `end` statement. In constrast to Python, tabs are only for clarity but do not impact functionality."

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

Looping using a for loop can be done by iterating over a list or range. Don't forget to end with an ènd`at the end.
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
pets = ["Hedwig", "Pig", "Crookhanks"]

# ╔═╡ 5ebe25c0-661d-11eb-389b-3d81570f7cf0
md"`zip` binds two or more iterators and yield tuples of the pairs."

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

Quite essential for scientific programming is the visualisation of the results. `Plots` is the Julia package that handles a lot of the visualisation. `rand(10)`, returns an array of 10 random floats between 0 and 1.
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

# ╔═╡ 094f35e3-6326-4f3c-853b-b8d08c96d483
hint(text) = Markdown.MD(Markdown.Admonition("hint", "Hint", [text]));

# ╔═╡ 060f0de1-1ced-4f10-9137-0267c3572ebc
hint(md"Check out `min`and `max`.")

# ╔═╡ e13aaaab-bc7f-48c7-98ed-be0d0180e151
almost(text) = Markdown.MD(Markdown.Admonition("warning", "Almost there!", [text]));

# ╔═╡ c70d8364-8a34-42ef-b0c2-96ec46242063
keep_working(text=md"The answer is not quite right.") = Markdown.MD(Markdown.Admonition("danger", "Keep working on it!", [text]));

# ╔═╡ e05273ae-d35b-44e5-b41a-f662410af88c
correct(text=md"Great! You got the right answer! Let's move on to the next section.") = Markdown.MD(Markdown.Admonition("correct", "Got it!", [text]));

# ╔═╡ 243d30cb-9072-4eb9-adc0-ff7dcecc2bb5
if !ismissing(clip(0.1)) 
	if clip(-1) == 0 && clip(0.25) ≈ 0.25 && clip(3.6) ≈ 1
		correct()
	else
		keep_working()
	end
end
    

# ╔═╡ b6b068c1-8f45-4f5d-b11d-8021d581708a
sol_stirling(n) = √(2π * n) * (n/exp(1))^n;

# ╔═╡ 9e211fbf-4c43-439a-95ba-86b48767c053
if !ismissing(stirling(5)) 
	if sol_stirling(20) ≈ stirling(20)
		correct()
	else
		keep_working()
	end
end
    

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Markdown = "d6f4376e-aef5-505a-96c1-9c027394607a"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
Plots = "~1.38.0"
PlutoUI = "~0.7.55"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "c278dfab760520b8bb7e9511b968bf4ba38b7acc"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.3"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BitFlags]]
git-tree-sha1 = "2dc09997850d68179b69dafb58ae806167a32b1b"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.8"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "59939d8a997469ee05c4b4944560a820f9ba0d73"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.4"

[[ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"

    [ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "75bd5b6fc5089df449b5d35fa501c846c9b6549b"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.12.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "8cfa272e8bdedfa88b6aefbbca7c19f1befac519"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.3.0"

[[Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "ac67408d9ddf207de5cfa9a97e114352430f01ed"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.16"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "ff38ba61beff76b8f4acad8ab0c97ef73bb670cb"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.9+0"

[[GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "27442171f28c952804dede8ff72828a96f2bfc1f"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.10"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "025d171a2847f616becc0f84c8dc62fe18f0f6dd"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.10+0"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "e94c92c7bf4819685eb80186d51c43e71d4afa17"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.76.5+0"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "abbbb9ec3afd783a7cbd82ef01dcd088ea051398"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.1"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "8b72179abc660bfab5e28472e019392b97d0985c"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.4"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "a53ebe394b71470c7f97c2e7e170d51df21b17af"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.7"

[[JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "60b1194df0a3298f460063de985eae7b01bc011a"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.1+0"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "f428ae552340899a935973270b8d98e5a31c49fe"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.1"

    [Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"

    [LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "51901a49222b09e3743c65b8847687ae5fc78eb2"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.1"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "60e3045590bd104a16fefb12836c00c0ef8c7f8c"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.13+0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "862942baf5663da528f66d24996eb6da85218e76"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.0"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "9f8675a55b37a70aa23177ec110f6e3f4dd68466"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.38.17"

    [Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "68723afdb616445c6caaef6255067a8339f91325"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.55"

[[PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "37b7bb7aabf9a085e0044307e1717436117f2b3b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.5.3+1"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "1d77abd07f617c4868c33d4f5b9e1dbb2643c9cf"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.2"

[[SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TranscodingStreams]]
git-tree-sha1 = "54194d92959d8ebaa8e26227dbe3cdefcdcd594f"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.3"
weakdeps = ["Random", "Test"]

    [TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "3c793be6df9dd77a0cf49d80984ef9ff996948fa"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.19.0"

    [Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "801cbe47eae69adc50f36c3caec4758d2650741b"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.2+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "522b8414d40c4cbbab8dee346ac3a09f9768f25d"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.5+0"

[[Xorg_libICE_jll]]
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "e5becd4411063bdcac16be8b66fc2f9f6f1e8fe5"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.0.10+1"

[[Xorg_libSM_jll]]
deps = ["Libdl", "Pkg", "Xorg_libICE_jll"]
git-tree-sha1 = "4a9d9e4c180e1e8119b5ffc224a7b59d3a7f7e18"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.3+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a68c9655fbe6dfcab3d972808f1aafec151ce3f8"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.43.0+0"

[[gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "93284c28274d9e75218a416c65ec49d0e0fcdf3d"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.40+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╟─f089cbaa-4ab9-11eb-09d1-05f49911487f
# ╠═9a1cca86-0fff-4f8a-a033-8cf6db337503
# ╠═e97e5984-4ab9-11eb-3efb-9f54c6c307dd
# ╟─fd21a9fa-4ab9-11eb-05e9-0d0963826b9f
# ╟─2f210d77-c1fe-4854-b8cb-2c33dcf64292
# ╠═8dc1e292-43a1-459e-b69d-9ef620e121a1
# ╠═e4f67601-03c5-4920-a1fa-de4ec5947868
# ╠═fe9353dd-f2a5-4597-962d-df267c3c70b8
# ╟─f7c7d372-53eb-402e-a0a9-f7a9a5bdb4de
# ╟─fdc8ef8f-6a52-4374-a343-4f44a419639f
# ╟─27e432a0-cf43-4ab4-b51d-bed405e9f791
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
# ╠═2a5fca7c-4ac0-11eb-33a3-23d972ca27b8
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
# ╠═094f35e3-6326-4f3c-853b-b8d08c96d483
# ╠═e13aaaab-bc7f-48c7-98ed-be0d0180e151
# ╠═c70d8364-8a34-42ef-b0c2-96ec46242063
# ╠═e05273ae-d35b-44e5-b41a-f662410af88c
# ╠═b6b068c1-8f45-4f5d-b11d-8021d581708a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
