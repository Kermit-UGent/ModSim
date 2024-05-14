#======================================================================================
=======================================================================================
  Ecosystem Model
=======================================================================================
======================================================================================#

using DifferentialEquations
using Plots

#======================================================================================
  Solving using DifferentialEquations
======================================================================================#

#======================================================================================
  Model function
======================================================================================#

function ecosystem!(du, u, p, t)
    Pc, γ, Lc, σ, pin = p
    c = u[1]
    v = u[2]
    L = c / σ
    λ = L / (L + Lc)
    cdot = (1 - λ) * Pc - γ * c + pin
    vdot = λ * (1 - v) * v * Pc / c - γ *v
    du[1] = cdot
    du[2] = vdot
end

#======================================================================================
  Leaf growth
======================================================================================#

# parameters
Pc = 20
γ = 0.1
Lc = 5
σ = 5

# constant input value
pin = 0

# initial values
c₀ = 0.5
v₀ = 0.1

u0 = [c₀, v₀]
params = Pc, γ, Lc, σ, pin
tspan = (0.0, 150.0)
oprob = ODEProblem(ecosystem!, u0, tspan, params)

sol = solve(oprob, abstol=1e-9, reltol=1e-9)

# convert 92-element Vector{Vector{Float64}} into 92×2 Matrix{Float64}
solmat = mapreduce(permutedims, vcat, sol.u)

t = sol.t
c = solmat[:,1]
v = solmat[:,2]
λ = (c ./ σ) ./ (c ./ σ .+ Lc)

Pc = plot(t, c, linecolor=:black, linewidth=1.5, xaxis="time [s]", yaxis="c [kg C / m²]", label="c")
Pv = plot(t, v, linecolor=:blue, linewidth=1.5, xaxis="time [s]", yaxis="v [m² / m²]", label="v")
Pλ = plot(t, λ , linecolor=:green, linewidth=1.5, xaxis="time [s]", yaxis="λ", label="λ")
plot(Pc, Pv, Pλ, layout=(1, 3))

c[end]
v[end]


#======================================================================================
  Finding steady state
======================================================================================#

u_guess = [c[end], v[end]]
v_oper, c_oper = solve(SteadyStateProblem(ODEProblem(ecosystem!, u_guess, tspan, params)))
# u: 2-element Vector{Float64}:
#  59.30703308172534
#   0.5784648345913733

#======================================================================================
  Effect of initial vegetation area
======================================================================================#

Pcv = plot(xaxis="c", yaxis="v")
v₀_values = collect(0:0.1:1)
for v₀ in v₀_values
    println(v₀)
    oprob = ODEProblem(ecosystem!, [c₀, v₀], tspan, params)
    sol = solve(oprob, abstol=1e-9, reltol=1e-9)
    solmat = mapreduce(permutedims, vcat, sol.u)
    t = sol.t
    c = solmat[:,1]
    v = solmat[:,2]
    plot!(Pcv, c, v, linewidth=1.5, label="v₀ = "*string(v₀))
end
plot(Pcv)


