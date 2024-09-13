using ModelingToolkit
# using ModelingToolkit: t_nounits as t, D_nounits as D
using DifferentialEquations
# using Symbolics
using Plots

# @mtkmodel FOL begin
#     @parameters begin
#         τ # parameters
#     end
#     @variables begin
#         x(t)   # dependent variables
#     end
#     @equations begin
#         D(x) ~ (1 - x) / τ
#     end
# end

# using DifferentialEquations: solve

# @mtkbuild fol = FOL()
# equations(fol)
# states(fol)
# parameters(fol)

# prob = ODEProblem(fol, [fol.x => 0.0], (0.0, 10.0), [fol.τ => 3.0])
# sol = solve(prob)

# plot(sol)

# # matrices, simplified_sys = linearize(equations(fol), [fol.x], [fol.y])





@variables t x(t)=0 y(t)=0 u(t)=0 r(t)=0
D = Differential(t)

@parameters kp = 1

eqs = [u ~ kp * (r - y) # P controller
       D(x) ~ -x + u    # First-order plant
       y ~ x]           # Output equation

@named sys = ODESystem(eqs, t)
op = Dict(x => 1, r => 2)
matrices, simplified_sys = linearize(sys, [r], [y], op=op) # Linearize from r to y
matrices

using ModelingToolkit: inputs, outputs, unknown_states
[unknown_states(simplified_sys); inputs(simplified_sys); outputs(simplified_sys)]

[unknown_states(sys); inputs(sys); outputs(sys)]


# ###############################################################################

# Step function, go from v0 to v1 at t = t0
function heaviside(t, t0, v0, v1)
    (v1 - v0) .* (t .> t0) .+ v0
end



@variables t x(t)=0 y(t)=0 u(t)=0 r(t)=2
D = Differential(t)

@parameters K=1.5 τ=300

eqs = [u ~ r
       D(x) ~ K/τ*u - 1/τ*x
       y ~ x]           # Output equation

@named sys = ODESystem(eqs, t)
equations(sys)
states(sys)
op = Dict(x => 0)
matrices, simplified_sys = linearize(sys, [u], [y], op=op)
matrices

equations(simplified_sys)
states(simplified_sys)
outputs(simplified_sys)

u0 = [:x => 2.5]
tspan = (0, 3600)
oprob = ODEProblem(simplified_sys, u0, tspan)
sol = solve(oprob)

plot(sol)