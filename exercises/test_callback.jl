using Catalyst
using DifferentialEquations, Plots

rn = @reaction_network begin
    (k,1), X1 <--> X2
end
u0 = [:X1 => 10.0,:X2 => 0.0]
tspan = (0.0, 20.0)
p = [:k => 1.0]
oprob = ODEProblem(rn, u0, tspan, p)

parameters(rn)

condition = [5.0]
affect!(integrator) = integrator[:X1] += 5.0
affect!(integrator) = integrator.p[1] = 5.0
# affect!(integrator) = setp(oprob, :k)(integrator, 5.0)
ps_cb = PresetTimeCallback(condition, affect!)

sol = solve(deepcopy(oprob); callback = ps_cb)
plot(sol)