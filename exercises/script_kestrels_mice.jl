
using Catalyst
using DifferentialEquations, Plots

# valken_konijnen = @reaction_network begin
#     @species V(t)=3 K(t)=12
#     @parameters rv=0.6 Kv=10 dv=0.06 rk=4.2 Kk=60 dk=1.9 p=5.0e-2
#     rv*(1 - V/Kv)*V, 0 --> V
#     dv, V --> 0
#     rk*(1 - K/Kk)*K, 0 --> K
#     dk, K --> 0
#     p, K + V --> V
# end

valken_konijnen = @reaction_network begin
    @species V(t)=3 K(t)=12
    @parameters rv=0.6 Kv=10.0 dv=0.06 rk=2.2 Kk=60.0 dk=1.9 p=0.05
    rv*(1 - V/Kv), V + K --> 2V   # 2V, must be an integer for JumpProblem !!!! so, not (1+0.01)*V
    dv, V --> 0
    rk*(1 - K/Kk), K --> 2*K
    dk, K --> 0
end


species(valken_konijnen)
parameters(valken_konijnen)
reactions(valken_konijnen)

osys  = convert(ODESystem, valken_konijnen)   # convert rn to a ODESystem
equations(osys)

u0 = [:V => 3, :K => 12]
tspan = (0.0, 10.0)
# params = [:rv=>0.6, :Kv=>10, :dv=>0.06, :rk=>4.2, :Kk=>60, :dk=>1.9, :p=>5.0e-2]
params = [:rv=>0.6, :Kv=>15, :dv=>0.06, :rk=>4.2, :Kk=>60, :dk=>0.9, :p=>0.10]

oprob = ODEProblem(valken_konijnen, u0, tspan, params)
osol = solve(oprob, Tsit5(), saveat=0.1)

plot(osol)


dprob = DiscreteProblem(valken_konijnen, u0, tspan, params)
jdprob = JumpProblem(valken_konijnen, dprob, Direct())
jdsol = solve(jdprob, SSAStepper())
plot(jdsol)
