# A julia script to run the notebooks and save the figures

using Plots

fontscale = 1.5

Plots.scalefontsizes(fontscale)

notebooks = Dict(
     "modelling_distributions.jl" => "probmod",
    "MCMC.jl" => "MCMC",
    "simulation_tools.jl" => "simulation_tools",
)

for (nb, dir) in notebooks
    # run notebook and generate a plots dict
    include(joinpath("scripts", nb))
    for (name, pl) in plots
        savefig(pl, joinpath("figures", dir, name*".pdf"))
    end
end

# reset font
Plots.scalefontsizes(1/fontscale)