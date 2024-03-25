# A julia script to run the notebooks and save the figures

using Plots

notebooks = Dict(
    "probmod" => "modelling_distributions.jl",
    "simulation_tools" => "simulation_tools.jl"
)



for (dir, nb) in notebooks
    # run notebook and generate a plots dict
    include(joinpath("scripts", nb))
    for (name, pl) in plots
        savefig(pl, joinpath("figures", dir, name*".pdf"))
    end
end