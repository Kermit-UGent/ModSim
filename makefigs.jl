# A julia script to run the notebooks and save the figures

using Plots

fontscale = 1.5

notebooks = Dict(
    #"modelling_distributions.jl" => "probmod",
    #"MCMC.jl" => "MCMC",
    #"simulation_tools.jl" => "simulation_tools",
    #"modelling_ODEs.jl" => "ODEs",
    #"optimization.jl" => "optimization",
    "calibration.jl" => "calibration",  # LV very unstable, might need to rerun several times
)


for (nb, dir) in notebooks
    # run notebook and generate a plots dict
    let
        include(joinpath("scripts", nb))
        for (name, pl) in plots
            Plots.scalefontsizes(fontscale)
            savefig(pl, joinpath("figures", dir, name * "_sb.pdf"))
            Plots.scalefontsizes(1/fontscale)
            savefig(pl, joinpath("figures", dir, name * ".pdf"))
        end
    end
end

# reset font
