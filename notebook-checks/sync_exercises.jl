#!/usr/bin/env julia
# Generate src/exercises/*.jl (the website copies) from the student notebooks in
# exercises/student_notebooks/P*/*.jl. The student notebooks are the source of truth;
# NEVER edit src/exercises/*.jl by hand.
#
#   julia notebook-checks/sync_exercises.jl          # (re)write src/exercises/*.jl
#   julia notebook-checks/sync_exercises.jl --check  # exit 1 if src/exercises is stale
#
# stdlib only, no --project needed.
#
# What the generated copy adds on top of the student notebook:
#   1. a `#> [frontmatter]` block (order/title/description/author from the table below)
#      right after the Pluto header, followed by a blank line;
#   2. in the (unique) `Pkg.activate("..")` cell: two advice comment lines and the
#      body normalised to `using Pkg; Pkg.activate("../../pluto-deployment-environment")`.
# Everything else (including the Pluto version stamp on line 2) is copied verbatim.

const ROOT        = normpath(joinpath(@__DIR__, ".."))
const STUDENT_DIR = joinpath(ROOT, "exercises", "student_notebooks")
const SRC_DIR     = joinpath(ROOT, "src", "exercises")

const PDE_ACTIVATE = "using Pkg; Pkg.activate(\"../../pluto-deployment-environment\")"
const ADVICE = [
    "# Running this yourself? Point this at your own environment —",
    "# we advise one shared project in the parent folder: Pkg.activate(\"..\")",
]
const TAGS   = "[\"exercises\"]"
const LAYOUT = "\"layout.jlhtml\""

# filename => (order, title, description, author). `author = nothing` omits the author block.
# Order decides the position in the website sidebar.
const FRONTMATTER = Dict{String,NamedTuple}(
    # P1_mtk
    "ode_model_mtk_intro.jl"                  => (order = "1",  title = "1. ODE_model_MTK_intro",                   description = "Introduction to ModelingToolkit",                                       author = nothing),
    "ode_model_irrigation_mtk.jl"             => (order = "2",  title = "1. ODE_model_irrigation",                  description = "Modeling of irrigation with MTK",                                       author = "Gauthier Vanhaelewyn"),
    "ode_model_diver_mtk.jl"                  => (order = "3",  title = "1. ODE_model_diver",                       description = "modeling of pressure on diver with MTK",                                author = "Gauthier Vanhaelewyn"),
    "ode_model_tank_h_mtk.jl"                 => (order = "4",  title = "1. ODE_model_tank_h_mtk",                  description = "modeling the height of the water in a tank",                            author = nothing),
    "ode_model_tractor_seat_mtk.jl"           => (order = "5",  title = "1. ODE_model_tractor_seat",                description = "modeling the movement of a oscilatory tractor seat",                    author = nothing),
    "ode_model_XTRA_tank_T_h_mtk.jl"          => (order = "6",  title = "1. ODE_model_Xtra_tank",                   description = "extra exercises on modeling water height in a tank",                    author = nothing),
    "ode_model_XTRA_temp_reactors_mtk.jl"     => (order = "7",  title = "1. ODE_model_temp_reactor",                description = "modeling temperature in a CSTR",                                        author = nothing),
    "ode_model_XTRA_water_evap_infil_mtk.jl"  => (order = "8",  title = "1. ODE_model_Xtra_evaporation",            description = "modeling evaporation and infiltration in ground",                       author = nothing),
    # P2_ode
    "ode_model_catalyst_intro.jl"             => (order = "9",  title = "2. ODE_model_Catalyst_intro",              description = "Introduction to Catalyst as an alternative to ModelingToolkit",         author = "Gauthier Vanhaelewyn"),
    "ode_model_birth_death.jl"                => (order = "10", title = "2. ODE_model_birth_death",                 description = "Simple birth-death model for a mice population",                        author = "Gauthier Vanhaelewyn"),
    "ode_model_fermenter_monod.jl"            => (order = "11", title = "2. ODE_model_fermenter_monod",             description = "Fermenter with biomass growing on substrate through Monod kinetics",    author = "Gauthier Vanhaelewyn"),
    "ode_model_infection.jl"                  => (order = "12", title = "2. ODE_model_infection",                   description = "Infection model built as a reaction network",                           author = "Gauthier Vanhaelewyn"),
    "ode_model_XTRA_fermenter_firstorder.jl"  => (order = "13", title = "2. ODE_model_Xtra_fermenter_firstorder",   description = "Extra exercise on a fermenter with first-order kinetics",               author = "Gauthier Vanhaelewyn"),
    "ode_model_XTRA_anaerobic_fermentation.jl"=> (order = "14", title = "2. ODE_model_Xtra_anaerobic_fermentation", description = "Extra exercise on anaerobic fermentation",                              author = "Gauthier Vanhaelewyn"),
    "ode_model_XTRA_soil_cont_plant_uptake.jl"=> (order = "15", title = "2. ODE_model_Xtra_soil_contamination",     description = "Extra exercise on soil contamination with plant uptake",                author = "Gauthier Vanhaelewyn"),
    "ode_model_XTRA_water_evap_infil.jl"      => (order = "16", title = "2. ODE_model_Xtra_evaporation",            description = "Extra exercise on water evaporation and infiltration in a reservoir",   author = "Gauthier Vanhaelewyn"),
    # P3_sde
    "sde_model_heston_mtk_intro.jl"           => (order = "17", title = "3. SDE_model_Heston_intro",                description = "Introduction to solving SDE problems with ModelingToolkit",             author = "Gauthier Vanhaelewyn"),
    "sde_model_aging_mtk.jl"                  => (order = "18", title = "3. SDE_model_aging",                       description = "Aging with saturated repair, modelled as an SDE",                       author = "Gauthier Vanhaelewyn"),
    "sde_model_fermenter_secondorder_mtk.jl"  => (order = "19", title = "3. SDE_model_fermenter_secondorder",       description = "Fermenter with second-order kinetics, modelled as an SDE",              author = "Gauthier Vanhaelewyn"),
    "dje_model_catalyst_intro.jl"             => (order = "20", title = "3. DJE_model_Catalyst_intro",              description = "Introduction to solving discrete jump problems with Catalyst",          author = "Gauthier Vanhaelewyn"),
    "dje_model_bike_sharing.jl"               => (order = "21", title = "3. DJE_model_bike_sharing",                description = "Discrete jump model of a simple bike sharing system",                   author = "Gauthier Vanhaelewyn"),
    "dje_model_festival_toilet.jl"            => (order = "22", title = "3. DJE_model_festival_toilet",             description = "Discrete jump model of a festival toilet queue",                        author = "Gauthier Vanhaelewyn"),
    # P4_probmod
    "probmod_1-intro.jl"                      => (order = "23", title = "4. ProbMod intro",                         description = "Introduction to the sampling practicals",                               author = "Bram Spanoghe"),
    "probmod_2-basics.jl"                     => (order = "24", title = "4. ProbMod basics",                        description = "Basic sampling exercises",                                              author = "Bram Spanoghe"),
    "probmod_3-advanced.jl"                   => (order = "25", title = "4. ProbMod advanced",                      description = "Advanced sampling exercises",                                           author = "Bram Spanoghe"),
    "probmod_4-review.jl"                     => (order = "26", title = "4. ProbMod review",                        description = "Review sampling exercise",                                              author = "Bram Spanoghe"),
    # P5_mcmc
    "MCMC_1-intro.jl"                         => (order = "27", title = "5. MCMC intro",                            description = "MCMC intro",                                                            author = "Bram Spanoghe"),
    "MCMC_2-basics.jl"                        => (order = "28", title = "5. MCMC basics",                           description = "MCMC basics",                                                           author = "Bram Spanoghe"),
    "MCMC_3-advanced.jl"                      => (order = "29", title = "5. MCMC advanced",                         description = "MCMC advanced",                                                         author = "Bram Spanoghe"),
    "MCMC_4-review.jl"                        => (order = "30", title = "5. MCMC review",                           description = "MCMC review",                                                           author = "Bram Spanoghe"),
    # P6_calib
    "calib_intro.jl"                          => (order = "31", title = "6. Calibration intro",                     description = "Calibration intro",                                                     author = "Gauthier Vanhaelewyn"),
    "calib_fermenter_monod.jl"                => (order = "32", title = "6. Calibration fermenter monod",           description = "Calibration fermenter monod",                                           author = "Gauthier Vanhaelewyn"),
    "calib_irrigation.jl"                     => (order = "33", title = "6. Calibration irrigation",                description = "Calibration irrigation",                                                author = "Gauthier Vanhaelewyn"),
    "optim_wastewater_treatment.jl"           => (order = "34", title = "6. Optimisation wastewater treatment",     description = "Optimisation wastewater treatment",                                     author = "Gauthier Vanhaelewyn"),
    # P7_sens_uncert
    "sens_intro.jl"                           => (order = "35", title = "7. Sensitivity intro",                     description = "Sensitivity intro",                                                     author = "Gauthier Vanhaelewyn"),
    "sens_fermenter_monod.jl"                 => (order = "36", title = "7. Sensitivity fermenter monod",           description = "Sensitivity fermenter monod",                                           author = "Gauthier Vanhaelewyn"),
    "sens_bitrophic_model.jl"                 => (order = "37", title = "7. Sensitivity bitrophic model",           description = "Sensitivity bitrophic model",                                           author = "Gauthier Vanhaelewyn"),
    "sens_insuline.jl"                        => (order = "38", title = "7. Sensitivity insuline",                  description = "Sensitivity insuline",                                                  author = "Gauthier Vanhaelewyn"),
    "uncert_intro.jl"                         => (order = "39", title = "7. Uncertainty intro",                     description = "Uncertainty intro",                                                     author = "Gauthier Vanhaelewyn"),
    "uncert_fermenter_monod.jl"               => (order = "40", title = "7. Uncertainty fermenter monod",           description = "Uncertainty fermenter monod",                                           author = "Gauthier Vanhaelewyn"),
    "uncert_bitrophic_model.jl"               => (order = "41", title = "7. Uncertainty bitrophic model",           description = "Uncertainty bitrophic model",                                           author = "Gauthier Vanhaelewyn"),
    # P8_modselect
    "model_selection_intro.jl"                => (order = "42", title = "8. Model selection intro",                 description = "Model selection intro",                                                 author = nothing),
    "probabilistic_selection.jl"              => (order = "43", title = "8. Probability selection",                 description = "Probability selection",                                                 author = nothing),
)

frontmatter_lines(meta) = begin
    lines = [
        "#> [frontmatter]",
        "#> order = \"$(meta.order)\"",
        "#> title = \"$(meta.title)\"",
        "#> tags = $TAGS",
        "#> layout = $LAYOUT",
        "#> description = \"$(meta.description)\"",
    ]
    if meta.author !== nothing
        append!(lines, ["#> ", "#>     [[frontmatter.author]]", "#>     name = \"$(meta.author)\""])
    end
    lines
end

"""
    render(student_path, meta) -> String

Build the website copy of a student notebook: insert the frontmatter block and rewrite
the `Pkg.activate("..")` cell. Throws if the student notebook does not have the shape
we expect (so a broken notebook fails loudly instead of producing garbage).
"""
function render(student_path::AbstractString, meta)
    src = read(student_path, String)
    occursin('\r', src) && error("$student_path: CRLF line endings are not supported")
    lines = split(src, '\n')
    length(lines) ≥ 4 || error("$student_path: file too short to be a Pluto notebook")
    lines[1] == "### A Pluto.jl notebook ###" || error("$student_path: line 1 is not the Pluto header")
    startswith(lines[2], "# v") || error("$student_path: line 2 is not the Pluto version stamp")
    lines[3] == "" || error("$student_path: expected an empty line 3")
    any(startswith(l, "#> ") for l in lines) && error("$student_path: student notebook already has a `#>` frontmatter block")

    # Locate the unique activate cell.
    act = findall(l -> occursin("Pkg.activate(\"..\")", l), lines)
    length(act) == 1 || error("$student_path: expected exactly one `Pkg.activate(\"..\")`, found $(length(act))")
    i = act[1]
    header = findlast(k -> startswith(lines[k], "# ╔═╡ "), 1:i)
    header === nothing && error("$student_path: `Pkg.activate` is not inside a Pluto cell")
    body_start = header + 1
    while body_start ≤ length(lines) && startswith(lines[body_start], "# ╠═╡ ")
        body_start += 1  # keep cell metadata lines (e.g. `# ╠═╡ show_logs = false`)
    end
    body_end = findnext(isempty, lines, body_start)  # cell body runs until the next empty line
    body_end === nothing && error("$student_path: activate cell has no terminating empty line")
    body_end > i || error("$student_path: `Pkg.activate` found outside its cell body")

    out = String[]
    append!(out, lines[1:3])
    append!(out, frontmatter_lines(meta))
    push!(out, "")
    append!(out, lines[4:body_start-1])
    append!(out, ADVICE)
    push!(out, PDE_ACTIVATE)
    append!(out, lines[body_end:end])
    join(out, '\n')
end

student_files() = begin
    files = Pair{String,String}[]  # basename => path
    for p in sort(readdir(STUDENT_DIR))
        (startswith(p, "P") && isdir(joinpath(STUDENT_DIR, p))) || continue
        for f in sort(readdir(joinpath(STUDENT_DIR, p)))
            (endswith(f, ".jl") && !endswith(f, "_sol.jl")) || continue
            push!(files, f => joinpath(STUDENT_DIR, p, f))
        end
    end
    files
end

function main(args)
    check = "--check" in args
    fail = false
    files = student_files()
    produced = Set{String}()

    seen = Set{String}()
    for (name, _) in files
        name in seen && (println("ERROR   $name appears in more than one P*/ folder"); fail = true)
        push!(seen, name)
    end
    for k in sort(collect(keys(FRONTMATTER)))
        k in seen || (println("ERROR   table entry $k has no student notebook"); fail = true)
    end

    for (name, path) in files
        rel = relpath(path, ROOT)
        haskey(FRONTMATTER, name) || (println("ERROR   $rel is not in the FRONTMATTER table of $(relpath(@__FILE__, ROOT))"); fail = true; continue)
        rendered = try
            render(path, FRONTMATTER[name])
        catch e
            println("ERROR   $rel: ", sprint(showerror, e)); fail = true; continue
        end
        push!(produced, name)
        target = joinpath(SRC_DIR, name)
        current = isfile(target) ? read(target, String) : nothing
        if current == rendered
            println("OK      src/exercises/$name")
        elseif check
            println("DIFF    src/exercises/$name is out of date (run: julia notebook-checks/sync_exercises.jl)"); fail = true
        else
            write(target, rendered)
            println(current === nothing ? "NEW     " : "WROTE   ", "src/exercises/$name")
        end
    end

    for f in sort(readdir(SRC_DIR))
        (endswith(f, ".jl") && !(f in produced)) || continue
        println("STALE   src/exercises/$f has no student notebook (remove it by hand)"); fail = true
    end

    if fail
        println(check ? "\nsrc/exercises is NOT in sync with exercises/student_notebooks." : "\nsync finished with errors.")
        exit(1)
    end
    println(check ? "\nsrc/exercises is in sync." : "\nsync done.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
