using Pkg; Pkg.activate("./exercises/solved_notebooks")
using PlutoStaticHTML

cd(@__DIR__)

notebookdir = "../exercises/solved_notebooks"
failed_notebooks_dict = Dict{String, Vector{String}}() # keep track of notebooks that failed compilation

subdirs = readdir(notebookdir, join = true) |> x -> filter(isdir, x)
for subdir in subdirs
    bopts = BuildOptions(subdir, write_files = false, use_distributed = false)
    for file in readdir(subdir)
        try # prevent error in one practical from stopping the rest
            build_notebooks(bopts, [file])
        catch e
            @warn e.msg # give error message as a warning instead
            cd(@__DIR__) # errors in notebook building changes working directory sometimes, for some reason

            # log what notebooks failed (per practical)
            global failed_notebooks_dict # github actions likes explicit scoping
            subdirname = split(subdir, ['\\', '/'])[end] # get name of subdir from its path
            failed_notebooks_in_subdir = get!(failed_notebooks_dict, subdirname, String[])

            # get essence of error
            error_lines = split(e.msg, "\n")
            startidx = findfirst(x -> !isnothing(match(r"^Error:", x)), error_lines) + 1
            stopidx = findfirst(x -> !isnothing(match(r"Stacktrace:", x)), error_lines) - 1
            if !(isnothing(startidx) || isnothing(stopidx))
                push!(failed_notebooks_in_subdir, file * *("\n" .* error_lines[startidx:stopidx]...))
            else
                push!(failed_notebooks_in_subdir, file)
            end
        end
    end
end

if !isempty(failed_notebooks_dict)
    error_msg = [
            subdirname * "\n" * *(("\t - " .* failed_notebooks_dict[subdirname] .* "\n")...) 
            for subdirname in sort(collect(keys(failed_notebooks_dict)))
        ] |> 
        x -> reduce(vcat, x) |>
        x -> *(x...)

    error("The following notebook(s) failed to precompile:" * "\n" * error_msg)
end
