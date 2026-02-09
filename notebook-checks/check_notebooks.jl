using Pkg; Pkg.activate("./exercises/solved_notebooks")
using PlutoStaticHTML

cd(@__DIR__)

notebookdir = "../exercises/solved_notebooks"
any_notebook_failed = false # check if any notebook has failed compilation

subdirs = readdir(notebookdir, join = true) |> x -> filter(isdir, x)
for subdir in subdirs
    bopts = BuildOptions(subdir, write_files = false, use_distributed = false)
    try # prevent error in one practical from stopping the rest
        build_notebooks(bopts)
    catch e
        global any_notebook_failed = true
        @warn e.msg
    end
end

any_notebook_failed && error("One or more notebooks failed to compile.")