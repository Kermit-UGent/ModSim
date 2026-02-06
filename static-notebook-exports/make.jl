# use julia 1.10!
using Pkg; Pkg.activate("./exercises/solved_notebooks")
using PlutoStaticHTML
using Documenter

cd(@__DIR__)


# Export notebooks to markdown using PlutoStaticHTML.jl

notebookdir = "../exercises/solved_notebooks"
outputdir = "./src"
write("log.txt", "") # clean log

subdirs = readdir(notebookdir, join = true) |> x -> filter(isdir, x)
for subdir in subdirs
    # generate documenter-style markdown files from all pluto notebooks
    bopts = BuildOptions(subdir, output_format = documenter_output)
    oopts = OutputOptions(
        show_output_above_code = true, # show output above code rather than below, to mirror OG Pluto.jl
        convert_admonitions = false # dont change html style of admonitions like `!!! note Hello I am note`, looks weird
    )
    try # prevent error in one practical from stopping the rest
        build_notebooks(bopts)
    catch e
        # log error
        open("log.txt", "a") do f
            write(f, e.msg, "\n\n\n")
        end
    end

    # move the markdown files to their designated folder (they are generated in the folder of the pluto notebooks)
    md_files = readdir(subdir) |> x -> filter(name -> name[end-2:end] == ".md", x) # get all .md files
    subdirname = split(subdir, ['\\', '/'])[end] # get name of subdir from its path
    for md_file in md_files
        mv(subdir * "/" * md_file, outputdir * "/" * subdirname * "/" * md_file, force = true)
    end
end


# Generate HTML from markdown files using Documenter.jl

function get_title(filepath::String)
    lines = readlines(filepath)
    matches = [match(r"(?<=<h1).*?(?=</h1>)", line) for line in lines] |>
        x -> filter(!isnothing, x)
    firstmatch = first(matches).match
    title = split(firstmatch, '>')[2] |> string
    return title
end

subdirnames = split.(subdirs, ['\\', '/']) .|> last .|> string

exercise_pages = [
    "Practicum $i" => [
        get_title("./src/" * subdirname * "/" * file) => subdirname * "/" * file
        for file in readdir("./src/" * subdirname)
    ]
    for (i, subdirname) in enumerate(subdirnames)
]

pages = ["Introduction" => "index.md"; exercise_pages]
write("./src/index.md", "# ModSim exercises\n\nProof of concept for static notebook generation.")

makedocs(;
    sitename = "ModSim",
    pages,
    format = Documenter.HTML(
        size_threshold = 2000 * 1024, # Allow generation of large notebooks
        mathengine = MathJax3() # for LaTeX equations in notebooks
    ),
    warnonly = true # Only warn, don't error if notebooks are large
)


# Throw error at end if any notebook failed to compile

if !isempty(readlines("log.txt"))
    error("Notebook compilation failed for one or more practicals.")
end

using LiveServer
serve(dir = "./build", launch_browser = true)