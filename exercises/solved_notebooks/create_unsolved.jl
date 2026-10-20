"""
    create_unsolved(file::String; output_path::String = homedir() * "/Downloads")

Make the unsolved version of a notebook.

Cells that contain a line that starts with `#` and contains a "missing" will be changed: 
all uncommented lines are removed and all commented lines will be then uncommented.
"""
function create_unsolved(file::String; output_path::String = homedir() * "/Downloads")
    lines = readlines(file)
    cell_start_idxs = findall(occursin.("# ╔═╡", lines))
    delete_idxs = Int64[]
    for cell_nr in 1:(length(cell_start_idxs)-1)
        cell_idxs = (cell_start_idxs[cell_nr]+1):(cell_start_idxs[cell_nr+1]-1)
        needs_replacement = occursin.(r"^(\s*)\#.*missing", lines[cell_idxs]) |> any
        
        needs_replacement && for cell_idx in cell_idxs
            isempty(lines[cell_idx]) && continue
            if lines[cell_idx][1] == '#'
                println(lines[cell_idx])
                lines[cell_idx] = lines[cell_idx][3:end]
            else
                push!(delete_idxs, cell_idx)
            end
        end
    end
    deleteat!(lines, delete_idxs)
    
    output_filename = replace(splitpath(file)[end], "_sol" => "")
    write(output_path * '/' * output_filename, *(lines .* '\n'...))

    return nothing
end

file = "P7_sens_uncert/sens_intro_sol.jl"
create_unsolved(file)