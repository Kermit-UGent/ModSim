# Exercises

- `student_notebooks/P*/` is the **source of truth**: edit the exercise notebooks here.
- `solved_notebooks/P*/` holds the solutions; CI builds every one of them (`.github/workflows/CheckNotebooks.yml`).
- `src/exercises/*.jl` (the website copies) are **generated** from `student_notebooks/` by `julia notebook-checks/sync_exercises.jl` — never edit them by hand; rerun the script after changing a student notebook, CI (`--check`) fails otherwise.
- To work on a notebook, start Pluto with the top-level `launch_pluto.jl` (`julia launch_pluto.jl` in the repo root).
