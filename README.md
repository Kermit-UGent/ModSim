# Modeling and Simulation

Course Modelling and Simulation for Bioscience Engineers at UGent. 

🚧 Work in progress till 2025. 

## launching Pluto

To launch pluto, either:
- open a terminal in the top level folder, and run `julia launch_pluto.jl`
- open a julia terminal in the top level folder, and run `include("launch_pluto.jl")`

## gh-pages / website

all the info and files are in the `src` directory. There will be some redundancy between files/scripts in other folders and the ones in the `src` folder. 

see `website_maintenance.md`

Also don't forget to update the toml and manifest files in `pluto-deployment-environment`. Make sure to respect the Julia version, currently this is `1.11.2`.

If you need a more in-depth example of the different pages, usage of tags, etc. Check out release [v2425.1](https://github.com/Kermit-UGent/ModSim/releases/tag/v2425.1) en run the server on that code, alternatively check out the [original github source of the template](https://github.com/JuliaPluto/computational-thinking-template) or (the current website accompanying that repo)[https://juliapluto.github.io/computational-thinking-template/].