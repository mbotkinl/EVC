# README for Julia EVC

Masters project for electric vehicle charging


## General File/Folder Structure

- Scripts are labeled by the problem (EVC or Hub) and formulation method (central, ADMM, ALAD, etc)
- Wrapper scripts are in the Julia directory
- Wrappers call functions from the files in the Julia/functions directory
- Other scripts used to create supplimental analysis plots are in "Analysis Scripts" directory

### Main Function Structure

The main functions that run the EVC PWL problem are in Julia/functions/funEVCpwl.jl 

Example description of ALADIN functions: 

- `EVCdecentral_ALADIN.jl` wrapper script calls `pwlEValad` wrapper function
- `pwlEValad` loops over MPC timesteps and calls `runEVALADStep` function
- `runEVALADStep` loops over iterations and calls `runEVALADIt` function
- `runEVALADIt` is the wrapper for one iteration of the ALADIN algorithm and calls `localEVALAD`, `localXFRMALAD`, and `coordALAD` functions

## Running Code

### Prerequisites

Main Julia Packages needed
```
JuMP
Gurobi
```

### Initialization 

Run `EVC_init.jl` to load packages, set parameters, and load or create scenario. Setting `datafile="n"` will create new scenario and `datafile="jld2"` will read saved scenario file from `path` and `file`. 

### Create/Read Central Solution

Run `EVCcentral.jl` will read or run central solution depending on `loadResults` parameter. 

### Create/Read Decentralized Solution

Similarly running scripts `EVCdecentral.jl`, `EVCdecentral_ADMM.jl`, or `EVCdecentral_ALADIN.jl` will run dual decomposition, ADMM, and ALADIN distributed methods, respectively. 


## Plotting
Plots use the `Plots` and `PyPlot` packages. Main plots exist at the end of each EVC or Hub wrapper script in the main Julia directory. Other visualization can be created using the scripts in the "Analysis Scripts" or `funEVCvis.jl` but are not fully built out and are still in development.
