# README for Julia EVC

Masters project for electric vehicle charging


## General Code/Folder Structure

## Running Code

### Prerequisites

Main Julia Packages needed
```
JuMP
Gurobi
```

### Initialization 

Run `EVC_init.jl` to load packages, set parameters, and load or create scenario. Setting `datafile="n"` will create new scenario and `datafile="jld2"` will read saved scenario file. 

### Create/Read Central Solution

Run `EVCcentral.jl` will read or run central solution depending on `loadResults` parameter. 

### Create/Read Decentralized Solution

Similarly running scripts `EVCdecentral.jl`, `EVCdecentral_ADMM.jl`, or `EVCdecentral_ALADIN.jl` will run dual decomposition, ADMM, and ALADIN distributed methods, respectively. 
