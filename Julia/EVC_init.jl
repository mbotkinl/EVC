# Wrapper Code to load scenario and most packages to run EVC code
# Micah Botkin-Levy

println("Loading Packages...")
runParallel=false   #run local EV problems in parallel

using Distributed

if runParallel
	addprocs(3)
	@everywhere using Suppressor
	@everywhere using JuMP
	@everywhere using Printf
	@everywhere using Parameters
	@everywhere using SharedArrays
	@everywhere include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//structEVC.jl")
else
	using Suppressor
	using JuMP
	using Printf
	using Parameters
	using SharedArrays
	include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//structEVC.jl")
end

using LinearAlgebra
using Plots;pyplot()
using Statistics
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVChelpers.jl")

N=100 # number of EVs

# path to working directory
#path="C:\\Users\\micah\\Documents\\uvm\\Research\\Results\\N$(N)_new\\"
path="C:\\Users\\micah\\Documents\\uvm\\Research\\Results\\N$(N)_largeQ\\"
#path="C:\\Users\\micah\\Documents\\uvm\\Research\\Results\\N$(N)\\PWL\\"
#path="C:\\Users\\micah\\Documents\\uvm\\Research\\Results\\N$(N)_K\\"

# Main algorithm parameters
datafile="jld2"    # "jld2" for reading in scenario file, "n" for creating new scenario
updateMethod="dualAscent" #dual decompostion method: "dualAscent" or "fastAscent"
maxIt=500           # max number of iteration per timestep
dualChk = 5e-4      # stopping criteria for ||lam^p-lam^(p-1)||_2
primChk = 5e-4      # stopping criteria for ||coupling constraint||_1
saveLogInd=[1,2,71,141,210] # time step index: used for comparing hot start timesteps between central and decentralized methods
noTlimit=false      # false for coordinator, true for uncooridated (no temperature limit)
forecastError=false # adds forecast to iD (see below, not fully implemented)
relaxedMode=0       # controls convex relaxation (0=NL, 1=QCQP, 2=SOCP) for NL scripts (not fully implemented)
slack=false         # adds slack to temperature constraint (not fully implemented)
eqForm=false        # ALADIN equality vs inequality form
tempAugment=false   # add temperature augmentation objective term (not fully implemented)
ψ=-0                # tempAugment weight for objective function

# Parameters that control I/O
drawFig=false      # draws output figures
saveResults=false  # saves result file
saveS=false        # saves scenario file
loadResults=false  # loads results
silent=true       # prevents output to console
solverSilent=true  # prevents feedback from JuMP solver

file="EVCscenarioN$(N)."*datafile
if datafile=="jld2"  # load scenario
	using FileIO
	println("Reading in Scenario Data...")
	loadF=load(path*file)
	evS=loadF["evScenario"]
else               #create scenario
	println("Creating EV Scenario...")
	Tmax=100            # Temperature limit (Celsius)
	num_homes=1000      # affects background demand magnitude
	Tamb_amplitude = 18 # affects ambeint temperature magnitude
	Dload_error=0       # legacy

	include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCscenario.jl")
    if saveS using FileIO end
	using Distributions
	evS=setupScenario(N;Tmax=Tmax,num_homes=num_homes,Dload_error=Dload_error,saveS=saveS,path=path)
end

# add forcast error to iD_actual
if forecastError
	using Distributions
	L=Int(round(evS.K/2))
	noise = pdf.(Normal(0,1),range(-5,5,length=L))
	noise_norm = (noise.-minimum(noise))/(maximum(noise)-minimum(noise))
	noisePerc=5
	iD_error=zeros(length(evS.iD_pred),1)
	iD_error[1:L]=round.(noise_norm*noisePerc,digits=3)
	iD_error[iD_error.<0.1].=0
	iD_actual=evS.iD_pred .+iD_error
else
	iD_actual=evS.iD_pred
end

# xticks for graphing
stT1=Time(20,0)
endT1=Time(23,59)
stT2=Time(0,0)
endT2=Time(10,0)
Xlabels=vcat(collect(stT1:Dates.Second(round(evS.Ts)):endT1),collect(stT2:Dates.Second(round(evS.Ts)):endT2))
xticks=(1:40:evS.K,Dates.format.(Xlabels[1:40:evS.K],"HH:MM"))
