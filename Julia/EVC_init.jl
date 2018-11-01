# C:\Users\micah\AppData\Local\JuliaPro-0.6.2.2\Julia-0.6.2\bin\julia
# C:\Users\micah\AppData\Local\Julia-0.7.0\bin\julia
# C:\Users\micah\AppData\Local\Julia-1.0.1\bin\julia

println("Loading Packages...")
using JuMP
using Parameters
using Printf
using Distributed
using LinearAlgebra
using Plots
using Statistics
pyplot()


include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//structEVC.jl")
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVChelpers.jl")

N=10
path="C:\\Users\\micah\\Documents\\uvm\\Research\\Results\\N$(N)\\"
datafile="jld2" #"mat" #"jld" #"n"
file="EVCscenarioN$(N)."*datafile

updateMethod="dualAscent" #dualAscent #fastAscent
maxIt=100
noTlimit=false
relaxed=false
slack=false

drawFig=false
saveResults=false
loadResults=false
#verbose=false

if datafile=="jld2"
	using JLD2
	println("Reading in Data...")
	loadF=load(path*file)
	evS=loadF["evScenario"]
else #create scenario

	Tmax=371
	Dload_amplitude=0
	saveS=true

	include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCscenario.jl")
    if saveS==true using JLD2 end
	using Distributions
	println("Creating EV Scenario...")
	evS=setupScenario(N;Tmax=Tmax,Dload_amplitude=Dload_amplitude,saveS=saveS,path=path)
end


#run comparison
path = clips()
path=path*"\\"
cRun,runs=readRuns(path)
compareRunsGraph(runs, cRun)

# @time runALADit(1)
#@time testALAD(1)

# Profile.clear()
# runALADit(1)
# @profile runALADit(1)
# Juno.profiler()


# Profile.clear()
# testDual(1)
# @profile testDual(1)
# Juno.profiler()
