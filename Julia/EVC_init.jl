# C:\Users\micah\AppData\Local\JuliaPro-0.6.2.2\Julia-0.6.2\bin\julia
# C:\Users\micah\AppData\Local\Julia-0.7.0\bin\julia
# C:\Users\micah\AppData\Local\Julia-1.0.1\bin\julia
# C:\Users\micah\AppData\Local\Julia-1.0.2\bin\julia


println("Loading Packages...")
using JuMP
using Parameters
using SharedArrays
using Printf
using Distributed
using LinearAlgebra
using Plots
using Statistics
pyplot()

include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//structEVC.jl")
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVChelpers.jl")

N=80
path="C:\\Users\\micah\\Documents\\uvm\\Research\\Results\\N$(N)\\"
datafile="jld2" #"mat" #"jld" #"n"
file="EVCscenarioN$(N)."*datafile

updateMethod="dualAscent" #dualAscent #fastAscent
maxIt=20
noTlimit=false
forecastError=false
relaxedMode=2
slack=false
eqForm=true

drawFig=false
saveResults=false
saveS=false
loadResults=false
silent=true

if datafile=="jld2"
	using FileIO
	println("Reading in Scenario Data...")
	loadF=load(path*file)
	evS=loadF["evScenario"]

	# if isassigned(evS.iDnoise)==false
	# 	evS.iDnoise=zeros(length(iD),1)
	# end

else #create scenario
	println("Creating EV Scenario...")
	Tmax=400/1000
	Dload_amplitude=20
	Dload_error=0

	include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCscenario.jl")
    if saveS using FileIO end
	using Distributions
	evS=setupScenario(N;Tmax=Tmax,Dload_amplitude=Dload_amplitude,Dload_error=Dload_error,saveS=saveS,path=path)
end




#run comparison
#path = clips()
# path=path*"Relax\\"
# cRun,runs, noLim=readRuns(path);
# lowRes=true
# resPlot, convPlot=compareRunsGraph(runs, cRun, noLim, saveResults,lowRes)
# cTable=compareRunsTable(runs)

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
