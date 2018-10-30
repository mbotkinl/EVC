# C:\Users\micah\AppData\Local\JuliaPro-0.6.2.2\Julia-0.6.2\bin\julia
# C:\Users\micah\AppData\Local\Julia-0.7.0\bin\julia

println("Loading Packages...")
using JuMP
using Parameters
using Printf


include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//structEVC.jl")
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVChelpers.jl")

N=10
path="C:\\Users\\micah\\Documents\\uvm\\Research\\Results\\N$(N)\\"
datafile="n" #"mat" #"jld" #"n"
file="EVCscenarioN$(N)."*datafile

updateMethod="dualAscent" #dualAscent #fastAscent
drawFig=0
saveResults=0
loadResults=0
noTlimit=false
maxIt=100
relaxed=false
slack=false
verbose=false

if datafile=="jld"
	using JLD2
	println("Reading in Data...")
	loadF=load(path*file)
	evS=loadF["evScenario"]
else #create scenario

	Tmax=393
	Dload_amplitude=10
	saveS=false

	include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCscenario.jl")
    if saveS==true using JLD2 end
	using Distributions
	println("Creating EV Scenario...")
	evS=setupScenario(N;Tmax=Tmax,Dload_amplitude=Dload_amplitude,saveS=saveS)
end

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
