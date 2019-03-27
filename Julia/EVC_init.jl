# C:\Users\micah\AppData\Local\JuliaPro-0.6.2.2\Julia-0.6.2\bin\julia
# C:\Users\micah\AppData\Local\Julia-0.7.0\bin\julia
# C:\Users\micah\AppData\Local\Julia-1.0.1\bin\julia
# C:\Users\micah\AppData\Local\Julia-1.0.2\bin\julia
println("Loading Packages...")
runParallel=false

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
#using Plots;gr()
using Statistics
#using Dates
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVChelpers.jl")

N=100
path="C:\\Users\\micah\\Documents\\uvm\\Research\\Results\\N$(N)_new\\"
#path="C:\\Users\\micah\\Documents\\uvm\\Research\\Results\\N$(N)\\PWL\\"
#path="C:\\Users\\micah\\Documents\\uvm\\Research\\Results\\N$(N)_K\\"

datafile="jld2" #"mat" #"jld" #"n"
file="EVCscenarioN$(N)."*datafile

updateMethod="dualAscent" #dualAscent #fastAscent
maxIt=500 #500 for Dual Ascent
dualChk = 5e-4 #lamIt=0
primChk = 5e-4 # Ax-B=0
saveLogInd=[1,2,71,141,210]
noTlimit=false
forecastError=false
relaxedMode=0
slack=false
eqForm=false
tempAugment=false
ψ=-0

drawFig=false
saveResults=false
saveS=false
loadResults=false
silent=true
solverSilent=true

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
	#Tmax=500/1000
	Tmax=100 # Celsius
	num_homes=1000
	Tamb_amplitude = 18 #C
	Dload_error=0

	include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCscenario.jl")
    if saveS using FileIO end
	using Distributions
	evS=setupScenario(N;Tmax=Tmax,num_homes=num_homes,Dload_error=Dload_error,saveS=saveS,path=path)
end

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

stT1=Time(20,0)
endT1=Time(23,59)
stT2=Time(0,0)
endT2=Time(10,0)
Xlabels=vcat(collect(stT1:Dates.Second(round(evS.Ts)):endT1),collect(stT2:Dates.Second(round(evS.Ts)):endT2))
xticks=(1:40:evS.K,Dates.format.(Xlabels[1:40:evS.K],"HH:MM"))
