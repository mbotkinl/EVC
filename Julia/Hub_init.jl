# Wrapper Code to load scenario and most packages to run HUb code
# Micah Botkin-Levy

using JuMP
using Statistics
using Plots;pyplot()
using Printf
using Parameters
using SharedArrays
using Distributed
using LinearAlgebra
using Gurobi
using Suppressor

H=4                  # number of hubs
maxIt=1000           # maximum number of iterations per time step
auxChk = 1e3         # auxillary gap check for ALADIN
primChk = 1e-1       # coupling constraint gap check
datafile="jld2"      # "jld2" for reading in scenario file, "n" for creating new scenario
mode="PWL"           # PWL, NL, relax1 (not fully implemented)
eqForm=false         # ALADIN equality vs inequality form
coordinated=true     # true for tempature coordinated, false for coordinated (no temperature limit)
runParallel=false    # run local EV problems in parallel

# Parameters that control I/O
silent=true          # prevents output to console
solverSilent=true    # prevents feedback from JuMP solver
saveS=false          # saves scenario file
saveF=false          # draws output figures
loadResults=false    # loads results
saveResults=true     # saves result file

include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//structEVC.jl")
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVChelpers.jl")

path="C:\\Users\\micah\\Documents\\uvm\\Research\\Results\\hub4_new\\"
# path="C:\\Users\\micah\\Documents\\uvm\\Research\\Results\\hub4_100\\"
#path="C:\\Users\\micah\\Documents\\uvm\\Research\\Results\\hub4_100_old\\"

file="HubscenarioH$(H).jld2"
if datafile=="jld2"
	using FileIO
	println("Reading in Hub Scenario...")
	loadF=load(path*file)
	hubS=loadF["hubS"]
else #create scenario
	println("Creating Hub Scenario...")
	Nh=100
	Tmax=100
	include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCscenario.jl")
    if saveS using FileIO end
	using Distributions
	hubS=setupHubScenario(H,Nh,Tmax=Tmax,saveS=saveS,path=path)
end

# xticks for graphing
stT1=Time(20,0)
endT1=Time(23,59)
stT2=Time(0,0)
endT2=Time(10,0)
Xlabels=vcat(collect(stT1:Dates.Second(round(hubS.Ts)):endT1),collect(stT2:Dates.Second(round(hubS.Ts)):endT2))
xticks=(1:40:hubS.K,Dates.format.(Xlabels[1:40:hubS.K],"HH:MM"))
