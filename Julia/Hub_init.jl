#test hub model
using JuMP
using Statistics
using Plots;pyplot()
using Printf
using Parameters
using SharedArrays
using Distributed
using LinearAlgebra
#using Ipopt
using Gurobi
using Suppressor

# check all indexing????
#especially stepI +k for iD and Tamb
H=4
maxIt=5000
auxChk = 1e1 #lamIt=0
primChk = 1e-1 # Ax-B=0s
datafile="jld2"
mode="PWL"
silent=true
solverSilent=true
saveS=false
saveF=false
loadResults=false
saveResults=true
eqForm=false
coordinated=true
runParallel=false
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
	#Dload_amplitude=20

	include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCscenario.jl")
    if saveS using FileIO end
	using Distributions
	hubS=setupHubScenario(H,Nh,Tmax=Tmax,saveS=saveS,path=path)
end

stT1=Time(20,0)
endT1=Time(23,59)
stT2=Time(0,0)
endT2=Time(10,0)
Xlabels=vcat(collect(stT1:Dates.Second(round(hubS.Ts)):endT1),collect(stT2:Dates.Second(round(hubS.Ts)):endT2))
xticks=(1:40:hubS.K,Dates.format.(Xlabels[1:40:hubS.K],"HH:MM"))
