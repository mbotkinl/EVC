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

# check all indexing????
#especially stepI +k for iD and Tamb
H=4
maxIt=1000
datafile="n"
mode="PWL"
silent=true
saveS=false
saveF=false
loadResults=false
saveResults=false
eqForm=false
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//structEVC.jl")
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVChelpers.jl")

path="C:\\Users\\micah\\Documents\\uvm\\Research\\Results\\hub4\\"
file="HubscenarioH$(H).jld2"
if datafile=="jld2"
	using FileIO
	println("Reading in Hub Scenario...")
	loadF=load(path*file)
	hubS=loadF["hubS"]
else #create scenario
	println("Creating Hub Scenario...")
	Nh=40
	Tmax=.393
	#Dload_amplitude=20

	include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCscenario.jl")
    if saveS using FileIO end
	using Distributions
	hubS=setupHubScenario(H,Nh,Tmax=Tmax,saveS=saveS,path=path)
end
