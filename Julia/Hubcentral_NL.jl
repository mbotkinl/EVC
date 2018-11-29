#test hub model
using JuMP
using Ipopt
using Statistics
using Plots;pyplot()
using Printf
using Parameters
using SharedArrays

# check all indexing????
#especially stepI +k for iD and Tamb
# fix for multiple hubs ????

N=3
silent=true

include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//structEVC.jl")
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCscenario.jl")
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funHub.jl")

println("Creating Hub Scenario...")
hubS=setupHubScenario(N)
H=hubS.H
K=hubS.K
N=hubS.N

#initialize
t0=hubS.t0
e0=hubS.e0
cSol=centralHubSolutionStruct()

for stepI=1:K
    runHubCentralStep(stepI,hubS,cSol,silent)
end


p1nl=plot(cSol.E,xlabel="Time",ylabel="Hub Energy",legend=false,xlims=(1,K))
p2nl=plot(cSol.U,xlabel="Time",ylabel="Hub Current (kA)",legend=false,xlims=(1,K))

p3nl=plot(cSol.T,label="XFRM Temp",xlims=(1,K),xlabel="Time",ylabel="Temp (K)")
plot!(p3nl,hubS.Tmax*ones(K),label="XFRM Limit",line=(:dash,:red))

p4nl=plot(cSol.Lam,label="Time",ylabel=raw"Lambda ($/kA)",xlims=(1,K),legend=false)


t=false
