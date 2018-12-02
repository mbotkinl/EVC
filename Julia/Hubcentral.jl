#test hub model
using JuMP
using Statistics
using Plots;pyplot()
using Printf
using Parameters
using SharedArrays
#using Ipopt
using Gurobi

# check all indexing????
#especially stepI +k for iD and Tamb

H=3
Nh=20
Tmax=.393
mode="PWL"
silent=true
saveS=false
saveF=false

include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//structEVC.jl")
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCscenario.jl")
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funHub.jl")
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVChelpers.jl")

path="C:\\Users\\micah\\Documents\\uvm\\Research\\Results\\hubSimple\\"

println("Creating Hub Scenario...")
hubS=setupHubScenario(H,Nh,Tmax=Tmax,saveS=saveS,path=path)
H=hubS.H
K=hubS.K
Nh=hubS.Nh

#initialize
t0=hubS.t0
e0=hubS.e0
cSol=centralHubSolutionStruct()

for stepI=1:K
    runHubCentralStep(stepI,hubS,cSol,mode,silent)
end

stT1=Time(20,0)
endT1=Time(23,0)
stT2=Time(0,0)
endT2=Time(10,0)
Xlabels=vcat(collect(stT1:Dates.Second(round(hubS.Ts)):endT1),collect(stT2:Dates.Second(round(hubS.Ts)):endT2))
#Xlabels=vcat(collect(stT1:Dates.Minute(3):endT1),collect(stT2:Dates.Minute(3):endT2))
xticks=(1:40:K,Dates.format.(Xlabels[1:40:K],"HH:MM"))


p1nl=plot(1:K,cSol.E,xlabel="Time",ylabel="Energy (kWh)",label="Hub Energy",xlims=(0,K),seriestype=:bar)
plot!(p1nl,1:K,cSol.E_depart,label="Depart Energy",seriestype=:scatter,markersize=10)
plot!(p1nl,1:K,cSol.E_arrive,label="Arrive Energy",seriestype=:scatter,markersize=10)
plot!(twinx(),1:K,cSol.U,label="Hub Current",seriestype=:line,seriescolor=:red,legend=false,ylabel="Current (kA)",
        xlims=(0,K))

p2nl=plot(cSol.U,xlabel="Time",ylabel="Hub Current (kA)",legend=false,xlims=(0,K))

p3nl=plot(cSol.T*1000,label="XFRM Temp",xlims=(0,K),xlabel="Time",ylabel="Temp (K)")
plot!(p3nl,hubS.Tmax*ones(K)*1000,label="XFRM Limit",line=(:dash,:red))

p4nl=plot(cSol.Lam,label="Time",ylabel=raw"Lambda ($/kA)",xlims=(0,K),legend=false)


h1=plot(p1nl,p3nl,p4nl,layout=(3,1))
lowRes=true
if lowRes
    pubPlot(h1,thickscale=0.4,sizeWH=(400,300),dpi=40)
else
    pubPlot(h1,thickscale=0.8,sizeWH=(800,600),dpi=100)
en
if saveF savefig(h1,path*"hubPlot.png") end
