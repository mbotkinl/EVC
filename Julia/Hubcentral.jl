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


mode="PWL"
silent=true
saveS=false
saveF=false

include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//structEVC.jl")
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funHub.jl")
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVChelpers.jl")

path="C:\\Users\\micah\\Documents\\uvm\\Research\\Results\\hub4\\"
file="HubscenarioH4.jld2"
if datafile=="jld2"
	using FileIO
	println("Reading in Hub Scenario...")
	loadF=load(path*file)
	hubS=loadF["hubS"]

else #create scenario
	println("Creating Hub Scenario...")
	H=4
	Nh=25
	Tmax=.393
	#Dload_amplitude=20

	include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCscenario.jl")
    if saveS using FileIO end
	using Distributions
	hubS=setupHubScenario(H,Nh,Tmax=Tmax,saveS=saveS,path=path)
end


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
endT1=Time(23,59)
stT2=Time(0,0)
endT2=Time(10,0)
Xlabels=vcat(collect(stT1:Dates.Second(round(hubS.Ts)):endT1),collect(stT2:Dates.Second(round(hubS.Ts)):endT2))
#Xlabels=vcat(collect(stT1:Dates.Minute(3):endT1),collect(stT2:Dates.Minute(3):endT2))
xticks=(1:40:K,Dates.format.(Xlabels[1:40:K],"HH:MM"))



hubLabels=permutedims(["Hub $(h)" for h=1:H])
p1=plot(1:K,cSol.E,xlabel="",ylabel="Energy (kWh)",xlims=(0,K),seriestype=:line,labels=hubLabels,xticks=xticks)

sumPlot=plot(1:K,sum(cSol.E,dims=2),xlabel="",ylabel="Energy (kWh)",label="Hub Energy",xlims=(0,K),seriestype=:bar,xticks=xticks)
plot!(sumPlot,1:K,sum(cSol.E_depart,dims=2),label="Depart Energy",seriestype=:scatter,markersize=10)
plot!(sumPlot,1:K,sum(cSol.E_arrive,dims=2),label="Arrive Energy",seriestype=:scatter,markersize=10)
plot!(twinx(),1:K,sum(cSol.U,dims=2),label="Hub Current",seriestype=:line,seriescolor=:red,legend=false,ylabel="Current (kA)",
        xlims=(0,K),xticks=xticks)

p2=plot(cSol.U,xlabel="",ylabel="Hub Current (kA)",legend=false,xlims=(0,K),xticks=xticks)

p3=plot(cSol.T*1000,label="XFRM Temp",xlims=(0,K),xlabel="",ylabel="Temp (K)")
plot!(p3,hubS.Tmax*ones(K)*1000,label="XFRM Limit",line=(:dash,:red),xticks=xticks)

p4=plot(cSol.Lam,label="Time",ylabel=raw"Lambda ($/kA)",xlims=(0,K),legend=false,xticks=xticks)


h1=plot(p1,p2,p3,p4,layout=(4,1))
lowRes=true
if lowRes
    pubPlot(h1,thickscale=0.4,sizeWH=(400,300),dpi=40)
else
    pubPlot(h1,thickscale=0.8,sizeWH=(800,600),dpi=100)
end
if saveF savefig(h1,path*"hubPlot1.png") end

h2=plot(sumPlot,p3,p4,layout=(3,1))
lowRes=true
if lowRes
    pubPlot(h2,thickscale=0.4,sizeWH=(400,300),dpi=40)
else
    pubPlot(h2,thickscale=0.8,sizeWH=(800,600),dpi=100)
end
if saveF savefig(h2,path*"hubPlot2.png") end
