using Distributions
using Gurobi

#pull out a few key variables
H=hubS.H
S=hubS.S

packLen=2 #number of time steps
mttr=hubS.Ts*packLen #300
setSOC=0.1
roundSigFigs=12

include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCpem.jl")

fname = "d_PEM_H$(H)"

if loadResults
	println("Reading in PEM Sim")
	loadF=load(path*fname*".jld2")
	pemSol=loadF["solution"]
else
	println("Running PEM Sim")
	timeT=@elapsed pemSol,Req=pemHub(hubS,silent)
	if saveResults saveRun(path,fname,timeT,pemSol) end
end
println("plotting....")

p1pem=plot(pemSol.E,xlabel="Time",ylabel="Hub Energy",legend=false,xlims=(0,hubS.K))
if drawFig savefig(p1pem,path*"J_PEM_SoC.png") end

p2pem=plot(pemSol.U,xlabel="Time",ylabel="HUb Current (kA)",legend=false,xlims=(0,hubS.K))
if drawFig savefig(p2pem,path*"J_PEM_Curr.png") end

p3pem=plot(hcat(cSol.Tactual,pemSol.Tactual),label=["Central" "PEM"],xlims=(0,hubS.K),xlabel="Time",ylabel="Temp (K)")
plot!(p3pem,hubS.Tmax*ones(hubS.K),label="XFRM Limit",line=(:dash,:red))
if drawFig savefig(p3pem,path*"J_PEM_Temp.png") end

requests = [count(Req[k,:].>0) for k=1:hubS.K]
reqPlot=plot(requests,xlims=(0,hubS.K),xlabel="Time",ylabel="Number of Requests")

aggUpem=plot(hcat(cSol.uSum,pemSol.uSum),label=["Central" "PEM"],
			xlims=(0,hubS.K),xlabel="Time",ylabel="PEV Current (kA)")
if drawFig savefig(aggUpem,path*"J_PEM_uSum.png") end


aggU=plot(hcat(cSol.uSum .+ hubS.iD_actual,pemSol.Itotal,hubS.iD_actual),label=["Central" "Total" "iD"],xlims=(0,hubS.K),xlabel="Time",ylabel="Current (kA)")
plot!(aggU,1:hubS.K,hubS.ItotalMax*ones(hubS.K),label="XFRM Limit",line=(:dash,:red))
