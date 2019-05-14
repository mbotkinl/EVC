# wrapper to run PWL centralized problem (must run EVC_init.jl first)
#Micah Botkin-Levy
#4/8/18

if runParallel
	@everywhere using Gurobi
	@everywhere include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCpwl.jl")
else
	using Gurobi
	include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCpwl.jl")
end

fname = "central_N$(N)"
if noTlimit
	fname = fname*"_noLim"
end

if loadResults
	println("Reading in Central Sim")
	loadF=load(path*fname*".jld2")
	cSol=loadF["solution"]
	cSave=loadF["centralLog"]
else
	#initial states
    t0=evS.t0
    s0=evS.s0

	roundSigFigs=12

	println("Running Central Sim")
	timeT=@elapsed cSol, cSave=pwlEVcentral(evS,slack,silent)
	if saveResults saveRun(path,fname,timeT,cSol,cSave=cSave) end
end


println("plotting....")

p1=plot(cSol.Sn,xlabel="Time",ylabel="PEV SoC",legend=false,xlims=(1,evS.K))
if drawFig savefig(p1,path*"J_central_SoC.png") end

p2=plot(cSol.Un,xlabel="Time",ylabel="PEV Current (kA)",legend=false,xlims=(1,evS.K))
if drawFig savefig(p2,path*"J_central_Curr.png") end

p3=plot(hcat(cSol.Tpred,cSol.Tactual),label=["Pred Temp" "Actual Temp"],xlabel="Time",ylabel="Temp (C)")
plot!(p3,1:evS.K,evS.Tmax*ones(evS.K),label="XFRM Limit",line=(:dash,:red),xticks=xticks)
if drawFig savefig(p3,path*"J_central_Temp.png") end

p4b=plot(cSol.lamTemp,xlabel="Time",ylabel=raw"Lambda ($/K)",legend=false)
p4=plot(cSol.lamCoupl,xlabel="Time",ylabel=raw"Lambda",legend=false,xticks=xticks)
if drawFig savefig(p4,path*"J_central_Lam.png") end


aggU=plot(hcat(cSol.uSum,cSol.Itotal,evS.iD_actual),label=["Central" "Total" "iD"],xlabel="Time",ylabel="Current (kA)")
plot!(aggU,1:evS.K,evS.ItotalMax*ones(evS.K),label="XFRM Limit",line=(:dash,:red),xticks=xticks)

timePlot = plot(cSol.timeT,xlabel="Time",ylabel="Computational Time (s)",legend=false, )

checkDesiredStates(cSol.Sn,evS.Kn,evS.Snmin)
