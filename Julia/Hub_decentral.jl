
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funHub.jl")
fname = "dDual_H$(H)"

if loadResults
	println("Reading in Dual Hub Sim")
	loadF=load(path*fname*".jld2")
	dSol=loadF["solution"]
else
	t0=hubS.t0
	e0=hubS.e0
	#prevLam=ones(hubS.K1+1,1)
	prevLam=3e3*ones(hubS.K1+1,1)
	#alpha0=3e2
	# alpha0 =1 #for kA
	roundSigFigs=16

	println("Running Dual Hub Sim")
	timeT=@elapsed dSol=hubDual(maxIt,hubS,cSol,mode,silent)
	if saveResults saveRun(path,fname,timeT, dSol) end
end

println("plotting....")
pd1alad=plot(dSol.Sn,xlabel="Time",ylabel="PEV SOC",legend=false,xlims=(0,hubS.K),ylims=(0,1))
if drawFig savefig(pd1alad,path*"J_decentral_ALADIN_SOC.png") end

pd2alad=plot(dSol.U,xlabel="Time",ylabel="PEV Current (kA)",legend=false,xlims=(0,hubS.K))
if drawFig savefig(pd2alad,path*"J_decentral_ALADIN_Curr.png") end

pd3alad=plot(hcat(dSol.Tactual[:,1],dSol.Tpred[:,1]),label=["Actual Temp" "PWL Temp"],xlims=(0,hubS.K),xlabel="Time",ylabel="Temp (K)")
plot!(pd3alad,1:hubS.K,hubS.Tmax*ones(hubS.K),label="XFRM Limit",line=(:dash,:red))
if drawFig savefig(pd3alad,path*"J_decentral_ALADIN_Temp.png") end

pd4alad=plot(hcat(cSol.Lam,dSol.Lam),xlabel="Time",ylabel=raw"Lambda ($/kA)",
             xlims=(0,hubS.K),labels=["Central" "Dual"])
if drawFig savefig(pd4alad,path*"J_decentral_ALADIN_Lam.png") end

aggU=plot(hcat(cSol.uSum,dSol.uSum),label=["Central" "ALAD"],
			xlims=(0,hubS.K),xlabel="Time",ylabel="PEV Current (kA)")

checkDesiredStates(dSol.Sn,evS.Kn,evS.Snmin)
