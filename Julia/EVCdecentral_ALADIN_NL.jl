
#EVC with ALADIN for PWL convex relaxation
#Micah Botkin-Levy
#5/22/18

#formulation try 2
#u, sn, xt, and z are all in "y", v is the auxilliary variable corresponding to x in literature
#current constraint is coupling
@everywhere include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCnl.jl")

eqString=if eqForm "_eq" else "_ineq" end
relaxString= "_R$(relaxedMode)"
fname = "dALADIN_NL_N$(N)"*eqString*relaxString

if loadResults
	println("Reading in NL AlAD Sim")
	loadF=load(path*fname*".jld2")
	dSolaladnl=loadF["solution"]
	dCMaladnl=loadF["convMetrics"]
else
	t0=evS.t0
	s0=evS.s0

	#prevLam=cSolnl.lamCoupl[1:evS.K1+1]
	prevLam=2e3*ones(evS.K1+1,1)
    prevVt=evS.Tmax*ones(evS.K1+1,1)
    prevVi=evS.ItotalMax*ones(evS.K1+1,1)
    prevVu=.01*ones(evS.N*(evS.K1+1),1)
    prevVs=.5*ones(evS.N*(evS.K1+1),1)
	ρALADp = 1

	roundSigFigs=12

	println("Running NL AlAD Sim")
	timeT=@elapsed dSolaladnl,dCMaladnl=nlEValad(maxIt,evS,cSavenl,slack,eqForm,roundSigFigs,silent)
	if saveResults saveRun(path,fname,timeT,dSolaladnl,convMetrics=dCMaladnl) end
end

println("plotting....")

pd1aladNL=plot(dSolaladnl.Sn,xlabel="Time",ylabel="PEV SoC",legend=false,xlims=(0,evS.K),ylims=(0,1))
if drawFig==1 savefig(pd1aladNL,path*"J_decentralNL_ALADIN_SOC.png") end

pd2aladNL=plot(dSolaladnl.Un,xlabel="Time",ylabel="PEV Current (kA)",legend=false,xlims=(0,evS.K))
if drawFig==1 savefig(pd2aladNL,path*"J_decentralNL_ALADIN_Curr.png") end

pd3aladNL=plot(hcat(dSolaladnl.Tactual[:,1],dSolaladnl.Tpred[:,1])*1000,label=["Actual Temp" "Pred Temp"],xlims=(0,evS.K),xlabel="Time",ylabel="Temp (K)")
plot!(pd3aladNL,1:evS.K,evS.Tmax*ones(evS.K)*1000,label="XFRM Limit",line=(:dash,:red))
if drawFig==1 savefig(pd3aladNL,path*"J_decentralNL_ALADIN_Temp.png") end

pd4aladNL=plot(hcat(cSolnl.lamCoupl,dSolaladnl.lamCoupl[:,convIt]),xlabel="Time",ylabel=raw"Lambda ($/kA)",labels=["Central" "ALADIN"])
if drawFig==1 savefig(pd4aladNL,path*"J_decentralNL_ALADIN_Lam.png") end

aggU=plot(hcat(cSolnl.uSum,dSolaladnl.uSum),label=["Central" "ALAD"],
			xlims=(0,evS.K),xlabel="Time",ylabel="PEV Current (kA)")

checkDesiredStates(dSolaladnl.Sn,evS.Kn,evS.Snmin)
