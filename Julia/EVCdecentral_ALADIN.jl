#EVC with ALADIN for PWL convex relaxation
#Micah Botkin-Levy
#5/22/18

#formulation try 2
#u, sn, xt, and z are all in "y", v is the auxilliary variable corresponding to x in literature
#current constraint is coupling

@everywhere include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCpwl.jl")
eqString=if eqForm "_eq" else "_ineq" end
fname = "dALADIN_N$(N)"*eqString

if loadResults
	println("Reading in ALAD Sim")
	loadF=load(path*fname*".jld2")
	dLogalad=loadF["solution"]
	dCMalad=loadF["convMetrics"]
	convIt=loadF["convIt"]
else
	t0=evS.t0
	s0=evS.s0
	prevLam=2e3*ones(evS.K1+1,1)
    prevVt=evS.Tmax*ones(evS.K1+1,1)
    prevVz=.01*ones(evS.S*(evS.K1+1),1)
    prevVu=.01*ones(evS.N*(evS.K1+1),1)
    prevVs=.5*ones(evS.N*(evS.K1+1),1)
	ρALADp = 1e3

	roundSigFigs=8

	println("Running ALAD Sim")
	timeT=@elapsed dSolalad,dCMalad=pwlEValad(maxIt,evS,cSave,slack,eqForm,silent)
	if saveResults saveRun(path,fname,timeT,dSolalad,convMetrics=dCMalad) end
end

println("plotting....")
pd1alad=plot(dSolalad.Sn,xlabel="Time",ylabel="PEV SOC",legend=false,xlims=(0,evS.K),ylims=(0,1))
if drawFig savefig(pd1alad,path*"J_decentral_ALADIN_SOC.png") end

pd2alad=plot(dSolalad.Un,xlabel="Time",ylabel="PEV Current (kA)",legend=false,xlims=(0,evS.K))
if drawFig savefig(pd2alad,path*"J_decentral_ALADIN_Curr.png") end

pd3alad=plot(hcat(dSolalad.Tactual[:,1],dSolalad.Tpwl[:,1])*1000,label=["Actual Temp" "PWL Temp"],xlims=(0,evS.K),xlabel="Time",ylabel="Temp (K)")
plot!(pd3alad,1:evS.K,evS.Tmax*ones(evS.K)*1000,label="XFRM Limit",line=(:dash,:red))
if drawFig savefig(pd3alad,path*"J_decentral_ALADIN_Temp.png") end

pd4alad=plot(hcat(cSol.lamCoupl,dSolalad.lamCoupl),xlabel="Time",ylabel=raw"Lambda ($/kA)",
             xlims=(0,evS.K),labels=["Central" "ALADIN"])
if drawFig savefig(pd4alad,path*"J_decentral_ALADIN_Lam.png") end

aggU=plot(hcat(cSol.uSum,dSolalad.uSum),label=["Central" "ALAD"],
			xlims=(0,evS.K),xlabel="Time",ylabel="PEV Current (kA)")

checkDesiredStates(dSolalad.Sn,evS.Kn,evS.Snmin)
