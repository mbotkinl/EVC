
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funHub.jl")
eqString=if eqForm "_eq" else "_ineq" end
fname = "dALADIN_H$(H)"*eqString

if loadResults
	println("Reading in ALAD Sim")
	loadF=load(path*fname*".jld2")
	hubS=loadF["scenario"]
	dSol=loadF["solution"]
else
	t0=hubS.t0
	e0=hubS.e0
	prevLam=2e3*ones(hubS.K1+1,1)
    prevVt=hubS.Tmax*ones(hubS.K1+1,1)
    prevVz=.01*ones((hubS.K1+1),hubS.S)
    prevVu=.01*ones((hubS.K1+1),hubS.H)
    prevVe=ones((hubS.K1+1),hubS.H)
	prevVd=ones((hubS.K1+1),hubS.H)
	ρALADp = 1e3

	println("Running ALAD Sim")
	timeT=@elapsed dSolalad=hubALAD(maxIt,hubS,cSol,mode,eqForm,silent)
	if saveResults saveRun(path,fname,timeT, evS, dSolalad) end
end
#
# println("plotting....")
# pd1alad=plot(dSolalad.Sn,xlabel="Time",ylabel="PEV SOC",legend=false,xlims=(0,evS.K),ylims=(0,1))
# if drawFig savefig(pd1alad,path*"J_decentral_ALADIN_SOC.png") end
#
# pd2alad=plot(dSolalad.Un,xlabel="Time",ylabel="PEV Current (kA)",legend=false,xlims=(0,evS.K))
# if drawFig savefig(pd2alad,path*"J_decentral_ALADIN_Curr.png") end
#
# pd3alad=plot(hcat(dSolalad.Tactual[:,1],dSolalad.Tpwl[:,1])*1000,label=["Actual Temp" "PWL Temp"],xlims=(0,evS.K),xlabel="Time",ylabel="Temp (K)")
# plot!(pd3alad,1:evS.K,evS.Tmax*ones(evS.K)*1000,label="XFRM Limit",line=(:dash,:red))
# if drawFig savefig(pd3alad,path*"J_decentral_ALADIN_Temp.png") end
#
# pd4alad=plot(hcat(cSol.lamCoupl,dSolalad.lamCoupl),xlabel="Time",ylabel=raw"Lambda ($/kA)",
#              xlims=(0,evS.K),labels=["Central" "ALADIN"])
# if drawFig savefig(pd4alad,path*"J_decentral_ALADIN_Lam.png") end
#
# aggU=plot(hcat(cSol.uSum,dSolalad.uSum),label=["Central" "ALAD"],
# 			xlims=(0,evS.K),xlabel="Time",ylabel="PEV Current (kA)")
#
# checkDesiredStates(dSolalad.Sn,evS.Kn,evS.Snmin)
