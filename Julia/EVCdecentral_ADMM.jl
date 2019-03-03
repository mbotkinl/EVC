#EVC with ADMM for PWL convex relaxation
#Micah Botkin-Levy
#5/22/18

#formulation try 2
#u, sn, xt, and z are all in "x" v is the auxilliary variable corresponding to z in literature
#current constraint is coupling

@everywhere include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCpwl.jl")
#include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCpwl.jl")
fname = "dADMM_N$(N)"

if loadResults
	println("Reading in ADMM Sim")
	loadF=load(path*fname*".jld2")
	dSoladmm=loadF["solution"]
	dCMadmm=loadF["convMetrics"]
	convIt=loadF["convIt"]
else
	#initialize
    t0=evS.t0
    s0=evS.s0
	prevLam=10*ones(evS.K1+1,1)
	prevVz=-evS.deltaI/2*ones(evS.S*(evS.K1+1),1)
	prevVu=.02*ones(evS.N*(evS.K1+1),1)
	# prevVz=-0.01*ones(evS.S*(evS.K1+1),1)
	# prevVu=.01*ones(evS.N*(evS.K1+1),1)

	# ρADMMp = 1
	# ρDivRate=1.1
	#
	# ρADMMp = 10
	# ρDivRate=1.05

	ρADMMp = 100
	ρDivRate=1.01
	#149
	#lastGap  4.437284e-02 after 149 iterations
	#constGap 6.319717e-02 after 149 iterations


	roundSigFigs=12


	println("Running ADMM Sim")
	timeT=@elapsed dSoladmm,dCMadmm=pwlEVadmm(maxIt,evS,cSave,slack,roundSigFigs,silent)
	if saveResults saveRun(path,fname,timeT,dSoladmm,convMetrics=dCMadmm) end
end


println("plotting....")

pd1admm=plot(dSoladmm.Sn,xlabel="Time",ylabel="PEV SoC",legend=false,xlims=(0,evS.K),ylims=(0,1))
if drawFig savefig(pd1admm,path*"J_decentral_ADMM_SoC.png") end

pd2admm=plot(dSoladmm.Un,xlabel="Time",ylabel="PEV Current (kA)",legend=false,xlims=(0,evS.K))
if drawFig savefig(pd2admm,path*"J_decentral_ADMM_Curr.png") end

pd3admm=plot(hcat(dSoladmm.Tactual[:,1],dSoladmm.Tpred[:,1]),label=["Actual Temp" "PWL Temp"],xlims=(0,evS.K),xlabel="Time",ylabel="Temp (K)")
plot!(pd3admm,evS.Tmax*ones(evS.K),label="XFRM Limit",line=(:dash,:red))
if drawFig savefig(pd3admm,path*"J_decentral_ADMM_Temp.png") end

pd4admm=plot(hcat(cSol.lamCoupl,dSoladmm.lamCoupl),xlabel="Time",ylabel=raw"Lambda ($/kA)",
             xlims=(0,evS.K),labels=["Central" "ADMM"])
if drawFig savefig(pd4admm,path*"J_decentral_ADMM_Lam.png") end



#compare central and decentral current agg
aggU=plot(hcat(cSol.uSum,dSoladmm.uSum),label=["Central" "ADMM"],
			xlims=(0,evS.K),xlabel="Time",ylabel="PEV Current (kA)")


checkDesiredStates(dSoladmm.Sn,evS.Kn,evS.Snmin)
