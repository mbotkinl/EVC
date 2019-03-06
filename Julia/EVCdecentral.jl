


#Micah Botkin-Levy
#4/10/18
@everywhere include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCpwl.jl")
fname = "d_$(updateMethod)_N$(N)"

if loadResults
	println("Reading in Dual Decomp Sim")
	loadF=load(path*fname*".jld2")
	dLog=loadF["solution"]
	dCM=loadF["convMetrics"]
	convIt=loadF["convIt"]
else
	#initialize
    t0=evS.t0
    s0=evS.s0
	prevLam=1e5*ones(evS.K1+1,1)

	alpha0 = 1e4 #for kA
	#prevLam=round.(cSol.lamCoupl[1:evS.K1+1],digits=4)
	roundSigFigs=30

	println("Running Dual Decomp Sim")
	timeT=@elapsed dSol,dCM=pwlEVdual(maxIt,updateMethod,evS,cSave,slack,roundSigFigs,silent)
	# s=Symbol(@sprintf("dCM_%s",updateMethod))
	# v=Symbol(@sprintf("dCM"))
	# @eval(($s)=($v))
	if saveResults saveRun(path,fname,timeT,dSol,convMetrics=dCM) end
end

println("plotting....")
#solution plots
pd1=plot(dSol.Sn,xlabel="Time",ylabel="PEV SoC",legend=false,xlims=(0,evS.K),ylims=(0,1))
if drawFig savefig(p1d,path*"J_"*updateMethod*"_SoC.png") end

pd2=plot(dSol.Un,xlabel="Time",ylabel="PEV Current (kA)",legend=false,xlims=(0,evS.K))
if drawFig savefig(p2d,path*"J_"*updateMethod*"_Curr.png") end

pd3=plot(dSol.Tactual[:,1],label="Actual Temp",xlims=(0,evS.K),xlabel="Time",ylabel="Temp (K)")
plot!(pd3,1:evS.K,evS.Tmax*ones(evS.K),label="XFRM Limit",line=(:dash,:red))
if updateMethod=="dualAscent" plot!(pd3,1:evS.K,dSol.Tpred[:,1]*1000,label="PWL Temp") end
if drawFig savefig(pd3,path*"J_"*updateMethod*"_Temp.png") end

if updateMethod=="fastAscent"
	lamLabel=raw"Lambda ($/K)"
else
	lamLabel=raw"Lambda ($/kA)"
end

pd4=plot(hcat(cSol.lamCoupl,dSol.lamCoupl[:,1]),xlabel="Time",xlims=(0,evS.K),labels=["Central" "Dual"])
if drawFig savefig(pd4,path*"J_"*updateMethod*"_Lam.png") end

#compare central and decentral current agg
aggU=plot(hcat(cSol.uSum,dSol.uSum),label=["Central" "Decentral"],
			xlims=(0,evS.K),xlabel="Time",ylabel="PEV Current (kA)")

checkDesiredStates(dSol.Sn,evS.Kn,evS.Snmin)


#draw(PNG(path*"aggPlot_fast.png", 13inch, 8inch), aggU)
