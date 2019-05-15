# Run file for hub ALADIN simulation
# Micah Botkin-Levy

include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funHub.jl")
eqString=if eqForm "_eq" else "_ineq" end
fname = "dALADIN_H$(H)"*eqString

if loadResults
	println("Reading in ALAD Sim")
	loadF=load(path*fname*".jld2")
	hubS=loadF["scenario"]
	dSolalad=loadF["solution"]
else
	# initialize
	t0=hubS.t0
	e0=hubS.e0
	prevLam=8e3*ones(hubS.K1+1,1)
    prevVt=hubS.Tmax*ones(hubS.K1+1,1)
    prevVz=hubS.deltaI*ones((hubS.K1+1),hubS.S)
    prevVu=repeat(maximum(hubS.uMax,dims=1),outer=[(hubS.K1+1),1])
    prevVe=repeat(maximum(hubS.eMax,dims=1),outer=[(hubS.K1+1),1])
	prevVd=repeat(maximum(hubS.slackMax,dims=1),outer=[(hubS.K1+1),1])

	ρALADp = .01
	ρRate=1.1
	ρALADmax=1e6
	roundSigFigs=16

	println("Running ALAD Sim")
	timeT=@elapsed dSolalad=hubALAD(maxIt,hubS,cSol,mode,eqForm,silent)
	if saveResults saveRun(path,fname,timeT, dSolalad) end
end

hubLabels=permutedims(["Hub $(h)" for h=1:hubS.H])

pd1alad=plot(dSolalad.E,xlabel="",ylabel="Energy (MWh)",seriestype=:line,labels=hubLabels,xticks=xticks, xlims=(0,hubS.K))
plot!(pd1alad,hubS.eMax,label="Hub Max",line=(:dash))

pd2alad=plot(dSolalad.U,xlabel="",ylabel="Hub Current (kA)",legend=false,xticks=xticks, xlims=(0,hubS.K))

pd3alad=plot(dSolalad.Tactual,label="ALAD",xlabel="",ylabel="Temp (K)", xlims=(0,hubS.K))
plot!(pd3alad,cSol.Tactual,label="Central")
plot!(pd3alad,hubS.Tmax*ones(hubS.K),label="XFRM Limit",line=(:dash,:red),xticks=xticks)

pd4alad=plot(dSolalad.Lam,xlabel="Time",ylabel=raw"Lambda",label="ALAD",xticks=xticks, xlims=(0,hubS.K))
plot!(pd4alad,cSol.Lam,label="Central")

aggU=plot(hcat(cSol.uSum,dSolalad.uSum),label=["Central" "ALAD"],
			xlims=(0,hubS.K),xlabel="",ylabel="Current (kA)")

aggE=plot(hcat(sum(cSol.E,dims=2),sum(dSolalad.E,dims=2)),label=["Central" "ALAD"], xlims=(0,hubS.K),xlabel="",ylabel="Energy (MWh)")
plot!(aggE,sum(hubS.eMax,dims=2),label="Max Energy",line=(:dash,:red),xticks=xticks)

h1alad=plot(aggE,aggU,pd3alad,pd4alad,layout=(4,1))
lowRes=true
if lowRes
    pubPlot(h1alad,thickscale=0.4,sizeWH=(400,300),dpi=40)
else
    pubPlot(h1alad,thickscale=0.8,sizeWH=(1000,600),dpi=100)
end
if saveF savefig(h1alad,path*"hubPlot1ALAD.png") end
