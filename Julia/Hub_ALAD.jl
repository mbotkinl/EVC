
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
	prevLam=ones(hubS.K1+1,1)
    prevVt=hubS.Tmax*ones(hubS.K1+1,1)
    prevVz=hubS.deltaI*ones((hubS.K1+1),hubS.S)
    prevVu=repeat(maximum(hubS.uMax,dims=1),outer=[(hubS.K1+1),1])
    prevVe=repeat(maximum(hubS.eMax,dims=1),outer=[(hubS.K1+1),1])
	prevVd=repeat(maximum(hubS.slackMax,dims=1),outer=[(hubS.K1+1),1])
	ρALADp = 1e-2


	println("Running ALAD Sim")
	timeT=@elapsed dSolalad=hubALAD(maxIt,hubS,cSol,mode,eqForm,silent)
	if saveResults saveRun(path,fname,timeT, hubS, dSolalad) end
end

stT1=Time(20,0)
endT1=Time(23,59)
stT2=Time(0,0)
endT2=Time(10,0)
Xlabels=vcat(collect(stT1:Dates.Second(round(hubS.Ts)):endT1),collect(stT2:Dates.Second(round(hubS.Ts)):endT2))
xticks=(1:40:hubS.K,Dates.format.(Xlabels[1:40:hubS.K],"HH:MM"))
hubLabels=permutedims(["Hub $(h)" for h=1:hubS.H])

p1=plot(dSolalad.E,xlabel="",ylabel="Energy (MWh)",seriestype=:line,labels=hubLabels,xticks=xticks)
plot!(hubS.eMax,label="Hub Max",line=(:dash))

p2=plot(dSolalad.U,xlabel="",ylabel="Hub Current (kA)",legend=false,xticks=xticks)

p3=plot(dSolalad.Tactual*1000,label="XFRM Temp",xlabel="",ylabel="Temp (K)")
plot!(p3,hubS.Tmax*ones(hubS.K)*1000,label="XFRM Limit",line=(:dash,:red),xticks=xticks)

p4=plot(dSolalad.Lam,label="Time",ylabel=raw"Lambda ($/kA)",legend=false,xticks=xticks)
plot!(p4,cSol.Lam,label="Central")

aggU=plot(hcat(cSol.uSum,dSolalad.uSum),label=["Central" "ALAD"],
			xlims=(0,hubS.K),xlabel="Time",ylabel="PEV Current (kA)")

#checkDesiredStates(dSolalad.Sn,evS.Kn,evS.Snmin)
