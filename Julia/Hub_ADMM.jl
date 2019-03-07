
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funHub.jl")
fname = "dADMM_H$(H)"

if loadResults
	println("Reading in ADMM Hub Sim")
	loadF=load(path*fname*".jld2")
	hubS=loadF["scenario"]
	dSoladmm=loadF["solution"]
else
	t0=hubS.t0
	e0=hubS.e0
	prevLam=ones(hubS.K1+1,1)
	prevVz=-hubS.deltaI*ones((hubS.K1+1),hubS.S)
	#prevVu=repeat(maximum(hubS.uMax,dims=1),outer=[(hubS.K1+1),1])
	prevVu=.01*ones((hubS.K1+1),hubS.H)
	#ρADMMp = 1e2

	ρADMMp = 1e2
	ρDivRate=1.02
	#ρDivRate=1.02


	println("Running ADMM Hub Sim")
	timeT=@elapsed dSoladmm=hubADMM(maxIt,hubS,cSol,mode,silent)
	if saveResults saveRun(path,fname,timeT, dSoladmm) end
end


stT1=Time(20,0)
endT1=Time(23,59)
stT2=Time(0,0)
endT2=Time(10,0)
Xlabels=vcat(collect(stT1:Dates.Second(round(hubS.Ts)):endT1),collect(stT2:Dates.Second(round(hubS.Ts)):endT2))
xticks=(1:40:hubS.K,Dates.format.(Xlabels[1:40:hubS.K],"HH:MM"))
hubLabels=permutedims(["Hub $(h)" for h=1:hubS.H])

pd1admm=plot(dSoladmm.E,xlabel="",ylabel="Energy (MWh)",seriestype=:line,labels=hubLabels,xticks=xticks, xlims=(0,hubS.K))
plot!(pd1admm,hubS.eMax,label="Hub Max",line=(:dash))

pd2admm=plot(dSoladmm.U,xlabel="",ylabel="Hub Current (kA)",legend=false,xticks=xticks, xlims=(0,hubS.K))

pd3admm=plot(dSoladmm.Tactual*1000,label="ADMM",xlabel="",ylabel="Temp (K)", xlims=(0,hubS.K))
plot!(pd3admm,cSol.Tactual*1000,label="Central")
plot!(pd3admm,hubS.Tmax*ones(hubS.K)*1000,label="XFRM Limit",line=(:dash,:red),xticks=xticks)

pd4admm=plot(dSoladmm.Lam,xlabel="Time",ylabel=raw"Lambda ($/kA)",label="ADMM",xticks=xticks, xlims=(0,hubS.K))
plot!(pd4admm,cSol.Lam,label="Central")

aggU=plot(hcat(cSol.uSum,dSoladmm.uSum),label=["Central" "ADMM"],
			xlims=(0,hubS.K),xlabel="",ylabel="PEV Current (kA)")
# plot!(aggU,sum(hubS.uMax,dims=2),label="Max Current",line=(:dash,:red),xticks=xticks)

aggE=plot(hcat(sum(cSol.E,dims=2),sum(dSoladmm.E,dims=2)),label=["Central" "ADMM"], xlims=(0,hubS.K),xlabel="",ylabel="Energy (MWh)")
plot!(aggE,sum(hubS.eMax,dims=2),label="Max Energy",line=(:dash,:red),xticks=xticks)
