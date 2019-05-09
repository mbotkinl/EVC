
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funHub.jl")
fname = "dADMM_H$(H)"

if loadResults
	println("Reading in ADMM Hub Sim")
	loadF=load(path*fname*".jld2")
	dSoladmm=loadF["solution"]
else
	t0=hubS.t0
	e0=hubS.e0
	prevLam=8e3*ones(hubS.K1+1,1)
	prevVz=-hubS.deltaI*ones((hubS.K1+1),hubS.S)
	#prevVu=repeat(maximum(hubS.uMax,dims=1),outer=[(hubS.K1+1),1])
	prevVu=.01*ones((hubS.K1+1),hubS.H)

	# ρADMMp = 5e1
	# ρDivRate=1.1
	# maxRho=1e7
	#
	# ρADMMp = 8e3
	# ρDivRate=1
	# maxRho=1e10
	#
	# ρADMMp = 5e3
	# ρDivRate=1.01
	# maxRho=1e10

	# ρDivRate=1.01
	# maxRho=1e10

	# ρADMMp = 5e3
	# ρDivRate=1.05
	# maxRho=1e8

	ρADMMp = 1e3
	ρDivRate=1.001
	maxRho=1e8


	roundSigFigs=15


	println("Running ADMM Hub Sim")
	timeT=@elapsed dSoladmm=hubADMM(maxIt,hubS,cSol,mode,silent)
	if saveResults saveRun(path,fname,timeT, dSoladmm) end
end

hubLabels=permutedims(["Hub $(h)" for h=1:hubS.H])

pd1admm=plot(dSoladmm.E,xlabel="",ylabel="Energy (MWh)",seriestype=:line,labels=hubLabels,xticks=xticks, xlims=(0,hubS.K))
plot!(pd1admm,hubS.eMax,label="Hub Max",line=(:dash))

pd2admm=plot(dSoladmm.U,xlabel="",ylabel="Hub Current (kA)",legend=false,xticks=xticks, xlims=(0,hubS.K))

pd3admm=plot(dSoladmm.Tactual,label="ADMM",xlabel="",ylabel="Temp (K)", xlims=(0,hubS.K))
plot!(pd3admm,cSol.Tactual,label="Central")
plot!(pd3admm,hubS.Tmax*ones(hubS.K),label="XFRM Limit",line=(:dash,:red),xticks=xticks)

pd4admm=plot(dSoladmm.Lam,xlabel="Time",ylabel=raw"Lambda",label="ADMM",xticks=xticks)
plot!(pd4admm,cSol.Lam,label="Central")


aggU=plot(hcat(cSol.uSum,dSoladmm.uSum),label=["Central" "ADMM"],
			xlims=(0,hubS.K),xlabel="",ylabel="PEV Current (kA)")
# plot!(aggU,sum(hubS.uMax,dims=2),label="Max Current",line=(:dash,:red),xticks=xticks)

aggE=plot(hcat(sum(cSol.E,dims=2),sum(dSoladmm.E,dims=2)),label=["Central" "ADMM"], xlims=(0,hubS.K),xlabel="",ylabel="Energy (MWh)")
plot!(aggE,sum(hubS.eMax,dims=2),label="Max Energy",line=(:dash,:red),xticks=xticks)
