# Run file for hub dual decompotions simulation
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funHub.jl")
fname = "dDual_H$(H)"

if loadResults
	println("Reading in Dual Hub Sim")
	loadF=load(path*fname*".jld2")
	dSol=loadF["solution"]
else
	t0=hubS.t0
	e0=hubS.e0
	prevLam=8e3*ones(hubS.K1+1,1)
	roundSigFigs=16

	println("Running Dual Hub Sim")
	timeT=@elapsed dSol=hubDual(maxIt,hubS,cSol,mode,silent)
	if saveResults saveRun(path,fname,timeT, dSol) end
end

hubLabels=permutedims(["Hub $(h)" for h=1:hubS.H])
allColors=get_color_palette(:auto, plot_color(:white), H)
plotColors=allColors'

println("plotting....")
pd1alad=plot(dSol.E,xlabel="Time",ylabel="PEV SOC",legend=false,xlims=(0,hubS.K))
plot!(hubS.eMax,label=hubLabels.*" Max",line=(:dash),seriescolor=plotColors)

pd2alad=plot(dSol.U,xlabel="Time",ylabel="PEV Current (kA)",legend=false,xlims=(0,hubS.K))

pd3alad=plot(hcat(dSol.Tactual[:,1],dSol.Tpred[:,1]),label=["Actual Temp" "PWL Temp"],xlims=(0,hubS.K),xlabel="Time",ylabel="Temp (K)")
plot!(pd3alad,1:hubS.K,hubS.Tmax*ones(hubS.K),label="XFRM Limit",line=(:dash,:red))

pd4alad=plot(hcat(cSol.Lam,dSol.Lam),xlabel="Time",ylabel=raw"Lambda ($/kA)",
             xlims=(0,hubS.K),labels=["Central" "Dual"])

aggU=plot(hcat(cSol.uSum,dSol.uSum),label=["Central" "Dual"],
			xlims=(0,hubS.K),xlabel="Time",ylabel="PEV Current (kA)")

aggZ=plot(hcat(cSol.zSum,dSol.zSum),label=["Central" "Dual"],
			xlims=(0,hubS.K),xlabel="Time",ylabel="PEV Current (kA)")
