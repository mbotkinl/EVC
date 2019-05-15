# Run file for hub central simulation
# Micah Botkin-Levy

using Gurobi
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funHub.jl")
fname = "central_H$(H)"
if coordinated==false
	fname = fname*"_noLim"
end

if loadResults
	println("Reading in Central Sim")
	loadF=load(path*fname*".jld2")
	cSol=loadF["solution"]
else
	#initialize
	t0=hubS.t0
	e0=hubS.e0

	timeT=@elapsed cSol=hubCentral(hubS,mode,silent)
	if saveResults saveRun(path,fname,timeT,cSol) end
end

hubLabels=permutedims(["Hub $(h)" for h=1:hubS.H])
allColors=get_color_palette(:auto, plot_color(:white), H)
plotColors=allColors'

p1=plot(cSol.E,xlabel="",ylabel="Energy (MWh)",seriestype=:line,labels=hubLabels,xticks=xticks, xlims=(0,hubS.K))
plot!(hubS.eMax,label=hubLabels.*" Max",line=(:dash),seriescolor=plotColors)

eD=copy(cSol.E_depart)
eD[eD.==0].=NaN
eDmin=copy(hubS.eDepart_min)
eDmin[eDmin.==0].=NaN
eDmax=copy(hubS.eDepart_min .+ hubS.slackMax)
eDmax[eDmax.==0].=NaN
eA=copy(cSol.E_arrive)
eA[eA.==0].=NaN
sumPlot=plot(sum(cSol.E,dims=2),xlabel="",ylabel="Energy (MWh)",label="Hub Energy",seriestype=:bar,xticks=xticks)
plot!(sumPlot,sum(eD,dims=2),label="Depart Energy",seriestype=:scatter,markersize=10)
plot!(sumPlot,sum(eA,dims=2),label="Arrive Energy",seriestype=:scatter,markersize=10)
plot!(twinx(),sum(cSol.U,dims=2),label="Hub Current",seriestype=:line,seriescolor=:red,legend=false,ylabel="Current (kA)",xticks=xticks)

h=3
sumPlot=plot(cSol.E[:,h],xlabel="",ylabel="Energy (MWh)",label="Hub Energy",seriestype=:bar,xticks=xticks)
plot!(sumPlot,eD[:,h],label="Depart Energy",seriestype=:scatter,markersize=10)
plot!(sumPlot,eA[:,h],label="Arrive Energy",seriestype=:scatter,markersize=10)
plot!(twinx(),cSol.U[:,h],label="Hub Current",seriestype=:line,seriescolor=:red,legend=false,ylabel="Current (kA)",xticks=xticks)

#departPlot=plot(eD[:,h],label="Depart Actual Energy",seriestype=:bar,xlims=(200,hubS.K),yerror=eDmin[:,h] .- eD[:,h])

stT=Time(7,00)
endT=Time(10,0)
XlabelsAM=vcat(collect(stT:Dates.Second(Int(round(hubS.Ts))):endT))
xticksAM=(220:10:hubS.K,Dates.format.(XlabelsAM[1:10:61],"HH:MM"))

departPlot=plot(eD[:,h],label="Depart Actual Energy",seriestype=:bar,xlims=(220,hubS.K),ylabel="Energy (MWh)",
				xlabel="Time",xticks=xticksAM,size=(1200,800))
plot!(departPlot,eDmin[:,h],label="Depart Min Energy",seriestype=:scatter,markersize=12)
plot!(departPlot,eDmax[:,h],label="Depart Max Energy",seriestype=:scatter,markersize=12)

p2=plot(cSol.U,xlabel="",ylabel="Current (kA)",legend=false,xticks=xticks, xlims=(0,hubS.K))

p3=plot(hcat(cSol.Tactual,cSol.Tpred),label=["Actual Temp" "PWL Temp"],xlabel="",ylabel="Temp (C)", xlims=(0,hubS.K))
plot!(p3,hubS.Tmax*ones(hubS.K),label="XFRM Limit",line=(:dash,:red),xticks=xticks)

p4=plot(cSol.Lam,xlabel="Time",ylabel=raw"Lambda",legend=false,xticks=xticks, xlims=(0,hubS.K))

aggU=plot((cSol.uSum+hubS.iD_actual)/Ntf*1000,label=["Central" "ADMM"],xlabel="Time",ylabel="Total Hub Current (kA)")
plot!(aggU,hubS.iD_actual/Ntf*1000,label="Background Demand",line=(:dash),linewidth=2)

# # test to make sure all hubs meet demand
# epsilon=.01
# compMin=cSol.E_depart.>=hubS.eDepart_min
# all(compMin) # this should be true
# compDept=cSol.E_depart.-(hubS.eDepart_min.+hubS.slackMax)
# testDept=abs.(compDept).<=epsilon
# all(testDept) # this can be false
# t=findfirst(testDept.==false)
