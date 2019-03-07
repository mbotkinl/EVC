using Gurobi
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funHub.jl")
fname = "central_H$(H)"

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

sumPlot=plot(sum(cSol.E,dims=2),xlabel="",ylabel="Energy (MWh)",label="Hub Energy",seriestype=:bar,xticks=xticks)
plot!(sumPlot,sum(cSol.E_depart,dims=2),label="Depart Energy",seriestype=:scatter,markersize=10)
plot!(sumPlot,sum(cSol.E_arrive,dims=2),label="Arrive Energy",seriestype=:scatter,markersize=10)
plot!(twinx(),sum(cSol.U,dims=2),label="Hub Current",seriestype=:line,seriescolor=:red,legend=false,ylabel="Current (kA)",xticks=xticks)

p2=plot(cSol.U,xlabel="",ylabel="Current (kA)",legend=false,xticks=xticks, xlims=(0,hubS.K))

p3=plot(hcat(cSol.Tactual,cSol.Tpred),label=["Actual Temp" "PWL Temp"],xlabel="",ylabel="Temp (K)", xlims=(0,hubS.K))
plot!(p3,hubS.Tmax*ones(hubS.K),label="XFRM Limit",line=(:dash,:red),xticks=xticks)

p4=plot(cSol.Lam/1000,xlabel="Time",ylabel=raw"Lambda ($/A)",legend=false,xticks=xticks, xlims=(0,hubS.K))


# epsilon=.01
# compMin=cSol.E_depart.>=hubS.eDepart_min
# all(compMin)
# compDept=cSol.E_depart.-(hubS.eDepart_min.+hubS.slackMax)
# testDept=abs.(compDept).<=epsilon
# all(testDept)
# findfirst(testDept.==false)


#nolim plots
# currComp=plot(hcat(cSol.uSum,noLim.uSum),label=["Central" "Uncoordinated"],xlabel="",ylabel="Current (kA)", xlims=(0,hubS.K),xticks=xticks)
#
# tempComp=plot(hcat(cSol.Tactual*1000,noLim.Tactual*1000),label=["Central" "Uncoordinated"],xlabel="",ylabel="Temp (K)", xlims=(0,hubS.K),xticks=xticks)
# plot!(tempComp,hubS.Tmax*ones(hubS.K)*1000,label="XFRM Limit",line=(:dash,:red))
#
# compP=plot(currComp,tempComp,p4,layout=(3,1))
# pubPlot(compP,thickscale=0.8,sizeWH=(1000,600),dpi=100)
# savefig(compP,path*"centralPlot1.png")


# h1=plot(p1,p2,p3,p4,layout=(4,1))
# lowRes=true
# if lowRes
#     pubPlot(h1,thickscale=0.4,sizeWH=(400,300),dpi=40)
# else
#     pubPlot(h1,thickscale=0.8,sizeWH=(1000,600),dpi=100)
# end
# if saveF savefig(h1,path*"hubPlot1.png") end
#
# h2=plot(sumPlot,p3,p4,layout=(3,1))
# lowRes=true
# if lowRes
#     pubPlot(h2,thickscale=0.4,sizeWH=(400,300),dpi=40)
# else
#     pubPlot(h2,thickscale=0.8,sizeWH=(800,600),dpi=100)
# end
# if saveF savefig(h2,path*"hubPlot2.png") end
