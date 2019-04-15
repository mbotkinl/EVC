# test temperature augmented objective expressions

using Gurobi
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCpwl.jl")

tempAugment=true
slack=true
silent=true
#ψs=-range(maxPsi/numPsi,maxPsi,length=numPsi)
# setRange=range(2, stop=6, length=6)
# ψs=- vcat(0,round.(exp10.(setRange)))
ψs=[0 -1e6 -5e6 -1e7]
#ψs=[0 -1e4 -1e5 -1e6 -1e7]
#ψs=[0 -1e1 -1e3 -1e5 -1e6 -1e7]

lenCals=length(ψs)
uSum=zeros(evS.K,lenCals)
Tactual=zeros(evS.K,lenCals)
lamCoupl=zeros(evS.K,lenCals)
plotLabels=Array{String}(undef, lenCals,1)
for i=1:lenCals
    @printf "%s: calc %g of %g....\n" Dates.format(Dates.now(),"HH:MM:SS") i lenCals
    global t0=evS.t0
    global s0=evS.s0
    global ψ=ψs[i]
    timeT=@elapsed cSol, cSave=pwlEVcentral(evS,slack,silent)
    @printf "SOC satisfied %s \n" checkDesiredStates(cSol.Sn,evS.Kn,evS.Snmin)
    uSum[:,i]=cSol.uSum
    lamCoupl[:,i]=cSol.lamCoupl
    Tactual[:,i]=cSol.Tactual
    plotLabels[i]=@sprintf "ψ=%s" ψ
end



allColors=get_color_palette(:auto, plot_color(:white), lenCals+1)
plotColors=allColors[1:lenCals]'

stT1=Time(20,0)
endT1=Time(23,59)
stT2=Time(0,0)
endT2=Time(10,0)
Xlabels=vcat(collect(stT1:Dates.Second(round(evS.Ts)):endT1),collect(stT2:Dates.Second(round(evS.Ts)):endT2))
#Xlabels=vcat(collect(stT1:Dates.Minute(3):endT1),collect(stT2:Dates.Minute(3):endT2))
xticks=(1:40:evS.K,Dates.format.(Xlabels[1:40:evS.K],"HH:MM"))

iD=evS.iD_actual[1:evS.K]#.+10
loadPlot=plot(1:evS.K,iD,label="Background Demand",line=(:dash),ylabel="Total Load (kA)",xlims=(0,evS.K),xticks=xticks)
plot!(loadPlot,uSum.+iD,label="",seriescolor=plotColors)

tempPlot=plot(1:evS.K,evS.Tmax*ones(evS.K)*1000,label="XFRM Limit",line=(:dash,:red),ylabel="Temperature (K)",xlims=(0,evS.K),xticks=xticks)
plot!(tempPlot,Tactual*1000,label="",seriescolor=plotColors)

lamPlot =plot(1:evS.K,lamCoupl,ylabel="Lambda",label=permutedims(plotLabels),xlims=(0,evS.K),xticks=xticks,seriescolor=plotColors)


resPlot=plot(tempPlot,loadPlot,lamPlot,layout=(3,1))
pubPlot(resPlot,thickscale=0.7,sizeWH=(600,400),dpi=60)
savefig(resPlot,path*"tempObjPlot.png")
