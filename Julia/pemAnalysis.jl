# test different PEM parameters
using Gurobi
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCpem.jl")

N=evS.N
K=evS.K

numRuns=10
packLen=2 #number of time steps
mttr=evS.Ts*packLen
#setSOCs=[0.01 0.25 0.5 0.75 .99]
setSOCs=[0.35 0.5 0.65]

# mapping ratio to probability of request for n
rLen=1000
ratio=range(0,1,length=rLen)
setP = zeros(rLen,length(setSOCs))
plotLabels=Array{String}(undef, length(setSOCs),1)
for i=1:length(setSOCs)
    mu = 1 / mttr*((ratio.-0)./(1 .- ratio))*((1 - setSOCs[i])/(setSOCs[i] - 0))
    setP[:,i] = min.(max.(1 .- exp.(-mu*evS.Ts),0),1)
    plotLabels[i]=@sprintf "setSOC=%s" setSOCs[i]
end
setPlot=plot(ratio,setP,xlabel="Ratio",ylabel="Req Probability",label=permutedims(plotLabels))
pubPlot(setPlot,thickscale=1,sizeWH=(600,400),dpi=60)
# savefig(setPlot,path*"setPlot.png")


### changing setSOC
lenCalcs=length(setSOCs)
Sn=zeros(K,N,lenCalcs)
Tactual=zeros(K,1,lenCalcs)
uSum=zeros(K,1,lenCalcs)
plotLabels=Array{String}(undef, lenCalcs,1)
for i=1:lenCalcs
	for ii=1:numRuns
	    @printf "%s: calc %g of %g....\n" Dates.format(Dates.now(),"HH:MM:SS") (i-1)*numRuns+ii lenCalcs*numRuns
	    global setSOC = setSOCs[i]
	    pemSol=pemEVC(evS,slack,silent)
		Sn[:,:,i]=Sn[:,:,i].+pemSol.Sn
		Tactual[:,1,i]=Tactual[:,1,i].+pemSol.Tactual
		uSum[:,1,i]=uSum[:,1,i].+pemSol.uSum
	end
	Sn[:,:,i]=Sn[:,:,i]./numRuns
	Tactual[:,1,i]=Tactual[:,1,i]./numRuns
	uSum[:,1,i]=uSum[:,1,i]./numRuns
	plotLabels[i]=@sprintf "setSOC=%s" setSOC
end

if true
	saveFile=path*"pemAnalysis"*".jld2"
	save(saveFile,"packLen",packLen,"numRuns",numRuns,"Sn", Sn, "Tactual", Tactual,"uSum",uSum)
end


allColors=get_color_palette(:auto, plot_color(:white), lenCalcs+1)
plotColors=allColors[2:lenCalcs+1]'

stT1=Time(20,0)
endT1=Time(23,59)
stT2=Time(0,0)
endT2=Time(10,0)
Xlabels=vcat(collect(stT1:Dates.Second(round(evS.Ts)):endT1),collect(stT2:Dates.Second(round(evS.Ts)):endT2))
#Xlabels=vcat(collect(stT1:Dates.Minute(3):endT1),collect(stT2:Dates.Minute(3):endT2))
xticks=(1:40:evS.K,Dates.format.(Xlabels[1:40:evS.K],"HH:MM"))

# normal plots
p3pem=plot(Tactual[:,1,:]*1000,xlims=(0,evS.K),xlabel="Time",ylabel="Temp (K)",
			label="",xticks=xticks,seriescolor=plotColors)
plot!(p3pem,evS.Tmax*ones(evS.K)*1000,label="XFRM Limit",line=(:dash,:red))
plot!(p3pem,cSol.Tactual*1000,label="",seriescolor=allColors[1])
aggUpem=plot(hcat(cSol.uSum,uSum[:,1,:]),label="",xticks=xticks, #label=["Central" permutedims(plotLabels)]
			xlims=(0,evS.K),xlabel="Time",ylabel="PEV Current (kA)",seriescolor=allColors')


#plot parallelogram for n
n=2
parPlotn=plot(1:K,Sn[:,n,:],xlabel="Time",ylabel="PEV SOC",legend=false,xlims=(0,evS.K),ylims=(0,1),label=permutedims(plotLabels))
plot!(parPlotn,evS.s0[n]*ones(evS.Kn[n]),line=(:dash,:red))
plot!(parPlotn,evS.Snmin[n]*ones(K),line=(:dash,:red))
slope=evS.ηP[n]*evS.imax[n]
stop1=ceil((evS.Snmin[n]-evS.s0[n])/slope)
stop2=evS.Kn[n]-stop1
plot!(parPlotn,slope*range(1,stop1).+evS.s0[n],line=(:dash,:red))
endLine=zeros(evS.Kn[n])
endLine[Int(stop2)+1:evS.Kn[n]]=slope*range(1,stop1).+evS.s0[n]
plot!(parPlotn,endLine,line=(:dash,:red))


#plot parallelogram for all n
parPlot=plot(1:K,mean(Sn,dims=2)[:,1,:],xlabel="Time",ylabel="Avg PEV SOC",xlims=(0,evS.K),ylims=(0,1),
			label=permutedims(plotLabels),xticks=xticks,seriescolor=plotColors)
plot!(parPlot,mean(cSol.Sn,dims=2),label="Central",seriescolor=allColors[1])

resPlot=plot(p3pem,aggUpem,parPlot,layout=(3,1))
pubPlot(resPlot,thickscale=0.7,sizeWH=(600,400),dpi=60)
savefig(resPlot,path*"setSOCAnalaysisPlot.png")
