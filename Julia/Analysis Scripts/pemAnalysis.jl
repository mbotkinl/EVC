# test different PEM parameters
using Gurobi
using LaTeXStrings
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCpem.jl")

N=evS.N
K=evS.K

runAnalysis=false

numRuns=10
packLen=2 #number of time steps
mttr=evS.Ts*packLen
setSOCs=[0.01 0.1 0.5]

# mapping ratio to probability of request for n
rLen=1000
ratio=range(0,1,length=rLen)
setP = zeros(rLen,length(setSOCs))
plotLabels=Array{String}(undef, length(setSOCs),1)
for i=1:length(setSOCs)
    mu = 1 / mttr*((ratio.-0)./(1 .- ratio))*((1 - setSOCs[i])/(setSOCs[i] - 0))
    setP[:,i] = min.(max.(1 .- exp.(-mu*evS.Ts),0),1)
    temp=@sprintf "=%s" setSOCs[i]
	# plotLabels[i]=L"\hat r_{\text{set},n}"*temp
	plotLabels[i]=L"\hat r_{set,n}"*temp
end
setPlot=plot(ratio,setP,xlabel="Ratio",ylabel="Request Probability",label=permutedims(plotLabels))
pubPlot(setPlot,thickscale=1,sizeWH=(600,400),dpi=60)
#savefig(setPlot,path*"setSOCPlot.png")

saveFile=path*"Analysis\\pemAnalysis"*".jld2"

lenCalcs=length(setSOCs)

if runAnalysis

	### changing setSOC
	Sn=zeros(K,N,lenCalcs)
	Tactual=zeros(K,1,lenCalcs)
	uSum=zeros(K,1,lenCalcs)
	requests=zeros(K,1,lenCalcs)
	for i=1:lenCalcs
		for ii=1:numRuns
		    @printf "%s: calc %g of %g....\n" Dates.format(Dates.now(),"HH:MM:SS") (i-1)*numRuns+ii lenCalcs*numRuns
		    global setSOC = setSOCs[i]
		    pemSol,Req=pemEVC(evS,slack,silent)
			Sn[:,:,i]=Sn[:,:,i].+pemSol.Sn
			Tactual[:,1,i]=Tactual[:,1,i].+pemSol.Tactual
			uSum[:,1,i]=uSum[:,1,i].+pemSol.uSum
			requests[:,1,i] = requests[:,1,i] .+ [count(Req[k,:].>0)  for k=1:evS.K]
		end
		Sn[:,:,i]=Sn[:,:,i]./numRuns
		Tactual[:,1,i]=Tactual[:,1,i]./numRuns
		uSum[:,1,i]=uSum[:,1,i]./numRuns
		requests[:,1,i]=requests[:,1,i]./numRuns
	end

	save(saveFile,"packLen",packLen,"numRuns",numRuns,"Sn", Sn, "Tactual", Tactual,"uSum",uSum,"requests",requests)
else
	loadF=load(saveFile)
	numRuns=loadF["numRuns"]
	packLen=loadF["packLen"]
	Tactual=loadF["Tactual"]
	requests=loadF["requests"]
	uSum=loadF["uSum"]
	Sn=loadF["Sn"]
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
p3pem=plot(Tactual[:,1,:],xlabel="Time",ylabel="Temp (C)",
			label="",xticks=xticks,seriescolor=plotColors)
plot!(p3pem,evS.Tmax*ones(evS.K),label="XFRM Limit",line=(:dash,:red))
plot!(p3pem,cSol.Tactual,label="",seriescolor=allColors[1])
aggUpem=plot(hcat(cSol.uSum/Ntf*1000,uSum[:,1,:]/Ntf*1000),label="",xticks=xticks, #label=["Central" permutedims(plotLabels)]
			xlabel="Time",ylabel="PEV Current (A)",seriescolor=allColors')


#plot parallelogram for n
n= rand(1:N)
parPlotn=plot(1:K,Sn[:,n,:],xlabel="Time",ylabel="PEV SoC",xlims=(0,evS.K),labels=permutedims(plotLabels),ylims=(0,1),xticks=xticks)
plot!(parPlotn,evS.s0[n]*ones(evS.Kn[n]),line=(:dash,:red),label="")
plot!(parPlotn,evS.Snmin[n]*ones(K),line=(:dash,:red),label="")
slope=evS.ηP[n]*evS.imax[n]
stop1=ceil((evS.Snmin[n]-evS.s0[n])/slope)
stop2=evS.Kn[n]-stop1
plot!(parPlotn,slope*range(1,stop1).+evS.s0[n],line=(:dash,:red),label="")
endLine=zeros(evS.Kn[n])
endLine[Int(stop2)+1:evS.Kn[n]]=slope*range(1,stop1).+evS.s0[n]
plot!(parPlotn,endLine,line=(:dash,:red),label="")

p1=parPlotn
p2=parPlotn
p3=parPlotn
parPlot3=plot(p1,p2,p3,layout=(3,1))
pubPlot(parPlot3,thickscale=0.7,sizeWH=(800,400),dpi=60)
savefig(parPlot3,path*"setSOCPlot3par.png")


#plot requests per timestep
reqPlot=plot(requests[:,1,:],xlabel="Time",ylabel="Number of Requests",xticks=xticks,seriescolor=plotColors,legend=false)

#plot parallelogram for all n
parPlot=plot(mean(Sn,dims=2)[:,1,:],xlabel="Time",ylabel="Avg PEV SoC",ylims=(0.35,1),
			label=permutedims(plotLabels),xticks=xticks,seriescolor=plotColors)
plot!(parPlot,mean(cSol.Sn,dims=2),label="Central",seriescolor=allColors[1])

# relativeSoC=zeros(evS.K,lenCalcs)
# for i =1:lenCalcs
# 	for n=1:N
# 		relativeSoC[:,i]= sum(mapLin(Sn[:,n,i],evS.s0[n],evS.Snmin[n],0,1) for n=1:N)/N
# 	end
# end

resPlot=plot(p3pem,aggUpem,reqPlot,parPlot,layout=(2,2))
pubPlot(resPlot,thickscale=1.8,sizeWH=(1400,800),dpi=100)
savefig(resPlot,path*"setSOCAnalaysisPlot.png")
