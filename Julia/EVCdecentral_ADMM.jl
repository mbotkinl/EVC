#EVC with ADMM for PWL convex relaxation
#Micah Botkin-Levy
#5/22/18

#formulation try 2
#u, sn, xt, and z are all in "x" v is the auxilliary variable corresponding to z in literature
#current constraint is coupling

include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCpwl.jl")
errorString= if forecastError "_Error" else "" end
fname = "dADMM_N$(N)"*errorString

if loadResults
	println("Reading in ADMM Sim")
	loadF=load(path*fname*".jld2")
	evS=loadF["scenario"]
	dLogadmm=loadF["solution"]
	dCMadmm=loadF["convMetrics"]
	convIt=loadF["convIt"]
else
	println("Running ADMM Sim")
	timeT=@elapsed dLogadmm,dCMadmm,convIt=pwlEVadmm(N,S,horzLen,maxIt,evS,cSol,forecastError,slack)
	if saveResults saveRun(path,fname,timeT, evS,dLogadmm, dCMadmm, convIt) end
end


println("plotting....")
xPlot=zeros(horzLen+1,N)
uPlotd=zeros(horzLen+1,N)
for ii= 1:N
	xPlot[:,ii]=dLogadmm.Sn[collect(ii:N:length(dLogadmm.Sn[:,convIt])),convIt]
	uPlotd[:,ii]=dLogadmm.Un[collect(ii:N:length(dLogadmm.Un[:,convIt])),convIt]
end

pd1admm=plot(xPlot,xlabel="Time",ylabel="PEV SOC",legend=false,xlims=(0,horzLen+1),ylims=(0,1))
if drawFig savefig(pd1admm,path*"J_decentral_ADMM_SOC.png") end

pd2admm=plot(uPlotd,xlabel="Time",ylabel="PEV Current (kA)",legend=false,xlims=(0,horzLen+1))
if drawFig savefig(pd2admm,path*"J_decentral_ADMM_Curr.png") end

pd3admm=plot(1:horzLen+1,hcat(dLogadmm.Tactual[:,convIt],dLogadmm.Xt[:,convIt]),label=["Actual Temp" "PWL Temp"],xlims=(0,horzLen+1),xlabel="Time",ylabel="Temp (K)")
plot!(pd3admm,1:horzLen+1,evS.Tmax*ones(horzLen+1),label="XFRM Limit",line=(:dash,:red))
if drawFig savefig(pd3admm,path*"J_decentral_ADMM_Temp.png") end

pd4admm=plot(1:horzLen+1,hcat(cSol.lamCoupl,dLogadmm.Lam[:,convIt]),xlabel="Time",ylabel=raw"Lambda ($/kA)",
             xlims=(0,horzLen+1),labels=["Central" "ADMM"])
if drawFig savefig(pd4admm,path*"J_decentral_ADMM_Lam.png") end


#convergence plots
halfCI=Int(floor(convIt/2))
CList=reshape([range(colorant"blue", stop=colorant"yellow",length=halfCI);
               range(colorant"yellow", stop=colorant"red",length=convIt-halfCI)], 1, convIt);

uSumPlotadmm=plot(dLogadmm.uSum[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Current Sum (kA)",xlims=(0,horzLen+1),legend=false)
plot!(uSumPlotadmm,1:horzLen+1,cSol.uSum,seriescolor=:black,linewidth=2,linealpha=0.8)

zSumPlotadmm=plot(dLogadmm.zSum[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Z sum",xlims=(0,horzLen+1),legend=false)
plot!(zSumPlotadmm,1:horzLen+1,cSol.zSum,seriescolor=:black,linewidth=2,linealpha=0.8)

constPlotadmm2=plot(dLogadmm.couplConst[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="curr constraint diff",xlims=(0,horzLen+1),legend=false)

lamPlotadmm=plot(dLogadmm.Lam[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Lambda",xlims=(0,horzLen+1),legend=false)
plot!(lamPlotadmm,1:horzLen+1,cSol.lamCoupl,seriescolor=:black,linewidth=2,linealpha=0.8)

fPlotadmm=plot(1:convIt,dCMadmm.obj[1:convIt,1],xlabel="Iteration",ylabel="obj function gap",xlims=(1,convIt),legend=false,yscale=:log10)
convItPlotadmm=plot(1:convIt,dCMadmm.lamIt[1:convIt,1],xlabel="Iteration",ylabel="2-Norm Lambda Gap",xlims=(1,convIt),legend=false,yscale=:log10)
convPlotadmm=plot(1:convIt,dCMadmm.lam[1:convIt,1],xlabel="Iteration",ylabel="central lambda gap",xlims=(1,convIt),legend=false,yscale=:log10)
constPlotadmm=plot(1:convIt,dCMadmm.couplConst[1:convIt,1],xlabel="Iteration",ylabel="curr constraint Gap",xlims=(1,convIt),legend=false,yscale=:log10)


#compare central and decentral current agg
aggU=plot(1:horzLen+1,hcat(sum(uPlot[:,i] for i=1:N),sum(uPlotd[:,i] for i=1:N)),label=["Central" "ADMM"],
			xlims=(0,horzLen+1),xlabel="Time",ylabel="PEV Current (kA)")


checkDesiredStates(dLogadmm.Sn[:,convIt],evS.Kn,evS.Snmin)
