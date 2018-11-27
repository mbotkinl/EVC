#Micah Botkin-Levy
#4/10/18
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCpwl.jl")
errorString= if forecastError "_Error" else "" end
fname = "d_$(updateMethod)_N$(N)"*errorString

if loadResults
	println("Reading in Dual Decomp Sim")
	loadF=load(path*fname*".jld2")
	evS=loadF["scenario"]
	dLog=loadF["solution"]
	dCM=loadF["convMetrics"]
	convIt=loadF["convIt"]
else
	println("Running Dual Decomp Sim")
	timeT=@elapsed dLog,dCM,convIt=pwlEVdual(N,S,horzLen,maxIt,updateMethod,evS,cSol,forecastError,slack)

	s=Symbol(@sprintf("dCM_%s",updateMethod))
	v=Symbol(@sprintf("dCM"))
	@eval(($s)=($v))
	if saveResults saveRun(path,fname,timeT, evS,dLog, dCM, convIt) end
end



println("plotting....")
snPlotd=zeros(horzLen+1,N)
uPlotd=zeros(horzLen+1,N)
for ii= 1:N
	snPlotd[:,ii]=dLog.Sn[collect(ii:N:length(dLog.Sn[:,convIt])),convIt]
	uPlotd[:,ii]=dLog.Un[collect(ii:N:length(dLog.Un[:,convIt])),convIt]
end

#solution plots
pd1=plot(snPlotd,xlabel="Time",ylabel="PEV SOC",legend=false,xlims=(0,horzLen+1),ylims=(0,1))
if drawFig savefig(p1d,path*"J_"*updateMethod*"_SOC.png") end

pd2=plot(uPlotd,xlabel="Time",ylabel="PEV Current (kA)",legend=false,xlims=(0,horzLen+1))
if drawFig savefig(p2d,path*"J_"*updateMethod*"_Curr.png") end

pd3=plot(1:horzLen+1,dLog.Tactual[:,convIt],label="Actual Temp",xlims=(0,horzLen+1),xlabel="Time",ylabel="Temp (K)")
plot!(pd3,1:horzLen+1,evS.Tmax*ones(horzLen+1),label="XFRM Limit",line=(:dash,:red))
if updateMethod=="dualAscent" plot!(pd3,1:horzLen+1,dLog.Xt[:,convIt],label="PWL Temp") end
if drawFig savefig(pd3,path*"J_"*updateMethod*"_Temp.png") end

if updateMethod=="fastAscent"
	lamLabel=raw"Lambda ($/K)"
else
	lamLabel=raw"Lambda ($/kA)"
end

pd4=plot(1:horzLen+1,hcat(cSol.lamCoupl,dLog.Lam[:,convIt]),xlabel="Time",ylabel=lamLabel,xlims=(0,horzLen+1),labels=["Central" "ALADIN"])
if drawFig savefig(pd4,path*"J_"*updateMethod*"_Lam.png") end

#convergence plots
halfCI=Int(floor(convIt/2))
CList=reshape([range(colorant"blue", stop=colorant"yellow",length=halfCI);
               range(colorant"yellow", stop=colorant"red",length=convIt-halfCI)], 1, convIt);

uSumPlotd=plot(dLog.uSum[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Current Sum (kA)",xlims=(0,horzLen+1),legend=false)
plot!(uSumPlotd,1:horzLen+1,cSol.uSum,seriescolor=:black,linewidth=2,linealpha=0.8)

# uSumPlotd=plot(dLog.uSum[:,1:convIt],palette=:greens,line_z=(1:convIt)',legend=false,colorbar=:right,colorbar_title="Iteration",
#      xlabel="Time",ylabel="Current Sum (kA)",xlims=(0,horzLen+1))
# plot!(uSumPlotd,1:horzLen+1,cSol.uSum,seriescolor=:black,linewidth=2,linealpha=0.8)


zSumPlotd=plot(dLog.zSum[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Z sum",xlims=(0,horzLen+1),legend=false)
plot!(zSumPlotd,1:horzLen+1,cSol.zSum,seriescolor=:black,linewidth=2,linealpha=0.8)

lamPlotd=plot(dLog.Lam[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Lambda",xlims=(0,horzLen+1),legend=false)
plot!(lamPlotd,1:horzLen+1,cSol.lamCoupl,seriescolor=:black,linewidth=2,linealpha=0.8)

#plot(dLog.Lam[:,1:convIt],color=:RdYlBu,line_z=(1:convIt)')

constPlot2=plot(dLog.couplConst[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="curr constraint diff",xlims=(0,horzLen+1),legend=false)

#convergence metric plots
fPlot=plot(1:convIt,dCM.obj[1:convIt,1],xlabel="Iteration",ylabel="obj function gap",xlims=(1,convIt),legend=false,yscale=:log10)
convItPlot=plot(1:convIt,dCM.lamIt[1:convIt,1],xlabel="Iteration",ylabel="2-Norm Lambda Gap",xlims=(1,convIt),legend=false,yscale=:log10)
convPlot=plot(1:convIt,dCM.lam[1:convIt,1],xlabel="Iteration",ylabel="central lambda gap",xlims=(2,convIt),legend=false,yscale=:log10)
constPlot=plot(1:convIt,dCM.couplConst[1:convIt,1],xlabel="Iteration",ylabel="curr constraint Gap",xlims=(2,convIt),legend=false,yscale=:log10)

#compare central and decentral current agg
aggU=plot(1:horzLen+1,hcat(sum(uPlot[:,i] for i=1:N),sum(uPlotd[:,i] for i=1:N)),label=["Central" "Decentral"],
			xlims=(0,horzLen+1),xlabel="Time",ylabel="PEV Current (kA)")

checkDesiredStates(dLog.Sn[:,convIt],evS.Kn,evS.Snmin)


#draw(PNG(path*"aggPlot_fast.png", 13inch, 8inch), aggU)
