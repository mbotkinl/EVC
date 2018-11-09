#Micah Botkin-Levy
#4/10/18
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCnl.jl")

relaxString= if relaxed==true "_relax"else "" end
fname = "d_$(updateMethod)_NL_N$(N)"*relaxString

if loadResults
	loadF=load(path*fname*".jld2")
	evS=loadF["scenario"]
	dLognl=loadF["solution"]
	dCMnl=loadF["convMetrics"]
	convIt=loadF["convIt"]
else
	timeT=@elapsed dLognl,dCMnl,convIt=nlEVdual(N,S,horzLen,maxIt,updateMethod,evS,cSolnl,relaxed,slack)
	s=Symbol(@sprintf("dCMnl_%s",updateMethod))
	v=Symbol(@sprintf("dCMnl"))
	@eval(($s)=($v))
	if saveResults saveRun(path,fname,timeT, evS,dLognl, dCMnl, convIt) end
end


println("plotting....")
xPlotd=zeros(horzLen+1,N)
uPlotd=zeros(horzLen+1,N)
for ii= 1:N
	xPlotd[:,ii]=dLognl.Sn[collect(ii:N:length(dLognl.Sn[:,convIt])),convIt]
	uPlotd[:,ii]=dLognl.Un[collect(ii:N:length(dLognl.Un[:,convIt])),convIt]
end
pd1nl=plot(xPlotd,xlabel="Time",ylabel="PEV SOC",legend=false,xlims=(0,horzLen+1),ylims=(0,1))
if drawFig==1 savefig(pd1nl,path*"J_NL_"*updateMethod*"_SOC.png") end

pd2nl=plot(uPlotd,xlabel="Time",ylabel="PEV Current (kA)",legend=false,xlims=(0,horzLen+1))
if drawFig==1 savefig(pd2nl,path*"J_NL_"*updateMethod*"_Curr.png") end

pd3nl=plot(1:horzLen+1,evS.Tmax*ones(horzLen+1),label="XFRM Limit",line=(:dash,:red),xlims=(0,horzLen+1),xlabel="Time",ylabel="Temp (K)")
if updateMethod=="dualAscent" plot!(pd3nl,1:horzLen+1,dLognl.Xt[:,convIt],label="XFRM Temp") end
if drawFig==1 savefig(pd3nl,path*"J_NL_"*updateMethod*"_Temp.png") end


if updateMethod=="fastAscent"
	lamLabel=raw"Lambda ($/K)"
else
	lamLabel=raw"Lambda ($/kA)"
end
pd4nl=plot(1:horzLen+1,dLognl.Lam[:,convIt],xlabel="Time",ylabel=lamLabel,xlims=(0,horzLen+1),legend=false)
if drawFig==1 savefig(pd4nl,path*"J_NL_"*updateMethod*"_Lam.png") end

#convergence plots
halfCI=Int(floor(convIt/2))
CList=reshape([range(colorant"blue", stop=colorant"yellow",length=halfCI);
               range(colorant"yellow", stop=colorant"red",length=convIt-halfCI)], 1, convIt);


uSumPlotnl=plot(dLognl.uSum[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Current Sum (kA)",xlims=(0,horzLen+1),legend=false)
plot!(uSumPlotnl,1:horzLen+1,cSolnl.uSum,seriescolor=:black,linewidth=2,linealpha=0.8)

constPlotnl2=plot(dLognl.couplConst[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="curr constraint diff",xlims=(0,horzLen+1),legend=false)

lamPlotnl=plot(dLognl.Lam[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Lambda",xlims=(0,horzLen+1),legend=false)
plot!(lamPlotnl,1:horzLen+1,cSolnl.lamCoupl,seriescolor=:black,linewidth=2,linealpha=0.8)


fPlotNL=plot(1:convIt,dCMnl.obj[1:convIt,1],xlabel="Iteration",ylabel="obj function gap",xlims=(1,convIt),legend=false,yscale=:log10)
convItPlotNL=plot(1:convIt,dCMnl.lamIt[1:convIt,1],xlabel="Iteration",ylabel="2-Norm Lambda Gap",xlims=(1,convIt),legend=false,yscale=:log10)
convPlotNL=plot(1:convIt,dCMnl.lam[1:convIt,1],xlabel="Iteration",ylabel="central lambda gap",xlims=(1,convIt),legend=false,yscale=:log10)
constPlotNL=plot(1:convIt,dCMnl.couplConst[1:convIt,1],xlabel="Iteration",ylabel="curr constraint Gap",xlims=(1,convIt),legend=false,yscale=:log10)

checkDesiredStates(dLognl.Sn[:,convIt],evS.Kn,evS.Snmin)
