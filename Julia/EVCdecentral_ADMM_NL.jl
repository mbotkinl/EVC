#EVC with ADMM for PWL convex relaxation
#Micah Botkin-Levy
#5/22/18

#formulation try 2
#u, sn, xt, and z are all in "x" v is the auxilliary variable corresponding to z in literature
#current constraint is coupling
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCnl.jl")

relaxString= if relaxed==true "_relax"else "" end
fname = "dADMM_NL_N$(N)"*relaxString

if loadResults
	loadF=load(path*fname*".jld2")
	evS=loadF["scenario"]
	dLognladmm=loadF["solution"]
	dCMnladmm=loadF["convMetrics"]
	convIt=loadF["convIt"]
else
	timeT=@elapsed dLognladmm,dCMnladmm,convIt=nlEVadmm(N,S,horzLen,maxIt,evS,cSolnl,relaxed,slack)
	if saveResults saveRun(path,fname,timeT, evS,dLognladmm, dCMnladmm, convIt) end
end

println("plotting....")
xPlot=zeros(horzLen+1,N)
uPlot=zeros(horzLen+1,N)
for ii= 1:N
	xPlot[:,ii]=dLognladmm.Sn[collect(ii:N:length(dLognladmm.Sn[:,convIt])),convIt]
    uPlot[:,ii]=dLognladmm.Un[collect(ii:N:length(dLognladmm.Un[:,convIt])),convIt]
end

pd1NLadmm=plot(xPlot,xlabel="Time",ylabel="PEV SOC",legend=false,xlims=(0,horzLen+1),ylims=(0,1))
if drawFig==1 savefig(pd1NLadmm,path*"J_decentralNL_ADMM_SOC.png") end

pd2NLadmm=plot(uPlot,xlabel="Time",ylabel="PEV Current (kA)",legend=false,xlims=(0,horzLen+1))
if drawFig==1 savefig(pd2NLadmm,path*"J_decentralNL_ADMM_Curr.png") end

pd3NLadmm=plot(1:horzLen+1,dLognladmm.Xt[:,convIt],label="XFRM Temp",xlims=(0,horzLen+1),xlabel="Time",ylabel="Temp (K)")
plot!(pd3NLadmm,1:horzLen+1,evS.Tmax*ones(horzLen+1),label="XFRM Limit",line=(:dash,:red))
if drawFig==1 savefig(pd3NLadmm,path*"J_decentralNL_ADMM_Temp.png") end


pd4NLadmm=plot(1:horzLen+1,hcat(cSolnl.lamCoupl,dLognladmm.Lam[:,convIt]),xlabel="Time",ylabel=raw"Lambda ($/kA)",
             xlims=(0,horzLen+1),labels=["Central" "ADMM"])
if drawFig==1 savefig(pd4NLadmm,path*"J_decentralNL_ADMM_Lam.png") end

fName="J_Central.png"

#convergence plots
halfCI=Int(floor(convIt/2))
CList=reshape([range(colorant"blue", stop=colorant"yellow",length=halfCI);
               range(colorant"yellow", stop=colorant"red",length=convIt-halfCI)], 1, convIt);


uSumPlotnladmm=plot(dLognladmm.uSum[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Current Sum (kA)",xlims=(0,horzLen+1),legend=false)
plot!(uSumPlotnladmm,1:horzLen+1,cSolnl.uSum,seriescolor=:black,linewidth=2,linealpha=0.8)

constPlotnladmm=plot(dLognladmm.couplConst[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="curr constraint diff",xlims=(0,horzLen+1),legend=false)

lamPlotnladmm=plot(dLognladmm.Lam[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Lambda",xlims=(0,horzLen+1),legend=false)
plot!(lamPlotnladmm,1:horzLen+1,cSolnl.lamCoupl,seriescolor=:black,linewidth=2,linealpha=0.8)

fPlotadmm=plot(1:convIt,dCMnladmm.obj[1:convIt,1],xlabel="Iteration",ylabel="obj function gap",xlims=(1,convIt),legend=false,yscale=:log10)
convItPlotadmm=plot(1:convIt,dCMnladmm.lamIt[1:convIt,1],xlabel="Iteration",ylabel="2-Norm Lambda Gap",xlims=(1,convIt),legend=false,yscale=:log10)
convPlotadmm=plot(1:convIt,dCMnladmm.lam[1:convIt,1],xlabel="Iteration",ylabel="central lambda gap",xlims=(1,convIt),legend=false,yscale=:log10)
constPlotadmm=plot(1:convIt,dCMnladmm.couplConst[1:convIt,1],xlabel="Iteration",ylabel="curr constraint Gap",xlims=(1,convIt),legend=false,yscale=:log10)

checkDesiredStates(dLognladmm.Sn[:,convIt],evS.Kn,evS.Snmin)
