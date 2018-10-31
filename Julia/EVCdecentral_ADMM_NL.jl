#EVC with ADMM for PWL convex relaxation
#Micah Botkin-Levy
#5/22/18

#formulation try 2
#u, sn, xt, and z are all in "x" v is the auxilliary variable corresponding to z in literature
#current constraint is coupling

relaxString= if relaxed==true "_relax"else "" end
filename = "dADMM_NL_N$(N)"*relaxString

if loadResults
	loadF=load(path*filename*".jld2")
	evS=loadF["scenario"]
	dLognladmm=loadF["solution"]
	dCMnladmm=loadF["convMetrics"]
	convIt=loadF["convIt"]
else
	timeT=@elapsed dLognladmm,dCMnladmm,convIt=nlEVadmm(N,S,horzLen,maxIt,evS,cSolnl,relaxed,slack)
	if saveResults saveRun(path,filename,timeT, evS,dLognladmm, dCMnladmm, convIt) end
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


# lamPlotNLadmm=plot(dLognladmm.Lam[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
#             layer(x=1:horzLen+1,y=cSolnl.lamCoupl,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
# 			Guide.xlabel("Time"), Guide.ylabel("Lambda",orientation=:vertical),Guide.ColorKey(title="Iteration"),
# 			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
# 			minor_label_font_size=26pt,key_label_font_size=26pt))
# constPlotNLadmm2=plot(dLognladmm.couplConst[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
# 			Guide.xlabel("Time"), Guide.ylabel("curr constraint diff",orientation=:vertical),Guide.ColorKey(title="Iteration"),
# 			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
# 			minor_label_font_size=26pt,key_label_font_size=26pt))
# if drawFig==1 draw(PNG(path*"J_ADMM_LamConv.png", 36inch, 12inch), lamPlotadmm) end



fPlotadmm=plot(1:horzLen+1,dCMnladmm.obj[1:convIt-1,1],xlabel="Iteration",ylabel="obj function gap",xlims=(0,convIt),legend=false,yscale=:log10)
convItPlotadmm=plot(1:horzLen+1,dCMnladmm.lamIt[1:convIt-1,1],xlabel="Iteration",ylabel="2-Norm Lambda Gap",xlims=(0,convIt),legend=false,yscale=:log10)
convPlotadmm=plot(1:horzLen+1,dCMnladmm.lam[1:convIt-1,1],xlabel="Iteration",ylabel="central lambda gap",xlims=(0,convIt),legend=false,yscale=:log10)
constPlotadmm=plot(1:horzLen+1,dCMnladmm.couplConst[1:convIt-1,1],xlabel="Iteration",ylabel="curr constraint Gap",xlims=(0,convIt),legend=false,yscale=:log10)

#checkDesiredStates(dLognladmm.Sn,evS.Kn,evS.Snmin)
