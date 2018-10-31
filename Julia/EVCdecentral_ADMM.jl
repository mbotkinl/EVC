#EVC with ADMM for PWL convex relaxation
#Micah Botkin-Levy
#5/22/18

#formulation try 2
#u, sn, xt, and z are all in "x" v is the auxilliary variable corresponding to z in literature
#current constraint is coupling


filename = "dADMM_N$(N)"

if loadResults
	loadF=load(path*filename*".jld2")
	evS=loadF["scenario"]
	dLogadmm=loadF["solution"]
	dCMadmm=loadF["convMetrics"]
	convIt=loadF["convIt"]
else
	timeT=@elapsed dLogadmm,dCMadmm,convIt=pwlEVadmm(N,S,horzLen,maxIt,evS,cSol,slack)
	if saveResults saveRun(path,filename,timeT, evS,dLogadmm, dCMadmm, convIt) end
end


println("plotting....")
xPlot=zeros(horzLen+1,N)
uPlot=zeros(horzLen+1,N)
for ii= 1:N
	xPlot[:,ii]=dLogadmm.Sn[collect(ii:N:length(dLogadmm.Sn[:,convIt])),convIt]
	uPlot[:,ii]=dLogadmm.Un[collect(ii:N:length(dLogadmm.Un[:,convIt])),convIt]
end

pd1admm=plot(xPlot,xlabel="Time",ylabel="PEV SOC",legend=false,xlims=(0,horzLen+1),ylims=(0,1))
if drawFig==1 savefig(pd1admm,path*"J_decentral_ADMM_SOC.png") end

pd2admm=plot(uPlot,xlabel="Time",ylabel="PEV Current (kA)",legend=false,xlims=(0,horzLen+1))
if drawFig==1 savefig(pd2admm,path*"J_decentral_ADMM_Curr.png") end

pd3admm=plot(1:horzLen+1,hcat(dLogadmm.Tactual[:,convIt],dLogadmm.Xt[:,convIt]),label=["Actual Temp" "PWL Temp"],xlims=(0,horzLen+1),xlabel="Time",ylabel="Temp (K)")
plot!(pd3admm,1:horzLen+1,evS.Tmax*ones(horzLen+1),label="XFRM Limit",line=(:dash,:red))
if drawFig==1 savefig(pd3admm,path*"J_decentral_ADMM_Temp.png") end

pd4admm=plot(1:horzLen+1,hcat(cSol.lamCoupl,dLogadmm.Lam[:,convIt]),xlabel="Time",ylabel=raw"Lambda ($/kA)",
             xlims=(0,horzLen+1),labels=["Central" "ADMM"])
if drawFig==1 savefig(pd4admm,path*"J_decentral_ADMM_Lam.png") end
#
# lamPlotadmm=plot(dLogadmm.Lam[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
# 			layer(x=1:horzLen+1,y=cSol.lamCoupl,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
# 			Guide.xlabel("Time"), Guide.ylabel("Lambda"),Guide.ColorKey(title="Iteration"),
# 			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
# 			minor_label_font_size=26pt,key_label_font_size=26pt))
# constPlotadmm2=plot(dLogadmm.couplConst[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
# 			Guide.xlabel("Time"), Guide.ylabel("curr constraint diff"),Guide.ColorKey(title="Iteration"),
# 			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
# 			minor_label_font_size=26pt,key_label_font_size=26pt))
# zSumPlotadmm=plot(dLogadmm.zSum[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
# 			layer(x=1:horzLen+1,y=cSol.zSum,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
# 			Guide.xlabel("Time"), Guide.ylabel("Z sum"),Guide.ColorKey(title="Iteration"),
# 			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
# 			minor_label_font_size=26pt,key_label_font_size=26pt))
# uSumPlotadmm=plot(dLogadmm.uSum[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
# 			layer(x=1:horzLen+1,y=cSol.uSum,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
# 			Guide.xlabel("Time"), Guide.ylabel("U sum"),Guide.ColorKey(title="Iteration"),
# 			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
# 			minor_label_font_size=26pt,key_label_font_size=26pt))
# if drawFig==1 draw(PNG(path*"J_ADMM_LamConv.png", 36inch, 12inch), lamPlotadmm) end



fPlotadmm=plot(1:convIt,dCMadmm.obj[1:convIt,1],xlabel="Iteration",ylabel="obj function gap",xlims=(1,convIt),legend=false,yscale=:log10)
convItPlotadmm=plot(1:convIt,dCMadmm.lamIt[1:convIt,1],xlabel="Iteration",ylabel="2-Norm Lambda Gap",xlims=(1,convIt),legend=false,yscale=:log10)
convPlotadmm=plot(1:convIt,dCMadmm.lam[1:convIt,1],xlabel="Iteration",ylabel="central lambda gap",xlims=(1,convIt),legend=false,yscale=:log10)
constPlotadmm=plot(1:convIt,dCMadmm.couplConst[1:convIt,1],xlabel="Iteration",ylabel="curr constraint Gap",xlims=(1,convIt),legend=false,yscale=:log10)

#checkDesiredStates(dLogadmm.Sn,evS.Kn,evS.Snmin)
