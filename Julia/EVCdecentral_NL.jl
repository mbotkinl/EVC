#Micah Botkin-Levy
#4/10/18

relaxString= if relaxed==true "_relax"else "" end
filename = "d_$(updateMethod)_NL_N$(N)"*relaxString

if loadResults
	loadF=load(path*filename*".jld2")
	evS=loadF["scenario"]
	dLognl=loadF["solution"]
	dCMnl=loadF["convMetrics"]
	convIt=loadF["convIt"]
else
	timeT=@elapsed dLognl,dCMnl,convIt=nlEVdual(N,S,horzLen,maxIt,updateMethod,evS,cSolnl,relaxed,slack)
	s=Symbol(@sprintf("dCMnl_%s",updateMethod))
	v=Symbol(@sprintf("dCMnl"))
	@eval(($s)=($v))
	if saveResults saveRun(path,filename,timeT, evS,dLognl, dCMnl, convIt) end
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

#fName="J_Decentral_notfast.png"
#fName="J_Decentral_fast.png"
#draw(PNG(path*fName, 13inch, 14inch), vstack(pd1,pd2,pd3,pd4))
#
# uSumPlotNL=plot(dLognl.uSum[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
# 			layer(x=1:horzLen+1,y=cSolnl.uSum,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
# 			Guide.xlabel("Time"), Guide.ylabel("U sum"),Guide.ColorKey(title="Iteration"),
# 			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
# 			minor_label_font_size=26pt,key_label_font_size=26pt))
# iPlotNL=plot(dLognl.Itotal[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
# 			layer(x=1:horzLen+1,y=cSolnl.Itotal,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
# 			Guide.xlabel("Time"), Guide.ylabel("Z sum"),Guide.ColorKey(title="Iteration"),
# 			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
# 			minor_label_font_size=26pt,key_label_font_size=26pt))
# lamPlotNL=plot(dLognl.Lam[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
# 			layer(x=1:horzLen+1,y=cSolnl.lamCoupl,Geom.line,Theme(default_color=colorant"black",line_width=4pt)),
# 			Guide.xlabel("Time"), Guide.ylabel("Lambda"),Guide.ColorKey(title="Iteration"),
# 			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
# 			minor_label_font_size=26pt,key_label_font_size=26pt))
# if drawFig==1 draw(PNG(path*"J_"*updateMethod*"_LamConv.png", 36inch, 12inch), lamPlot) end


fPlotNL=plot(1:convIt,dCMnl.obj[1:convIt,1],xlabel="Iteration",ylabel="obj function gap",xlims=(1,convIt),legend=false,yscale=:log10)
convItPlotNL=plot(1:convIt,dCMnl.lamIt[1:convIt,1],xlabel="Iteration",ylabel="2-Norm Lambda Gap",xlims=(1,convIt),legend=false,yscale=:log10)
convPlotNL=plot(1:convIt,dCMnl.lam[1:convIt,1],xlabel="Iteration",ylabel="central lambda gap",xlims=(1,convIt),legend=false,yscale=:log10)
constPlotNL=plot(1:convIt,dCMnl.couplConst[1:convIt,1],xlabel="Iteration",ylabel="curr constraint Gap",xlims=(1,convIt),legend=false,yscale=:log10)

checkDesiredStates(dLognl.Sn,evS.Kn,evS.Snmin)
