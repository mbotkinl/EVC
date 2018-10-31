#Micah Botkin-Levy
#4/10/18

timeT=@elapsed dLog,dCM,convIt=pwlEVdual(N,S,horzLen,maxIt,updateMethod,evS,cSol,slack)

s=Symbol(@sprintf("dCM_%s",updateMethod))
v=Symbol(@sprintf("dCM"))
@eval(($s)=($v))

fname = "d_$(updateMethod)_N$(N)"
# save
if saveResults==1 saveRun(path,fname,timeT, evS,dLog, dCM, convIt) end
# load
if loadResults==1
	loadF=load(path*fname*".jld2")
	evS=loadF["scenario"]
	dLog=loadF["solution"]
	dCM=loadF["convMetrics"]
	convIt=loadF["convIt"]
end


println("plotting....")
snPlotd=zeros(horzLen+1,N)
uPlotd=zeros(horzLen+1,N)
for ii= 1:N
	snPlotd[:,ii]=dLog.Sn[collect(ii:N:length(dLog.Sn[:,convIt])),convIt]
	uPlotd[:,ii]=dLog.Un[collect(ii:N:length(dLog.Un[:,convIt])),convIt]
end

pd1=plot(snPlotd,xlabel="Time",ylabel="PEV SOC",legend=false,xlims=(0,horzLen+1),ylims=(0,1))
if drawFig==1 savefig(p1d,path*"J_"*updateMethod*"_SOC.png") end

pd2=plot(uPlotd,xlabel="Time",ylabel="PEV Current (kA)",legend=false,xlims=(0,horzLen+1))
if drawFig==1 savefig(p2d,path*"J_"*updateMethod*"_Curr.png") end

pd3=plot(1:horzLen+1,dLog.Tactual[:,convIt],label="Actual Temp",xlims=(0,horzLen+1),xlabel="Time",ylabel="Temp (K)")
plot!(pd3,1:horzLen+1,evS.Tmax*ones(horzLen+1),label="XFRM Limit",line=(:dash,:red))
if updateMethod=="dualAscent" plot!(pd3,1:horzLen+1,dLog.Xt[:,convIt],label="PWL Temp") end
if drawFig==1 savefig(pd3,path*"J_"*updateMethod*"_Temp.png") end

if updateMethod=="fastAscent"
	lamLabel=raw"Lambda ($/K)"
else
	lamLabel=raw"Lambda ($/kA)"
end

pd4=plot(1:horzLen+1,dLog.Lam[:,convIt],xlabel="Time",ylabel=lamLabel,xlims=(0,horzLen+1),legend=false)
if drawFig==1 savefig(pd4,path*"J_"*updateMethod*"_Lam.png") end


#fName="J_Decentral_notfast.png"
#fName="J_Decentral_fast.png"
#draw(PNG(path*fName, 13inch, 14inch), vstack(pd1,pd2,pd3,pd4))

#
# A = [i for i=1:100, j=1:100]
# heatmap(A, c=ColorGradient([:red,:yellow,:blue]))
#
# C(g::ColorGradient) = RGB[g[z] for z=range(1,step=1,length=10)]
# g = :inferno
# cgrad(g) |> C
#
# colors=[:red]
# plot(1:horzLen+1,dLog.uSum[:,1:4], linecolor=cgrad(:inferno))
# plot(1:horzLen+1,dLog.uSum[:,1:10],	color=cgrad(:grays))
# plot(1:horzLen+1,dLog.uSum[:,1:10],	zcolor=1:10)
#
# plot(1:horzLen+1,dLog.uSum[:,1:10], palette=:blues)
#

#
# uSumPlotd=plot(dLog.uSum[:,1:convIt], line_z=1:convIt,xlabel="Time",ylabel="Current Sum",xlims=(0,horzLen+1),legend=false)
#
# if drawFig==1 savefig(pd4,path*"J_"*updateMethod*"_Lam.png") end
#
# uSumPlotd=plot(dLog.uSum[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
# 			layer(x=1:horzLen+1,y=cSol.uSum,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
# 			Guide.xlabel("Time"), Guide.ylabel("U sum"),Guide.ColorKey(title="Iteration"),
# 			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
# 			minor_label_font_size=26pt,key_label_font_size=26pt))
# zSumPlotd=plot(dLog.zSum[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
# 			layer(x=1:horzLen+1,y=cSol.zSum,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
# 			Guide.xlabel("Time"), Guide.ylabel("Z sum"),Guide.ColorKey(title="Iteration"),
# 			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
# 			minor_label_font_size=26pt,key_label_font_size=26pt))
# lamPlot=plot(dLog.Lam[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
# 			layer(x=1:horzLen+1,y=cSol.lamCoupl,Geom.line,Theme(default_color=colorant"black",line_width=4pt)),
# 			Guide.xlabel("Time"), Guide.ylabel("Lambda"),Guide.ColorKey(title="Iteration"),
# 			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
# 			minor_label_font_size=26pt,key_label_font_size=26pt))
# if drawFig==1 draw(PNG(path*"J_"*updateMethod*"_LamConv.png", 36inch, 12inch), lamPlot) end

fPlot=plot(1:horzLen+1,dCM.obj[1:convIt-1,1],xlabel="Iteration",xlims=(0,convIt),legend=false)
yaxis!(fPlot,"obj function gap",:log10)
convItPlot=plot(1:horzLen+1,dCM.lamIt[1:convIt-1,1],xlabel="Iteration",xlims=(0,convIt),legend=false)
yaxis!(convItPlot,"2-Norm Lambda Gap",:log10)
convPlot=plot(1:horzLen+1,dCM.lam[1:convIt-1,1],xlabel="Iteration",xlims=(0,convIt),legend=false)
yaxis!(convPlot,"central lambda gap",:log10)

#checkDesiredStates(dLog.Sn,evS.Kn,evS.Snmin)


#compare central and decentral current agg
aggU=plot(1:horzLen+1,hcat(sum(uPlot[:,i] for i=1:N),sum(uPlotd[:,i] for i=1:N)),label=["Central" "Decentral"],
			xlims=(0,horzLen+1),xlabel="Time",ylabel="PEV Current (kA)")

#draw(PNG(path*"aggPlot_fast.png", 13inch, 8inch), aggU)
