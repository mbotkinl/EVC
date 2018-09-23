#Micah Botkin-Levy
#4/10/18

tic()
dLognl,dCMnl,convIt=nlEVdual(N,S,horzLen,maxIt,updateMethod,evS,cSolnl)
dLognl.timeT=toc()

s=Symbol(@sprintf("dCMnl_%s",updateMethod))
v=Symbol(@sprintf("dCMnl"))
@eval(($s)=($v))


filename = "d_$updateMethod_NL_N$(N)"
# save
if saveResults==1 saveRun(path,filename,timeT, evS,dLognl, dCMnl, convIt) end
# load
if loadResults==1
	loadF=JLD.load(path*filename*".jld")
	evS=loadF["scenario"]
	dLognl=loadF["solution"]
	convIt=loadF["convIt"]
end


println("plotting....")
xPlot=zeros(horzLen+1,N)
uPlotd=zeros(horzLen+1,N)
for ii= 1:N
	xPlot[:,ii]=dLognl.Sn[collect(ii:N:length(dLognl.Sn[:,convIt])),convIt]
	uPlotd[:,ii]=dLognl.Un[collect(ii:N:length(dLognl.Un[:,convIt])),convIt]
end

pd1nl=plot(xPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV SOC"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=24pt,line_width=3pt,
		minor_label_font_size=20pt,key_label_font_size=20pt))
if drawFig==1 draw(PNG(path*"J_"*updateMethod*"_SOC.png", 24inch, 12inch), pd1) end

pd2nl=plot(uPlotd,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV Current (kA)"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1,ymin=0),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=24pt,line_width=3pt,
		minor_label_font_size=20pt,key_label_font_size=20pt))
if drawFig==1 draw(PNG(path*"J_"*updateMethod*"_Curr.png", 24inch, 12inch), pd2) end

pd3nl=plot(layer(yintercept=[evS.Tmax],Geom.hline(color=["red"],style=:dot),Theme(line_width=3pt)),
		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=24pt,
		minor_label_font_size=20pt,key_label_font_size=20pt))
if updateMethod=="dualAscent"
	push!(pd3nl,layer(x=1:horzLen+1,y=dLognl.Xt[:,convIt],Geom.line,Theme(default_color=colorant"blue",line_width=3pt)))
	push!(pd3nl,Guide.manual_color_key("", ["PWL Temp", "Actual Temp"], ["blue", "green"]))
	push!(pd3nl,Theme(key_position = :top,background_color=colorant"white",major_label_font_size=24pt,line_width=3pt,
	minor_label_font_size=20pt,key_label_font_size=20pt))
end
if drawFig==1 draw(PNG(path*"J_"*updateMethod*"_Temp.png", 24inch, 12inch), pd3) end

if updateMethod=="fastAscent"
	lamLabel=raw"Lambda ($/K)"
else
	lamLabel=raw"Lambda ($/kA)"
end
pd4nl=plot(layer(x=1:horzLen+1,y=dLognl.Lam[:,convIt],Geom.line,Theme(line_width=3pt)),
		Guide.xlabel("Time"), Guide.ylabel(lamLabel),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=24pt,
		minor_label_font_size=20pt,key_label_font_size=20pt))
if drawFig==1 draw(PNG(path*"J_"*updateMethod*"_Lam.png", 24inch, 12inch), pd4) end

#fName="J_Decentral_notfast.png"
#fName="J_Decentral_fast.png"
#draw(PNG(path*fName, 13inch, 14inch), vstack(pd1,pd2,pd3,pd4))

uSumPlotNL=plot(dLognl.uSum[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=cSolnl.uSum,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("U sum"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
iPlotNL=plot(dLognl.Itotal[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=cSolnl.itotal,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("Z sum"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
lamPlotNL=plot(dLognl.Lam[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=cSolnl.lamCoupl,Geom.line,Theme(default_color=colorant"black",line_width=4pt)),
			Guide.xlabel("Time"), Guide.ylabel("Lambda"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_"*updateMethod*"_LamConv.png", 36inch, 12inch), lamPlot) end

fPlotNL=plot(x=1:convIt-1,y=dCMnl.objVal[1:convIt-1,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("obj function gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
convItPlotNL=plot(x=1:convIt,y=dCMnl.lamIt[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap",orientation=:vertical),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
convPlotNL=plot(x=1:convIt,y=dCMnl.lam[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("Lambda Star Gap",orientation=:vertical),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
constPlotNL=plot(x=1:convIt,y=dCMnl.couplConst[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("Const Gap",orientation=:vertical),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_"*updateMethod*"_Conv.png", 36inch, 12inch), convPlot) end
