#EVC with ADMM for PWL convex relaxation
#Micah Botkin-Levy
#5/22/18

#formulation try 2
#u, sn, xt, and z are all in "x" v is the auxilliary variable corresponding to z in literature
#current constraint is coupling


tic()
dLognladmm,dCMnladmm,convIt=nlEVadmm(N,S,horzLen,maxIt,evS,cSolnl,relaxed)
timeT=toc()

relaxString= if relaxed==true "_relax"else "" end
filename = "dADMM_NL_N$(N)"*relaxString
# save
if saveResults==1 saveRun(path,filename,timeT, evS,dLognladmm, dCMnladmm, convIt) end
# load
if loadResults==1
	loadF=JLD.load(path*filename*".jld")
	evS=loadF["scenario"]
	dLognladmm=loadF["solution"]
	dCMnladmm=loadF["convMetrics"]
	convIt=loadF["convIt"]
end

println("plotting....")
xPlot=zeros(horzLen+1,N)
uPlot=zeros(horzLen+1,N)
for ii= 1:N
	xPlot[:,ii]=dLognladmm.Sn[collect(ii:N:length(dLognladmm.Sn[:,convIt])),convIt]
    uPlot[:,ii]=dLognladmm.Un[collect(ii:N:length(dLognladmm.Un[:,convIt])),convIt]
end

pd1NLadmm=plot(xPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV SOC"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_centralNL_ADMM_SOC.png", 24inch, 12inch), pd1NLadmm) end

pd2NLadmm=plot(uPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV Current (A)"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_centralNL_ADMM_Curr.png", 24inch, 12inch), pd2NLadmm) end

pd3NLadmm=plot(layer(x=1:horzLen+1,y=dLognladmm.Xt[:,convIt],Geom.line,Theme(default_color=colorant"blue")),
		#layer(x=1:horzLen+1,y=Tactual,Geom.line,Theme(default_color=colorant"green")),
		yintercept=[evS.Tmax],Geom.hline(color=["red"],style=:dot),
		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",key_position = :top,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
		#Guide.manual_color_key("", ["PWL Temp", "Actual Temp"], ["blue", "green"]))
if drawFig==1 draw(PNG(path*"J_centralNL_ADMM_Temp.png", 24inch, 12inch), pd3NLadmm) end

pd4NLadmm=plot(x=1:horzLen+1,y=dLognladmm.Lam[:,convIt],Geom.line,
		Guide.xlabel("Time"), Guide.ylabel(raw"Lambda ($/A)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_centralNL_ADMM_Lam.png", 24inch, 12inch), pd4NLadmm) end

fName="J_Central.png"


lamPlotNLadmm=plot(dLognladmm.Lam[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
            layer(x=1:horzLen+1,y=cSolnl.lamCoupl,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("Lambda",orientation=:vertical),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
constPlotNLadmm2=plot(dLognladmm.couplConst[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel("curr constraint diff",orientation=:vertical),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_ADMM_LamConv.png", 36inch, 12inch), lamPlotadmm) end

convItPlotNLadmm=plot(x=1:convIt,y=dCMnladmm.lamIt[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap",orientation=:vertical),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
convPlotNLadmm=plot(x=1:convIt,y=dCMnladmm.lam[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap",orientation=:vertical),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
fPlotadmm=plot(x=1:convIt-1,y=dCMnladmm.objVal[1:convIt-1,1],Geom.line,#Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap",orientation=:vertical),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
constPlotadmm=plot(x=1:convIt,y=dCMnladmm.couplConst[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap",orientation=:vertical),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_ADMM_Conv.png", 36inch, 12inch), convPlotadmm) end
