#Micah Botkin-Levy
#4/8/18

using Ipopt

#pull out a few key variables
N=evS.N
S=evS.S
horzLen=evS.K1

include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCnl.jl")

timeT=@elapsed cSolnl=nlEVcentral(N,S,horzLen,evS,relaxed)

relaxString= if relaxed==true "_relax"else "" end
filename = "central_NL_N$(N)"*relaxString
# save
if saveResults==1 saveRun(path,filename,timeT, evS,cSolnl) end
# load
if loadResults==1
	loadF=JLD.load(path*filename*".jld")
	evS=loadF["scenario"]
	cSolnl=loadF["solution"]
end

println("plotting....")
xPlot=zeros(horzLen+1,N)
xPlot2=zeros(horzLen+1,N)
uPlot=zeros(horzLen+1,N)
for ii= 1:N
	xPlot[:,ii]=cSolnl.sn[collect(ii:N:length(cSolnl.sn))]
	xPlot2[:,ii]=(evS.Snmin[ii,1].-xPlot[:,ii])./(evS.Kn[ii,1].-(1:1:length(xPlot[:,ii])))
    uPlot[:,ii]=cSolnl.un[collect(ii:N:length(cSolnl.un))]
end

p1nl=plot(xPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV SOC"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1,ymax=1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
p1bnl=plot(xPlot2,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV SOC"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_centralNL_SOC.png", 24inch, 12inch), p1nl) end

p2nl=plot(uPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV Current (kA)"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1,ymin=0),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_centralNL_Curr.png", 24inch, 12inch), p2nl) end

p3nl=plot(layer(x=1:horzLen+1,y=cSolnl.xt,Geom.line,Theme(default_color=colorant"blue")),
		yintercept=[evS.Tmax],Geom.hline(color=["red"],style=:dot),
		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",key_position = :top,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_centralNL_Temp.png", 24inch, 12inch), p3nl) end

# p4nlb=plot(x=1:horzLen+1,y=kappaUpperT,Geom.line,
# 		Guide.xlabel("Time"), Guide.ylabel(raw"Lambda ($/K)",orientation=:vertical),
# 		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=18pt,
# 		minor_label_font_size=16pt,key_label_font_size=16pt))
p4nl=plot(x=1:horzLen+1,y=cSolnl.lamCoupl,Geom.line,
        Guide.xlabel("Time"), Guide.ylabel(raw"Lambda ($/kA)",orientation=:vertical),
        Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=18pt,
        minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_centralNL_Lam.png", 24inch, 12inch), p4nl) end

fName="J_Central.png"

#draw(PNG(path*fName, 13inch, 14inch), vstack(p1,p2,p3,p4))



#uPlotNoLim=uPlot
#do this more elegantly
# aggU=plot(layer(x=1:horzLen+1,y=sum(uPlot[:,i] for i=1:N),Geom.line,Theme(default_color=colorant"blue")),
# 		layer(x=1:horzLen+1,y=sum(uPlotNoLim[:,i] for i=1:N),Geom.line,Theme(default_color=colorant"red")),
# 		Guide.xlabel("Time"), Guide.ylabel("PEV Current (kA)"),
# 		Coord.Cartesian(xmin=0,xmax=horzLen+1),
# 		Theme(background_color=colorant"white"),
# 		Guide.manual_color_key("", ["Aggregate Current", "Unconstrained Aggregate Current"], ["blue","red"]))
# fName="aggCentral_noLim.png"
#  draw(PNG(path*fName, 13inch, 14inch), aggU)


# p3nolimit=plot(layer(x=1:horzLen+1,y=xtRaw,Geom.line,Theme(default_color=colorant"blue")),
# 		layer(x=1:horzLen+1,y=Tactual,Geom.line,Theme(default_color=colorant"green")),
# 		layer(x=1:horzLen+1,y=noLimitt,Geom.line,Theme(default_color=colorant"red")),
# 		yintercept=[Tmax],Geom.hline(color=["red"],style=:dot),
# 		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)",orientation=:vertical),
# 		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",key_position = :top),
# 		Guide.manual_color_key("", ["PWL Temp", "Actual Temp","Unconstrained Temp"], ["blue", "green","red"]))
# fName="Temp_w_nolimit.png"
# draw(PNG(path*fName, 13inch, 14inch), p3nolimit)
