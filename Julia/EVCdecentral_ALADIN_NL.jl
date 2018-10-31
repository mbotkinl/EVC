
#EVC with ALADIN for PWL convex relaxation
#Micah Botkin-Levy
#5/22/18

#formulation try 2
#u, sn, xt, and z are all in "y", v is the auxilliary variable corresponding to x in literature
#current constraint is coupling


timeT=@elapsed dLognlalad,dCMnlalad,convIt,ΔY,convCheck=nlEValad(N,S,horzLen,maxIt,evS,cSolnl,relaxed,slack)

relaxString= if relaxed==true "_relax"else "" end
filename = "dALADIN_NL_N$(N)"*relaxString
# save
if saveResults==1 saveRun(path,filename,timeT, evS,dLognlalad, dCMnlalad, convIt) end
# load
if loadResults==1
	loadF=load(path*filename*".jld2")
	evS=loadF["scenario"]
	dLognlalad=loadF["solution"]
	dCMnlalad=loadF["convMetrics"]
	convIt=loadF["convIt"]
end

println("plotting....")
xPlot=zeros(horzLen+1,N)
uPlot=zeros(horzLen+1,N)
for ii= 1:N
	xPlot[:,ii]=dLognlalad.Sn[collect(ii:N:length(dLognlalad.Sn[:,convIt])),convIt]
    uPlot[:,ii]=dLognlalad.Un[collect(ii:N:length(dLognlalad.Un[:,convIt])),convIt]
end

pd1aladNL=plot(xPlot,xlabel="Time",ylabel="PEV SOC",legend=false,xlims=(0,horzLen+1),ylims=(0,1))
if drawFig==1 savefig(pd1aladNL,path*"J_decentralNL_ALADIN_SOC.png") end

pd2aladNL=plot(uPlot,xlabel="Time",ylabel="PEV Current (kA)",legend=false,xlims=(0,horzLen+1))
if drawFig==1 savefig(pd2aladNL,path*"J_decentralNL_ALADIN_Curr.png") end

pd3aladNL=plot(1:horzLen+1,dLognlalad.Xt[:,convIt),label="XFRM Temp",xlims=(0,horzLen+1),xlabel="Time",ylabel="Temp (K)")
plot!(pd3aladNL,1:horzLen+1,evS.Tmax*ones(horzLen+1),label="XFRM Limit",line=(:dash,:red))
if drawFig==1 savefig(pd3aladNL,path*"J_decentralNL_ALADIN_Temp.png") end

pd4aladNL=plot(1:horzLen+1,hcat(cSolnl.lamCoupl,dLognlalad.Lam[:,convIt]),xlabel="Time",ylabel=raw"Lambda ($/kA)",
             xlims=(0,horzLen+1),labels=["Central" "ADMM"])
if drawFig==1 savefig(pd4aladNL,path*"J_decentralNL_ALADIN_Lam.png") end

#
#
# lamPlotaladNL=plot(dLognlalad.Lam[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
# 			layer(x=1:horzLen+1,y=cSolnl.lamCoupl,Geom.line,Theme(default_color=colorant"black",line_width=4pt)),
# 			Guide.xlabel("Time"), Guide.ylabel("Lambda"),Guide.ColorKey(title="Iteration"),
# 			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
# 			minor_label_font_size=26pt,key_label_font_size=26pt))
# uSumPlotaladNL=plot(dLognlalad.uSum[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
# 			layer(x=1:horzLen+1,y=cSolnl.uSum,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
# 			Guide.xlabel("Time"), Guide.ylabel("U sum"),Guide.ColorKey(title="Iteration"),
# 			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
# 			minor_label_font_size=26pt,key_label_font_size=26pt))
# iPlotaladNL=plot(dLognlalad.Itotal[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
# 			layer(x=1:horzLen+1,y=cSolnl.Itotal,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
# 			Guide.xlabel("Time"), Guide.ylabel("Z sum"),Guide.ColorKey(title="Iteration"),
# 			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
# 			minor_label_font_size=26pt,key_label_font_size=26pt))
# xtPlotaladNL=plot(dLognlalad.Xt[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
# 			layer(x=1:horzLen+1,y=cSolnl.Xt,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
# 			Guide.xlabel("Time"), Guide.ylabel("Z sum"),Guide.ColorKey(title="Iteration"),
# 			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
# 			minor_label_font_size=26pt,key_label_font_size=26pt))
# constPlotaladNL2=plot(dLognlalad.couplConst,x=Row.index,y=Col.value,color=Col.index,Geom.line,
# 			Guide.xlabel("Time"), Guide.ylabel("curr constraint diff"),Guide.ColorKey(title="Iteration"),
# 			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
# 			minor_label_font_size=26pt,key_label_font_size=26pt))
# if drawFig==1 draw(PNG(path*"J_ALADIN_LamConv.png", 36inch, 12inch), lamPlotalad) end

activeSet=zeros(convIt,1)
setChanges=zeros(convIt,1)
for ii=1:convIt
    activeSet[ii,1]=sum(abs.(dLognlalad.Cs[:,ii]))+sum(abs.(dLognlalad.Ct[:,ii]))+
              sum(abs.(dLognlalad.Cu[:,ii]))+sum(abs.(dLognlalad.Ci[:,ii]))
    setChanges[ii,1]=sum(abs.(dLognlalad.Cs[:,ii]-dLognlalad.Cs[:,ii-1]))+sum(abs.(dLognlalad.Ct[:,ii]-dLognlalad.Ct[:,ii-1]))+
                     sum(abs.(dLognlalad.Cu[:,ii]-dLognlalad.Cu[:,ii-1]))+sum(abs.(dLognlalad.Ci[:,ii]-dLognlalad.Ci[:,ii-1]))
end

activeSetPlot=plot(1:convIt,activeSet[1:convIt],xlabel="Iteration",ylabel="Total Active inequality constraints",
                   legend=false,xlims=(1,convIt))
setChangesPlot=plot(1:convIt,setChanges[1:convIt],xlabel="Iteration",ylabel="Changes in Active inequality constraints",
                  legend=false,xlims=(1,convIt))
solChangesplot=plot(1:convIt,hcat(ΔY[1:convIt],convCheck[1:convIt]),xlabel="Iteration",labels=["ΔY" "y-x"],xlims=(1,convIt))

fPlotaladNL=plot(1:horzLen+1,dCMnlalad.obj[1:convIt-1,1],xlabel="Iteration",xlims=(0,convIt),legend=false)
yaxis!(fPlotaladNL,"obj function gap",:log10)
convItPlotaladNL=plot(1:horzLen+1,dCMnlalad.lamIt[1:convIt-1,1],xlabel="Iteration",xlims=(0,convIt),legend=false)
yaxis!(convItPlotaladNL,"2-Norm Lambda Gap",:log10)
convPlotaladNL=plot(1:horzLen+1,dCMnlalad.lam[1:convIt-1,1],xlabel="Iteration",xlims=(0,convIt),legend=false)
yaxis!(convPlotaladNL,"central lambda gap",:log10)
constPlotaladNL=plot(1:horzLen+1,dCMnlalad.couplConst[1:convIt-1,1],xlabel="Iteration",xlims=(0,convIt),legend=false)
yaxis!(constPlotaladNL,"curr constraint Gap",:log10)


convItPlotaladNL=plot(x=1:convIt,y=dCMnlalad.lamIt[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
convPlotaladNL=plot(x=1:convIt,y=dCMnlalad.lam[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("central lambda gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
constPlotaladNL=plot(x=1:convIt,y=dCMnlalad.couplConst[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("consensus gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
fPlotaladNL=plot(x=1:convIt-1,y=dCMnlalad.obj[1:convIt-1,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("obj function gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_ALADIN_Conv.png", 36inch, 12inch), vstack(convItPlotalad,convPlotalad,fPlotalad)) end
