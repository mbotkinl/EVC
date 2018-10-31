#EVC with ALADIN for PWL convex relaxation
#Micah Botkin-Levy
#5/22/18

#formulation try 2
#u, sn, xt, and z are all in "y", v is the auxilliary variable corresponding to x in literature
#current constraint is coupling


filename = "dALADIN_N$(N)"

if loadResults
	loadF=load(path*filename*".jld2")
	evS=loadF["scenario"]
	dLogalad=loadF["solution"]
	dCMalad=loadF["convMetrics"]
	convIt=loadF["convIt"]
else
	timeT=@elapsed dLogalad,dCMalad,convIt,ΔY,convCheck=pwlEValad(N,S,horzLen,maxIt,evS,cSol,slack)
	if saveResults saveRun(path,filename,timeT, evS,dLogalad, dCMalad, convIt) end
end


println("plotting....")
xPlot=zeros(horzLen+1,N)
uPlot=zeros(horzLen+1,N)
for ii= 1:N
	xPlot[:,ii]=dLogalad.Sn[collect(ii:N:length(dLogalad.Sn[:,convIt])),convIt]
    uPlot[:,ii]=dLogalad.Un[collect(ii:N:length(dLogalad.Un[:,convIt])),convIt]
end

pd1alad=plot(xPlot,xlabel="Time",ylabel="PEV SOC",legend=false,xlims=(0,horzLen+1),ylims=(0,1))
if drawFig==1 savefig(pd1alad,path*"J_decentral_ALADIN_SOC.png") end

pd2alad=plot(uPlot,xlabel="Time",ylabel="PEV Current (kA)",legend=false,xlims=(0,horzLen+1))
if drawFig==1 savefig(pd2alad,path*"J_decentral_ALADIN_Curr.png") end

pd3alad=plot(1:horzLen+1,hcat(dLogalad.Tactual[:,convIt],dLogalad.Xt[:,convIt]),label=["Actual Temp" "PWL Temp"],xlims=(0,horzLen+1),xlabel="Time",ylabel="Temp (K)")
plot!(pd3alad,1:horzLen+1,evS.Tmax*ones(horzLen+1),label="XFRM Limit",line=(:dash,:red))
if drawFig==1 savefig(pd3alad,path*"J_decentral_ALADIN_Temp.png") end

pd4alad=plot(1:horzLen+1,hcat(cSol.lamCoupl,dLogalad.Lam[:,convIt]),xlabel="Time",ylabel=raw"Lambda ($/kA)",
             xlims=(0,horzLen+1),labels=["Central" "ADMM"])
if drawFig==1 savefig(pd4alad,path*"J_decentral_ALADIN_Lam.png") end

#
# lamPlotalad=plot(dLogalad.Lam[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
# 			layer(x=1:horzLen+1,y=cSol.lamCoupl,Geom.line,Theme(default_color=colorant"black",line_width=4pt)),
# 			Guide.xlabel("Time"), Guide.ylabel("Lambda"),Guide.ColorKey(title="Iteration"),
# 			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
# 			minor_label_font_size=26pt,key_label_font_size=26pt))
# uSumPlotalad=plot(dLogalad.uSum[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
# 			layer(x=1:horzLen+1,y=cSol.uSum,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
# 			Guide.xlabel("Time"), Guide.ylabel("U sum"),Guide.ColorKey(title="Iteration"),
# 			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
# 			minor_label_font_size=26pt,key_label_font_size=26pt))
# zSumPlotalad=plot(dLogalad.zSum[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
# 			layer(x=1:horzLen+1,y=cSol.zSum,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
# 			Guide.xlabel("Time"), Guide.ylabel("Z sum"),Guide.ColorKey(title="Iteration"),
# 			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
# 			minor_label_font_size=26pt,key_label_font_size=26pt))
# constPlotalad2=plot(dLogalad.couplConst[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
# 			Guide.xlabel("Time"), Guide.ylabel("curr constraint diff",orientation=:vertical),Guide.ColorKey(title="Iteration"),
# 			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
# 			minor_label_font_size=26pt,key_label_font_size=26pt))
# if drawFig==1 draw(PNG(path*"J_ALADIN_LamConv.png", 36inch, 12inch), lamPlotalad) end

activeSet=zeros(convIt,1)
setChanges=zeros(convIt,1)
for ii=2:convIt
    activeSet[ii,1]=sum(abs.(dLogalad.Cs[:,ii]))+sum(abs.(dLogalad.Ct[:,ii]))+
              sum(abs.(dLogalad.Cu[:,ii]))+sum(abs.(dLogalad.Cz[:,ii]))
    setChanges[ii,1]=sum(abs.(dLogalad.Cs[:,ii]-dLogalad.Cs[:,ii-1]))+sum(abs.(dLogalad.Ct[:,ii]-dLogalad.Ct[:,ii-1]))+
                     sum(abs.(dLogalad.Cu[:,ii]-dLogalad.Cu[:,ii-1]))+sum(abs.(dLogalad.Cz[:,ii]-dLogalad.Cz[:,ii-1]))
end
activeSetPlot=plot(2:convIt,activeSet[2:convIt],xlabel="Iteration",ylabel="Total Active inequality constraints",
                   legend=false,xlims=(2,convIt))
setChangesPlot=plot(2:convIt,setChanges[2:convIt],xlabel="Iteration",ylabel="Changes in Active inequality constraints",
                  legend=false,xlims=(2,convIt))
solChangesplot=plot(2:convIt,hcat(ΔY[2:convIt],convCheck[2:convIt]),xlabel="Iteration",labels=["ΔY" "y-x"],xlims=(1,convIt))

fPlotalad=plot(1:horzLen+1,dCMalad.obj[1:convIt-1,1],xlabel="Iteration",ylabel="obj function gap",xlims=(0,convIt),legend=false,yscale=:log10)
convItPlotalad=plot(1:horzLen+1,dCMalad.lamIt[1:convIt-1,1],xlabel="Iteration",ylabel="2-Norm Lambda Gap",xlims=(0,convIt),legend=false,yscale=:log10)
convPlotalad=plot(1:horzLen+1,dCMalad.lam[1:convIt-1,1],xlabel="Iteration",ylabel="central lambda gap",xlims=(0,convIt),legend=false,yscale=:log10)
constPlotalad=plot(1:horzLen+1,dCMalad.couplConst[1:convIt-1,1],xlabel="Iteration",ylabel="curr constraint Gap",xlims=(0,convIt),legend=false,yscale=:log10)

#checkDesiredStates(dLogalad.Sn,evS.Kn,evS.Snmin)
