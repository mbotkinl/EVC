#Micah Botkin-Levy
#4/8/18
if relaxedMode==2 #SOCP
	using Mosek
elseif relaxedMode==1 #QCQP
	using Gurobi
else
	using Ipopt
end
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCnl.jl")

#pull out a few key variables
N=evS.N
S=evS.S
horzLen=evS.K1

relaxString= "_R$(relaxedMode)"
fname = "central_NL_N$(N)"*relaxString

if loadResults
	println("Reading in NL Central Sim")
	loadF=load(path*fname*".jld2")
	evS=loadF["scenario"]
	cSolnl=loadF["solution"]
else
	println("Running NL Central Sim")
	timeT=@elapsed cSolnl=nlEVcentral(N,S,horzLen,evS,relaxedMode,slack)
	if saveResults saveRun(path,fname,timeT, evS,cSolnl) end
end

println("plotting....")
snPlot=zeros(horzLen+1,N)
#xPlot2=zeros(horzLen+1,N)
uPlot=zeros(horzLen+1,N)
for ii= 1:N
	snPlot[:,ii]=cSolnl.Sn[collect(ii:N:length(cSolnl.Sn))]
	#xPlot2[:,ii]=(evS.Snmin[ii,1]-xPlot[:,ii])./(evS.Kn[ii,1]-(1:1:length(xPlot[:,ii])))
    uPlot[:,ii]=cSolnl.Un[collect(ii:N:length(cSolnl.Un))]
end


p1nl=plot(snPlot,xlabel="Time",ylabel="PEV SOC",legend=false,xlims=(1,horzLen+1),ylims=(0,1))
if drawFig==1 savefig(p1nl,path*"J_centralNL_SOC.png") end
p2nl=plot(uPlot,xlabel="Time",ylabel="PEV Current (kA)",legend=false,xlims=(1,horzLen+1))
if drawFig==1 savefig(p2nl,path*"J_centralNL_Curr.png") end

p3nl=plot(1:horzLen+1,cSolnl.Xt,label="XFRM Temp",xlims=(1,horzLen+1),xlabel="Time",ylabel="Temp (K)")
plot!(p3nl,1:horzLen+1,evS.Tmax*ones(horzLen+1),label="XFRM Limit",line=(:dash,:red))
if drawFig==1 savefig(p3nl,path*"J_centralNL_Temp.png") end

p4nl=plot(1:horzLen+1,cSolnl.lamCoupl,xlabel="Time",ylabel=raw"Lambda ($/kA)",xlims=(1,horzLen+1),legend=false)
if drawFig==1 savefig(p4nl,path*"J_centralNL_Lam.png") end

fName="J_Central.png"
#draw(PNG(path*fName, 13inch, 14inch), vstack(p1,p2,p3,p4))

checkDesiredStates(cSolnl.Sn,evS.Kn,evS.Snmin)


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
