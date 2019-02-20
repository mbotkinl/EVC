#Micah Botkin-Levy
#4/8/18
if relaxedMode==2 #SOCP
	using Mosek
elseif relaxedMode==1 #QCQP
	using Gurobi
else
	@everywhere using Gurobi
	using Ipopt
end
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCnl.jl")

relaxString= "_R$(relaxedMode)"
errorString= if forecastError "_Error" else "" end
fname = "central_NL_N$(N)"*relaxString*errorString

if loadResults
	println("Reading in NL Central Sim")
	loadF=load(path*fname*".jld2")
	cSavenl=loadF["centralLog"]
	cSolnl=loadF["solution"]
else
	#initialize
	t0=evS.t0
	s0=evS.s0

	println("Running NL Central Sim")
	timeT=@elapsed cSolnl,cSave=nlEVcentral(evS,slack,relaxedMode,silent)
	if saveResults saveRun(path,fname,timeT, cSolnl,cSave=cSavenl) end
end

println("plotting....")

p1nl=plot(cSolnl.Sn,xlabel="Time",ylabel="PEV SOC",legend=false,xlims=(1,evS.K),ylims=(0,1))
if drawFig==1 savefig(p1nl,path*"J_centralNL_SOC.png") end
p2nl=plot(cSolnl.Un,xlabel="Time",ylabel="PEV Current (kA)",legend=false,xlims=(1,evS.K))
if drawFig==1 savefig(p2nl,path*"J_centralNL_Curr.png") end

p3nl=plot(cSolnl.Tactual*1000,label="XFRM Temp",xlims=(1,evS.K),xlabel="Time",ylabel="Temp (K)")
plot!(p3nl,evS.Tmax*ones(evS.K)*1000,label="XFRM Limit",line=(:dash,:red))
plot!(p3nl,cSolnl.Tpred*1000,label="Predict Temp")
if drawFig==1 savefig(p3nl,path*"J_centralNL_Temp.png") end

p4nl=plot(cSolnl.lamCoupl,xlabel="Time",ylabel=raw"Lambda ($/kA)",xlims=(1,evS.K),legend=false)
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
