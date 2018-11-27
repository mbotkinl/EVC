#Micah Botkin-Levy
#4/8/18
using Gurobi
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCpwl.jl")

#pull out a few key variables
N=evS.N
S=evS.S
horzLen=evS.K1
errorString= if forecastError "_Error" else "" end
fname = "central_N$(N)"*errorString

if loadResults
	println("Reading in Central Sim")
	loadF=load(path*fname*".jld2")
	evS=loadF["scenario"]
	cSol=loadF["solution"]
else
	println("Running Central Sim")
	timeT=@elapsed cSol=pwlEVcentral(N,S,horzLen,evS,forecastError,slack)
	if saveResults saveRun(path,fname,timeT, evS,cSol) end
end


println("plotting....")
snPlot=zeros(horzLen+1,N)
#xPlot2=zeros(horzLen+1,N)
uPlot=zeros(horzLen+1,N)
for ii= 1:N
	snPlot[:,ii]=cSol.Sn[collect(ii:N:length(cSol.Sn))]
	#xPlot2[:,ii]=(evS.Snmin[ii,1]-xPlot[:,ii])./(evS.Kn[ii,1]-(1:1:length(xPlot[:,ii])))
    uPlot[:,ii]=cSol.Un[collect(ii:N:length(cSol.Un))]
end


p1=plot(snPlot,xlabel="Time",ylabel="PEV SOC",legend=false,xlims=(1,horzLen+1))
if drawFig savefig(p1,path*"J_central_SOC.png") end

p2=plot(uPlot,xlabel="Time",ylabel="PEV Current (kA)",legend=false,xlims=(1,horzLen+1))
if drawFig savefig(p2,path*"J_central_Curr.png") end

p3=plot(1:horzLen+1,hcat(cSol.Xt*1000,cSol.Tactual*1000),label=["PWL Temp" "Actual Temp"],xlims=(1,horzLen+1),xlabel="Time",ylabel="Temp (K)")
plot!(p3,1:horzLen+1,evS.Tmax*ones(horzLen+1)*1000,label="XFRM Limit",line=(:dash,:red))
if drawFig savefig(p3,path*"J_central_Temp.png") end

p4b=plot(1:horzLen+1,cSol.lamTemp,xlabel="Time",ylabel=raw"Lambda ($/K)",xlims=(1,horzLen+1),legend=false)
p4=plot(1:horzLen+1,cSol.lamCoupl,xlabel="Time",ylabel=raw"Lambda ($/kA)",xlims=(1,horzLen+1),legend=false)
if drawFig savefig(p4,path*"J_central_Lam.png") end

fName="J_Central.png"

#draw(PNG(path*fName, 13inch, 14inch), vstack(p1,p2,p3,p4))

checkDesiredStates(cSol.Sn,evS.Kn,evS.Snmin)


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
