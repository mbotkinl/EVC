
#EVC with ALADIN for PWL convex relaxation
#Micah Botkin-Levy
#5/22/18

#formulation try 2
#u, sn, xt, and z are all in "y", v is the auxilliary variable corresponding to x in literature
#current constraint is coupling
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCnl.jl")

eqString=if eqForm "_eq" else "_ineq" end
relaxString= "_R$(relaxedMode)"
fname = "dALADIN_NL_N$(N)"*eqString*relaxString

if loadResults
	println("Reading in NL AlAD Sim")
	loadF=load(path*fname*".jld2")
	evS=loadF["scenario"]
	dLognlalad=loadF["solution"]
	dCMnlalad=loadF["convMetrics"]
	convIt=loadF["convIt"]
else
	t0=evS.t0
	s0=evS.s0

	#prevLam=cSolnl.lamCoupl[1:evS.K1+1]
	prevLam=2e3*ones(evS.K1+1,1)
    prevVt=evS.Tmax*ones(evS.K1+1,1)
    prevVi=evS.ItotalMax*ones(evS.K1+1,1)
    prevVu=.01*ones(evS.N*(evS.K1+1),1)
    prevVs=.5*ones(evS.N*(evS.K1+1),1)
	ρALADp = 1e3


	println("Running NL AlAD Sim")
	timeT=@elapsed dSolaladnl,dCMaladnl=nlEValad(maxIt,evS,cSavenl,slack,eqForm,silent)
	if saveResults saveRun(path,fname,timeT,dSolaladnl,convMetrics=dCMaladnl) end
end

println("plotting....")

pd1aladNL=plot(xPlot,xlabel="Time",ylabel="PEV SOC",legend=false,xlims=(0,horzLen+1),ylims=(0,1))
if drawFig==1 savefig(pd1aladNL,path*"J_decentralNL_ALADIN_SOC.png") end

pd2aladNL=plot(uPlot,xlabel="Time",ylabel="PEV Current (kA)",legend=false,xlims=(0,horzLen+1))
if drawFig==1 savefig(pd2aladNL,path*"J_decentralNL_ALADIN_Curr.png") end

pd3aladNL=plot(1:horzLen+1,dLognlalad.Xt[:,convIt]*1000,label="XFRM Temp",xlims=(0,horzLen+1),xlabel="Time",ylabel="Temp (K)")
plot!(pd3aladNL,1:horzLen+1,evS.Tmax*ones(horzLen+1)*1000,label="XFRM Limit",line=(:dash,:red))
if drawFig==1 savefig(pd3aladNL,path*"J_decentralNL_ALADIN_Temp.png") end

pd4aladNL=plot(1:horzLen+1,hcat(cSolnl.lamCoupl,dLognlalad.Lam[:,convIt]),xlabel="Time",ylabel=raw"Lambda ($/kA)",
             xlims=(0,horzLen+1),labels=["Central" "ALADIN"])
if drawFig==1 savefig(pd4aladNL,path*"J_decentralNL_ALADIN_Lam.png") end


#convergence plots
halfCI=Int(floor(convIt/2))
CList=reshape([range(colorant"blue", stop=colorant"yellow",length=halfCI);
               range(colorant"yellow", stop=colorant"red",length=convIt-halfCI)], 1, convIt);

uSumPlotnlalad=plot(dLognlalad.uSum[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Current Sum (kA)",xlims=(0,horzLen+1),legend=false)
plot!(uSumPlotnlalad,1:horzLen+1,cSolnl.uSum,seriescolor=:black,linewidth=2,linealpha=0.8)

constPlotnlalad=plot(dLognlalad.couplConst[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="curr constraint diff",xlims=(0,horzLen+1),legend=false)

lamPlotnlalad=plot(dLognlalad.Lam[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Lambda",xlims=(0,horzLen+1),legend=false)
plot!(lamPlotnlalad,1:horzLen+1,cSolnl.lamCoupl,seriescolor=:black,linewidth=2,linealpha=0.8)

gradPlots=plot(uSumPlotnlalad, constPlotnlalad, lamPlotnlalad,layout=(3,1))
pubPlot(gradPlots,thickscale=0.8,sizeWH=(1000,600),dpi=300)
if drawFig savefig(gradPlots,path*"J_decentral_ALADIN_gradPlots.png") end


activeSet=zeros(convIt,1)
setChanges=zeros(convIt,1)
for ii=2:convIt
    activeSet[ii,1]=sum(abs.(dLognlalad.Csu[:,ii]))+sum(abs.(dLognlalad.Ctu[:,ii]))+
              		sum(abs.(dLognlalad.Cuu[:,ii]))+sum(abs.(dLognlalad.Ciu[:,ii]))+
					sum(abs.(dLognlalad.Csl[:,ii]))+sum(abs.(dLognlalad.Ctl[:,ii]))+
				    sum(abs.(dLognlalad.Cul[:,ii]))+sum(abs.(dLognlalad.Cil[:,ii]))
    setChanges[ii,1]=sum(abs.(dLognlalad.Csu[:,ii]-dLognlalad.Csu[:,ii-1]))+sum(abs.(dLognlalad.Ctu[:,ii]-dLognlalad.Ctu[:,ii-1]))+
                     sum(abs.(dLognlalad.Cuu[:,ii]-dLognlalad.Cuu[:,ii-1]))+sum(abs.(dLognlalad.Ciu[:,ii]-dLognlalad.Ciu[:,ii-1]))+
					 sum(abs.(dLognlalad.Csl[:,ii]-dLognlalad.Csl[:,ii-1]))+sum(abs.(dLognlalad.Ctl[:,ii]-dLognlalad.Ctl[:,ii-1]))+
				     sum(abs.(dLognlalad.Cul[:,ii]-dLognlalad.Cul[:,ii-1]))+sum(abs.(dLognlalad.Cil[:,ii]-dLognlalad.Cil[:,ii-1]))
end

activeSetPlot=plot(2:convIt,activeSet[2:convIt],xlabel="Iteration",ylabel="Total Active inequality constraints",
                   legend=false,xlims=(2,convIt))
setChangesPlot=plot(2:convIt,setChanges[2:convIt],xlabel="Iteration",ylabel="Changes in Active inequality constraints",
                  legend=false,xlims=(2,convIt))
solChangesplot=plot(2:convIt,hcat(ΔY[1:convIt],convCheck[1:convIt]),xlabel="Iteration",labels=["ΔY" "y-x"],xlims=(1,convIt))

setChangePlots=plot(activeSetPlot, setChangesPlot,layout=(2,1))
pubPlot(setChangePlots,thickscale=0.5,sizeWH=(1000,600),dpi=300)
if drawFig savefig(setChangePlots,path*"J_decentralNL_ALADIN_setPlots.png") end


fPlotaladNL=plot(1:convIt,dCMnlalad.obj[1:convIt,1],xlabel="Iteration",ylabel="obj function gap",xlims=(1,convIt),legend=false,yscale=:log10)
convItPlotaladNL=plot(1:convIt,dCMnlalad.lamIt[1:convIt,1],xlabel="Iteration",ylabel="2-Norm Lambda Gap",xlims=(1,convIt),legend=false,yscale=:log10)
convPlotaladNL=plot(1:convIt,dCMnlalad.lam[1:convIt,1],xlabel="Iteration",ylabel="central lambda gap",xlims=(1,convIt),legend=false,yscale=:log10)
constPlotaladNL=plot(1:convIt,dCMnlalad.couplConst[1:convIt,1],xlabel="Iteration",ylabel="curr constraint Gap",xlims=(1,convIt),legend=false,yscale=:log10)

convPlotsNL=plot(fPlotaladNL, convItPlotaladNL,convPlotaladNL, constPlotaladNL,layout=(2,2))
pubPlot(convPlotsNL,thickscale=0.8,sizeWH=(1000,600),dpi=300)
if drawFig savefig(convPlotsNL,path*"J_decentralNL_ALADIN_convPlots.png") end

checkDesiredStates(dLognlalad.Sn[:,convIt],evS.Kn,evS.Snmin)
