#EVC with ALADIN for PWL convex relaxation
#Micah Botkin-Levy
#5/22/18

#formulation try 2
#u, sn, xt, and z are all in "y", v is the auxilliary variable corresponding to x in literature
#current constraint is coupling

include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCpwl.jl")
fname = "dALADIN_N$(N)"

if loadResults
	println("Reading in ALAD Sim")
	loadF=load(path*fname*".jld2")
	evS=loadF["scenario"]
	dLogalad=loadF["solution"]
	dCMalad=loadF["convMetrics"]
	convIt=loadF["convIt"]
else
	println("Running ALAD Sim")
	timeT=@elapsed dLogalad,dCMalad,convIt,ΔY,convCheck=pwlEValad(N,S,horzLen,maxIt,evS,cSol,slack)
	if saveResults saveRun(path,fname,timeT, evS,dLogalad, dCMalad, convIt) end
end


println("plotting....")
xPlot=zeros(horzLen+1,N)
uPlot=zeros(horzLen+1,N)
for ii= 1:N
	xPlot[:,ii]=dLogalad.Sn[collect(ii:N:length(dLogalad.Sn[:,convIt])),convIt]
    uPlot[:,ii]=dLogalad.Un[collect(ii:N:length(dLogalad.Un[:,convIt])),convIt]
end

pd1alad=plot(xPlot,xlabel="Time",ylabel="PEV SOC",legend=false,xlims=(0,horzLen+1),ylims=(0,1))
if drawFig savefig(pd1alad,path*"J_decentral_ALADIN_SOC.png") end

pd2alad=plot(uPlot,xlabel="Time",ylabel="PEV Current (kA)",legend=false,xlims=(0,horzLen+1))
if drawFig savefig(pd2alad,path*"J_decentral_ALADIN_Curr.png") end

pd3alad=plot(1:horzLen+1,hcat(dLogalad.Tactual[:,convIt],dLogalad.Xt[:,convIt]),label=["Actual Temp" "PWL Temp"],xlims=(0,horzLen+1),xlabel="Time",ylabel="Temp (K)")
plot!(pd3alad,1:horzLen+1,evS.Tmax*ones(horzLen+1),label="XFRM Limit",line=(:dash,:red))
if drawFig savefig(pd3alad,path*"J_decentral_ALADIN_Temp.png") end

pd4alad=plot(1:horzLen+1,hcat(cSol.lamCoupl,dLogalad.Lam[:,convIt]),xlabel="Time",ylabel=raw"Lambda ($/kA)",
             xlims=(0,horzLen+1),labels=["Central" "ALADIN"])
if drawFig savefig(pd4alad,path*"J_decentral_ALADIN_Lam.png") end




#convergence plots
halfCI=Int(floor(convIt/2))
CList=reshape([range(colorant"blue", stop=colorant"yellow",length=halfCI);
               range(colorant"yellow", stop=colorant"red",length=convIt-halfCI)], 1, convIt);

uSumPlotalad=plot(dLogalad.uSum[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Current Sum (kA)",xlims=(0,horzLen+1),legend=false)
plot!(uSumPlotalad,1:horzLen+1,cSol.uSum,seriescolor=:black,linewidth=2,linealpha=0.8)

zSumPlotalad=plot(dLogalad.zSum[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Z sum",xlims=(0,horzLen+1),legend=false)
plot!(zSumPlotalad,1:horzLen+1,cSol.zSum,seriescolor=:black,linewidth=2,linealpha=0.8)

constPlotalad2=plot(dLogalad.couplConst[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="curr constraint diff",xlims=(0,horzLen+1),legend=false)

lamPlotalad=plot(dLogalad.Lam[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Lambda",xlims=(0,horzLen+1),legend=false)
plot!(lamPlotalad,1:horzLen+1,cSol.lamCoupl,seriescolor=:black,linewidth=2,linealpha=0.8)

gradPlots=plot(uSumPlotalad, zSumPlotalad,constPlotalad2, lamPlotalad,layout=(2,2))
pubPlot(gradPlots,thickscale=0.8,sizeWH=(1000,600),dpi=300)
if drawFig savefig(gradPlots,path*"J_decentral_ALADIN_gradPlots.png") end

activeSet=zeros(convIt,1)
setChanges=zeros(convIt,1)
for ii=2:convIt
    activeSet[ii,1]=sum(abs.(dLogalad.Csu[:,ii]))+sum(abs.(dLogalad.Ctu[:,ii]))+
              		sum(abs.(dLogalad.Cuu[:,ii]))+sum(abs.(dLogalad.Czu[:,ii]))+
					sum(abs.(dLogalad.Csl[:,ii]))+sum(abs.(dLogalad.Ctl[:,ii]))+
				    sum(abs.(dLogalad.Cul[:,ii]))+sum(abs.(dLogalad.Czl[:,ii]))
    setChanges[ii,1]=sum(abs.(dLogalad.Csu[:,ii]-dLogalad.Csu[:,ii-1]))+sum(abs.(dLogalad.Ctu[:,ii]-dLogalad.Ctu[:,ii-1]))+
                     sum(abs.(dLogalad.Cuu[:,ii]-dLogalad.Cuu[:,ii-1]))+sum(abs.(dLogalad.Czu[:,ii]-dLogalad.Czu[:,ii-1]))+
					 sum(abs.(dLogalad.Csl[:,ii]-dLogalad.Csl[:,ii-1]))+sum(abs.(dLogalad.Ctl[:,ii]-dLogalad.Ctl[:,ii-1]))+
				     sum(abs.(dLogalad.Cul[:,ii]-dLogalad.Cul[:,ii-1]))+sum(abs.(dLogalad.Czl[:,ii]-dLogalad.Czl[:,ii-1]))
end

activeSetPlot=plot(2:convIt,activeSet[2:convIt],xlabel="Iteration",ylabel="Total Active inequality constraints",
                   legend=false,xlims=(2,convIt))
setChangesPlot=plot(2:convIt,setChanges[2:convIt],xlabel="Iteration",ylabel="Changes in Active inequality constraints",
                  legend=false,xlims=(2,convIt))
solChangesplot=plot(2:convIt,hcat(ΔY[2:convIt],convCheck[2:convIt]),xlabel="Iteration",labels=["ΔY" "y-x"],xlims=(1,convIt))

setChangePlots=plot(activeSetPlot, setChangesPlot,layout=(2,1))
pubPlot(setChangePlots,thickscale=0.5,sizeWH=(1000,600),dpi=300)
if drawFig savefig(setChangePlots,path*"J_decentral_ALADIN_setPlots.png") end


fPlotalad=plot(1:convIt,dCMalad.obj[1:convIt,1],xlabel="Iteration",ylabel="obj function gap",xlims=(1,convIt),legend=false,yscale=:log10)
convItPlotalad=plot(1:convIt,dCMalad.lamIt[1:convIt,1],xlabel="Iteration",ylabel="2-Norm Lambda Gap",xlims=(1,convIt),legend=false) #,yscale=:log10
convPlotalad=plot(1:convIt,dCMalad.lam[1:convIt,1],xlabel="Iteration",ylabel="central lambda gap",xlims=(1,convIt),legend=false,yscale=:log10)
constPlotalad=plot(1:convIt,dCMalad.couplConst[1:convIt,1],xlabel="Iteration",ylabel="curr constraint Gap",xlims=(1,convIt),legend=false,yscale=:log10)

convPlots=plot(fPlotalad, convItPlotalad,convPlotalad, constPlotalad,layout=(2,2))
pubPlot(convPlots,thickscale=0.8,sizeWH=(1000,600),dpi=300)
if drawFig savefig(convPlots,path*"J_decentral_ALADIN_convPlots.png") end

checkDesiredStates(dLogalad.Sn[:,convIt],evS.Kn,evS.Snmin)
