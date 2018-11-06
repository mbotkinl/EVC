
#EVC with ALADIN for PWL convex relaxation
#Micah Botkin-Levy
#5/22/18

#formulation try 2
#u, sn, xt, and z are all in "y", v is the auxilliary variable corresponding to x in literature
#current constraint is coupling

relaxString= if relaxed==true "_relax"else "" end
fname = "dALADIN_NL_N$(N)"*relaxString

if loadResults
	println("Reading in NL AlAD Sim")
	loadF=load(path*fname*".jld2")
	evS=loadF["scenario"]
	dLognlalad=loadF["solution"]
	dCMnlalad=loadF["convMetrics"]
	convIt=loadF["convIt"]
else
	println("Running NL AlAD Sim")
	timeT=@elapsed dLognlalad,dCMnlalad,convIt,ΔY,convCheck=nlEValad(N,S,horzLen,maxIt,evS,cSolnl,relaxed,slack)
	if saveResults==1 saveRun(path,fname,timeT, evS,dLognlalad, dCMnlalad, convIt) end
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

pd3aladNL=plot(1:horzLen+1,dLognlalad.Xt[:,convIt],label="XFRM Temp",xlims=(0,horzLen+1),xlabel="Time",ylabel="Temp (K)")
plot!(pd3aladNL,1:horzLen+1,evS.Tmax*ones(horzLen+1),label="XFRM Limit",line=(:dash,:red))
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



activeSet=zeros(convIt,1)
setChanges=zeros(convIt,1)
for ii=2:convIt
    activeSet[ii,1]=sum(abs.(dLognlalad.Cs[:,ii]))+sum(abs.(dLognlalad.Ct[:,ii]))+
              sum(abs.(dLognlalad.Cu[:,ii]))+sum(abs.(dLognlalad.Ci[:,ii]))
    setChanges[ii,1]=sum(abs.(dLognlalad.Cs[:,ii]-dLognlalad.Cs[:,ii-1]))+sum(abs.(dLognlalad.Ct[:,ii]-dLognlalad.Ct[:,ii-1]))+
                     sum(abs.(dLognlalad.Cu[:,ii]-dLognlalad.Cu[:,ii-1]))+sum(abs.(dLognlalad.Ci[:,ii]-dLognlalad.Ci[:,ii-1]))
end

activeSetPlot=plot(2:convIt,activeSet[2:convIt],xlabel="Iteration",ylabel="Total Active inequality constraints",
                   legend=false,xlims=(2,convIt))
setChangesPlot=plot(2:convIt,setChanges[2:convIt],xlabel="Iteration",ylabel="Changes in Active inequality constraints",
                  legend=false,xlims=(2,convIt))
solChangesplot=plot(2:convIt,hcat(ΔY[1:convIt],convCheck[1:convIt]),xlabel="Iteration",labels=["ΔY" "y-x"],xlims=(1,convIt))

fPlotaladNL=plot(1:convIt,dCMnlalad.obj[1:convIt,1],xlabel="Iteration",ylabel="obj function gap",xlims=(1,convIt),legend=false,yscale=:log10)
convItPlotaladNL=plot(1:convIt,dCMnlalad.lamIt[1:convIt,1],xlabel="Iteration",ylabel="2-Norm Lambda Gap",xlims=(1,convIt),legend=false,yscale=:log10)
convPlotaladNL=plot(1:convIt,dCMnlalad.lam[1:convIt,1],xlabel="Iteration",ylabel="central lambda gap",xlims=(1,convIt),legend=false,yscale=:log10)
constPlotaladNL=plot(1:convIt,dCMnlalad.couplConst[1:convIt,1],xlabel="Iteration",ylabel="curr constraint Gap",xlims=(1,convIt),legend=false,yscale=:log10)

checkDesiredStates(dLognlalad.Sn[:,convIt],evS.Kn,evS.Snmin)
