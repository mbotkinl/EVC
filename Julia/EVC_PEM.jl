using Distributions
using Gurobi

#pull out a few key variables
N=evS.N
S=evS.S
horzLen=evS.K1

include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCpem.jl")

fname = "d_PEM_N$(N)"

if loadResults
	println("Reading in PEM Sim")
	loadF=load(path*fname*".jld2")
	evS=loadF["scenario"]
	pemSol=loadF["solution"]
else
	println("Running PEM Sim")
	timeT=@elapsed pemSol,ratio=pemEVC(N,S,horzLen,evS,slack)
	timeT=timeT/horzLen
	if saveResults saveRun(path,fname,timeT, evS,pemSol) end
end

p1pem=plot(pemSol.Sn,xlabel="Time",ylabel="PEV SOC",legend=false,xlims=(0,horzLen+1),ylims=(0,1))
if drawFig savefig(p1pem,path*"J_PEM_SOC.png") end

p2pem=plot(pemSol.Un,xlabel="Time",ylabel="PEV Current (kA)",legend=false,xlims=(0,horzLen+1))
if drawFig savefig(p2pem,path*"J_PEM_Curr.png") end

p3pem=plot(1:horzLen+1,pemSol.Xt,label="XFRM Temp",xlims=(0,horzLen+1),xlabel="Time",ylabel="Temp (K)")
plot!(p3pem,1:horzLen+1,evS.Tmax*ones(horzLen+1),label="XFRM Limit",line=(:dash,:red))
if drawFig savefig(p3pem,path*"J_PEM_Temp.png") end

# pRpem=plot(ratio,xlabel="Time",ylabel="Ratio",xlims=(0,horzLen+1),legend=false)
# if drawFig savefig(pRpem,path*"J_PEM_Ratio.png") end


aggUpem=plot(1:horzLen+1,hcat(sum(uPlot[:,i] for i=1:N),pemSol.uSum),label=["Central" "PEM"],
			xlims=(0,horzLen+1),xlabel="Time",ylabel="PEV Current (kA)")
if drawFig savefig(aggUpem,path*"J_PEM_uSum.png") end


checkDesiredStates(pemSol.Sn,evS.Kn,evS.Snmin)
checkDesiredStates(pemSol.Sn,evS.Kn,ones(N))
