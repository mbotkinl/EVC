using Distributions
using Gurobi

#pull out a few key variables
N=evS.N
S=evS.S

include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCpem.jl")

fname = "d_PEM_N$(N)"

if loadResults
	println("Reading in PEM Sim")
	loadF=load(path*fname*".jld2")
	evS=loadF["scenario"]
	pemSol=loadF["solution"]
else
	println("Running PEM Sim")
	timeT=@elapsed pemSol=pemEVC(evS,slack,silent)
	if saveResults saveRun(path,fname,timeT, evS,pemSol) end
end

println("plotting....")

p1pem=plot(pemSol.Sn,xlabel="Time",ylabel="PEV SOC",legend=false,xlims=(0,evS.K),ylims=(0,1))
if drawFig savefig(p1pem,path*"J_PEM_SOC.png") end

p2pem=plot(pemSol.Un,xlabel="Time",ylabel="PEV Current (kA)",legend=false,xlims=(0,evS.K))
if drawFig savefig(p2pem,path*"J_PEM_Curr.png") end

p3pem=plot(pemSol.Tactual*1000,label="XFRM Temp",xlims=(0,evS.K),xlabel="Time",ylabel="Temp (K)")
plot!(p3pem,evS.Tmax*ones(evS.K)*1000,label="XFRM Limit",line=(:dash,:red))
if drawFig savefig(p3pem,path*"J_PEM_Temp.png") end

# pRpem=plot(ratio,xlabel="Time",ylabel="Ratio",xlims=(0,horzLen+1),legend=false)
# if drawFig savefig(pRpem,path*"J_PEM_Ratio.png") end


aggUpem=plot(hcat(cSol.uSum,pemSol.uSum),label=["Central" "PEM"],
			xlims=(0,evS.K),xlabel="Time",ylabel="PEV Current (kA)")
if drawFig savefig(aggUpem,path*"J_PEM_uSum.png") end


checkDesiredStates(pemSol.Sn,evS.Kn,evS.Snmin)
checkDesiredStates(pemSol.Sn,evS.Kn,ones(N))
