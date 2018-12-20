
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funHub.jl")
fname = "dADMM_H$(H)"

if loadResults
	println("Reading in ADMM Hub Sim")
	loadF=load(path*fname*".jld2")
	hubS=loadF["scenario"]
	dSoladmm=loadF["solution"]
else
	t0=hubS.t0
	e0=hubS.e0
	prevLam=ones(hubS.K1+1,1)
	prevVz=-hubS.deltaI*ones((hubS.K1+1),hubS.S)
	#prevVu=repeat(maximum(hubS.uMax,dims=1),outer=[(hubS.K1+1),1])
	prevVu=.01*ones((hubS.K1+1),hubS.H)
	ρADMMp = 1e3

	println("Running ADMM Hub Sim")
	timeT=@elapsed dSoladmm=hubADMM(maxIt,hubS,cSol,mode,silent)
	if saveResults saveRun(path,fname,timeT, evS, dSoladmm) end
end
