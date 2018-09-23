#C:\Users\micah\AppData\Local\JuliaPro-0.6.2.2\Julia-0.6.2\bin\julia

include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//structEVC.jl")

N=10
datafile="jld" #"mat" #"jld" #"n"
updateMethod="dualAscent" #dualAscent #fastAscent
drawFig=0
noTlimit=0
maxIt=10

Tmax=371
Dload_amplitude=0
saveS=false

println("Loading Packages...")

using JuMP
using Gadfly
using Cairo #for png output
using Fontconfig
using Parameters

if datafile=="mat"
	using MAT #to read in scenarios from matlab

	function string_as_varname(s::String,v::Any)
		 s=Symbol(s)
		 @eval (($s) = ($v))
	end

	#read in mat scenario
	path="C:\\Users\\micah\\Documents\\uvm\\Research\\EVC code\\N$(N)\\"
	file="EVCscenarioN$(N)."*datafile
	if datafile=="mat"
		vars = matread(path*file)
	elseif datafile=="jld"
		vars=load(path*file)
	end
	varnames=keys(vars)
	varNum=length(varnames)
	varKeys=collect(varnames)
	varValues=collect(values(vars))

	for i =1:varNum
		n=varKeys[i]
		v=varValues[i]
		if datafile=="mat"
			if n in ["N" "K" "S"]
				v=convert(Int, v)
			end
		end
		string_as_varname(n,v)
	end
	println("done reading in")

	if datafile=="mat"
		Kn=convert(Array{Int,2},Kn)
	end
elseif datafile=="jld"
	using JLD
	println("Reading in Data...")
	path="C:\\Users\\micah\\Documents\\uvm\\Research\\EVC code\\Julia\\PWLvNL\\"
	file="EVCscenarioN$(N)."*datafile
	loadF=JLD.load(path*file)
	evS=loadF["evScenario"]
else #create scenatio
	include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCscenario.jl")
    if saveS==true using JLD end
	using Distributions
	println("Creating EV Scenario...")
	evS=setupScenario(N;Tmax=Tmax,Dload_amplitude=Dload_amplitude,saveS=saveS)
end



#temp backwards capability
if isassigned(evS.ηP)==false
	evS.ηP=eta
	evS.γP=gamma
	evS.ρP=rho
	evS.τP=tau
end

if isassigned(evS.Snmin)==false
	evS.Snmin=Sn
end


# @time runALADit(1)
#@time testALAD(1)

# Profile.clear()
# runALADit(1)
# @profile runALADit(1)
# Juno.profiler()


# Profile.clear()
# testDual(1)
# @profile testDual(1)
# Juno.profiler()
