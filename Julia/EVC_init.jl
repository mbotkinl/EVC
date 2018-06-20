
N=5
datafile="n" #"mat" #"jld" #"n"
updateMethod="dualAscent" #dualAscent #fastAscent
drawFig=0
noTlimit=0
#if datafile in ["mat" "jld"]; N=5 end

println("Loading Packages...")

using JuMP
using Gadfly
using Cairo #for png output
using Fontconfig

if datafile=="mat"
	using MAT #to read in scenarios from matlab
elseif datafile=="jld"
	using JLD
else #create scenatio
	using Distributions
end

if datafile in ["mat" "jld"]
	println("Reading in Data...")

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
else
	println("Creating EV Scenario...")
	include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCscenario.jl")
	evS=setupScenario(5)
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
