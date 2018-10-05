#C:\Users\micah\AppData\Local\JuliaPro-0.6.2.2\Julia-0.6.2\bin\julia
# C:\Users\micah\AppData\Local\Julia-0.7.0\bin\julia

include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//structEVC.jl")
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVChelpers.jl")

N=50
path="C:\\Users\\micah\\Documents\\uvm\\Research\\Results\\PWLvNL\\N$(N)\\"
datafile="n" #"mat" #"jld" #"n"
file="EVCscenarioN$(N)."*datafile

updateMethod="dualAscent" #dualAscent #fastAscent
drawFig=0
saveResults=0
loadResults=0
noTlimit=false
maxIt=100
relaxed=false
verbose=false

println("Loading Packages...")

using JuMP
using Gadfly

# these cause pangolayout error
# using Cairo #for png output
# using Fontconfig

using Parameters
using Printf
using Distributed

if datafile=="mat"
	using MAT #to read in scenarios from matlab

	function string_as_varname(s::String,v::Any)
		 s=Symbol(s)
		 @eval (($s) = ($v))
	end

	#read in mat scenario
	vars = matread(path*file)

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
	using JLD2
	println("Reading in Data...")
	loadF=@load(path*file)
	evS=loadF["evScenario"]
else #create scenario

	Tmax=400
	Dload_amplitude=0
	saveS=false

	include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCscenario.jl")
    if saveS==true using JLD2 end
	using Distributions
	println("Creating EV Scenario...")
	evS=setupScenario(N;Tmax=Tmax,Dload_amplitude=Dload_amplitude,saveS=saveS)
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
