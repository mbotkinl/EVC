#remove converts once julia scenario


using Gadfly
using DataFrames
using JuMP
using Gurobi
#using Ipopt
#using Clp
using MAT #to read in scenarios from matlab


function string_as_varname(s::String,v::Any)
	 s=Symbol(s)
	 @eval (($s) = ($v))
end

#read in mat scenario
file = matopen("EVCscenarioN10.mat")
#file = matopen("test.mat")
varnames=names(file)
varNum=length(varnames)
varnamesS=collect(varnames)

#println(varnamesS)

for i =1:varNum
	n=varnamesS[i]
	#println(n)
	if n!="MCOS"
		string_as_varname(n,read(file,n))
		#if isa(n,Array)
		#	n=convert(DataFrame, n)
		#end
	end
end
close(file)
println("done reading in")

N=convert(Int,N)
K=convert(Int,K)
S=convert(Int,S)


#initialize
xi=x0

#setup model
#n1=convert(Int,N*(K+1))
#n2=convert(Int,(N+1)*(K+1))
#n3=convert(Int,S*(K+1))

println("setting up model")
m = Model(solver = GurobiSolver())
#m = Model(solver = ClpSolver())
#m = Model(solver = IpoptSolver())

@variable(m,u[1:N*(K+1)])
@variable(m,x[1:(N+1)*(K+1)])
@variable(m,z[1:S*(K+1)])

println("obj")
#@objective(m,Min,sum(u)) #fix???
#@objective(m,Min,transpose(u)*Rt*u+transpose(x)*Qt*x-2*ones(1,n2)*Qt*x)

#@objective(m,Min,sum(u'*Rt*u+x'*Qt*x-2*ones(1,n2)*Qt*x))
#@NLobjective(m,Min,(u'*Rt*u+x'*Qt*x-2*ones(1,n2)*Qt*x))


objFun(x,u)=sum(sum((x[(k-1)*(N+1)+n,1]-1)^2*Qsi[n] for n=1:N) for k=1:K)+sum(sum((u[(k-1)*N+n,1])^2*Ri[n] for n=1:N) for k=1:(K-1))

@objective(m,Min, objFun(x,u))

println("constraints")
@constraint(m,(eye(n2)-Ahat)*x.==A0hat*xi+Bhat*u+Vhat*w+Ehat*z)
@constraint(m,0.==Hhat*u+Ghat*w+Fhat*z)
@constraint(m,x.<=repmat([ones(N,1);Tmax],K+1,1))
@constraint(m,x.>=0) #add target???
@constraint(m,u.<=repmat(imax,K+1,1))
@constraint(m,u.>=repmat(imin,K+1,1))
@constraint(m,z.>=0)
@constraint(m,z.<=deltaI)

println("solving....")
status = solve(m)



xRaw=getvalue(x)
xPlot=DataFrame()
for ii= 1:N
	xPlot[:,ii]=xRaw[ii:N+1:length(xRaw)]
end

#p1=
plot(xPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("SOC"))
#gui()
#img=SVG("SOC.svg",6inch,4inch)
#draw(img,p1)


uRaw=getvalue(u)
uPlot=DataFrame()
for ii= 1:N
	uPlot[:,ii]=uRaw[ii:N:length(uRaw)]
end

#p2=
plot(uPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV Current"))

#img=SVG("Current.svg",6inch,4inch)
#draw(img,p2)



#getvalue(u)

#pwd()
#cd(homedir())
#cd("Documents\\uvm\\Research\\EVC code\\Julia")
#include("EVC_central_J.jl") to run
