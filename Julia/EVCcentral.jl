#Micah Botkin-Levy
#4/8/18

println("Loading Packages...")

using Gadfly
using DataFrames
using JuMP
using Gurobi
#using Ipopt
#using Clp
using MAT #to read in scenarios from matlab
using Cairo #for png output
using Fontconfig



println("Reading in Data...")

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
		#	n=convert(DataFrame, n) #fix this???
		#end
	end
end
close(file)
println("done reading in")

N=convert(Int,N)
K=convert(Int,K)
S=convert(Int,S)
Kn=convert(Array{Int,2},Kn)

#initialize
xi=x0

#add mpc loop here ??
stepI=1

println("setting up model")
m = Model(solver = GurobiSolver())
#m = Model(solver = ClpSolver())
#m = Model(solver = IpoptSolver())

@variable(m,u[1:N*(K+1)])
@variable(m,x[1:(N+1)*(K+1)])
@variable(m,z[1:S*(K+1)])

#desired SOC
target=zeros((N+1)*(K+1),1);
for ii=1:N
   cur=Kn[ii]-(stepI-1);
   ind=max(0,(cur-1)*(N+1))+ii:N+1:length(target);
   target[ind]=Sn[ii];
end

println("obj")
#@objective(m,Min,sum(u'*Rt*u+x'*Qt*x-2*ones(1,n2)*Qt*x))
objFun(x,u)=sum(sum((x[(k-1)*(N+1)+n,1]-1)^2*Qsi[n] for n=1:N+1) for k=1:K+1)+
			sum(sum((u[(k-1)*N+n,1])^2*Ri[n]        for n=1:N) for k=1:K+1)
@objective(m,Min, objFun(x,u))

println("constraints")
@constraint(m,(eye((N+1)*(K+1))-Ahat)*x.==A0hat*xi+Bhat*u+Vhat*w+Ehat*z) #this is slow???
@constraint(m,0.==Hhat*u+Ghat*w+Fhat*z)
@constraint(m,x.<=repmat([ones(N,1);Tmax],K+1,1))
@constraint(m,x.>=target)
@constraint(m,u.<=repmat(imax,K+1,1))
@constraint(m,u.>=repmat(imin,K+1,1))
@constraint(m,z.>=0)
@constraint(m,z.<=deltaI)

println("solving....")
status = solve(m)


println("plotting....")
xRaw=getvalue(x)
xPlot=DataFrame()
for ii= 1:N
	xPlot[:,ii]=xRaw[collect(ii:N+1:length(xRaw))]
end

p1=plot(xPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("SOC"))
display(p1)
#draw(PNG("SOC.png", 4inch, 3inch), p1)


uRaw=getvalue(u)
uPlot=DataFrame()
for ii= 1:N
	uPlot[:,ii]=uRaw[collect(ii:N:length(uRaw))]
end

p2=plot(uPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV Current"))
display(p2)
#draw(PNG("Current.png", 4inch, 3inch), p2)



p3=plot(x=1:K+1,y=xRaw[N+1:N+1:length(xRaw)],Geom.line,
	Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)"),)
display(p3)
#draw(PNG("Temp.png", 4inch, 3inch), p3)
