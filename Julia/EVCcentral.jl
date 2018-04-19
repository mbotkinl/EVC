#Micah Botkin-Levy
#4/8/18
N=20
datafile="jld" #mat #"jld" #"n"

println("Loading Packages...")

using Gadfly
using JuMP
using Gurobi
using Cairo #for png output
using Fontconfig

if datafile=="mat"
	using MAT #to read in scenarios from matlab
elseif datafile=="jld"
	using JLD
end

if datafile in ["mat" "jld" "n"]
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
		#if isa(v,Array)
		#	v=convert(DataFrame, v)
		#end
		string_as_varname(n,v)
	end
	println("done reading in")

	if datafile=="mat"
		Kn=convert(Array{Int,2},Kn)
	end
end


#initialize
xi=x0

#add mpc loop here ??
stepI=1

println("setting up model")
m = Model(solver = GurobiSolver())
#m = Model(solver = ClpSolver())
#m = Model(solver = IpoptSolver())

#u w and z are one index ahead of x. i.e the x[k+1]=x[k]+eta*u[k+1]
@variable(m,u[1:N*(K+1)])
@variable(m,x[1:(N+1)*(K+1)])
@variable(m,z[1:S*(K+1)])

#desired SOC
target=zeros((N+1)*(K+1),1);
for ii=1:N
   cur=Kn[ii]-(stepI-1);
   ind=max(0,(cur-1)*(N+1))+ii:N+1:length(target);
   target[ind]=Sn[ii,1];
end

println("obj")
objFun(x,u)=sum(sum((x[(k-1)*(N+1)+n,1]-1)^2*Qsi[n,1] for n=1:N+1) for k=1:K+1) +
			sum(sum((u[(k-1)*N+n,1])^2*Ri[n,1]        for n=1:N)   for k=1:K+1)
@objective(m,Min, objFun(x,u))

println("constraints")
@constraint(m,stateCon1,x[1:N,1].==xi[1:N,1]+eta[:,1].*u[1:N,1])
@constraint(m,stateCon2[k=1:K,n=1:N],x[n+(k)*(N+1),1]==x[n+(k-1)*(N+1),1]+eta[n,1]*u[n+(k)*(N),1])
@constraint(m,tempCon1,x[N+1,1]==tau*T0+sum(Et*z[1:S,1])+rho*w[2,1]) #fix for MPC loop???
@constraint(m,tempCon2[k=1:K],x[(N+1)*(k+1),1]==tau*x[(N+1)*(k),1]+sum(Et*z[(k)*S+(1:1:S),1])+rho*w[k*2+2,1])
@constraint(m,currCon[k=1:K+1],0==-sum(u[(k-1)*(N)+(1:N)])-w[(k-1)*2+1]+sum(z[(k-1)*(S)+(1:S)]))
@constraint(m,upperTCon,x.<=repmat([ones(N,1);Tmax],K+1,1))
@constraint(m,x.>=target)
#@constraint(m,x.>=0)
@constraint(m,upperCCon,u.<=repmat(imax,K+1,1))
@constraint(m,u.>=repmat(imin,K+1,1))
@constraint(m,z.>=0)
@constraint(m,z.<=deltaI)

println("solving....")
status = solve(m)
if status!=:Optimal
    return
end

uRaw=getvalue(u)
xRaw=getvalue(x)

#calculate actual temp
Tactual=zeros(K+1,1)
ztotal=zeros(K+1,1)
for k=1:K+1
	ztotal[k,1]=sum((uRaw[(k-1)*N+n,1]) for n=1:N) + w[(k-1)*2+1,1]
end
Tactual[1,1]=tau*T0+gamma*ztotal[1,1]^2+rho*w[2,1]#fix for MPC loop
for k=1:K
	Tactual[k+1,1]=tau*Tactual[k,1]+gamma*ztotal[k+1,1]^2+rho*w[k*2+2,1]
end

#getobjectivevalue(m)
#lambdaCurr=getdual(currCon)
#lambdaTemp=[getdual(tempCon1);getdual(tempCon2)]
#lambdaState=[getdual(stateCon1)';getdual(stateCon2)]
lambdaUpperT=-getdual(upperTCon)
#lambdaUpperC=-getdual(upperCCon)

println("plotting....")
xPlot=zeros(K+1,N)
for ii= 1:N
	xPlot[:,ii]=xRaw[collect(ii:N+1:length(xRaw))]
end

p1=plot(xPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV SOC"),
		Coord.Cartesian(xmin=0,xmax=K+1),
		Theme(background_color=colorant"white",key_position = :none))
#display(p1)

uPlot=zeros(K+1,N)
for ii= 1:N
	uPlot[:,ii]=uRaw[collect(ii:N:length(uRaw))]
end

p2=plot(uPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV Current"),
		Coord.Cartesian(xmin=0,xmax=K+1),
		Theme(background_color=colorant"white",key_position = :none))
#display(p2)

p3=plot(layer(x=1:K+1,y=xRaw[N+1:N+1:length(xRaw)],Geom.line,Theme(default_color=colorant"blue")),
		layer(x=1:K+1,y=Tactual,Geom.line,Theme(default_color=colorant"green")),
		yintercept=[Tmax],Geom.hline(color=["red"],style=:dot),
		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=K+1),Theme(background_color=colorant"white",key_position = :top),
		Guide.manual_color_key("", ["PWL Temp", "Actual Temp"], ["blue", "green"]))
#display(p3)

p4=plot(x=1:K+1,y=lambdaUpperT[collect(N+1:N+1:length(lambdaUpperT))],Geom.line,
		Guide.xlabel("Time"), Guide.ylabel(raw"Lambda ($/A)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=K+1),Theme(background_color=colorant"white"))
#display(p4)

fName="J_Central.png"

#draw(PNG(path*fName, 13inch, 14inch), vstack(p1,p2,p3,p4))
