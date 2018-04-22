#Micah Botkin-Levy
#4/8/18
datafile="n" #mat #"jld" #"n"

if datafile in ["mat" "jld"]; N=20 end

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

#initialize with current states
xn0=s0
xt0=T0

#save matricies
Xtmodel=zeros(K,1) #rows are time
Xtactual=zeros(K,1) #rows are time
Xn=zeros(K,N) #row are time, column are EV
Un=zeros(K,N) #row are time, column are EV
Lambda=zeros(K,1)

for stepI=1:K
	#stepI=1

	@printf "time step %g of %g....\n" stepI K

	horzLen=min(K1,K-stepI)

	#println("setting up model")
	m = Model(solver = GurobiSolver(OutputFlag=0))

	#u w and z are one index ahead of x. i.e the x[k+1]=x[k]+eta*u[k+1]
	@variable(m,u[1:N*(horzLen+1)])
	@variable(m,xn[1:(N)*(horzLen+1)])
	@variable(m,xt[1:(horzLen+1)])
	@variable(m,z[1:S*(horzLen+1)])

	#desired SOC
	target=zeros(N*(horzLen+1),1);
	for ii=1:N
	   cur=Kn[ii]-(stepI-1)
	   ind=max(0,(cur-1)*N)+ii:N:length(target)
	   target[ind]=Sn[ii,1]
	end

	#println("obj")
	objFun(xn,xt,u)=sum(sum((xn[(k-1)*(N)+n,1]-1)^2*Qsi[n,1]     for n=1:N) for k=1:horzLen+1) +
					sum((xt[k,1]-1)^2*Qsi[N+1,1]                 for k=1:horzLen+1) +
					sum(sum((u[(k-1)*N+n,1])^2*Ri[n,1]           for n=1:N) for k=1:horzLen+1)
	@objective(m,Min, objFun(xn,xt,u))

	#println("constraints")
	@constraint(m,stateCon1,xn[1:N,1].==xn0[1:N,1]+eta[:,1].*u[1:N,1])
	@constraint(m,stateCon2[k=1:horzLen,n=1:N],xn[n+(k)*(N),1]==xn[n+(k-1)*(N),1]+eta[n,1]*u[n+(k)*(N),1])
	@constraint(m,tempCon1,xt[1,1]==tau*xt0+gamma*deltaI*sum((2*m+1)*z[m+1,1] for m=0:S-1)+rho*w[stepI*2,1])
	@constraint(m,tempCon2[k=1:horzLen],xt[k+1,1]==tau*xt[k,1]+gamma*deltaI*sum((2*m+1)*z[k*S+(m+1),1] for m=0:S-1)+rho*w[stepI*2+k*2,1])
	@constraint(m,currCon[k=1:horzLen+1],0==-sum(u[(k-1)*(N)+(1:N)])-w[(k-1)*2+(stepI*2-1)]+sum(z[(k-1)*(S)+(1:S)]))
	@constraint(m,xn.<=1)
	@constraint(m,xn.>=target)
	@constraint(m,upperTCon,xt.<=Tmax)
	@constraint(m,xt.>=0)
	@constraint(m,upperCCon,u.<=repmat(imax,horzLen+1,1))
	@constraint(m,u.>=repmat(imin,horzLen+1,1))
	@constraint(m,z.>=0)
	@constraint(m,z.<=deltaI)

	#println("solving....")

	TT = STDOUT # save original STDOUT stream
	redirect_stdout()
	status = solve(m)
	redirect_stdout(TT)
	if status!=:Optimal
	    return
	end

	#getobjectivevalue(m)
	#lambdaCurr=getdual(currCon)
	#lambdaTemp=[getdual(tempCon1);getdual(tempCon2)]
	#lambdaState=[getdual(stateCon1)';getdual(stateCon2)]
	lambdaUpperT=-getdual(upperTCon)
	#lambdaUpperC=-getdual(upperCCon)

	uRaw=getvalue(u)
	xnRaw=getvalue(xn)
	xtRaw=getvalue(xt)

	#calculate actual temp
	Tactual=zeros(horzLen+1,1)
	ztotal=zeros(horzLen+1,1)
	for k=1:horzLen+1
		ztotal[k,1]=sum((uRaw[(k-1)*N+n,1]) for n=1:N) + w[(k-1)*2+(stepI*2-1),1]
	end
	Tactual[1,1]=tau*xt0+gamma*ztotal[1,1]^2+rho*w[stepI*2,1]

	#do we need this for anything???
	for k=1:horzLen
		Tactual[k+1,1]=tau*Tactual[k,1]+gamma*ztotal[k+1,1]^2+rho*w[k*2+stepI*2,1]
	end


	#save states
	Xtmodel[stepI,1]=xtRaw[1,1]
	Xtactual[stepI,1]=Tactual[1,1]
	Xn[stepI,:]=xnRaw[1:N,1]
	Un[stepI,:]=uRaw[1:N,1]
	Lambda[stepI,1]=lambdaUpperT[1,1]

	#apply currents and assume new states
	xn0=xnRaw[1:N,1]
	xt0=Tactual[1,1]

end


println("plotting....")
# xPlot=zeros(K+1,N)
# for ii= 1:N
# 	xPlot[:,ii]=xRaw[collect(ii:N+1:length(xRaw))]
# end

p1mpc=plot(Xn,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV SOC"),
		Coord.Cartesian(xmin=0,xmax=K+1),
		Theme(background_color=colorant"white",key_position = :none))
#display(p1)

# uPlot=zeros(K+1,N)
# for ii= 1:N
# 	uPlot[:,ii]=uRaw[collect(ii:N:length(uRaw))]
# end

p2mpc=plot(Un,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV Current"),
		Coord.Cartesian(xmin=0,xmax=K+1),
		Theme(background_color=colorant"white",key_position = :none))
#display(p2)

p3mpc=plot(layer(x=1:K,y=Xtmodel,Geom.line,Theme(default_color=colorant"blue")),
		layer(x=1:K,y=Xtactual,Geom.line,Theme(default_color=colorant"green")),
		yintercept=[Tmax],Geom.hline(color=["red"],style=:dot),
		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=K+1),Theme(background_color=colorant"white",key_position = :top),
		Guide.manual_color_key("", ["PWL Temp", "Actual Temp"], ["blue", "green"]))
#display(p3)

p4mpc=plot(x=1:K,y=Lambda,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel(raw"Lambda ($/A)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=K+1),Theme(background_color=colorant"white"))
#display(p4)

fName="J_Central.png"
#draw(PNG(path*fName, 13inch, 14inch), vstack(p1,p2,p3,p4))

#lambda0=lambdaUpperT[collect(N+1:N+1:length(lambdaUpperT))]


#do this more elegantly
aggPmpc=plot(x=1:K,y=sum(Un[:,i] for i=1:N),Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV Current"),
		Coord.Cartesian(xmin=0,xmax=K+1),
		Theme(background_color=colorant"white"))
