#Micah Botkin-Levy
#4/10/18


datafile="n" #mat #"jld" #"n"

if datafile in ["mat" "jld"]; N=20 end

updateMethod="fastAscent" #dualAscent #fastAscent

println("Loading Packages...")

using Gadfly
using JuMP
using Gurobi
using Cairo #for png output
using Fontconfig
using Distributions

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



#initialize
xn0=s0
xt0=T0
d = Truncated(Normal(0), 0, 0.5)
lambda0=rand(d, K1+1)
#lambdaGuess=2
#lambda0=ones(K+1,1)*lambdaGuess
#lambda0=[0;linspace(7,0,K)]
#lambda0=rand(K+1)/2
lambda=lambda0

#convergence criteria
alpha=.2
convChk = 1e-6
numIteration=100
convIt=numIteration


#save matricies
Xtactual=zeros(K,1) #rows are time
Xn=zeros(K,N) #row are time, column are EV
Un=zeros(K,N) #row are time, column are EV
Lambda=zeros(K,1)
if updateMethod=="dualAscent"
	Xtmodel=zeros(K,1) #rows are time
end

for stepI=1:K
	#stepI = 1
	@printf "time step %g of %g....\n" stepI K

	horzLen=min(K1,K-stepI)

	Conv=zeros(numIteration,1)
	Lam=zeros(horzLen+1,numIteration) #(rows are time, columns are iteration)
	Lam[:,1]=lambda0
	Xtit=zeros(horzLen+1,1) #rows are time
	Tactual=zeros(horzLen+1,1) #rows are time
	Xnit=zeros(horzLen+1,N) #row are time, column are EV
	Unit=zeros(horzLen+1,N) #row are time, column are EV

	#iterate at each time step until convergence
	for p=2:numIteration

	    #@printf "iteration step %g of %g....\n" p numIteration

	    #solve subproblem for each EV
	    for evInd=1:N
	        target=zeros(horzLen+1,1)
	        #target[max(1,Kn[evInd]-(stepI-1)*Ts):1:length(target),1]=s0[evInd] #fix Ts for time loop???
	        target[Kn[evInd,1]:1:length(target),1]=xn0[evInd,1]
	        evM=Model(solver = GurobiSolver(OutputFlag=0))
	        @variable(evM,un[1:horzLen+1])
	        @variable(evM,xn[1:horzLen+1])
	        objFun(x,u)=sum((x[k,1]-1)^2*Qsi[evInd,1] for k=1:horzLen+1) +
	        			sum((u[k,1])^2*Ri[evInd,1]    for k=1:horzLen+1) +
	                    sum(lambda[k,1]*u[k,1]        for k=1:horzLen+1)
	        @objective(evM,Min, objFun(xn,un))
			@constraint(evM,xn[1,1]==xn0[evInd,1]+eta[evInd,1]*un[1,1]) #fix for MPC loop
			@constraint(evM,[k=1:horzLen],xn[k+1,1]==xn[k,1]+eta[evInd,1]*un[k+1,1]) #check K+1
	        @constraint(evM,xn.<=1)
	        @constraint(evM,xn.>=target)
	        @constraint(evM,un.<=imax[evInd,1])
	        @constraint(evM,un.>=imin[evInd,1])

			TT = STDOUT # save original STDOUT stream
			redirect_stdout()
	        status = solve(evM)
			redirect_stdout(TT)

			if status!=:Optimal
	            break
	        else
	            Xnit[:,evInd]=getvalue(xn) #solved state goes in next time slot
	            Unit[:,evInd]=getvalue(un) #current go
	        end
	    end


		if updateMethod=="dualAscent"
		    #solve coordinator problem
		    coorM=Model(solver = GurobiSolver(OutputFlag=0))
		    @variable(coorM,z[1:S*(horzLen+1)])
		    @variable(coorM,xt[1:horzLen+1])
		    @objective(coorM,Min,-sum(lambda[k,1]*sum(z[(k-1)*S+ii,1] for ii=1:S) for k=1:horzLen+1))
			@constraint(coorM,xt[1,1]==tau*xt0+gamma*deltaI*sum((2*m+1)*z[m+1,1] for m=0:S-1)+rho*w[stepI*2,1]) #fix for MPC loop
			@constraint(coorM,[k=1:horzLen],xt[k+1,1]==tau*xt[k,1]+gamma*deltaI*sum((2*m+1)*z[k*S+(m+1),1] for m=0:S-1)+rho*w[k*2+stepI*2,1])
		    @constraint(coorM,xt.<=Tmax)
		    @constraint(coorM,xt.>=0)
		    @constraint(coorM,z.<=deltaI)
		    @constraint(coorM,z.>=0)
			TT = STDOUT # save original STDOUT stream
			redirect_stdout()
		    status = solve(coorM)
			redirect_stdout(TT)
		    if status!=:Optimal
		        break
			else
				 Xtit=getvalue(xt);
			end

		    #grad of lagragian
			zSum=zeros(horzLen+1,1)
			for k=1:horzLen+1
				zSum[k,1]=sum(getvalue(z)[(k-1)*(S)+(1:S)])
			end
			gradL=sum(Unit[:,ii] for ii=1:N)+w[stepI*2-1:2:stepI*2+horzLen*2,1]-zSum
		end

		#calculate actual temperature from nonlinear model of XFRM
		ztotal=sum(Unit[:,ii] for ii=1:N) + w[stepI*2-1:2:stepI*2+horzLen*2,1]
		Tactual[1,1]=tau*xt0+gamma*ztotal[1,1]^2+rho*w[stepI*2,1]#fix for MPC loop
		for k=1:horzLen
			Tactual[k+1,1]=tau*Tactual[k,1]+gamma*ztotal[k+1,1]^2+rho*w[k*2+stepI*2,1] #check this????
		end

		if updateMethod=="fastAscent"
			#fast ascent
			gradL=Tactual-Tmax*ones(horzLen+1,1)
		end

	    #update lambda
	    alpha_p = alpha/ceil(p/2)
		#alpha_p = alpha/(p*5)

	    lambda_new=max.(lambda+alpha_p*gradL,0)
		#lambda_new=lambda+alpha_p*gradL
	    Lam[:,p]=lambda_new
	    lambda=lambda_new

		#check convergence
		convGap = norm(Lam[:,p]-Lam[:,p-1],2)
		Conv[p,1]=convGap
		if(convGap <= convChk )
			@printf "Converged after %g iterations\n" p
			convIt=p
			break
		else
			#@printf "convGap %e after %g iterations\n" convGap p
		end
	end


	#save states
	Xn[stepI,:]=Xnit[1,:]
	Un[stepI,:]=Unit[1,:]
	Xtactual[stepI,1]=Tactual[1,1]
	Lambda[stepI,1]=Lam[1,convIt]
	if updateMethod=="dualAscent"
		Xtmodel[stepI,1]=Xtit[1,1]
	end

	#apply current and set new states
	xn0=Xnit[1,:]
	xt0=Tactual[1,1]
	lambda0=Lam[1:horzLen+1,convIt]

end


println("plotting....")

pd1mpc=plot(Xn,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV SOC"),
		Coord.Cartesian(xmin=0,xmax=K+1),
		Theme(background_color=colorant"white",key_position = :none))
#display(pd1)

pd2mpc=plot(Un,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV Current"),
		Coord.Cartesian(xmin=0,xmax=K+1),
		Theme(background_color=colorant"white",key_position = :none))
#display(pd2)

pd3mpc=plot(layer(x=1:K+1,y=Xtactual,Geom.line,Theme(default_color=colorant"green")),
		yintercept=[Tmax],Geom.hline(color=["red"],style=:dot),
		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)"),
		Coord.Cartesian(xmin=0,xmax=K+1),Theme(background_color=colorant"white"))
if updateMethod=="dualAscent"
	push!(pd3,layer(x=1:K+1,y=Xtmodel,Geom.line,Theme(default_color=colorant"blue")))
	push!(pd3,Guide.manual_color_key("", ["PWL Temp", "Actual Temp"], ["blue", "green"]))
	push!(pd3,Theme(key_position = :top,background_color=colorant"white"))
end

#display(pd3)

pd4mpc=plot(x=1:K+1,y=Lambda,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel(raw"Lambda ($/A)"),
		Coord.Cartesian(xmin=0,xmax=K+1),Theme(background_color=colorant"white"))
#display(pd4)


#fName="J_Decentral_notfast.png"
fName="J_Decentral_fast.png"
# draw(PNG(path*fName, 13inch, 14inch), vstack(pd1,pd2,pd3,pd4))



lamPlot=plot(Lam[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel("Lambda"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white"))
convPlot=plot(x=1:convIt,y=Conv[1:convIt,1],Geom.line,
			Guide.xlabel("Iteration"), Guide.ylabel("Convergance Gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white"))
# draw(PNG(path*"distPlots_fast.png", 13inch, 8inch), vstack(lamPlot,convPlot))


#do this more elegantly
aggPlot=plot(x=1:K+1,y=sum(Un[:,i] for i=1:N),Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV Current"),
		Coord.Cartesian(xmin=0,xmax=K+1),
		Theme(background_color=colorant"white"))
