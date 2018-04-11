#Micah Botkin-Levy
#4/10/18

println("Loading Packages...")

using Gadfly
using DataFrames
using JuMP
using Gurobi
using MAT #to read in scenarios from matlab


println("Reading in Data...")
function string_as_varname(s::String,v::Any)
	 s=Symbol(s)
	 @eval (($s) = ($v))
end

#read in mat scenario
path="C:\\Users\\micah\\Documents\\uvm\\Research\\EVC code\\N20\\EVCscenarioN20.mat"
vars = matread(path)
varnames=keys(vars)
varNum=length(varnames)
varKeys=collect(varnames)
varValues=collect(values(vars))

for i =1:varNum
	n=varKeys[i]
	v=varValues[i]
	if n in ["N" "K" "S"]
		v=convert(Int, v)
	end
	#if isa(v,Array)
	#	v=convert(DataFrame, v)
	#end
	string_as_varname(n,v)
end
println("done reading in")

Kn=convert(Array{Int,2},Kn)

#initialize
lambdaGuess=10
lambda0=ones(K+1,1)*lambdaGuess
lambda=lambda0
alpha=.05
convChk = 0.1
numIteration=2000
steps=K+1



#MPC here???
stepI = 1;

Lam=zeros((K+1),numIteration) #(rows are time, columns are iteration)
Lam[:,1]=lambda0
Xt=zeros((K+1),1) #rows are time
Xt[1,1]=T0
Xn=zeros((K+2),N) #row are time, column are EV
Xn[1,:]=s0'
Un=zeros((K+1),N) #row are time, column are EV

#iterate at each time step until convergence
for p=2:numIteration #change this to a while
#p=2
    @printf "iteration step %g of %g....\n" p numIteration

    #solve subproblem for each EV
    for evInd=1:N
        target=zeros((K+1),1)
        #target[max(1,Kn[evInd]-(stepI-1)*Ts):1:length(target),1]=s0[evInd] #fix Ts for time loop???
        target[Kn[evInd,1]:1:length(target),1]=s0[evInd,1]
        evM=Model(solver = GurobiSolver(OutputFlag=0))
        @variable(evM,un[1:K+1])
        @variable(evM,xn[1:K+1])
        objFun(x,u)=sum((x[k,1]-1)^2*Qsi[evInd,1] for k=1:K+1) +
        			sum((u[k,1])^2*Ri[evInd,1]    for k=1:K+1) +
                    sum(lambda[k,1]*u[k,1]      for k=1:K+1)
        @objective(evM,Min, objFun(xn,un))

        #@constraint(evM,(eye(K+1)-Ahats)*xn.==Ahats0*s0[evInd,1]+Bhats[evInd,1]*un) #this is slow???
		@constraint(evM,xn[1,1]==s0[evInd,1]) #fix for MPC loop
		@constraint(evM,[k=1:K],xn[k+1,1]==xn[k,1]+eta[evInd,1]*un[k,1]) #check K+1

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
            Xn[2:K+2,evInd]=getvalue(xn) #solved state goes in next time slot

            Un[:,evInd]=getvalue(un) #current go
        end
    end

    #solve coordinator problem
    coorM=Model(solver = GurobiSolver(OutputFlag=0))
    @variable(coorM,z[1:S*(K+1)])
    @variable(coorM,xt[1:K+1])
    @objective(coorM,Min,-sum(lambda[k,1]*sum(z[(k-1)*S+ii,1] for ii=1:S) for k=1:K+1))

    #@constraint(coorM,(eye(K+1)-AhatT)*xt.==AhatT0*T0+VhatT*w+EhatT*z)
	@constraint(coorM,xt[1,1]==T0) #fix for MPC loop
	@constraint(coorM,[k=1:K],xt[(k+1),1]==tau*xt[(k),1]+sum(Et*z[(k-1)*S+(1:1:S),1])+rho*w[k*2,1]) #check Z index

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
		 Xt=getvalue(xt);
	end

    #grad of lagragian
    gradL=sum(Un[:,ii] for ii=1:N)+Ghat*w-Fhat*getvalue(z); #need to sum U to get one value for each time step???

    #update lambda
    alpha_p = alpha/ceil(p/10);
    lambda_new=lambda+alpha_p*gradL;
    Lam[:,p]=lambda_new;
    lambda=lambda_new;

	#check convergence
	convGap = norm(Lam[:,p]-Lam[:,p-1],2)
	if(convGap < convChk )
		@printf "Converged after %g iterations\n" p
		break
	else
		@printf "convGap %f after %g iterations\n" convGap p
	end
end


p1=plot(Xn,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("SOC"))
display(p1)

p2=plot(Un,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV Current"))
display(p2)

p3=plot(x=1:K+1,y=Xt,Geom.line,
	Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)"),)
display(p3)
