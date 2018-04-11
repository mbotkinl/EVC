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
file = matopen("EVCscenarioN10.mat")
varnames=names(file)
varNum=length(varnames)
varnamesS=collect(varnames)

for i =1:varNum
	n=varnamesS[i]
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
lambdaGuess=10
lambda0=ones(K+1,1)*lambdaGuess
lambda=lambda0
alpha=.05
convChk = 1
numIteration=300
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
        target[Kn[evInd]:1:length(target),1]=s0[evInd]
        evM=Model(solver = GurobiSolver(OutputFlag=0))
        @variable(evM,un[1:K+1])
        @variable(evM,xn[1:K+1])
        objFun(x,u)=sum((x[k,1]-1)^2*Qsi[evInd] for k=1:K+1) +
        			sum((u[k,1])^2*Ri[evInd]    for k=1:K+1) +
                    sum(lambda[k,1]*u[k,1]      for k=1:K+1)
        @objective(evM,Min, objFun(xn,un))
        @constraint(evM,(eye(K+1)-Ahats)*xn.==Ahats0*s0[evInd]+Bhats[evInd,1]*un) #this is slow???
        @constraint(evM,xn.<=1)
        @constraint(evM,xn.>=target)
        @constraint(evM,un.<=imax[evInd])
        @constraint(evM,un.>=imin[evInd])
        status = solve(evM)
        if status!=:Optimal
            #break
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
    @constraint(coorM,(eye(K+1)-AhatT)*xt.==AhatT0*T0+VhatT*w+EhatT*z)
    @constraint(coorM,xt.<=Tmax)
    @constraint(coorM,xt.>=0)
    @constraint(coorM,z.<=deltaI)
    @constraint(coorM,z.>=0)
    status = solve(coorM)
    if status!=:Optimal
        #break
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
