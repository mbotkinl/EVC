# test of ALADIN with simple example
#Micah Botkin-Levy
#5/12/18

using JuMP
using Gurobi
using Gadfly

#minimize sum_i=1^N x_i^2*q_ii
#s.t. sum_i=1^N a_i*x_i=b

#set up Problem
N=1000
Q=ones(N,1) #change to different weights??
A=round.(rand(N),4)
b=round.(rand(1),4)

#central solution
m=Model(solver = GurobiSolver(OutputFlag=0))
@variable(m,xC[1:N])
@objective(m,Min,sum(0.5*xC[i,1]^2*Q[i,1] for i=1:N))
@constraint(m,lambdaP,sum(A[i,1]*xC[i,1] for i=1:N).==b)
status = solve(m)
if status!=:Optimal
    return
end

#ALADIN
#initialize
maxIt=50
convIt=maxIt
X=zeros(N,maxIt) #(rows are variable, columns are iteration)
X0=ones(N,1)
x=X0
Lambda=zeros(maxIt,1)
lambda0=0
lambda=lambda0
H=Q+round.(rand(N)/10,4)
epsilon=1e-6
avgN=zeros(maxIt,1)

for k=1:maxIt
    @printf "iteration step %g of %g....\n" k maxIt
    G=zeros(N,1)
    Y=zeros(N,1)
    #solve NLP
    for i=1:N
        m=Model(solver = GurobiSolver(OutputFlag=0))
        @variable(m,y)
        @objective(m,Min,0.5*y^2*Q[i,1]+lambda*A[i,1]*y+0.5*(y-x[i,1])^2*H[i,1])
        TT = STDOUT # save original STDOUT stream
		redirect_stdout()
        status = solve(m)
		redirect_stdout(TT)
        if status!=:Optimal
            return
        else
            Y[i,1]=getvalue(y)
            G[i,1]=H[i,1]*(x[i,1]-Y[i,1])-A[i,1]*lambda
        end
    end

    #update H here


    avgN[k,1]=mean(norm.(x-Y))
    if all(norm.(x-Y).<=epsilon)
        @printf "Converged after %g iterations\n" k
        convIt=k
        return
    end

    #coupled QP
    m=Model(solver = GurobiSolver(OutputFlag=0))
    @variable(m,dy[1:N])
    @objective(m,Min,sum(0.5*dy[i,1]^2*H[i,1]+G[i,1]*dy[i,1] for i=1:N))
    @constraint(m,lambdaP,sum(A[i,1]*(Y[i,1]+dy[i,1]) for i=1:N).==b)
    TT = STDOUT # save original STDOUT stream
    redirect_stdout()
    status = solve(m)
    redirect_stdout(TT)
    if status!=:Optimal
        return
    end

    #update
    x=Y+getvalue(dy)
    X[:,k]=x
    lambda=-getdual(lambdaP[1])
    Lambda[k,1]=lambda
end


compPlot=plot(layer(x=1:N,y=getvalue(xC),Geom.point),
              layer(x=1:N,y=X[:,convIt-1],Geom.point))

lamPlot=plot(x=1:(convIt-1),y=Lambda[1:convIt-1,1],Geom.point)
convPlot=plot(x=1:(convIt-1),y=avgN[1:convIt-1,1],Geom.point)
