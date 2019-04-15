# test of ALADIN with simple example
#Micah Botkin-Levy
#5/12/18

using JuMP
using Gurobi
using Gadfly

#minimize sum_i=1^N x_i^2*q_ii
#s.t. sum_i=1^N a_i*x_i=b


#set up Problem
N=10
Q=ones(N,1) #change to different weights??
A=round.(2.*rand(N) -1,4)
b=100*round.(rand(),4)

#central solution
m=Model(solver = GurobiSolver(OutputFlag=0))
@variable(m,xC[1:N])
@objective(m,Min,sum(0.5*xC[i,1]^2*Q[i,1] for i=1:N))
@constraint(m,lambdaP,sum(A[i,1]*xC[i,1] for i=1:N)==b)
status = solve(m)
if status!=:Optimal
    return
end
fStar=getobjectivevalue(m)
xStar=getvalue(xC)
lamStar=-getdual(lambdaP)

#ALADIN
#initialize
maxIt=50
epsilon=1e-10
rhoALAD=1

X0=ones(N,1)
#X0=xStar
lambda0=1
#lambda0=lamStar

convIt=maxIt
X=zeros(N,maxIt) #(rows are variable, columns are iteration)
G=zeros(N,maxIt)
Y=zeros(N,maxIt)
X[:,1]=X0
Lambda=zeros(maxIt,1)
Lambda[1,1]=lambda0
H=Q+round.(rand(N)/10,4)
avgN=zeros(maxIt,1)
fGap=zeros(maxIt,1)
itGap=zeros(maxIt,1)
constGap=zeros(maxIt,1)
convGap=zeros(maxIt,1)

for p=1:maxIt-1

    #solve NLP
    for i=1:N
        m=Model(solver = GurobiSolver(OutputFlag=0))
        @variable(m,y)
        @objective(m,Min,0.5*y^2*Q[i,1]+Lambda[p,1]*A[i,1]*y+rhoALAD/2*(y-X[i,p])^2*H[i,1])
        TT = STDOUT # save original STDOUT stream
		redirect_stdout()
        status = solve(m)
		redirect_stdout(TT)
        if status!=:Optimal
            return
        else
            Y[i,p+1]=getvalue(y)
            #G[i,p+1]=H[i,1]*(X[i,p]-Y[i,p+1])-A[i,1]*Lambda[p,1]
            G[i,p+1]=Q[i,1]*Y[i,p+1]
        end
    end

    #update H here
    #check for convergence
    objFun(xC)=sum(0.5*xC[i,1]^2*Q[i,1] for i=1:N)
    fGap[p,1]= abs(objFun(Y[:,p+1])-fStar)
    constGap[p,1]=norm(A'*Y[:,p+1]-b)
    itGap[p,1] = norm(Lambda[p,1]-Lambda[max(p-1,1),1],2)
    convGap[p,1] = norm(Lambda[p,1]-lamStar,2)
    avgN[p,1]=mean(norm.(X[:,p]-Y[:,p+1]))
    if  all(norm.(X[:,p]-Y[:,p+1]).<=epsilon)
        @printf "Converged after %g iterations\n" p
        convIt=p+1
        break
    else
        @printf "avgN       %e after %g iterations\n" avgN[p,1] p
        @printf "lastGap    %e after %g iterations\n" itGap[p,1] p
        @printf "convLamGap %e after %g iterations\n" convGap[p,1] p
        @printf "constGap   %e after %g iterations\n" constGap[p,1] p
        @printf("fGap       %e after %g iterations\n\n",fGap[p,1],p)
    end

    #coupled QP
    m=Model(solver = GurobiSolver(OutputFlag=0))
    @variable(m,dy[1:N])
    @objective(m,Min,sum(0.5*dy[i,1]^2*H[i,1]+G[i,p+1]*dy[i,1] for i=1:N))
    @constraint(m,lambdaP,sum(A[i,1]*(Y[i,p+1]+dy[i,1]) for i=1:N)==b)
    TT = STDOUT # save original STDOUT stream
    redirect_stdout()
    status = solve(m)
    redirect_stdout(TT)
    if status!=:Optimal
        return
    end

    #update
    alpha=1
    #X[:,p+1]=Y[:,p+1]+getvalue(dy)
    X[:,p+1]=X[:,p]+alpha*(Y[:,p+1]-X[:,p])+alpha*getvalue(dy)
    #Lambda[p+1,1]=-getdual(lambdaP)
    Lambda[p+1,1]=Lambda[p,1]+alpha*(-getdual(lambdaP)-Lambda[p,1])
end


compPlot=plot(layer(x=1:N,y=getvalue(xC),Geom.point),
              layer(x=1:N,y=Y[:,convIt],Geom.point))

convGapPlot=plot(x=1:(convIt-1),y=convGap[1:convIt-1,1],Geom.point,Scale.y_log10)
itGapPlot=plot(x=1:(convIt-1),y=itGap[1:convIt-1,1],Geom.point,Scale.y_log10)
fGapPlot=plot(x=1:(convIt-1),y=fGap[1:convIt-1,1],Geom.point,Scale.y_log10)
constGapPlot=plot(x=1:(convIt-1),y=constGap[1:convIt-1,1],Geom.point,Scale.y_log10)
convPlot=plot(x=1:(convIt-1),y=avgN[1:convIt-1,1],Geom.point,Scale.y_log10)

xPlot=plot(X[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel("X"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=N),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))

yPlot=plot(Y[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel("Y"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=N),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
