# test of ALADIN with simple example
#Micah Botkin-Levy
#5/12/18

using JuMP
using Gurobi
using Ipopt
using Gadfly

#minimize sum_i=1^N x_i^2*q_ii+sum_i=1^N (z_i-1)^2*r_ii
#s.t. sum_i=1^N a_i*x_i=b
# xmin<=x_i<=xmax
# zmin<=z_i<=zmax
# x-z=d


#set up Problem

#function testALAD()
N=10
Q=ones(N,1) #change to different weights??
R=10*ones(N,1) #change to different weights??
#A=round.(2.*rand(N) -1,4)
A=[-0.3567;0.27;-0.1726;0.4832;-0.3183;-0.3668;0.954;0.0975;0.9486;0.3887]
#b=round.(20*rand()+10,4)
b=24.7813
d=-10
xmin=-8
xmax=8
zmin=2
zmax=16
#central solution
m=Model(solver = GurobiSolver(OutputFlag=0))
@variable(m,xC[1:N])
@variable(m,zC[1:N])
@objective(m,Min,sum(xC[i,1]^2*Q[i,1] for i=1:N)+sum((zC[i,1]-1)^2*R[i,1] for i=1:N))
@constraint(m,lambdaP,sum(A[i,1]*xC[i,1] for i=1:N)==b)
@constraint(m,xC-zC.==d)
@constraint(m,kapMax,xC.<=xmax)
@constraint(m,kapMin,xC.>=xmin)
@constraint(m,kapMaxZ,zC.<=zmax)
@constraint(m,kapMinZ,zC.>=zmin)
status = solve(m)


if status!=:Optimal
    return
end

fStar=getobjectivevalue(m)
xStar=getvalue(xC)
zStar=getvalue(zC)
lamStar=-getdual(lambdaP)


#ALADIN
#initialize
maxIt=50
epsilon=1e-8
tolC=1e-3
rhoALAD=1
muALAD=10^8

#x0=ones(N,1)
#z0=ones(N,1)
x0=xStar
z0=zStar
#lambda0=1
lambda0=lamStar

convIt=maxIt
X=zeros(N,maxIt) #(rows are variable, columns are iteration)
V=zeros(N,maxIt) #(rows are variable, columns are iteration)
Z=zeros(N,maxIt)
Y=zeros(N,maxIt)
Cy=zeros(N,maxIt)
Cv=zeros(N,maxIt)
Gy=zeros(N,maxIt)
Gv=zeros(N,maxIt)


X[:,1]=x0
Z[:,1]=z0
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
        nlp=Model(solver = GurobiSolver(OutputFlag=0))
        @variable(nlp,y)
        @variable(nlp,v)
        @objective(nlp,Min,y^2*Q[i,1]+(v-1)^2*R[i,1]+Lambda[p,1]*A[i,1]*y+
                            rhoALAD/2*(y-X[i,p])^2*Q[i,1]+rhoALAD/2*(v-Z[i,p])^2*R[i,1])
        @constraint(nlp,y-v==d)
        @constraint(nlp,kapMax,y<=xmax)
        @constraint(nlp,kapMin,y>=xmin)
        @constraint(nlp,kapMaxZ,v<=zmax)
        @constraint(nlp,kapMinZ,v>=zmin)
        TT = STDOUT # save original STDOUT stream
		redirect_stdout()
        status = solve(nlp)
		redirect_stdout(TT)
        if status!=:Optimal
            return
        else
            Y[i,p+1]=getvalue(y)
            V[i,p+1]=getvalue(v)


            cValMax=abs.(Y[i,p+1]-xmax).<tolC
            cValMin=abs.(Y[i,p+1]-xmin).<tolC
            Cy[i,p+1]=1cValMax-1cValMin

            cValMax=abs.(V[i,p+1]-zmax).<tolC
            cValMin=abs.(V[i,p+1]-zmin).<tolC
            Cv[i,p+1]=1cValMax-1cValMin

            #Gy[i,p+1]=H[i,1]*(X[i,p]-Y[i,p+1])-A[i,1]*Lambda[p,1]
            Gy[i,p+1]=2*Q[i,1]*Y[i,p+1]

            Gv[i,p+1]=2*R[i,1]*V[i,p+1]-2*R[i,1]
        end
    end

    #update H here
    #check for convergence
    objFun(xC,zC)=sum(xC[i,1]^2*Q[i,1] for i=1:N)+sum((zC[i,1]-1)^2*R[i,1] for i=1:N)
    fGap[p,1]= abs(objFun(Y[:,p+1],V[:,p+1])-fStar)
    constGap[p,1]=norm(A'*Y[:,p+1]-b)
    itGap[p,1] = norm(Lambda[p,1]-Lambda[max(p-1,1),1],2)
    convGap[p,1] = norm(Lambda[p,1]-lamStar,2)
    avgN[p,1]=mean(norm.(X[:,p]-Y[:,p+1])) #fix this
    if  all(norm.(X[:,p]-Y[:,p+1]).<=epsilon) && constGap[p,1]<=epsilon
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
    #mC=Model(solver = GurobiSolver(OutputFlag=0))
    mC=Model(solver = IpoptSolver())
    @variable(mC,dy[1:N])
    @variable(mC,dv[1:N])
    #@variable(mC,relaxS)
    objExp=sum(0.5*dy[i,1]^2*Q[i,1]+Gy[i,p+1]*dy[i,1] for i=1:N)+sum(0.5*dv[i,1]^2*R[i,1]+Gv[i,p+1]*dv[i,1] for i=1:N)
    #objExp=objExp+Lambda[p,1]*relaxS+muALAD/2*relaxS^2
    @objective(mC,Min,objExp)
    @constraint(mC,lambdaP,sum(A[i,1]*(Y[i,p+1]+dy[i,1]) for i=1:N)==b)#+relaxS)
    @constraint(mC,dy-dv.==0)
    @constraint(mC,Cy[:,p+1].*dy.<=0)
    @constraint(mC,Cv[:,p+1].*dv.<=0)
    TT = STDOUT # save original STDOUT stream
    redirect_stdout()
    status = solve(mC)
    redirect_stdout(TT)
    if status!=:Optimal
        return
    end

    #update
    alpha=1
    #X[:,p+1]=Y[:,p+1]+getvalue(dy)
    X[:,p+1]=X[:,p]+alpha*(Y[:,p+1]-X[:,p])+alpha*getvalue(dy)
    Z[:,p+1]=Z[:,p]+alpha*(V[:,p+1]-Z[:,p])+alpha*getvalue(dv)

    #Lambda[p+1,1]=-getdual(lambdaP)
    Lambda[p+1,1]=Lambda[p,1]+alpha*(-getdual(lambdaP)-Lambda[p,1])
end
#end

#testALAD()

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
