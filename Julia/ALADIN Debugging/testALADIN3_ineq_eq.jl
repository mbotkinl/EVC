# test of ALADIN with simple example
#Micah Botkin-Levy
#5/12/18

using JuMP
using Gurobi
using Ipopt
using Gadfly

#minimize sum_i=1^N+1 x_i^2*q_ii+sum_i=1^N+1 (z_i-1)^2*r_ii (q_ii=r_ii=0)
#s.t. sum_i=1^N+1 a_i*x_i=b (a_N+1=-1)
# xmin<=x_i<=xmax
# zmin<=z_i<=zmax
# zi=z0+eta*xi (eta_N+1=gamma)


#set up Problem

#function testALAD()
N=10
Q=[ones(N,1);0] #change to different weights??
R=[10*ones(N,1);0] #change to different weights??
A=[round.(rand(N),4);-1]
#A=[-0.3567;0.27;-0.1726;0.4832;-0.3183;-0.3668;0.954;0.0975;0.9486;0.3887]
#b=-round.(20*rand()+10,4)
b=0
#b=24.7813
eta=[.2*ones(N,1);5]
z0=[0.2*rand(N,1)+.5 ;90]
xmin=zeros(N+1,1)
xmax=[8*ones(N,1);8*N]
zmin=zeros(N+1,1)
zmax=[ones(N,1);95]


#central solution
m=Model(solver = IpoptSolver())
#m=Model(solver = GurobiSolver())
@variable(m,xC[1:N+1])
@variable(m,zC[1:N+1])
@objective(m,Min,sum(xC[i,1]^2*Q[i,1] for i=1:N+1)+sum((zC[i,1]-1)^2*R[i,1] for i=1:N+1))
@constraint(m,lambdaP,sum(A[i]*xC[i] for i=1:N+1)==b)
@constraint(m,socConst[i=1:N+1],zC[i]==z0[i]+eta[i]*xC[i])
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
epsilon=1e-4
tolC=1e-1
rhoALAD=1
muALAD=10^8
sigmaX=[ones(N,1);1]
sigmaZ=[10*ones(N,1);1]
tolX=[.01*ones(N,1);.1]
tolZ=[.01*ones(N,1);.1]

vx0=ones(N+1,1)
vz0=ones(N+1,1)
lambda0=1
# vx0=xStar
# vz0=zStar
# lambda0=lamStar

convIt=maxIt
X=zeros(N+1,maxIt) #(rows are variable, columns are iteration)
Z=zeros(N+1,maxIt)
Vx=zeros(N+1,maxIt) #(rows are variable, columns are iteration)
Vz=zeros(N+1,maxIt) #(rows are variable, columns are iteration)
Cz=zeros(N+1,maxIt)
Cx=zeros(N+1,maxIt)
Gx=zeros(N+1,maxIt)
Gz=zeros(N+1,maxIt)


Vx[:,1]=vx0
Vz[:,1]=vz0
Lambda=zeros(1,maxIt)
Lambda[1,1]=lambda0
avgN=zeros(maxIt,1)
fGap=zeros(maxIt,1)
itGap=zeros(maxIt,1)
constGap=zeros(maxIt,1)
convGap=zeros(maxIt,1)

for p=1:maxIt-1

    #solve NLP
    for i=1:N+1
        nlp=Model(solver = GurobiSolver(OutputFlag=0))
        @variable(nlp,x)
        @variable(nlp,z)
        @objective(nlp,Min,x^2*Q[i,1]+(z-1)^2*R[i,1]+Lambda[1,p]*A[i,1]*x)#+
                            #rhoALAD/2*(x-Vx[i,p])^2*sigmaX[i,1]+rhoALAD/2*(z-Vz[i,p])^2*sigmaZ[i,1])
        @constraint(nlp,z==z0[i,1]+eta[i,1]*x)
        @constraint(nlp,kapMax,x<=xmax[i,1])
        @constraint(nlp,kapMin,x>=xmin[i,1])
        @constraint(nlp,kapMaxZ,z<=zmax[i,1])
        @constraint(nlp,kapMinZ,z>=zmin[i,1])
        TT = STDOUT # save original STDOUT stream
		redirect_stdout()
        status = solve(nlp)
		redirect_stdout(TT)
        if status!=:Optimal
            return
        else
            X[i,p+1]=getvalue(x)
            Z[i,p+1]=getvalue(z)


            cValMax=abs.(X[i,p+1]-xmax[i,1]).<tolX[i,1]
            cValMin=abs.(X[i,p+1]-xmin[i,1]).<tolX[i,1]
            Cx[i,p+1]=1cValMax-1cValMin

            cValMax=abs.(Z[i,p+1]-zmax[i,1]).<tolZ[i,1]
            cValMin=abs.(Z[i,p+1]-zmin[i,1]).<tolZ[i,1]
            Cz[i,p+1]=1cValMax-1cValMin

            #Gy[i,p+1]=H[i,1]*(X[i,p]-Y[i,p+1])-A[i,1]*Lambda[p,1]
            Gx[i,p+1]=2*Q[i,1]*X[i,p+1]

            Gz[i,p+1]=2*R[i,1]*Z[i,p+1]-2*R[i,1]
        end
    end

    #update H here
    #check for convergence
    objFun(xC,zC)=sum(xC[i,1]^2*Q[i,1] for i=1:N+1)+sum((zC[i,1]-1)^2*R[i,1] for i=1:N+1)
    fGap[p,1]= abs(objFun(X[:,p+1],Z[:,p+1])-fStar)
    constGap[p,1]=norm(A'*X[:,p+1]-b)
    itGap[p,1] = norm(Lambda[1,p]-Lambda[1,max(p-1,1)],2)
    convGap[p,1] = norm(Lambda[1,p]-lamStar,2)
    testC=norm(vcat(Vx[:,p]-X[:,p+1],Vz[:,p]-Z[:,p+1])) #add sigma???
    if  testC<=epsilon && constGap[p,1]<=epsilon
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
    @variable(mC,dx[1:N+1])
    @variable(mC,dz[1:N+1])
    #@variable(mC,relaxS)
    objExp=sum(0.5*dx[i,1]^2*Q[i,1]+Gx[i,p+1]*dx[i,1]+0.5*dz[i,1]^2*R[i,1]+Gz[i,p+1]*dz[i,1] for i=1:N+1)
    #objExp=objExp+Lambda[1,p]*relaxS+muALAD/2*relaxS^2
    @objective(mC,Min,objExp)
    @constraint(mC,lambdaP,sum(A[i,1]*(X[i,p+1]+dx[i,1]) for i=1:N+1)==b)#+relaxS)
    @constraint(mC,eqC[i=1:N+1],dz[i]-eta[i,1]*dx[i]==0)
    @constraint(mC,Cz[:,p+1].*dz.<=0)
    @constraint(mC,Cx[:,p+1].*dx.<=0)
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
    Vx[:,p+1]=Vx[:,p]+alpha*(X[:,p+1]-Vx[:,p])+alpha*getvalue(dx)
    Vz[:,p+1]=Vz[:,p]+alpha*(Z[:,p+1]-Vz[:,p])+alpha*getvalue(dz)

    #Lambda[p+1,1]=-getdual(lambdaP)
    Lambda[1,p+1]=Lambda[1,p]+alpha*(-getdual(lambdaP)-Lambda[1,p])
end
#end

#testALAD()



convGapPlot=plot(x=1:(convIt-1),y=convGap[1:convIt-1,1],Geom.point,Scale.y_log10)
itGapPlot=plot(x=1:(convIt-1),y=itGap[1:convIt-1,1],Geom.point,Scale.y_log10)
fGapPlot=plot(x=1:(convIt-1),y=fGap[1:convIt-1,1],Geom.point,Scale.y_log10)
constGapPlot=plot(x=1:(convIt-1),y=constGap[1:convIt-1,1],Geom.point,Scale.y_log10)

xPlot=plot(X[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel("X"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=N+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))

zPlot=plot(Z[1:N,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel("Y"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=N),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
