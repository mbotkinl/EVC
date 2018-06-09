# test of ALADIN with simple example
#Micah Botkin-Levy
#5/12/18

using JuMP
using Gurobi
using Ipopt
using Gadfly

sn0=s0
xt0=T0
target=0
imax=imax*100

etaP=.2*ones(N,1)
gammaP=.1
Tmax=80000

#set up Problem
m = Model(solver = IpoptSolver())

#u w and z are one index ahead of x. i.e the x[k+1]=x[k]+eta*u[k+1]
@variable(m,u[1:N])
@variable(m,sn[1:N])
@variable(m,xt)
@variable(m,z[1:S])
objFun(sn,u)=sum((sn[n,1]-1)^2*Qsi[n,1]     for n=1:N) +
			 sum((u[n,1])^2*Ri[n,1]           for n=1:N)
@objective(m,Min, objFun(sn,u))
@constraint(m,stateCon1[n=1:N],sn[n,1]/etaP[n,1]==sn0[n,1]/etaP[n,1]+u[n,1])
#@constraint(m,stateCon1,sn[1:N,1].==sn0[1:N,1]+etaP[:,1].*u[1:N,1])
@constraint(m,tempCon1,xt/gammaP==tauP*xt0/gammaP+deltaI*sum((2*m+1)*z[m+1,1] for m=0:S-1)+rhoP*w[2,1]/gammaP)
#@constraint(m,tempCon1,xt==tauP*xt0+gammaP*deltaI*sum((2*m+1)*z[m+1,1] for m=0:S-1)+rhoP*w[2,1])
@constraint(m,currCon,sum(u[n] for n=1:N)-sum(z[s] for s=1:S)==-w[1])
@constraint(m,sn.<=1)
@constraint(m,sn.>=target)
@constraint(m,upperTCon,xt<=Tmax)
@constraint(m,xt>=0)
@constraint(m,upperCCon,u.<=imax)
@constraint(m,u.>=imin)
@constraint(m,z.>=0)
@constraint(m,z.<=deltaI)

println("solving....")
statusM = solve(m)
if status!=:Optimal
    return
end

xtStar=getvalue(xt)
snStar=getvalue(sn)
uStar=getvalue(u)
zStar=getvalue(z)
fStar=getobjectivevalue(m)
lamTempStar=-getdual(upperTCon)
lamCurrStar=-getdual(currCon)



#ALADIN
#initialize
maxIt=50
convIt=maxIt
epsilon=1e-6
rhoALAD=1
muALAD=10^8
sigmaU=10*ones(N,1)
sigmaS=100*ones(N,1)
sigmaZ=10
sigmaT=1
tolU=1e-3
tolS=1e-6
tolT=1e-1
tolZ=1e-3

#guesses
#vs0=max.(snStar + rand(Truncated(Normal(0), -0.02, 0.02), N),0)
#vu0=max.(uStar + rand(Truncated(Normal(0), -0.1, 0.1), N),0)
#vz0=max.(zStar-2*rand(Truncated(Normal(0), 0, 1), S),0)
#vt0=max.(xtStar-10*rand(Truncated(Normal(0), 0, 1),1),0)
#lambda0=max.(lamCurrStar-rand(Truncated(Normal(0), 0, 1), 1),0)
#lambda0=zeros(horzLen+1,1)

# vs0=snStar
# vu0=uStar
# vz0=zStar
# vt0=xtStar
# lambda0=lamCurrStar

vs0=rand(Truncated(Normal(0), 0, 1), N)
vu0=imax[1,1]*0.8*rand(Truncated(Normal(0), 0, 1), N)
vz0=ItotalMax*rand(Truncated(Normal(0), 0, 1), S)
vt0=Tmax*rand(Truncated(Normal(0), 0, 1), 1)
lambda0=5*rand(Truncated(Normal(0), 0, 1), 1)

#save matrices
Un=SharedArray{Float64}(N,maxIt) #row are time (N states for k=1, them N states for k=2),  columns are iteration
Sn=SharedArray{Float64}(N,maxIt)  #row are time,  columns are iteration
T=zeros(1,maxIt)  #row are time,  columns are iteration
Z=zeros(S,maxIt)  #row are time,  columns are iteration
Lambda=zeros(1,maxIt) #(rows are time, columns are iteration)
Lambda[1,1]=lambda0[1]

#auxillary variables
Vu=zeros(N,maxIt) #row are time,  columns are iteration
Vu[:,1]=vu0 #initial guess goes in column 1
Vs=zeros(N,maxIt) #row are time,  columns are iteration
Vs[:,1]=vs0 #initial guess goes in column 1
Vz=zeros(S,maxIt) #row are time,  columns are iteration
Vz[:,1]=vz0
Vt=zeros(1,maxIt) #row are time,  columns are iteration
Vt[1,1]=vt0[1]

#Gradian Vectors
Gu=SharedArray{Float64}(N,maxIt) #row are time (N states for k=1, them N states for k=2),  columns are iteration
Gs=SharedArray{Float64}(N,maxIt) #row are time (N states for k=1, them N states for k=2),  columns are iteration
Gz=zeros(S,maxIt) #row are time (N states for k=1, them N states for k=2),  columns are iteration
Gt=zeros(1,maxIt) #row are time

#Jacobian C Vectors
Cs=SharedArray{Float64}(N,maxIt)  #row are time,  columns are iteration
Cu=SharedArray{Float64}(N,maxIt)  #row are time,  columns are iteration
Cz=zeros(S,maxIt)  #row are time,  columns are iteration
Ct=zeros(1,maxIt)  #row are time,  columns are iteration

#conv metrics
convGap=zeros(maxIt,1)
constGap=zeros(maxIt,1)
itGap=zeros(maxIt,1)
fGap=zeros(maxIt,1)


for p=1:maxIt-1

    #solve NLP
    for i=1:N
        #nlp=Model(solver = GurobiSolver(OutputFlag=0))
        nlp=Model(solver = IpoptSolver())
        @variable(nlp,un)
        @variable(nlp,sn)
        @objective(nlp,Min,un^2*Ri[i,1]+(sn-1)^2*Qsi[i,1]+Lambda[1,p]*un+
                            rhoALAD/2*(un-Vu[i,p])^2*sigmaU[i,1]+rhoALAD/2*(sn-Vs[i,p])^2*sigmaS[i,1])
        @constraint(nlp,sn==sn0[i,1]+etaP[i,1]*un)
        @constraint(nlp,kapMax,sn<=1)
        @constraint(nlp,kapMin,sn>=0)
        @constraint(nlp,kapMaxU,un<=imax[i,1])
        @constraint(nlp,kapMinU,un>=0)
        TT = STDOUT # save original STDOUT stream
		redirect_stdout()
        status = solve(nlp)
		redirect_stdout(TT)
        if status!=:Optimal
            return
        else
            Sn[i,p+1]=getvalue(sn)
            Un[i,p+1]=getvalue(un)


            cValMax=abs.(Sn[i,p+1]-1).<tolS
            cValMin=abs.(Sn[i,p+1]-0).<tolS
            Cs[i,p+1]=1cValMax-1cValMin

            cValMax=abs.(Un[i,p+1]-imax[i,1]).<tolU
            cValMin=abs.(Un[i,p+1]-0).<tolU
            Cu[i,p+1]=1cValMax-1cValMin

            #Gy[i,p+1]=H[i,1]*(X[i,p]-Y[i,p+1])-A[i,1]*Lambda[p,1]
            Gu[i,p+1]=2*Ri[i,1]*Un[i,p+1]

            Gs[i,p+1]=2*Qsi[i,1]*Sn[i,p+1]-2*Qsi[i,1]
        end
    end

    #N+1 decupled
    nlp=Model(solver = IpoptSolver())
    @variable(nlp,z[1:S])
    @variable(nlp,t)
    @objective(nlp,Min,-Lambda[1,p]*sum(z[i] for i=1:S)+
                        rhoALAD/2*(t-Vt[1,p])^2*sigmaT+rhoALAD/2*sigmaZ*sum((z[s]-Vz[s,p])^2 for s=1:S))#rhoALAD/2*norm(z-Vz[:,p],2)^2*sigmaZ)
    #@constraint(nlp,t==tauP*T0+gammaP*deltaI*sum((2*s-1)*z[s] for s=1:S)+rhoP*w[2,1])
	@constraint(nlp,t/gammaP==tauP*T0/gammaP+deltaI*sum((2*s-1)*z[s] for s=1:S)+rhoP*w[2,1]/gammaP)
    @constraint(nlp,kapMax,t<=Tmax)
    @constraint(nlp,kapMin,t>=0)
    @constraint(nlp,kapMaxU,z.<=deltaI)
    @constraint(nlp,kapMinU,z.>=0)
    TT = STDOUT # save original STDOUT stream
    redirect_stdout()
    status = solve(nlp)
    redirect_stdout(TT)
    if status!=:Optimal
        return
    else
        Z[:,p+1]=getvalue(z)
        T[1,p+1]=getvalue(t)


        cValMax=abs.(Z[:,p+1]-deltaI).<tolZ
        cValMin=abs.(Z[:,p+1]-0).<tolZ
        Cz[:,p+1]=1cValMax-1cValMin

        cValMax=abs(T[1,p+1]-Tmax)<tolT
        cValMin=abs(T[1,p+1]-0)<tolT
        Ct[1,p+1]=1cValMax-1cValMin

        Gt[1,p+1]=0
        Gz[:,p+1]=0
        #Gz[:,p+1]=Lambda[1,p]
    end


    #update H here
    #check for convergence
    objFun(sn,u)=sum((sn[n,1]-1)^2*Qsi[n,1]     for n=1:N) +
    			 sum((u[n,1])^2*Ri[n,1]           for n=1:N)
    fGap[p,1]= abs(objFun(Sn[:,p+1],Un[:,p+1])-fStar)
    constGap[p,1]=norm(sum(Un[i,p+1] for i=1:N)-sum(Z[s,p+1] for s=1:S)+w[1,1])
    itGap[p,1] = norm(Lambda[1,p]-Lambda[1,max(p-1,1)],2)
    convGap[p,1] = norm(Lambda[1,p]-lamCurrStar,2)
    testC=norm(vcat(Vs[:,p]-Sn[:,p+1],Vu[:,p]-Un[:,p+1]))
    if  all(testC<=epsilon)
        @printf "Converged after %g iterations\n" p
        convIt=p+1
        break
    else
		@printf "testC      %e after %g iterations\n" testC p
        @printf "lastGap    %e after %g iterations\n" itGap[p,1] p
        @printf "convLamGap %e after %g iterations\n" convGap[p,1] p
        @printf "constGap   %e after %g iterations\n" constGap[p,1] p
        @printf("fGap       %e after %g iterations\n\n",fGap[p,1],p)
    end

    #coupled QP
    #mC=Model(solver = GurobiSolver(OutputFlag=0))
    mC=Model(solver = IpoptSolver())
    @variable(mC,dUn[1:N])
    @variable(mC,dSn[1:N])
    @variable(mC,dT)
    @variable(mC,dZ[1:S])
    #@variable(mC,relaxS)
    objExp=sum(0.5*dUn[i,1]^2*Ri[i,1]+Gu[i,p+1]*dUn[i,1] for i=1:N)+sum(0.5*dSn[i,1]^2*Qsi[i,1]+Gs[i,p+1]*dSn[i,1] for i=1:N)
    #objExp=objExp+Lambda[p,1]*relaxS+muALAD/2*relaxS^2
    @objective(mC,Min,objExp)
    @constraint(mC,lambdaP,sum((Un[i,p+1]+dUn[i,1]) for i=1:N)+sum(-(Z[s,p+1]+dZ[s,1]) for s=1:S)==-w[1])#+relaxS)
    @constraint(mC,[i=1:N],dSn[i]/etaP[i,1]-dUn[i]==0)
	#@constraint(mC,dSn-etaP[:,1].*dUn.==0)
    @constraint(mC,dT/gammaP-deltaI*sum((2*s-1)*dZ[s] for s=1:S)==0)
	#@constraint(mC,dT-gammaP*deltaI*sum((2*s-1)*dZ[s] for s=1:S)==0)
    @constraint(mC,Cu[:,p+1].*dUn.<=0)
    @constraint(mC,Cs[:,p+1].*dSn.<=0)
    @constraint(mC,Cz[:,p+1].*dZ.<=0)
    @constraint(mC,Ct[1,p+1]*dT<=0)

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
    Vs[:,p+1]=Vs[:,p]+alpha*(Sn[:,p+1]-Vs[:,p])+alpha*getvalue(dSn)
    Vz[:,p+1]=Vz[:,p]+alpha*(Z[:,p+1]-Vz[:,p])+alpha*getvalue(dZ)
    Vu[:,p+1]=Vu[:,p]+alpha*(Un[:,p+1]-Vu[:,p])+alpha*getvalue(dUn)
    Vt[1,p+1]=Vt[1,p]+alpha*(T[1,p+1]-Vt[1,p])+alpha*getvalue(dT)

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
			Coord.Cartesian(xmin=0,xmax=N),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))

yPlot=plot(Y[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel("Y"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=N),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
