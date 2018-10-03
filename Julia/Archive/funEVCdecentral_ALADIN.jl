
#function for setting up ALADIN variables

#functions for initializing opt models
function evModel_init(;evInd=1,p=1)#horzLen,ηP,σU,σS,target,imax,imin,Qsi,Ri)
    evM = Model(solver = GurobiSolver(NumericFocus=3))
    @variable(evM,sn[1:(horzLen+1)])
    @variable(evM,u[1:(horzLen+1)])
    @NLparameter(evM,evVu[1:(horzLen+1)]==1)
    @NLparameter(evM,evVs[1:(horzLen+1)]==1)
    @NLparameter(evM,lambda[1:(horzLen+1)]==1)
    @NLparameter(evM,target[1:(horzLen+1)]==1)
    @NLobjective(evM,Min,sum((sn[k,1]-1)^2*Qsi[evInd,1]+(u[k,1])^2*Ri[evInd,1]+
                            lambda[k,1]*(u[k,1])+
                            ρALADp[1,p]/2*(u[k,1]-evVu[k,1])*σU[evInd,1]*(u[k,1]-evVu[k,1])+
                            ρALADp[1,p]/2*(sn[k,1]-evVs[k,1])*σS[evInd,1]*(sn[k,1]-evVs[k,1]) for k=1:horzLen+1))
    @constraint(evM,sn[1,1]==sn0[evInd,1]+ηP[evInd,1]*u[1,1])
    @constraint(evM,[k=1:horzLen],sn[k+1,1]==sn[k,1]+ηP[evInd,1]*u[k+1,1])
    @constraint(evM,socKappaMax,sn.<=1)
    @constraint(evM,socKappaMin,sn.>=target)
    @constraint(evM,curKappaMax,u.<=imax[evInd,1])
    @constraint(evM,curKappaMin,u.>=imin[evInd,1])
    TT = STDOUT # save original STDOUT stream
    redirect_stdout()
    status = solve(evM)
    redirect_stdout(TT)
    return evM
end

#XFRM nlp
function xfrmModel_init()
    tM = Model(solver = IpoptSolver())

    #optimization variables
    @variable(tM,z[1:(S)*(horzLen+1)])
    @variable(tM,xt[1:(horzLen+1)])

    #input variables
    @variable(tM,lambda[1:(horzLen+1)]==1)
    @variable(tM,tempVz[1:S*(horzLen+1)]==1)
    @variable(tM,tempVt[1:(horzLen+1)]==1)
    @variable(tM,ρALADp==1)


    @NLobjective(tM,Min, sum(-lambda[k,1]*sum(z[(k-1)*(S)+s] for s=1:S)+
                        ρALADp/2*σZ*(sum(z[(k-1)*(S)+s] for s=1:S)-sum(Vz[(k-1)*(S)+s] for s=1:S))^2+
                        ρALADp/2*σT*(xt[k]-Vt[k,1])^2  for k=1:(horzLen+1)))
    @constraint(tM,tempCon1,xt[1]-τP*xt0-γP*deltaI*sum((2*s-1)*z[s] for s=1:S)-ρP*w[stepI*2,1]==0) #need to add w to variables for MPC???
    @constraint(tM,tempCon2[k=1:horzLen],xt[k+1]-τP*xt[k]-γP*deltaI*sum((2*s-1)*z[(k)*(S)+s] for s=1:S)-ρP*w[stepI*2+k*2,1]==0)
    if noTlimit==0
        @constraint(tM,upperTCon,xt.<=Tmax)
    end
    @constraint(tM,lowerTCon,xt.>=0)
    @constraint(tM,pwlKappaMin,z.>=0)
    @constraint(tM,pwlKappaMax,z.<=deltaI)
    TT = STDOUT # save original STDOUT stream
    redirect_stdout()
    statusTM = solve(tM)
    redirect_stdout(TT)

    return tM,lambda,tempVz,tempVt,ρALADp
end

#Central QP
function cModel_init()
end


#functions used in loop
function nlpEVfun(p,evInd,Cu,Cs,Sn,Un,Gs,Gu)

    ind=[evInd]
    for k=1:horzLen
        append!(ind,k*N+evInd)
    end
    target=zeros((horzLen+1),1)
    target[(Kn[evInd,1]-(stepI-1)):1:length(target),1]=Snmin[evInd,1]

    #assume new values
    setvalue(evVu,Vu[ind,p])
    setvalue(evVs,Vs[ind,p])
    setvalue(lambdaEV, Lam[:,p])
    setvalue(targetEV, target)
    setvalue(ρALADev,ρALADp[1,p])
    setvalue(Qn,Qsi[evInd,1])
    setvalue(Rn,Ri[evInd,1])
    setvalue(σUn,σU[evInd,1])
    setvalue(σSn,σS[evInd,1])
    setvalue(ηPn,ηP[evInd,1])
    setvalue(sn0n,sn0[evInd,1])
    setvalue(imaxn,imax[evInd,1])
    setvalue(iminn,imin[evInd,1])

    TT = STDOUT # save original STDOUT stream
    redirect_stdout()
    statusEVM = solve(evM)
    redirect_stdout(TT)

    @assert statusEVM==:Optimal "solver issues with EV NLP"

    # kappaMax=-getdual(curKappaMax)
    # kappaMin=-getdual(curKappaMin)
    # socMax=-getdual(socKappaMax)
    # socMin=-getdual(socKappaMin)
    uVal=getvalue(u)
    snVal=getvalue(sn)

    cValMax=abs.(uVal-imax[evInd,1]).<tolU
    cValMin=abs.(uVal-imin[evInd,1]).<tolU
    # cVal=kappaMax
    # cVal[cVal.>0]=1
    # cVal=kappaMin
    # cVal[cVal.<0]=-1
    Cu[ind,p+1]=1cValMax-1cValMin


    cValMax=abs.(snVal-1).<tolS
    cValMin=abs.(snVal-target).<tolS
    # cVal=socMax
    # cVal[cVal.>0]=1
    # cVal=socMin
    # cVal[cVal.<0]=-1
    Cs[ind,p+1]=1cValMax-1cValMin

    Sn[ind,p+1]=snVal
    Un[ind,p+1]=uVal
    Gu[ind,p+1]=2*Ri[evInd,1]*uVal
    Gs[ind,p+1]=2*Qsi[evInd,1]*snVal-2*Qsi[evInd,1]
    #Gu[collect(evInd:N:length(Gu[:,p+1])),p+1]=σU[evInd,1]*(evVu-uVal)+lambda
    #Gs[collect(evInd:N:length(Gs[:,p+1])),p+1]=σN[evInd,1]*(evVs-snVal)-lambda
    return
end
function nlpXFRMfun(p,Ct,Cz,Xt,Z,Gz,Gt)

    #assume new values
    setvalue(tempVz,Vz[:,p])
    setvalue(tempVt,Vt[:,p])
    setvalue(lambdaXFRM,Lam[:,p])
    setvalue(ρALADXFRM,ρALADp[1,p])

    TT = STDOUT # save original STDOUT stream
    redirect_stdout()
    statusTM = solve(tM)
    redirect_stdout(TT)

    @assert statusTM==:Optimal "solver issues with XFRM NLP"

    #kappaMax=-getdual(pwlKappaMax)
    #kappaMin=-getdual(pwlKappaMin)
    #tMax=-getdual(upperTCon)
    #tMin=-getdual(lowerTCon)
    lambdaTemp=[-getdual(tempCon1);-getdual(tempCon2)]
    zVal=getvalue(z)
    xtVal=getvalue(xt)

    cValMax=abs.(zVal-deltaI).<tolZ
    cValMin=abs.(zVal-0).<tolZ
    # cVal=kappaMax
    # cVal=kappaMin
    # cVal[cVal.<0]=-1
    # cVal[cVal.>0]=1
    Cz[:,p+1]=1cValMax-1cValMin

    cValMax=abs.(xtVal-Tmax).<tolT
    cValMin=abs.(xtVal-0).<tolT
    # cVal=tMin
    # cVal[cVal.<0]=-1
    # cVal=tMax
    # cVal[cVal.>0]=1
    Ct[:,p+1]=1cValMax-1cValMin

    Xt[:,p+1]=xtVal
    Z[:,p+1]=zVal
    Gz[:,p+1]=0
    #Gz[:,p+1]=σZ*(Vz[:,p]-zVal)-repeat(-Lam[:,p],inner=S)
    Gt[:,p+1]=0
    return
end
function centralQPfun(p,lamQP,save_dUn,save_dZ,save_dSn,save_dXt)

    #assume new values
    setvalue(tempCt,Ct[:,p+1])
    setvalue(tempCu,Cu[:,p+1])
    setvalue(tempCs,Cs[:,p+1])
    setvalue(tempUn,Un[:,p+1])
    setvalue(tempGu,Gu[:,p+1])
    setvalue(tempGs,Gs[:,p+1])
    setvalue(tempCz,Cz[:,p+1])
    setvalue(tempZ,Z[:,p+1])

    TT = STDOUT # save original STDOUT stream
    redirect_stdout()
    statusCM = solve(cM)
    redirect_stdout(TT)


    @assert statusCM==:Optimal "solver issues with Central QP"

    save_dXt[:,p+1]=getvalue(dXt)
    save_dUn[:,p+1]=getvalue(dUn)
    save_dSn[:,p+1]=getvalue(dSn)
    save_dZ[:,p+1]=getvalue(dZ)
    lamQP[:,p+1]=-getdual(currCon)
    return
end

#wrapper for each iteration
function runALADit(p)

    #solve decoupled
    @sync @parallel for evInd=1:N
        nlpEVfun(p,evInd,Cu,Cs,Sn,Un,Gs,Gu)
    end

    #N+1 decoupled problem aka transformer current
    nlpXFRMfun(p,Ct,Cz,Xt,Z,Gz,Gt)

    for k=1:horzLen+1
        uSum[k,p+1]=sum(Un[(k-1)*N+n,p+1] for n=1:N)
        zSum[k,p+1]=sum(Z[(k-1)*(S)+s,p+1] for s=1:S)
        currConst[k,p+1]=uSum[k,p+1] + w[(k-1)*2+(stepI*2-1),1] - zSum[k,p+1]
    end

    #check for convergence
    constGap=norm(currConst[:,p+1],1)
    cc=norm(vcat((Vu[:,p]-Un[:,p+1]),(Vz[:,p]-Z[:,p+1])),1)
    #convCheck=ρALAD*norm(vcat(repmat(σU,horzLen+1,1).*(Vu[:,p]-Un[:,p+1]),σZ*(Vz[:,p]-Z[:,p+1])),1)
    objFun(sn,xt,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*Qsi[n,1]     for n=1:N) for k=1:horzLen+1) +
                    sum((xt[k,1]-1)^2*Qsi[N+1,1]                 for k=1:horzLen+1) +
                    sum(sum((u[(k-1)*N+n,1])^2*Ri[n,1]           for n=1:N) for k=1:horzLen+1)
    fGap= abs(objFun(Sn[:,p+1],Xt[:,p+1],Un[:,p+1])-fStar)
    snGap=norm((Sn[:,p+1]-snStar),2)
    unGap=norm((Un[:,p+1]-uStar),2)
    itGap = norm(Lam[:,p]-Lam[:,max(p-1,1)],2)
    convGap = norm(Lam[:,p]-lamCurrStar,2)
    fConvALAD[p,1]=fGap
    snConvALAD[p,1]=snGap
    unConvALAD[p,1]=unGap
    itConvALAD[p,1]=itGap
    constConvALAD[p,1]=constGap
    ConvALAD[p,1]=convGap
    convCheck[p,1]=cc
    if  constGap<=epsilon && convCheck<=epsilon
        @printf "Converged after %g iterations\n" p
        convIt=p+1
        #break
    else
        @printf "lastGap    %e after %g iterations\n" itGap p
        @printf "convLamGap %e after %g iterations\n" convGap p
        @printf "convCheck  %e after %g iterations\n" cc p
        @printf "constGap   %e after %g iterations\n" constGap p
        @printf "snGap      %e after %g iterations\n" snGap p
        @printf("fGap       %e after %g iterations\n\n",fGap,p)
    end

    centralQPfun(p,lamQP,save_dUn,save_dZ,save_dSn,save_dXt)

    #update step
    # Lam[:,p+1]=-getdual(currCon)
    α1=1
    α2=1
    α3=1
    #alpha3=alpha3/ceil(p/2)
    #Lam[:,p+1]=Lam[:,p]+alpha3*(-getdual(currCon)-Lam[:,p])
    Lam[:,p+1]=max.(Lam[:,p]+α3*(lamQP[:,p+1]-Lam[:,p]),0)

    Vu[:,p+1]=Vu[:,p]+α1*(Un[:,p+1]-Vu[:,p])+α2*save_dUn[:,p+1]
    Vz[:,p+1]=Vz[:,p]+α1*(Z[:,p+1]-Vz[:,p])+α2*save_dZ[:,p+1]
    Vs[:,p+1]=Vs[:,p]+α1*(Sn[:,p+1]-Vs[:,p])+α2*save_dSn[:,p+1]
    Vt[:,p+1]=Vt[:,p]+α1*(Xt[:,p+1]-Vt[:,p])+α2*save_dXt[:,p+1]

    ρALADp[1,p+1]=ρALADp[1,p]*ρRate #increase ρ every iteration
    deltaY[1,p+1]=norm(vcat(save_dUn[:,p+1],save_dZ[:,p+1],save_dSn[:,p+1],save_dXt[:,p+1]),Inf)
    return
end
