#functions to run PWL EVC problems

#central
function runEVCCentralStep(stepI,evS,cSol,cSave,silent)

    K=evS.K
    N=evS.N
    S=evS.S
    #deltaI=evS.deltaI
    #S=8
    deltaI=evS.ItotalMax/S
    horzLen=min(evS.K1,K-stepI)

    #initialize with current states
    global s0
    global t0

    #desired SOC
    target=zeros(N*(horzLen+1),1);
    for ii=1:N
       cur=evS.Kn[ii]-(stepI-1)
       ind=(max(0,(cur-1)*N)+ii):N:length(target)
       target[ind].=evS.Snmin[ii,1]
    end

    if !silent println("setting up model") end

    centralModel = Model(solver = GurobiSolver())
    #centralModel = Model(solver = GurobiSolver(Presolve=0,BarHomogeneous=1,NumericFocus=3))

    #u iD and z are one index ahead of sn and T. i.e the x[k+1]=x[k]+eta*u[k+1]
    @variable(centralModel,sn[1:(N)*(horzLen+1)])
    @variable(centralModel,u[1:(N)*(horzLen+1)])
    @variable(centralModel,t[1:(horzLen+1)])
    @variable(centralModel,z[1:S*(horzLen+1)])
    if slack
        @variable(centralModel,slackSn[1:N])
        #@variable(centralModel,slackT[1:(horzLen+1)])
    end
    if !silent println("obj") end
    #objExp =sum((sn[n,1]-1)^2*evS.Qsi[n,stepI]+(u[n,1])^2*evS.Ri[n,stepI] for n=1:N)
    objExp=0*u[1,1]^2
    for k=1:horzLen+1
        append!(objExp,sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,stepI+k-1]+(u[(k-1)*N+n,1])^2*evS.Ri[n,stepI+k-1]  for n=1:N))
    end
    if slack append!(objExp,sum(1e3*evS.β[n]*slackSn[n]^2 for n=1:N)) end
    if tempAugment append!(objExp,ψ*sum((evS.Tmax-t[k]) for k=1:horzLen+1)) end
    @objective(centralModel,Min,objExp)

    if !silent println("constraints") end
    @constraint(centralModel,stateCon1,sn[1:N,1].==s0[1:N,1]+evS.ηP[:,1].*u[1:N,1])
    @constraint(centralModel,stateCon2[k=1:horzLen,n=1:N],sn[n+(k)*(N),1]==sn[n+(k-1)*(N),1]+evS.ηP[n,1]*u[n+(k)*(N),1])
    @constraint(centralModel,tempCon1,t[1,1]==evS.τP*t0+evS.γP*deltaI*sum((2*s-1)*z[s,1] for s=1:S)+evS.ρP*evS.Tamb[stepI,1])
    @constraint(centralModel,tempCon2[k=1:horzLen],t[k+1,1]==evS.τP*t[k,1]+evS.γP*deltaI*sum((2*s-1)*z[k*S+s,1] for s=1:S)+evS.ρP*evS.Tamb[stepI+k,1])
    @constraint(centralModel,currCon[k=1:horzLen+1],0==-sum(u[(k-1)*(N)+n] for n=1:N)-evS.iD_pred[stepI+(k-1)]+sum(z[(k-1)*(S)+s] for s=1:S))
    @constraint(centralModel,sn.<=1)
    if slack
        @constraint(centralModel,sn.>=target.*(1-repeat(slackSn,horzLen+1,1)))
        @constraint(centralModel,slackSn.>=0)
    else
        @constraint(centralModel,sn.>=target)
    end
    if noTlimit==false
    	@constraint(centralModel,upperTCon,t.<=evS.Tmax)
        # @constraint(centralModel,upperTCon,t.<=evS.Tmax+slackT)
        # @constraint(centralModel,slackT.>=0)
    end
    @constraint(centralModel,t.>=0)
    @constraint(centralModel,upperCCon,u.<=repeat(evS.imax,horzLen+1,1))
    @constraint(centralModel,u.>=repeat(evS.imin,horzLen+1,1))
    @constraint(centralModel,z.>=0)
    @constraint(centralModel,z.<=deltaI)

    if !silent println("solving....") end

    if solverSilent
        @suppress_out begin
            status = solve(centralModel)
        end
    else
        status = solve(centralModel)
    end
    @assert status==:Optimal "Central optimization not solved to optimality"

    uRaw=getvalue(u)
    snRaw=getvalue(sn)
    tRaw=getvalue(t)
    zRaw=getvalue(z)
    if slack
        slackSnRaw=getvalue(slackSn)
    else
        slackSnRaw=zeros(N)
    end

    #calculate actual temp
    Tactual=zeros(horzLen+1,1)
    Tactual_background=zeros(horzLen+1,1)
    itotal=zeros(horzLen+1,1)
    itotal_background=zeros(horzLen+1,1)
    for k=1:horzLen+1
        itotal[k,1]=sum((uRaw[(k-1)*N+n,1]) for n=1:N) + iD_actual[stepI+(k-1),1]
        itotal_background[k,1]= iD_actual[stepI+(k-1),1]
    end
    Tactual[1,1]=evS.τP*t0+evS.γP*itotal[1,1]^2+evS.ρP*evS.Tamb[stepI,1]
    Tactual_background[1,1]=evS.τP*t0+evS.γP*itotal_background[1,1]^2+evS.ρP*evS.Tamb[stepI,1]
    for k=1:horzLen
        Tactual[k+1,1]=evS.τP*Tactual[k,1]+evS.γP*itotal[k+1,1]^2+evS.ρP*evS.Tamb[stepI+k,1]
        Tactual_background[k+1,1]=evS.τP*Tactual_background[k,1]+evS.γP*itotal_background[k+1,1]^2+evS.ρP*evS.Tamb[stepI+k,1]
    end

    lambdaCurr=-getdual(currCon)
    if noTlimit==false
    	lambdaUpperT=-getdual(upperTCon)
    else
    	lambdaUpperT=zeros(horzLen+1,1)
    end

    uSum=zeros(horzLen+1,1)
    zSum=zeros(horzLen+1,1)
    for k=1:horzLen+1
        zSum[k,1]=sum(zRaw[(k-1)*(S)+s,1] for s=1:S)
        uSum[k,1]=sum(uRaw[(k-1)*N+n,1] for n=1:N)
    end

    cSol.Tpred[stepI,1]=tRaw[1]
    cSol.Un[stepI,:]=uRaw[1:N]
    cSol.Sn[stepI,:]=snRaw[1:N]
    cSol.Z[stepI,:]=zRaw[1:S]
    cSol.uSum[stepI,1]=uSum[1]
    cSol.zSum[stepI,1]=zSum[1]
    cSol.Itotal[stepI,1]=itotal[1]
    cSol.Tactual[stepI,1]=Tactual[1]
    cSol.lamCoupl[stepI,1]=lambdaCurr[1]
    cSol.lamTemp[stepI,1]=lambdaUpperT[1]
    cSol.timeSolve[stepI,1]=getsolvetime(centralModel)
    #cSol.objVal[1,1,stepI]=getobjectivevalue(centralModel)

    if stepI in saveLogInd
        ind=findall(x->x==stepI,saveLogInd)[1]
        cSave.Obj[1,1,ind]=getobjectivevalue(centralModel)

        zReshape=zeros(horzLen+1,S)
        uReshape=zeros(horzLen+1,N)
        for ii= 1:N
            uReshape[:,ii]=uRaw[collect(ii:N:length(uRaw))]
        end
        for ii= 1:S
            zReshape[:,ii]=zRaw[collect(ii:S:length(zRaw))]
        end

        cSave.Un[1:(horzLen+1),:,ind]=uReshape
        #cSave.uSum[1:(horzLen+1),:,ind]=uSum
        #cSave.Itotal[1:(horzLen+1),:,ind]=itotal
        cSave.Lam[1:(horzLen+1),:,ind]=lambdaCurr
        cSave.Z[1:(horzLen+1),:,ind]=zReshape
        #cSave.zSum[1:(horzLen+1),:,ind]=zSum
        #cSave.Sn[1:(horzLen+1),:,ind]=snRaw
        cSave.Tactual[1:(horzLen+1),:,ind]=Tactual
        #cSave.Tpred[1:(horzLen+1),:,ind]=tRaw
    end

    # new states
    t0=round.(cSol.Tactual[stepI,1],sigdigits=roundSigFigs)
    s0=round.(cSol.Sn[stepI,:],sigdigits=roundSigFigs)
    return nothing
end

function pwlEVcentral(evS::scenarioStruct,slack::Bool,silent::Bool)

    horzLen=evS.K1
    K=evS.K
    N=evS.N
    S=evS.S
    #S=8
    cSol=solutionStruct(K=K,N=N,S=S)
    cSave=centralLogStruct(logLength=length(saveLogInd),horzLen=horzLen,N=N,S=S)

    stepI=1
    # using Profile
    # Profile.clear()
    # runEVCCentralStep(stepI,evS,cSol,silent)
    # cSol=centralSolutionStruct(horzLen=horzLen,K=K,N=N,S=S)
    # @profile  runEVCCentralStep(stepI,evS,cSol,silent)
    # Juno.profiler()

    Juno.progress() do id
        for stepI=1:K
            @info "$(Dates.format(Dates.now(),"HH:MM:SS")): $(stepI) of $(K)" progress=stepI/K _id=id
            @printf "%s: time step %g of %g....\n" Dates.format(Dates.now(),"HH:MM:SS") stepI K
            try
                cSol.timeT[stepI]=@elapsed runEVCCentralStep(stepI,evS,cSol,cSave,silent)
            catch e
                @printf "error: %s" e
                break
            end
        end
    end

    objFun(sn,u)=sum(sum((sn[k,n]-1)^2*evS.Qsi[n,k] for n=1:N) +sum((u[k,n])^2*evS.Ri[n,k] for n=1:N) for k=1:evS.K)
    cSol.objVal[1,1]=objFun(cSol.Sn,cSol.Un)

    return cSol, cSave
end

#dual
function localEVDual(evInd::Int,p::Int,stepI::Int,evS::scenarioStruct,dLog::itLogPWL,itLam,s0,slack,roundSigFigs,solverSilent)
    N=evS.N
    S=evS.S
    horzLen=min(evS.K1,evS.K-stepI)

    target=zeros((horzLen+1),1)
    target[max(1,(evS.Kn[evInd,1]-(stepI-1))):1:length(target),1].=evS.Snmin[evInd,1]
    # evM=Model(solver = GurobiSolver(NumericFocus=3))
    evM=Model(solver = GurobiSolver())
    @variable(evM,un[1:horzLen+1])
    @variable(evM,sn[1:horzLen+1])
    if slack @variable(evM,slackSn) end
    objExp=sum((sn[k,1]-1)^2*evS.Qsi[evInd,stepI+k-1]+(un[k,1])^2*evS.Ri[evInd,stepI+k-1]+itLam[k,1]*un[k,1] for k=1:horzLen+1)
    if slack
        append!(objExp,sum(evS.β[n]*slackSn^2 for n=1:N))
    end
    @objective(evM,Min,objExp)
    @constraint(evM,sn[1,1]==s0[evInd,1]+evS.ηP[evInd,1]*un[1,1]) #fix for MPC loop
    @constraint(evM,[k=1:horzLen],sn[k+1,1]==sn[k,1]+evS.ηP[evInd,1]*un[k+1,1]) #check K+1
    @constraint(evM,sn.<=1)
    if slack
        @constraint(evM,sn.>=target.*(1-slackSn))
        @constraint(evM,slackSn>=0)
    else
        @constraint(evM,sn.>=target)
    end
    @constraint(evM,un.<=evS.imax[evInd,1])
    @constraint(evM,un.>=evS.imin[evInd,1])

    if solverSilent
        @suppress_out begin
            statusEVM = solve(evM)
        end
    else
        statusEVM = solve(evM)
    end

    @assert statusEVM==:Optimal "EV NLP optimization not solved to optimality"

    dLog.Sn[collect(evInd:N:N*(horzLen+1)),p]=round.(getvalue(sn),sigdigits=roundSigFigs) #solved state goes in next time slot
    dLog.Un[collect(evInd:N:N*(horzLen+1)),p]=round.(getvalue(un),sigdigits=roundSigFigs) #current go
    dLog.slackSn[evInd]= if slack getvalue(slackSn) else 0 end
    dLog.timeSolve[1,p]=max(getsolvetime(evM),dLog.timeSolve[1,p])
    return nothing
end

function localXFRMDual(p::Int,stepI::Int,evS::scenarioStruct,dLog::itLogPWL,itLam,roundSigFigs)
    N=evS.N
    S=evS.S
    horzLen=min(evS.K1,evS.K-stepI)

    if updateMethod=="dualAscent"
        #solve coordinator problem
        #coorM=Model(solver = GurobiSolver(Presolve=0,NumericFocus=1))
        coorM=Model(solver = GurobiSolver())
        #coorM=Model(solver = IpoptSolver())
        @variable(coorM,z[1:S*(horzLen+1)])
        @variable(coorM,t[1:horzLen+1])
        @objective(coorM,Min,sum(itLam[k,1]*sum(-z[(k-1)*S+s,1] for s=1:S) for k=1:(horzLen+1)))
        @constraint(coorM,t[1,1]==evS.τP*t0+evS.γP*evS.deltaI*sum((2*s-1)*z[s,1] for s=1:S)+evS.ρP*evS.Tamb[stepI,1]) #fix for MPC loop
        @constraint(coorM,[k=1:horzLen],t[k+1,1]==evS.τP*t[k,1]+evS.γP*evS.deltaI*sum((2*s-1)*z[k*S+s,1] for s=1:S)+evS.ρP*evS.Tamb[stepI+k,1])
        if noTlimit==false
            @constraint(coorM,upperTCon,t.<=evS.Tmax)
        end
        @constraint(coorM,t.>=0)
        @constraint(coorM,z.<=evS.deltaI)
        @constraint(coorM,z.>=0)

        if solverSilent
            @suppress_out begin
                statusC = solve(coorM)
            end
        else
            statusC = solve(coorM)
        end

        @assert statusC==:Optimal "Dual Ascent XFRM optimization not solved to optimality"

        dLog.Tpred[:,p]=round.(getvalue(t),sigdigits=roundSigFigs)
        dLog.Z[:,p]=round.(getvalue(z),sigdigits=roundSigFigs)
        dLog.timeSolve[1,p]=max(getsolvetime(coorM),dLog.timeSolve[1,p])

        #grad of lagragian
        for k=1:horzLen+1
            dLog.zSum[k,p]=sum(dLog.Z[(k-1)*(S)+s,p] for s=1:S)
            #dLog.uSum[k,p]=sum(dLog.Un[(k-1)*N+n,p] for n=1:N)
            #dLog.couplConst[k,p]=dLog.uSum[k,p] + evS.iD_pred[stepI+(k-1),1] - dLog.zSum[k,p]
        end
    end

    return nothing
end

function runEVDualIt(p,stepI,evS,itLam,dLog,dCM,dSol,cSave,roundSigFigs,silent)
    K=evS.K
    N=evS.N
    S=evS.S
    horzLen=min(evS.K1,K-stepI)

    #initialize with current states
    global s0
    global t0
    #global prevLam

    if updateMethod=="fastAscent"
        #alpha0 = 0.1  #for A
        #alpha0 = 5e4 #for kA
        # alphaDivRate=2
        # minAlpha=1e-6
    else
        #alpha0 = 3e-3 #for A
        # alpha = 3e4 #for kA
        # alphaDivRate=4
        #alpha0 = 2e6 #for kA
        #alphaDivRate=2
        # minAlpha=1e-6
    end
    #solve subproblem for each EV
    if runParallel
        @sync @distributed for evInd=1:N
            localEVDual(evInd,p,stepI,evS,dLog,itLam,s0,slack,roundSigFigs,solverSilent)
        end
    else
        for evInd=1:N
            localEVDual(evInd,p,stepI,evS,dLog,itLam,s0,slack,roundSigFigs,solverSilent)
        end
    end

    localXFRMDual(p,stepI,evS,dLog,itLam,roundSigFigs)

    for k=1:horzLen+1
        dLog.uSum[k,p]=sum(dLog.Un[(k-1)*N+n,p] for n=1:N)
    end
    #calculate actual temperature from nonlinear model of XFRM
    for k=1:horzLen+1
        dLog.Itotal[k,p]=dLog.uSum[k,p] + evS.iD_actual[stepI+(k-1),1]
    end
    dLog.Tactual[1,p]=evS.τP*t0+evS.γP*dLog.Itotal[1,p]^2+evS.ρP*evS.Tamb[stepI,1]
    for k=1:horzLen
        dLog.Tactual[k+1,p]=evS.τP*dLog.Tactual[k,p]+evS.γP*dLog.Itotal[k+1,p]^2+evS.ρP*evS.Tamb[stepI+k,1]  #fix for mpc
    end

    if updateMethod=="fastAscent"
        #fast ascent
        if noTlimit==false
            dLog.couplConst[:,p]=dLog.Tactual[:,p]-evS.Tmax*ones(horzLen+1,1)
        else
            dLog.couplConst[:,p]=zeros(horzLen+1,1)
        end
        #add some amount of future lambda
        for k=1:(horzLen+1-2)
            dLog.couplConst[k,p]=.6*dLog.couplConst[k,p]+.3*dLog.couplConst[k+1,p]+.1*dLog.couplConst[k+2,p]
            #gradL[k,1]=.5*gradL[k,1]+.2*gradL[k+1,1]+.2*gradL[k+2,1]+.1*gradL[k+3,1]+.1*gradL[k+4,1]
        end
    else
        for k=1:horzLen+1
            dLog.couplConst[k,p]=dLog.uSum[k,p] + evS.iD_pred[stepI+(k-1),1] - dLog.zSum[k,p]
        end
    end


    #update lambda
    alphaP= max(alpha0/ceil(p/alphaDivRate),minAlpha)
    #alphaP= max(alpha0/(p/alphaDivRate),minAlpha)
    dLog.itUpdate[1,p]=alphaP
    #alphaP= alphaP*alphaRate

    #dLog.Lam[:,p]=max.(itLam[:,1]+alphaP*dLog.couplConst[:,p],0)
    dLog.Lam[:,p]=round.(itLam[:,1]+alphaP*dLog.couplConst[:,p],sigdigits=roundSigFigs)

    #check convergence
    objFun(sn,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,stepI+k-1]     for n=1:N) for k=1:horzLen+1) +
                    sum(sum((u[(k-1)*N+n,1])^2*evS.Ri[n,stepI+k-1]           for n=1:N) for k=1:horzLen+1)
    dLog.objVal[1,p]=objFun(dLog.Sn[:,p],dLog.Un[:,p])
    itGap = norm(dLog.Lam[:,p]-itLam[:,1],2)
    #itGap = norm((dLog.Lam[:,p] .- itLam[:,1])./itLam[:,1],2)
    couplGap=norm(dLog.couplConst[:,p],1)

    if stepI in saveLogInd
        ind=findall(x->x==stepI,saveLogInd)[1]
        dCM.coupl1Norm[p,ind]=couplGap
        dCM.lamIt2Norm[p,ind]=itGap
        dCM.objAbs[p,ind]=abs(dLog.objVal[1,p]-cSave.Obj[1,1,ind])
        dCM.objPerc[p,ind]=abs(dLog.objVal[1,p]-cSave.Obj[1,1,ind])/cSave.Obj[1,1,ind]
        prevObj = if p==1 dLog.objVal[1,p] else dLog.objVal[1,p-1] end
        dCM.objItPerc[p,ind]=abs(dLog.objVal[1,p]-prevObj)/prevObj
        dCM.lam1Norm[p,ind]= norm(dLog.Lam[:,p]-cSave.Lam[1:(horzLen+1),:,ind],1)
        dCM.lam2Norm[p,ind]= norm(dLog.Lam[:,p]-cSave.Lam[1:(horzLen+1),:,ind],2)
        dCM.lamInfNorm[p,ind]= norm(dLog.Lam[:,p]-cSave.Lam[1:(horzLen+1),:,ind],Inf)
        dCM.t1Norm[p,ind]= norm(dLog.Tactual[:,p]-cSave.Tactual[1:(horzLen+1),:,ind],1)
        dCM.t2Norm[p,ind]= norm(dLog.Tactual[:,p]-cSave.Tactual[1:(horzLen+1),:,ind],2)
        dCM.tInfNorm[p,ind]= norm(dLog.Tactual[:,p]-cSave.Tactual[1:(horzLen+1),:,ind],Inf)
        zReshape=zeros(horzLen+1,S)
        uReshape=zeros(horzLen+1,N)
        for ii= 1:N
            uReshape[:,ii]=dLog.Un[collect(ii:N:length(dLog.Un[:,p])),p]
        end
        for ii= 1:S
            zReshape[:,ii]=dLog.Z[collect(ii:S:length(dLog.Z[:,p])),p]
        end
        dCM.un1Norm[p,ind]= norm(uReshape-cSave.Un[1:(horzLen+1),:,ind],1)
        dCM.un2Norm[p,ind]= norm(uReshape-cSave.Un[1:(horzLen+1),:,ind],2)
        dCM.unInfNorm[p,ind]= norm(uReshape-cSave.Un[1:(horzLen+1),:,ind],Inf)
        dCM.z1Norm[p,ind]= norm(zReshape-cSave.Z[1:(horzLen+1),:,ind],1)
        dCM.z2Norm[p,ind]= norm(zReshape-cSave.Z[1:(horzLen+1),:,ind],2)
        dCM.zInfNorm[p,ind]= norm(zReshape-cSave.Z[1:(horzLen+1),:,ind],Inf)
    end

    # if((couplGap<=primChk) && (itGap <= dualChk))
    if(couplGap<=primChk)
        if !silent @printf "Converged after %g iterations\n" p end
        @printf "Y"
        convIt=p
        return true
    else
        if !silent
            #@printf "lastGap  %e after %g iterations\n" itGap p
            @printf "couplGap %e after %g iterations\n" couplGap p
            #@printf "convGap %e after %g iterations\n\n" convGap p
            #@printf "snGap   %e after %g iterations\n" snGap p
            #@printf "unGap   %e after %g iterations\n" unGap p
            #@printf("fGap    %e after %g iterations\n\n",fGap,p)
        end
        #prevLam=dLog.Lam[:,p]
        return false
    end
end

function runEVDualStep(stepI,maxIt,evS,dSol,dCM,cSave,roundSigFigs,silent)

    K=evS.K
    N=evS.N
    S=evS.S
    horzLen=min(evS.K1,K-stepI)
    #u iD and z are one index ahead of sn and T. i.e the x[k+1]=x[k]+eta*u[k+1]
    dLog=itLogPWL(horzLen=horzLen,N=N,S=S)
    p=1
    timeStart=now()
    while (p<=maxIt && round(now()-timeStart,Second)<=Dates.Second(9/10*evS.Ts))
        #while (p<maxIt)
        #global p
        if p==1
            itLam=prevLam
            #global alpha0=max(min(maximum(prevLam)/5,1e6),1e-3)
        else
            itLam=round.(dLog.Lam[:,(p-1)],sigdigits=roundSigFigs)
        end
        cFlag=runEVDualIt(p,stepI,evS,itLam,dLog,dCM,dSol,cSave,roundSigFigs,silent)
        global convIt=p
        if cFlag
            break
        end
        p+=1
    end
    if stepI in saveLogInd
        ind=findall(x->x==stepI,saveLogInd)[1]
        dCM.convIt[1,1,ind]=convIt
    end

    # # # # convergence plots
    # halfCI=Int(floor(convIt/2))
    # CList=reshape([range(colorant"blue", stop=colorant"yellow",length=halfCI);
    #                range(colorant"yellow", stop=colorant"red",length=convIt-halfCI)], 1, convIt);
    #
    # uSumPlotd=plot(dLog.uSum[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Current Sum (kA)",xlims=(0,horzLen+1),legend=false)
    # plot!(uSumPlotd,sum(cSave.Un[:,:,ind],dims=2),seriescolor=:black,linewidth=2,linealpha=0.8)
    # plot!(uSumPlotd,cSol.uSum[stepI:stepI+horzLen],seriescolor=:black,linewidth=2,linealpha=0.8)
    # plot!(uSumPlotd,uSum,seriescolor=:black,linewidth=2,linealpha=0.8)
    #
    # zSumPlotd=plot(dLog.zSum[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Z sum",xlims=(0,horzLen+1),legend=false)
    # plot!(zSumPlotd,sum(cSave.Z[:,:,ind],dims=2),seriescolor=:black,linewidth=2,linealpha=0.8)
    #
    # lamPlotd=plot(dLog.Lam[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Lambda",xlims=(0,horzLen+1),legend=false)
    # plot!(lamPlotd,cSave.Lam[:,:,ind],seriescolor=:black,linewidth=2,linealpha=0.8)
    # plot!(lamPlotd,cSol.lamCoupl[stepI:stepI+horzLen],seriescolor=:black,linewidth=2,linealpha=0.8)
    # plot!(lamPlotd,lambdaCurr,seriescolor=:black,linewidth=2,linealpha=0.8)
    #
    # constPlot2=plot(dLog.couplConst[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="curr constraint diff",xlims=(0,horzLen+1),legend=false)
    #
    # tempPlot=plot(dLog.Tactual[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Temp",xlims=(0,horzLen+1),legend=false)
    #
    # #convergence metric plots
    # fAbsPlot=plot(1:convIt,dCM.objAbs[1:convIt,ind],xlabel="Iteration",ylabel="obj function gap",xlims=(1,convIt),legend=false,yscale=:log10)
    # fPlot=plot(1:convIt,dCM.objPerc[2:convIt-1,ind],xlabel="Iteration",ylabel="obj function gap",xlims=(2,convIt-1),legend=false,yscale=:log10)
    # fItPlot=plot(1:convIt,dCM.objItPerc[2:convIt-1,ind],xlabel="Iteration",ylabel="obj function gap",xlims=(2,convIt),legend=false,yscale=:log10)
    # convItPlot=plot(1:convIt,dCM.lamIt2Norm[1:convIt,1],xlabel="Iteration",ylabel="2-Norm Lambda Gap",xlims=(1,convIt),legend=false,yscale=:log10)
    # convLamPlot=plot(1:convIt,dCM.lam2Norm[1:convIt,1],xlabel="Iteration",ylabel="central lambda gap",xlims=(2,convIt),legend=false,yscale=:log10)
    # constPlot=plot(1:convIt,dCM.coupl1Norm[1:convIt,1],xlabel="Iteration",ylabel="curr constraint Gap",xlims=(2,convIt),legend=false,yscale=:log10)
    # #savefig(lamPlotd,path*"lamPlotd"*".png")
    # testRel2 = zeros(convIt-1)
    # testRelInf = zeros(convIt-1)
    # testV = zeros(convIt-1)
    # testV2 = zeros(convIt-1)
    # testZ = zeros(convIt-1)
    # for i in 2:convIt
    #     v1=dLog.Lam[:,i]
    #     v2=dLog.Lam[:,i-1]
    #     testV[i-1]=1/2*(var(v1-v2))/(var(v1)+var(v2))
    #     testV2[i-1]=sqrt(sum((v1[k]-v2[k])^2/(var(v1)+var(v2)) for k=1:horzLen+1))
    #     testRel2[i-1] = norm(v1-v2,2)/norm(v2,2)
    #     testRelInf[i-1] = norm(v1-v2,Inf)/norm(v2,Inf)
    #
    #     z1=(v1 .- mean(v1))./std(v1)
    #     z2=(v2 .- mean(v2))./std(v2)
    #     testZ[i-1]=norm(z1-z2,2)
    # end
    # plot(testRel2,yscale=:log10,legend=false,xlabel="Iteration")
    # plot(testRelInf,yscale=:log10,legend=false,xlabel="Iteration")
    # # plot(testV,yscale=:log10,legend=false,xlabel="Iteration")
    # # plot(testV2,yscale=:log10,legend=false,xlabel="Iteration")
    # # plot(testZ,yscale=:log10,legend=false,xlabel="Iteration")

    #save current state and update for next timeSteps
    dSol.Tpred[stepI,1]=dLog.Tpred[1,convIt]
    dSol.Un[stepI,:]=dLog.Un[1:N,convIt]
    dSol.Sn[stepI,:]=dLog.Sn[1:N,convIt]
    dSol.Z[stepI,:]=dLog.Z[1:S,convIt]
    dSol.uSum[stepI,1]=dLog.uSum[1,convIt]
    dSol.zSum[stepI,1]=dLog.zSum[1,convIt]
    dSol.Itotal[stepI,1]=dLog.Itotal[1,convIt]
    dSol.Tactual[stepI,1]=dLog.Tactual[1,convIt]
    dSol.convIt[stepI,1]=convIt
    dSol.timeSolve[stepI,1]=mean(dLog.timeSolve[1,1:convIt])

    # new states
    global t0=round.(dSol.Tactual[stepI,1],sigdigits=roundSigFigs)
    global s0=round.(dSol.Sn[stepI,:],sigdigits=roundSigFigs)

    if convIt==1
        dSol.lamCoupl[stepI,1]=prevLam[1,1]
        if stepI+horzLen==evS.K
            newLam=prevLam[2:horzLen+1,1]
        else
            newLam=vcat(prevLam[2:horzLen+1,1],prevLam[horzLen+1,1])
        end
    else
        dSol.lamCoupl[stepI,1]=dLog.Lam[1,convIt-1]
        if stepI+horzLen==evS.K
            newLam=dLog.Lam[2:horzLen+1,convIt-1]
        else
            newLam=vcat(dLog.Lam[2:horzLen+1,convIt-1],dLog.Lam[horzLen+1,convIt-1])
        end
    end
    global prevLam=round.(newLam,sigdigits=roundSigFigs)
    return nothing
end

function pwlEVdual(maxIt::Int,updateMethod::String,evS::scenarioStruct,cSave::centralLogStruct,slack::Bool,roundSigFigs::Int,silent::Bool)

    horzLen=evS.K1
    K=evS.K
    N=evS.N
    S=evS.S
    dSol=solutionStruct(K=K,N=N,S=S)
    dCM=convMetricsStruct(maxIt=maxIt,logLength=length(saveLogInd))

    Juno.progress() do id
        for stepI=1:K
            @info "$(Dates.format(Dates.now(),"HH:MM:SS")): $(stepI) of $(K)" progress=stepI/K _id=id
            @printf "%s: time step %g of %g...." Dates.format(Dates.now(),"HH:MM:SS") stepI K
            try
                dSol.timeT[stepI]=@elapsed runEVDualStep(stepI,maxIt,evS,dSol,dCM,cSave,roundSigFigs,silent)
                @printf "convIt: %g\n" dSol.convIt[stepI,1]
            catch e
                @printf "error: %s" e
                break
            end
        end
    end

    objFun(sn,u)=sum(sum((sn[k,n]-1)^2*evS.Qsi[n,k] for n=1:N) +sum((u[k,n])^2*evS.Ri[n,k] for n=1:N) for k=1:evS.K)
    dSol.objVal[1,1]=objFun(dSol.Sn,dSol.Un)

    return dSol, dCM
end

#ADMM
function localEVADMM(evInd::Int,p::Int,stepI::Int,evS,dLogadmm::itLogPWL,
    evV,itLam,s0,itρ,slack,roundSigFigs,solverSilent)
    N=evS.N
    S=evS.S
    horzLen=min(evS.K1,evS.K-stepI)

    #evV=prevVu[collect(evInd:N:length(prevVu)),1]
    target=zeros((horzLen+1),1)
    target[max(1,(evS.Kn[evInd,1]-(stepI-1))):1:length(target),1].=evS.Snmin[evInd,1]
    evM = Model(solver = GurobiSolver(NumericFocus=3))
    @variable(evM,sn[1:(horzLen+1)])
    @variable(evM,u[1:(horzLen+1)])
    if slack @variable(evM,slackSn) end
    objExp= sum((sn[k,1]-1)^2*evS.Qsi[evInd,stepI+k-1]+(u[k,1])^2*evS.Ri[evInd,stepI+k-1]+
                            itLam[k,1]*(u[k,1]-evV[k,1])+
                            itρ/2*(u[k,1]-evV[k,1])^2 for k=1:horzLen+1)
    # objExp= sum((sn[k,1]-1)^2*evS.Qsi[evInd,1]+(u[k,1])^2*evS.Ri[evInd,1]+
    #                         itLam[k,1]*(u[k,1]-evV[k,1])+
    #                         itρ/2*(u[k,1]-evV[k,1])^2 for k=1:horzLen+1)
    if slack
        append!(objExp,sum(evS.β[n]*slackSn^2 for n=1:N))
    end
    @objective(evM,Min,objExp)
    @constraint(evM,sn[1,1]==s0[evInd,1]+evS.ηP[evInd,1]*u[1,1])
    @constraint(evM,[k=1:horzLen],sn[k+1,1]==sn[k,1]+evS.ηP[evInd,1]*u[k+1,1])
    @constraint(evM,sn.<=1)
    if slack
        @constraint(evM,sn.>=target.*(1-slackSn))
        @constraint(evM,slackSn>=0)
    else
        @constraint(evM,sn.>=target)
    end
    @constraint(evM,u.<=evS.imax[evInd,1])
    @constraint(evM,u.>=evS.imin[evInd,1])
    if solverSilent
        @suppress_out begin
            statusEVM = solve(evM)
        end
    else
        statusEVM = solve(evM)
    end
    @assert statusEVM==:Optimal "ADMM EV NLP optimization not solved to optimality"

    dLogadmm.Sn[collect(evInd:N:length(dLogadmm.Sn[:,p])),p]=round.(getvalue(sn),sigdigits=roundSigFigs)
    dLogadmm.Un[collect(evInd:N:length(dLogadmm.Un[:,p])),p]=round.(getvalue(u),sigdigits=roundSigFigs)
    dLogadmm.slackSn[evInd]= if slack getvalue(slackSn) else 0 end
    dLogadmm.timeSolve[1,p]=max(getsolvetime(evM),dLogadmm.timeSolve[1,p])

    return nothing
end

function localXFRMADMM(p::Int,stepI::Int,evS,dLogadmm::itLogPWL,itLam,itVz,itρ,roundSigFigs,solverSilent)
    N=evS.N
    S=evS.S
    horzLen=min(evS.K1,evS.K-stepI)

    #N+1 decoupled problem aka transformer current
    tM = Model(solver = GurobiSolver(NumericFocus=3))
    @variable(tM,z[1:(S)*(horzLen+1)])
    @variable(tM,t[1:(horzLen+1)])
    # constFun1(u,v)=sum(Lam[k,p]*sum(u[(k-1)*(S)+s,1]-v[(k-1)*(S)+s,1] for s=1:S)  for k=1:(horzLen+1))
    # constFun2(u,v)=ρADMMp/2*sum(sum((u[(k-1)*(S)+s,1]-v[(k-1)*(S)+s,1])*(u[(k-1)*(S)+s,1]-v[(k-1)*(S)+s,1]) for s=1:S)  for k=1:(horzLen+1))
    # @objective(tM,Min, constFun1(-z,Vz[:,p])+constFun2(-z,Vz[:,p]))
    @objective(tM,Min,sum(itLam[k,1]*(sum(-z[(k-1)*(S)+s,1] for s=1:S)-sum(itVz[(k-1)*(S)+s,1] for s=1:S)) +
                        itρ/2*(sum(-z[(k-1)*(S)+s,1] for s=1:S)-sum(itVz[(k-1)*(S)+s,1] for s=1:S))^2  for k=1:(horzLen+1)))
    @constraint(tM,tempCon1,t[1,1]==evS.τP*t0+evS.γP*evS.deltaI*sum((2*m+1)*z[m+1,1] for m=0:S-1)+evS.ρP*evS.Tamb[stepI,1])
    @constraint(tM,tempCon2[k=1:horzLen],t[k+1,1]==evS.τP*t[k,1]+evS.γP*evS.deltaI*sum((2*m+1)*z[k*S+(m+1),1] for m=0:S-1)+evS.ρP*evS.Tamb[stepI+k,1])
    if noTlimit==false
            @constraint(tM,upperTCon,t.<=evS.Tmax)
    end
    @constraint(tM,t.>=0)
    @constraint(tM,z.>=0)
    @constraint(tM,z.<=evS.deltaI)
    #@constraint(tM,zC[k=1:horzLen+1],zSum[k,1]==sum(z[(k-1)*(S)+s] for s=1:S))

    if solverSilent
        @suppress_out begin
            statusC = solve(tM)
        end
    else
        statusC = solve(tM)
    end

    @assert statusC==:Optimal "ADMM XFRM NLP optimization not solved to optimality"

    dLogadmm.Tpred[:,p]=round.(getvalue(t),sigdigits=roundSigFigs)
    dLogadmm.Z[:,p]=round.(getvalue(z),sigdigits=roundSigFigs)
    dLogadmm.timeSolve[1,p]=max(getsolvetime(tM),dLogadmm.timeSolve[1,p])

    return nothing
end

function runEVADMMIt(p,stepI,evS,itLam,itVu,itVz,itρ,dLogadmm,dCM,dSol,cSave,roundSigFigs,silent)
    #@printf "a"
    K=evS.K
    N=evS.N
    S=evS.S
    horzLen=min(evS.K1,K-stepI)

    #initialize with current states
    global s0
    global t0
    # global prevLam
    # global prevVu
    # global prevVz
    # global ρADMMp

    #ρDivRate=1.1
    #ρRate=1.5

    #x minimization eq 7.66 in Bertsekas
    if runParallel
        @sync @distributed for evInd=1:N
            evV=itVu[collect(evInd:N:length(itVu)),1]
            localEVADMM(evInd,p,stepI,evS,dLogadmm,evV,itLam,s0,itρ,slack,roundSigFigs,solverSilent)
        end
    else
        for evInd=1:N
            evV=itVu[collect(evInd:N:length(itVu)),1]
            localEVADMM(evInd,p,stepI,evS,dLogadmm,evV,itLam,s0,itρ,slack,roundSigFigs,solverSilent)
        end
    end
    #@printf "b"


    localXFRMADMM(p,stepI,evS,dLogadmm,itLam,itVz,itρ,roundSigFigs,solverSilent)
    #@printf "c"

    #lambda update eq 7.68
    for k=1:horzLen+1
        dLogadmm.uSum[k,p]=sum(dLogadmm.Un[(k-1)*N+n,p] for n=1:N)
        dLogadmm.zSum[k,p]=sum(dLogadmm.Z[(k-1)*(S)+s,p] for s=1:S)
        dLogadmm.couplConst[k,p]= dLogadmm.uSum[k,p] + evS.iD_pred[stepI+(k-1),1] - dLogadmm.zSum[k,p]
        dLogadmm.Lam[k,p]=round.(itLam[k,1]+itρ/(max(horzLen+1,N))*(dLogadmm.couplConst[k,p]),sigdigits=roundSigFigs)
        #dLogadmm.Lam[k,p]=max.(round.(itLam[k,1]+itρ/(max(horzLen+1,N))*(dLogadmm.couplConst[k,p]),sigdigits=roundSigFigs),0)
    end

    #calculate actual temperature from nonlinear model of XFRM
    for k=1:horzLen+1
        dLogadmm.Itotal[k,p]=dLogadmm.uSum[k,p] + evS.iD_actual[stepI+(k-1),1]
    end
    dLogadmm.Tactual[1,p]=evS.τP*t0+evS.γP*dLogadmm.Itotal[1,p]^2+evS.ρP*evS.Tamb[stepI,1] #fix for mpc
    for k=1:horzLen
        dLogadmm.Tactual[k+1,p]=evS.τP*dLogadmm.Tactual[k,p]+evS.γP*dLogadmm.Itotal[k+1,p]^2+evS.ρP*evS.Tamb[stepI+k,1]  #fix for mpc
    end

    #v upate eq 7.67
    for k=1:horzLen+1
        dLogadmm.Vu[(k-1)*N.+collect(1:N),p]=round.(dLogadmm.Un[(k-1)*N.+collect(1:N),p].+(itLam[k,1]-dLogadmm.Lam[k,p])/(itρ/1),sigdigits=roundSigFigs)
        dLogadmm.Vz[(k-1)*(S).+collect(1:S),p]=round.(-dLogadmm.Z[(k-1)*(S).+collect(1:S),p].+(itLam[k,1]-dLogadmm.Lam[k,p])/(itρ/1),sigdigits=roundSigFigs)

        # dLogadmm.Vu[(k-1)*N.+collect(1:N),p]=min.(max.(dLogadmm.Un[(k-1)*N.+collect(1:N),p].+(itLam[k,1]-dLogadmm.Lam[k,p])/itρ,evS.imin),evS.imax)
        # dLogadmm.Vz[(k-1)*(S).+collect(1:S),p]=max.(min.(-dLogadmm.Z[(k-1)*(S).+collect(1:S),p].+(itLam[k,1]-dLogadmm.Lam[k,p])/itρ,0),-evS.deltaI)
    end

    #check convergence
    cc=norm(vcat((itVu[:,1]-dLogadmm.Vu[:,p]),(itVz[:,1]-dLogadmm.Vz[:,p])),1)
    objFun(sn,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,stepI+k-1]     for n=1:N) +
                     sum((u[(k-1)*N+n,1])^2*evS.Ri[n,stepI+k-1]           for n=1:N) for k=1:horzLen+1)
    dLogadmm.objVal[1,p]=objFun(dLogadmm.Sn[:,p],dLogadmm.Un[:,p])

    constGap=norm(dLogadmm.couplConst[:,p],1)
    itGap = norm(dLogadmm.Lam[:,p]-itLam[:,1],2)

    #update rho
    #ρADMMp = ρADMM/ceil(p/ρDivRate)
    dLogadmm.itUpdate[1,p]= min(itρ*ρDivRate,maxRho)

    #Boyds varying penalty parameter
    # if norm(dLogadmm.couplConst[:,p],2)>10*itGap
    #     dLogadmm.itUpdate[1,p]=min(itρ*ρRate,maxRho)
    # elseif 10*norm(dLogadmm.couplConst[:,p],2)<itGap
    #     dLogadmm.itUpdate[1,p]=max(itρ/ρRate,1)
    # else
    #     dLogadmm.itUpdate[1,p]=itρ
    # end

    #only if saving
    if stepI in saveLogInd
        ind=findall(x->x==stepI,saveLogInd)[1]
        dCM.coupl1Norm[p,ind]=constGap
        dCM.lamIt2Norm[p,ind]=itGap
        dCM.objAbs[p,ind]=abs(dLogadmm.objVal[1,p]-cSave.Obj[1,1,ind])
        dCM.objPerc[p,ind]=abs(dLogadmm.objVal[1,p]-cSave.Obj[1,1,ind])/cSave.Obj[1,1,ind]
        prevObj = if p==1 dLogadmm.objVal[1,p] else dLogadmm.objVal[1,p-1] end
        dCM.objItPerc[p,ind]=abs(dLogadmm.objVal[1,p]-prevObj)/prevObj
        dCM.lam1Norm[p,ind]= norm(dLogadmm.Lam[:,p]-cSave.Lam[1:(horzLen+1),:,ind],1)
        dCM.lam2Norm[p,ind]= norm(dLogadmm.Lam[:,p]-cSave.Lam[1:(horzLen+1),:,ind],2)
        dCM.lamInfNorm[p,ind]= norm(dLogadmm.Lam[:,p]-cSave.Lam[1:(horzLen+1),:,ind],Inf)
        dCM.t1Norm[p,ind]= norm(dLogadmm.Tactual[:,p]-cSave.Tactual[1:(horzLen+1),:,ind],1)
        dCM.t2Norm[p,ind]= norm(dLogadmm.Tactual[:,p]-cSave.Tactual[1:(horzLen+1),:,ind],2)
        dCM.tInfNorm[p,ind]= norm(dLogadmm.Tactual[:,p]-cSave.Tactual[1:(horzLen+1),:,ind],Inf)
        zReshape=zeros(horzLen+1,S)
        uReshape=zeros(horzLen+1,N)
        for ii= 1:N
            uReshape[:,ii]=dLogadmm.Un[collect(ii:N:length(dLogadmm.Un[:,p])),p]
        end
        for ii= 1:S
            zReshape[:,ii]=dLogadmm.Z[collect(ii:S:length(dLogadmm.Z[:,p])),p]
        end
        dCM.un1Norm[p,ind]= norm(uReshape-cSave.Un[1:(horzLen+1),:,ind],1)
        dCM.un2Norm[p,ind]= norm(uReshape-cSave.Un[1:(horzLen+1),:,ind],2)
        dCM.unInfNorm[p,ind]= norm(uReshape-cSave.Un[1:(horzLen+1),:,ind],Inf)
        dCM.z1Norm[p,ind]= norm(zReshape-cSave.Z[1:(horzLen+1),:,ind],1)
        dCM.z2Norm[p,ind]= norm(zReshape-cSave.Z[1:(horzLen+1),:,ind],2)
        dCM.zInfNorm[p,ind]= norm(zReshape-cSave.Z[1:(horzLen+1),:,ind],Inf)
    end

    #if(constGap <= primChk  && itGap <=dualChk)
    if(constGap <= primChk)
        if !silent @printf "Converged after %g iterations\n" p end
        @printf "Y"
        convIt=p
        return true
    else
        if !silent
            #@printf "convGap  %e after %g iterations\n" convGap p
            #@printf "dual residue  %e after %g iterations\n" cc p
            #@printf "lastGap  %e after %g iterations\n" itGap p
            @printf "constGap %e after %g iterations\n" constGap p
            #@printf "snGap    %e after %g iterations\n" snGap p
            #@printf "unGap    %e after %g iterations\n" unGap p
            #@printf("fGap     %e after %g iterations\n\n",fGap,p)
        end
        return false
    end
end

function runEVADMMStep(stepI::Int,maxIt::Int,evS,dSol::solutionStruct,dCM,cSave::centralLogStruct,roundSigFigs::Int,silent::Bool)
    K=evS.K
    N=evS.N
    S=evS.S
    horzLen=min(evS.K1,K-stepI)
    #u iD and z are one index ahead of sn and T. i.e the x[k+1]=x[k]+η*u[k+1]
    dLogadmm=itLogPWL(horzLen=horzLen,N=N,S=S)
    p=1
    timeStart=now()
    while (p<=maxIt && round(now()-timeStart,Second)<=Dates.Second(9/10*evS.Ts))
        #while (p<=maxIt)
        #global p
        if p==1
            itLam=prevLam
            itVu=prevVu
            itVz=prevVz
            itρ=ρADMMp
        else
            itLam=round.(dLogadmm.Lam[:,(p-1)],sigdigits=roundSigFigs)
            itVu=round.(dLogadmm.Vu[:,(p-1)],sigdigits=roundSigFigs)
            itVz=round.(dLogadmm.Vz[:,(p-1)],sigdigits=roundSigFigs)
            itρ=round.(dLogadmm.itUpdate[1,(p-1)],sigdigits=roundSigFigs)
        end
        cFlag=runEVADMMIt(p,stepI,evS,itLam,itVu,itVz,itρ,dLogadmm,dCM,dSol,cSave,roundSigFigs,silent)
        global convIt=p
        if cFlag
            break
        end
        p+=1
    end
    if stepI in saveLogInd
        ind=findall(x->x==stepI,saveLogInd)[1]
        dCM.convIt[1,1,ind]=convIt
    end

    # # # #convergence plots
    # halfCI=Int(floor(convIt/2))
    # CList=reshape([range(colorant"blue", stop=colorant"yellow",length=halfCI);
    #                range(colorant"yellow", stop=colorant"red",length=convIt-halfCI)], 1, convIt);
    #
    # uSumPlotadmm=plot(dLogadmm.uSum[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Current Sum (kA)",xlims=(0,horzLen+1),legend=false)
    # plot!(uSumPlotadmm,sum(cSave.Un[:,:,ind],dims=2),seriescolor=:black,linewidth=2,linealpha=0.8)
    #
    # zSumPlotadmm=plot(dLogadmm.zSum[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Z sum",xlims=(0,horzLen+1),legend=false)
    # plot!(zSumPlotadmm,sum(cSave.Z[:,:,ind],dims=2),seriescolor=:black,linewidth=2,linealpha=0.8)
    #
    # constPlotadmm2=plot(dLogadmm.couplConst[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="curr constraint diff",xlims=(0,horzLen+1),legend=false)
    #
    # lamPlotadmm=plot(dLogadmm.Lam[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Lambda",xlims=(0,horzLen+1),legend=false)
    # plot!(lamPlotadmm,cSave.Lam[:,:,ind],seriescolor=:black,linewidth=2,linealpha=0.8)
    #
    # fAbsPlot=plot(dCM.objAbs[1:convIt,ind],xlabel="Iteration",ylabel="obj function gap",xlims=(1,convIt),legend=false,yscale=:log10)
    # fPlot=plot(dCM.objPerc[2:convIt-1,ind],xlabel="Iteration",ylabel="obj function gap",legend=false,yscale=:log10)
    # fItPlot=plot(dCM.objItPerc[2:convIt-1,ind],xlabel="Iteration",ylabel="obj function gap",legend=false,yscale=:log10)
    #
    # convPlotalad=plot(dCM.lam2Norm[1:convIt,1],xlabel="Iteration",ylabel="central lambda gap",xlims=(1,convIt),legend=false,yscale=:log10)
    # convItPlotalad2=plot(dCM.lamIt2Norm[1:convIt,1],xlabel="Iteration",ylabel="2-Norm Dual",xlims=(1,convIt),legend=false,yscale=:log10)
    # constPlotalad1=plot(dCM.coupl1Norm[1:convIt,1],xlabel="Iteration",ylabel="1-Norm coupl",xlims=(1,convIt),legend=false,yscale=:log10)
    # #savefig(lamPlotadmm,path*"lamPlotadmm"*".png")

    #save current state and update for next timeSteps
    dSol.Tpred[stepI,1]=dLogadmm.Tpred[1,convIt]
    dSol.Un[stepI,:]=dLogadmm.Un[1:N,convIt]
    dSol.Sn[stepI,:]=dLogadmm.Sn[1:N,convIt]
    dSol.Z[stepI,:]=dLogadmm.Z[1:S,convIt]
    dSol.uSum[stepI,1]=dLogadmm.uSum[1,convIt]
    dSol.zSum[stepI,1]=dLogadmm.zSum[1,convIt]
    dSol.Itotal[stepI,1]=dLogadmm.Itotal[1,convIt]
    dSol.Tactual[stepI,1]=dLogadmm.Tactual[1,convIt]
    dSol.convIt[stepI,1]=convIt
    dSol.timeSolve[stepI,1]=mean(dLogadmm.timeSolve[1,1:convIt])

    # new states
    global t0=dSol.Tactual[stepI,1]
    global s0=dSol.Sn[stepI,:]

    if convIt==1
        dSol.lamCoupl[stepI,1]=prevLam[1,1]
        if stepI+horzLen==evS.K
            newLam=prevLam[2:horzLen+1,1]
            newVu=prevVu[(N+1):(N*(horzLen+1)),1]
            newVz=prevVz[(S+1):(S*(horzLen+1)),1]
        else
            newLam=vcat(prevLam[2:horzLen+1,1],prevLam[horzLen+1,1])
            newVu=vcat(prevVu[(N+1):(N*(horzLen+1)),1],prevVu[((N*horzLen)+1):(N*(horzLen+1)),1])
            newVz=vcat(prevVz[(S+1):(S*(horzLen+1)),1],prevVz[((S*horzLen)+1):(S*(horzLen+1)),1])
        end
    else
        dSol.lamCoupl[stepI,1]=dLogadmm.Lam[1,convIt-1]
        if stepI+horzLen==evS.K
            newLam=dLogadmm.Lam[2:horzLen+1,convIt-1]
            newVu=dLogadmm.Vu[(N+1):(N*(horzLen+1)),convIt-1]
            newVz=dLogadmm.Vz[(S+1):(S*(horzLen+1)),convIt-1]
        else
            newLam=vcat(dLogadmm.Lam[2:horzLen+1,convIt-1],dLogadmm.Lam[horzLen+1,convIt-1])
            newVu=vcat(dLogadmm.Vu[(N+1):(N*(horzLen+1)),convIt-1],dLogadmm.Vu[((N*horzLen)+1):(N*(horzLen+1)),convIt-1])
            newVz=vcat(dLogadmm.Vz[(S+1):(S*(horzLen+1)),convIt-1],dLogadmm.Vz[((S*horzLen)+1):(S*(horzLen+1)),convIt-1])
        end
    end
    global prevLam=round.(newLam,sigdigits=roundSigFigs)
    global prevVu=round.(newVu,sigdigits=roundSigFigs)
    global prevVz=round.(newVz,sigdigits=roundSigFigs)
    #global ρADMMp=round.(ogρ,digits=2)

    return nothing
end

function pwlEVadmm(maxIt::Int,evS,cSave::centralLogStruct,slack::Bool,roundSigFigs::Int, silent::Bool)

    horzLen=evS.K1
    K=evS.K
    N=evS.N
    S=evS.S
    dSol=solutionStruct(K=K,N=N,S=S)
    dCM=convMetricsStruct(maxIt=maxIt,logLength=length(saveLogInd))

    Juno.progress() do id
        for stepI=1:K
            @info "$(Dates.format(Dates.now(),"HH:MM:SS")): $(stepI) of $(K)" progress=stepI/K _id=id
            @printf "%s: time step %g of %g...." Dates.format(Dates.now(),"HH:MM:SS") stepI K
            try
                dSol.timeT[stepI]=@elapsed runEVADMMStep(stepI,maxIt,evS,dSol,dCM,cSave,roundSigFigs,silent)
                @printf "convIt: %g\n" dSol.convIt[stepI,1]
            catch e
                @printf "error: %s" e
                break
            end
        end
    end

    objFun(sn,u)=sum(sum((sn[k,n]-1)^2*evS.Qsi[n,k] for n=1:N) +sum((u[k,n])^2*evS.Ri[n,k] for n=1:N) for k=1:evS.K)
    dSol.objVal[1,1]=objFun(dSol.Sn,dSol.Un)

    return dSol, dCM
end

#ALADIN
function localEVALAD(evInd::Int,p::Int,stepI::Int,σU::Array{Float64,2},σS::Array{Float64,2},evS,dLogalad::itLogPWL,
    ind,evVu,evVs,itLam,s0,itρ,slack,solverSilent,roundSigFigs,silent)
    # local decoupled ALADIN EV problem

    N=evS.N
    S=evS.S
    horzLen=min(evS.K1,evS.K-stepI)

    tolU=1e-6
    tolS=1e-6

    target=zeros((horzLen+1),1)
    target[max(1,(evS.Kn[evInd,1]-(stepI-1))):1:length(target),1].=evS.Snmin[evInd,1]
    evM = Model(solver = GurobiSolver(NumericFocus=3))
    @variable(evM,sn[1:(horzLen+1)])
    @variable(evM,u[1:(horzLen+1)])
    if slack @variable(evM,slackSn) end
    objExp=sum((sn[k,1]-1)^2*evS.Qsi[evInd,stepI+k-1]+(u[k,1])^2*evS.Ri[evInd,stepI+k-1]+
                            itLam[k,1]*(u[k,1])+
                            itρ/2*(u[k,1]-evVu[k,1])*σU[evInd,1]*(u[k,1]-evVu[k,1])+
                            itρ/2*(sn[k,1]-evVs[k,1])*σS[evInd,1]*(sn[k,1]-evVs[k,1]) for k=1:horzLen+1)
    if slack
        append!(objExp,sum(evS.β[n]*slackSn^2 for n=1:N))
    end
    # if reg
    #     append!(objExp,reg_weight*sum(sn[k,1]^2+u[k,1]^2 for k=1:horzLen+1))
    # end
    @objective(evM,Min,objExp)
    @constraint(evM,sn[1,1]==s0[evInd,1]+evS.ηP[evInd,1]*u[1,1])
    @constraint(evM,[k=1:horzLen],sn[k+1,1]==sn[k,1]+evS.ηP[evInd,1]*u[k+1,1])
    @constraint(evM,socKappaMax,sn.<=1)
    if slack
        @constraint(evM,sn.>=target.*(1-slackSn))
        @constraint(evM,slackSn>=0)
    else
        @constraint(evM,socKappaMin,sn.>=target)
    end
    @constraint(evM,curKappaMax,u.<=evS.imax[evInd,1])
    @constraint(evM,curKappaMin,u.>=evS.imin[evInd,1])

    if solverSilent
        @suppress_out begin
            statusEVM = solve(evM)
        end
    else
        statusEVM = solve(evM)
    end

    @assert statusEVM==:Optimal "ALAD EV NLP optimization not solved to optimality"

    uVal=getvalue(u)
    snVal=getvalue(sn)

    cValMax=abs.(uVal.-evS.imax[evInd,1]).<tolU
    cValMin=abs.(uVal.-evS.imin[evInd,1]).<tolU
    dLogalad.Cuu[ind,p]=1cValMax
    dLogalad.Cul[ind,p]=-1cValMin

    cValMax=abs.(snVal.-1).<tolS
    if slack
        cValMin=abs.(snVal.-target*(1-getvalue(slackSn))).<tolS
    else
        cValMin=abs.(snVal.-target).<tolS
    end
    dLogalad.Csu[ind,p]=1cValMax
    dLogalad.Csl[ind,p]=-1cValMin

    dLogalad.slackSn[evInd]= if slack getvalue(slackSn) else 0 end
    dLogalad.Sn[ind,p]=round.(snVal,sigdigits=roundSigFigs)
    dLogalad.Un[ind,p]=round.(uVal,sigdigits=roundSigFigs)

    # if reg
    #     dLogalad.Gu[ind,p]=round.(2*evS.Ri[evInd,stepI:stepI+horzLen]*uVal+2*reg_weight*uVal,sigdigits=roundSigFigs)
    #     dLogalad.Gs[ind,p]=round.(2*evS.Qsi[evInd,stepI:stepI+horzLen]*snVal.-2*evS.Qsi[evInd,stepI:stepI+horzLen]+2*reg_weight*snVal,sigdigits=roundSigFigs)
    # else
    dLogalad.Gu[ind,p]=round.(2*evS.Ri[evInd,stepI:stepI+horzLen].*uVal,sigdigits=roundSigFigs)
    dLogalad.Gs[ind,p]=round.(2*evS.Qsi[evInd,stepI:stepI+horzLen].*snVal .- 2*evS.Qsi[evInd,stepI:stepI+horzLen],sigdigits=roundSigFigs)
    #end
    dLogalad.timeSolve[1,p]=max(getsolvetime(evM),dLogalad.timeSolve[1,p])

    #use convex ALADIN approach
    #dLogalad.Gu[ind,p]=σU[evInd,1]*(evVu-uVal)-prevLam[:,1]
    #dLogalad.Gs[ind,p]=σS[evInd,1]*(evVs-snVal)#-prevLam[:,1]
    return nothing
end

function localXFRMALAD(p::Int,stepI::Int,σZ::Float64,σT::Float64,evS,dLogalad::itLogPWL,
    itLam,itVz,itVt,itρ,roundSigFigs)
    # local transformer ALADIN problem

    N=evS.N
    S=evS.S
    deltaI=evS.deltaI
    horzLen=min(evS.K1,evS.K-stepI)

    tolT=1e-3
    tolZ=1e-4

    #N+1 decoupled problem aka transformer current
    tM = Model(solver = GurobiSolver(NumericFocus=3))
    @variable(tM,z[1:(S)*(horzLen+1)])
    @variable(tM,t[1:(horzLen+1)])
    objExp=sum(-itLam[k,1]*sum(z[(k-1)*(S)+s] for s=1:S)+
              itρ/2*σZ*(sum(z[(k-1)*(S)+s] for s=1:S)-sum(itVz[(k-1)*(S)+s,1] for s=1:S))^2+
              itρ/2*σT*(t[k]-itVt[k,1])^2  for k=1:(horzLen+1))
      # if reg
      #     append!(objExp,reg_weight*sum(sum(z[k,1]^2 for s=1:S)+t[k,1]^2 for k=1:horzLen+1))
      # end
    @objective(tM,Min, objExp)
    @constraint(tM,tempCon1,t[1]-evS.τP*t0-evS.γP*deltaI*sum((2*s-1)*z[s] for s=1:S)-evS.ρP*evS.Tamb[stepI,1]==0)
    @constraint(tM,tempCon2[k=1:horzLen],t[k+1]-evS.τP*t[k]-evS.γP*deltaI*sum((2*s-1)*z[(k)*(S)+s] for s=1:S)-evS.ρP*evS.Tamb[stepI+k,1]==0)
    if noTlimit==false
        @constraint(tM,upperTCon,t.<=evS.Tmax)
    end
    @constraint(tM,lowerTCon,t.>=0)
    @constraint(tM,pwlKappaMin,z.>=0)
    @constraint(tM,pwlKappaMax,z.<=deltaI)
    if solverSilent
        @suppress_out begin
            statusTM = solve(tM)
        end
    else
        statusTM = solve(tM)
    end
    @assert statusTM==:Optimal "ALAD XFRM NLP optimization not solved to optimality"

    zVal=getvalue(z)
    tVal=getvalue(t)

    cValMax=abs.(zVal.-deltaI).<tolZ
    cValMin=abs.(zVal.-0).<tolZ
    dLogalad.Czu[:,p]=1cValMax
    dLogalad.Czl[:,p]=-1cValMin

    cValMax=abs.(tVal.-evS.Tmax).<tolT
    cValMin=abs.(tVal.-0).<tolT
    dLogalad.Ctu[:,p]=1cValMax
    dLogalad.Ctl[:,p]=-1cValMin

    dLogalad.Tpred[:,p]=round.(tVal,sigdigits=roundSigFigs)
    dLogalad.Z[:,p]=round.(zVal,sigdigits=roundSigFigs)

    # if reg
    #     dLogalad.Gz[:,p]=round.(2*reg_weight*zVal,sigdigits=roundSigFigs)
    #     dLogalad.Gt[:,p]=round.(2*reg_weight*tVal,sigdigits=roundSigFigs)
    # else
        dLogalad.Gz[:,p].=0
        dLogalad.Gt[:,p].=0
    #end
    dLogalad.timeSolve[1,p]=max(getsolvetime(tM),dLogalad.timeSolve[1,p])

    #dLogalad.Gz[:,p]=σZ*(prevVz-zVal)-repeat(-prevLam[:,1],inner=S)
    #dLogalad.Gt[:,p]=σT*(prevVt-xtVal)
    return nothing
end

function coordALAD(p::Int,stepI::Int,itμ,evS,itLam,itVu,itVs,itVz,itVt,itρ,
    dLogalad::itLogPWL,roundSigFigs)
    # ALADIN coordinator problem

    N=evS.N
    S=evS.S
    deltaI=evS.deltaI
    horzLen=min(evS.K1,evS.K-stepI)

    # if reg
    #     Hu=2*evS.Ri.+2*reg_weight
    #     Hs=2*evS.Qsi.+2*reg_weight
    #     Hz=2*reg_weight
    #     Ht=2*reg_weight
    # else
        Hu=2*evS.Ri[:,stepI:stepI+horzLen]
        Hs=2*evS.Qsi[:,stepI:stepI+horzLen]
        # Hu=2*evS.Ri *((1.5-2.5)*rand()+2.5)
        # Hs=2*evS.Qsi *((1.5-2.5)*rand()+2.5)
        Hz=0
        Ht=0
        # Hz=1e-6
        # Ht=1e-6
    #end

    #μALADp=μALADp*p^2

    #coupled QP
    cM = Model(solver = GurobiSolver(NumericFocus=3))
    #cM = Model(solver = GurobiSolver(Presolve=1))
    #cM = Model(solver = IpoptSolver())
    #cM = Model(solver = MosekSolver())

    @variable(cM,dUn[1:(N)*(horzLen+1)])
    @variable(cM,dSn[1:(N)*(horzLen+1)])
    @variable(cM,dZ[1:(S)*(horzLen+1)])
    @variable(cM,dT[1:(horzLen+1)])
    @variable(cM,relaxS[1:(horzLen+1)])
    objExp=sum(sum(0.5*dUn[(k-1)*N+n,1]^2*Hu[n,k]+dLogalad.Gu[(k-1)*N+n,p]*dUn[(k-1)*N+n,1] +
                   0.5*dSn[(k-1)*N+n,1]^2*Hs[n,k]+dLogalad.Gs[(k-1)*N+n,p]*dSn[(k-1)*N+n,1] for n=1:N) +
               sum(0.5*dZ[(k-1)*(S)+s,1]^2*Hz for s=1:S)+
               0.5*dT[k,1]^2*Ht   for k=1:(horzLen+1))
    objExp=objExp+itLam[:,1]'*relaxS+itμ/2*sum(relaxS[k,1]^2 for k=1:horzLen+1)
    #objExp=objExp+itμ/2*sum(relaxS[k,1]^2 for k=1:horzLen+1)
    objExp=objExp+dot(dLogalad.Gz[:,p],dZ)+dot(dLogalad.Gt[:,p],dT)

    @objective(cM,Min,objExp)
    Unp=round.(dLogalad.Un[:,p],sigdigits=roundSigFigs)
    Zp=round.(dLogalad.Z[:,p],sigdigits=roundSigFigs)
    @constraint(cM,currCon[k=1:horzLen+1],sum(Unp[(k-1)*(N)+n,1]+dUn[(k-1)*(N)+n,1] for n=1:N)-
                                          sum(Zp[(k-1)*(S)+s,1]+dZ[(k-1)*(S)+s,1] for s=1:S)==
                                          -evS.iD_pred[stepI+(k-1)]+relaxS[k,1])
    #local equality constraints C*(X+deltaX)=0 is same as C*deltaX=0 since we already know CX=0
    @constraint(cM,stateCon1[n=1:N],dSn[n,1]==evS.ηP[n,1]*dUn[n,1])
    @constraint(cM,stateCon2[k=1:horzLen,n=1:N],dSn[n+(k)*(N),1]==dSn[n+(k-1)*(N),1]+evS.ηP[n,1]*dUn[n+(k)*(N),1])
    @constraint(cM,tempCon1,dT[1,1]==evS.γP*deltaI*sum((2*s-1)*dZ[s,1] for s=1:S))
    @constraint(cM,tempCon2[k=1:horzLen],dT[k+1,1]==evS.τP*dT[k,1]+evS.γP*deltaI*sum((2*s-1)*dZ[k*S+s,1] for s=1:S))

    #active local constraints
    if eqForm
        @constraint(cM,dLogalad.Czu[:,p].*dZ.==0)
        @constraint(cM,dLogalad.Cuu[:,p].*dUn.==0)
        @constraint(cM,dLogalad.Csu[:,p].*dSn.==0)
        @constraint(cM,dLogalad.Ctu[:,p].*dT.==0)
        @constraint(cM,dLogalad.Czl[:,p].*dZ.==0)
        @constraint(cM,dLogalad.Cul[:,p].*dUn.==0)
        @constraint(cM,dLogalad.Csl[:,p].*dSn.==0)
        @constraint(cM,dLogalad.Ctl[:,p].*dT.==0)
    else
        @constraint(cM,dLogalad.Czu[:,p].*dZ.<=0)
        @constraint(cM,dLogalad.Cuu[:,p].*dUn.<=0)
        @constraint(cM,dLogalad.Csu[:,p].*dSn.<=0)
        @constraint(cM,dLogalad.Ctu[:,p].*dT.<=0)
        @constraint(cM,dLogalad.Czl[:,p].*dZ.<=0)
        @constraint(cM,dLogalad.Cul[:,p].*dUn.<=0)
        @constraint(cM,dLogalad.Csl[:,p].*dSn.<=0)
        @constraint(cM,dLogalad.Ctl[:,p].*dT.<=0)
    end

    if solverSilent
        @suppress_out begin
            statusM = solve(cM)
        end
    else
        statusM = solve(cM)
    end
    @assert statusM==:Optimal "ALAD Central QP optimization not solved to optimality"

    #update step
    α1=1
    α2=1
    α3=1
    #α3=1/ceil(p/2)

    # dLogalad.Lam[:,p]=round.(itLam[:,1]+α3*(-getdual(currCon)-itLam[:,1]),sigdigits=roundSigFigs)
    dLogalad.Lam[:,p]=max.(round.(itLam[:,1]+α3*(-getdual(currCon)-itLam[:,1]),sigdigits=roundSigFigs),0)
    dLogalad.Vu[:,p]=round.(itVu[:,1]+α1*(dLogalad.Un[:,p]-itVu[:,1])+α2*getvalue(dUn),sigdigits=roundSigFigs)
    dLogalad.Vz[:,p]=round.(itVz[:,1]+α1*(dLogalad.Z[:,p]-itVz[:,1])+α2*getvalue(dZ),sigdigits=roundSigFigs)
    dLogalad.Vs[:,p]=round.(itVs[:,1]+α1*(dLogalad.Sn[:,p]-itVs[:,1])+α2*getvalue(dSn),sigdigits=roundSigFigs)
    dLogalad.Vt[:,p]=round.(itVt[:,1]+α1*(dLogalad.Tpred[:,p]-itVt[:,1])+α2*getvalue(dT),sigdigits=roundSigFigs)
    dLogalad.timeSolve[1,p]+=getsolvetime(cM)

    # lamPlot=plot(dLogalad.Lam[:,p],label="ineq")
    # vtPlot=plot(dLogalad.Vt[:,p],label="ineq")
    # plot!(lamPlot,dLogalad.Lam[:,p],label="eq")
    # plot!(vtPlot,dLogalad.Vt[:,p],label="eq")


    #dCM.lamIt[p,1]=norm(dLogalad.Lam[:,p]-itLam[:,1],2)
    #dCM.lam[p,1]=norm(dLogalad.Lam[:,p]-cSave.lamCoupl[stepI:(horzLen+stepI)],2)
    if !silent
        #@printf "lastGap    %e after %g iterations\n" dCM.lamIt[p,1] p
        #@printf "convLamGap %e after %g iterations\n\n" dCM.lam[p,1] p
    end

    dLogalad.itUpdate[1,p]=min(itρ*ρRate,ρALADmax) #increase ρ every iteration
    dLogalad.itμUpdate[1,p]=min(itμ*μRate,μALADmax) #increase μ every iteration
    #ΔY[1,p]=norm(vcat(getvalue(dUn),getvalue(dZ),getvalue(dSn),getvalue(dT)),Inf)
    return nothing
end

function runEVALADIt(p,stepI,evS,itLam,itVu,itVz,itVs,itVt,itρ,itμ,dLogalad,dCM,dSol,cSave,eqForm,roundSigFigs,silent)
    # Function for a single iteration of the ALADIN algorithm

    K=evS.K
    N=evS.N
    S=evS.S
    horzLen=min(evS.K1,K-stepI)

    #initialize with current states
    global s0
    global t0

    #ALADIN tuning
    # if eqForm
    #     #println("Running Eq ALADIN")
    #     scalingF=1e-4
    # else
    #     #println("Running ineq ALADIN")
    #     scalingF=1e-4
    # end
    # σZ=scalingF*1e1
    # σT=scalingF
    # σU=ones(N,1)
    # σS=ones(N,1)/10#for kA

    σZ=1.0/2.5
    σT=1/200
    σU=ones(N,1)/.05
    σS=ones(N,1)


    #μALADp=5e3
    #μALADp=1e8
    # μRate=1
    # μALADmax=2e9

    #solve decoupled local EV problems
    if runParallel
        @sync @distributed for evInd=1:N
            ind=[evInd]
            for k=1:horzLen
                append!(ind,k*N+evInd)
            end
            evVu=itVu[ind,1]
            evVs=itVs[ind,1]
            localEVALAD(evInd,p,stepI,σU,σS,evS,dLogalad,ind,evVu,evVs,itLam,s0,itρ,slack,solverSilent,roundSigFigs,silent)
        end
    else
        for evInd=1:N
            ind=[evInd]
            for k=1:horzLen
                append!(ind,k*N+evInd)
            end
            evVu=itVu[ind,1]
            evVs=itVs[ind,1]
            localEVALAD(evInd,p,stepI,σU,σS,evS,dLogalad,ind,evVu,evVs,itLam,s0,itρ,slack,solverSilent,roundSigFigs,silent)
        end
    end

    # solve local transformer problem
    localXFRMALAD(p,stepI,σZ,σT,evS,dLogalad,itLam,itVz,itVt,itρ,roundSigFigs)

    # save useful vectors
    for k=1:horzLen+1
        dLogalad.uSum[k,p]=sum(dLogalad.Un[(k-1)*N+n,p] for n=1:N)
        dLogalad.zSum[k,p]=sum(dLogalad.Z[(k-1)*(S)+s,p] for s=1:S)
        dLogalad.couplConst[k,p]=dLogalad.uSum[k,p] + evS.iD_pred[stepI+(k-1),1] - dLogalad.zSum[k,p]
    end

    # calculate actual temperature from nonlinear model of XFRM
    for k=1:horzLen+1
        dLogalad.Itotal[k,p]=dLogalad.uSum[k,p] + evS.iD_actual[stepI+(k-1),1]
    end
    dLogalad.Tactual[1,p]=evS.τP*t0+evS.γP*dLogalad.Itotal[1,p]^2+evS.ρP*evS.Tamb[stepI,1] #fix for mpc
    for k=1:horzLen
        dLogalad.Tactual[k+1,p]=evS.τP*dLogalad.Tactual[k,p]+evS.γP*dLogalad.Itotal[k,p]^2+evS.ρP*evS.Tamb[stepI+k,1]  #fix for mpc
    end

    # ALADIN coordination problem
    coordALAD(p,stepI,itμ,evS,itLam,itVu,itVs,itVz,itVt,itρ,dLogalad,roundSigFigs)

    #check for convergence
    constGap=norm(dLogalad.couplConst[:,p],1)
    itGap = norm(dLogalad.Lam[:,p]-itLam[:,1],2)
    objFun(sn,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,stepI+k-1] for n=1:N) +
                     sum((u[(k-1)*N+n,1])^2*evS.Ri[n,stepI+k-1]       for n=1:N) for k=1:horzLen+1)
    dLogalad.objVal[1,p]=objFun(dLogalad.Sn[:,p],dLogalad.Un[:,p])
    auxGap=norm(vcat(σU[1]*(itVu[:,1]-dLogalad.Un[:,p]),σZ*(itVz[:,1]-dLogalad.Z[:,p]),
                 σT*(itVt[:,1]-dLogalad.Tpred[:,p]),σS[1]*(itVs[:,1]-dLogalad.Sn[:,p])),1)

    #only if saving
    if stepI in saveLogInd
        ind=findall(x->x==stepI,saveLogInd)[1]
        dCM.coupl1Norm[p,ind]=constGap
        dCM.coupl2Norm[p,ind]=norm(dLogalad.couplConst[:,p],2)
        dCM.couplInfNorm[p,ind]=norm(dLogalad.couplConst[:,p],Inf)
        dCM.lamIt1Norm[p,ind]= norm(dLogalad.Lam[:,p]-itLam[:,1],1)
        dCM.lamIt2Norm[p,ind]=itGap
        dCM.lamItInfNorm[p,ind]= norm(dLogalad.Lam[:,p]-itLam[:,1],Inf)
        dCM.objAbs[p,ind]=abs(dLogalad.objVal[1,p]-cSave.Obj[1,1,ind])
        dCM.objPerc[p,ind]=abs(dLogalad.objVal[1,p]-cSave.Obj[1,1,ind])/cSave.Obj[1,1,ind]
        prevObj = if p==1 dLogalad.objVal[1,p] else dLogalad.objVal[1,p-1] end
        dCM.objItPerc[p,ind]=abs(dLogalad.objVal[1,p]-prevObj)/prevObj
        dCM.lam1Norm[p,ind]= norm(dLogalad.Lam[:,p]-cSave.Lam[1:(horzLen+1),:,ind],1)
        dCM.lam2Norm[p,ind]= norm(dLogalad.Lam[:,p]-cSave.Lam[1:(horzLen+1),:,ind],2)
        dCM.lamInfNorm[p,ind]= norm(dLogalad.Lam[:,p]-cSave.Lam[1:(horzLen+1),:,ind],Inf)
        dCM.t1Norm[p,ind]= norm(dLogalad.Tactual[:,p]-cSave.Tactual[1:(horzLen+1),:,ind],1)
        dCM.t2Norm[p,ind]= norm(dLogalad.Tactual[:,p]-cSave.Tactual[1:(horzLen+1),:,ind],2)
        dCM.tInfNorm[p,ind]= norm(dLogalad.Tactual[:,p]-cSave.Tactual[1:(horzLen+1),:,ind],Inf)
        zReshape=zeros(horzLen+1,S)
        uReshape=zeros(horzLen+1,N)
        for ii= 1:N
            uReshape[:,ii]=dLogalad.Un[collect(ii:N:length(dLogalad.Un[:,p])),p]
        end
        for ii= 1:S
            zReshape[:,ii]=dLogalad.Z[collect(ii:S:length(dLogalad.Z[:,p])),p]
        end
        dCM.un1Norm[p,ind]= norm(uReshape-cSave.Un[1:(horzLen+1),:,ind],1)
        dCM.un2Norm[p,ind]= norm(uReshape-cSave.Un[1:(horzLen+1),:,ind],2)
        dCM.unInfNorm[p,ind]= norm(uReshape-cSave.Un[1:(horzLen+1),:,ind],Inf)
        dCM.z1Norm[p,ind]= norm(zReshape-cSave.Z[1:(horzLen+1),:,ind],1)
        dCM.z2Norm[p,ind]= norm(zReshape-cSave.Z[1:(horzLen+1),:,ind],2)
        dCM.zInfNorm[p,ind]= norm(zReshape-cSave.Z[1:(horzLen+1),:,ind],Inf)
    end

    # if  constGap<=primChk && itGap<=dualChk
    if  constGap<=primChk && auxGap<=auxChk
        if !silent @printf "Converged after %g iterations" p end
        @printf "Y"
        convIt=p
        return true
    else
        if !silent
            #@printf "objItPerc  %e after %g iterations\n" itGap p
            # @printf "lamIt      %e after %g iterations\n" itGap p
            @printf "auxGap     %e after %g iterations\n" auxGap p
            @printf "constGap   %e after %g iterations\n\n" constGap p
        end
    end
    return false
end

function runEVALADStep(stepI,maxIt,evS,dSol,dCM,cSave,eqForm,roundSigFigs,silent)
    # Function for each MPC time step

    # # working for eqForm for simple example
    # # μALADp = 20
    # # μRate=1.01
    # # μALADmax=1e2
    #
    # # working for ineq form
    # eqForm=false
    # ρALADp = 1
    # ρRate=1.1
    # ρALADmax=1e4
    # μALADp = 1e3
    # μRate=1.15
    # μALADmax=1e4
    #
    # #slow convergence for equality
    # eqForm=true
    # ρALADp = .1
    # ρRate=1.1
    # ρALADmax=1e1
    # μALADp = .1
    # μRate=1.1
    # μALADmax=1e1
    #
    # # also ok
    # eqForm=true
    # ρALADp = 1
    # ρRate=1.1
    # ρALADmax=2e3
    # μALADp = 15
    # μRate=1.01
    # μALADmax=1e3
    #
    # eqForm=true
    # ρALADp = .001
    # ρRate=1.1
    # ρALADmax=2e3
    # μALADp = 1
    # μRate=1.1
    # μALADmax=1e2

    K=evS.K
    N=evS.N
    S=evS.S
    horzLen=min(evS.K1,K-stepI)
    dLogalad=itLogPWL(horzLen=horzLen,N=N,S=S)
    p=1
    timeStart=now()
    # keep iterating while not converged, under maximum number of iterations, and time left in timestep
    while (p<=maxIt && round(now()-timeStart,Second)<=Dates.Second(9/10*evS.Ts))
        # while (p<=maxIt)
        # global p
        if p==1 #first iteration: hot start with previous timestep solutions
            itLam=prevLam
            itVu=prevVu
            itVz=prevVz
            itVs=prevVs
            itVt=prevVt
            itρ=ρALADp
            itμ=μALADp
        else   # use last iteration solutions
            itLam=round.(dLogalad.Lam[:,(p-1)],sigdigits=roundSigFigs)
            itVu=round.(dLogalad.Vu[:,(p-1)],sigdigits=roundSigFigs)
            itVz=round.(dLogalad.Vz[:,(p-1)],sigdigits=roundSigFigs)
            itVs=round.(dLogalad.Vs[:,(p-1)],sigdigits=roundSigFigs)
            itVt=round.(dLogalad.Vt[:,(p-1)],sigdigits=roundSigFigs)
            itρ=round.(dLogalad.itUpdate[1,(p-1)],sigdigits=roundSigFigs)
            # itμ=round(min(itμ*μRate,μALADmax),sigdigits=4)
            itμ=round(dLogalad.itμUpdate[1,(p-1)],sigdigits=4)
        end
        #println(itμ)
        cFlag=runEVALADIt(p,stepI,evS,itLam,itVu,itVz,itVs,itVt,itρ,itμ,dLogalad,dCM,dSol,cSave,eqForm,roundSigFigs,silent)
        global convIt=p
        if cFlag
            break # if converged, end iterations
        end
        p+=1
    end
    if stepI in saveLogInd
        ind=findall(x->x==stepI,saveLogInd)[1]
        dCM.convIt[1,1,ind]=convIt
    end

    # # # #these are plots I use to debug a single time step
    # print(round(now()-timeStart,Second))
    # xPlot=zeros(horzLen+1,N)
    # uPlot=zeros(horzLen+1,N)
    # for ii= 1:N
    # 	xPlot[:,ii]=dLogalad.Sn[collect(ii:N:length(dLogalad.Sn[:,convIt])),convIt]
    #     uPlot[:,ii]=dLogalad.Un[collect(ii:N:length(dLogalad.Un[:,convIt])),convIt]
    # end
    #
    # p1=plot(dLogalad.uSum[:,convIt],xlabel="Time",ylabel="Current Sum (kA)",xlims=(0,horzLen+1),label="ALADIN Open Loop")
    # plot!(p1,sum(cSave.Un[:,:,ind],dims=2),seriescolor=:black,linewidth=2,linealpha=0.8,label="Central Open Loop")
    # plot(xPlot)
    # pd3alad=plot(hcat(dLogalad.Tactual[:,convIt],dLogalad.Tpred[:,convIt]),label=["Actual Temp" "PWL Temp"],xlims=(0,evS.K),xlabel="Time",ylabel="Temp (K)")
    # plot!(pd3alad,1:horzLen+1,evS.Tmax*ones(evS.K),label="XFRM Limit",line=(:dash,:red))
    #
    # #convergence plots
    # halfCI=Int(floor(convIt/2))
    # if halfCI>0
    #     CList=reshape([range(colorant"blue", stop=colorant"yellow",length=halfCI);
    #                    range(colorant"yellow", stop=colorant"red",length=convIt-halfCI)], 1, convIt);
    # else
    #     CList=colorant"red";
    # end
    # uSumPlotalad=plot(dLogalad.uSum[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Current Sum (kA)",xlims=(0,horzLen+1),legend=false)
    # plot!(uSumPlotalad,sum(cSave.Un[:,:,ind],dims=2),seriescolor=:black,linewidth=2,linealpha=0.8)
    #
    # zSumPlotalad=plot(dLogalad.zSum[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Z sum",xlims=(0,horzLen+1),legend=false)
    # plot!(zSumPlotalad,sum(cSave.Z[:,:,ind],dims=2),seriescolor=:black,linewidth=2,linealpha=0.8)
    #
    # constPlotalad2=plot(dLogalad.couplConst[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="curr constraint diff",xlims=(0,horzLen+1),legend=false)
    #
    # lamPlotalad=plot(dLogalad.Lam[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Lambda",legend=false)
    # plot!(lamPlotalad,cSave.Lam[:,:,ind],seriescolor=:black,linewidth=2,linealpha=0.8)
    #
    # activeSet=zeros(convIt,1)
    # setChanges=zeros(convIt,1)
    # for ii=2:convIt
    #     activeSet[ii,1]=sum(abs.(dLogalad.Csu[:,ii]))+sum(abs.(dLogalad.Ctu[:,ii]))+
    #               		sum(abs.(dLogalad.Cuu[:,ii]))+sum(abs.(dLogalad.Czu[:,ii]))+
    # 					sum(abs.(dLogalad.Csl[:,ii]))+sum(abs.(dLogalad.Ctl[:,ii]))+
    # 				    sum(abs.(dLogalad.Cul[:,ii]))+sum(abs.(dLogalad.Czl[:,ii]))
    #     setChanges[ii,1]=sum(abs.(dLogalad.Csu[:,ii]-dLogalad.Csu[:,ii-1]))+sum(abs.(dLogalad.Ctu[:,ii]-dLogalad.Ctu[:,ii-1]))+
    #                      sum(abs.(dLogalad.Cuu[:,ii]-dLogalad.Cuu[:,ii-1]))+sum(abs.(dLogalad.Czu[:,ii]-dLogalad.Czu[:,ii-1]))+
    # 					 sum(abs.(dLogalad.Csl[:,ii]-dLogalad.Csl[:,ii-1]))+sum(abs.(dLogalad.Ctl[:,ii]-dLogalad.Ctl[:,ii-1]))+
    # 				     sum(abs.(dLogalad.Cul[:,ii]-dLogalad.Cul[:,ii-1]))+sum(abs.(dLogalad.Czl[:,ii]-dLogalad.Czl[:,ii-1]))
    # end
    #
    # activeSetPlot=plot(2:convIt,activeSet[2:convIt],xlabel="Iteration",ylabel="Total Active inequality constraints",
    #                    legend=false,xlims=(2,convIt))
    # setChangesPlot=plot(10:convIt,setChanges[10:convIt],xlabel="Iteration",ylabel="Changes in Active inequality constraints",
    #                   legend=false,xlims=(2,convIt))
    # #solChangesplot=plot(2:convIt,hcat(ΔY[2:convIt],convCheck[2:convIt]),xlabel="Iteration",labels=["ΔY" "y-x"],xlims=(1,convIt))
    #
    # fAbsPlot=plot(dCM.objAbs[1:convIt,ind],xlabel="Iteration",ylabel="obj function gap",legend=false,yscale=:log10)
    # fPlot=plot(dCM.objPerc[2:convIt-1,ind],xlabel="Iteration",ylabel="obj function gap",legend=false,yscale=:log10)
    # fItPlot=plot(dCM.objItPerc[2:convIt-1,ind],xlabel="Iteration",ylabel="obj function gap",legend=false,yscale=:log10)
    # convPlotalad=plot(dCM.lam2Norm[1:convIt,1],xlabel="Iteration",ylabel="central lambda gap",xlims=(1,convIt),legend=false,yscale=:log10)
    # convItPlotalad1=plot(dCM.lamIt1Norm[1:convIt,1],xlabel="Iteration",ylabel="1-Norm Dual",xlims=(1,convIt),legend=false,yscale=:log10)
    # convItPlotalad2=plot(dCM.lamIt2Norm[1:convIt,1],xlabel="Iteration",ylabel="2-Norm Dual",xlims=(1,convIt),legend=false,yscale=:log10)
    # convItPlotaladInf=plot(dCM.lamItInfNorm[1:convIt,1],xlabel="Iteration",ylabel="Inf-Norm Dual",xlims=(1,convIt),legend=false,yscale=:log10)
    # constPlotalad1=plot(dCM.coupl1Norm[1:convIt,1],xlabel="Iteration",ylabel="1-Norm coupl",xlims=(1,convIt),legend=false,yscale=:log10)
    # constPlotalad2=plot(dCM.coupl2Norm[1:convIt,1],xlabel="Iteration",ylabel="2-Norm coupl",xlims=(1,convIt),legend=false,yscale=:log10)
    # constPlotaladInf=plot(dCM.couplInfNorm[1:convIt,1],xlabel="Iteration",ylabel="Inf-Norm coupl",xlims=(1,convIt),legend=false,yscale=:log10)
    # #savefig(uSumPlotalad,path*"uSumPlotalad"*Dates.format(Dates.now(),"_HHMMSS_") *".png")
    #
    # # checkPlot2=plot(convItPlotalad1,constPlotalad1,convItPlotalad2,
    # #                 constPlotalad2,convItPlotaladInf,constPlotaladInf,layout=(3,2))
    # # pubPlot(checkPlot2,thickscale=1,sizeWH=(600,400),dpi=60)
    # # savefig(checkPlot2,path*"checkPlot2"*Dates.format(Dates.now(),"_HHMMSS_") *".png")
    # #
    # # checkPlot=plot(convItPlotalad,constPlotalad,fPlotalad,convPlotalad,layout=(2,2))
    # # pubPlot(checkPlot,thickscale=1,sizeWH=(600,400),dpi=60)
    # # savefig(checkPlot,path*"checkPlot"*Dates.format(Dates.now(),"_HHMMSS_") *".png")

    #save current state and update for next timeSteps
    dSol.Tpred[stepI,1]=dLogalad.Tpred[1,convIt]
    dSol.Un[stepI,:]=dLogalad.Un[1:N,convIt]
    dSol.Sn[stepI,:]=dLogalad.Sn[1:N,convIt]
    dSol.Z[stepI,:]=dLogalad.Z[1:S,convIt]
    dSol.uSum[stepI,1]=dLogalad.uSum[1,convIt]
    dSol.zSum[stepI,1]=dLogalad.zSum[1,convIt]
    dSol.Itotal[stepI,1]=dLogalad.Itotal[1,convIt]
    dSol.Tactual[stepI,1]=dLogalad.Tactual[1,convIt]
    dSol.convIt[stepI,1]=convIt
    dSol.timeSolve[stepI,1]=mean(dLogalad.timeSolve[1,1:convIt])

    # new states
    global t0=round(dSol.Tactual[stepI,1],sigdigits=roundSigFigs)
    global s0=dSol.Sn[stepI,:]

    # prepare solutions for next time step
    #TODO: clean this up
    if convIt==1
        dSol.lamCoupl[stepI,1]=prevLam[1,1]
        if stepI+horzLen==evS.K
            newLam=prevLam[2:horzLen+1,1]
            newVu=prevVu[(N+1):(N*(horzLen+1)),1]
            newVz=prevVz[(S+1):(S*(horzLen+1)),1]
            newVt=prevVt[2:horzLen+1,1]
            newVs=prevVs[(N+1):(N*(horzLen+1)),1]
        else
            newLam=vcat(prevLam[2:horzLen+1,1],prevLam[horzLen+1,1])
            newVu=vcat(prevVu[(N+1):(N*(horzLen+1)),1],prevVu[((N*horzLen)+1):(N*(horzLen+1)),1])
            newVz=vcat(prevVz[(S+1):(S*(horzLen+1)),1],prevVz[((S*horzLen)+1):(S*(horzLen+1)),1])
            newVt=vcat(prevVt[2:horzLen+1,1],prevVt[horzLen+1,1])
            newVs=vcat(prevVs[(N+1):(N*(horzLen+1)),1],prevVs[((N*horzLen)+1):(N*(horzLen+1)),1])
        end
    else
        dSol.lamCoupl[stepI,1]=dLogalad.Lam[1,convIt]
        if stepI+horzLen==evS.K
            newLam=dLogalad.Lam[2:horzLen+1,convIt]
            newVu=dLogalad.Vu[(N+1):(N*(horzLen+1)),convIt]
            newVz=dLogalad.Vz[(S+1):(S*(horzLen+1)),convIt]
            newVt=dLogalad.Vt[2:horzLen+1,convIt]
            newVs=dLogalad.Vs[(N+1):(N*(horzLen+1)),convIt]
        else
            newLam=vcat(dLogalad.Lam[2:horzLen+1,convIt-1],dLogalad.Lam[horzLen+1,convIt])
            newVu=vcat(dLogalad.Vu[(N+1):(N*(horzLen+1)),convIt-1],dLogalad.Vu[((N*horzLen)+1):(N*(horzLen+1)),convIt])
            newVz=vcat(dLogalad.Vz[(S+1):(S*(horzLen+1)),convIt-1],dLogalad.Vz[((S*horzLen)+1):(S*(horzLen+1)),convIt])
            newVt=vcat(dLogalad.Vt[2:horzLen+1,convIt-1],dLogalad.Vt[horzLen+1,convIt])
            newVs=vcat(dLogalad.Vs[(N+1):(N*(horzLen+1)),convIt-1],dLogalad.Vs[((N*horzLen)+1):(N*(horzLen+1)),convIt])
        end
    end

    global prevLam=round.(newLam,sigdigits=roundSigFigs)
    global prevVu=round.(newVu,sigdigits=roundSigFigs)
    global prevVz=round.(newVz,sigdigits=roundSigFigs)
    global prevVt=round.(newVt,sigdigits=roundSigFigs)
    global prevVs=round.(newVs,sigdigits=roundSigFigs)
    #global ρALADp=round.(ogρ,digits=2)

    return nothing
end

function pwlEValad(maxIt::Int,evS,cSave::centralLogStruct,slack::Bool,eqForm::Bool,roundSigFigs::Int,silent::Bool)
    # wrapper function for MPC PWL ALADIN

    horzLen=evS.K1
    K=evS.K
    N=evS.N
    S=evS.S
    dSol=solutionStruct(K=K,N=N,S=S)
    dCM=convMetricsStruct(maxIt=maxIt,logLength=length(saveLogInd))
    stepI=1

    p = plot(2,label=["Central" "ALAD"])
    Juno.progress() do id
        for stepI=1:K
            @info "$(Dates.format(Dates.now(),"HH:MM:SS")): $(stepI) of $(K)" progress=stepI/K _id=id
            @printf "%s: time step %g of %g...." Dates.format(Dates.now(),"HH:MM:SS") stepI K
            try
                dSol.timeT[stepI]=@elapsed runEVALADStep(stepI,maxIt,evS,dSol,dCM,cSave,eqForm,roundSigFigs,silent)
                @printf "convIt: %g\n" dSol.convIt[stepI,1]
            catch e
                @printf "error: %s" e
                break
            end

            push!(p, stepI, [cSol.lamCoupl[stepI], dSol.lamCoupl[stepI]])
        	display(p)
        end
    end

    objFun(sn,u)=sum(sum((sn[k,n]-1)^2*evS.Qsi[n,k] for n=1:N) +sum((u[k,n])^2*evS.Ri[n,k] for n=1:N) for k=1:evS.K)
    dSol.objVal[1,1]=objFun(dSol.Sn,dSol.Un)

    return dSol, dCM
end
