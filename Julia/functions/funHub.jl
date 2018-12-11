# hub functions


#Central
function runHubCentralStep(stepI,hubS,cSol,mode,silent)

    global t0
    global e0
    H=hubS.H
    K=hubS.K
    Nh=hubS.Nh
    horzLen=min(hubS.K1,K-stepI)

    eMax=hubS.eMax[stepI:(stepI+horzLen),:]
    uMax=hubS.uMax[stepI:(stepI+horzLen),:]
    eDepart_min=hubS.eDepart_min[stepI:(stepI+horzLen),:]
    eArrive_pred=hubS.eArrive_pred[stepI:(stepI+horzLen),:]
    eArrive_actual=hubS.eArrive_actual[stepI:(stepI+horzLen),:]
    slackMax=hubS.slackMax[stepI:(stepI+horzLen),:]

    if silent
        #cModel = Model(solver = IpoptSolver(print_level=0))
        #cModel = Model(solver = GurobiSolver(OutputFlag=0,QCPDual=1))
        cModel = Model(solver = GurobiSolver(OutputFlag=0))
    else
        #cModel = Model(solver = IpoptSolver())
        #cModel = Model(solver = GurobiSolver(QCPDual=1))
        cModel = Model(solver = GurobiSolver())

    end

    @variable(cModel,u[1:(horzLen+1),1:H])
    @variable(cModel,t[1:(horzLen+1)])
    @variable(cModel,e[1:(horzLen+1),1:H])
    @variable(cModel,eΔ[1:(horzLen+1),1:H])

    #objective
    @objective(cModel,Min,sum(sum(hubS.Qh[h]*(e[k,h]-eMax[k,h])^2+hubS.Rh[h]*u[k,h]^2 for k=1:horzLen+1) for h=1:H))

    #transformer constraints
    @constraint(cModel,upperTCon,t.<=hubS.Tmax)
    @constraint(cModel,t.>=0)
    if mode=="NL"
        @variable(cModel,itotal[1:(horzLen+1)])
        @constraint(cModel,itotal.<=hubS.ItotalMax)
        @constraint(cModel,itotal.>=0)
        @NLconstraint(cModel,tempCon1,t[1,1]==hubS.τP*t0+hubS.γP*(itotal[1])^2+hubS.ρP*hubS.Tamb[stepI,1])
        @NLconstraint(cModel,tempCon2[k=1:horzLen],t[k+1,1]==hubS.τP*t[k,1]+hubS.γP*(itotal[k+1])^2+hubS.ρP*hubS.Tamb[stepI+k,1])
        #coupling constraint
        @constraint(cModel,currCon[k=1:horzLen+1],0==-sum(u[k,h] for h=1:H)-hubS.iD_pred[stepI+(k-1)]+itotal[k,1])
    elseif mode=="relax1"
        @variable(cModel,itotal[1:(horzLen+1)])
        @constraint(cModel,itotal.<=hubS.ItotalMax)
        @constraint(cModel,itotal.>=0)
        @constraint(cModel,tempCon1,t[1,1]>=hubS.τP*t0+hubS.γP*(itotal[1])^2+hubS.ρP*hubS.Tamb[stepI,1])
        @constraint(cModel,tempCon2[k=1:horzLen],t[k+1,1]>=hubS.τP*t[k,1]+hubS.γP*(itotal[k+1])^2+hubS.ρP*hubS.Tamb[stepI+k,1])
        #coupling constraint
        @constraint(cModel,currCon[k=1:horzLen+1],0==-sum(u[k,h] for h=1:H)-hubS.iD_pred[stepI+(k-1)]+itotal[k,1])
    elseif mode=="PWL"
        @variable(cModel,z[1:(horzLen+1),1:hubS.S])
        @constraint(cModel,tempCon1,t[1,1]==hubS.τP*t0+hubS.γP*hubS.deltaI*sum((2*s-1)*z[1,s] for s=1:hubS.S)+hubS.ρP*hubS.Tamb[stepI,1])
        @constraint(cModel,tempCon2[k=1:horzLen],t[k+1,1]==hubS.τP*t[k,1]+hubS.γP*hubS.deltaI*sum((2*s-1)*z[k+1,s] for s=1:hubS.S)+hubS.ρP*hubS.Tamb[stepI+k,1])
        @constraint(cModel,z.>=0)
        @constraint(cModel,z.<=hubS.deltaI)
        #coupling constraint
        @constraint(cModel,currCon[k=1:horzLen+1],0==-sum(u[k,h] for h=1:H)-hubS.iD_pred[stepI+(k-1)]+sum(z[k,s] for s=1:hubS.S))
    end

    #hub constraints
    @constraint(cModel,stateCon1[h=1:H],e[1,h]==e0[h]+hubS.ηP[h]*u[1,h]-(eDepart_min[1,h]+eΔ[1,h])+eArrive_pred[1,h])
    @constraint(cModel,stateCon[k=1:horzLen,h=1:H],e[k+1,h]==e[k,h]+hubS.ηP[h]*u[k+1,h]-(eDepart_min[k+1,h]+eΔ[k+1,h])+eArrive_pred[k+1,h])
    @constraint(cModel,e.>=0)
    @constraint(cModel,eMaxCon[k=1:horzLen+1,h=1:H],e[k,h]<=eMax[k,h])
    @constraint(cModel,uMaxCon,u.<=uMax)
    @constraint(cModel,u.>=0)
    @constraint(cModel,eΔ.<=slackMax)
    @constraint(cModel,eΔ.>=0)

    TT = stdout # save original stdout stream
    redirect_stdout()
    status = solve(cModel)
    redirect_stdout(TT)
    @assert status==:Optimal "Central Hub optimization not solved to optimality"

    eRaw=getvalue(e)
    uRaw=getvalue(u)
    tRaw=getvalue(t)
    lambdaCurr=-getdual(currCon)
    extraE=getvalue(eΔ)
    #
    # p1nl=plot(eRaw,xlabel="Time",ylabel="Energy",xlims=(1,horzLen+1),label="hub energy")
    # plot!(p1nl,eMax,label="hub max")
    # plot!(p1nl,eMax*.8,label="minimum departure")
    #
    # p2nl=plot(uRaw,xlabel="Time",ylabel="Hub Current (kA)",legend=false,xlims=(1,horzLen+1))
    # plot!(p2nl,uMax)
    #
    # p3nl=plot(1:horzLen+1,tRaw,label="XFRM Temp",xlims=(1,horzLen+1),xlabel="Time",ylabel="Temp (K)")
    # plot!(p3nl,1:horzLen+1,Tmax*ones(horzLen+1),label="XFRM Limit",line=(:dash,:red))
    # p4nl=plot(1:horzLen+1,lambdaCurr,label="Time",ylabel=raw"Lambda ($/kA)",xlims=(1,horzLen+1),legend=false)
    #
    # plot(sum(eMax,dims=2),label="max")
    # plot!(sum(eRaw,dims=2),label="e")
    # plot(eDepart_min,label="min")
    # plot!(eDepart_min+extraE,label="actual")
    # plot!(eDepart_min+slackMax,label="max")


    #apply current and actual departures and arrivals
    nextU=uRaw[1,:]
    cSol.U[stepI,:]=nextU
    cSol.Lam[stepI,1]=lambdaCurr[1,1]
    cSol.E_depart[stepI,:]=eDepart_min[1,:]+extraE[1,:]
    cSol.E_arrive[stepI,:]=eArrive_actual[1,:]
    cSol.E[stepI,:]=e0[:]+hubS.ηP[:].*nextU-(cSol.E_depart[stepI,:])+cSol.E_arrive[stepI,:]
    cSol.Tactual[stepI,1]=hubS.τP*t0+hubS.γP*(sum(nextU)+hubS.iD_actual[stepI,1])^2+hubS.ρP*hubS.Tamb[stepI,1]

    # new states
    t0=cSol.Tactual[stepI,1]
    e0=cSol.E[stepI,:]

    return nothing
end

function hubCentral(hubS::scenarioHubStruct,mode::String,silent::Bool)
    H=hubS.H
    K=hubS.K
    Nh=hubS.Nh
    cSol=hubSolutionStruct(K=K,H=H)

    for stepI=1:K
        @printf "%s: time step %g of %g....\n" Dates.format(Dates.now(),"HH:MM:SS") stepI K
        try
            runHubCentralStep(stepI,hubS,cSol,mode,silent)
        catch e
            @printf "error: %s" e
            break
        end
    end

    objFun(e,u)=sum(sum((e[k,h]-hubS.eMax[k,h])^2*hubS.Qh[1,h] for h=1:H) +sum((u[k,h])^2*hubS.Rh[1,h] for h=1:H) for k=1:K)
    cSol.objVal[1,1]=objFun(cSol.E,cSol.U)

    return cSol
end


#Dual Ascent
function localEVDual(hubInd::Int,p::Int,stepI::Int,hubS::scenarioHubStruct,dLog::hubItLogPWL)
    horzLen=min(hubS.K1,hubS.K-stepI)

    eMax=hubS.eMax[stepI:(stepI+horzLen),hubInd]
    uMax=hubS.uMax[stepI:(stepI+horzLen),hubInd]
    eDepart_min=hubS.eDepart_min[stepI:(stepI+horzLen),hubInd]
    eArrive_pred=hubS.eArrive_pred[stepI:(stepI+horzLen),hubInd]
    eArrive_actual=hubS.eArrive_actual[stepI:(stepI+horzLen),hubInd]
    slackMax=hubS.slackMax[stepI:(stepI+horzLen),hubInd]

    hubM=Model(solver = GurobiSolver())
    @variable(hubM,u[1:horzLen+1])
    @variable(hubM,e[1:horzLen+1])
    @variable(hubM,eΔ[1:(horzLen+1)])
    objExp=sum((e[k,1]-eMax[k,1])^2*hubS.Qh[1,hubInd]+(u[k,1])^2*hubS.Rh[1,hubInd]+prevLam[k,1]*u[k,1] for k=1:horzLen+1)
    @objective(hubM,Min,objExp)

    #hub constraints
    @constraint(hubM,stateCon1,e[1,1]==e0[hubInd]+hubS.ηP[1,hubInd]*u[1]-(eDepart_min[1]+eΔ[1])+eArrive_pred[1])
    @constraint(hubM,stateCon[k=1:horzLen],e[k+1,1]==e[k,1]+hubS.ηP[1,hubInd]*u[k+1]-(eDepart_min[k+1]+eΔ[k+1])+eArrive_pred[k+1])
    @constraint(hubM,e.>=0)
    @constraint(hubM,eMaxCon,e.<=eMax)
    @constraint(hubM,uMaxCon,u.<=uMax)
    @constraint(hubM,u.>=0)
    @constraint(hubM,eΔ.<=slackMax)
    @constraint(hubM,eΔ.>=0)

    TT = stdout # save original stdout stream
    redirect_stdout()
    statusM = solve(hubM)
    redirect_stdout(TT)
    @assert statusM==:Optimal "Hub optimization not solved to optimality"

    dLog.U[:,hubInd,p]=getvalue(u)
    dLog.E[:,hubInd,p]=getvalue(e)
    return nothing
end

function localXFRMDual(p::Int,stepI::Int,hubS::scenarioHubStruct,dLog::hubItLogPWL,mode::String)
    S=hubS.S
    horzLen=min(hubS.K1,hubS.K-stepI)

    coorM=Model(solver = GurobiSolver())
    @variable(coorM,t[1:(horzLen+1)])
    @constraint(coorM,upperTCon,t.<=hubS.Tmax)
    @constraint(coorM,t.>=0)
    if mode=="NL"
        @variable(coorM,itotal[1:(horzLen+1)])
        @objective(coorM,Min,sum(prevLam[k,1]*itotal[k] for k=1:(horzLen+1)))
        @constraint(coorM,itotal.<=hubS.ItotalMax)
        @constraint(coorM,itotal.>=0)
        @NLconstraint(coorM,tempCon1,t[1,1]==hubS.τP*t0+hubS.γP*(itotal[1])^2+hubS.ρP*hubS.Tamb[stepI,1])
        @NLconstraint(coorM,tempCon2[k=1:horzLen],t[k+1,1]==hubS.τP*t[k,1]+hubS.γP*(itotal[k+1])^2+hubS.ρP*hubS.Tamb[stepI+k,1])
    elseif mode=="relax1"
        @variable(coorM,itotal[1:(horzLen+1)])
        @objective(coorM,Min,sum(prevLam[k,1]*itotal[k] for k=1:(horzLen+1)))
        @constraint(coorM,itotal.<=hubS.ItotalMax)
        @constraint(coorM,itotal.>=0)
        @constraint(coorM,tempCon1,t[1,1]>=hubS.τP*t0+hubS.γP*(itotal[1])^2+hubS.ρP*hubS.Tamb[stepI,1])
        @constraint(coorM,tempCon2[k=1:horzLen],t[k+1,1]>=hubS.τP*t[k,1]+hubS.γP*(itotal[k+1])^2+hubS.ρP*hubS.Tamb[stepI+k,1])
    elseif mode=="PWL"
        @variable(coorM,z[1:(horzLen+1),1:hubS.S])
        @objective(coorM,Min,sum(prevLam[k,1]*sum(-z[k,s] for s=1:S) for k=1:(horzLen+1)))
        @constraint(coorM,tempCon1,t[1,1]==hubS.τP*t0+hubS.γP*hubS.deltaI*sum((2*s-1)*z[1,s] for s=1:hubS.S)+hubS.ρP*hubS.Tamb[stepI,1])
        @constraint(coorM,tempCon2[k=1:horzLen],t[k+1,1]==hubS.τP*t[k,1]+hubS.γP*hubS.deltaI*sum((2*s-1)*z[k+1,s] for s=1:hubS.S)+hubS.ρP*hubS.Tamb[stepI+k,1])
        @constraint(coorM,z.>=0)
        @constraint(coorM,z.<=hubS.deltaI)
    end

    TT = stdout # save original stdout stream
    redirect_stdout()
    statusC = solve(coorM)
    redirect_stdout(TT)
    @assert statusC==:Optimal "Dual Ascent XFRM optimization not solved to optimality"

     dLog.Tpwl[:,1,p]=getvalue(t)
     dLog.Z[:,:,p]=getvalue(z)

    return nothing
end

function runHubDualIt(p,stepI,hubS,dLog,dSol,cSol,mode,silent)

    global t0
    global e0
    global prevLam

    H=hubS.H
    K=hubS.K
    horzLen=min(hubS.K1,K-stepI)

    alpha0 = 5e5 #for kA
    alphaDivRate=2
    minAlpha=1e-6
    convChk = 1e-4

    #solve subproblem for each EV
    @sync @distributed for hubInd=1:H
        localEVDual(hubInd,p,stepI,hubS,dLog)
    end

    localXFRMDual(p,stepI,hubS,dLog,mode)

    #grad of lagragian
    for k=1:horzLen+1
        dLog.uSum[k,1,p]=sum(dLog.U[k,h,p] for h=1:H)
        dLog.zSum[k,1,p]=sum(dLog.Z[k,s,p] for s=1:S)
        dLog.couplConst[k,1,p]=dLog.uSum[k,1,p]+ hubS.iD_pred[stepI+(k-1),1] - dLog.zSum[k,1,p]
    end
    #dCM.couplConst[p,1]=norm(dLog.couplConst[:,p],2)

    #update lambda
    alphaP= max(alpha0/ceil(p/alphaDivRate),minAlpha)
    #dLog.itUpdate[1,p]=alphaP
    #alphaP= alphaP*alphaRate

    #dLog.Lam[:,p]=max.(prevLam[:,1]+alphaP*dLog.couplConst[:,p],0)
    dLog.Lam[:,1,p]=prevLam[:,1]+alphaP*dLog.couplConst[:,1,p]

    #calculate actual temperature from nonlinear model of XFRM
    for k=1:horzLen+1
        dLog.Itotal[k,1,p]=dLog.uSum[k,1,p] + hubS.iD_actual[stepI+(k-1),1]
    end
    dLog.Tactual[1,1,p]=hubS.τP*t0+hubS.γP*dLog.Itotal[1,1,p]^2+hubS.ρP*hubS.Tamb[stepI,1]
    for k=1:horzLen
        dLog.Tactual[k+1,1,p]=hubS.τP*dLog.Tactual[k,1,p]+hubS.γP*dLog.Itotal[k+1,1,p]^2+hubS.ρP*hubS.Tamb[stepI+k,1]
    end

    #check convergence
    # objFun(sn,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*hubS.Qsi[n,1]     for n=1:N) for k=1:horzLen+1) +
    #                 sum(sum((u[(k-1)*N+n,1])^2*hubS.Ri[n,1]           for n=1:N) for k=1:horzLen+1)
    # dLog.objVal[1,p]=objFun(dLog.Sn[:,p],dLog.Un[:,p])
    #fGap=abs(dLog.objVal[1,p] -cSol.objVal[1,1])
    #snGap=norm((dLog.Sn[:,p]-cSol.Sn),2)
    #unGap=norm((dLog.Un[:,p]-cSol.Un),2)
    itGap = norm(dLog.Lam[:,1,p]-prevLam[:,1],2)
    convGap = norm(dLog.Lam[:,1,p]-cSol.Lam[stepI:(horzLen+stepI)],2)
    #dCM.obj[p,1]=fGap
    #dCM.sn[p,1]=snGap
    #dCM.un[p,1]=unGap
    #dCM.lamIt[p,1]=itGap
    #dCM.lam[p,1]=convGap
    if(itGap <= convChk )
        if !silent @printf "Converged after %g iterations\n" p end
        convIt=p
        return true
    else
        if !silent
            @printf "lastGap %e after %g iterations\n" itGap p
            @printf "convGap %e after %g iterations\n\n" convGap p
            #@printf "snGap   %e after %g iterations\n" snGap p
            #@printf "unGap   %e after %g iterations\n" unGap p
            #@printf("fGap    %e after %g iterations\n\n",fGap,p)
        end
        prevLam=dLog.Lam[:,1,p]
        return false
    end
end

function runHubDualStep(stepI,maxIt,hubS,dSol,cSol,mode,silent)
    K=hubS.K
    S=hubS.S
    horzLen=min(hubS.K1,K-stepI)
    dLog=hubItLogPWL(horzLen=horzLen,H=hubS.H,S=hubS.S)
    for p=1:maxIt
        cFlag=runHubDualIt(p,stepI,hubS,dLog,dSol,cSol,mode,silent)
        global convIt=p
        if cFlag
            break
        end
    end

    # # convergence plots
    # halfCI=Int(floor(convIt/2))
    # CList=reshape([range(colorant"blue", stop=colorant"yellow",length=halfCI);
    #                range(colorant"yellow", stop=colorant"red",length=convIt-halfCI)], 1, convIt);
    #
    # uSumPlotd=plot(dLog.uSum[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Current Sum (kA)",xlims=(0,horzLen+1),legend=false)
    # plot!(uSumPlotd,1:horzLen+1,cSol.uSum,seriescolor=:black,linewidth=2,linealpha=0.8)
    #
    # # uSumPlotd=plot(dLog.uSum[:,1:convIt],palette=:greens,line_z=(1:convIt)',legend=false,colorbar=:right,colorbar_title="Iteration",
    # #      xlabel="Time",ylabel="Current Sum (kA)",xlims=(0,horzLen+1))
    # # plot!(uSumPlotd,1:horzLen+1,cSol.uSum,seriescolor=:black,linewidth=2,linealpha=0.8)
    #
    #
    # zSumPlotd=plot(dLog.zSum[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Z sum",xlims=(0,horzLen+1),legend=false)
    # plot!(zSumPlotd,1:horzLen+1,cSol.zSum,seriescolor=:black,linewidth=2,linealpha=0.8)
    #
    # lamPlotd=plot(dLog.Lam[:,1,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Lambda",xlims=(0,horzLen+1),legend=false)
    # plot!(lamPlotd,cSol.Lam,seriescolor=:black,linewidth=2,linealpha=0.8)
    #
    # #plot(dLog.Lam[:,1:convIt],color=:RdYlBu,line_z=(1:convIt)')
    #
    # constPlot2=plot(dLog.couplConst[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="curr constraint diff",xlims=(0,horzLen+1),legend=false)
    #
    # #convergence metric plots
    # fPlot=plot(1:convIt,dCM.obj[1:convIt,1],xlabel="Iteration",ylabel="obj function gap",xlims=(1,convIt),legend=false,yscale=:log10)
    # convItPlot=plot(1:convIt,dCM.lamIt[1:convIt,1],xlabel="Iteration",ylabel="2-Norm Lambda Gap",xlims=(1,convIt),legend=false,yscale=:log10)
    # convPlot=plot(1:convIt,dCM.lam[1:convIt,1],xlabel="Iteration",ylabel="central lambda gap",xlims=(2,convIt),legend=false,yscale=:log10)
    # constPlot=plot(1:convIt,dCM.couplConst[1:convIt,1],xlabel="Iteration",ylabel="curr constraint Gap",xlims=(2,convIt),legend=false,yscale=:log10)
    #

    #save current state and update for next timeSteps
    dSol.Tpwl[stepI,1]=dLog.Tpwl[1,1,convIt]
    dSol.U[stepI,:]=dLog.U[1,1:H,convIt]
    dSol.E[stepI,:]=dLog.E[1,1:H,convIt]
    #dSol.Z[stepI,:]=dLog.Z[1:S,convIt]
    dSol.Itotal[stepI,1]=dLog.Itotal[1,1,convIt]
    dSol.Tactual[stepI,1]=dLog.Tactual[1,1,convIt]
    dSol.convIt[stepI,1]=convIt

    # new states
    global t0=dSol.Tactual[stepI,1]
    global e0=dSol.E[stepI,:]

    if convIt==1
        dSol.Lam[stepI,1]=prevLam[1,1]
        if stepI+horzLen==hubS.K
            newLam=prevLam[2:horzLen+1,1]
        else
            newLam=vcat(prevLam[2:horzLen+1,1],prevLam[horzLen+1,1])
        end
    else
        dSol.Lam[stepI,1]=dLog.Lam[1,1,convIt-1]
        if stepI+horzLen==hubS.K
            newLam=dLog.Lam[2:horzLen+1,1,convIt-1]
        else
            newLam=vcat(dLog.Lam[2:horzLen+1,1,convIt-1],dLog.Lam[horzLen+1,1,convIt-1])
        end
    end
    global prevLam=newLam
    return nothing
end

function hubDual(maxIt::Int,hubS::scenarioHubStruct,cSol::hubSolutionStruct,mode::String,silent::Bool)
    H=hubS.H
    K=hubS.K
    Nh=hubS.Nh
    dSol=hubSolutionStruct(K=K,H=H)

    for stepI=1:K
        @printf "%s: time step %g of %g...." Dates.format(Dates.now(),"HH:MM:SS") stepI K
        try
            runHubDualStep(stepI,maxIt,hubS,dSol,cSol,mode,silent)
            @printf "convIt: %g\n" dSol.convIt[stepI,1]
        catch e
            @printf "error: %s" e
            break
        end
    end

    objFun(e,u)=sum(sum((e[k,h]-hubS.eMax[k,h])^2*hubS.Qh[1,h] for h=1:H) +sum((u[k,h])^2*hubS.Rh[1,h] for h=1:H) for k=1:K)
    dSol.objVal[1,1]=objFun(dSol.E,dSol.U)

    return dSol
end



#ALADIN
function localEVALAD(evInd::Int,p::Int,stepI::Int,σU::Array{Float64,2},σE::Array{Float64,2},hubS::scenarioHubStruct,dLogalad::hubItLogPWL)
    horzLen=min(hubS.K1,hubS.K-stepI)

    tolU=1e-6
    tolS=1e-8

    eMax=hubS.eMax[stepI:(stepI+horzLen),hubInd]
    uMax=hubS.uMax[stepI:(stepI+horzLen),hubInd]
    eDepart_min=hubS.eDepart_min[stepI:(stepI+horzLen),hubInd]
    eArrive_pred=hubS.eArrive_pred[stepI:(stepI+horzLen),hubInd]
    eArrive_actual=hubS.eArrive_actual[stepI:(stepI+horzLen),hubInd]
    slackMax=hubS.slackMax[stepI:(stepI+horzLen),hubInd]


    hubM=Model(solver = GurobiSolver())
    @variable(hubM,u[1:horzLen+1])
    @variable(hubM,e[1:horzLen+1])
    @variable(hubM,eΔ[1:(horzLen+1)])
    objExp=sum((e[k,1]-eMax[k,1])^2*hubS.Qh[1,hubInd]+(u[k,1])^2*hubS.Rh[1,hubInd]+
                            prevLam[k,1]*(u[k,1])+
                            ρALADp/2*(u[k,1]-prevVu[k,hubInd])*σU[hubInd,1]*(u[k,1]-prevVu[k,hubInd])+
                            ρALADp/2*(eΔ[k,1]-prevVd[k,hubInd])*σE[hubInd,1]*(eΔ[k,1]-prevVd[k,hubInd])+
                            ρALADp/2*(e[k,1]-prevVe[k,hubInd])*σE[hubInd,1]*(e[k,1]-prevVe[k,hubInd]) for k=1:horzLen+1)

    @objective(hubM,Min,objExp)

    #hub constraints
    @constraint(hubM,stateCon1,e[1,1]==e0[hubInd]+hubS.ηP[1,hubInd]*u[1]-(eDepart_min[1]+eΔ[1])+eArrive_pred[1])
    @constraint(hubM,stateCon[k=1:horzLen],e[k+1,1]==e[k,1]+hubS.ηP[1,hubInd]*u[k+1]-(eDepart_min[k+1]+eΔ[k+1])+eArrive_pred[k+1])
    @constraint(hubM,e.>=0)
    @constraint(hubM,eMaxCon,e.<=eMax)
    @constraint(hubM,uMaxCon,u.<=uMax)
    @constraint(hubM,u.>=0)
    @constraint(hubM,eΔ.<=slackMax)
    @constraint(hubM,eΔ.>=0)

    TT = stdout # save original stdout stream
    redirect_stdout()
    statusM = solve(hubM)
    redirect_stdout(TT)
    @assert statusM==:Optimal "Hub optimization not solved to optimality"


    uVal=getvalue(u)
    eVal=getvalue(e)
    eΔVal=getvalue(eΔ)

    cValMax=abs.(uVal.-uMax).<tolU
    cValMin=abs.(uVal.-0).<tolU
    dLogalad.Cuu[:,hubInd,p]=1cValMax
    dLogalad.Cul[:,hubInd,p]=-1cValMin

    cValMax=abs.(eVal.-eMax).<tolE
    cValMin=abs.(eVal.-0).<tolE
    dLogalad.Ceu[:,hubInd,p]=1cValMax
    dLogalad.Cel[:,hubInd,p]=-1cValMin

    cValMax=abs.(eΔVal.-slackMax).<tolE
    cValMin=abs.(eΔVal.-0).<tolE
    dLogalad.Cdu[:,hubInd,p]=1cValMax
    dLogalad.Cdl[:,hubInd,p]=-1cValMin

    #dLogalad.slackSn[evInd]= if slack getvalue(slackSn) else 0 end
    dLogalad.E[:,hubInd,p]=eVal
    dLogalad.U[:,hubInd,p]=uVal

    dLogalad.Gu[:,hubInd,p]=round.(2*hubS.Rh[hubInd,1]*uVal,digits=4)
    dLogalad.Ge[:,hubInd,p]=round.(2*hubS.Qh[hubInd,1]*eVal.-2*hubS.Qh[hubInd,1]*eMax,digits=4)
    dLogalad.GeΔ[:,hubInd,p]=0
    return nothing
end

function localXFRMALAD(p::Int,stepI::Int,σZ::Float64,σT::Float64,hubS::scenarioHubStruct,dLogalad::hubItLogPWL)
    S=hubS.S
    horzLen=min(hubS.K1,hubS.K-stepI)

    tolT=1e-6
    tolZ=1e-6

    #N+1 decoupled problem aka transformer current
    coorM=Model(solver = GurobiSolver())
    @variable(coorM,t[1:(horzLen+1)])
    @constraint(coorM,upperTCon,t.<=hubS.Tmax)
    @constraint(coorM,t.>=0)
    if mode=="NL"
        # @variable(coorM,itotal[1:(horzLen+1)])
        # @objective(coorM,Min,sum(prevLam[k,1]*itotal[k] for k=1:(horzLen+1)))
        # @constraint(coorM,itotal.<=hubS.ItotalMax)
        # @constraint(coorM,itotal.>=0)
        # @NLconstraint(coorM,tempCon1,t[1,1]==hubS.τP*t0+hubS.γP*(itotal[1])^2+hubS.ρP*hubS.Tamb[stepI,1])
        # @NLconstraint(coorM,tempCon2[k=1:horzLen],t[k+1,1]==hubS.τP*t[k,1]+hubS.γP*(itotal[k+1])^2+hubS.ρP*hubS.Tamb[stepI+k,1])
    elseif mode=="relax1"
        # @variable(coorM,itotal[1:(horzLen+1)])
        # @objective(coorM,Min,sum(prevLam[k,1]*itotal[k] for k=1:(horzLen+1)))
        # @constraint(coorM,itotal.<=hubS.ItotalMax)
        # @constraint(coorM,itotal.>=0)
        # @constraint(coorM,tempCon1,t[1,1]>=hubS.τP*t0+hubS.γP*(itotal[1])^2+hubS.ρP*hubS.Tamb[stepI,1])
        # @constraint(coorM,tempCon2[k=1:horzLen],t[k+1,1]>=hubS.τP*t[k,1]+hubS.γP*(itotal[k+1])^2+hubS.ρP*hubS.Tamb[stepI+k,1])
    elseif mode=="PWL"
        @variable(coorM,z[1:(horzLen+1),1:hubS.S])
        objExp=sum(-prevLam[k,1]*sum(z[k,s] for s=1:S)+
                  ρALADp/2*σZ*(sum(z[k,s] for s=1:S)-sum(prevVz[k,s] for s=1:S))^2+
                  ρALADp/2*σT*(t[k]-prevVt[k,1])^2  for k=1:(horzLen+1))
        @objective(coorM,Min, objExp)
        @constraint(coorM,tempCon1,t[1,1]==hubS.τP*t0+hubS.γP*hubS.deltaI*sum((2*s-1)*z[1,s] for s=1:hubS.S)+hubS.ρP*hubS.Tamb[stepI,1])
        @constraint(coorM,tempCon2[k=1:horzLen],t[k+1,1]==hubS.τP*t[k,1]+hubS.γP*hubS.deltaI*sum((2*s-1)*z[k+1,s] for s=1:hubS.S)+hubS.ρP*hubS.Tamb[stepI+k,1])
        @constraint(coorM,z.>=0)
        @constraint(coorM,z.<=hubS.deltaI)
    end

    TT = stdout # save original stdout stream
    redirect_stdout()
    statusC = solve(coorM)
    redirect_stdout(TT)
    @assert statusC==:Optimal "Dual Ascent XFRM optimization not solved to optimality"

    #kappaMax=-getdual(pwlKappaMax)
    #kappaMin=-getdual(pwlKappaMin)
    #tMax=-getdual(upperTCon)
    #tMin=-getdual(lowerTCon)
    #lambdaTemp=[-getdual(tempCon1);-getdual(tempCon2)]
    zVal=getvalue(z)
    tVal=getvalue(t)

    cValMax=abs.(zVal.-hubS.deltaI).<tolZ
    cValMin=abs.(zVal.-0).<tolZ
    dLogalad.Czu[:,:,p]=1cValMax
    dLogalad.Czl[:,:,p]=-1cValMin

    cValMax=abs.(tVal.-hubS.Tmax).<tolT
    cValMin=abs.(tVal.-0).<tolT
    dLogalad.Ctu[:,:,p]=1cValMax
    dLogalad.Ctl[:,:,p]=-1cValMin

    dLogalad.Tpwl[:,:,p]=tVal
    dLogalad.Z[:,:,p]=zVal
    dLogalad.Gz[:,:,p].=0
    dLogalad.Gt[:,:,p].=0

    return nothing
end

function coordALAD(p::Int,stepI::Int,μALADp::Float64,hubS::scenarioHubStruct,dLogalad::hubItLogPWL)
    H=hubS.H
    S=hubS.S
    horzLen=min(hubS.K1,hubS.K-stepI)

    Hu=2*hubS.Rh
    Hs=2*hubS.Qh
    # Hu=2*hubS.Ri *((1.5-2.5)*rand()+2.5)
    # Hs=2*hubS.Qsi *((1.5-2.5)*rand()+2.5)
    Hz=0
    Ht=0
    HeΔ=0
    # Hz=1e-6
    # Ht=1e-6
    ρRate=1.1
    ρALADmax=4e5

    #coupled QP
    cM = Model(solver = GurobiSolver(NumericFocus=3))
    @variable(cM,dU[1:(horzLen+1),1:H])
    @variable(cM,dE[1:(horzLen+1),1:H])
    @variable(cM,dEΔ[1:(horzLen+1),1:H])
    @variable(cM,dZ[1:(horzLen+1),1:S])
    @variable(cM,dT[1:(horzLen+1)])
    @variable(cM,relaxS[1:(horzLen+1)])

    objExp=sum(sum(0.5*dU[k,h]^2*Hu[h,1]+dLogalad.Gu[k,h,p]*dU[k,h] +
                   0.5*dE[k,h]^2*He[h,1]+dLogalad.Ge[k,h,p]*dE[k,h] +
                   0.5*dEΔ[k,h]^2*HeΔ[h,1]+dLogalad.GeΔ[k,h,p]*dEΔ[k,h] for h=1:H) +
               sum(0.5*dZ[k,s]^2*Hz for s=1:S)+
               0.5*dT[k,1]^2*Ht   for k=1:(horzLen+1))
    objExp=objExp+prevLam[:,1]'*relaxS+μALADp/2*sum(relaxS[k,1]^2 for k=1:horzLen+1)
    objExp=objExp+dot(dLogalad.Gz[:,1,p],dZ)+dot(dLogalad.Gt[:,1,p],dT)

    @objective(cM,Min,objExp)

    Up=round.(dLogalad.U[:,:,p],digits=8)
    Zp=round.(dLogalad.Z[:,:,p],digits=8)
    @constraint(cM,currCon[k=1:horzLen+1],sum(Up[k,h]+dU[k,h] for h=1:H)-
                                          sum(Zp[k,s]+dZ[k,s] for s=1:S)==
                                          -hubS.iD_pred[stepI+(k-1)]+relaxS[k,1])
    #local equality constraints C*(X+deltaX)=0 is same as C*deltaX=0 since we already know CX=0
    @constraint(cM,stateCon1[h=1:H],dE[1,h]==hubS.ηP[h]*dU[1,h]-dEΔ[1,h])
    @constraint(cM,stateCon2[k=1:horzLen,h=1:H],dE[k+1,h]==dE[k,h]+hubS.ηP[h]*dU[k+1,h]-dEΔ[k+1,h])
    @constraint(cM,tempCon1,dT[1,1]==hubS.γP*hubS.deltaI*sum((2*s-1)*dZ[1,s] for s=1:S))
    @constraint(cM,tempCon2[k=1:horzLen],dT[k+1,1]==hubS.τP*dT[k,1]+hubS.γP*hubS.deltaI*sum((2*s-1)*dZ[k+1,s] for s=1:S))

    #active local constraints
    if eqForm
        @constraint(cM,dLogalad.Czu[:,:,p].*dZ.==0)
        @constraint(cM,dLogalad.Cuu[:,:,p].*dU.==0)
        @constraint(cM,dLogalad.Ceu[:,:,p].*dE.==0)
        @constraint(cM,dLogalad.Cdu[:,:,p].*dEΔ.==0)
        @constraint(cM,dLogalad.Ctu[:,:,p].*dT.==0)
        @constraint(cM,dLogalad.Czl[:,:,p].*dZ.==0)
        @constraint(cM,dLogalad.Cul[:,:,p].*dU.==0)
        @constraint(cM,dLogalad.Cdl[:,:,p].*dEΔ.==0)
        @constraint(cM,dLogalad.Cel[:,:,p].*dE.==0)
        @constraint(cM,dLogalad.Ctl[:,:,p].*dT.==0)
    else
        @constraint(cM,dLogalad.Czu[:,:,p].*dZ.<=0)
        @constraint(cM,dLogalad.Cuu[:,:,p].*dU.<=0)
        @constraint(cM,dLogalad.Ceu[:,:,p].*dE.<=0)
        @constraint(cM,dLogalad.Cdu[:,:,p].*dEΔ.<=0)
        @constraint(cM,dLogalad.Ctu[:,:,p].*dT.<=0)
        @constraint(cM,dLogalad.Czl[:,:,p].*dZ.<=0)
        @constraint(cM,dLogalad.Cul[:,:,p].*dU.<=0)
        @constraint(cM,dLogalad.Cdl[:,:,p].*dEΔ.<=0)
        @constraint(cM,dLogalad.Cel[:,:,p].*dE.<=0)
        @constraint(cM,dLogalad.Ctl[:,:,p].*dT.<=0)
    end


    TT = stdout # save original stdout stream
    redirect_stdout()
    statusM = solve(cM)
    redirect_stdout(TT)
    @assert statusM==:Optimal "ALAD Central QP optimization not solved to optimality"

    #update step
    # Lam[:,p+1]=-getdual(currCon)
    α1=1
    α2=1
    α3=1
    #α3=1/ceil(p/2)

    dLogalad.Lam[:,p]=prevLam[:,1]+α3*(-getdual(currCon)-prevLam[:,1])
    #dLogalad.Lam[:,p]=max.(prevLam[:,1]+α3*(-getdual(currCon)-prevLam[:,1]),0)
    dLogalad.Vu[:,:,p]=prevVu[:,:,1]+α1*(dLogalad.U[:,:,p]-prevVu[:,:,1])+α2*getvalue(dU)
    dLogalad.Vz[:,:,p]=prevVz[:,:,1]+α1*(dLogalad.Z[:,:,p]-prevVz[:,:,1])+α2*getvalue(dZ)
    dLogalad.Ve[:,:,p]=prevVe[:,:,1]+α1*(dLogalad.E[:,:,p]-prevVe[:,:,1])+α2*getvalue(dE)
    dLogalad.VeΔ[:,:,p]=prevVd[:,:,1]+α1*(dLogalad.VeΔ[:,:,p]-prevVd[:,:,1])+α2*getvalue(dEΔ)
    dLogalad.Vt[:,:,p]=prevVt[:,:,1]+α1*(dLogalad.Tpwl[:,:,p]-prevVt[:,:,1])+α2*getvalue(dT)

    #dCMalad.lamIt[p,1]=norm(dLogalad.Lam[:,p]-prevLam[:,1],2)
    #dCMalad.lam[p,1]=norm(dLogalad.Lam[:,p]-cSol.lamCoupl[stepI:(horzLen+stepI)],2)
    # if !silent
    #     @printf "lastGap    %e after %g iterations\n" dCMalad.lamIt[p,1] p
    #     @printf "convLamGap %e after %g iterations\n\n" dCMalad.lam[p,1] p
    # end

    dLogalad.itUpdate[1,p]=min(ρALADp*ρRate,ρALADmax) #increase ρ every iteration
    #μALADp=min(μALADp*μRate,μALADmax) #increase μ every iteration
    #ΔY[1,p]=norm(vcat(getvalue(dUn),getvalue(dZ),getvalue(dSn),getvalue(dT)),Inf)
    return nothing
end

function runEVALADIt(p,stepI,hubS,dLogalad,dCMalad,dSol,cSol,eqForm,silent)
    H=hubS.H
    K=hubS.K
    S=hubS.S
    horzLen=min(hubS.K1,K-stepI)

    #initialize with current states
    global s0
    global t0
    global prevLam
    global prevVu
    global prevVz
    global prevVt
    global prevVs
    global ρALADp

    #other parameters
    epsilon = 1e-3

    #ALADIN tuning
    if eqForm
        #println("Running Eq ALADIN")
        scalingF=1
    else
        #println("Running ineq ALADIN")
        scalingF=1e-4
    end
    σZ=scalingF*1e1
    σT=scalingF
    σU=ones(H,1)
    σE=ones(H,1)/10#for kA

    μALADp=1e8
    # μALAD=1e8
    # μRate=1
    # μALADmax=2e9

    #solve decoupled
    @sync @distributed for hubInd=1:H
        localEVALAD(hubInd,p,stepI,σU,σE,hubS,dLogalad)
    end

    localXFRMALAD(p,stepI,σZ,σT,hubS,dLogalad)

    for k=1:horzLen+1
        dLogalad.uSum[k,1,p]=sum(dLogalad.U[k,h,p] for h=1:H)
        dLogalad.zSum[k,1,p]=sum(dLogalad.Z[k,s,p] for s=1:S)
        dLogalad.couplConst[k,1,p]=dLogalad.uSum[k,1,p] + hubS.iD_pred[stepI+(k-1),1] - dLogalad.zSum[k,1,p]
    end

    #calculate actual temperature from nonlinear model of XFRM
    for k=1:horzLen+1
        dLogalad.Itotal[k,1,p]=dLogalad.uSum[k,1,p] + hubS.iD_actual[stepI+(k-1),1]
    end
    dLogalad.Tactual[1,1,p]=hubS.τP*t0+hubS.γP*dLogalad.Itotal[1,1,p]^2+hubS.ρP*hubS.Tamb[stepI,1] #fix for mpc
    for k=1:horzLen
        dLogalad.Tactual[k+1,1,p]=hubS.τP*dLogalad.Tactual[k,1,p]+hubS.γP*dLogalad.Itotal[k,1,p]^2+hubS.ρP*hubS.Tamb[stepI+k,1]  #fix for mpc
    end

    #check for convergence
    constGap=norm(dLogalad.couplConst[:,1,p],1)
    cc=norm(vcat(σU[1]*(prevVu[:,:]-dLogalad.U[:,:,p]),σZ*(prevVz[:,:]-dLogalad.Z[:,:,p]),
                 σT*(prevVt[:,:]-dLogalad.Tpwl[:,:,p]),σS[1]*(prevVe[:,:]-dLogalad.E[:,:,p])),1)
    #cc=ρALAD*norm(vcat(repeat(σU,horzLen+1,1).*(Vu[:,p]-Un[:,p+1]),σZ*(Vz[:,p]-Z[:,p+1])),1)
    # objFun(sn,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*hubS.Qsi[n,1] for n=1:N) +
    #                  sum((u[(k-1)*N+n,1])^2*hubS.Ri[n,1]       for n=1:N) for k=1:horzLen+1)
    # dLogalad.objVal[1,p]=objFun(dLogalad.Sn[:,p],dLogalad.Un[:,p])
    #fGap= abs(dLogalad.objVal[1,p]-cSol.objVal[1,1])
    #snGap=norm((dLogalad.Sn[:,p]-cSol.Sn),2)
    #unGap=norm((dLogalad.Un[:,p]-cSol.Un),2)
    #dCMalad.obj[p,1]=fGap
    #dCMalad.sn[p,1]=snGap
    #dCMalad.un[p,1]=unGap
    #dCMalad.couplConst[p,1]=constGap
    #convCheck[p,1]=cc
    if  constGap<=epsilon && cc<=epsilon
        if !silent @printf "Converged after %g iterations\n" p end
        convIt=p
        #break
        return true
    else
        if !silent
            @printf "convCheck  %e after %g iterations\n" cc p
            @printf "constGap   %e after %g iterations\n" constGap p
            #@printf "snGap      %e after %g iterations\n" snGap p
            #@printf("fGap       %e after %g iterations\n",fGap,p)
        end
    end

    coordALAD(p,stepI,μALADp,hubS,dLogalad,dCMalad)

    #reset for next iteration
    prevVu=dLogalad.Vu[:,:,p]
    prevVs=dLogalad.Vs[:,:,p]
    prevVz=dLogalad.Vz[:,:,p]
    prevVt=dLogalad.Vt[:,:,p]
    prevLam=dLogalad.Lam[:,:,p]
    ρALADp=dLogalad.itUpdate[1,p]

    return false
end

function runHubALADStep(stepI,maxIt,hubS,dSol,cSol,mode,eqForm,silent)
    K=hubS.K
    S=hubS.S
    horzLen=min(hubS.K1,K-stepI)
    dLogalad=hubItLogPWL(horzLen=horzLen,H=hubS.H,S=hubS.S)
    for p=1:maxIt
        @printf "%git" p
        cFlag=runEVALADIt(p,stepI,hubS,dLogalad,dCMalad,dSol,cSol,mode,eqForm,silent)
        global convIt=p
        if cFlag
            break
        end
    end

    # xPlot=zeros(horzLen+1,N)
    # uPlot=zeros(horzLen+1,N)
    # for ii= 1:N
    # 	xPlot[:,ii]=dLogalad.Sn[collect(ii:N:length(dLogalad.Sn[:,convIt])),convIt]
    #     uPlot[:,ii]=dLogalad.Un[collect(ii:N:length(dLogalad.Un[:,convIt])),convIt]
    # end
    #
    # plot(uPlot)
    # plot(xPlot)
    # pd3alad=plot(hcat(dLogalad.Tactual[:,convIt],dLogalad.Tpwl[:,convIt])*1000,label=["Actual Temp" "PWL Temp"],xlims=(0,hubS.K),xlabel="Time",ylabel="Temp (K)")
    # plot!(pd3alad,1:horzLen+1,hubS.Tmax*ones(hubS.K)*1000,label="XFRM Limit",line=(:dash,:red))
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
    # plot!(uSumPlotalad,cSol.uSum,seriescolor=:black,linewidth=2,linealpha=0.8)
    #
    # zSumPlotalad=plot(dLogalad.zSum[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Z sum",xlims=(0,horzLen+1),legend=false)
    # plot!(zSumPlotalad,cSol.zSum,seriescolor=:black,linewidth=2,linealpha=0.8)
    #
    # constPlotalad2=plot(dLogalad.couplConst[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="curr constraint diff",xlims=(0,horzLen+1),legend=false)
    #
    # lamPlotalad=plot(dLogalad.Lam[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Lambda",xlims=(0,horzLen+1),legend=false)
    # plot!(lamPlotalad,cSol.lamCoupl,seriescolor=:black,linewidth=2,linealpha=0.8)
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
    # setChangesPlot=plot(2:convIt,setChanges[2:convIt],xlabel="Iteration",ylabel="Changes in Active inequality constraints",
    #                   legend=false,xlims=(2,convIt))
    # #solChangesplot=plot(2:convIt,hcat(ΔY[2:convIt],convCheck[2:convIt]),xlabel="Iteration",labels=["ΔY" "y-x"],xlims=(1,convIt))
    #
    # fPlotalad=plot(dCMalad.obj[1:convIt,1],xlabel="Iteration",ylabel="obj function gap",xlims=(1,convIt),legend=false,yscale=:log10)
    # convItPlotalad=plot(dCMalad.lamIt[1:convIt,1],xlabel="Iteration",ylabel="2-Norm Lambda Gap",xlims=(1,convIt),legend=false) #,yscale=:log10
    # convPlotalad=plot(dCMalad.lam[1:convIt,1],xlabel="Iteration",ylabel="central lambda gap",xlims=(1,convIt),legend=false,yscale=:log10)
    # constPlotalad=plot(dCMalad.couplConst[1:convIt,1],xlabel="Iteration",ylabel="curr constraint Gap",xlims=(1,convIt),legend=false,yscale=:log10)

    #save current state and update for next timeSteps
    dSol.Tpwl[stepI,1]=dLog.Tpwl[1,1,convIt]
    dSol.U[stepI,:]=dLog.U[1,1:H,convIt]
    dSol.E[stepI,:]=dLog.E[1,1:H,convIt]
    #dSol.Z[stepI,:]=dLogalad.Z[1:S,convIt]
    dSol.uSum[stepI,1]=dLogalad.uSum[1,1,convIt]
    dSol.zSum[stepI,1]=dLogalad.zSum[1,1,convIt]
    dSol.Itotal[stepI,1]=dLog.Itotal[1,1,convIt]
    dSol.Tactual[stepI,1]=dLog.Tactual[1,1,convIt]
    dSol.convIt[stepI,1]=convIt

    # new states
    global t0=dSol.Tactual[stepI,1]
    global s0=dSol.Sn[stepI,:]

    #function getAttr()
    #clean this up
    if convIt==1
        dSol.lamCoupl[stepI,1]=prevLam[1,1]
        if stepI+horzLen==hubS.K
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
        dSol.lamCoupl[stepI,1]=dLogalad.Lam[1,1,convIt-1]
        if stepI+horzLen==hubS.K
            newLam=dLogalad.Lam[2:horzLen+1,1,convIt-1]
            newVu=dLogalad.Vu[2:horzLen+1,:,convIt-1]
            newVz=dLogalad.Vz[2:horzLen+1,:,convIt-1]
            newVt=dLogalad.Vt[2:horzLen+1,1,convIt-1]
            newVe=dLogalad.Ve[2:horzLen+1,:,convIt-1]
            newVd=dLogalad.Vd[2:horzLen+1,:,convIt-1]
        else
            newLam=vcat(dLogalad.Lam[2:horzLen+1,:,convIt-1],dLogalad.Lam[2:horzLen+1,:,convIt-1])
            newVu=vcat(dLogalad.Vu[2:horzLen+1,:,convIt-1],dLogalad.Vu[2:horzLen+1,:,convIt-1])
            newVz=vcat(dLogalad.Vz[2:horzLen+1,:,convIt-1],dLogalad.Vz[2:horzLen+1,:,convIt-1])
            newVt=vcat(dLogalad.Vt[2:horzLen+1,:,convIt-1],dLogalad.Vt[2:horzLen+1,:,convIt-1])
            newVe=vcat(dLogalad.Ve[2:horzLen+1,:,convIt-1],dLogalad.Ve[2:horzLen+1,:,convIt-1])
            newVd=vcat(dLogalad.Vd[2:horzLen+1,:,convIt-1],dLogalad.Vd[2:horzLen+1,:,convIt-1])
        end
    end

    global prevLam=newLam
    global prevVu=newVu
    global prevVz=newVz
    global prevVt=newVt
    global prevVs=newVs

    return nothing
end

function hubALAD(maxIt::Int,hubS::scenarioHubStruct,cSol::hubSolutionStruct,mode::Bool,eqForm::Bool,silent::Bool)
    H=hubS.H
    K=hubS.K
    Nh=hubS.Nh
    dSol=hubSolutionStruct(K=K,H=H)

    for stepI=1:K
        @printf "%s: time step %g of %g...." Dates.format(Dates.now(),"HH:MM:SS") stepI K
        try
            runHubALADStep(stepI,maxIt,hubS,dSol,cSol,mode,eqForm,silent)
            @printf "convIt: %g\n" dSol.convIt[stepI,1]
        catch e
            @printf "error: %s" e
            break
        end
    end

    objFun(e,u)=sum(sum((e[k,h]-hubS.eMax[k,h])^2*hubS.Qh[1,h] for h=1:H) +sum((u[k,h])^2*hubS.Rh[1,h] for h=1:H) for k=1:K)
    dSol.objVal[1,1]=objFun(dSol.E,dSol.U)

    return dSol
end
