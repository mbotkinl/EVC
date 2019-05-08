# hub functions


#Central
function runHubCentralStep(stepI,hubS,cSol,mode,silent)

    global t0
    global e0
    H=hubS.H
    K=hubS.K
    horzLen=min(hubS.K1,K-stepI)

    eMax=hubS.eMax[stepI:(stepI+horzLen),:]
    uMax=hubS.uMax[stepI:(stepI+horzLen),:]
    eDepart_min=hubS.eDepart_min[stepI:(stepI+horzLen),:]
    eArrive_pred=hubS.eArrive_pred[stepI:(stepI+horzLen),:]
    eArrive_actual=hubS.eArrive_actual[stepI:(stepI+horzLen),:]
    slackMax=hubS.slackMax[stepI:(stepI+horzLen),:]

    #cModel = Model(solver = IpoptSolver())
    #cModel = Model(solver = GurobiSolver(OutputFlag=0,QCPDual=1))
    cModel = Model(solver = GurobiSolver(NumericFocus=3))

    @variable(cModel,u[1:(horzLen+1),1:H])
    @variable(cModel,t[1:(horzLen+1)])
    @variable(cModel,e[1:(horzLen+1),1:H])
    @variable(cModel,eΔ[1:(horzLen+1),1:H])

    #objective
    @objective(cModel,Min,sum(sum(hubS.Qh[h]*(e[k,h]-eMax[k,h])^2+hubS.Rh[h]*u[k,h]^2-hubS.Oh[h]*eΔ[k,h] for k=1:horzLen+1) for h=1:H))

    #transformer constraints
    if coordinated==true
        @constraint(cModel,upperTCon,t.<=hubS.Tmax)
    end
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
		if coordinated==true
    		@constraint(cModel,z.<=hubS.deltaI)
		end
    	#coupling constraint
    	@constraint(cModel,currCon[k=1:horzLen+1],0==-sum(u[k,h] for h=1:H)-hubS.iD_pred[stepI+(k-1)]+sum(z[k,s] for s=1:hubS.S))
    end

    #hub constraints
    @constraint(cModel,stateCon1[h=1:H],e[1,h]==e0[h]+hubS.ηP[h]*u[1,h]-(eDepart_min[1,h]+eΔ[1,h])+eArrive_pred[1,h])
    @constraint(cModel,stateCon[k=1:horzLen,h=1:H],e[k+1,h]==e[k,h]+hubS.ηP[h]*u[k+1,h]-(eDepart_min[k+1,h]+eΔ[k+1,h])+eArrive_pred[k+1,h])
    @constraint(cModel,e.>=0)
    @constraint(cModel,eMaxCon[k=1:horzLen+1,h=1:H],e[k,h]<=eMax[k,h])
    #@constraint(cModel,uMaxCon,u.<=uMax)
    @constraint(cModel,u.>=0)
    @constraint(cModel,eΔ.<=slackMax)
    @constraint(cModel,eΔ.>=0)

    if solverSilent
        @suppress_out begin
			status = solve(cModel)
        end
    else
		status = solve(cModel)
    end

    @assert status==:Optimal "Central Hub optimization not solved to optimality"

    eRaw=getvalue(e)
    uRaw=getvalue(u)
    tRaw=getvalue(t)
    lambdaCurr=-getdual(currCon)
    extraE=getvalue(eΔ)
	cSol.Tpred[stepI,1]=tRaw[1,1]
    if mode=="PWL"
		zRaw=getvalue(z)
        cSol.zSum[stepI,1]=sum(zRaw[1,:])
    end

    #apply current and actual departures and arrivals
    nextU=uRaw[1,:]
    cSol.U[stepI,:]=nextU
    cSol.uSum[stepI,1]=sum(nextU)
    cSol.Lam[stepI,1]=lambdaCurr[1,1]
	cSol.D[stepI,:]=extraE[1,:]
    cSol.E_depart[stepI,:]=eDepart_min[1,:]+extraE[1,:]
    cSol.E_arrive[stepI,:]=eArrive_actual[1,:]
    cSol.E[stepI,:]=e0[:]+hubS.ηP[:].*nextU-(cSol.E_depart[stepI,:])+cSol.E_arrive[stepI,:]
    cSol.Tactual[stepI,1]=hubS.τP*t0+hubS.γP*(cSol.uSum[stepI,1]+hubS.iD_actual[stepI,1])^2+hubS.ρP*hubS.Tamb[stepI,1]
	cSol.timeSolve[stepI,1]=getsolvetime(cModel)

	# p1nl=plot(eRaw,xlabel="Time",ylabel="Energy",xlims=(1,horzLen+1),label="hub energy")
    # plot!(p1nl,eMax,label="hub max")
    # plot!(p1nl,eMax*.8,label="minimum departure")
	#
    # p2nl=plot(uRaw,xlabel="Time",ylabel="Hub Current (kA)",legend=false,xlims=(1,horzLen+1))
    # plot!(p2nl,uMax)
	#
    # p3nl=plot(1:horzLen+1,tRaw,label="XFRM Temp",xlims=(1,horzLen+1),xlabel="Time",ylabel="Temp (C)")
    # plot!(p3nl,1:horzLen+1,hubS.Tmax*ones(horzLen+1),label="XFRM Limit",line=(:dash,:red))
	#
    # p4nl=plot(1:horzLen+1,lambdaCurr,label="Time",ylabel=raw"Lambda ($/kA)",xlims=(1,horzLen+1),legend=false)
	#
    # plot(sum(eMax,dims=2),label="max")
    # plot!(sum(eRaw,dims=2),label="e")
	# h=1
    # plot(eDepart_min[:,h],label="min",xlims=(200,horzLen+1))
    # plot!(eDepart_min[:,h]+extraE[:,h],label="actual")
    # plot!(eDepart_min[:,h]+slackMax[:,h],label="max")
	#
	# aggU=plot(hcat(sum(uRaw,dims=2),sum(uRaw,dims=2) .+ hubS.iD_actual[1:horzLen+1],hubS.iD_actual[1:horzLen+1]),label=["Central" "Total" "iD"],xlims=(0,horzLen+1),xlabel="Time",ylabel="Current (kA)")
	# plot!(aggU,1:hubS.K,hubS.ItotalMax*ones(hubS.K),label="XFRM Limit",line=(:dash,:red))

    # new states
    t0=cSol.Tactual[stepI,1]
    e0=cSol.E[stepI,:]

    return nothing
end

function hubCentral(hubS::scenarioHubStruct,mode::String,silent::Bool)
    H=hubS.H
    K=hubS.K
    cSol=hubSolutionStruct(K=K,H=H)

    Juno.progress() do id
        for stepI=1:K
            @info "$(Dates.format(Dates.now(),"HH:MM:SS")): $(stepI) of $(K)....\n" progress=stepI/K _id=id
            @printf "%s: time step %g of %g....\n" Dates.format(Dates.now(),"HH:MM:SS") stepI K
            try
                cSol.timeT[stepI]=@elapsed runHubCentralStep(stepI,hubS,cSol,mode,silent)
            catch e
                @printf "error: %s" e
                break
            end
        end
    end

	objFun(e,u,d)=sum(sum((e[k,h]-hubS.eMax[k,h])^2*hubS.Qh[1,h]+(u[k,h])^2*hubS.Rh[1,h]-hubS.Oh[h]*d[k,h] for h=1:H) for k=1:K)
    cSol.objVal[1,1]=objFun(cSol.E,cSol.U,cSol.D)

    return cSol
end


#Dual Ascent
function localEVDual(hubInd::Int,p::Int,stepI::Int,hubS::scenarioHubStruct,dLog::hubItLogPWL,itLam,e0,solverSilent)
    horzLen=min(hubS.K1,hubS.K-stepI)

    eMax=hubS.eMax[stepI:(stepI+horzLen),hubInd]
    uMax=hubS.uMax[stepI:(stepI+horzLen),hubInd]
    eDepart_min=hubS.eDepart_min[stepI:(stepI+horzLen),hubInd]
    eArrive_pred=hubS.eArrive_pred[stepI:(stepI+horzLen),hubInd]
    eArrive_actual=hubS.eArrive_actual[stepI:(stepI+horzLen),hubInd]
    slackMax=hubS.slackMax[stepI:(stepI+horzLen),hubInd]

    hubM=Model(solver = GurobiSolver(NumericFocus=3))
    @variable(hubM,u[1:horzLen+1])
    @variable(hubM,e[1:horzLen+1])
    @variable(hubM,eΔ[1:(horzLen+1)])
    objExp=sum((e[k,1]-eMax[k,1])^2*hubS.Qh[1,hubInd]+(u[k,1])^2*hubS.Rh[1,hubInd]-hubS.Oh[1,hubInd]*eΔ[k,1]+itLam[k,1]*u[k,1] for k=1:horzLen+1)
    @objective(hubM,Min,objExp)

    #hub constraints
    @constraint(hubM,stateCon1,e[1,1]==e0[hubInd]+hubS.ηP[1,hubInd]*u[1]-(eDepart_min[1]+eΔ[1])+eArrive_pred[1])
    @constraint(hubM,stateCon[k=1:horzLen],e[k+1,1]==e[k,1]+hubS.ηP[1,hubInd]*u[k+1]-(eDepart_min[k+1]+eΔ[k+1])+eArrive_pred[k+1])
    @constraint(hubM,e.>=0)
    @constraint(hubM,eMaxCon,e.<=eMax)
    #@constraint(hubM,uMaxCon,u.<=uMax)
    @constraint(hubM,u.>=0)
    @constraint(hubM,eΔ.<=slackMax)
    @constraint(hubM,eΔ.>=0)

	if solverSilent
        @suppress_out begin
			statusM = solve(hubM)
        end
    else
		statusM = solve(hubM)
    end

    @assert statusM==:Optimal "Hub optimization not solved to optimality"

    dLog.U[:,hubInd,p]=round.(getvalue(u),sigdigits=roundSigFigs)
    dLog.E[:,hubInd,p]=round.(getvalue(e),sigdigits=roundSigFigs)
	dLog.D[:,hubInd,p]=round.(getvalue(eΔ),sigdigits=roundSigFigs)
	dLog.timeSolve[1,p]=max(getsolvetime(hubM),dLog.timeSolve[1,p])

    return nothing
end

function localXFRMDual(p::Int,stepI::Int,hubS::scenarioHubStruct,dLog::hubItLogPWL,mode::String,itLam,t0,solverSilent)
    S=hubS.S
    horzLen=min(hubS.K1,hubS.K-stepI)

    coorM=Model(solver = GurobiSolver(NumericFocus=3))
    @variable(coorM,t[1:(horzLen+1)])
    @constraint(coorM,upperTCon,t.<=hubS.Tmax)
    @constraint(coorM,t.>=0)
    if mode=="NL"
        @variable(coorM,itotal[1:(horzLen+1)])
        @objective(coorM,Min,sum(itLam[k,1]*itotal[k] for k=1:(horzLen+1)))
        @constraint(coorM,itotal.<=hubS.ItotalMax)
        @constraint(coorM,itotal.>=0)
        @NLconstraint(coorM,tempCon1,t[1,1]==hubS.τP*t0+hubS.γP*(itotal[1])^2+hubS.ρP*hubS.Tamb[stepI,1])
        @NLconstraint(coorM,tempCon2[k=1:horzLen],t[k+1,1]==hubS.τP*t[k,1]+hubS.γP*(itotal[k+1])^2+hubS.ρP*hubS.Tamb[stepI+k,1])
    elseif mode=="relax1"
        @variable(coorM,itotal[1:(horzLen+1)])
        @objective(coorM,Min,sum(itLam[k,1]*itotal[k] for k=1:(horzLen+1)))
        @constraint(coorM,itotal.<=hubS.ItotalMax)
        @constraint(coorM,itotal.>=0)
        @constraint(coorM,tempCon1,t[1,1]>=hubS.τP*t0+hubS.γP*(itotal[1])^2+hubS.ρP*hubS.Tamb[stepI,1])
        @constraint(coorM,tempCon2[k=1:horzLen],t[k+1,1]>=hubS.τP*t[k,1]+hubS.γP*(itotal[k+1])^2+hubS.ρP*hubS.Tamb[stepI+k,1])
    elseif mode=="PWL"
        @variable(coorM,z[1:(horzLen+1),1:hubS.S])
        @objective(coorM,Min,sum(itLam[k,1]*sum(-z[k,s] for s=1:S) for k=1:(horzLen+1)))
        @constraint(coorM,tempCon1,t[1,1]==hubS.τP*t0+hubS.γP*hubS.deltaI*sum((2*s-1)*z[1,s] for s=1:hubS.S)+hubS.ρP*hubS.Tamb[stepI,1])
        @constraint(coorM,tempCon2[k=1:horzLen],t[k+1,1]==hubS.τP*t[k,1]+hubS.γP*hubS.deltaI*sum((2*s-1)*z[k+1,s] for s=1:hubS.S)+hubS.ρP*hubS.Tamb[stepI+k,1])
        @constraint(coorM,z.>=0)
        @constraint(coorM,z.<=hubS.deltaI)
    end

	if solverSilent
        @suppress_out begin
			statusC = solve(coorM)
        end
    else
		statusC = solve(coorM)
    end

    @assert statusC==:Optimal "Dual Ascent XFRM optimization not solved to optimality"

    dLog.Tpred[:,1,p]=round.(getvalue(t),sigdigits=roundSigFigs)
    dLog.Z[:,:,p]=round.(getvalue(z),sigdigits=roundSigFigs)
	dLog.timeSolve[1,p]=max(getsolvetime(coorM),dLog.timeSolve[1,p])

    return nothing
end

function runHubDualIt(p,stepI,hubS,itLam,dLog,dSol,cSol,mode,silent)

    global t0
    global e0
    #global prevLam

    H=hubS.H
    K=hubS.K
    S=hubS.S
    horzLen=min(hubS.K1,K-stepI)

    minAlpha=1e-9
	alphaDivRate=20


    #solve subproblem for each EV
	if runParallel
		@sync @distributed for hubInd=1:H
			localEVDual(hubInd,p,stepI,hubS,dLog,itLam,e0,solverSilent)
		end
	else
		for hubInd=1:H
			localEVDual(hubInd,p,stepI,hubS,dLog,itLam,e0,solverSilent)
		end
	end

    localXFRMDual(p,stepI,hubS,dLog,mode,itLam,t0,solverSilent)

    #grad of lagragian
    for k=1:horzLen+1
        dLog.uSum[k,1,p]=sum(dLog.U[k,h,p] for h=1:H)
        dLog.zSum[k,1,p]=sum(dLog.Z[k,s,p] for s=1:S)
        dLog.couplConst[k,1,p]=dLog.uSum[k,1,p]+ hubS.iD_pred[stepI+(k-1),1] - dLog.zSum[k,1,p]
    end
    #dCM.couplConst[p,1]=norm(dLog.couplConst[:,p],2)

    #update lambda
    alphaP= max(alpha0/ceil(p/alphaDivRate),minAlpha)
    dLog.itUpdate[1,1,p]=alphaP
    #alphaP= alphaP*alphaRate

	# dLog.Lam[:,1,p]=round.(max.(itLam[:,1]+alphaP*dLog.couplConst[:,1,p],0),sigdigits=roundSigFigs)
    dLog.Lam[:,1,p]=round.(itLam[:,1]+alphaP*dLog.couplConst[:,1,p],sigdigits=roundSigFigs)

    #calculate actual temperature from nonlinear model of XFRM
    for k=1:horzLen+1
        dLog.Itotal[k,1,p]=dLog.uSum[k,1,p] + hubS.iD_actual[stepI+(k-1),1]
    end
    dLog.Tactual[1,1,p]=hubS.τP*t0+hubS.γP*dLog.Itotal[1,1,p]^2+hubS.ρP*hubS.Tamb[stepI,1]
    for k=1:horzLen
        dLog.Tactual[k+1,1,p]=hubS.τP*dLog.Tactual[k,1,p]+hubS.γP*dLog.Itotal[k+1,1,p]^2+hubS.ρP*hubS.Tamb[stepI+k,1]
    end

    #check convergence
	constGap=norm(dLog.couplConst[:,1,p],1)
	itGap = norm(dLog.Lam[:,1,p]-itLam[:,1],2)
    #if((constGap<=primChk) && (itGap <= dualChk))
	if(constGap<=primChk)

        if !silent @printf "Converged after %g iterations\n" p end
		@printf "Y"
        convIt=p
        return true
    else
        if !silent
            #@printf "itGap    %e after %g iterations\n" itGap p
            @printf "constGap %e after %g iterations\n" constGap p
            #@printf "snGap   %e after %g iterations\n" snGap p
            #@printf "unGap   %e after %g iterations\n" unGap p
            #@printf("fGap    %e after %g iterations\n\n",fGap,p)
        end
        #prevLam=dLog.Lam[:,1,p]
        return false
    end
end

function runHubDualStep(stepI,maxIt,hubS,dSol,cSol,mode,silent)
    K=hubS.K
    S=hubS.S
    horzLen=min(hubS.K1,K-stepI)
    dLog=hubItLogPWL(horzLen=horzLen,H=hubS.H,S=hubS.S)
	p=1
    timeStart=now()
    while (p<=maxIt && round(now()-timeStart,Second)<=Dates.Second(9/10*hubS.Ts))
		# global p
		if p==1
			itLam=prevLam
			#alpha0 =1 #for kA
			global alpha0=max(min(maximum(prevLam)/100,1e6),1e-3)
		else
			itLam=round.(dLog.Lam[:,1,(p-1)],digits=8)
		end
        cFlag=runHubDualIt(p,stepI,hubS,itLam,dLog,dSol,cSol,mode,silent)
        global convIt=p
        if cFlag
            break
        end
		p+=1
    end

    # plot(dLog.U[:,:,convIt])
    # plot(dLog.E[:,:,convIt])
    # pd3alad=plot(hcat(dLog.Tactual[:,1,convIt],dLog.Tpred[:,1,convIt])*1000,label=["Actual Temp" "PWL Temp"],xlims=(0,hubS.K),xlabel="Time",ylabel="Temp (K)")
    # plot!(pd3alad,1:horzLen+1,hubS.Tmax*ones(hubS.K)*1000,label="XFRM Limit",line=(:dash,:red))
	#
    # #convergence plots
    # halfCI=Int(floor(convIt/2))
    # CList=reshape([range(colorant"blue", stop=colorant"yellow",length=halfCI);
    #                range(colorant"yellow", stop=colorant"red",length=convIt-halfCI)], 1, convIt);
	#
    # uSumPlotd=plot(dLog.uSum[:,1,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Current Sum (kA)",xlims=(0,horzLen+1),legend=false)
    # plot!(uSumPlotd,sum(uRaw,dims=2),seriescolor=:black,linewidth=2,linealpha=0.8)
	# plot!(uSumPlotd,cSol.uSum[stepI:stepI+horzLen],seriescolor=:black,linewidth=2,linealpha=0.8)
	#
    # zSumPlotd=plot(dLog.zSum[:,1,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Z sum",xlims=(0,horzLen+1),legend=false)
	# plot!(zSumPlotd,sum(zRaw,dims=2),seriescolor=:black,linewidth=2,linealpha=0.8)
	#
    # lamPlotd=plot(dLog.Lam[:,1,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Lambda",xlims=(0,horzLen+1),legend=false)
    # plot!(lamPlotd,lambdaCurr,seriescolor=:black,linewidth=2,linealpha=0.8)
	# plot!(lamPlotd,cSol.Lam[stepI:stepI+horzLen],seriescolor=:black,linewidth=2,linealpha=0.8)
    # #plot(dLog.Lam[:,1:convIt],color=:RdYlBu,line_z=(1:convIt)')
	#
    # constPlot2=plot(dLog.couplConst[:,1,1:convIt],seriescolor=CList,xlabel="Time",ylabel="curr constraint diff",xlims=(0,horzLen+1),legend=false)
	#
    # #convergence metric plots
    # fPlot=plot(1:convIt,dCM.obj[1:convIt,1],xlabel="Iteration",ylabel="obj function gap",xlims=(1,convIt),legend=false,yscale=:log10)
    # convItPlot=plot(1:convIt,dCM.lamIt[1:convIt,1],xlabel="Iteration",ylabel="2-Norm Lambda Gap",xlims=(1,convIt),legend=false,yscale=:log10)
    # convPlot=plot(1:convIt,dCM.lam[1:convIt,1],xlabel="Iteration",ylabel="central lambda gap",xlims=(2,convIt),legend=false,yscale=:log10)
    # constPlot=plot(1:convIt,dCM.couplConst[1:convIt,1],xlabel="Iteration",ylabel="curr constraint Gap",xlims=(2,convIt),legend=false,yscale=:log10)


    #save current state and update for next timeSteps
    dSol.Tpred[stepI,1]=dLog.Tpred[1,1,convIt]
    dSol.U[stepI,:]=dLog.U[1,1:H,convIt]
	dSol.uSum[stepI,1]=sum(dLog.U[1,1:H,convIt])
    dSol.E[stepI,:]=dLog.E[1,1:H,convIt]
	dSol.D[stepI,:]=dLog.D[1,1:H,convIt]
	dSol.zSum[stepI,1]=dLog.zSum[1,1,convIt]
    dSol.Itotal[stepI,1]=dLog.Itotal[1,1,convIt]
    dSol.Tactual[stepI,1]=dLog.Tactual[1,1,convIt]
    dSol.convIt[stepI,1]=convIt
	dSol.timeSolve[stepI,1]=mean(dLog.timeSolve[1,1:convIt])

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
        dSol.Lam[stepI,1]=dLog.Lam[1,1,convIt]
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

    Juno.progress() do id
        for stepI=1:K
            @info "$(Dates.format(Dates.now(),"HH:MM:SS")): $(stepI) of $(K)....\n" progress=stepI/K _id=id
            @printf "%s: time step %g of %g...." Dates.format(Dates.now(),"HH:MM:SS") stepI K
            try
                dSol.timeT[stepI]=@elapsed runHubDualStep(stepI,maxIt,hubS,dSol,cSol,mode,silent)
                @printf "convIt: %g\n" dSol.convIt[stepI,1]
            catch e
                @printf "error: %s" e
                break
            end
        end
    end

	objFun(e,u,d)=sum(sum((e[k,h]-hubS.eMax[k,h])^2*hubS.Qh[1,h]+(u[k,h])^2*hubS.Rh[1,h]-hubS.Oh[h]*d[k,h] for h=1:H) for k=1:K)
    dSol.objVal[1,1]=objFun(dSol.E,dSol.U,dSol.D)

    return dSol
end

#ADMM
function localEVADMM(hubInd::Int,p::Int,stepI::Int,hubS::scenarioHubStruct,dLogadmm::hubItLogPWL,itLam,itVu,e0,itρ,solverSilent)

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
    objExp=sum((e[k,1]-eMax[k,1])^2*hubS.Qh[1,hubInd]+(u[k,1])^2*hubS.Rh[1,hubInd]-hubS.Oh[1,hubInd]*eΔ[k,1]+
                            itLam[k,1]*(u[k,1]-itVu[k,hubInd])+itρ/2*(u[k,1]-itVu[k,hubInd])^2 for k=1:horzLen+1)

    @objective(hubM,Min,objExp)

    #hub constraints
    @constraint(hubM,stateCon1,e[1,1]==e0[hubInd]+hubS.ηP[1,hubInd]*u[1]-(eDepart_min[1]+eΔ[1])+eArrive_pred[1])
    @constraint(hubM,stateCon[k=1:horzLen],e[k+1,1]==e[k,1]+hubS.ηP[1,hubInd]*u[k+1]-(eDepart_min[k+1]+eΔ[k+1])+eArrive_pred[k+1])
    @constraint(hubM,e.>=0)
    @constraint(hubM,eMaxCon,e.<=eMax)
    #@constraint(hubM,uMaxCon,u.<=uMax)
    @constraint(hubM,u.>=0)
    @constraint(hubM,eΔ.<=slackMax)
    @constraint(hubM,eΔ.>=0)

	if solverSilent
        @suppress_out begin
			statusM = solve(hubM)
        end
    else
		statusM = solve(hubM)
    end

    @assert statusM==:Optimal "Hub optimization not solved to optimality"

    uVal=getvalue(u)
    eVal=getvalue(e)
    eΔVal=getvalue(eΔ)

    dLogadmm.E[:,hubInd,p]=round.(eVal,sigdigits=roundSigFigs)
    dLogadmm.U[:,hubInd,p]=round.(uVal,sigdigits=roundSigFigs)
	dLogadmm.D[:,hubInd,p]=round.(eΔVal,sigdigits=roundSigFigs)
	dLogadmm.timeSolve[1,p]=max(getsolvetime(hubM),dLogadmm.timeSolve[1,p])

    return nothing
end

function localXFRMADMM(p::Int,stepI::Int,hubS::scenarioHubStruct,dLogadmm::hubItLogPWL,mode::String,
	t0,itLam,itVz,itρ,solverSilent)

    S=hubS.S
    horzLen=min(hubS.K1,hubS.K-stepI)

    #N+1 decoupled problem aka transformer current
    coorM=Model(solver = GurobiSolver())
    @variable(coorM,t[1:(horzLen+1)])
    @constraint(coorM,upperTCon,t.<=hubS.Tmax)
    @constraint(coorM,t.>=0)
    if mode=="NL"
        # @variable(coorM,itotal[1:(horzLen+1)])
        # @objective(coorM,Min,sum(itLam[k,1]*itotal[k] for k=1:(horzLen+1)))
        # @constraint(coorM,itotal.<=hubS.ItotalMax)
        # @constraint(coorM,itotal.>=0)
        # @NLconstraint(coorM,tempCon1,t[1,1]==hubS.τP*t0+hubS.γP*(itotal[1])^2+hubS.ρP*hubS.Tamb[stepI,1])
        # @NLconstraint(coorM,tempCon2[k=1:horzLen],t[k+1,1]==hubS.τP*t[k,1]+hubS.γP*(itotal[k+1])^2+hubS.ρP*hubS.Tamb[stepI+k,1])
    elseif mode=="relax1"
        # @variable(coorM,itotal[1:(horzLen+1)])
        # @objective(coorM,Min,sum(itLam[k,1]*itotal[k] for k=1:(horzLen+1)))
        # @constraint(coorM,itotal.<=hubS.ItotalMax)
        # @constraint(coorM,itotal.>=0)
        # @constraint(coorM,tempCon1,t[1,1]>=hubS.τP*t0+hubS.γP*(itotal[1])^2+hubS.ρP*hubS.Tamb[stepI,1])
        # @constraint(coorM,tempCon2[k=1:horzLen],t[k+1,1]>=hubS.τP*t[k,1]+hubS.γP*(itotal[k+1])^2+hubS.ρP*hubS.Tamb[stepI+k,1])
    elseif mode=="PWL"
        @variable(coorM,z[1:(horzLen+1),1:hubS.S])
        objExp=sum(itLam[k,1]*(sum(-z[k,s] for s=1:S)-sum(itVz[k,s] for s=1:S))+
                  itρ/2*(sum(-z[k,s] for s=1:S)-sum(itVz[k,s] for s=1:S))^2 for k=1:(horzLen+1))
        @objective(coorM,Min, objExp)
        @constraint(coorM,tempCon1,t[1,1]==hubS.τP*t0+hubS.γP*hubS.deltaI*sum((2*s-1)*z[1,s] for s=1:hubS.S)+hubS.ρP*hubS.Tamb[stepI,1])
        @constraint(coorM,tempCon2[k=1:horzLen],t[k+1,1]==hubS.τP*t[k,1]+hubS.γP*hubS.deltaI*sum((2*s-1)*z[k+1,s] for s=1:hubS.S)+hubS.ρP*hubS.Tamb[stepI+k,1])
        @constraint(coorM,z.>=0)
        @constraint(coorM,z.<=hubS.deltaI)
    end

	if solverSilent
        @suppress_out begin
			statusC = solve(coorM)
        end
    else
		statusC = solve(coorM)
    end

    @assert statusC==:Optimal "Dual Ascent XFRM optimization not solved to optimality"

    zVal=getvalue(z)
    tVal=getvalue(t)

    dLogadmm.Tpred[:,:,p]=round.(tVal,sigdigits=roundSigFigs)
    dLogadmm.Z[:,:,p]=round.(zVal,sigdigits=roundSigFigs)
	dLogadmm.timeSolve[1,p]=max(getsolvetime(coorM),dLogadmm.timeSolve[1,p])
    return nothing
end

function runEVADMMIt(p,stepI,hubS,itLam,itVu,itVz,itρ,dLogadmm,dSol,cSol,mode,eqForm,silent)
    H=hubS.H
    K=hubS.K
    S=hubS.S
    horzLen=min(hubS.K1,K-stepI)

    #initialize with current states
    global e0
    global t0

    #solve decoupled
    if runParallel
        @sync @distributed for hubInd=1:H
			localEVADMM(hubInd,p,stepI,hubS,dLogadmm,itLam,itVu,e0,itρ,solverSilent)
		end
    else
        for hubInd=1:H
			localEVADMM(hubInd,p,stepI,hubS,dLogadmm,itLam,itVu,e0,itρ,solverSilent)
		end
    end

	localXFRMADMM(p,stepI,hubS,dLogadmm,mode,t0,itLam,itVz,itρ,solverSilent)

	#lambda update eq 7.68
    for k=1:horzLen+1
		dLogadmm.uSum[k,1,p]=sum(dLogadmm.U[k,h,p] for h=1:H)
		dLogadmm.zSum[k,1,p]=sum(dLogadmm.Z[k,s,p] for s=1:S)
		dLogadmm.couplConst[k,1,p]=dLogadmm.uSum[k,1,p] + hubS.iD_pred[stepI+(k-1),1] - dLogadmm.zSum[k,1,p]
        dLogadmm.Lam[k,1,p]=itLam[k,1]+itρ/(max(horzLen+1,H))*(dLogadmm.couplConst[k,1,p])
    end

    #calculate actual temperature from nonlinear model of XFRM
    for k=1:horzLen+1
        dLogadmm.Itotal[k,1,p]=dLogadmm.uSum[k,1,p] + hubS.iD_actual[stepI+(k-1),1]
    end
    dLogadmm.Tactual[1,1,p]=hubS.τP*t0+hubS.γP*dLogadmm.Itotal[1,1,p]^2+hubS.ρP*hubS.Tamb[stepI,1] #fix for mpc
    for k=1:horzLen
        dLogadmm.Tactual[k+1,1,p]=hubS.τP*dLogadmm.Tactual[k,1,p]+hubS.γP*dLogadmm.Itotal[k,1,p]^2+hubS.ρP*hubS.Tamb[stepI+k,1]  #fix for mpc
    end

	#v upate eq 7.67
	for k=1:horzLen+1
		dLogadmm.Vu[k,:,p]=dLogadmm.U[k,:,p].+(itLam[k,1]-dLogadmm.Lam[k,1,p])/(itρ/1)
		dLogadmm.Vz[k,:,p]=-dLogadmm.Z[k,:,p].+(itLam[k,1]-dLogadmm.Lam[k,1,p])/(itρ/1)
	end

	#update rho
	#ρADMMp = ρADMM/ceil(p/ρDivRate)
	dLogadmm.itUpdate[1,1,p]= min(itρ*ρDivRate,maxRho)


    #check for convergence
	constGap=norm(dLogadmm.couplConst[:,1,p],1)
	itGap = norm(dLogadmm.Lam[:,1,p]-itLam[:,1],2)
	#cc=norm(hcat((itVu[:,:]-dLogadmm.Vu[:,:,p]),(itVz[:,:]-dLogadmm.Vz[:,:,p])),2)
	# if(constGap <= primChk  && itGap <=dualChk)
	if(constGap <= primChk)
        if !silent @printf "Converged after %g iterations\n" p end
		@printf "Y"
        convIt=p
        #break
        return true
    else
        if !silent
            #@printf "itGap  %e after %g iterations\n" itGap p
            @printf "constGap   %e after %g iterations\n" constGap p
            #@printf "snGap      %e after %g iterations\n" snGap p
            #@printf("fGap       %e after %g iterations\n",fGap,p)
        end
    end

    return false
end

function runHubADMMStep(stepI,maxIt,hubS,dSol,cSol,mode,eqForm,silent)
    K=hubS.K
    S=hubS.S
    horzLen=min(hubS.K1,K-stepI)
    #ogρ=ρALADp #save to reset later
    dLogadmm=hubItLogPWL(horzLen=horzLen,H=hubS.H,S=hubS.S)
	timeStart=now()
	p=1
	while (p<=maxIt && round(now()-timeStart,Second)<=Dates.Second(9/10*hubS.Ts))
		#global p
        #@printf "%git" p
		if p==1
			itLam=prevLam
			itVu=prevVu
			itVz=prevVz
			itρ=ρADMMp
			#itρ=max(min(maximum(prevLam),maxRho),1e-3)
		else
			itLam=round.(dLogadmm.Lam[:,1,(p-1)],sigdigits=roundSigFigs)
			itVu=round.(dLogadmm.Vu[:,:,(p-1)],sigdigits=roundSigFigs)
			itVz=round.(dLogadmm.Vz[:,:,(p-1)],sigdigits=roundSigFigs)
			itρ=round.(dLogadmm.itUpdate[1,1,(p-1)],sigdigits=roundSigFigs)
		end
		#@printf "itρ: %f\n" itρ

        cFlag=runEVADMMIt(p,stepI,hubS,itLam,itVu,itVz,itρ,dLogadmm,dSol,cSol,mode,eqForm,silent)
        global convIt=p
        if cFlag
            break
        end
		p+=1
    end

    # plot(dLogadmm.U[:,:,convIt])
    # plot(dLogadmm.E[:,:,convIt])
    # pd3alad=plot(hcat(dLogadmm.Tactual[:,1,convIt],dLogadmm.Tpred[:,1,convIt])*1000,label=["Actual Temp" "PWL Temp"],xlims=(0,hubS.K),xlabel="Time",ylabel="Temp (K)")
    # plot!(pd3alad,1:horzLen+1,hubS.Tmax*ones(hubS.K)*1000,label="XFRM Limit",line=(:dash,:red))
	# pd4admm=plot(dLogadmm.Lam[:,1,convIt],xlabel="Time",ylabel=raw"Lambda ($/kA)",label="ADMM", xlims=(0,hubS.K))
	# plot!(pd4admm,cSol.Lam,label="Central")
	#
    # #convergence plots
    # halfCI=Int(floor(convIt/2))
    # if halfCI>0
    #     CList=reshape([range(colorant"blue", stop=colorant"yellow",length=halfCI);
    #                    range(colorant"yellow", stop=colorant"red",length=convIt-halfCI)], 1, convIt);
    # else
    #     CList=colorant"red";
    # end
    # uSumPlotalad=plot(dLogadmm.uSum[:,1,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Current Sum (kA)",xlims=(0,horzLen+1),legend=false)
    # #plot!(uSumPlotalad,cSol.uSum,seriescolor=:black,linewidth=2,linealpha=0.8)
	# plot!(uSumPlotalad,sum(uRaw,dims=2),seriescolor=:black,linewidth=2,linealpha=0.8)
	#
    # zSumPlotalad=plot(dLogadmm.zSum[:,1,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Z sum",xlims=(0,horzLen+1),legend=false)
    # #plot!(zSumPlotalad,cSol.zSum,seriescolor=:black,linewidth=2,linealpha=0.8)
	# plot!(zSumPlotalad,sum(zRaw,dims=2),seriescolor=:black,linewidth=2,linealpha=0.8)
	#
    # constPlotalad2=plot(dLogadmm.couplConst[:,1,1:convIt],seriescolor=CList,xlabel="Time",ylabel="curr constraint diff",xlims=(0,horzLen+1),legend=false)
	#
    # lamPlotalad=plot(dLogadmm.Lam[:,1,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Lambda",xlims=(0,horzLen+1),legend=false)
    # #plot!(lamPlotalad,cSol.Lam,seriescolor=:black,linewidth=2,linealpha=0.8)
	# plot!(lamPlotalad,lambdaCurr,seriescolor=:black,linewidth=2,linealpha=0.8)
	#
    # fPlotalad=plot(dCMalad.obj[1:convIt,1],xlabel="Iteration",ylabel="obj function gap",xlims=(1,convIt),legend=false,yscale=:log10)
    # convItPlotalad=plot(dCMalad.lamIt[1:convIt,1],xlabel="Iteration",ylabel="2-Norm Lambda Gap",xlims=(1,convIt),legend=false) #,yscale=:log10
    # convPlotalad=plot(dCMalad.lam[1:convIt,1],xlabel="Iteration",ylabel="central lambda gap",xlims=(1,convIt),legend=false,yscale=:log10)
    # constPlotalad=plot(dCMalad.couplConst[1:convIt,1],xlabel="Iteration",ylabel="curr constraint Gap",xlims=(1,convIt),legend=false,yscale=:log10)

    #save current state and update for next timeSteps
    dSol.Tpred[stepI,1]=dLogadmm.Tpred[1,1,convIt]
    dSol.U[stepI,:]=dLogadmm.U[1,1:H,convIt]
    dSol.E[stepI,:]=dLogadmm.E[1,1:H,convIt]
	dSol.D[stepI,:]=dLogadmm.D[1,1:H,convIt]
    dSol.zSum[stepI,1]=dLogadmm.zSum[1,1,convIt]
    dSol.uSum[stepI,1]=dLogadmm.uSum[1,1,convIt]
    dSol.zSum[stepI,1]=dLogadmm.zSum[1,1,convIt]
    dSol.Itotal[stepI,1]=dLogadmm.Itotal[1,1,convIt]
    dSol.Tactual[stepI,1]=dLogadmm.Tactual[1,1,convIt]
    dSol.convIt[stepI,1]=convIt
	dSol.timeSolve[stepI,1]=mean(dLogadmm.timeSolve[1,1:convIt])

    # new states
    global t0=round.(dSol.Tactual[stepI,1],sigdigits=roundSigFigs)
    global e0=round.(dSol.E[stepI,:],sigdigits=roundSigFigs)

    #function getAttr()
    #clean this up
    if convIt==1
        dSol.Lam[stepI,1]=prevLam[1,1]
        if stepI+horzLen==hubS.K
            newLam=prevLam[2:horzLen+1,1]
            newVu=prevVu[2:horzLen+1,:]
            newVz=prevVz[2:horzLen+1,:]
        else
            newLam=vcat(prevLam[2:horzLen+1,1],prevLam[horzLen+1,1])
            newVu=vcat(prevVu[2:horzLen+1,:],prevVu[horzLen+1,:]')
            newVz=vcat(prevVz[2:horzLen+1,:],prevVz[horzLen+1,:]')
        end
    else
        dSol.Lam[stepI,1]=dLogadmm.Lam[1,1,convIt]
        if stepI+horzLen==hubS.K
            newLam=dLogadmm.Lam[2:horzLen+1,1,convIt-1]
            newVu=dLogadmm.Vu[2:horzLen+1,:,convIt-1]
            newVz=dLogadmm.Vz[2:horzLen+1,:,convIt-1]
        else
            newLam=vcat(dLogadmm.Lam[2:horzLen+1,:,convIt-1],dLogadmm.Lam[horzLen+1,:,convIt-1])
            newVu=vcat(dLogadmm.Vu[2:horzLen+1,:,convIt-1],dLogadmm.Vu[horzLen+1,:,convIt-1]')
            newVz=vcat(dLogadmm.Vz[2:horzLen+1,:,convIt-1],dLogadmm.Vz[horzLen+1,:,convIt-1]')
        end
    end

    global prevLam=round.(newLam,sigdigits=roundSigFigs)
    global prevVu=round.(newVu,sigdigits=roundSigFigs)
    global prevVz=round.(newVz,sigdigits=roundSigFigs)
    #global ρALADp=ogρ

    return nothing
end

function hubADMM(maxIt::Int,hubS::scenarioHubStruct,cSol::hubSolutionStruct,mode::String,silent::Bool)
    H=hubS.H
    K=hubS.K
    dSol=hubSolutionStruct(K=K,H=H)

    Juno.progress() do id
        for stepI=1:K
            @info "$(Dates.format(Dates.now(),"HH:MM:SS")): $(stepI) of $(K)....\n" progress=stepI/K _id=id
            @printf "%s: time step %g of %g...." Dates.format(Dates.now(),"HH:MM:SS") stepI K
            try
                dSol.timeT[stepI]=@elapsed runHubADMMStep(stepI,maxIt,hubS,dSol,cSol,mode,eqForm,silent)
                @printf "convIt: %g\n" dSol.convIt[stepI,1]
            catch e
                @printf "error: %s" e
                break
            end
        end
    end

	objFun(e,u,d)=sum(sum((e[k,h]-hubS.eMax[k,h])^2*hubS.Qh[1,h]+(u[k,h])^2*hubS.Rh[1,h]-hubS.Oh[h]*d[k,h] for h=1:H) for k=1:K)
    dSol.objVal[1,1]=objFun(dSol.E,dSol.U,dSol.D)

    return dSol
end

#ALADIN
function localEVALAD(hubInd::Int,p::Int,stepI::Int,σU::Array{Float64,2},σE::Array{Float64,2},σD::Array{Float64,2},
	hubS::scenarioHubStruct,dLogalad::hubItLogPWL,itLam,itVu,itVe,itVd,e0,itρ,solverSilent)

    horzLen=min(hubS.K1,hubS.K-stepI)

    tolU=1e-6
    tolE=1e-6

    eMax=hubS.eMax[stepI:(stepI+horzLen),hubInd]
    uMax=hubS.uMax[stepI:(stepI+horzLen),hubInd]
    eDepart_min=hubS.eDepart_min[stepI:(stepI+horzLen),hubInd]
    eArrive_pred=hubS.eArrive_pred[stepI:(stepI+horzLen),hubInd]
    eArrive_actual=hubS.eArrive_actual[stepI:(stepI+horzLen),hubInd]
    slackMax=hubS.slackMax[stepI:(stepI+horzLen),hubInd]

    hubM=Model(solver = GurobiSolver(NumericFocus=3))
    @variable(hubM,u[1:horzLen+1])
    @variable(hubM,e[1:horzLen+1])
    @variable(hubM,eΔ[1:(horzLen+1)])
    objExp=sum((e[k,1]-eMax[k,1])^2*hubS.Qh[1,hubInd]+(u[k,1])^2*hubS.Rh[1,hubInd]-hubS.Oh[1,hubInd]*eΔ[k,1]+
                            itLam[k,1]*(u[k,1])+
                            ρALADp/2*(u[k,1]-itVu[k,hubInd])*σU[hubInd,1]*(u[k,1]-itVu[k,hubInd])+
                            ρALADp/2*(eΔ[k,1]-itVd[k,hubInd])*σD[hubInd,1]*(eΔ[k,1]-itVd[k,hubInd])+
                            ρALADp/2*(e[k,1]-itVe[k,hubInd])*σE[hubInd,1]*(e[k,1]-itVe[k,hubInd]) for k=1:horzLen+1)

    @objective(hubM,Min,objExp)

    #hub constraints
    @constraint(hubM,stateCon1,e[1,1]==e0[hubInd]+hubS.ηP[1,hubInd]*u[1]-(eDepart_min[1]+eΔ[1])+eArrive_pred[1])
    @constraint(hubM,stateCon[k=1:horzLen],e[k+1,1]==e[k,1]+hubS.ηP[1,hubInd]*u[k+1]-(eDepart_min[k+1]+eΔ[k+1])+eArrive_pred[k+1])
    @constraint(hubM,e.>=0)
    @constraint(hubM,eMaxCon,e.<=eMax)
    #@constraint(hubM,uMaxCon,u.<=uMax)
    @constraint(hubM,u.>=0)
    @constraint(hubM,eΔ.<=slackMax)
    @constraint(hubM,eΔ.>=0)

	if solverSilent
        @suppress_out begin
			statusM = solve(hubM)
        end
    else
		statusM = solve(hubM)
    end

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
    dLogalad.E[:,hubInd,p]=round.(eVal,sigdigits=roundSigFigs)
    dLogalad.U[:,hubInd,p]=round.(uVal,sigdigits=roundSigFigs)
	dLogalad.D[:,hubInd,p]=round.(eΔVal,sigdigits=roundSigFigs)

    dLogalad.Gu[:,hubInd,p]=round.(2*hubS.Rh[hubInd]*uVal,sigdigits=roundSigFigs)
    dLogalad.Ge[:,hubInd,p]=round.(2*hubS.Qh[hubInd]*eVal.-2*hubS.Qh[hubInd]*eMax,sigdigits=roundSigFigs)
    dLogalad.GeΔ[:,hubInd,p].=-hubS.Oh[1,hubInd]
	dLogalad.timeSolve[1,p]=max(getsolvetime(hubM),dLogalad.timeSolve[1,p])
    return nothing
end

function localXFRMALAD(p::Int,stepI::Int,σZ::Float64,σT::Float64,hubS::scenarioHubStruct,dLogalad::hubItLogPWL,
	mode::String,itLam,itVz,itVt,itρ)
    S=hubS.S
    horzLen=min(hubS.K1,hubS.K-stepI)

    tolT=1e-3
    tolZ=1e-4

    #N+1 decoupled problem aka transformer current
    coorM=Model(solver = GurobiSolver(NumericFocus=3))
    @variable(coorM,t[1:(horzLen+1)])
    @constraint(coorM,upperTCon,t.<=hubS.Tmax)
    @constraint(coorM,t.>=0)
    if mode=="NL"
        # @variable(coorM,itotal[1:(horzLen+1)])
        # @objective(coorM,Min,sum(itLam[k,1]*itotal[k] for k=1:(horzLen+1)))
        # @constraint(coorM,itotal.<=hubS.ItotalMax)
        # @constraint(coorM,itotal.>=0)
        # @NLconstraint(coorM,tempCon1,t[1,1]==hubS.τP*t0+hubS.γP*(itotal[1])^2+hubS.ρP*hubS.Tamb[stepI,1])
        # @NLconstraint(coorM,tempCon2[k=1:horzLen],t[k+1,1]==hubS.τP*t[k,1]+hubS.γP*(itotal[k+1])^2+hubS.ρP*hubS.Tamb[stepI+k,1])
    elseif mode=="relax1"
        # @variable(coorM,itotal[1:(horzLen+1)])
        # @objective(coorM,Min,sum(itLam[k,1]*itotal[k] for k=1:(horzLen+1)))
        # @constraint(coorM,itotal.<=hubS.ItotalMax)
        # @constraint(coorM,itotal.>=0)
        # @constraint(coorM,tempCon1,t[1,1]>=hubS.τP*t0+hubS.γP*(itotal[1])^2+hubS.ρP*hubS.Tamb[stepI,1])
        # @constraint(coorM,tempCon2[k=1:horzLen],t[k+1,1]>=hubS.τP*t[k,1]+hubS.γP*(itotal[k+1])^2+hubS.ρP*hubS.Tamb[stepI+k,1])
    elseif mode=="PWL"
        @variable(coorM,z[1:(horzLen+1),1:hubS.S])
        objExp=sum(-itLam[k,1]*sum(z[k,s] for s=1:S)+
                  ρALADp/2*σZ*(sum(z[k,s] for s=1:S)-sum(itVz[k,s] for s=1:S))^2+
                  ρALADp/2*σT*(t[k]-itVt[k,1])^2  for k=1:(horzLen+1))
        @objective(coorM,Min, objExp)
        @constraint(coorM,tempCon1,t[1,1]==hubS.τP*t0+hubS.γP*hubS.deltaI*sum((2*s-1)*z[1,s] for s=1:hubS.S)+hubS.ρP*hubS.Tamb[stepI,1])
        @constraint(coorM,tempCon2[k=1:horzLen],t[k+1,1]==hubS.τP*t[k,1]+hubS.γP*hubS.deltaI*sum((2*s-1)*z[k+1,s] for s=1:hubS.S)+hubS.ρP*hubS.Tamb[stepI+k,1])
        @constraint(coorM,z.>=0)
        @constraint(coorM,z.<=hubS.deltaI)
    end

	if solverSilent
        @suppress_out begin
			statusC = solve(coorM)
        end
    else
		statusC = solve(coorM)
    end

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

    dLogalad.Tpred[:,:,p]=round.(tVal,sigdigits=roundSigFigs)
    dLogalad.Z[:,:,p]=round.(zVal,sigdigits=roundSigFigs)
    dLogalad.Gz[:,:,p].=0
    dLogalad.Gt[:,:,p].=0
	dLogalad.timeSolve[1,p]=max(getsolvetime(coorM),dLogalad.timeSolve[1,p])
    return nothing
end

function coordALAD(p::Int,stepI::Int,μALADp::Float64,hubS::scenarioHubStruct,dLogalad::hubItLogPWL,
	itLam,itVu,itVe,itVd,itVz,itVt,itρ)

    H=hubS.H
    S=hubS.S
    horzLen=min(hubS.K1,hubS.K-stepI)

    Hu=2*hubS.Rh
    He=2*hubS.Qh
    # Hu=2*hubS.Ri *((1.5-2.5)*rand()+2.5)
    # Hs=2*hubS.Qsi *((1.5-2.5)*rand()+2.5)
    Hz=0
    Ht=0
    HeΔ=zeros(H)
    # Hz=1e-6
    # Ht=1e-6
    ρRate=1.1
    ρALADmax=1e6

    #coupled QP
    cM = Model(solver = GurobiSolver(NumericFocus=3))
    @variable(cM,dU[1:(horzLen+1),1:H])
    @variable(cM,dE[1:(horzLen+1),1:H])
    @variable(cM,dEΔ[1:(horzLen+1),1:H])
    @variable(cM,dZ[1:(horzLen+1),1:S])
    @variable(cM,dT[1:(horzLen+1)])
    @variable(cM,relaxS[1:(horzLen+1)])

    objExp=sum(sum(0.5*dU[k,h]^2*Hu[h]+dLogalad.Gu[k,h,p]*dU[k,h] +
                   0.5*dE[k,h]^2*He[h]+dLogalad.Ge[k,h,p]*dE[k,h] +
                   0.5*dEΔ[k,h]^2*HeΔ[h]+dLogalad.GeΔ[k,h,p]*dEΔ[k,h] for h=1:H) +
               sum(0.5*dZ[k,s]^2*Hz for s=1:S)+
               0.5*dT[k,1]^2*Ht   for k=1:(horzLen+1))
    objExp=objExp+itLam[:,1]'*relaxS+μALADp/2*sum(relaxS[k,1]^2 for k=1:horzLen+1)
    objExp=objExp+dot(dLogalad.Gz[:,:,p],dZ)+dot(dLogalad.Gt[:,1,p],dT)

    @objective(cM,Min,objExp)

    Up=round.(dLogalad.U[:,:,p],digits=6)
    Zp=round.(dLogalad.Z[:,:,p],digits=6)
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

	if solverSilent
        @suppress_out begin
			statusM = solve(cM)
        end
    else
		statusM = solve(cM)
    end

    @assert statusM==:Optimal "ALAD Central QP optimization not solved to optimality"

    #update step
    # Lam[:,p+1]=-getdual(currCon)
    α1=1
    α2=1
    α3=1
    #α3=1/ceil(p/2)

    dLogalad.Lam[:,1,p]=round.(itLam[:,1]+α3*(-getdual(currCon)-itLam[:,1]),sigdigits=roundSigFigs)
    #dLogalad.Lam[:,p]=max.(itLam[:,1]+α3*(-getdual(currCon)-itLam[:,1]),0)
    dLogalad.Vu[:,:,p]=round.(itVu[:,:,1]+α1*(dLogalad.U[:,:,p]-itVu[:,:,1])+α2*getvalue(dU),sigdigits=roundSigFigs)
    dLogalad.Vz[:,:,p]=round.(itVz[:,:,1]+α1*(dLogalad.Z[:,:,p]-itVz[:,:,1])+α2*getvalue(dZ),sigdigits=roundSigFigs)
    dLogalad.Ve[:,:,p]=round.(itVe[:,:,1]+α1*(dLogalad.E[:,:,p]-itVe[:,:,1])+α2*getvalue(dE),sigdigits=roundSigFigs)
    dLogalad.Vd[:,:,p]=round.(itVd[:,:,1]+α1*(dLogalad.Vd[:,:,p]-itVd[:,:,1])+α2*getvalue(dEΔ),sigdigits=roundSigFigs)
    dLogalad.Vt[:,:,p]=round.(itVt[:,:,1]+α1*(dLogalad.Tpred[:,:,p]-itVt[:,:,1])+α2*getvalue(dT),sigdigits=roundSigFigs)
	dLogalad.timeSolve[1,p]+=getsolvetime(cM)

    #dCMalad.lamIt[p,1]=norm(dLogalad.Lam[:,p]-itLam[:,1],2)
    #dCMalad.lam[p,1]=norm(dLogalad.Lam[:,p]-cSol.lamCoupl[stepI:(horzLen+stepI)],2)
    # if !silent
    #     @printf "lastGap    %e after %g iterations\n" dCMalad.lamIt[p,1] p
    #     @printf "convLamGap %e after %g iterations\n\n" dCMalad.lam[p,1] p
    # end

    dLogalad.itUpdate[1,1,p]=min(ρALADp*ρRate,ρALADmax) #increase ρ every iteration
    #μALADp=min(μALADp*μRate,μALADmax) #increase μ every iteration
    #ΔY[1,p]=norm(vcat(getvalue(dUn),getvalue(dZ),getvalue(dSn),getvalue(dT)),Inf)
    return nothing
end

function runEVALADIt(p,stepI,hubS,itLam,itVu,itVz,itVe,itVd,itVt,itρ,dLogalad,dSol,cSol,mode,eqForm,silent)
    H=hubS.H
    K=hubS.K
    S=hubS.S
    horzLen=min(hubS.K1,K-stepI)

    #initialize with current states
    global e0
    global t0
    # global prevLam
    # global prevVu
    # global prevVz
    # global prevVt
    # global prevVe
    # global prevVd
    # global ρALADp

    #other parameters
    epsilon = 1e-2

    #ALADIN tuning
    # if eqForm
    #     #println("Running Eq ALADIN")
    #     scalingF=1.0
    # else
    #     #println("Running ineq ALADIN")
    #     scalingF=1e-2
    # end
    # σZ=scalingF*10
    # σT=scalingF
    # σU=ones(H,1)*scalingF*10
    # σE=ones(H,1)*scalingF
    # σD=ones(H,1)*scalingF*10
	σZ=1/8
    σT=1/200
    σU=ones(H,1)/50
    σE=ones(H,1)
    σD=ones(H,1)*5
    μALADp=1e8
    # μALAD=1e8
    # μRate=1
    # μALADmax=2e9

    #solve decoupled
    if runParallel
        @sync @distributed for hubInd=1:H
			localEVALAD(hubInd,p,stepI,σU,σE,σD,hubS,dLogalad,itLam,itVu,itVe,itVd,e0,itρ,solverSilent)
		end
    else
        for hubInd=1:H
			localEVALAD(hubInd,p,stepI,σU,σE,σD,hubS,dLogalad,itLam,itVu,itVe,itVd,e0,itρ,solverSilent)
		end
    end

    localXFRMALAD(p,stepI,σZ,σT,hubS,dLogalad,mode,itLam,itVz,itVt,itρ)

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

	coordALAD(p,stepI,μALADp,hubS,dLogalad,itLam,itVu,itVe,itVd,itVz,itVt,itρ)

    #check for convergence
    constGap=norm(dLogalad.couplConst[:,1,p],1)
	itGap = norm(dLogalad.Lam[:,1,p]-itLam[:,1],2)
	#dCM.coupl1Norm[p,1]=constGap
	#dCM.lamIt2Norm[p,1]=itGap
    auxGap=norm(hcat(σU[1]*(itVu[:,:]-dLogalad.U[:,:,p]),σZ*(itVz[:,:]-dLogalad.Z[:,:,p]),
                 σT*(itVt[:,:]-dLogalad.Tpred[:,:,p]),σE[1]*(itVe[:,:]-dLogalad.E[:,:,p])),1)
    #cc=ρALAD*norm(vcat(repeat(σU,horzLen+1,1).*(Vu[:,p]-Un[:,p+1]),σZ*(Vz[:,p]-Z[:,p+1])),1)
    # objFun(sn,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*hubS.Qsi[n,1] for n=1:N) +
    #                  sum((u[(k-1)*N+n,1])^2*hubS.Ri[n,1]       for n=1:N) for k=1:horzLen+1)
    # dLogalad.objVal[1,p]=objFun(dLogalad.Sn[:,p],dLogalad.Un[:,p])
    #fGap= abs(dLogalad.objVal[1,p]-cSol.objVal[1,1])
    # snGap=norm((dLogalad.Sn[:,p]-cSol.Sn),2)
    # unGap=norm((dLogalad.Un[:,p]-cSol.Un),2)
    #dCMalad.obj[p,1]=fGap
    # dCMalad.sn[p,1]=snGap
    # dCMalad.un[p,1]=unGap
    # dCMalad.couplConst[p,1]=constGap
    #convCheck[p,1]=cc
    if  constGap<=primChk && auxGap<=auxChk
        if !silent @printf "Converged after %g iterations\n" p end
        convIt=p
        #break
        return true
    else
        if !silent
            @printf "auxGap     %e after %g iterations\n" auxGap p
            @printf "constGap   %e after %g iterations\n" constGap p
            #@printf "snGap      %e after %g iterations\n" snGap p
            #@printf("fGap       %e after %g iterations\n",fGap,p)
        end
    end

    return false
end

function runHubALADStep(stepI,maxIt,hubS,dSol,cSol,mode,eqForm,silent)
    K=hubS.K
    S=hubS.S
    horzLen=min(hubS.K1,K-stepI)
    #ogρ=ρALADp #save to reset later
    dLogalad=hubItLogPWL(horzLen=horzLen,H=hubS.H,S=hubS.S)
	timeStart=now()
	p=1
    while (p<=maxIt && round(now()-timeStart,Second)<=Dates.Second(9/10*hubS.Ts))
		# global p
		if p==1
			itLam=prevLam
			itVu=prevVu
			itVz=prevVz
			itVe=prevVe
			itVd=prevVd
			itVt=prevVt
			itρ=ρALADp
		else
			itLam=round.(dLogalad.Lam[:,1,(p-1)],digits=8)
			itVu=round.(dLogalad.Vu[:,:,(p-1)],digits=8)
			itVz=round.(dLogalad.Vz[:,:,(p-1)],digits=8)
			itVd=round.(dLogalad.Vd[:,:,(p-1)],digits=6)
			itVe=round.(dLogalad.Ve[:,:,(p-1)],digits=6)
			itVt=round.(dLogalad.Vt[:,1,(p-1)],digits=6)
			itρ=round.(dLogalad.itUpdate[1,1,(p-1)],digits=2)
		end
        cFlag=runEVALADIt(p,stepI,hubS,itLam,itVu,itVz,itVe,itVd,itVt,itρ,dLogalad,dSol,cSol,mode,eqForm,silent)
		global convIt=p
        if cFlag
            break
        end
        p+=1
    end

    # plot(dLogalad.U[:,:,convIt])
    # plot(dLogalad.E[:,:,convIt])
    # pd3alad=plot(hcat(dLogalad.Tactual[:,1,convIt],dLogalad.Tpred[:,1,convIt])*1000,label=["Actual Temp" "PWL Temp"],xlims=(0,hubS.K),xlabel="Time",ylabel="Temp (K)")
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
    # uSumPlotalad=plot(dLogalad.uSum[:,1,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Current Sum (kA)",xlims=(0,horzLen+1),legend=false)
    # plot!(uSumPlotalad,cSol.uSum,seriescolor=:black,linewidth=2,linealpha=0.8)
	#
    # zSumPlotalad=plot(dLogalad.zSum[:,1,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Z sum",xlims=(0,horzLen+1),legend=false)
    # plot!(zSumPlotalad,cSol.zSum,seriescolor=:black,linewidth=2,linealpha=0.8)
	#
    # constPlotalad2=plot(dLogalad.couplConst[:,1,1:convIt],seriescolor=CList,xlabel="Time",ylabel="curr constraint diff",xlims=(0,horzLen+1),legend=false)
	#
    # lamPlotalad=plot(dLogalad.Lam[:,1,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Lambda",xlims=(0,horzLen+1),legend=false)
    # plot!(lamPlotalad,cSol.Lam,seriescolor=:black,linewidth=2,linealpha=0.8)
	# plot!(lamPlotalad,lambdaCurr,seriescolor=:black,linewidth=2,linealpha=0.8)
	#
    # activeSet=zeros(convIt,1)
    # setChanges=zeros(convIt,1)
    # for ii=2:convIt
    #     activeSet[ii,1]=sum(abs.(dLogalad.Ceu[:,:,ii]))+sum(abs.(dLogalad.Ctu[:,:,ii]))+sum(abs.(dLogalad.Cdu[:,:,ii]))+
    #               		sum(abs.(dLogalad.Cuu[:,:,ii]))+sum(abs.(dLogalad.Czu[:,:,ii]))+
    # 					sum(abs.(dLogalad.Cel[:,:,ii]))+sum(abs.(dLogalad.Ctl[:,:,ii]))+sum(abs.(dLogalad.Cdl[:,:,ii]))+
    # 				    sum(abs.(dLogalad.Cul[:,:,ii]))+sum(abs.(dLogalad.Czl[:,:,ii]))
    #     setChanges[ii,1]=sum(abs.(dLogalad.Ceu[:,:,ii]-dLogalad.Ceu[:,:,ii-1]))+sum(abs.(dLogalad.Ctu[:,:,ii]-dLogalad.Ctu[:,:,ii-1]))+
    #                      sum(abs.(dLogalad.Cuu[:,:,ii]-dLogalad.Cuu[:,:,ii-1]))+sum(abs.(dLogalad.Czu[:,:,ii]-dLogalad.Czu[:,:,ii-1]))+
    # 					 sum(abs.(dLogalad.Cel[:,:,ii]-dLogalad.Cel[:,:,ii-1]))+sum(abs.(dLogalad.Ctl[:,:,ii]-dLogalad.Ctl[:,:,ii-1]))+
    #                      sum(abs.(dLogalad.Cdu[:,:,ii]-dLogalad.Cdu[:,:,ii-1]))+sum(abs.(dLogalad.Cdl[:,:,ii]-dLogalad.Cdl[:,:,ii-1]))+
    # 				     sum(abs.(dLogalad.Cul[:,:,ii]-dLogalad.Cul[:,:,ii-1]))+sum(abs.(dLogalad.Czl[:,:,ii]-dLogalad.Czl[:,:,ii-1]))
    # end
	#
    # activeSetPlot=plot(2:convIt,activeSet[2:convIt],xlabel="Iteration",ylabel="Total Active inequality constraints",
    #                    legend=false,xlims=(2,convIt))
    # setChangesPlot=plot(10:convIt,setChanges[10:convIt],xlabel="Iteration",ylabel="Changes in Active inequality constraints",
    #                   legend=false,xlims=(2,convIt))
    # #solChangesplot=plot(2:convIt,hcat(ΔY[2:convIt],convCheck[2:convIt]),xlabel="Iteration",labels=["ΔY" "y-x"],xlims=(1,convIt))
	#
    # #fPlotalad=plot(dCMalad.obj[1:convIt,1],xlabel="Iteration",ylabel="obj function gap",xlims=(1,convIt),legend=false,yscale=:log10)
    # convItPlotalad=plot(dCM.lamIt2Norm[1:convIt,1],xlabel="Iteration",ylabel="2-Norm Lambda Gap",xlims=(1,convIt),legend=false,yscale=:log10) #,yscale=:log10
    # #convPlotalad=plot(dCMalad.lam[1:convIt,1],xlabel="Iteration",ylabel="central lambda gap",xlims=(1,convIt),legend=false,yscale=:log10)
    # constPlotalad=plot(dCM.coupl1Norm[1:convIt,1],xlabel="Iteration",ylabel="curr constraint Gap",xlims=(1,convIt),legend=false,yscale=:log10)

    #save current state and update for next timeSteps
    dSol.Tpred[stepI,1]=dLogalad.Tpred[1,1,convIt]
    dSol.U[stepI,:]=dLogalad.U[1,1:H,convIt]
    dSol.E[stepI,:]=dLogalad.E[1,1:H,convIt]
	dSol.D[stepI,:]=dLogalad.D[1,1:H,convIt]
    #dSol.Z[stepI,:]=dLogalad.Z[1:S,convIt]
    dSol.uSum[stepI,1]=dLogalad.uSum[1,1,convIt]
    dSol.zSum[stepI,1]=dLogalad.zSum[1,1,convIt]
    dSol.Itotal[stepI,1]=dLogalad.Itotal[1,1,convIt]
    dSol.Tactual[stepI,1]=dLogalad.Tactual[1,1,convIt]
    dSol.convIt[stepI,1]=convIt
	dSol.timeSolve[stepI,1]=mean(dLogalad.timeSolve[1,1:convIt])

    # new states
    global t0=dSol.Tactual[stepI,1]
    global e0=dSol.E[stepI,:]

    #function getAttr()
    #clean this up
    if convIt==1
        dSol.Lam[stepI,1]=prevLam[1,1]
        if stepI+horzLen==hubS.K
            newLam=prevLam[2:horzLen+1,1]
            newVu=prevVu[2:horzLen+1,:]
            newVz=prevVz[2:horzLen+1,:]
            newVt=prevVt[2:horzLen+1,1]
            newVe=prevVe[2:horzLen+1,:]
            newVd=prevVd[2:horzLen+1,:]
        else
            newLam=vcat(prevLam[2:horzLen+1,1],prevLam[horzLen+1,1])
            newVu=vcat(prevVu[2:horzLen+1,:],prevVu[horzLen+1,:]')
            newVz=vcat(prevVz[2:horzLen+1,:],prevVz[horzLen+1,:]')
            newVt=vcat(prevVt[2:horzLen+1,:],prevVt[horzLen+1,:])
            newVe=vcat(prevVe[2:horzLen+1,:],prevVe[horzLen+1,:]')
            newVd=vcat(prevVd[2:horzLen+1,:],prevVd[horzLen+1,:]')
        end
    else
        dSol.Lam[stepI,1]=dLogalad.Lam[1,1,convIt]
        if stepI+horzLen==hubS.K
            newLam=dLogalad.Lam[2:horzLen+1,1,convIt]
            newVu=dLogalad.Vu[2:horzLen+1,:,convIt]
            newVz=dLogalad.Vz[2:horzLen+1,:,convIt]
            newVt=dLogalad.Vt[2:horzLen+1,1,convIt]
            newVe=dLogalad.Ve[2:horzLen+1,:,convIt]
            newVd=dLogalad.Vd[2:horzLen+1,:,convIt]
        else
            newLam=vcat(dLogalad.Lam[2:horzLen+1,:,convIt],dLogalad.Lam[horzLen+1,:,convIt])
            newVu=vcat(dLogalad.Vu[2:horzLen+1,:,convIt],dLogalad.Vu[horzLen+1,:,convIt]')
            newVz=vcat(dLogalad.Vz[2:horzLen+1,:,convIt],dLogalad.Vz[horzLen+1,:,convIt]')
            newVt=vcat(dLogalad.Vt[2:horzLen+1,:,convIt],dLogalad.Vt[horzLen+1,:,convIt])
            newVe=vcat(dLogalad.Ve[2:horzLen+1,:,convIt],dLogalad.Ve[horzLen+1,:,convIt]')
            newVd=vcat(dLogalad.Vd[2:horzLen+1,:,convIt],dLogalad.Vd[horzLen+1,:,convIt]')
        end
    end

    global prevLam=newLam
    global prevVu=newVu
    global prevVz=newVz
    global prevVt=newVt
    global prevVe=newVe
    global prevVd=newVd
    #global ρALADp=ogρ

    return nothing
end

function hubALAD(maxIt::Int,hubS::scenarioHubStruct,cSol::hubSolutionStruct,mode::String,eqForm::Bool,silent::Bool)
    H=hubS.H
    K=hubS.K
    dSol=hubSolutionStruct(K=K,H=H)
	dCM=convMetricsStruct(maxIt=maxIt,logLength=1)

	p = plot(2,label=["Central" "ALAD"])
    Juno.progress() do id
        for stepI=1:K
            @info "$(Dates.format(Dates.now(),"HH:MM:SS")): $(stepI) of $(K)....\n" progress=stepI/K _id=id
            @printf "%s: time step %g of %g...." Dates.format(Dates.now(),"HH:MM:SS") stepI K
            try
                dSol.timeT[stepI]=@elapsed runHubALADStep(stepI,maxIt,hubS,dSol,cSol,mode,eqForm,silent)
                @printf "convIt: %g\n" dSol.convIt[stepI,1]
            catch e
                @printf "error: %s" e
                break
            end
			push!(p, stepI, [cSol.Lam[stepI], dSol.Lam[stepI]])
			display(p)
        end
    end

	objFun(e,u,d)=sum(sum((e[k,h]-hubS.eMax[k,h])^2*hubS.Qh[1,h]+(u[k,h])^2*hubS.Rh[1,h]-hubS.Oh[h]*d[k,h] for h=1:H) for k=1:K)
    dSol.objVal[1,1]=objFun(dSol.E,dSol.U,dSol.D)

    return dSol
end
