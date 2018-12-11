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


#ALADIN
