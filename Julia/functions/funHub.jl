# hub functions

function runHubCentralStep(stepI,hubS,cSol,mode,silent)

    global t0
    global e0
    H=hubS.H
    K=hubS.K
    Nh=hubS.Nh
    @printf "time step %g of %g....\n" stepI K

    horzLen=min(hubS.K1,K-stepI)

    #prepare predicted values for optimization

    #ηH=zeros(horzLen+1,1)
    eMax=zeros(horzLen+1,H)
    uMax=zeros(horzLen+1,H)
    eDepart=zeros(horzLen+1,H)
    eArrive_pred=zeros(horzLen+1,H)
    eArrive_actual=zeros(horzLen+1,H)
    slackMax=zeros(horzLen+1,H)
    for (i,k) in enumerate(stepI:(stepI+horzLen))
        for h=1:H # find a way to do this without another loop
            depart=[n for n=1:Nh if k==hubS.K_depart_pred[n,h]]
            if length(depart)!=0
                eDepart[i,h] = sum(hubS.Sn_depart_min[n,h]*hubS.EVcap[n,h] for n in depart)
                slackMax[i,h] = sum(hubS.EVcap[n,h]-hubS.Sn_depart_min[n,h]*hubS.EVcap[n,h] for n in depart)
            end

            # arrive=[n for n=1:Nh,h=1:H if k==hubS.K_arrive_pred[n,h]]
            arrive=[n for n=1:Nh if k==hubS.K_arrive_pred[n,h]]
            if length(arrive)!=0
                eArrive_pred[i,h] = sum(hubS.Sn_arrive_pred[n,h]*hubS.EVcap[n,h] for n in arrive)
                eArrive_actual[i,h] = sum(hubS.Sn_arrive_actual[n,h]*hubS.EVcap[n,h] for n in arrive)
            end

            parked=[n for n=1:Nh  if hubS.K_arrive_pred[n,h]<=k<=hubS.K_depart_pred[n,h]]
            if length(parked)!=0
                #ηH[k]=mean(eta[n] for n in parked)
                uMax[i,h]=sum(hubS.imax[n,h] for n in parked)
                eMax[i,h]=sum(hubS.EVcap[n,h] for n in parked)
            end
        end
    end

    if silent
        #cModel = Model(solver = IpoptSolver(print_level=0))
        cModel = Model(solver = GurobiSolver(OutputFlag=0,QCPDual=1))
        cModel = Model(solver = GurobiSolver(OutputFlag=0))
    else
        #cModel = Model(solver = IpoptSolver())
        cModel = Model(solver = GurobiSolver(QCPDual=1))
        cModel = Model(solver = GurobiSolver())

    end

    @variable(cModel,u[1:(horzLen+1),1:H])
    @variable(cModel,t[1:(horzLen+1)])
    @variable(cModel,e[1:(horzLen+1),1:H])
    @variable(cModel,slackE[1:(horzLen+1),1:H])

    #objective
    @objective(cModel,Min,sum(sum(hubS.Q[h]*(e[k,h]-eMax[k,h])^2+hubS.R[h]*u[k,h]^2 for k=1:horzLen+1) for h=1:H))

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
        @variable(cModel,z[1:(horzLen+1)*hubS.S])
        @constraint(cModel,tempCon1,t[1,1]==hubS.τP*t0+hubS.γP*hubS.deltaI*sum((2*s-1)*z[s,1] for s=1:hubS.S)+hubS.ρP*hubS.Tamb[stepI,1])
        @constraint(cModel,tempCon2[k=1:horzLen],t[k+1,1]==hubS.τP*t[k,1]+hubS.γP*hubS.deltaI*sum((2*s-1)*z[k*hubS.S+s,1] for s=1:hubS.S)+hubS.ρP*hubS.Tamb[stepI+k,1])
        @constraint(cModel,z.>=0)
        @constraint(cModel,z.<=hubS.deltaI)
        #coupling constraint
        @constraint(cModel,currCon[k=1:horzLen+1],0==-sum(u[k,h] for h=1:H)-hubS.iD_pred[stepI+(k-1)]+sum(z[(k-1)*(hubS.S)+s] for s=1:hubS.S))
    end

    #hub constraints
    @constraint(cModel,stateCon1[h=1:H],e[1,h]==e0[h]+hubS.ηP[h]*u[1,h]-(eDepart[1,h]+slackE[1,h])+eArrive_pred[1,h])
    @constraint(cModel,stateCon[k=1:horzLen,h=1:H],e[k+1,h]==e[k,h]+hubS.ηP[h]*u[k+1,h]-(eDepart[k+1,h]+slackE[k+1,h])+eArrive_pred[k+1,h])
    @constraint(cModel,e.>=0)
    @constraint(cModel,eMaxCon[k=1:horzLen+1,h=1:H],e[k,h]<=eMax[k,h])
    @constraint(cModel,uMaxCon,u.<=uMax)
    @constraint(cModel,u.>=0)
    @constraint(cModel,slackE.<=slackMax)
    @constraint(cModel,slackE.>=0)

    status=solve(cModel)
    @assert status==:Optimal "Central Hub optimization not solved to optimality"

    eRaw=getvalue(e)
    uRaw=getvalue(u)
    tRaw=getvalue(t)
    lambdaCurr=-getdual(currCon)
    extraE=getvalue(slackE)
    #
    # p1nl=plot(eRaw,xlabel="Time",ylabel="Hub Energy",legend=false,xlims=(1,horzLen+1))
    # p2nl=plot(uRaw,xlabel="Time",ylabel="Hub Current (kA)",legend=false,xlims=(1,horzLen+1))
    # p3nl=plot(1:horzLen+1,tRaw,label="XFRM Temp",xlims=(1,horzLen+1),xlabel="Time",ylabel="Temp (K)")
    # plot!(p3nl,1:horzLen+1,Tmax*ones(horzLen+1),label="XFRM Limit",line=(:dash,:red))
    # p4nl=plot(1:horzLen+1,lambdaCurr,label="Time",ylabel=raw"Lambda ($/kA)",xlims=(1,horzLen+1),legend=false)
    #
    # plot(sum(eMax,dims=2),label="max")
    # plot!(sum(eRaw,dims=2),label="e")

    #apply current and actual departures and arrivals
    nextU=uRaw[1,:]
    cSol.U[stepI,:]=nextU
    cSol.Lam[stepI,1]=lambdaCurr[1,1]
    cSol.E_depart[stepI,:]=eDepart[1,:]+extraE[1,:]
    cSol.E_arrive[stepI,:]=eArrive_actual[1,:]
    cSol.E[stepI,:]=e0[:]+hubS.ηP[:].*nextU-(cSol.E_depart[stepI,:])+cSol.E_arrive[stepI,:]
    cSol.T[stepI,1]=hubS.τP*t0+hubS.γP*(sum(nextU)+hubS.iD_actual[stepI,1])^2+hubS.ρP*hubS.Tamb[stepI,1]

    # new states
    t0=cSol.T[stepI,1]
    e0=cSol.E[stepI,:]

    return nothing
end
