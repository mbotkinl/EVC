#functions to run PWL EVC problems


#central
function pwlEVcentral(N::Int,S::Int,horzLen::Int,evS::scenarioStruct,slack::Bool)
    #initialize with current states
    sn0=evS.s0
    xt0=evS.t0
    stepI=1
    #desired SOC
    target=zeros(N*(horzLen+1),1);
    for ii=1:N
       cur=evS.Kn[ii]-(stepI-1)
       ind=(max(0,(cur-1)*N)+ii):N:length(target)
       target[ind].=evS.Snmin[ii,1]
    end

    println("setting up model")
    #centralModel = Model(solver = GurobiSolver(Presolve=0,BarHomogeneous=1,NumericFocus=3))
    centralModel = Model(solver = GurobiSolver())

    #u iD and z are one index ahead of sn and T. i.e the x[k+1]=x[k]+eta*u[k+1]
    @variable(centralModel,u[1:N*(horzLen+1)])
    @variable(centralModel,sn[1:(N)*(horzLen+1)])
    @variable(centralModel,xt[1:(horzLen+1)])
    @variable(centralModel,z[1:evS.S*(horzLen+1)])
    if slack
        @variable(centralModel,slackSn[1:N])
    end
    println("obj")

    objExp =sum((sn[n,1]-1)^2*evS.Qsi[n,1]+(u[n,1])^2*evS.Ri[n,1] for n=1:N)
    for k=2:horzLen+1
        append!(objExp,sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,1]+(u[(k-1)*N+n,1])^2*evS.Ri[n,1]  for n=1:N))
    end
    if slack append!(objExp,sum(evS.β[n]*slackSn[n]^2 for n=1:N)) end

    #objExp=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,1]+(u[(k-1)*N+n,1])^2*evS.Ri[n,1]  for n=1:N) for k=1:horzLen+1)

    @objective(centralModel,Min,objExp)
    println("constraints")
    @constraint(centralModel,stateCon1,sn[1:N,1].==sn0[1:N,1]+evS.ηP[:,1].*u[1:N,1])
    @constraint(centralModel,stateCon2[k=1:horzLen,n=1:N],sn[n+(k)*(N),1]==sn[n+(k-1)*(N),1]+evS.ηP[n,1]*u[n+(k)*(N),1])
    @constraint(centralModel,tempCon1,xt[1,1]==evS.τP*xt0+evS.γP*evS.deltaI*sum((2*s-1)*z[s,1] for s=1:S)+evS.ρP*evS.Tamb[stepI,1])
    @constraint(centralModel,tempCon2[k=1:horzLen],xt[k+1,1]==evS.τP*xt[k,1]+evS.γP*evS.deltaI*sum((2*s-1)*z[k*S+s,1] for s=1:S)+evS.ρP*evS.Tamb[stepI+k,1])
    @constraint(centralModel,currCon[k=1:horzLen+1],0==-sum(u[(k-1)*(N)+n] for n=1:N)-evS.iD[stepI+(k-1)]+sum(z[(k-1)*(S)+s] for s=1:S))
    @constraint(centralModel,sn.<=1)
    if slack
        @constraint(centralModel,sn.>=target.*(1-repeat(slackSn,horzLen+1,1)))
        @constraint(centralModel,slackSn.>=0)
    else
        @constraint(centralModel,sn.>=target)
    end
    if noTlimit==false
    	@constraint(centralModel,upperTCon,xt.<=evS.Tmax)
    end
    @constraint(centralModel,xt.>=0)
    @constraint(centralModel,upperCCon,u.<=repeat(evS.imax,horzLen+1,1))
    @constraint(centralModel,u.>=repeat(evS.imin,horzLen+1,1))
    @constraint(centralModel,z.>=0)
    @constraint(centralModel,z.<=evS.deltaI)

    println("solving....")
    statusM = solve(centralModel)
    @assert statusM==:Optimal "Central optimization not solved to optimality"

    uRaw=getvalue(u)
    snRaw=getvalue(sn)
    xtRaw=getvalue(xt)
    zRaw=getvalue(z)
    if slack
        slackSnRaw=getvalue(slackSn)
    else
        slackSnRaw=zeros(N)
    end

    #calculate actual temp
    Tactual=zeros(horzLen+1,1)
    itotal=zeros(horzLen+1,1)
    for k=1:horzLen+1
        itotal[k,1]=sum((uRaw[(k-1)*N+n,1]) for n=1:N) + evS.iD[stepI+(k-1),1]
    end
    Tactual[1,1]=evS.τP*xt0+evS.γP*itotal[1,1]^2+evS.ρP*evS.Tamb[stepI,1]
    for k=1:horzLen
        Tactual[k+1,1]=evS.τP*Tactual[k,1]+evS.γP*itotal[k+1,1]^2+evS.ρP*evS.Tamb[stepI+(k-1),1]
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

    cSol=centralSolutionStruct(Xt=xtRaw,Un=uRaw,Sn=snRaw,z=zRaw,
                        Itotal=itotal,uSum=uSum,zSum=zSum,
                        objVal=getobjectivevalue(centralModel),
                        lamTemp=lambdaUpperT,lamCoupl=lambdaCurr,
                        Tactual=Tactual,slackSn=slackSnRaw)
    return cSol
end

#dual
function pwlEVdual(N::Int,S::Int,horzLen::Int,maxIt::Int,updateMethod::String,evS::scenarioStruct,
    cSol::centralSolutionStruct,slack::Bool)

    #initialize
    sn0=evS.s0
    xt0=evS.t0

    dCM=convMetricsStruct()
    dLog=itLogPWL()

    if updateMethod=="fastAscent"
    	#alpha = 0.1  #for A
    	alpha = 5e4 #for kA
    	alphaDivRate=2
    	minAlpha=1e-6
    else
    	#alpha = 3e-3 #for A
    	alpha = 2e3 #for kA
    	alphaDivRate=8
    	minAlpha=1e-6
    end

    stepI = 1
    convChk = 1e-16
    convIt=maxIt

    #u iD and z are one index ahead of sn and T. i.e the x[k+1]=x[k]+eta*u[k+1]
    alphaP=alpha*ones(maxIt+1,1)

    #initialize with guess
    lambda0=1000*ones(horzLen+1,1)
    #lambda0=lamCurrStar

    # dLog.Lam[:,1]=lambda0
    prevLam=lambda0

    #iterate at each time step until convergence
    for p=1:maxIt
        #solve subproblem for each EV
    	@sync @distributed for evInd=1:N
            target=zeros((horzLen+1),1)
    		target[(evS.Kn[evInd,1]-(stepI-1)):1:length(target),1].=evS.Snmin[evInd,1]
            evM=Model(solver = GurobiSolver(NumericFocus=1))
            @variable(evM,un[1:horzLen+1])
            @variable(evM,sn[1:horzLen+1])
            if slack @variable(evM,slackSn) end
            objExp=sum((sn[k,1]-1)^2*evS.Qsi[evInd,1]+(un[k,1])^2*evS.Ri[evInd,1]+prevLam[k,1]*un[k,1] for k=1:horzLen+1)
            if slack
                append!(objExp,sum(evS.β[n]*slackSn^2 for n=1:N))
            end
            @objective(evM,Min,objExp)
    		@constraint(evM,sn[1,1]==sn0[evInd,1]+evS.ηP[evInd,1]*un[1,1]) #fix for MPC loop
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

    		TT = stdout # save original stdout stream
    		redirect_stdout()
            statusEVM = solve(evM)
    		redirect_stdout(TT)

    		@assert statusEVM==:Optimal "EV NLP optimization not solved to optimality"

            dLog.Sn[collect(evInd:N:N*(horzLen+1)),p]=getvalue(sn) #solved state goes in next time slot
            dLog.Un[collect(evInd:N:N*(horzLen+1)),p]=getvalue(un) #current go
            dLog.slackSn[evInd]= if slack getvalue(slackSn) else 0 end
        end

    	if updateMethod=="dualAscent"
    	    #solve coordinator problem
    	    #coorM=Model(solver = GurobiSolver(Presolve=0,NumericFocus=1))
    		coorM=Model(solver = GurobiSolver())
    		#coorM=Model(solver = IpoptSolver())
    	    @variable(coorM,z[1:S*(horzLen+1)])
    	    @variable(coorM,xt[1:horzLen+1])
    	    @objective(coorM,Min,sum(prevLam[k,1]*sum(-z[(k-1)*S+s,1] for s=1:S) for k=1:(horzLen+1)))
    		@constraint(coorM,xt[1,1]==evS.τP*xt0+evS.γP*evS.deltaI*sum((2*s-1)*z[s,1] for s=1:S)+evS.ρP*evS.Tamb[stepI,1]) #fix for MPC loop
    		@constraint(coorM,[k=1:horzLen],xt[k+1,1]==evS.τP*xt[k,1]+evS.γP*evS.deltaI*sum((2*s-1)*z[k*S+s,1] for s=1:S)+evS.ρP*evS.Tamb[stepI+k,1])
    		if noTlimit==false
    			@constraint(coorM,upperTCon,xt.<=evS.Tmax)
    		end
    	    @constraint(coorM,xt.>=0)
    	    @constraint(coorM,z.<=evS.deltaI)
    	    @constraint(coorM,z.>=0)
    		TT = stdout # save original stdout stream
    		redirect_stdout()
    	    statusC = solve(coorM)
    		redirect_stdout(TT)

    		@assert statusC==:Optimal "Dual Ascent central optimization not solved to optimality"

    		 dLog.Xt[:,p]=getvalue(xt)
    		 dLog.Z[:,p]=getvalue(z)

    	    #grad of lagragian
    		for k=1:horzLen+1
    			dLog.zSum[k,p]=sum(dLog.Z[(k-1)*(S)+s,p] for s=1:S)
    			dLog.uSum[k,p]=sum(dLog.Un[(k-1)*N+n,p] for n=1:N)
    			dLog.couplConst[k,p]=dLog.uSum[k,p] + evS.iD[stepI+(k-1),1] - dLog.zSum[k,p]
    		end
    		dCM.couplConst[p,1]=norm(dLog.couplConst[:,p],2)
    	end

    	#calculate actual temperature from nonlinear model of XFRM
    	itotal=zeros(horzLen+1,1)
    	for k=1:horzLen+1
    		itotal[k,1]=dLog.uSum[k,p] + evS.iD[stepI+(k-1),1]
    	end
    	dLog.Tactual[1,p]=evS.τP*xt0+evS.γP*itotal[1,1]^2+evS.ρP*evS.Tamb[stepI,1]
    	for k=1:horzLen
    		dLog.Tactual[k+1,p]=evS.τP*dLog.Tactual[k,p]+evS.γP*itotal[k+1,1]^2+evS.ρP*evS.Tamb[stepI+k,1]  #fix for mpc
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
    	end

        #update lambda
    	alphaP[p,1] = max(alpha/ceil(p/alphaDivRate),minAlpha)
    	#alphaP[p+1,1] = alphaP[p,1]*alphaRate

        dLog.Lam[:,p]=max.(prevLam[:,1]+alphaP[p,1]*dLog.couplConst[:,p],0)
        #dLog.Lam[:,p]=prevLam[:,1]+alphaP[p,1]*dLog.couplConst[:,p]

    	#check convergence
    	objFun(sn,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,1]     for n=1:N) for k=1:horzLen+1) +
    					sum(sum((u[(k-1)*N+n,1])^2*evS.Ri[n,1]           for n=1:N) for k=1:horzLen+1)
        dLog.objVal[1,p]=objFun(dLog.Sn[:,p],dLog.Un[:,p])
    	fGap=abs(dLog.objVal[1,p] -cSol.objVal)
    	snGap=norm((dLog.Sn[:,p]-cSol.Sn),2)
    	unGap=norm((dLog.Un[:,p]-cSol.Un),2)
    	itGap = norm(dLog.Lam[:,p]-prevLam[:,1],2)
    	if updateMethod=="fastAscent"
    		convGap = norm(dLog.Lam[:,p]-cSol.lamTemp,2)
    	else
    		convGap = norm(dLog.Lam[:,p]-cSol.lamCoupl,2)
    	end
    	dCM.obj[p,1]=fGap
    	dCM.sn[p,1]=snGap
    	dCM.un[p,1]=unGap
    	dCM.lamIt[p,1]=itGap
    	dCM.lam[p,1]=convGap
    	if(itGap <= convChk )
    		@printf "Converged after %g iterations\n" p
    		convIt=p
    		break
    	else
    		@printf "lastGap %e after %g iterations\n" itGap p
    		@printf "convGap %e after %g iterations\n" convGap p
            @printf "snGap   %e after %g iterations\n" snGap p
    		@printf "unGap   %e after %g iterations\n" unGap p
    		@printf("fGap    %e after %g iterations\n\n",fGap,p)
            prevLam=dLog.Lam[:,p]
    	end
    end
    #
    # dLog=itLogPWL(Lam=Lam,Un=Un,Sn=Sn,Z=Z,Xt=Xt,Tactual=Tactual,
    #               uSum=uSum,zSum=zSum,couplConst=couplConst)
    # dCM=covnMetrics(objVal=objVal,couplConst=couplConst,lam=lam,
    #                 sn=sn,un=un,lamIt=lamIt,snIt=snIt,unIt=unIt)

    return (dLog,dCM,convIt)
end

#ADMM
function pwlEVadmm(N::Int,S::Int,horzLen::Int,maxIt::Int,evS::scenarioStruct,cSol::centralSolutionStruct,slack::Bool)
    #initialize
    sn0=evS.s0
    xt0=evS.t0

    stepI = 1
    convChk = 1e-16
    convIt=maxIt

    #admm  initial parameters and guesses
    #ρADMM=1    #for A
    ρADMM=1e6 #for kA
    ρDivRate=10


    #u iD and z are one index ahead of sn and T. i.e the x[k+1]=x[k]+η*u[k+1]
    dCMadmm=convMetricsStruct()
    dLogadmm=itLogPWL()

    #initialize with guess
    lambda0=1000*ones(horzLen+1,1)
    vz0=-ones(S*(horzLen+1),1)
    vu0=.01*ones(N*(horzLen+1),1)
    #vz0=-zStar
    #vu0=uStar
    #lambda0=lamCurrStar

    # dLogadmm.Lam[:,1]=lambda0
    # dLogadmm.Vz[:,1]=vz0
    # dLogadmm.Vu[:,1]=vu0
    prevLam=lambda0
    prevVz=vz0
    prevVu=vu0

    for p in 1:maxIt
    	ρI = ρADMM/ceil(p/ρDivRate)
    	#ρI = ρADMM

        #x minimization eq 7.66 in Bertsekas
        @sync @distributed for evInd=1:N
            evV=prevVu[collect(evInd:N:length(prevVu)),1]
            target=zeros((horzLen+1),1)
    		target[(evS.Kn[evInd,1]-(stepI-1)):1:length(target),1].=evS.Snmin[evInd,1]
        	evM = Model(solver = GurobiSolver())
        	@variable(evM,sn[1:(horzLen+1)])
        	@variable(evM,u[1:(horzLen+1)])
            if slack @variable(evM,slackSn) end
            objExp= sum((sn[k,1]-1)^2*evS.Qsi[evInd,1]+(u[k,1])^2*evS.Ri[evInd,1]+
    								prevLam[k,1]*(u[k,1]-evV[k,1])+
    								ρI/2*(u[k,1]-evV[k,1])^2 for k=1:horzLen+1)
            if slack
                append!(objExp,sum(evS.β[n]*slackSn^2 for n=1:N))
            end
    		@objective(evM,Min,objExp)
            @constraint(evM,sn[1,1]==sn0[evInd,1]+evS.ηP[evInd,1]*u[1,1])
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
        	TT = stdout # save original stdout stream
        	redirect_stdout()
        	statusEVM = solve(evM)
        	redirect_stdout(TT)
    		@assert statusEVM==:Optimal "ADMM EV NLP optimization not solved to optimality"

    		dLogadmm.Sn[collect(evInd:N:length(dLogadmm.Sn[:,p])),p]=getvalue(sn)
    		dLogadmm.Un[collect(evInd:N:length(dLogadmm.Un[:,p])),p]=getvalue(u)
            dLogadmm.slackSn[evInd]= if slack getvalue(slackSn) else 0 end
        end

        #N+1 decoupled problem aka transformer current
        tM = Model(solver = GurobiSolver())
        @variable(tM,z[1:(S)*(horzLen+1)])
        @variable(tM,xt[1:(horzLen+1)])
        # constFun1(u,v)=sum(Lam[k,p]*sum(u[(k-1)*(S)+s,1]-v[(k-1)*(S)+s,1] for s=1:S)  for k=1:(horzLen+1))
        # constFun2(u,v)=ρI/2*sum(sum((u[(k-1)*(S)+s,1]-v[(k-1)*(S)+s,1])*(u[(k-1)*(S)+s,1]-v[(k-1)*(S)+s,1]) for s=1:S)  for k=1:(horzLen+1))
        # @objective(tM,Min, constFun1(-z,Vz[:,p])+constFun2(-z,Vz[:,p]))
    	@objective(tM,Min,sum(prevLam[k,1]*(sum(-z[(k-1)*(S)+s,1] for s=1:S)-sum(prevVz[(k-1)*(S)+s,1] for s=1:S)) +
    						ρI/2*(sum(-z[(k-1)*(S)+s,1] for s=1:S)-sum(prevVz[(k-1)*(S)+s,1] for s=1:S))^2  for k=1:(horzLen+1)))
        @constraint(tM,tempCon1,xt[1,1]==evS.τP*xt0+evS.γP*evS.deltaI*sum((2*m+1)*z[m+1,1] for m=0:S-1)+evS.ρP*evS.Tamb[stepI,1])
        @constraint(tM,tempCon2[k=1:horzLen],xt[k+1,1]==evS.τP*xt[k,1]+evS.γP*evS.deltaI*sum((2*m+1)*z[k*S+(m+1),1] for m=0:S-1)+evS.ρP*evS.Tamb[stepI+k,1])
        if noTlimit==false
        	@constraint(tM,upperTCon,xt.<=evS.Tmax)
        end
        @constraint(tM,xt.>=0)
        @constraint(tM,z.>=0)
        @constraint(tM,z.<=evS.deltaI)
        #@constraint(tM,zC[k=1:horzLen+1],zSum[k,1]==sum(z[(k-1)*(S)+s] for s=1:S))

        TT = stdout # save original stdout stream
        redirect_stdout()
        statusC = solve(tM)
        redirect_stdout(TT)
    	@assert statusC==:Optimal "ADMM XFRM NLP optimization not solved to optimality"

        dLogadmm.Xt[:,p]=getvalue(xt)
        dLogadmm.Z[:,p]=getvalue(z)

        #lambda update eq 7.68
    	for k=1:horzLen+1
    		dLogadmm.uSum[k,p]=sum(dLogadmm.Un[(k-1)*N+n,p] for n=1:N)
    		dLogadmm.zSum[k,p]=sum(dLogadmm.Z[(k-1)*(S)+s,p] for s=1:S)
    		dLogadmm.couplConst[k,p]= dLogadmm.uSum[k,p] + evS.iD[stepI+(k-1),1] - dLogadmm.zSum[k,p]
            dLogadmm.Lam[k,p]=max.(prevLam[k,1]+ρI/(S*(N))*(dLogadmm.couplConst[k,p]),0)
    		#dLogadmm.Lam[k,p]=dLogadmm.Lam[k,p]+ρI/(S*(N))*(dLogadmm.couplConst[k,p])
    	end

    	#calculate actual temperature from nonlinear model of XFRM
        itotal=zeros(horzLen+1,1)
        for k=1:horzLen+1
            itotal[k,1]=dLogadmm.uSum[k,p] + evS.iD[stepI+(k-1),1]
        end
    	dLogadmm.Tactual[1,p]=evS.τP*xt0+evS.γP*itotal[1,1]^2+evS.ρP*evS.Tamb[stepI,1] #fix for mpc
    	for k=1:horzLen
    		dLogadmm.Tactual[k+1,p]=evS.τP*dLogadmm.Tactual[k,p]+evS.γP*itotal[k+1,1]^2+evS.ρP*evS.Tamb[stepI+k,1]  #fix for mpc
    	end

        #v upate eq 7.67
        for k=1:horzLen+1
            dLogadmm.Vu[(k-1)*N.+collect(1:N),p]=dLogadmm.Un[(k-1)*N.+collect(1:N),p].+(prevLam[k,1]-dLogadmm.Lam[k,p])/ρI
            dLogadmm.Vz[(k-1)*(S).+collect(1:S),p]=-dLogadmm.Z[(k-1)*(S).+collect(1:S),p].+(prevLam[k,1]-dLogadmm.Lam[k,p])/ρI
            #
            # dLogadmm.Vu[(k-1)*N+collect(1:N),p]=min.(max.(dLogadmm.Un[(k-1)*N+collect(1:N),p]+(prevLam[k,1]-dLogadmm.Lam[k,p])/ρI,evS.imin),evS.imax)
            # dLogadmm.Vz[(k-1)*(S)+collect(1:S),p]=max.(min.(-dLogadmm.Z[(k-1)*(S)+collect(1:S),p]+(prevLam[k,1]-dLogadmm.Lam[k,p])/ρI,0),-evS.deltaI)
        end

        #check convergence
    	objFun(sn,xt,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,1]     for n=1:N) for k=1:horzLen+1) +
    					sum((xt[k,1]-1)^2*evS.Qsi[N+1,1]                 for k=1:horzLen+1) +
    					sum(sum((u[(k-1)*N+n,1])^2*evS.Ri[n,1]           for n=1:N) for k=1:horzLen+1)
        dLogadmm.objVal[1,p]=objFun(dLogadmm.Sn[:,p],dLogadmm.Xt[:,p],dLogadmm.Un[:,p])
    	fGap= abs(dLogadmm.objVal[1,p]-cSol.objVal)
    	snGap=norm((dLogadmm.Sn[:,p]-cSol.Sn),2)
    	unGap=norm((dLogadmm.Un[:,p]-cSol.Un),2)
    	constGap=norm(dLogadmm.couplConst[:,p],2)
    	itGap = norm(dLogadmm.Lam[:,p]-prevLam[:,1],2)
    	convGap = norm(dLogadmm.Lam[:,p]-cSol.lamCoupl,2)
    	dCMadmm.obj[p,1]=fGap
    	dCMadmm.sn[p,1]=snGap
    	dCMadmm.un[p,1]=unGap
    	dCMadmm.couplConst[p,1]=constGap
    	dCMadmm.lamIt[p,1]=itGap
    	dCMadmm.lam[p,1]=convGap
    	if(convGap <= convChk )
    		@printf "Converged after %g iterations\n" p
    		convIt=p
    		break
    	else
    		@printf "lastGap  %e after %g iterations\n" itGap p
    		@printf "convGap  %e after %g iterations\n" convGap p
    		@printf "constGap %e after %g iterations\n" constGap p
            @printf "snGap    %e after %g iterations\n" snGap p
    		@printf "unGap    %e after %g iterations\n" unGap p
    		@printf("fGap     %e after %g iterations\n\n",fGap,p)
            prevLam=dLogadmm.Lam[:,p]
            prevVz=dLogadmm.Vz[:,p]
            prevVu=dLogadmm.Vu[:,p]
    	end
    end

    return (dLogadmm,dCMadmm,convIt)
end

#ALADIN
function pwlEValad(N::Int,S::Int,horzLen::Int,maxIt::Int,evS::scenarioStruct,cSol::centralSolutionStruct,slack::Bool)
    #initialize
    sn0=evS.s0
    xt0=evS.t0

    stepI = 1
    epsilon = 1e-8
    tolU=1e-4
    tolS=1e-8
    tolT=1e-4
    tolZ=1e-6

    convIt=maxIt

    #ALADIN tuning and initial guess
    σU=1*ones(N,1)
    σS=ones(N,1)/10 #for kA
    #σS=ones(N,1)*100 #for A
    σZ=1/N
    σT=1/10000 #for kA
    #σT=1/10  #for A

    Hu=2*evS.Ri
    Hs=2*evS.Qsi
    # Hu=2*evS.Ri *((1.5-2.5)*rand()+2.5)
    # Hs=2*evS.Qsi *((1.5-2.5)*rand()+2.5)
    Hz=1e-6
    Ht=1e-6

    ρALAD=1
    ρRate=1.15
    muALAD=10^8

    dCMalad=convMetricsStruct()
    dLogalad=itLogPWL()
    convCheck=zeros(maxIt+1,1)

    # lambda0=5*rand(Truncated(Normal(0), 0, 1), horzLen+1)
    # vt0=Tmax*rand(Truncated(Normal(0), 0, 1), (horzLen+1))
    # vz0=ItotalMax/1000*rand(Truncated(Normal(0), 0, 1), S*(horzLen+1))
    # vu0=imax[1,1]*0.8*rand(Truncated(Normal(0), 0, 1), N*(horzLen+1))
    # vs0=rand(Truncated(Normal(0), 0, 1), N*(horzLen+1))
    lambda0=1000*ones(horzLen+1,1)
    vt0=ones(horzLen+1,1)
    vz0=ones(S*(horzLen+1),1)
    vu0=.01*ones(N*(horzLen+1),1)
    vs0=.5*ones(N*(horzLen+1),1)
    # lambda0=lamCurrStar
    # vt0=xtStar
    # vz0=zStar
    # vu0=uStar
    # vs0=snStar



    # dLogalad.Vu[:,1]=vu0
    # dLogalad.Vs[:,1]=vs0
    # dLogalad.Vz[:,1]=vz0
    # dLogalad.Vt[:,1]=vt0
    # dLogalad.Lam[:,1]=lambda0
    prevVu=vu0
    prevVs=vs0
    prevVz=vz0
    prevVt=vt0
    prevLam=lambda0

    ΔY=zeros(1,maxIt+1)
    ρALADp=ρALAD*ones(1,maxIt+1)

    for p=1:maxIt

        #solve decoupled
        @sync @distributed for evInd=1:N
            ind=[evInd]
            for k=1:horzLen
                append!(ind,k*N+evInd)
            end
            evVu=prevVu[ind,1]
            evVs=prevVs[ind,1]
            #evV=zeros(horzLen+1,1)
            target=zeros((horzLen+1),1)
            target[(evS.Kn[evInd,1]-(stepI-1)):1:length(target),1].=evS.Snmin[evInd,1]
            evM = Model(solver = GurobiSolver())
            @variable(evM,sn[1:(horzLen+1)])
            @variable(evM,u[1:(horzLen+1)])
            if slack @variable(evM,slackSn) end
            objExp=sum((sn[k,1]-1)^2*evS.Qsi[evInd,1]+(u[k,1])^2*evS.Ri[evInd,1]+
                                    prevLam[k,1]*(u[k,1])+
                                    ρALADp[1,p]/2*(u[k,1]-evVu[k,1])*σU[evInd,1]*(u[k,1]-evVu[k,1])+
                                    ρALADp[1,p]/2*(sn[k,1]-evVs[k,1])*σS[evInd,1]*(sn[k,1]-evVs[k,1]) for k=1:horzLen+1)
            if slack
                append!(objExp,sum(evS.β[n]*slackSn^2 for n=1:N))
            end
            @objective(evM,Min,objExp)
            @constraint(evM,sn[1,1]==sn0[evInd,1]+evS.ηP[evInd,1]*u[1,1])
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


            TT = stdout # save original stdout stream
            redirect_stdout()
            statusEVM = solve(evM)
            redirect_stdout(TT)
            @assert statusEVM==:Optimal "ALAD EV NLP optimization not solved to optimality"

    		# kappaMax=-getdual(curKappaMax)
    		# kappaMin=-getdual(curKappaMin)
            # socMax=-getdual(socKappaMax)
            # socMin=-getdual(socKappaMin)
            uVal=getvalue(u)
            snVal=getvalue(sn)

            cValMax=abs.(uVal.-evS.imax[evInd,1]).<tolU
            cValMin=abs.(uVal.-evS.imin[evInd,1]).<tolU
            # cVal=kappaMax
            # cVal[cVal.>0]=1
            # cVal=kappaMin
            # cVal[cVal.<0]=-1
            dLogalad.Cu[ind,p]=1cValMax-1cValMin


            cValMax=abs.(snVal.-1).<tolS
            if slack
                cValMin=abs.(snVal.-target*(1-getvalue(slackSn))).<tolS
            else
                cValMin=abs.(snVal.-target).<tolS
            end
            # cVal=socMax
            # cVal[cVal.>0]=1
            # cVal=socMin
            # cVal[cVal.<0]=-1
            dLogalad.Cs[ind,p]=1cValMax-1cValMin

            dLogalad.slackSn[evInd]= if slack getvalue(slackSn) else 0 end
            dLogalad.Sn[ind,p]=snVal
    		dLogalad.Un[ind,p]=uVal
            dLogalad.Gu[ind,p]=2*evS.Ri[evInd,1]*uVal
            dLogalad.Gs[ind,p]=2*evS.Qsi[evInd,1]*snVal.-2*evS.Qsi[evInd,1]
            #Gu[collect(evInd:N:length(Gu[:,p+1])),p+1]=σU[evInd,1]*(evVu-uVal)+lambda
            #Gs[collect(evInd:N:length(Gs[:,p+1])),p+1]=σN[evInd,1]*(evVs-snVal)-lambda
        end

        #N+1 decoupled problem aka transformer current
        tM = Model(solver = GurobiSolver())
        #tM = Model(solver = IpoptSolver())
        @variable(tM,z[1:(S)*(horzLen+1)])
        @variable(tM,xt[1:(horzLen+1)])
        @objective(tM,Min, sum(-prevLam[k,1]*sum(z[(k-1)*(S)+s] for s=1:S)+
                  ρALADp[1,p]/2*σZ*(sum(z[(k-1)*(S)+s] for s=1:S)-sum(prevVz[(k-1)*(S)+s,1] for s=1:S))^2+
                  ρALADp[1,p]/2*σT*(xt[k]-prevVt[k,1])^2  for k=1:(horzLen+1)))
        @constraint(tM,tempCon1,xt[1]-evS.τP*xt0-evS.γP*evS.deltaI*sum((2*s-1)*z[s] for s=1:S)-evS.ρP*evS.Tamb[stepI,1]==0)
        @constraint(tM,tempCon2[k=1:horzLen],xt[k+1]-evS.τP*xt[k]-evS.γP*evS.deltaI*sum((2*s-1)*z[(k)*(S)+s] for s=1:S)-evS.ρP*evS.Tamb[stepI+k,1]==0)
        if noTlimit==false
        	@constraint(tM,upperTCon,xt.<=evS.Tmax)
        end
        @constraint(tM,lowerTCon,xt.>=0)
        @constraint(tM,pwlKappaMin,z.>=0)
        @constraint(tM,pwlKappaMax,z.<=evS.deltaI)
        TT = stdout # save original stdout stream
        redirect_stdout()
        statusTM = solve(tM)
        redirect_stdout(TT)
        @assert statusTM==:Optimal "ALAD XFRM NLP optimization not solved to optimality"

    	#kappaMax=-getdual(pwlKappaMax)
    	#kappaMin=-getdual(pwlKappaMin)
        #tMax=-getdual(upperTCon)
        #tMin=-getdual(lowerTCon)
        lambdaTemp=[-getdual(tempCon1);-getdual(tempCon2)]
        zVal=getvalue(z)
        xtVal=getvalue(xt)

        cValMax=abs.(zVal.-evS.deltaI).<tolZ
        cValMin=abs.(zVal.-0).<tolZ
        # cVal=kappaMax
        # cVal=kappaMin
        # cVal[cVal.<0]=-1
        # cVal[cVal.>0]=1
        dLogalad.Cz[:,p]=1cValMax-1cValMin

        cValMax=abs.(xtVal.-evS.Tmax).<tolT
        cValMin=abs.(xtVal.-0).<tolT
        # cVal=tMin
        # cVal[cVal.<0]=-1
        # cVal=tMax
        # cVal[cVal.>0]=1
        dLogalad.Ct[:,p]=1cValMax-1cValMin

        dLogalad.Xt[:,p]=xtVal
        dLogalad.Z[:,p]=zVal
        #dLogalad.Gz[:,p].=0
        #Gz[:,p+1]=σZ*(Vz[:,p]-zVal)-repeat(-Lam[:,p],inner=S)

        for k=1:horzLen+1
            dLogalad.uSum[k,p]=sum(dLogalad.Un[(k-1)*N+n,p] for n=1:N)
            dLogalad.zSum[k,p]=sum(dLogalad.Z[(k-1)*(S)+s,p] for s=1:S)
            dLogalad.couplConst[k,p]=dLogalad.uSum[k,p] + evS.iD[stepI+(k-1),1] - dLogalad.zSum[k,p]
        end

        #calculate actual temperature from nonlinear model of XFRM
        itotal=zeros(horzLen+1,1)
        for k=1:horzLen+1
            itotal[k,1]=dLogalad.uSum[k,p] + evS.iD[stepI+(k-1),1]
        end
        dLogalad.Tactual[1,p]=evS.τP*xt0+evS.γP*dLogalad.zSum[1,p]^2+evS.ρP*evS.Tamb[stepI,1] #fix for mpc
        for k=1:horzLen
            dLogalad.Tactual[k+1,p]=evS.τP*dLogalad.Tactual[k,p]+evS.γP*dLogalad.zSum[k+1,p]^2+evS.ρP*evS.Tamb[stepI+k,1]  #fix for mpc
        end

        #check for convergence
        constGap=norm(dLogalad.couplConst[:,p],1)
        cc=norm(vcat((prevVu[:,1]-dLogalad.Un[:,p]),(prevVz[:,1]-dLogalad.Z[:,p])),1)
        #convCheck=ρALAD*norm(vcat(repeat(σU,horzLen+1,1).*(Vu[:,p]-Un[:,p+1]),σZ*(Vz[:,p]-Z[:,p+1])),1)
        objFun(sn,xt,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,1]     for n=1:N) for k=1:horzLen+1) +
                        sum((xt[k,1]-1)^2*evS.Qsi[N+1,1]                 for k=1:horzLen+1) +
                        sum(sum((u[(k-1)*N+n,1])^2*evS.Ri[n,1]           for n=1:N) for k=1:horzLen+1)
        dLogalad.objVal[1,p]=objFun(dLogalad.Sn[:,p],dLogalad.Xt[:,p],dLogalad.Un[:,p])
        fGap= abs(dLogalad.objVal[1,p]-cSol.objVal)
        snGap=norm((dLogalad.Sn[:,p]-cSol.Sn),2)
        unGap=norm((dLogalad.Un[:,p]-cSol.Un),2)
        dCMalad.obj[p,1]=fGap
        dCMalad.sn[p,1]=snGap
        dCMalad.un[p,1]=unGap
        dCMalad.couplConst[p,1]=constGap
        convCheck[p,1]=cc
        if  constGap<=epsilon && cc<=epsilon
            @printf "Converged after %g iterations\n" p
            convIt=p
            break
        else
            @printf "convCheck  %e after %g iterations\n" cc p
            @printf "constGap   %e after %g iterations\n" constGap p
            @printf "snGap      %e after %g iterations\n" snGap p
            @printf("fGap       %e after %g iterations\n",fGap,p)
        end

        #coupled QP
        cM = Model(solver = GurobiSolver())
        #cM = Model(solver = IpoptSolver())
        @variable(cM,dUn[1:(N)*(horzLen+1)])
        @variable(cM,dSn[1:(N)*(horzLen+1)])
        @variable(cM,dZ[1:(S)*(horzLen+1)])
        @variable(cM,dXt[1:(horzLen+1)])
        #@variable(cM,relaxS[1:(horzLen+1)])
        # coupledObj(deltaY,Hi,gi)=1/2*deltaY'*Hi*deltaY+gi'*deltaY
    	# objExp=coupledObj(dZ,Hz,Gz[:,p+1])
        # objExp=objExp+coupledObj(dXt,Ht,zeros(length(dXt)))
    	# for n=1:N
        #     objExp=objExp+coupledObj(dUn[collect(n:N:(N)*(horzLen+1)),1],Hu[n,1],Gu[collect(n:N:(N)*(horzLen+1)),p+1])+
        #                   coupledObj(dSn[collect(n:N:(N)*(horzLen+1)),1],Hn[n,1],Gs[collect(n:N:(N)*(horzLen+1)),p+1])
    	# end

        #objExp=objExp+Lam[:,p]'*relaxS+muALAD/2*sum(relaxS[k,1]^2 for k=1:horzLen+1)
    	@objective(cM,Min, sum(sum(0.5*dUn[(k-1)*N+i,1]^2*Hu[i,1]+dLogalad.Gu[(k-1)*N+i,p]*dUn[(k-1)*N+i,1] +
                                   0.5*dSn[(k-1)*N+i,1]^2*Hs[i,1]+dLogalad.Gs[(k-1)*N+i,p]*dSn[(k-1)*N+i,1] for i=1:N) +
                                   0.5*dZ[k,1]^2*Hz+
                                   0.5*dXt[k,1]^2*Ht   for k=1:(horzLen+1)))
        # @constraint(cM,currCon[k=1:horzLen+1],iD[stepI+(k-1)]+relaxS[k,1]==-sum(Un[(k-1)*(N)+n,p+1]+dUn[(k-1)*(N)+n,1] for n=1:N)+
        #                                          sum(Z[(k-1)*(S)+s,p+1]+dZ[(k-1)*(S)+s,1] for s=1:S))
        @constraint(cM,currCon[k=1:horzLen+1],sum(dLogalad.Un[(k-1)*(N)+n,p]+dUn[(k-1)*(N)+n,1] for n=1:N)-
                                                 sum(dLogalad.Z[(k-1)*(S)+s,p]+dZ[(k-1)*(S)+s,1] for s=1:S)==-evS.iD[stepI+(k-1)])#+relaxS[k,1])
        #local equality constraints C*(X+deltaX)=0 is same as C*deltaX=0 since we already know CX=0
        @constraint(cM,stateCon1[n=1:N],dSn[n,1]==evS.ηP[n,1]*dUn[n,1])
        @constraint(cM,stateCon2[k=1:horzLen,n=1:N],dSn[n+(k)*(N),1]==dSn[n+(k-1)*(N),1]+evS.ηP[n,1]*dUn[n+(k)*(N),1])
        @constraint(cM,tempCon1,dXt[1,1]==evS.γP*evS.deltaI*sum((2*s-1)*dZ[s,1] for s=1:S))
        @constraint(cM,tempCon2[k=1:horzLen],dXt[k+1,1]==evS.τP*dXt[k,1]+evS.γP*evS.deltaI*sum((2*s-1)*dZ[k*S+s,1] for s=1:S))

        #local inequality constraints
        # @constraint(cM,(Z[:,p+1]+dZ).>=0)
        # @constraint(cM,(Z[:,p+1]+dZ).<=deltaI)
     	# @constraint(cM,(Un[:,p+1]+dUn).<=repeat(imax,horzLen+1,1))
     	# @constraint(cM,(Un[:,p+1]+dUn).>=repeat(imin,horzLen+1,1))

        #these shouldnt be elementwise?????
        # @constraint(cM,Cz[:,p+1].*dZ.==0)
        # @constraint(cM,Cu[:,p+1].*dUn.==0)
        # @constraint(cM,Cn[:,p+1].*dSn.==0)
        # @constraint(cM,Ct[:,p+1].*dXt.==0)

        @constraint(cM,dLogalad.Cz[:,p].*dZ.<=0)
        @constraint(cM,dLogalad.Cu[:,p].*dUn.<=0)
        @constraint(cM,dLogalad.Cs[:,p].*dSn.<=0)
        @constraint(cM,dLogalad.Ct[:,p].*dXt.<=0)


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
        #alpha3=alpha3/ceil(p/2)
        #Lam[:,p+1]=Lam[:,p]+alpha3*(-getdual(currCon)-Lam[:,p])

        dLogalad.Lam[:,p]=max.(prevLam[:,1]+α3*(-getdual(currCon)-prevLam[:,1]),0)
        dLogalad.Vu[:,p]=prevVu[:,1]+α1*(dLogalad.Un[:,p]-prevVu[:,1])+α2*getvalue(dUn)
        dLogalad.Vz[:,p]=prevVz[:,1]+α1*(dLogalad.Z[:,p]-prevVz[:,1])+α2*getvalue(dZ)
        dLogalad.Vs[:,p]=prevVs[:,1]+α1*(dLogalad.Sn[:,p]-prevVs[:,1])+α2*getvalue(dSn)
        dLogalad.Vt[:,p]=prevVt[:,1]+α1*(dLogalad.Xt[:,p]-prevVt[:,1])+α2*getvalue(dXt)

        dCMalad.lamIt[p,1]=norm(dLogalad.Lam[:,p]-prevLam[:,1],2)
        dCMalad.lam[p,1]=norm(dLogalad.Lam[:,p]-cSol.lamCoupl,2)
        @printf "lastGap    %e after %g iterations\n" dCMalad.lamIt[p,1] p
        @printf "convLamGap %e after %g iterations\n\n" dCMalad.lam[p,1] p

        ρALADp[1,p+1]=min(ρALADp[1,p]*ρRate,1e6) #increase ρ every iteration
        ΔY[1,p]=norm(vcat(getvalue(dUn),getvalue(dZ),getvalue(dSn),getvalue(dXt)),Inf)

        #reset for next iteration
        prevVu=dLogalad.Vu[:,p]
        prevVs=dLogalad.Vs[:,p]
        prevVz=dLogalad.Vz[:,p]
        prevVt=dLogalad.Vt[:,p]
        prevLam=dLogalad.Lam[:,p]
    end

    return dLogalad,dCMalad,convIt,ΔY,convCheck
end
