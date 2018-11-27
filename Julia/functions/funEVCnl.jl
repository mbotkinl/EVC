#functions to run NL EVC problems


#central
function nlEVcentral(N::Int,S::Int,horzLen::Int,evS::scenarioStruct,forecastError::Bool,relaxedMode=2,slack=false)

    #initialize with current states
    sn0=evS.s0
    xt0=evS.t0
    stepI=1
	if forecastError
		iD=evS.iDnoise
	else
		iD=evS.iD
	end

    #desired SOC
    target=zeros(N*(horzLen+1),1)
    for ii=1:N
       cur=evS.Kn[ii]-(stepI-1)
       ind=(max(0,(cur-1)*N)+ii):N:length(target)
       target[ind].=evS.Snmin[ii,1]
    end

    println("setting up model")
    if relaxedMode==2
		cModel = Model(solver = MosekSolver())
	elseif relaxedMode==1
		cModel = Model(solver = GurobiSolver(QCPDual=1))
    else
        cModel = Model(solver = IpoptSolver())
    end

    #u iD and z are one index ahead of sn and T. i.e the x[k+1]=x[k]+eta*u[k+1]
    @variable(cModel,u[1:N*(horzLen+1)])
    @variable(cModel,sn[1:(N)*(horzLen+1)])
    @variable(cModel,xt[1:(horzLen+1)])
    @variable(cModel,itotal[1:(horzLen+1)])

	#objExp=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,1]+(u[(k-1)*N+n,1])^2*evS.Ri[n,1] for n=1:N) for k=1:(horzLen+1))

    println("constraints")

	if relaxedMode==1
		objExp =sum((sn[n,1]-1)^2*evS.Qsi[n,1]+(u[n,1])^2*evS.Ri[n,1] for n=1:N)
		for k=2:horzLen+1
			append!(objExp,sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,1]+(u[(k-1)*N+n,1])^2*evS.Ri[n,1]  for n=1:N))
		end
		@constraint(cModel,tempCon1,xt[1,1]>=evS.τP*xt0+evS.γP*(itotal[1])^2+evS.ρP*evS.Tamb[stepI,1])
		@constraint(cModel,tempCon2[k=1:horzLen],xt[k+1,1]>=evS.τP*xt[k,1]+evS.γP*(itotal[k+1])^2+evS.ρP*evS.Tamb[stepI+k,1])
	elseif relaxedMode==2
		@variable(cModel,e[1:(horzLen+1)])
		@variable(cModel,t)
		objExp=t+sum(-2*evS.Qsi[n,1]*sn[n,1]+evS.Qsi[n,1] for n=1:N)
		for k=2:horzLen+1
			append!(objExp,sum(-2*evS.Qsi[n,1]*sn[(k-1)*(N)+n,1]+evS.Qsi[n,1]  for n=1:N))
		end
		@constraint(cModel,tempCon1,xt[1,1]>=evS.τP*xt0+evS.γP*e[1]+evS.ρP*evS.Tamb[stepI,1])
		@constraint(cModel,tempCon2[k=1:horzLen],xt[k+1,1]>=evS.τP*xt[k,1]+evS.γP*e[k+1]+evS.ρP*evS.Tamb[stepI+k,1])
        #@constraint(cModel,eCon[k=1:horzLen+1],(itotal[k])^2-e[k]<=0)

		#@constraint(cModel,eCon[k=1:horzLen+1],norm([1/2-1/2*e[k] itotal[k] 0])<=1/2+1/2*e[k])
		@constraint(cModel,eCon[k=1:horzLen+1],norm([2itotal[k] e[k]-1])<=e[k]+1)

		#@constraint(cModel,eCon[k=1:horzLen+1],norm(itotal[k])<=e[k])

		#reformulated objective function (do this smarter***)
		objExpCon=0*sn[1:2]
		for k=1:horzLen+1
			for n=1:N
				tt=[2*sqrt(evS.Qsi[n,1])*sn[(k-1)*(N)+n,1];2*sqrt(evS.Ri[n,1])*u[(k-1)*(N)+n,1]]
				append!(objExpCon,tt)
			end
		end
		append!(objExpCon,[t-1])

		# objExpCon1 = [sqrt(evS.Qsi[n,1])*sn[(k-1)*(N)+n,1]-evS.Qsi[n,1]/sqrt(evS.Qsi[n,1])  for n=1:N for k=1:(horzLen+1)]
		# objExpCon2 = [sqrt(evS.Ri[n,1])*u[(k-1)*(N)+n,1]  for n=1:N for k=1:(horzLen+1)]
		# objExpCon=vcat(objExpCon1,objExpCon2)

		@constraint(cModel,objCon,norm(objExpCon)<=t+1)
    else
		objExp =sum((sn[n,1]-1)^2*evS.Qsi[n,1]+(u[n,1])^2*evS.Ri[n,1] for n=1:N)
		for k=2:horzLen+1
			append!(objExp,sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,1]+(u[(k-1)*N+n,1])^2*evS.Ri[n,1]  for n=1:N))
		end
        @NLconstraint(cModel,tempCon1,xt[1,1]==evS.τP*xt0+evS.γP*(itotal[1])^2+evS.ρP*evS.Tamb[stepI,1])
        @NLconstraint(cModel,tempCon2[k=1:horzLen],xt[k+1,1]==evS.τP*xt[k,1]+evS.γP*(itotal[k+1])^2+evS.ρP*evS.Tamb[stepI+k,1])
    end
	if slack
		@variable(cModel,slackSn[1:N])
		append!(objExp,sum(evS.β[n]*slackSn[n]^2 for n=1:N))
	end

	@objective(cModel,Min, objExp)

    @constraint(cModel,stateCon1,sn[1:N,1].==sn0[1:N,1]+evS.ηP[:,1].*u[1:N,1])
    @constraint(cModel,stateCon2[k=1:horzLen,n=1:N],sn[n+(k)*(N),1]==sn[n+(k-1)*(N),1]+evS.ηP[n,1]*u[n+(k)*(N),1])
    @constraint(cModel,currCon[k=1:horzLen+1],0==-sum(u[(k-1)*(N)+n,1] for n=1:N)-iD[stepI+(k-1)]+itotal[k])
    @constraint(cModel,sn.<=1)
	if slack
		@constraint(cModel,sn.>=target.*(1-repeat(slackSn,horzLen+1,1)))
		@constraint(cModel,slackSn.>=0)
	else
		@constraint(cModel,sn.>=target)
	end
    if noTlimit==false
    	@constraint(cModel,upperTCon,xt.<=evS.Tmax)
    end
    @constraint(cModel,xt.>=0)
    @constraint(cModel,upperCCon,u.<=repeat(evS.imax,horzLen+1,1))
    @constraint(cModel,u.>=repeat(evS.imin,horzLen+1,1))
    @constraint(cModel,itotal.<=evS.ItotalMax)
    @constraint(cModel,itotal.>=0)

    println("solving....")
    statusC = solve(cModel)
    @assert statusC==:Optimal "Central NL optimization not solved to optimality"

    uRaw=getvalue(u)
    snRaw=getvalue(sn)
    xtRaw=getvalue(xt)
    itotalRaw=getvalue(itotal)
	if slack
        slackSnRaw=getvalue(slackSn)
    else
        slackSnRaw=zeros(N)
    end

    if noTlimit==false
    	kappaUpperT=-getdual(upperTCon)
    else
    	kappaUpperT=zeros(horzLen+1,1)
    end
    lambdaCurr=-getdual(currCon)
    if relaxedMode==1
    	lambdaTemp=-Gurobi.get_dblattrarray(getrawsolver(cModel),"QCPi",1,horzLen+1)
    else
        lambdaTemp=[-getdual(tempCon1);-getdual(tempCon2)]
    end

    uSum=zeros(horzLen+1,1)
    for k=1:horzLen+1
        uSum[k,1]=sum(uRaw[(k-1)*N+n,1] for n=1:N)
    end

    cSol=centralSolutionStruct(Xt=xtRaw,Un=uRaw,Sn=snRaw,
                        Itotal=itotalRaw,uSum=uSum,objVal=getobjectivevalue(cModel),
                        lamTemp=lambdaTemp,lamCoupl=lambdaCurr,slackSn=slackSnRaw)
    return cSol
end

#dual
function nlEVdual(N::Int,S::Int,horzLen::Int,maxIt::Int,updateMethod::String,
    evS::scenarioStruct,cSolnl::centralSolutionStruct,forecastError::Bool, relaxedMode=2,slack=false)

    #initialize with current states
    sn0=evS.s0
    xt0=evS.t0
	if forecastError
		iD=evS.iDnoise
	else
		iD=evS.iD
	end


    if updateMethod=="fastAscent"
    	#alpha = 0.1 #for A
    	alpha = 5e3 #for kA
    	alphaDivRate=2
    	minAlpha=1e-6
    	#alphaRate=.99
    else
    	#alpha = .01 #for A
    	alpha = 5e3 #for kA
    	alphaDivRate=2
    	minAlpha=1e-6
    	#alphaRate=.99
    end

    stepI = 1
    convChk = 1e-8
    convIt=maxIt

    alphaP=alpha*ones(maxIt+1,1)

    dCM=convMetricsStruct()
    dLog=itLogNL()

    #u iD and z are one index ahead of sn and T. i.e the x[k+1]=x[k]+eta*u[k+1]
    lambda0=2e3*ones(horzLen+1,1)
    #lambda0=lamCurrStarNL
    #lambda0=max.(lamCurrStarNL,0)

    #dLog.Lam[:,1]=lambda0
	prevLam=lambda0

    #iterate at each time step until convergence
    for p=1:maxIt
        #solve subproblem for each EV
    	@sync @distributed for evInd=1:N
    		ind=[evInd]
    		for k=1:horzLen
    			append!(ind,k*N+evInd)
    		end
            target=zeros((horzLen+1),1)
    		target[(evS.Kn[evInd,1]-(stepI-1)):1:length(target),1].=evS.Snmin[evInd,1]
			if relaxedMode==2
				evM = Model(solver = MosekSolver())
			elseif relaxedMode==1
				evM = Model(solver = GurobiSolver())
		    else
		        evM = Model(solver = IpoptSolver())
		    end
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

    		@assert statusEVM==:Optimal "EV NL optimization not solved to optimality"

            dLog.Sn[ind,p]=getvalue(sn)
            dLog.Un[ind,p]=getvalue(un)
			dLog.slackSn[evInd]= if slack getvalue(slackSn) else 0 end
        end

    	if updateMethod=="dualAscent"
    	    #solve coordinator problem

			if relaxedMode==2
				coorM = Model(solver = MosekSolver())
			elseif relaxedMode==1
				coorM = Model(solver = GurobiSolver())
		    else
		        coorM = Model(solver = IpoptSolver())
		    end
	        @variable(coorM,itotal[1:(horzLen+1)])
    	    @variable(coorM,xt[1:(horzLen+1)])
    	    @objective(coorM,Min,-sum(prevLam[k,1]*itotal[k,1] for k=1:(horzLen+1)))
			if relaxedMode ==1
                @constraint(coorM,xt[1,1]>=evS.τP*xt0+evS.γP*(itotal[1])^2+evS.ρP*evS.Tamb[stepI,1]) #fix for MPC loop
        		@constraint(coorM,[k=1:horzLen],xt[k+1,1]>=evS.τP*xt[k,1]+evS.γP*(itotal[k+1,1])^2+evS.ρP*evS.Tamb[stepI+k,1])
			elseif relaxedMode==2
				@variable(coorM,e[1:(horzLen+1)])
				@constraint(coorM,tempCon1,xt[1,1]>=evS.τP*xt0+evS.γP*e[1]+evS.ρP*evS.Tamb[stepI,1])
				@constraint(coorM,tempCon2[k=1:horzLen],xt[k+1,1]>=evS.τP*xt[k,1]+evS.γP*e[k+1]+evS.ρP*evS.Tamb[stepI+k,1])
				#@constraint(coorM,eCon[k=1:horzLen+1],norm([1/2-1/2*e[k] itotal[k]])<=1/2+1/2*e[k])
				@constraint(coorM,eCon[k=1:horzLen+1],norm([2itotal[k] e[k]-1])<=e[k]+1)
            else
                @NLconstraint(coorM,xt[1,1]==evS.τP*xt0+evS.γP*(itotal[1])^2+evS.ρP*evS.Tamb[stepI,1]) #fix for MPC loop
        		@NLconstraint(coorM,[k=1:horzLen],xt[k+1,1]==evS.τP*xt[k,1]+evS.γP*(itotal[k+1,1])^2+evS.ρP*evS.Tamb[stepI+k,1])
            end
    		if noTlimit==false
    			@constraint(coorM,upperTCon,xt.<=evS.Tmax)
    		end
    	    @constraint(coorM,xt.>=0)
    	    @constraint(coorM,itotal.<=evS.ItotalMax)
    	    @constraint(coorM,itotal.>=0)
    		TT = stdout # save original stdout strea
    		redirect_stdout()
    	    statusC = solve(coorM)
    		redirect_stdout(TT)

    		@assert statusC==:Optimal "XFRM NL optimization not solved to optimality"

    		dLog.Xt[:,p]=getvalue(xt)
    		dLog.Itotal[:,p]=getvalue(itotal)

    	    #grad of lagragian
    		gradL=zeros(horzLen+1,1)
    		for k=1:horzLen+1
    			dLog.uSum[k,p]=sum(dLog.Un[(k-1)*N+n,p] for n=1:N)
    			dLog.couplConst[k,p]=dLog.uSum[k,p] + iD[stepI+(k-1),1] - dLog.Itotal[k,p]
    		end
    		dCM.couplConst[p,1]=norm(dLog.couplConst[:,p],2)
    	end

    	if updateMethod=="fastAscent"
            ztotal=zeros(horzLen+1,1)
            for k=1:horzLen+1
                ztotal[k,1]=sum(dLog.Un[(k-1)*N+n,p] for n=1:N) + iD[stepI+(k-1),1]
            end
            dLog.Xt[1,p]=evS.τP*xt0+evS.γP*ztotal[1,1]^2+evS.ρP*evS.Tamb[stepI,1]
            for k=1:horzLen
                dLog.Xt[k+1,p]=evS.τP*dLog.Xt[k,p]+evS.γP*ztotal[k+1,1]^2+evS.ρP*evS.Tamb[stepI+k,1]
            end

    		#fast ascent
    		if noTlimit==false
    			dLog.couplConst[:,p]=dLog.Xt[:,p]-evS.Tmax*ones(horzLen+1,1)
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
    	alphaP = max(alpha/ceil(p/alphaDivRate),minAlpha)
    	#alphaP = alphaP*alphaRate
		dLog.itUpdate[1,p]=alphaP

		dLog.Lam[:,p]=prevLam[:,1]+alphaP*dLog.couplConst[:,p]
        #dLog.Lam[:,p]=max.(prevLam[:,1]+alphaP[p,1]*dLog.couplConst[:,p],0)

    	#check convergence
    	objFun(sn,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,1]     for n=1:N) for k=1:horzLen+1) +
    					sum(sum((u[(k-1)*N+n,1])^2*evS.Ri[n,1]           for n=1:N) for k=1:horzLen+1)
        dLog.objVal[1,p]=objFun(dLog.Sn[:,p],dLog.Un[:,p])
    	fGap=abs(dLog.objVal[1,p]-cSolnl.objVal)
    	snGap=norm((dLog.Sn[:,p]-cSolnl.Sn),2)
    	unGap=norm((dLog.Un[:,p]-cSolnl.Un),2)
    	itGap = norm(dLog.Lam[:,p]-prevLam[:,1],2)
    	if updateMethod=="fastAscent"
    		convGap = norm(dLog.Lam[:,p]-cSolnl.lamTemp,2)
    	else
    		convGap = norm(dLog.Lam[:,p]-cSolnl.lamCoupl,2)
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
    return (dLog,dCM,convIt)
end

#admm
function nlEVadmm(N::Int,S::Int,horzLen::Int,maxIt::Int,evS::scenarioStruct,cSolnl::centralSolutionStruct,
	forecastError::Bool,relaxedMode=false, slack=false)
    #initialize with current states
    sn0=evS.s0
    xt0=evS.t0
	if forecastError
		iD=evS.iDnoise
	else
		iD=evS.iD
	end


    stepI = 1
    convChk = 1e-8
    convIt=maxIt

    #admm  initial parameters and guesses
    #ρADMM=10.0^(0)
    ρADMM=1e5
    # ρDivRate=10
	ρDivRate=1.1
	maxRho=1e7

    #u iD and z are one index ahead of sn and T. i.e the x[k+1]=x[k]+eta*u[k+1]
    dCMadmm=convMetricsStruct()
    dLogadmm=itLogNL()

    # lambda0=lamCurrStarNL
    # vi0=-itotalStarNL
    # vu0=uStarNL
    lambda0=2e3*ones(horzLen+1,1)
    vi0=-ones(horzLen+1,1)
    vu0=.01*ones(N*(horzLen+1),1)

    # dLogadmm.Vi[:,1]=vi0
    # dLogadmm.Vu[:,1]=vu0
    # dLogadmm.Lam[:,1]=lambda0
	prevLam=lambda0
	prevVu=vu0
	prevVi=vi0
	ρADMMp = ρADMM

    for p in 1:maxIt
        #x minimization eq 7.66 in Bertsekas
        @sync @distributed for evInd=1:N
            evV=prevVu[collect(evInd:N:length(prevVu[:,1])),1]
            target=zeros((horzLen+1),1)
    		target[(evS.Kn[evInd,1]-(stepI-1)):1:length(target),1].=evS.Snmin[evInd,1]
			if relaxedMode==2
				evM = Model(solver = MosekSolver())
			elseif relaxedMode==1
				evM = Model(solver = GurobiSolver(NumericalFocus=3))
		    else
		        evM = Model(solver = IpoptSolver())
		    end
        	@variable(evM,sn[1:(horzLen+1)])
        	@variable(evM,u[1:(horzLen+1)])
			if slack @variable(evM,slackSn) end
			objExp=sum((sn[k,1]-1)^2*evS.Qsi[evInd,1]+(u[k,1])^2*evS.Ri[evInd,1]+
                                    prevLam[k,1]*(u[k,1]-evV[k,1])+
                                    ρADMMp/2*(u[k,1]-evV[k,1])^2 for k=1:horzLen+1)
			if slack
                append!(objExp,sum(evS.β[n]*slackSn^2 for n=1:N))
            end
			@objective(evM,Min, objExp)
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
            @assert statusEVM==:Optimal "EV NLP NL optimization not solved to optimality"

    		dLogadmm.Sn[collect(evInd:N:length(dLogadmm.Sn[:,p])),p]=getvalue(sn)
    		dLogadmm.Un[collect(evInd:N:length(dLogadmm.Un[:,p])),p]=getvalue(u)
			dLogadmm.slackSn[evInd]= if slack getvalue(slackSn) else 0 end
		end

        #N+1 decoupled problem aka transformer current
		if relaxedMode==2
			tM = Model(solver = MosekSolver())
		elseif relaxedMode==1
			tM = Model(solver = GurobiSolver())
	    else
	        tM = Model(solver = IpoptSolver())
	    end
        #@variable(tM,z[1:(S)*(horzLen+1)])
        @variable(tM,xt[1:(horzLen+1)])
        @variable(tM,itotal[1:(horzLen+1)])
    	# constFun1(u,v)=sum(Lam[k,p]*(u[k,1]-v[k,1])  for k=1:(horzLen+1))
    	# constFun2(u,v)=ρADMM/2*sum((u[k,1]-v[k,1])^2  for k=1:(horzLen+1))
        # @objective(tM,Min, constFun1(-itotal,Vi[:,p])+constFun2(-itotal,Vi[:,p]))


		if relaxedMode ==1
			@objective(tM,Min,sum(prevLam[k,1]*(-itotal[k,1]-prevVi[k,1])+
								  ρADMM/2*(-itotal[k,1]-prevVi[k,1])^2  for k=1:(horzLen+1)))
			@constraint(tM,tempCon1,xt[1,1]>=evS.τP*xt0+evS.γP*(itotal[1,1])^2+evS.ρP*evS.Tamb[stepI,1])
			@constraint(tM,tempCon2[k=1:horzLen],xt[k+1,1]>=evS.τP*xt[k,1]+evS.γP*(itotal[k+1,1])^2+evS.ρP*evS.Tamb[stepI+k,1])
		elseif relaxedMode==2
			@variable(tM,e[1:(horzLen+1)])
			@variable(tM,t)

			objExp=t+sum(prevLam[k,1]*(-itotal[k,1]-prevVi[k,1])+
								  ρADMM/2*(2*prevVi[k,1]*itotal[k,1]+(prevVi[k,1])^2)  for k=1:(horzLen+1))
			@objective(tM,Min,objExp)

			@constraint(tM,tempCon1,xt[1,1]>=evS.τP*xt0+evS.γP*e[1]+evS.ρP*evS.Tamb[stepI,1])
			@constraint(tM,tempCon2[k=1:horzLen],xt[k+1,1]>=evS.τP*xt[k,1]+evS.γP*e[k+1]+evS.ρP*evS.Tamb[stepI+k,1])
			#@constraint(tM,eCon[k=1:horzLen+1],norm([1/2-1/2*e[k] itotal[k]])<=1/2+1/2*e[k])
			@constraint(tM,eCon[k=1:horzLen+1],norm([2itotal[k] e[k]-1])<=e[k]+1)

			objExpCon = [2*sqrt(ρADMM/2)*itotal[k,1] for k=1:(horzLen+1)]
			append!(objExpCon,[t-1])
			@constraint(tM,objCon,norm(objExpCon)<=t+1)
        else
			@objective(tM,Min,sum(prevLam[k,1]*(-itotal[k,1]-prevVi[k,1])+
								  ρADMM/2*(-itotal[k,1]-prevVi[k,1])^2  for k=1:(horzLen+1)))
            @NLconstraint(tM,tempCon1,xt[1,1]==evS.τP*xt0+evS.γP*(itotal[1,1])^2+evS.ρP*evS.Tamb[stepI,1])
            @NLconstraint(tM,tempCon2[k=1:horzLen],xt[k+1,1]==evS.τP*xt[k,1]+evS.γP*(itotal[k+1,1])^2+evS.ρP*evS.Tamb[stepI+k,1])
        end
        if noTlimit==false
        	@constraint(tM,upperTCon,xt.<=evS.Tmax)
        end
        @constraint(tM,xt.>=0)
        @constraint(tM,itotal.>=0)
        @constraint(tM,itotal.<=evS.ItotalMax)

        TT = stdout # save original stdout stream
        redirect_stdout()
        statusTM = solve(tM)
        redirect_stdout(TT)
        @assert statusTM==:Optimal "XFRM NLP NL optimization not solved to optimality"

        dLogadmm.Xt[:,p]=getvalue(xt)
        dLogadmm.Itotal[:,p]=getvalue(itotal)

        #lambda update eq 7.68
    	for k=1:horzLen+1
    		dLogadmm.uSum[k,p]=sum(dLogadmm.Un[(k-1)*N+n,p] for n=1:N)
    		dLogadmm.couplConst[k,p]=dLogadmm.uSum[k,p] + iD[stepI+(k-1),1] - dLogadmm.Itotal[k,p]
            #dLogadmm.Lam[k,p]=max.(prevLam[k,1]+ρADMMp/(N+1)*(dLogadmm.couplConst[k,p]),0)
    		dLogadmm.Lam[k,p]=prevLam[k,1]+ρADMMp/(N+1)*(dLogadmm.couplConst[k,p])
    	end

        #v upate eq 7.67
        for k=1:horzLen+1
            # dLogadmm.Vu[(k-1)*N.+collect(1:N),p]=min.(max.(dLogadmm.Un[(k-1)*N.+collect(1:N),p].+(prevLam[k,1]-dLogadmm.Lam[k,p])/ρADMMp,evS.imin),evS.imax)
    		# dLogadmm.Vi[k,p]=max.(min.(-dLogadmm.Itotal[k,p]+(prevLam[k,1]-dLogadmm.Lam[k,p])/ρADMMp,0),-evS.ItotalMax)
			dLogadmm.Vu[(k-1)*N.+collect(1:N),p]=dLogadmm.Un[(k-1)*N.+collect(1:N),p].+(prevLam[k,1]-dLogadmm.Lam[k,p])/ρADMMp
			dLogadmm.Vi[k,p]=-dLogadmm.Itotal[k,p]+(prevLam[k,1]-dLogadmm.Lam[k,p])/ρADMMp
        end

		#update rho
		#ρADMMp = ρADMM/ceil(p/ρDivRate)
		dLogadmm.itUpdate[1,p]= min(ρADMMp*ρDivRate,maxRho)

        #check convergence
    	objFun(sn,xt,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,1]     for n=1:N) for k=1:horzLen+1) +
    					sum((xt[k,1]-1)^2*evS.Qsi[N+1,1]                 for k=1:horzLen+1) +
    					sum(sum((u[(k-1)*N+n,1])^2*evS.Ri[n,1]           for n=1:N) for k=1:horzLen+1)
        dLogadmm.objVal[1,p]=objFun(dLogadmm.Sn[:,p],dLogadmm.Xt[:,p],dLogadmm.Un[:,p])
    	fGap= dLogadmm.objVal[1,p]-cSolnl.objVal
    	#fGap= objFun(Sn[:,p],Xt[:,p],Un[:,p])-fStar
    	snGap=norm((dLogadmm.Sn[:,p]-cSolnl.Sn),2)
        unGap=norm((dLogadmm.Un[:,p]-cSolnl.Un),2)
    	constGap=norm(dLogadmm.couplConst[:,p],2)
    	itGap = norm(dLogadmm.Lam[:,p]-prevLam[:,1],2)
    	convGap = norm(dLogadmm.Lam[:,p]-cSolnl.lamCoupl,2)
    	dCMadmm.obj[p,1]=abs(fGap)
    	dCMadmm.sn[p,1]=snGap
        dCMadmm.un[p,1]=unGap
    	dCMadmm.couplConst[p,1]=constGap
    	dCMadmm.lamIt[p,1]=itGap
    	dCMadmm.lam[p,1]=convGap
    	if(itGap <= convChk )
    		@printf "Converged after %g iterations\n" p
    		convIt=p
    		break
    	else
    		@printf "lastGap  %e after %g iterations\n" itGap p
    		@printf "convGap  %e after %g iterations\n" convGap p
    		@printf "constGap %e after %g iterations\n" constGap p
            @printf "snGap    %e after %g iterations\n" snGap p
    		@printf("fGap     %e after %g iterations\n\n",fGap,p)
			prevLam=dLogadmm.Lam[:,p]
			prevVu=dLogadmm.Vu[:,p]
			prevVi=dLogadmm.Vi[:,p]
			ρADMMp=dLogadmm.itUpdate[1,p]
    	end
    end

    return (dLogadmm,dCMadmm,convIt)
end

#aladin
function nlEValad(N::Int,S::Int,horzLen::Int,maxIt::Int,evS::scenarioStruct,cSolnl::centralSolutionStruct,
	forecastError::Bool,relaxedMode::Int, slack=false,eqForm=true)
    #initialize with current states
    sn0=evS.s0
    xt0=evS.t0
	if forecastError
		iD=evS.iDnoise
	else
		iD=evS.iD
	end


    stepI = 1
    epsilon = 1e-12
    tolU=1e-4
    tolS=1e-8
    tolT=1e-4
    tolI=1e-6
    convIt=maxIt

    #ALADIN tuning and initial guess
    ##σs are tuned to i_n [.010kA]
    σU=1*ones(N,1)
    σS=ones(N,1)/1e2
    σI=1
    σT=1e-3

	# σU=1*ones(N,1)
	# σS=ones(N,1)/10
	# σI=1e-2
	# σT=1e-4


    Hu=2*evS.Ri#*(1+rand())
    Hs=2*evS.Qsi#*(1+rand())
    # Hi=1e-6
    # Ht=1e-6
	Hi=0
	Ht=0

	ρALAD=1e3
    ρRate=1.1
    ρALADmax=1e7

    μALAD=1e8
    μRate=1
    μALADmax=2e9

    #.1/.1 looks good expect for const
    #.01/2 increasing ρ helps const
    #.1/1.1 best so far???

    dCMalad=convMetricsStruct()
    dLogalad=itLogNL()
    convCheck=zeros(maxIt+1,1)

    lambda0=2e3*ones(horzLen+1,1)
    vt0=ones(horzLen+1,1)
    vi0=ones((horzLen+1),1)
    vu0=.01*ones(N*(horzLen+1),1)
    vs0=.5*ones(N*(horzLen+1),1)
    # lambda0=lamCurrStarNL
    # vt0=xtStarNL
    # vi0=itotalStarNL
    # vu0=uStarNL
    # vs0=snStarNL

	prevVu=vu0
    prevVs=vs0
    prevVi=vi0
    prevVt=vt0
    prevLam=lambda0
	ρALADp=ρALAD
    μALADp=μALAD

    ΔY=zeros(1,maxIt+1)

    for p=1:maxIt
        @printf "Starting iteration %g \n" p

        #solve decoupled
        @sync @distributed for evInd=1:N
            ind=[evInd]
            for k=1:horzLen
                append!(ind,k*N+evInd)
            end
            evVu=prevVu[ind,1]
            evVs=prevVs[ind,1]
            target=zeros((horzLen+1),1)
            target[(evS.Kn[evInd,1]-(stepI-1)):1:length(target),1].=evS.Snmin[evInd,1]
            #evM = Model(solver = GurobiSolver(NumericFocus=3))
			if relaxedMode==2
				evM = Model(solver = MosekSolver())
			elseif relaxedMode==1
				evM = Model(solver = GurobiSolver())
		    else
		        evM = Model(solver = IpoptSolver())
		    end
            @variable(evM,sn[1:(horzLen+1)])
            @variable(evM,u[1:(horzLen+1)])
			if slack @variable(evM,slackSn) end
			objExp=sum((sn[k,1]-1)^2*evS.Qsi[evInd,1]+(u[k,1])^2*evS.Ri[evInd,1]+
                                    prevLam[k,1]*(u[k,1])+
                                    ρALADp/2*(u[k,1]-evVu[k,1])*σU[evInd,1]*(u[k,1]-evVu[k,1])+
                                    ρALADp/2*(sn[k,1]-evVs[k,1])*σS[evInd,1]*(sn[k,1]-evVs[k,1]) for k=1:horzLen+1)
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
            @assert statusEVM==:Optimal "ALAD EV NLP NL optimization not solved to optimality"

    		# kappaMax=-getdual(curKappaMax)
    		# kappaMin=-getdual(curKappaMin)
            # socMax=-getdual(socKappaMax)
            # socMin=-getdual(socKappaMin)
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
            dLogalad.Csu[ind,p+1]=1cValMax
			dLogalad.Csl[ind,p+1]=-1cValMin

			dLogalad.slackSn[evInd]= if slack getvalue(slackSn) else 0 end
            dLogalad.Sn[ind,p]=snVal
    		dLogalad.Un[ind,p]=uVal

			# dLogalad.Gu[ind,p]=2*evS.Ri[evInd,1]*uVal
			# dLogalad.Gs[ind,p]=2*evS.Qsi[evInd,1]*snVal.-2*evS.Qsi[evInd,1]

			dLogalad.Gu[ind,p]=round.(2*evS.Ri[evInd,1]*uVal,digits=4) # these might matter....
			dLogalad.Gs[ind,p]=round.(2*evS.Qsi[evInd,1]*snVal.-2*evS.Qsi[evInd,1],digits=8)

            #Gu[collect(evInd:N:length(Gu[:,p+1])),p+1]=σU[evInd,1]*(evVu-uVal)+lambda
            #Gs[collect(evInd:N:length(Gs[:,p+1])),p+1]=σN[evInd,1]*(evVs-snVal)-lambda
        end

        #N+1 decoupled problem aka transformer current
		if relaxedMode==2
			tM = Model(solver = MosekSolver())
		elseif relaxedMode==1
			tM = Model(solver = GurobiSolver(QCPDual=1))
	    else
	        tM = Model(solver = IpoptSolver())
	    end
        @variable(tM,itotal[1:(horzLen+1)])
        @variable(tM,xt[1:(horzLen+1)])
		if relaxedMode ==1
			@objective(tM,Min, sum(-prevLam[k,1]*itotal[k]+
					  ρALADp/2*σI*(itotal[k]-prevVi[k,1])^2+
					  ρALADp/2*σT*(xt[k]-prevVt[k,1])^2  for k=1:(horzLen+1)))
		    @constraint(tM,tempCon1,xt[1]-evS.τP*xt0-evS.γP*(itotal[1])^2-evS.ρP*evS.Tamb[stepI,1]>=0)
		  	@constraint(tM,tempCon2[k=1:horzLen],xt[k+1]-evS.τP*xt[k]-evS.γP*(itotal[k+1])^2-evS.ρP*evS.Tamb[stepI+k,1]>=0)
		elseif relaxedMode==2
			@variable(tM,e[1:(horzLen+1)])
			@variable(tM,t)

			@objective(tM,Min, sum(-prevLam[k,1]*itotal[k]+
					  ρALADp/2*σI*(-2itotal[k]*prevVi[k,1]+(prevVi[k,1])^2)+
					  ρALADp/2*σT*(-2xt[k]*prevVt[k,1]+(prevVt[k,1])^2)  for k=1:(horzLen+1)))

			@constraint(tM,tempCon1,xt[1,1]>=evS.τP*xt0+evS.γP*e[1]+evS.ρP*evS.Tamb[stepI,1])
			@constraint(tM,tempCon2[k=1:horzLen],xt[k+1,1]>=evS.τP*xt[k,1]+evS.γP*e[k+1]+evS.ρP*evS.Tamb[stepI+k,1])
			#@constraint(tM,eCon[k=1:horzLen+1],norm([1/2-1/2*e[k] itotal[k]])<=1/2+1/2*e[k])
			@constraint(tM,eCon[k=1:horzLen+1],norm([2itotal[k] e[k]-1])<=e[k]+1)

			objExpCon=0*xt[1:2]
			for k=1:horzLen+1
				tt=[2sqrt(ρALADp/2*σI)*itotal[k];2sqrt(ρALADp/2*σT)*xt[k]]
				append!(objExpCon,tt)
			end
			append!(objExpCon,[t-1])
			@constraint(tM,objCon,norm(objExpCon)<=t+1)
        else
			@objective(tM,Min, sum(-prevLam[k,1]*itotal[k]+
					  ρALADp/2*σI*(itotal[k]-prevVi[k,1])^2+
					  ρALADp/2*σT*(xt[k]-prevVt[k,1])^2  for k=1:(horzLen+1)))
            @NLconstraint(tM,tempCon1,xt[1]-evS.τP*xt0-evS.γP*(itotal[1])^2-evS.ρP*evS.Tamb[stepI,1]==0)
            @NLconstraint(tM,tempCon2[k=1:horzLen],xt[k+1]-evS.τP*xt[k]-evS.γP*(itotal[k+1])^2-evS.ρP*evS.Tamb[stepI+k,1]==0)
        end
        if noTlimit==false
        	@constraint(tM,upperTCon,xt.<=evS.Tmax)
        end
        @constraint(tM,lowerTCon,xt.>=0)
        @constraint(tM,KappaMin,itotal.>=0)
        @constraint(tM,KappaMax,itotal.<=evS.ItotalMax)
        TT = stdout # save original stdout stream
        redirect_stdout()
        statusTM = solve(tM)
        redirect_stdout(TT)
        #@assert statusTM==:Optimal "ALAD XFRM NL optimization not solved to optimality"

        # kappaMax=-getdual(KappaMin)
        # kappaMin=-getdual(KappaMax)
        # tMax=-getdual(upperTCon)
        # tMin=-getdual(lowerTCon)
		if relaxedMode==1
			lambdaTemp=-Gurobi.get_dblattrarray(getrawsolver(tM),"QCPi",1,horzLen+1)
		else
        	lambdaTemp=[-getdual(tempCon1);-getdual(tempCon2)]
		end
        iVal=getvalue(itotal)
        xtVal=getvalue(xt)

        cValMax=abs.(iVal.-evS.ItotalMax).<tolI
        cValMin=abs.(iVal.-0).<tolI
        dLogalad.Ciu[:,p]=1cValMax
		dLogalad.Cil[:,p]=-1cValMin

        cValMax=abs.(xtVal.-evS.Tmax).<tolT
        cValMin=abs.(xtVal.-0).<tolT
        dLogalad.Ctu[:,p]=1cValMax
		dLogalad.Ctl[:,p]=-1cValMin

        dLogalad.Xt[:,p]=xtVal
        dLogalad.Itotal[:,p]=iVal
        dLogalad.Gi[:,p].=0
		dLogalad.Gt[:,p].=0

        #Gz[:,p+1]=σZ*(Vz[:,p]-zVal)-repeat(-Lam[:,p],inner=S)

        for k=1:horzLen+1
            dLogalad.uSum[k,p]=sum(dLogalad.Un[(k-1)*N+n,p] for n=1:N)
            dLogalad.couplConst[k,p]=dLogalad.uSum[k,p] + iD[stepI+(k-1),1] - dLogalad.Itotal[k,p]
        end

        #check for convergence
        constGap=norm(dLogalad.couplConst[:,p],1)
        cc=norm(vcat((prevVu[:,1]-dLogalad.Un[:,p]),(prevVi[:,1]-dLogalad.Itotal[:,p])),1)
        #convCheck=ρALADp*norm(vcat(repeat(σU,horzLen+1,1).*(Vu[:,p]-Un[:,p+1]),σZ*(Vz[:,p]-Z[:,p+1])),1)
        objFun(sn,xt,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,1]     for n=1:N) for k=1:horzLen+1) +
                        sum((xt[k,1]-1)^2*evS.Qsi[N+1,1]                 for k=1:horzLen+1) +
                        sum(sum((u[(k-1)*N+n,1])^2*evS.Ri[n,1]           for n=1:N) for k=1:horzLen+1)
        dLogalad.objVal[1,p]=objFun(dLogalad.Sn[:,p],dLogalad.Xt[:,p],dLogalad.Un[:,p])
        fGap= abs(dLogalad.objVal[1,p]-cSolnl.objVal)
        fGap2= abs((dLogalad.objVal[1,p]-cSolnl.objVal)/cSolnl.objVal)
        snGap=norm((dLogalad.Sn[:,p]-cSolnl.Sn),2)
        unGap=norm((dLogalad.Un[:,p]-cSolnl.Un),2)
        convGap2 = norm(round.(dLogalad.Lam[:,p]-cSolnl.lamCoupl;digits=10)./cSolnl.lamCoupl,2) #need some rounding here if both are 1e-8

        dCMalad.obj[p,1]=fGap
        dCMalad.sn[p,1]=snGap
        dCMalad.un[p,1]=unGap
        dCMalad.couplConst[p,1]=constGap
        convCheck[p,1]=cc
        if  constGap<=epsilon && convCheck<=epsilon
            @printf "Converged after %g iterations\n" p
            convIt=p
            #break
        else
            @printf "convCheck  %e after %g iterations\n" cc p
            @printf "constGap   %e after %g iterations\n" constGap p
			@printf "snGap      %e after %g iterations\n" snGap p
            @printf("fGap       %e after %g iterations\n",fGap,p)
        end

        #coupled QP
		if relaxedMode==2
		 	cM= Model(solver = MosekSolver())
		elseif relaxedMode==1
			cM = Model(solver = GurobiSolver(QCPDual=1))
	    else
	        cM = Model(solver = IpoptSolver())
	    end
        @variable(cM,dUn[1:(N)*(horzLen+1)])
        @variable(cM,dSn[1:(N)*(horzLen+1)])
        @variable(cM,dI[1:(horzLen+1)])
        @variable(cM,dXt[1:(horzLen+1)])
        @variable(cM,relaxS[1:(horzLen+1)])
		objExp=sum(sum(0.5*dUn[(k-1)*N+n,1]^2*Hu[n,1]+dLogalad.Gu[(k-1)*N+n,p]*dUn[(k-1)*N+n,1]+
                       0.5*dSn[(k-1)*N+n,1]^2*Hs[n,1]+dLogalad.Gs[(k-1)*N+n,p]*dSn[(k-1)*N+n,1] for n=1:N)+
                   0.5*dI[k,1]^2*(Hi-lambdaTemp[k,1]*2*evS.γP)+
                   0.5*dXt[k,1]^2*Ht   for k=1:(horzLen+1))
	    objExp=objExp+dot(dLogalad.Gi[:,p],dI)+dot(dLogalad.Gt[:,p],dXt)
        objExp=objExp+prevLam[:,1]'*relaxS+μALADp/2*sum(relaxS[k,1]^2 for k=1:horzLen+1)
    	@objective(cM,Min,objExp)

		Unp=round.(dLogalad.Un[:,p],digits=8)
		Ip=round.(dLogalad.Itotal[:,p],digits=8)
        @constraint(cM,currCon[k=1:horzLen+1],sum(Unp[(k-1)*(N)+n,1]+dUn[(k-1)*(N)+n,1] for n=1:N)-
                                                    (Ip[k,1]+dI[k])==-iD[stepI+(k-1)]+relaxS[k,1])
        @constraint(cM,stateCon1[n=1:N],dSn[n,1]==evS.ηP[n,1]*dUn[n,1])
        @constraint(cM,stateCon2[k=1:horzLen,n=1:N],dSn[n+(k)*(N),1]==dSn[n+(k-1)*(N),1]+evS.ηP[n,1]*dUn[n+(k)*(N),1])
        @constraint(cM,tempCon1,dXt[1,1]==2*evS.γP*dLogalad.Itotal[1,p]*dI[1])
        @constraint(cM,tempCon2[k=1:horzLen],dXt[k+1,1]==evS.τP*dXt[k,1]+2*evS.γP*dLogalad.Itotal[k+1,p]*dI[k+1,1])

		#active local constraints
		if eqForm
			@constraint(cM,dLogalad.Ciu[:,p].*dI.==0)
			@constraint(cM,dLogalad.Cuu[:,p].*dUn.==0)
			@constraint(cM,dLogalad.Csu[:,p].*dSn.==0)
			@constraint(cM,dLogalad.Ctu[:,p].*dXt.==0)
			@constraint(cM,dLogalad.Cil[:,p].*dI.==0)
			@constraint(cM,dLogalad.Cul[:,p].*dUn.==0)
			@constraint(cM,dLogalad.Csl[:,p].*dSn.==0)
			@constraint(cM,dLogalad.Ctl[:,p].*dXt.==0)
		else
			@constraint(cM,dLogalad.Ciu[:,p].*dI.<=0)
	        @constraint(cM,dLogalad.Cuu[:,p].*dUn.<=0)
	        @constraint(cM,dLogalad.Csu[:,p].*dSn.<=0)
	        @constraint(cM,dLogalad.Ctu[:,p].*dXt.<=0)
			@constraint(cM,dLogalad.Cil[:,p].*dI.<=0)
			@constraint(cM,dLogalad.Cul[:,p].*dUn.<=0)
			@constraint(cM,dLogalad.Csl[:,p].*dSn.<=0)
			@constraint(cM,dLogalad.Ctl[:,p].*dXt.<=0)
		end

    	TT = stdout # save original stdout stream
        redirect_stdout()
        statusM = solve(cM)
        redirect_stdout(TT)
        @assert statusM==:Optimal "ALAD Central optimization not solved to optimality"

        #update step
        # Lam[:,p]=-getdual(currCon)
        α1=1
        α2=1
        α3=1
        #α1=α1/ceil(p/2)

		dLogalad.Lam[:,p]=prevLam[:,1]+α3*(-getdual(currCon)-prevLam[:,1])
        #dLogalad.Lam[:,p]=max.(prevLam[:,1]+α3*(-getdual(currCon)-prevLam[:,1]),0)
        dLogalad.Vu[:,p]=prevVu[:,1]+α1*(dLogalad.Un[:,p]-prevVu[:,1])+α2*getvalue(dUn)
        dLogalad.Vi[:,p]=prevVi[:,1]+α1*(dLogalad.Itotal[:,p]-prevVi[:,1])+α2*getvalue(dI)
        dLogalad.Vs[:,p]=prevVs[:,1]+α1*(dLogalad.Sn[:,p]-prevVs[:,1])+α2*getvalue(dSn)
        dLogalad.Vt[:,p]=prevVt[:,1]+α1*(dLogalad.Xt[:,p]-prevVt[:,1])+α2*getvalue(dXt)

		dCMalad.lamIt[p,1]=norm(dLogalad.Lam[:,p]-prevLam[:,1],2)
        dCMalad.lam[p,1]=norm(dLogalad.Lam[:,p]-cSolnl.lamCoupl,2)
        @printf "lastGap    %e after %g iterations\n" dCMalad.lamIt[p,1] p
        @printf "convLamGap %e after %g iterations\n\n" dCMalad.lam[p,1] p

		dLogalad.itUpdate[1,p]=min(ρALADp*ρRate,ρALADmax) #increase ρ every iteration
		μALADp=min(μALADp*μRate,μALADmax) #increase μ every iteration

		ΔY[1,p]=norm(vcat(getvalue(dUn),getvalue(dI),getvalue(dSn),getvalue(dXt)),Inf)

		#reset for next iteration
		prevVu=dLogalad.Vu[:,p]
		prevVs=dLogalad.Vs[:,p]
		prevVi=dLogalad.Vi[:,p]
		prevVt=dLogalad.Vt[:,p]
		prevLam=dLogalad.Lam[:,p]
		ρALADp=dLogalad.itUpdate[1,p]
    end

    return dLogalad,dCMalad,convIt,ΔY,convCheck
end
