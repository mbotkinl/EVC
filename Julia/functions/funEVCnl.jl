#functions to run NL EVC problems


#central
function runNLCentralStep(stepI,evS,cSolnl,cSavenl,relaxedMode,silent)
	K=evS.K
    N=evS.N
    horzLen=min(evS.K1,K-stepI)

	#initialize with current states
	global s0
	global t0

	#desired SOC
	target=zeros(N*(horzLen+1),1)
	for ii=1:N
	   cur=evS.Kn[ii]-(stepI-1)
	   ind=(max(0,(cur-1)*N)+ii):N:length(target)
	   target[ind].=evS.Snmin[ii,1]
	end
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
	@variable(cModel,t[1:(horzLen+1)])
	@variable(cModel,itotal[1:(horzLen+1)])

	#objExp=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,1]+(u[(k-1)*N+n,1])^2*evS.Ri[n,1] for n=1:N) for k=1:(horzLen+1))

	if relaxedMode==1
		objExp =sum((sn[n,1]-1)^2*evS.Qsi[n,1]+(u[n,1])^2*evS.Ri[n,1] for n=1:N)
		for k=2:horzLen+1
			append!(objExp,sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,1]+(u[(k-1)*N+n,1])^2*evS.Ri[n,1]  for n=1:N))
		end
		@constraint(cModel,tempCon1,t[1,1]>=evS.τP*t0+evS.γP*(itotal[1])^2+evS.ρP*evS.Tamb[stepI,1])
		@constraint(cModel,tempCon2[k=1:horzLen],t[k+1,1]>=evS.τP*t[k,1]+evS.γP*(itotal[k+1])^2+evS.ρP*evS.Tamb[stepI+k,1])
	elseif relaxedMode==2
		@variable(cModel,e[1:(horzLen+1)])
		@variable(cModel,t1)
		objExp=t1+sum(-2*evS.Qsi[n,1]*sn[n,1]+evS.Qsi[n,1] for n=1:N)
		for k=2:horzLen+1
			append!(objExp,sum(-2*evS.Qsi[n,1]*sn[(k-1)*(N)+n,1]+evS.Qsi[n,1]  for n=1:N))
		end
		@constraint(cModel,tempCon1,t[1,1]>=evS.τP*t0+evS.γP*e[1]+evS.ρP*evS.Tamb[stepI,1])
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
		append!(objExpCon,[t1-1])

		# objExpCon1 = [sqrt(evS.Qsi[n,1])*sn[(k-1)*(N)+n,1]-evS.Qsi[n,1]/sqrt(evS.Qsi[n,1])  for n=1:N for k=1:(horzLen+1)]
		# objExpCon2 = [sqrt(evS.Ri[n,1])*u[(k-1)*(N)+n,1]  for n=1:N for k=1:(horzLen+1)]
		# objExpCon=vcat(objExpCon1,objExpCon2)

		@constraint(cModel,objCon,norm(objExpCon)<=t1+1)
	else
		objExp =sum((sn[n,1]-1)^2*evS.Qsi[n,1]+(u[n,1])^2*evS.Ri[n,1] for n=1:N)
		for k=2:horzLen+1
			append!(objExp,sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,1]+(u[(k-1)*N+n,1])^2*evS.Ri[n,1]  for n=1:N))
		end
		@NLconstraint(cModel,tempCon1,t[1,1]==evS.τP*t0+evS.γP*(itotal[1])^2+evS.ρP*evS.Tamb[stepI,1])
		@NLconstraint(cModel,tempCon2[k=1:horzLen],t[k+1,1]==evS.τP*t[k,1]+evS.γP*(itotal[k+1])^2+evS.ρP*evS.Tamb[stepI+k,1])
	end
	if slack
		@variable(cModel,slackSn[1:N])
		append!(objExp,sum(evS.β[n]*slackSn[n]^2 for n=1:N))
	end

	@objective(cModel,Min, objExp)

	@constraint(cModel,stateCon1,sn[1:N,1].==s0[1:N,1]+evS.ηP[:,1].*u[1:N,1])
	@constraint(cModel,stateCon2[k=1:horzLen,n=1:N],sn[n+(k)*(N),1]==sn[n+(k-1)*(N),1]+evS.ηP[n,1]*u[n+(k)*(N),1])
	@constraint(cModel,currCon[k=1:horzLen+1],0==-sum(u[(k-1)*(N)+n,1] for n=1:N)-evS.iD_pred[stepI+(k-1)]+itotal[k])
	@constraint(cModel,sn.<=1)
	if slack
		@constraint(cModel,sn.>=target.*(1-repeat(slackSn,horzLen+1,1)))
		@constraint(cModel,slackSn.>=0)
	else
		@constraint(cModel,sn.>=target)
	end
	if noTlimit==false
		@constraint(cModel,upperTCon,t.<=evS.Tmax)
	end
	@constraint(cModel,t.>=0)
	@constraint(cModel,upperCCon,u.<=repeat(evS.imax,horzLen+1,1))
	@constraint(cModel,u.>=repeat(evS.imin,horzLen+1,1))
	@constraint(cModel,itotal.<=evS.ItotalMax)
	@constraint(cModel,itotal.>=0)

	if solverSilent
		@suppress_out begin
			statusC = solve(cModel)
		end
	else
		statusC = solve(cModel)
	end
	@assert statusC==:Optimal "Central NL optimization not solved to optimality"

	uRaw=getvalue(u)
	snRaw=getvalue(sn)
	tRaw=getvalue(t)
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

	#calculate actual temp
	Tactual=zeros(horzLen+1,1)
	Iactual=zeros(horzLen+1,1)
	for k=1:horzLen+1
		Iactual[k,1]=sum((uRaw[(k-1)*N+n,1]) for n=1:N) + evS.iD_actual[stepI+(k-1),1]
	end
	Tactual[1,1]=evS.τP*t0+evS.γP*Iactual[1,1]^2+evS.ρP*evS.Tamb[stepI,1]
	for k=1:horzLen
		Tactual[k+1,1]=evS.τP*Tactual[k,1]+evS.γP*Iactual[k+1,1]^2+evS.ρP*evS.Tamb[stepI+(k-1),1]
	end

	cSolnl.Tpred[stepI,1]=tRaw[1]
	cSolnl.Un[stepI,:]=uRaw[1:N]
	cSolnl.Sn[stepI,:]=snRaw[1:N]
	cSolnl.uSum[stepI,1]=uSum[1]
	cSolnl.Ipred[stepI,1]=itotalRaw[1]
	cSolnl.Tactual[stepI,1]=Tactual[1]
	cSolnl.lamCoupl[stepI,1]=lambdaCurr[1]
	cSolnl.lamTemp[stepI,1]=lambdaTemp[1]

	#cSol.objVal[1,1,stepI]=getobjectivevalue(centralModel)

	if stepI in saveLogInd
		ind=findall(x->x==stepI,saveLogInd)[1]
		cSavenl.Obj[1,1,ind]=getobjectivevalue(cModel)
		uReshape=zeros(horzLen+1,N)
		for ii= 1:N
			uReshape[:,ii]=uRaw[collect(ii:N:length(uRaw))]
		end
		cSavenl.Un[1:(horzLen+1),:,ind]=uReshape
		cSavenl.Lam[1:(horzLen+1),:,ind]=lambdaCurr
		cSavenl.Tactual[1:(horzLen+1),:,ind]=Tactual
	end

	# new states
	t0=round.(cSolnl.Tactual[stepI,1],digits=6)
	s0=round.(cSolnl.Sn[stepI,:],digits=6)
	return nothing
end

function nlEVcentral(evS::scenarioStruct,slack::Bool,relaxedMode::Int,silent::Bool)

    horzLen=evS.K1
    K=evS.K
    N=evS.N
    S=evS.S
    cSolnl=solutionStruct(K=K,N=N,S=S)
    cSavenl=centralLogStruct(logLength=length(saveLogInd),horzLen=horzLen,N=N,S=S)

    Juno.progress() do id
        for stepI=1:K
            @info "$(Dates.format(Dates.now(),"HH:MM:SS")): $(stepI) of $(K)" progress=stepI/K _id=id
            @printf "%s: time step %g of %g....\n" Dates.format(Dates.now(),"HH:MM:SS") stepI K
            try
                runNLCentralStep(stepI,evS,cSolnl,cSavenl,relaxedMode,silent)
            catch e
                @printf "error: %s" e
                break
            end
        end
    end

    objFun(sn,u)=sum(sum((sn[k,n]-1)^2*evS.Qsi[n,1] for n=1:N) +sum((u[k,n])^2*evS.Ri[n,1] for n=1:N) for k=1:evS.K)
    cSolnl.objVal[1,1]=objFun(cSolnl.Sn,cSolnl.Un)

    return cSolnl, cSavenl
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

		#calculate actual temperature from nonlinear model of XFRM
    	itotal=zeros(horzLen+1,1)
    	for k=1:horzLen+1
    		itotal[k,1]=dLog.uSum[k,p] + evS.iD[stepI+(k-1),1]
    	end
    	dLog.Tactual[1,p]=evS.τP*xt0+evS.γP*itotal[1,1]^2+evS.ρP*evS.Tamb[stepI,1]
    	for k=1:horzLen
    		dLog.Tactual[k+1,p]=evS.τP*dLog.Tactual[k,p]+evS.γP*itotal[k+1,1]^2+evS.ρP*evS.Tamb[stepI+k,1]  #fix for mpc
    	end


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
function localEVALAD(evInd::Int,p::Int,stepI::Int,σU::Array{Float64,2},σS::Array{Float64,2},evS::scenarioStruct,dLogaladnl::itLogNL,
    ind,evVu,evVs,itLam,s0,itρ,slack,relaxedMode,solverSilent,roundSigFigs,silent)
	N=evS.N
    horzLen=min(evS.K1,evS.K-stepI)

    tolU=1e-6
    tolS=1e-8

    #evV=zeros(horzLen+1,1)
    target=zeros((horzLen+1),1)
    target[max(1,(evS.Kn[evInd,1]-(stepI-1))):1:length(target),1].=evS.Snmin[evInd,1]
	if relaxedMode==2
		evM = Model(solver = MosekSolver())
	elseif relaxedMode==1
		evM = Model(solver = GurobiSolver(NumericFocus=3))
	else
		#evM = Model(solver = IpoptSolver())
		evM = Model(solver = GurobiSolver(NumericFocus=3))
	end
    @variable(evM,sn[1:(horzLen+1)])
    @variable(evM,u[1:(horzLen+1)])
    if slack @variable(evM,slackSn) end
    objExp=sum((sn[k,1]-1)^2*evS.Qsi[evInd,1]+(u[k,1])^2*evS.Ri[evInd,1]+
                            itLam[k,1]*(u[k,1])+
                            itρ/2*(u[k,1]-evVu[k,1])*σU[evInd,1]*(u[k,1]-evVu[k,1])+
                            itρ/2*(sn[k,1]-evVs[k,1])*σS[evInd,1]*(sn[k,1]-evVs[k,1]) for k=1:horzLen+1)
    if slack
        append!(objExp,sum(evS.β[n]*slackSn^2 for n=1:N))
    end
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
    dLogaladnl.Cuu[ind,p]=1cValMax
    dLogaladnl.Cul[ind,p]=-1cValMin

    cValMax=abs.(snVal.-1).<tolS
    if slack
        cValMin=abs.(snVal.-target*(1-getvalue(slackSn))).<tolS
    else
        cValMin=abs.(snVal.-target).<tolS
    end
    dLogaladnl.Csu[ind,p]=1cValMax
    dLogaladnl.Csl[ind,p]=-1cValMin

    dLogaladnl.slackSn[evInd]= if slack getvalue(slackSn) else 0 end
    dLogaladnl.Sn[ind,p]=round.(snVal,sigdigits=roundSigFigs)
    dLogaladnl.Un[ind,p]=round.(uVal,sigdigits=roundSigFigs)

    dLogaladnl.Gu[ind,p]=round.(2*evS.Ri[evInd,1]*uVal,sigdigits=roundSigFigs)
    dLogaladnl.Gs[ind,p]=round.(2*evS.Qsi[evInd,1]*snVal.-2*evS.Qsi[evInd,1],sigdigits=roundSigFigs)

    return nothing
end

function localXFRMALAD(p::Int,stepI::Int,σI::Float64,σT::Float64,evS::scenarioStruct,dLogaladnl::itLogNL,
    itLam,itVt,itVi,itρ,roundSigFigs)
	#N+1 decoupled problem aka transformer current
	tolT=1e-6
    tolI=1e-6
	horzLen=min(evS.K1,evS.K-stepI)

	if relaxedMode==2
		tM = Model(solver = MosekSolver())
	elseif relaxedMode==1
		tM = Model(solver = GurobiSolver(QCPDual=1))
	else
		tM = Model(solver = IpoptSolver())
	end
	@variable(tM,itotal[1:(horzLen+1)])
	@variable(tM,t[1:(horzLen+1)])
	if relaxedMode ==1
		@objective(tM,Min, sum(-itLam[k,1]*itotal[k]+
				  itρ/2*σI*(itotal[k]-itVi[k,1])^2+
				  itρ/2*σT*(t[k]-itVt[k,1])^2  for k=1:(horzLen+1)))
		@constraint(tM,tempCon1,t[1]-evS.τP*t0-evS.γP*(itotal[1])^2-evS.ρP*evS.Tamb[stepI,1]>=0)
		@constraint(tM,tempCon2[k=1:horzLen],t[k+1]-evS.τP*t[k]-evS.γP*(itotal[k+1])^2-evS.ρP*evS.Tamb[stepI+k,1]>=0)
	elseif relaxedMode==2
		@variable(tM,e[1:(horzLen+1)])
		@variable(tM,t1)

		@objective(tM,Min, sum(-itLam[k,1]*itotal[k]+
				  itρ/2*σI*(-2itotal[k]*itVi[k,1]+(itVi[k,1])^2)+
				  itρ/2*σT*(-2xt[k]*itVt[k,1]+(itVt[k,1])^2)  for k=1:(horzLen+1)))

		@constraint(tM,tempCon1,t[1,1]>=evS.τP*t0+evS.γP*e[1]+evS.ρP*evS.Tamb[stepI,1])
		@constraint(tM,tempCon2[k=1:horzLen],t[k+1,1]>=evS.τP*t[k,1]+evS.γP*e[k+1]+evS.ρP*evS.Tamb[stepI+k,1])
		#@constraint(tM,eCon[k=1:horzLen+1],norm([1/2-1/2*e[k] itotal[k]])<=1/2+1/2*e[k])
		@constraint(tM,eCon[k=1:horzLen+1],norm([2itotal[k] e[k]-1])<=e[k]+1)

		objExpCon=0*xt[1:2]
		for k=1:horzLen+1
			tt=[2sqrt(ρALADp/2*σI)*itotal[k];2sqrt(ρALADp/2*σT)*t[k]]
			append!(objExpCon,tt)
		end
		append!(objExpCon,[t1-1])
		@constraint(tM,objCon,norm(objExpCon)<=t1+1)
	else
		@objective(tM,Min, sum(-itLam[k,1]*itotal[k]+
				  itρ/2*σI*(itotal[k]-itVi[k,1])^2+
				  itρ/2*σT*(t[k]-itVt[k,1])^2  for k=1:(horzLen+1)))
		@NLconstraint(tM,tempCon1,t[1]-evS.τP*t0-evS.γP*(itotal[1])^2-evS.ρP*evS.Tamb[stepI,1]==0)
		@NLconstraint(tM,tempCon2[k=1:horzLen],t[k+1]-evS.τP*t[k]-evS.γP*(itotal[k+1])^2-evS.ρP*evS.Tamb[stepI+k,1]==0)
	end
	if noTlimit==false
		@constraint(tM,upperTCon,t.<=evS.Tmax)
	end
	@constraint(tM,lowerTCon,t.>=0)
	@constraint(tM,KappaMin,itotal.>=0)
	@constraint(tM,KappaMax,itotal.<=evS.ItotalMax)

	if solverSilent
        @suppress_out begin
            statusTM = solve(tM)
        end
    else
        statusTM = solve(tM)
    end
	@assert statusTM==:Optimal "ALAD XFRM NLP optimization not solved to optimality"


	if relaxedMode==1
		lambdaTemp=-Gurobi.get_dblattrarray(getrawsolver(tM),"QCPi",1,horzLen+1)
	else
		lambdaTemp=[-getdual(tempCon1);-getdual(tempCon2)]
	end
	lambdaTemp=round.(min.(lambdaTemp,0),sigdigits=roundSigFigs)

	iVal=getvalue(itotal)
	tVal=getvalue(t)

	cValMax=abs.(iVal.-evS.ItotalMax).<tolI
	cValMin=abs.(iVal.-0).<tolI
	dLogaladnl.Ciu[:,p]=1cValMax
	dLogaladnl.Cil[:,p]=-1cValMin

	cValMax=abs.(tVal.-evS.Tmax).<tolT
	cValMin=abs.(tVal.-0).<tolT
	dLogaladnl.Ctu[:,p]=1cValMax
	dLogaladnl.Ctl[:,p]=-1cValMin

	dLogaladnl.Tpred[:,p]=round.(tVal,sigdigits=roundSigFigs)
	dLogaladnl.Ipred[:,p]=round.(iVal,sigdigits=roundSigFigs)
	dLogaladnl.Gi[:,p].=0
	dLogaladnl.Gt[:,p].=0
	return lambdaTemp
end

function coordALAD(p::Int,stepI::Int,μALADp::Float64,evS::scenarioStruct,itLam,itVu,itVs,itVt,itVi,itρ,
    lambdaTemp,dLogaladnl::itLogNL,roundSigFigs)

	horzLen=min(evS.K1,evS.K-stepI)

	Hu=2*evS.Ri#*(1+rand())
	Hs=2*evS.Qsi#*(1+rand())
	Hi=0
	Ht=0
	ρALAD=1e3
	ρRate=1.1
	ρALADmax=1e6

	#coupled QP
	if relaxedMode==2
		cM= Model(solver = MosekSolver())
	elseif relaxedMode==1
		cM = Model(solver = GurobiSolver(QCPDual=1))
	else
		#cM = Model(solver = IpoptSolver())
		cM = Model(solver = GurobiSolver(NumericFocus=3))
	end
	@variable(cM,dUn[1:(N)*(horzLen+1)])
	@variable(cM,dSn[1:(N)*(horzLen+1)])
	@variable(cM,dI[1:(horzLen+1)])
	@variable(cM,dT[1:(horzLen+1)])
	@variable(cM,relaxS[1:(horzLen+1)])
	objExp=sum(sum(0.5*dUn[(k-1)*N+n,1]^2*Hu[n,1]+dLogaladnl.Gu[(k-1)*N+n,p]*dUn[(k-1)*N+n,1]+
				   0.5*dSn[(k-1)*N+n,1]^2*Hs[n,1]+dLogaladnl.Gs[(k-1)*N+n,p]*dSn[(k-1)*N+n,1] for n=1:N)+
			   0.5*dI[k,1]^2*(Hi-lambdaTemp[k,1]*2*evS.γP)+
			   0.5*dT[k,1]^2*Ht   for k=1:(horzLen+1))
	objExp=objExp+dot(dLogaladnl.Gi[:,p],dI)+dot(dLogaladnl.Gt[:,p],dT)
	objExp=objExp+itLam[:,1]'*relaxS+μALADp/2*sum(relaxS[k,1]^2 for k=1:horzLen+1)
	@objective(cM,Min,objExp)

	Unp=round.(dLogaladnl.Un[:,p],digits=8)
	Ip=round.(dLogaladnl.Ipred[:,p],digits=8)
	@constraint(cM,currCon[k=1:horzLen+1],sum(Unp[(k-1)*(N)+n,1]+dUn[(k-1)*(N)+n,1] for n=1:N)-
												(Ip[k,1]+dI[k])==-evS.iD_pred[stepI+(k-1)]+relaxS[k,1])
	@constraint(cM,stateCon1[n=1:N],dSn[n,1]==evS.ηP[n,1]*dUn[n,1])
	@constraint(cM,stateCon2[k=1:horzLen,n=1:N],dSn[n+(k)*(N),1]==dSn[n+(k-1)*(N),1]+evS.ηP[n,1]*dUn[n+(k)*(N),1])
	@constraint(cM,tempCon1,dT[1,1]==2*evS.γP*Ip[1]*dI[1])
	@constraint(cM,tempCon2[k=1:horzLen],dT[k+1,1]==evS.τP*dT[k,1]+2*evS.γP*Ip[k+1,1]*dI[k+1,1])

	#active local constraints
	if eqForm
		@constraint(cM,dLogaladnl.Ciu[:,p].*dI.==0)
		@constraint(cM,dLogaladnl.Cuu[:,p].*dUn.==0)
		@constraint(cM,dLogaladnl.Csu[:,p].*dSn.==0)
		@constraint(cM,dLogaladnl.Ctu[:,p].*dT.==0)
		@constraint(cM,dLogaladnl.Cil[:,p].*dI.==0)
		@constraint(cM,dLogaladnl.Cul[:,p].*dUn.==0)
		@constraint(cM,dLogaladnl.Csl[:,p].*dSn.==0)
		@constraint(cM,dLogaladnl.Ctl[:,p].*dT.==0)
	else
		@constraint(cM,dLogaladnl.Ciu[:,p].*dI.<=0)
		@constraint(cM,dLogaladnl.Cuu[:,p].*dUn.<=0)
		@constraint(cM,dLogaladnl.Csu[:,p].*dSn.<=0)
		@constraint(cM,dLogaladnl.Ctu[:,p].*dT.<=0)
		@constraint(cM,dLogaladnl.Cil[:,p].*dI.<=0)
		@constraint(cM,dLogaladnl.Cul[:,p].*dUn.<=0)
		@constraint(cM,dLogaladnl.Csl[:,p].*dSn.<=0)
		@constraint(cM,dLogaladnl.Ctl[:,p].*dT.<=0)
	end

	if solverSilent
        @suppress_out begin
            statusM = solve(cM)
        end
    else
        statusM = solve(cM)
    end
	@assert statusM==:Optimal "ALAD Central optimization not solved to optimality"

	#update step
	# Lam[:,p]=-getdual(currCon)
	α1=1
	α2=1
	α3=1
	#α1=α1/ceil(p/2)

	#dLogaladnl.Lam[:,p]=round.(itLam[:,1]+α3*(-getdual(currCon)-itLam[:,1]),sigdigits=roundSigFigs)
	dLogaladnl.Lam[:,p]=round.(max.(itLam[:,1]+α3*(-getdual(currCon)-itLam[:,1]),0),sigdigits=roundSigFigs)
	dLogaladnl.Vu[:,p]=round.(itVu[:,1]+α1*(dLogaladnl.Un[:,p]-itVu[:,1])+α2*getvalue(dUn),sigdigits=roundSigFigs)
	dLogaladnl.Vi[:,p]=round.(itVi[:,1]+α1*(dLogaladnl.Ipred[:,p]-itVi[:,1])+α2*getvalue(dI),sigdigits=roundSigFigs)
	dLogaladnl.Vs[:,p]=round.(itVs[:,1]+α1*(dLogaladnl.Sn[:,p]-itVs[:,1])+α2*getvalue(dSn),sigdigits=roundSigFigs)
	dLogaladnl.Vt[:,p]=round.(itVt[:,1]+α1*(dLogaladnl.Tpred[:,p]-itVt[:,1])+α2*getvalue(dT),sigdigits=roundSigFigs)

	# dCMalad.lamIt[p,1]=norm(dLogaladnl.Lam[:,p]-itLam[:,1],2)
	# dCMalad.lam[p,1]=norm(dLogaladnl.Lam[:,p]-cSolnl.lamCoupl,2)
    if !silent
		#@printf "lastGap    %e after %g iterations\n" dCMnl.lamIt[p,1] p
		#@printf "convLamGap %e after %g iterations\n\n" dCMnl.lam[p,1] p
	end
	dLogaladnl.itUpdate[1,p]=min(itρ*ρRate,ρALADmax) #increase ρ every iteration
	#μALADp=min(μALADp*μRate,μALADmax) #increase μ every iteration
	#ΔY[1,p]=norm(vcat(getvalue(dUn),getvalue(dI),getvalue(dSn),getvalue(dXt)),Inf)
	return nothing
end

function runNLALADIt(p,stepI,evS,itLam,itVu,itVs,itVt,itVi,itρ,dLogaladnl,dCMnl,dSolnl,cSave,eqForm,roundSigFigs,silent)
	K=evS.K
    N=evS.N
    horzLen=min(evS.K1,K-stepI)

    global s0
    global t0

    #ALADIN tuning and initial guess
    ##σs are tuned to i_n [.010kA]
	# if eqForm
	# 	#println("Running Eq ALADIN")
	# 	scalingF=1
	# else
	# 	#println("Running ineq ALADIN")
	# 	scalingF=1e-4
	# end
    # σU=1*ones(N,1)
    # σS=ones(N,1)/1e1
    # σI=1.0*scalingF
    # σT=1.5*scalingF
	μALADp=1e8

	σI=0.5
    σT=2.0
    σU=ones(N,1)/.02
    σS=ones(N,1)

    # μALAD=1e8
    # μRate=1
    # μALADmax=2e9


	#solve decouple
	if runParallel
		@sync @distributed for evInd=1:N
			ind=[evInd]
			for k=1:horzLen
				append!(ind,k*N+evInd)
			end
			evVu=itVu[ind,1]
			evVs=itVs[ind,1]
			localEVALAD(evInd,p,stepI,σU,σS,evS,dLogaladnl,ind,evVu,evVs,itLam,s0,itρ,
						slack,relaxedMode,solverSilent,roundSigFigs,silent)
		end
	else
		for evInd=1:N
			#print(evInd)
			ind=[evInd]
			for k=1:horzLen
				append!(ind,k*N+evInd)
			end
			evVu=itVu[ind,1]
			evVs=itVs[ind,1]
			localEVALAD(evInd,p,stepI,σU,σS,evS,dLogaladnl,ind,evVu,evVs,itLam,s0,itρ,
						slack,relaxedMode,solverSilent,roundSigFigs,silent)
		end
	end

	lambdaTemp=localXFRMALAD(p,stepI,σI,σT,evS,dLogaladnl,itLam,itVt,itVi,itρ,roundSigFigs)

	for k=1:horzLen+1
		dLogaladnl.uSum[k,p]=sum(dLogaladnl.Un[(k-1)*N+n,p] for n=1:N)
		dLogaladnl.couplConst[k,p]=dLogaladnl.uSum[k,p] + evS.iD_pred[stepI+(k-1),1] - dLogaladnl.Ipred[k,p]
	end

	#calculate actual temperature from nonlinear model of XFRM
	for k=1:horzLen+1
		dLogaladnl.Iactual[k,p]=dLogaladnl.uSum[k,p] + evS.iD_actual[stepI+(k-1),1]
	end
	dLogaladnl.Tactual[1,p]=evS.τP*t0+evS.γP*dLogaladnl.Iactual[1,p]^2+evS.ρP*evS.Tamb[stepI,1] #fix for mpc
	for k=1:horzLen
		dLogaladnl.Tactual[k+1,p]=evS.τP*dLogaladnl.Tactual[k,p]+evS.γP*dLogaladnl.Iactual[k,p]^2+evS.ρP*evS.Tamb[stepI+k,1]  #fix for mpc
	end

	coordALAD(p,stepI,μALADp,evS,itLam,itVu,itVs,itVt,itVi,itρ,lambdaTemp,dLogaladnl,roundSigFigs)

	#check for convergence
	constGap=norm(dLogaladnl.couplConst[:,p],1)
	#cc=norm(vcat((itVu[:,1]-dLogaladnl.Un[:,p]),(itVi[:,1]-dLogaladnl.Ipred[:,p])),1)
	#convCheck=ρALADp*norm(vcat(repeat(σU,horzLen+1,1).*(Vu[:,p]-Un[:,p+1]),σZ*(Vz[:,p]-Z[:,p+1])),1)
	objFun(sn,xt,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,1]     for n=1:N) for k=1:horzLen+1) +
					sum((xt[k,1]-1)^2*evS.Qsi[N+1,1]                 for k=1:horzLen+1) +
					sum(sum((u[(k-1)*N+n,1])^2*evS.Ri[n,1]           for n=1:N) for k=1:horzLen+1)
	dLogaladnl.objVal[1,p]=objFun(dLogaladnl.Sn[:,p],dLogaladnl.Tpred[:,p],dLogaladnl.Un[:,p])
	itGap = norm(dLogaladnl.Lam[:,p]-itLam[:,1],2)

	#only if saving
    if stepI in saveLogInd
        ind=findall(x->x==stepI,saveLogInd)[1]
        dCMnl.coupl1Norm[p,ind]=constGap
        dCMnl.lamIt2Norm[p,ind]=itGap
        dCMnl.objAbs[p,ind]=abs(dLogaladnl.objVal[1,p]-cSavenl.Obj[1,1,ind])
        dCMnl.objPerc[p,ind]=abs(dLogaladnl.objVal[1,p]-cSavenl.Obj[1,1,ind])/cSavenl.Obj[1,1,ind]*100
        dCMnl.lam1Norm[p,ind]= norm(dLogaladnl.Lam[:,p]-cSavenl.Lam[1:(horzLen+1),:,ind],1)
        dCMnl.lam2Norm[p,ind]= norm(dLogaladnl.Lam[:,p]-cSavenl.Lam[1:(horzLen+1),:,ind],2)
        dCMnl.lamInfNorm[p,ind]= norm(dLogaladnl.Lam[:,p]-cSavenl.Lam[1:(horzLen+1),:,ind],Inf)
        dCMnl.t1Norm[p,ind]= norm(dLogaladnl.Tactual[:,p]-cSavenl.Tactual[1:(horzLen+1),:,ind],1)
        dCMnl.t2Norm[p,ind]= norm(dLogaladnl.Tactual[:,p]-cSavenl.Tactual[1:(horzLen+1),:,ind],2)
        dCMnl.tInfNorm[p,ind]= norm(dLogaladnl.Tactual[:,p]-cSavenl.Tactual[1:(horzLen+1),:,ind],Inf)
        uReshape=zeros(horzLen+1,N)
        for ii= 1:N
            uReshape[:,ii]=dLogaladnl.Un[collect(ii:N:length(dLogaladnl.Un[:,p])),p]
        end
        dCMnl.un1Norm[p,ind]= norm(uReshape-cSavenl.Un[1:(horzLen+1),:,ind],1)
        dCMnl.un2Norm[p,ind]= norm(uReshape-cSavenl.Un[1:(horzLen+1),:,ind],2)
        dCMnl.unInfNorm[p,ind]= norm(uReshape-cSavenl.Un[1:(horzLen+1),:,ind],Inf)
    end

	if  constGap<=primChk && itGap<=dualChk
		@printf "Converged after %g iterations\n" p
		convIt=p
		return true
	else
		if !silent
			#@printf "convCheck  %e after %g iterations\n" cc p
			@printf "lamIt      %e after %g iterations\n" itGap p
			@printf "constGap   %e after %g iterations\n\n" constGap p
		end
	end

	return false
end

function runNLALADStep(stepI,maxIt,evS,dSolnl,dCMnl,cSave,eqForm,roundSigFigs,silent)
    K=evS.K
    N=evS.N
    horzLen=min(evS.K1,K-stepI)
    #ogρ=ρALADp #save to reset later
    #convCheck=zeros(maxIt+1,1)
    #ΔY=zeros(1,maxIt+1)
    dLogaladnl=itLogNL(horzLen=horzLen,N=N)
    p=1
    timeStart=now()
    while (p<=maxIt && round(now()-timeStart,Second)<=Dates.Second(9/10*evS.Ts))
    #while (p<=maxIt)
        #global p
        #@printf "%git" p
        if p==1
            itLam=prevLam
            itVu=prevVu
            itVs=prevVs
			itVi=prevVi
            itVt=prevVt
            itρ=ρALADp
        else
            itLam=round.(dLogaladnl.Lam[:,(p-1)],sigdigits=roundSigFigs)
            itVu=round.(dLogaladnl.Vu[:,(p-1)],sigdigits=roundSigFigs)
            itVs=round.(dLogaladnl.Vs[:,(p-1)],sigdigits=roundSigFigs)
			itVi=round.(dLogaladnl.Vi[:,(p-1)],sigdigits=roundSigFigs)
            itVt=round.(dLogaladnl.Vt[:,(p-1)],sigdigits=roundSigFigs)
            itρ=round.(dLogaladnl.itUpdate[1,(p-1)],sigdigits=roundSigFigs)
        end
        cFlag=runNLALADIt(p,stepI,evS,itLam,itVu,itVs,itVt,itVi,itρ,dLogaladnl,dCMnl,dSolnl,cSavenl,eqForm,roundSigFigs,silent)
        global convIt=p
        if cFlag
            break
        end
        p+=1
    end
    if stepI in saveLogInd
        ind=findall(x->x==stepI,saveLogInd)[1]
        dCMnl.convIt[1,1,ind]=convIt
    end
	# print(round(now()-timeStart,Second))
    # #
	#
	# indCheck=2
	# plot(hcat(dLogaladnl.uSum[:,indCheck],dLogalad.uSum[:,indCheck]),label=["NL" "PWL"])
	# plot(hcat(dLogaladnl.Tactual[:,indCheck],dLogalad.Tactual[:,indCheck]),label=["NL" "PWL"])
	# plot(hcat(dLogaladnl.Lam[:,indCheck],dLogalad.Lam[:,indCheck]),label=["NL" "PWL"])
	# plot(hcat(dLogaladnl.Ipred[:,indCheck],dLogalad.Itotal[:,indCheck]),label=["NL" "PWL"])
	# plot(hcat(dLogaladnl.couplConst[:,indCheck],dLogalad.couplConst[:,indCheck]),label=["NL" "PWL"])
	#
	# plot(hcat(dLogaladnl.Gu[:,indCheck],dLogalad.Gu[:,indCheck]),label=["NL" "PWL"])
	# plot(hcat(dLogaladnl.Gs[:,indCheck],dLogalad.Gs[:,indCheck]),label=["NL" "PWL"])
	# plot(hcat(dLogaladnl.Ctu[:,indCheck],dLogalad.Ctu[:,indCheck]),label=["NL" "PWL"])
	#
	#
	# plot(hcat(dLogaladnl.Vt[:,indCheck],dLogalad.Vt[:,indCheck]),label=["NL" "PWL"])
	# plot(hcat(dLogaladnl.Vu[:,indCheck],dLogalad.Vu[:,indCheck]),label=["NL" "PWL"])
	# plot(hcat(dLogaladnl.Vs[:,indCheck],dLogalad.Vs[:,indCheck]),label=["NL" "PWL"])
	#
	#
	#
    # xPlotnl=zeros(horzLen+1,N)
    # uPlotnl=zeros(horzLen+1,N)
    # for ii= 1:N
    # 	xPlotnl[:,ii]=dLogaladnl.Sn[collect(ii:N:length(dLogaladnl.Sn[:,convIt])),convIt]
    #     uPlotnl[:,ii]=dLogaladnl.Un[collect(ii:N:length(dLogaladnl.Un[:,convIt])),convIt]
    # end
	# p1=plot(dLogaladnl.uSum[:,convIt],xlabel="Time",ylabel="Current Sum (kA)",xlims=(0,horzLen+1),label="ALADIN Open Loop")
	# plot!(p1,sum(cSavenl.Un[:,:,ind],dims=2),seriescolor=:black,linewidth=2,linealpha=0.8,label="Central Open Loop")
	#
    # plot(xPlotnl)
    # pd3alad=plot(hcat(dLogaladnl.Tactual[:,convIt],dLogaladnl.Tpred[:,convIt])*1000,label=["Actual Temp" "pred Temp"],xlims=(0,evS.K),xlabel="Time",ylabel="Temp (K)")
    # plot!(pd3alad,1:horzLen+1,evS.Tmax*ones(evS.K)*1000,label="XFRM Limit",line=(:dash,:red))
	#
    # #convergence plots
    # halfCI=Int(floor(convIt/2))
    # if halfCI>0
    #     CList=reshape([range(colorant"blue", stop=colorant"yellow",length=halfCI);
    #                    range(colorant"yellow", stop=colorant"red",length=convIt-halfCI)], 1, convIt);
    # else
    #     CList=colorant"red";
    # end
    # uSumPlotalad=plot(dLogaladnl.uSum[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Current Sum (kA)",xlims=(0,horzLen+1),legend=false)
    # plot!(uSumPlotalad,sum(cSavenl.Un[:,:,ind],dims=2),seriescolor=:black,linewidth=2,linealpha=0.8)
	#
    # constPlotalad2=plot(dLogaladnl.couplConst[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="curr constraint diff",xlims=(0,horzLen+1),legend=false)
	#
    # lamPlotalad=plot(dLogaladnl.Lam[:,1:convIt],seriescolor=CList,xlabel="Time",ylabel="Lambda",legend=false)
    # plot!(lamPlotalad,cSavenl.Lam[:,:,ind],seriescolor=:black,linewidth=2,linealpha=0.8)
	#
    # activeSet=zeros(convIt,1)
    # setChanges=zeros(convIt,1)
    # for ii=2:convIt
    #     activeSet[ii,1]=sum(abs.(dLogaladnl.Csu[:,ii]))+sum(abs.(dLogaladnl.Ctu[:,ii]))+
    #               		sum(abs.(dLogaladnl.Cuu[:,ii]))+sum(abs.(dLogaladnl.Ciu[:,ii]))+
    # 					sum(abs.(dLogaladnl.Csl[:,ii]))+sum(abs.(dLogaladnl.Ctl[:,ii]))+
    # 				    sum(abs.(dLogaladnl.Cul[:,ii]))+sum(abs.(dLogaladnl.Cil[:,ii]))
    #     setChanges[ii,1]=sum(abs.(dLogaladnl.Csu[:,ii]-dLogaladnl.Csu[:,ii-1]))+sum(abs.(dLogaladnl.Ctu[:,ii]-dLogaladnl.Ctu[:,ii-1]))+
    #                      sum(abs.(dLogaladnl.Cuu[:,ii]-dLogaladnl.Cuu[:,ii-1]))+sum(abs.(dLogaladnl.Ciu[:,ii]-dLogaladnl.Ciu[:,ii-1]))+
    # 					 sum(abs.(dLogaladnl.Csl[:,ii]-dLogaladnl.Csl[:,ii-1]))+sum(abs.(dLogaladnl.Ctl[:,ii]-dLogaladnl.Ctl[:,ii-1]))+
    # 				     sum(abs.(dLogaladnl.Cul[:,ii]-dLogaladnl.Cul[:,ii-1]))+sum(abs.(dLogaladnl.Cil[:,ii]-dLogaladnl.Cil[:,ii-1]))
    # end
	#
    # activeSetPlot=plot(2:convIt,activeSet[2:convIt],xlabel="Iteration",ylabel="Total Active inequality constraints",
    #                    legend=false,xlims=(2,convIt))
    # setChangesPlot=plot(10:convIt,setChanges[10:convIt],xlabel="Iteration",ylabel="Changes in Active inequality constraints",
    #                   legend=false,xlims=(2,convIt))
    # #solChangesplot=plot(2:convIt,hcat(ΔY[2:convIt],convCheck[2:convIt]),xlabel="Iteration",labels=["ΔY" "y-x"],xlims=(1,convIt))
	#
    # fPlotalad=plot(dCMnl.objAbs[1:convIt,1],xlabel="Iteration",ylabel="obj function gap",xlims=(1,convIt),legend=false,yscale=:log10)
    # convItPlotalad=plot(dCMnl.lamIt2Norm[1:convIt,1],xlabel="Iteration",ylabel="2-Norm It Lambda Gap",xlims=(1,convIt),legend=false,yscale=:log10)
    # convPlotalad=plot(dCMnl.lam2Norm[1:convIt,1],xlabel="Iteration",ylabel="central lambda gap",xlims=(1,convIt),legend=false,yscale=:log10)
    # constPlotalad=plot(dCMnl.coupl1Norm[1:convIt,1],xlabel="Iteration",ylabel="curr constraint Gap",xlims=(1,convIt),legend=false,yscale=:log10)

    #save current state and update for next timeSteps
    dSolnl.Tpred[stepI,1]=dLogaladnl.Tpred[1,convIt]
    dSolnl.Un[stepI,:]=dLogaladnl.Un[1:N,convIt]
    dSolnl.Sn[stepI,:]=dLogaladnl.Sn[1:N,convIt]
    dSolnl.uSum[stepI,1]=dLogaladnl.uSum[1,convIt]
    #dSolnl.Ipred[stepI,1]=dLogaladnl.Ipred[1,convIt]
	dSolnl.Itotal[stepI,1]=dLogaladnl.Iactual[1,convIt]
    dSolnl.Tactual[stepI,1]=dLogaladnl.Tactual[1,convIt]
    dSolnl.convIt[stepI,1]=convIt

    # new states
    global t0=round(dSolnl.Tactual[stepI,1],sigdigits=roundSigFigs)
    global s0=dSolnl.Sn[stepI,:]

    #function getAttr()
    #clean this up
    if convIt==1
        dSolnl.lamCoupl[stepI,1]=prevLam[1,1]
        if stepI+horzLen==evS.K
            newLam=prevLam[2:horzLen+1,1]
            newVu=prevVu[(N+1):(N*(horzLen+1)),1]
            newVt=prevVt[2:horzLen+1,1]
            newVs=prevVs[(N+1):(N*(horzLen+1)),1]
        else
            newLam=vcat(prevLam[2:horzLen+1,1],prevLam[horzLen+1,1])
            newVu=vcat(prevVu[(N+1):(N*(horzLen+1)),1],prevVu[((N*horzLen)+1):(N*(horzLen+1)),1])
            newVt=vcat(prevVt[2:horzLen+1,1],prevVt[horzLen+1,1])
            newVs=vcat(prevVs[(N+1):(N*(horzLen+1)),1],prevVs[((N*horzLen)+1):(N*(horzLen+1)),1])
        end
    else
        dSolnl.lamCoupl[stepI,1]=dLogaladnl.Lam[1,convIt-1]
        if stepI+horzLen==evS.K
            newLam=dLogaladnl.Lam[2:horzLen+1,convIt-1]
            newVu=dLogaladnl.Vu[(N+1):(N*(horzLen+1)),convIt-1]
            newVt=dLogaladnl.Vt[2:horzLen+1,convIt-1]
            newVs=dLogaladnl.Vs[(N+1):(N*(horzLen+1)),convIt-1]
        else
            newLam=vcat(dLogaladnl.Lam[2:horzLen+1,convIt-1],dLogaladnl.Lam[horzLen+1,convIt-1])
            newVu=vcat(dLogaladnl.Vu[(N+1):(N*(horzLen+1)),convIt-1],dLogaladnl.Vu[((N*horzLen)+1):(N*(horzLen+1)),convIt-1])
            newVt=vcat(dLogaladnl.Vt[2:horzLen+1,convIt-1],dLogaladnl.Vt[horzLen+1,convIt-1])
            newVs=vcat(dLogaladnl.Vs[(N+1):(N*(horzLen+1)),convIt-1],dLogaladnl.Vs[((N*horzLen)+1):(N*(horzLen+1)),convIt-1])
        end
    end

    global prevLam=round.(newLam,sigdigits=roundSigFigs)
    global prevVu=round.(newVu,sigdigits=roundSigFigs)
    global prevVt=round.(newVt,sigdigits=roundSigFigs)
    global prevVs=round.(newVs,sigdigits=roundSigFigs)
    #global ρALADp=round.(ogρ,digits=2)

    return nothing
end

function nlEValad(maxIt::Int,evS::scenarioStruct,cSave::centralLogStruct,slack::Bool,eqForm::Bool,roundSigFigs::Int,silent::Bool)

    horzLen=evS.K1
    K=evS.K
    N=evS.N
    dSolnl=solutionStruct(K=K,N=N,S=evS.S)
    dCMnl=convMetricsStruct(maxIt=maxIt,logLength=length(saveLogInd))

    Juno.progress() do id
        for stepI=1:K
            @info "$(Dates.format(Dates.now(),"HH:MM:SS")): $(stepI) of $(K)" progress=stepI/K _id=id
            @printf "%s: time step %g of %g...." Dates.format(Dates.now(),"HH:MM:SS") stepI K
            try
                runNLALADStep(stepI,maxIt,evS,dSolnl,dCMnl,cSave,eqForm,roundSigFigs,silent)
                @printf "convIt: %g\n" dSolnl.convIt[stepI,1]
            catch e
                @printf "error: %s" e
                break
            end
        end
    end

    objFun(sn,u)=sum(sum((sn[k,n]-1)^2*evS.Qsi[n,1] for n=1:N) +sum((u[k,n])^2*evS.Ri[n,1] for n=1:N) for k=1:evS.K)
    dSolnl.objVal[1,1]=objFun(dSolnl.Sn,dSolnl.Un)

    return dSolnl, dCMnl
end
