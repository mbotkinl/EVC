#functions to run PWL EVC problems


#central
function pwlEVcentral(N::Int,S::Int,horzLen::Int,evS::scenarioStruct)
    #initialize with current states
    sn0=evS.s0
    xt0=evS.t0
    stepI=1
    #desired SOC
    target=zeros(N*(horzLen+1),1);
    for ii=1:N
       cur=evS.Kn[ii]-(stepI-1)
       ind=max(0,(cur-1)*N)+ii:N:length(target)
       target[ind]=evS.Snmin[ii,1]
    end

    println("setting up model")
    centralModel = Model(solver = GurobiSolver(Presolve=0,BarHomogeneous=1,NumericFocus=3))

    #u w and z are one index ahead of x. i.e the x[k+1]=x[k]+eta*u[k+1]
    @variable(centralModel,u[1:N*(horzLen+1)])
    @variable(centralModel,sn[1:(N)*(horzLen+1)])
    @variable(centralModel,xt[1:(horzLen+1)])
    @variable(centralModel,z[1:evS.S*(horzLen+1)])
    println("obj")
    @objective(centralModel,Min,sum(sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,1]+(u[(k-1)*N+n,1])^2*evS.Ri[n,1]  for n=1:N) for k=1:horzLen+1))
    println("constraints")
    @constraint(centralModel,stateCon1,sn[1:N,1].==sn0[1:N,1]+evS.ηP[:,1].*u[1:N,1])
    @constraint(centralModel,stateCon2[k=1:horzLen,n=1:N],sn[n+(k)*(N),1]==sn[n+(k-1)*(N),1]+evS.ηP[n,1]*u[n+(k)*(N),1])
    @constraint(centralModel,tempCon1,xt[1,1]==evS.τP*xt0+evS.γP*evS.deltaI*sum((2*s-1)*z[s,1] for s=1:S)+evS.ρP*evS.w[stepI*2,1])
    @constraint(centralModel,tempCon2[k=1:horzLen],xt[k+1,1]==evS.τP*xt[k,1]+evS.γP*evS.deltaI*sum((2*s-1)*z[k*S+s,1] for s=1:S)+evS.ρP*evS.w[stepI*2+k*2,1])
    @constraint(centralModel,currCon[k=1:horzLen+1],0==-sum(u[(k-1)*(N)+n] for n=1:N)-evS.w[(k-1)*2+1]+sum(z[(k-1)*(S)+s] for s=1:S))
    @constraint(centralModel,sn.<=1)
    @constraint(centralModel,sn.>=target)
    if noTlimit==0
    	@constraint(centralModel,upperTCon,xt.<=evS.Tmax)
    end
    @constraint(centralModel,xt.>=0)
    @constraint(centralModel,upperCCon,u.<=repmat(evS.imax,horzLen+1,1))
    @constraint(centralModel,u.>=repmat(evS.imin,horzLen+1,1))
    @constraint(centralModel,z.>=0)
    @constraint(centralModel,z.<=evS.deltaI)

    println("solving....")
    statusM = solve(centralModel)
    @assert statusM==:Optimal "Central optimization not solved to optimality"

    uRaw=getvalue(u)
    snRaw=getvalue(sn)
    xtRaw=getvalue(xt)
    zRaw=getvalue(z)

    #calculate actual temp
    Tactual=zeros(horzLen+1,1)
    itotal=zeros(horzLen+1,1)
    for k=1:horzLen+1
        itotal[k,1]=sum((uRaw[(k-1)*N+n,1]) for n=1:N) + evS.w[(k-1)*2+1,1]
    end
    Tactual[1,1]=evS.τP*xt0+evS.γP*itotal[1,1]^2+evS.ρP*evS.w[2,1]
    for k=1:horzLen
        Tactual[k+1,1]=evS.τP*Tactual[k,1]+evS.γP*itotal[k+1,1]^2+evS.ρP*evS.w[k*2+2,1]
    end
    lambdaCurr=-getdual(currCon)
    if noTlimit==0
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

    cSol=centralSolutionStruct(xt=xtRaw,un=uRaw,sn=snRaw,z=zRaw,
                        itotal=itotal,uSum=uSum,zSum=zSum,
                        objVal=getobjectivevalue(centralModel),
                        lamTemp=lambdaUpperT,lamCoupl=lambdaCurr,
                        Tactual=Tactual)
    return cSol
end

#dual
function pwlEVdual(N::Int,S::Int,horzLen::Int,maxIt::Int,updateMethod::String,evS::scenarioStruct,cSol::centralSolutionStruct)

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
    	alpha = 5e5 #for kA
    	alphaDivRate=2
    	#alphaRate=.99
    	minAlpha=1e-6
    end

    stepI = 1;
    convChk = 1e-16
    convIt=maxIt

    #u w and z are one index ahead of x. i.e the x[k+1]=x[k]+eta*u[k+1]
    alphaP=alpha*ones(maxIt+1,1)

    #initialize with guess
    lambda0=1000*ones(horzLen+1,1)
    #lambda0=lamCurrStar
    dLog.Lam[:,1]=lambda0

    #iterate at each time step until convergence
    for p=1:maxIt
        #solve subproblem for each EV
    	@sync @parallel for evInd=1:N
            target=zeros((horzLen+1),1)
    		target[(evS.Kn[evInd,1]-(stepI-1)):1:length(target),1]=evS.Snmin[evInd,1]
            evM=Model(solver = GurobiSolver(NumericFocus=1))
            @variable(evM,un[1:horzLen+1])
            @variable(evM,sn[1:horzLen+1])
            @objective(evM,Min,sum((sn[k,1]-1)^2*evS.Qsi[evInd,1]+(un[k,1])^2*evS.Ri[evInd,1]+dLog.Lam[k,p]*un[k,1] for k=1:horzLen+1))
    		@constraint(evM,sn[1,1]==sn0[evInd,1]+evS.ηP[evInd,1]*un[1,1]) #fix for MPC loop
    		@constraint(evM,[k=1:horzLen],sn[k+1,1]==sn[k,1]+evS.ηP[evInd,1]*un[k+1,1]) #check K+1
            @constraint(evM,sn.<=1)
            @constraint(evM,sn.>=target)
            @constraint(evM,un.<=evS.imax[evInd,1])
            @constraint(evM,un.>=evS.imin[evInd,1])

    		TT = STDOUT # save original STDOUT stream
    		redirect_stdout()
            statusEVM = solve(evM)
    		redirect_stdout(TT)

    		@assert statusEVM==:Optimal "EV NLP optimization not solved to optimality"

            dLog.Sn[collect(evInd:N:N*(horzLen+1)),p+1]=getvalue(sn) #solved state goes in next time slot
            dLog.Un[collect(evInd:N:N*(horzLen+1)),p+1]=getvalue(un) #current go
        end

    	if updateMethod=="dualAscent"
    	    #solve coordinator problem
    	    #coorM=Model(solver = GurobiSolver(Presolve=0,NumericFocus=1))
    		coorM=Model(solver = GurobiSolver())
    		#coorM=Model(solver = IpoptSolver())
    	    @variable(coorM,z[1:S*(horzLen+1)])
    	    @variable(coorM,xt[1:horzLen+1])
    	    @objective(coorM,Min,sum(dLog.Lam[k,p]*sum(-z[(k-1)*S+s,1] for s=1:S) for k=1:(horzLen+1)))
    		@constraint(coorM,xt[1,1]==evS.τP*xt0+evS.γP*evS.deltaI*sum((2*s-1)*z[s,1] for s=1:S)+evS.ρP*evS.w[2,1]) #fix for MPC loop
    		@constraint(coorM,[k=1:horzLen],xt[k+1,1]==evS.τP*xt[k,1]+evS.γP*evS.deltaI*sum((2*s-1)*z[k*S+s,1] for s=1:S)+evS.ρP*evS.w[k*2+2,1])
    		if noTlimit==0
    			@constraint(coorM,upperTCon,xt.<=evS.Tmax)
    		end
    	    @constraint(coorM,xt.>=0)
    	    @constraint(coorM,z.<=evS.deltaI)
    	    @constraint(coorM,z.>=0)
    		TT = STDOUT # save original STDOUT stream
    		redirect_stdout()
    	    statusC = solve(coorM)
    		redirect_stdout(TT)

    		@assert statusC==:Optimal "Dual Ascent central optimization not solved to optimality"

    		 dLog.Xt[:,p+1]=getvalue(xt)
    		 dLog.Z[:,p+1]=getvalue(z)

    	    #grad of lagragian
    		gradL=zeros(horzLen+1,1)
    		for k=1:horzLen+1
    			dLog.zSum[k,p+1]=sum(dLog.Z[(k-1)*(S)+s,p+1] for s=1:S)
    			dLog.uSum[k,p+1]=sum(dLog.Un[(k-1)*N+n,p+1] for n=1:N)
    			gradL[k,1]=dLog.uSum[k,p+1] + evS.w[(k-1)*2+(stepI*2-1),1] - dLog.zSum[k,p+1]
    		end
    		dCM.couplConst[p,1]=norm(gradL,2)
    	end

    	#calculate actual temperature from nonlinear model of XFRM
    	ztotal=zeros(horzLen+1,1)
    	for k=1:horzLen+1
    		ztotal[k,1]=sum(dLog.Un[(k-1)*N+n,p+1]    for n=1:N) + evS.w[(k-1)*2+(stepI*2-1),1]
    	end
    	dLog.Tactual[1,p+1]=evS.τP*xt0+evS.γP*ztotal[1,1]^2+evS.ρP*evS.w[2,1] #fix for mpc
    	for k=1:horzLen
    		dLog.Tactual[k+1,p+1]=evS.τP*dLog.Tactual[k,p+1]+evS.γP*ztotal[k+1,1]^2+evS.ρP*evS.w[k*2+2,1]  #fix for mpc
    	end

    	if updateMethod=="fastAscent"
    		#fast ascent
    		if noTlimit==0
    			gradL=dLog.Tactual[:,p+1]-evS.Tmax*ones(horzLen+1,1)
    		else
    			gradL=zeros(horzLen+1,1)
    		end
    		#add some amount of future lambda
    		for k=1:(horzLen+1-2)
    			gradL[k,1]=.6*gradL[k,1]+.3*gradL[k+1,1]+.1*gradL[k+2,1]
    			#gradL[k,1]=.5*gradL[k,1]+.2*gradL[k+1,1]+.2*gradL[k+2,1]+.1*gradL[k+3,1]+.1*gradL[k+4,1]
    		end
    	end

        #update lambda
    	alphaP[p+1,1] = max(alpha/ceil(p/alphaDivRate),minAlpha)
    	#alphaP[p+1,1] = alphaP[p,1]*alphaRate

    	#lambda_new=lambda+alpha_p*gradL
        dLog.Lam[:,p+1]=max.(dLog.Lam[:,p]+alphaP[p+1,1]*gradL,0)

    	#check convergence
    	objFun(sn,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,1]     for n=1:N) for k=1:horzLen+1) +
    					sum(sum((u[(k-1)*N+n,1])^2*evS.Ri[n,1]           for n=1:N) for k=1:horzLen+1)
    	fGap= objFun(dLog.Sn[:,p+1],dLog.Un[:,p+1])-cSol.objVal
    	snGap=norm((dLog.Sn[:,p+1]-cSol.sn),2)
    	unGap=norm((dLog.Un[:,p+1]-cSol.un),2)
    	itGap = norm(dLog.Lam[:,p+1]-dLog.Lam[:,p],2)
    	if updateMethod=="fastAscent"
    		convGap = norm(dLog.Lam[:,p+1]-cSol.lamTemp,2)
    	else
    		convGap = norm(dLog.Lam[:,p+1]-cSol.lamCoupl,2)
    	end
    	dCM.objVal[p,1]=abs(fGap)
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
function pwlEVadmm(N::Int,S::Int,horzLen::Int,maxIt::Int,evS::scenarioStruct,cSol::centralSolutionStruct)
    #initialize
    sn0=evS.s0
    xt0=evS.t0

    stepI = 1;
    convChk = 1e-16
    convIt=maxIt

    #admm  initial parameters and guesses
    #ρADMM=10.0^(0)
    ρADMM=10^6 #for kA
    #ρADMM=1    #for A

    #u w and z are one index ahead of x. i.e the x[k+1]=x[k]+η*u[k+1]
    dCMadmm=convMetricsStruct()
    dLogadmm=itLogPWL()

    #initialize with guess
    lambda0=1000*ones(horzLen+1,1)
    vz0=-ones(S*(horzLen+1),1)
    vu0=.01*ones(N*(horzLen+1),1)
    #vz0=-zStar
    #vu0=uStar
    #lambda0=lamCurrStar
    dLogadmm.Lam[:,1]=lambda0
    dLogadmm.Vz[:,1]=vz0
    dLogadmm.Vu[:,1]=vu0

    for p in 1:maxIt
    	#ρ_p = ρADMM/ceil(p/2)
    	ρI = ρADMM
        #x minimization eq 7.66 in Bertsekas
        @sync @parallel for evInd=1:N
    		lambda=dLogadmm.Lam[:,p]
            evV=dLogadmm.Vu[collect(evInd:N:length(dLogadmm.Vu[:,p])),p]
            target=zeros((horzLen+1),1)
    		target[(evS.Kn[evInd,1]-(stepI-1)):1:length(target),1]=evS.Snmin[evInd,1]
        	evM = Model(solver = GurobiSolver())
        	@variable(evM,sn[1:(horzLen+1)])
        	@variable(evM,u[1:(horzLen+1)])
    		@objective(evM,Min, sum((sn[k,1]-1)^2*evS.Qsi[evInd,1]+(u[k,1])^2*evS.Ri[evInd,1]+
    								lambda[k,1]*(u[k,1]-evV[k,1])+
    								ρI/2*(u[k,1]-evV[k,1])^2 for k=1:horzLen+1))
            @constraint(evM,sn[1,1]==sn0[evInd,1]+evS.ηP[evInd,1]*u[1,1])
            @constraint(evM,[k=1:horzLen],sn[k+1,1]==sn[k,1]+evS.ηP[evInd,1]*u[k+1,1])
        	@constraint(evM,sn.<=1)
        	@constraint(evM,sn.>=target)
            @constraint(evM,u.<=evS.imax[evInd,1])
            @constraint(evM,u.>=evS.imin[evInd,1])
        	TT = STDOUT # save original STDOUT stream
        	redirect_stdout()
        	statusEVM = solve(evM)
        	redirect_stdout(TT)
    		@assert statusEVM==:Optimal "ADMM EV NLP optimization not solved to optimality"

    		dLogadmm.Sn[collect(evInd:N:length(dLogadmm.Sn[:,p+1])),p+1]=getvalue(sn)
    		dLogadmm.Un[collect(evInd:N:length(dLogadmm.Un[:,p+1])),p+1]=getvalue(u)
        end

        #N+1 decoupled problem aka transformer current
        tM = Model(solver = GurobiSolver())
        @variable(tM,z[1:(S)*(horzLen+1)])
        @variable(tM,xt[1:(horzLen+1)])
        # constFun1(u,v)=sum(Lam[k,p]*sum(u[(k-1)*(S)+s,1]-v[(k-1)*(S)+s,1] for s=1:S)  for k=1:(horzLen+1))
        # constFun2(u,v)=ρI/2*sum(sum((u[(k-1)*(S)+s,1]-v[(k-1)*(S)+s,1])*(u[(k-1)*(S)+s,1]-v[(k-1)*(S)+s,1]) for s=1:S)  for k=1:(horzLen+1))
        # @objective(tM,Min, constFun1(-z,Vz[:,p])+constFun2(-z,Vz[:,p]))
    	@objective(tM,Min,sum(dLogadmm.Lam[k,p]*(sum(-z[(k-1)*(S)+s,1] for s=1:S)-sum(dLogadmm.Vz[(k-1)*(S)+s,p] for s=1:S)) +
    						ρI/2*(sum(-z[(k-1)*(S)+s,1] for s=1:S)-sum(dLogadmm.Vz[(k-1)*(S)+s,p] for s=1:S))^2  for k=1:(horzLen+1)))
        @constraint(tM,tempCon1,xt[1,1]==evS.τP*xt0+evS.γP*evS.deltaI*sum((2*m+1)*z[m+1,1] for m=0:S-1)+evS.ρP*evS.w[stepI*2,1])
        @constraint(tM,tempCon2[k=1:horzLen],xt[k+1,1]==evS.τP*xt[k,1]+evS.γP*evS.deltaI*sum((2*m+1)*z[k*S+(m+1),1] for m=0:S-1)+evS.ρP*evS.w[stepI*2+k*2,1])
        if noTlimit==0
        	@constraint(tM,upperTCon,xt.<=evS.Tmax)
        end
        @constraint(tM,xt.>=0)
        @constraint(tM,z.>=0)
        @constraint(tM,z.<=evS.deltaI)
        #@constraint(tM,zC[k=1:horzLen+1],zSum[k,1]==sum(z[(k-1)*(S)+s] for s=1:S))

        TT = STDOUT # save original STDOUT stream
        redirect_stdout()
        statusC = solve(tM)
        redirect_stdout(TT)
    	@assert statusC==:Optimal "ADMM XFRM NLP optimization not solved to optimality"

        dLogadmm.Xt[:,p+1]=getvalue(xt)
        dLogadmm.Z[:,p+1]=getvalue(z)

        #lambda update eq 7.68
    	for k=1:horzLen+1
    		dLogadmm.uSum[k,p+1]=sum(dLogadmm.Un[(k-1)*N+n,p+1] for n=1:N)
    		dLogadmm.zSum[k,p+1]=sum(dLogadmm.Z[(k-1)*(S)+s,p+1] for s=1:S)
    		dLogadmm.couplConst[k,p+1]= dLogadmm.uSum[k,p+1] + evS.w[(k-1)*2+(stepI*2-1),1] - dLogadmm.zSum[k,p+1]
    		#Lam[k,p+1]=max.(Lam[k,p]+ρADMM/(horzLen+1)*(currConst[k,1]),0)
    		dLogadmm.Lam[k,p+1]=dLogadmm.Lam[k,p]+ρI/(S*(N))*(dLogadmm.couplConst[k,p+1])
    	end

    	#calculate actual temperature from nonlinear model of XFRM
    	dLogadmm.Tactual[1,p+1]=evS.τP*xt0+evS.γP*dLogadmm.zSum[1,p+1]^2+evS.ρP*evS.w[2,1] #fix for mpc
    	for k=1:horzLen
    		dLogadmm.Tactual[k+1,p+1]=evS.τP*dLogadmm.Tactual[k,p+1]+evS.γP*dLogadmm.zSum[k+1,p+1]^2+evS.ρP*evS.w[k*2+2,1]  #fix for mpc
    	end

        #v upate eq 7.67
        for k=1:horzLen+1
            dLogadmm.Vu[(k-1)*N+collect(1:N),p+1]=min.(max.(dLogadmm.Un[(k-1)*N+collect(1:N),p+1]+(dLogadmm.Lam[k,p]-dLogadmm.Lam[k,p+1])/ρI,evS.imin),evS.imax)
            dLogadmm.Vz[(k-1)*(S)+collect(1:S),p+1]=max.(min.(-dLogadmm.Z[(k-1)*(S)+collect(1:S),p+1]+(dLogadmm.Lam[k,p]-dLogadmm.Lam[k,p+1])/ρI,0),-evS.deltaI)
        end

        #check convergence
    	objFun(sn,xt,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,1]     for n=1:N) for k=1:horzLen+1) +
    					sum((xt[k,1]-1)^2*evS.Qsi[N+1,1]                 for k=1:horzLen+1) +
    					sum(sum((u[(k-1)*N+n,1])^2*evS.Ri[n,1]           for n=1:N) for k=1:horzLen+1)
    	fGap= abs(objFun(dLogadmm.Sn[:,p+1],dLogadmm.Xt[:,p+1],dLogadmm.Un[:,p+1])-cSol.objVal)
    	#fGap= objFun(Sn[:,p],Xt[:,p],Un[:,p])-fStar
    	snGap=norm((dLogadmm.Sn[:,p+1]-cSol.sn),2)
    	unGap=norm((dLogadmm.Un[:,p+1]-cSol.un),2)
    	constGap=norm(dLogadmm.couplConst[:,p+1],2)
    	itGap = norm(dLogadmm.Lam[:,p+1]-dLogadmm.Lam[:,p],2)
    	convGap = norm(dLogadmm.Lam[:,p+1]-cSol.lamCoupl,2)
    	dCMadmm.objVal[p,1]=fGap
    	dCMadmm.sn[p,1]=snGap
    	dCMadmm.un[p,1]=unGap
    	dCMadmm.couplConst[p,1]=constGap
    	dCMadmm.lamIt[p,1]=itGap
    	dCMadmm.lam[p,1]=convGap
    	if(itGap <= convChk )
    		@printf "Converged after %g iterations\n" p
    		convIt=p+1
    		break
    	else
    		@printf "lastGap  %e after %g iterations\n" itGap p
    		@printf "convGap  %e after %g iterations\n" convGap p
    		@printf "constGap %e after %g iterations\n" constGap p
            @printf "snGap    %e after %g iterations\n" snGap p
    		@printf "unGap    %e after %g iterations\n" unGap p
    		@printf("fGap     %e after %g iterations\n\n",fGap,p)
    	end
    end

    return (dLogadmm,dCMadmm,convIt)
end

#ALADIN
function pwlEValad(N::Int,S::Int,horzLen::Int,maxIt::Int,evS::scenarioStruct,cSol::centralSolutionStruct)
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
    dLogalad.Vu[:,1]=vu0
    dLogalad.Vs[:,1]=vs0
    dLogalad.Vz[:,1]=vz0
    dLogalad.Vt[:,1]=vt0
    dLogalad.Lam[:,1]=lambda0

    ΔY=zeros(1,maxIt+1)
    ρALADp=ρALAD*ones(1,maxIt+1)

    for p=1:maxIt

        #solve decoupled
        @sync @parallel for evInd=1:N
            ind=[evInd]
            for k=1:horzLen
                append!(ind,k*N+evInd)
            end
            lambda=dLogalad.Lam[:,p]
            evVu=dLogalad.Vu[ind,p]
            evVs=dLogalad.Vs[ind,p]
            #evV=zeros(horzLen+1,1)
            target=zeros((horzLen+1),1)
            target[(evS.Kn[evInd,1]-(stepI-1)):1:length(target),1]=evS.Snmin[evInd,1]
            evM = Model(solver = GurobiSolver())
            @variable(evM,sn[1:(horzLen+1)])
            @variable(evM,u[1:(horzLen+1)])
            @objective(evM,Min,sum((sn[k,1]-1)^2*evS.Qsi[evInd,1]+(u[k,1])^2*evS.Ri[evInd,1]+
                                    lambda[k,1]*(u[k,1])+
                                    ρALADp[1,p]/2*(u[k,1]-evVu[k,1])*σU[evInd,1]*(u[k,1]-evVu[k,1])+
                                    ρALADp[1,p]/2*(sn[k,1]-evVs[k,1])*σS[evInd,1]*(sn[k,1]-evVs[k,1]) for k=1:horzLen+1))
            @constraint(evM,sn[1,1]==sn0[evInd,1]+evS.ηP[evInd,1]*u[1,1])
            @constraint(evM,[k=1:horzLen],sn[k+1,1]==sn[k,1]+evS.ηP[evInd,1]*u[k+1,1])
            @constraint(evM,socKappaMax,sn.<=1)
            @constraint(evM,socKappaMin,sn.>=target)
            @constraint(evM,curKappaMax,u.<=evS.imax[evInd,1])
            @constraint(evM,curKappaMin,u.>=evS.imin[evInd,1])


            TT = STDOUT # save original STDOUT stream
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

            cValMax=abs.(uVal-evS.imax[evInd,1]).<tolU
            cValMin=abs.(uVal-evS.imin[evInd,1]).<tolU
            # cVal=kappaMax
            # cVal[cVal.>0]=1
            # cVal=kappaMin
            # cVal[cVal.<0]=-1
            dLogalad.Cu[ind,p+1]=1cValMax-1cValMin


            cValMax=abs.(snVal-1).<tolS
            cValMin=abs.(snVal-target).<tolS
            # cVal=socMax
            # cVal[cVal.>0]=1
            # cVal=socMin
            # cVal[cVal.<0]=-1
            dLogalad.Cs[ind,p+1]=1cValMax-1cValMin

            dLogalad.Sn[ind,p+1]=snVal
    		dLogalad.Un[ind,p+1]=uVal
            dLogalad.Gu[ind,p+1]=2*evS.Ri[evInd,1]*uVal
            dLogalad.Gs[ind,p+1]=2*evS.Qsi[evInd,1]*snVal-2*evS.Qsi[evInd,1]
            #Gu[collect(evInd:N:length(Gu[:,p+1])),p+1]=σU[evInd,1]*(evVu-uVal)+lambda
            #Gs[collect(evInd:N:length(Gs[:,p+1])),p+1]=σN[evInd,1]*(evVs-snVal)-lambda
        end

        #N+1 decoupled problem aka transformer current
        tM = Model(solver = GurobiSolver())
        #tM = Model(solver = IpoptSolver())
        @variable(tM,z[1:(S)*(horzLen+1)])
        @variable(tM,xt[1:(horzLen+1)])
        @objective(tM,Min, sum(-dLogalad.Lam[k,p]*sum(z[(k-1)*(S)+s] for s=1:S)+
                  ρALADp[1,p]/2*σZ*(sum(z[(k-1)*(S)+s] for s=1:S)-sum(dLogalad.Vz[(k-1)*(S)+s,p] for s=1:S))^2+
                  ρALADp[1,p]/2*σT*(xt[k]-dLogalad.Vt[k,p])^2  for k=1:(horzLen+1)))
        @constraint(tM,tempCon1,xt[1]-evS.τP*xt0-evS.γP*evS.deltaI*sum((2*s-1)*z[s] for s=1:S)-evS.ρP*evS.w[stepI*2,1]==0)
        @constraint(tM,tempCon2[k=1:horzLen],xt[k+1]-evS.τP*xt[k]-evS.γP*evS.deltaI*sum((2*s-1)*z[(k)*(S)+s] for s=1:S)-evS.ρP*evS.w[stepI*2+k*2,1]==0)
        if noTlimit==0
        	@constraint(tM,upperTCon,xt.<=evS.Tmax)
        end
        @constraint(tM,lowerTCon,xt.>=0)
        @constraint(tM,pwlKappaMin,z.>=0)
        @constraint(tM,pwlKappaMax,z.<=evS.deltaI)
        TT = STDOUT # save original STDOUT stream
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

        cValMax=abs.(zVal-evS.deltaI).<tolZ
        cValMin=abs.(zVal-0).<tolZ
        # cVal=kappaMax
        # cVal=kappaMin
        # cVal[cVal.<0]=-1
        # cVal[cVal.>0]=1
        dLogalad.Cz[:,p+1]=1cValMax-1cValMin

        cValMax=abs.(xtVal-evS.Tmax).<tolT
        cValMin=abs.(xtVal-0).<tolT
        # cVal=tMin
        # cVal[cVal.<0]=-1
        # cVal=tMax
        # cVal[cVal.>0]=1
        dLogalad.Ct[:,p+1]=1cValMax-1cValMin

        dLogalad.Xt[:,p+1]=xtVal
        dLogalad.Z[:,p+1]=zVal
        dLogalad.Gz[:,p+1]=0
        #Gz[:,p+1]=σZ*(Vz[:,p]-zVal)-repeat(-Lam[:,p],inner=S)

        for k=1:horzLen+1
            dLogalad.uSum[k,p+1]=sum(dLogalad.Un[(k-1)*N+n,p+1] for n=1:N)
            dLogalad.zSum[k,p+1]=sum(dLogalad.Z[(k-1)*(S)+s,p+1] for s=1:S)
            dLogalad.couplConst[k,p+1]=dLogalad.uSum[k,p+1] + evS.w[(k-1)*2+(stepI*2-1),1] - dLogalad.zSum[k,p+1]
        end

        #check for convergence
        constGap=norm(dLogalad.couplConst[:,p+1],1)
        cc=norm(vcat((dLogalad.Vu[:,p]-dLogalad.Un[:,p+1]),(dLogalad.Vz[:,p]-dLogalad.Z[:,p+1])),1)
        #convCheck=ρALAD*norm(vcat(repmat(σU,horzLen+1,1).*(Vu[:,p]-Un[:,p+1]),σZ*(Vz[:,p]-Z[:,p+1])),1)
        objFun(sn,xt,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,1]     for n=1:N) for k=1:horzLen+1) +
                        sum((xt[k,1]-1)^2*evS.Qsi[N+1,1]                 for k=1:horzLen+1) +
                        sum(sum((u[(k-1)*N+n,1])^2*evS.Ri[n,1]           for n=1:N) for k=1:horzLen+1)
        fGap= abs(objFun(dLogalad.Sn[:,p+1],dLogalad.Xt[:,p+1],dLogalad.Un[:,p+1])-cSol.objVal)
        snGap=norm((dLogalad.Sn[:,p+1]-cSol.sn),2)
        unGap=norm((dLogalad.Un[:,p+1]-cSol.un),2)
        itGap = norm(dLogalad.Lam[:,p]-dLogalad.Lam[:,max(p-1,1)],2)
        convGap = norm(dLogalad.Lam[:,p]-cSol.lamCoupl,2)
        dCMalad.objVal[p,1]=fGap
        dCMalad.sn[p,1]=snGap
        dCMalad.un[p,1]=unGap
        dCMalad.lamIt[p,1]=itGap
        dCMalad.couplConst[p,1]=constGap
        dCMalad.lam[p,1]=convGap
        convCheck[p,1]=cc
        if  constGap<=epsilon && cc<=epsilon
            @printf "Converged after %g iterations\n" p
            convIt=p+1
            break
        else
            @printf "lastGap    %e after %g iterations\n" itGap p
            @printf "convLamGap %e after %g iterations\n" convGap p
            @printf "convCheck  %e after %g iterations\n" cc p
            @printf "constGap   %e after %g iterations\n" constGap p
            @printf "snGap      %e after %g iterations\n" snGap p
            @printf("fGap       %e after %g iterations\n\n",fGap,p)
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
    	@objective(cM,Min, sum(sum(0.5*dUn[(k-1)*N+i,1]^2*Hu[i,1]+dLogalad.Gu[(k-1)*N+i,p+1]*dUn[(k-1)*N+i,1] +
                                   0.5*dSn[(k-1)*N+i,1]^2*Hs[i,1]+dLogalad.Gs[(k-1)*N+i,p+1]*dSn[(k-1)*N+i,1] for i=1:N) +
                                   0.5*dZ[k,1]^2*Hz+
                                   0.5*dXt[k,1]^2*Ht   for k=1:(horzLen+1)))
        # @constraint(cM,currCon[k=1:horzLen+1],w[(k-1)*2+1]+relaxS[k,1]==-sum(Un[(k-1)*(N)+n,p+1]+dUn[(k-1)*(N)+n,1] for n=1:N)+
        #                                          sum(Z[(k-1)*(S)+s,p+1]+dZ[(k-1)*(S)+s,1] for s=1:S))
        @constraint(cM,currCon[k=1:horzLen+1],sum(dLogalad.Un[(k-1)*(N)+n,p+1]+dUn[(k-1)*(N)+n,1] for n=1:N)-
                                                 sum(dLogalad.Z[(k-1)*(S)+s,p+1]+dZ[(k-1)*(S)+s,1] for s=1:S)==-evS.w[(k-1)*2+1])#+relaxS[k,1])
        #local equality constraints C*(X+deltaX)=0 is same as C*deltaX=0 since we already know CX=0
        @constraint(cM,stateCon1[n=1:N],dSn[n,1]==evS.ηP[n,1]*dUn[n,1])
        @constraint(cM,stateCon2[k=1:horzLen,n=1:N],dSn[n+(k)*(N),1]==dSn[n+(k-1)*(N),1]+evS.ηP[n,1]*dUn[n+(k)*(N),1])
        @constraint(cM,tempCon1,dXt[1,1]==evS.γP*evS.deltaI*sum((2*s-1)*dZ[s,1] for s=1:S))
        @constraint(cM,tempCon2[k=1:horzLen],dXt[k+1,1]==evS.τP*dXt[k,1]+evS.γP*evS.deltaI*sum((2*s-1)*dZ[k*S+s,1] for s=1:S))

        #local inequality constraints
        # @constraint(cM,(Z[:,p+1]+dZ).>=0)
        # @constraint(cM,(Z[:,p+1]+dZ).<=deltaI)
     	# @constraint(cM,(Un[:,p+1]+dUn).<=repmat(imax,horzLen+1,1))
     	# @constraint(cM,(Un[:,p+1]+dUn).>=repmat(imin,horzLen+1,1))

        #these shouldnt be elementwise?????
        # @constraint(cM,Cz[:,p+1].*dZ.==0)
        # @constraint(cM,Cu[:,p+1].*dUn.==0)
        # @constraint(cM,Cn[:,p+1].*dSn.==0)
        # @constraint(cM,Ct[:,p+1].*dXt.==0)

        @constraint(cM,dLogalad.Cz[:,p+1].*dZ.<=0)
        @constraint(cM,dLogalad.Cu[:,p+1].*dUn.<=0)
        @constraint(cM,dLogalad.Cs[:,p+1].*dSn.<=0)
        @constraint(cM,dLogalad.Ct[:,p+1].*dXt.<=0)


    	TT = STDOUT # save original STDOUT stream
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
        dLogalad.Lam[:,p+1]=max.(dLogalad.Lam[:,p]+α3*(-getdual(currCon)-dLogalad.Lam[:,p]),0)

        dLogalad.Vu[:,p+1]=dLogalad.Vu[:,p]+α1*(dLogalad.Un[:,p+1]-dLogalad.Vu[:,p])+α2*getvalue(dUn)
        dLogalad.Vz[:,p+1]=dLogalad.Vz[:,p]+α1*(dLogalad.Z[:,p+1]-dLogalad.Vz[:,p])+α2*getvalue(dZ)
        dLogalad.Vs[:,p+1]=dLogalad.Vs[:,p]+α1*(dLogalad.Sn[:,p+1]-dLogalad.Vs[:,p])+α2*getvalue(dSn)
        dLogalad.Vt[:,p+1]=dLogalad.Vt[:,p]+α1*(dLogalad.Xt[:,p+1]-dLogalad.Vt[:,p])+α2*getvalue(dXt)

        # Vu[:,p+1]=Un[:,p+1]+getvalue(dUn)
        # Vz[:,p+1]=Z[:,p+1]+getvalue(dZ)
        # Vs[:,p+1]=Sn[:,p+1]+getvalue(dSn)
        # Vt[:,p+1]=Xt[:,p+1]+getvalue(dXt)

        ρALADp[1,p+1]=min(ρALADp[1,p]*ρRate,1e6) #increase ρ every iteration
        ΔY[1,p+1]=norm(vcat(getvalue(dUn),getvalue(dZ),getvalue(dSn),getvalue(dXt)),Inf)
    end

    return dLogalad,dCMalad,convIt,ΔY,convCheck
end
