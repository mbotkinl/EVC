#functions to run NL EVC problems


#central
function nlEVcentral(N::Int,S::Int,horzLen::Int,evS::scenarioStruct, relaxed=false)

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
    cModel = Model(solver = IpoptSolver())
    #u w and z are one index ahead of x. i.e the x[k+1]=x[k]+eta*u[k+1]
    @variable(cModel,u[1:N*(horzLen+1)])
    @variable(cModel,sn[1:(N)*(horzLen+1)])
    @variable(cModel,xt[1:(horzLen+1)])
    @variable(cModel,itotal[1:(horzLen+1)])

    println("obj")
    @objective(cModel,Min, sum(sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,1]+
                                   (u[(k-1)*N+n,1])^2*evS.Ri[n,1] for n=1:N) for k=1:(horzLen+1)))

    println("constraints")
    @constraint(cModel,stateCon1,sn[1:N,1].==sn0[1:N,1]+evS.ηP[:,1].*u[1:N,1])
    @constraint(cModel,stateCon2[k=1:horzLen,n=1:N],sn[n+(k)*(N),1]==sn[n+(k-1)*(N),1]+evS.ηP[n,1]*u[n+(k)*(N),1])
    if relaxed==true
        @constraint(cModel,tempCon1,xt[1,1]>=evS.τP*xt0+evS.γP*(itotal[1])^2+evS.ρP*evS.w[stepI*2,1])
        @constraint(cModel,tempCon2[k=1:horzLen],xt[k+1,1]>=evS.τP*xt[k,1]+evS.γP*(itotal[k+1])^2+evS.ρP*evS.w[stepI*2+k*2,1]) #check id index???
    else
        @NLconstraint(cModel,tempCon1,xt[1,1]==evS.τP*xt0+evS.γP*(itotal[1])^2+evS.ρP*evS.w[stepI*2,1])
        @NLconstraint(cModel,tempCon2[k=1:horzLen],xt[k+1,1]==evS.τP*xt[k,1]+evS.γP*(itotal[k+1])^2+evS.ρP*evS.w[stepI*2+k*2,1]) #check id index???
    end
    @constraint(cModel,currCon[k=1:horzLen+1],0==-sum(u[(k-1)*(N)+n,1] for n=1:N)-evS.w[(k-1)*2+1]+itotal[k]) #fix for MPC
    @constraint(cModel,sn.<=1)
    @constraint(cModel,sn.>=target)
    if noTlimit==false
    	@constraint(cModel,upperTCon,xt.<=evS.Tmax)
    end
    @constraint(cModel,xt.>=0)
    @constraint(cModel,upperCCon,u.<=repmat(evS.imax,horzLen+1,1))
    @constraint(cModel,u.>=repmat(evS.imin,horzLen+1,1))
    @constraint(cModel,itotal.<=evS.ItotalMax)
    @constraint(cModel,itotal.>=0)

    println("solving....")
    statusC = solve(cModel)
    @assert statusC==:Optimal "Central NL optimization not solved to optimality"


    uRaw=getvalue(u)
    snRaw=getvalue(sn)
    xtRaw=getvalue(xt)
    itotalRaw=getvalue(itotal)

    if noTlimit==false
    	kappaUpperT=-getdual(upperTCon)
    else
    	kappaUpperT=zeros(horzLen+1,1)
    end
    lambdaCurr=-getdual(currCon)
    if relaxed==true
        lambdaTemp=zeros(horzLen+1,1)
    else
        lambdaTemp=[-getdual(tempCon1);-getdual(tempCon2)]
    end

    uSum=zeros(horzLen+1,1)
    for k=1:horzLen+1
        uSum[k,1]=sum(uRaw[(k-1)*N+n,1] for n=1:N)
    end

    cSol=centralSolutionStruct(xt=xtRaw,un=uRaw,sn=snRaw,
                        itotal=itotalRaw,uSum=uSum,
                        objVal=getobjectivevalue(cModel),
                        lamTemp=lambdaTemp,lamCoupl=lambdaCurr)
    return cSol
end

#dual
function nlEVdual(N::Int,S::Int,horzLen::Int,maxIt::Int,updateMethod::String,
    evS::scenarioStruct,cSolnl::centralSolutionStruct, relaxed=false)

    #initialize with current states
    sn0=evS.s0
    xt0=evS.t0

    if updateMethod=="fastAscent"
    	#alpha = 0.1 #for A
    	alpha = 5e3 #for kA
    	alphaDivRate=2
    	minAlpha=1e-6
    	#alphaRate=.99
    else
    	#alpha = .01 #for A
    	alpha = 1e4 #for kA
    	alphaDivRate=4
    	minAlpha=1e-6
    	#alphaRate=.99
    end

    stepI = 1
    convChk = 1e-8
    convIt=maxIt

    alphaP=alpha*ones(maxIt+1,1)

    dCM=convMetricsStruct()
    dLog=itLogNL()

    #u w and z are one index ahead of x. i.e the x[k+1]=x[k]+eta*u[k+1]
    lambda0=1000*ones(horzLen+1,1)
    #lambda0=lamCurrStarNL
    #lambda0=max.(lamCurrStarNL,0)
    dLog.Lam[:,1]=lambda0

    #iterate at each time step until convergence
    for p=1:maxIt
        #solve subproblem for each EV
    	@sync @parallel for evInd=1:N
    		ind=[evInd]
    		for k=1:horzLen
    			append!(ind,k*N+evInd)
    		end
            target=zeros((horzLen+1),1)
    		target[(evS.Kn[evInd,1]-(stepI-1)):1:length(target),1]=evS.Snmin[evInd,1]
            evM=Model(solver = IpoptSolver())
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

    		@assert statusEVM==:Optimal "EV NL optimization not solved to optimality"

            dLog.Sn[ind,p+1]=getvalue(sn) #solved state goes in next time slot
            dLog.Un[ind,p+1]=getvalue(un) #current go
        end

    	if updateMethod=="dualAscent"
    	    #solve coordinator problem

            coorM = Model(solver = IpoptSolver())
    	    @variable(coorM,itotal[1:(horzLen+1)])
    	    @variable(coorM,xt[1:(horzLen+1)])
    	    @objective(coorM,Min,-sum(dLog.Lam[k,p]*itotal[k,1] for k=1:(horzLen+1)))
            if relaxed==true
                @constraint(coorM,xt[1,1]>=evS.τP*xt0+evS.γP*(itotal[1])^2+evS.ρP*evS.w[2,1]) #fix for MPC loop
        		@constraint(coorM,[k=1:horzLen],xt[k+1,1]>=evS.τP*xt[k,1]+evS.γP*(itotal[k+1,1])^2+evS.ρP*evS.w[k*2+2,1])
            else
                @NLconstraint(coorM,xt[1,1]==evS.τP*xt0+evS.γP*(itotal[1])^2+evS.ρP*evS.w[2,1]) #fix for MPC loop
        		@NLconstraint(coorM,[k=1:horzLen],xt[k+1,1]==evS.τP*xt[k,1]+evS.γP*(itotal[k+1,1])^2+evS.ρP*evS.w[k*2+2,1])
            end
    		if noTlimit==false
    			@constraint(coorM,upperTCon,xt.<=evS.Tmax)
    		end
    	    @constraint(coorM,xt.>=0)
    	    @constraint(coorM,itotal.<=evS.ItotalMax)
    	    @constraint(coorM,itotal.>=0)
    		TT = STDOUT # save original STDOUT strea
    		redirect_stdout()
    	    statusC = solve(coorM)
    		redirect_stdout(TT)

    		@assert statusC==:Optimal "XFRM NL optimization not solved to optimality"

    		 dLog.Xt[:,p+1]=getvalue(xt)
    		 dLog.Itotal[:,p+1]=getvalue(itotal)

    	    #grad of lagragian
    		gradL=zeros(horzLen+1,1)
    		for k=1:horzLen+1
    			dLog.uSum[k,p+1]=sum(dLog.Un[(k-1)*N+n,p+1] for n=1:N)
    			gradL[k,1]=dLog.uSum[k,p+1] + evS.w[(k-1)*2+(stepI*2-1),1] - dLog.Itotal[k,p+1]
    		end
    		dLog.couplConst[p,1]=norm(gradL,2)
    	end

    	if updateMethod=="fastAscent"
            ztotal=zeros(horzLen+1,1)
            for k=1:horzLen+1
                ztotal[k,1]=sum(dLog.Un[(k-1)*N+n,p+1]    for n=1:N) + evS.w[(k-1)*2+(stepI*2-1),1]
            end
            dLog.Xt[1,p+1]=evS.τP*xt0+evS.γP*ztotal[1,1]^2+evS.ρP*evS.w[2,1] #fix for mpc
            for k=1:horzLen
                dLog.Xt[k+1,p+1]=evS.τP*dLog.Xt[k,p+1]+evS.γP*ztotal[k+1,1]^2+evS.ρP*evS.w[k*2+2,1]  #fix for mpc
            end

    		#fast ascent
    		if noTlimit==false
    			gradL=dLog.Xt[:,p+1]-evS.Tmax*ones(horzLen+1,1)
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

    	#Lam[:,p+1]=Lam[:,p]+alpha_p*gradL
        dLog.Lam[:,p+1]=max.(dLog.Lam[:,p]+alphaP[p+1,1]*gradL,0)

    	#check convergence
    	objFun(sn,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,1]     for n=1:N) for k=1:horzLen+1) +
    					sum(sum((u[(k-1)*N+n,1])^2*evS.Ri[n,1]           for n=1:N) for k=1:horzLen+1)
    	fGap= objFun(dLog.Sn[:,p+1],dLog.Un[:,p+1])-cSolnl.objVal
    	snGap=norm((dLog.Sn[:,p+1]-cSolnl.sn),2)
    	unGap=norm((dLog.Un[:,p+1]-cSolnl.un),2)
    	itGap = norm(dLog.Lam[:,p+1]-dLog.Lam[:,p],2)
    	if updateMethod=="fastAscent"
    		convGap = norm(dLog.Lam[:,p+1]-cSolnl.lamTemp,2)
    	else
    		convGap = norm(dLog.Lam[:,p+1]-cSolnl.lamCoupl,2)
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
    return (dLog,dCM,convIt)
end

#admm
function nlEVadmm(N::Int,S::Int,horzLen::Int,maxIt::Int,evS::scenarioStruct,cSolnl::centralSolutionStruct, relaxed=false)
    #initialize with current states
    sn0=evS.s0
    xt0=evS.t0

    stepI = 1
    convChk = 1e-8
    convIt=maxIt

    #admm  initial parameters and guesses
    #ρADMM=10.0^(0)
    ρADMM=5e5
    ρDivRate=10

    #u w and z are one index ahead of x. i.e the x[k+1]=x[k]+eta*u[k+1]
    dCMadmm=convMetricsStruct()
    dLogadmm=itLogNL()

    # lambda0=lamCurrStarNL
    # vi0=-itotalStarNL
    # vu0=uStarNL
    lambda0=1000*ones(horzLen+1,1)
    vi0=-ones(horzLen+1,1)
    vu0=.01*ones(N*(horzLen+1),1)
    dLogadmm.Vi[:,1]=vi0
    dLogadmm.Vu[:,1]=vu0
    dLogadmm.Lam[:,1]=lambda0

    for p in 1:maxIt
        ρADMMp = ρADMM/ceil(p/ρDivRate)
        #ρADMMp = ρADMM

        #x minimization eq 7.66 in Bertsekas
        @sync @parallel for evInd=1:N
            evV=dLogadmm.Vu[collect(evInd:N:length(dLogadmm.Vu[:,p])),p]
            target=zeros((horzLen+1),1)
    		target[(evS.Kn[evInd,1]-(stepI-1)):1:length(target),1]=evS.Snmin[evInd,1]
            evM = Model(solver = IpoptSolver())
        	@variable(evM,sn[1:(horzLen+1)])
        	@variable(evM,u[1:(horzLen+1)])
    		@objective(evM,Min, sum((sn[k,1]-1)^2*evS.Qsi[evInd,1]+(u[k,1])^2*evS.Ri[evInd,1]+
                                    dLogadmm.Lam[k,p]*(u[k,1]-evV[k,1])+
                                    ρADMMp/2*(u[k,1]-evV[k,1])^2 for k=1:horzLen+1))
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
            @assert statusEVM==:Optimal "EV NLP NL optimization not solved to optimality"

    		dLogadmm.Sn[collect(evInd:N:length(dLogadmm.Sn[:,p+1])),p+1]=getvalue(sn)
    		dLogadmm.Un[collect(evInd:N:length(dLogadmm.Un[:,p+1])),p+1]=getvalue(u)
        end

        #N+1 decoupled problem aka transformer current
        tM = Model(solver = IpoptSolver())
        #@variable(tM,z[1:(S)*(horzLen+1)])
        @variable(tM,xt[1:(horzLen+1)])
        @variable(tM,itotal[1:(horzLen+1)])
    	# constFun1(u,v)=sum(Lam[k,p]*(u[k,1]-v[k,1])  for k=1:(horzLen+1))
    	# constFun2(u,v)=ρADMM/2*sum((u[k,1]-v[k,1])^2  for k=1:(horzLen+1))
        # @objective(tM,Min, constFun1(-itotal,Vi[:,p])+constFun2(-itotal,Vi[:,p]))
        @objective(tM,Min,sum(dLogadmm.Lam[k,p]*(-itotal[k,1]-dLogadmm.Vi[k,p])+
                              ρADMM/2*(-itotal[k,1]-dLogadmm.Vi[k,p])^2  for k=1:(horzLen+1)))
        if relaxed==true
            @constraint(tM,tempCon1,xt[1,1]>=evS.τP*xt0+evS.γP*(itotal[1,1])^2+evS.ρP*evS.w[stepI*2,1])
            @constraint(tM,tempCon2[k=1:horzLen],xt[k+1,1]>=evS.τP*xt[k,1]+evS.γP*(itotal[k+1,1])^2+evS.ρP*evS.w[stepI*2+k*2,1])
        else
            @NLconstraint(tM,tempCon1,xt[1,1]==evS.τP*xt0+evS.γP*(itotal[1,1])^2+evS.ρP*evS.w[stepI*2,1])
            @NLconstraint(tM,tempCon2[k=1:horzLen],xt[k+1,1]==evS.τP*xt[k,1]+evS.γP*(itotal[k+1,1])^2+evS.ρP*evS.w[stepI*2+k*2,1])
        end
        if noTlimit==false
        	@constraint(tM,upperTCon,xt.<=evS.Tmax)
        end
        @constraint(tM,xt.>=0)
        @constraint(tM,itotal.>=0)
        @constraint(tM,itotal.<=evS.ItotalMax)

        TT = STDOUT # save original STDOUT stream
        redirect_stdout()
        statusTM = solve(tM)
        redirect_stdout(TT)
        @assert statusTM==:Optimal "XFRM NLP NL optimization not solved to optimality"

        dLogadmm.Xt[:,p+1]=getvalue(xt)
        dLogadmm.Itotal[:,p+1]=getvalue(itotal)

        #lambda update eq 7.68
    	for k=1:horzLen+1
    		dLogadmm.uSum[k,p+1]=sum(dLogadmm.Un[(k-1)*N+n,p+1] for n=1:N)
    		dLogadmm.couplConst[k,p+1]=dLogadmm.uSum[k,p+1] + evS.w[(k-1)*2+(stepI*2-1),1] - dLogadmm.Itotal[k,p+1]
            dLogadmm.Lam[k,p+1]=max.(dLogadmm.Lam[k,p]+ρADMMp/(N+1)*(dLogadmm.couplConst[k,p+1]),0)
    		# dLogadmm.Lam[k,p+1]=dLogadmm.Lam[k,p]+ρADMMp/(N+1)*(dLogadmm.couplConst[k,p+1])
    	end

        #v upate eq 7.67
        for k=1:horzLen+1
            dLogadmm.Vu[(k-1)*N+collect(1:N),p+1]=min.(max.(dLogadmm.Un[(k-1)*N+collect(1:N),p+1]+(dLogadmm.Lam[k,p]-dLogadmm.Lam[k,p+1])/ρADMMp,evS.imin),evS.imax)
    		dLogadmm.Vi[k,p+1]=max.(min.(-dLogadmm.Itotal[k,p+1]+(dLogadmm.Lam[k,p]-dLogadmm.Lam[k,p+1])/ρADMMp,0),-evS.ItotalMax)
        end

        #check convergence
    	objFun(sn,xt,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,1]     for n=1:N) for k=1:horzLen+1) +
    					sum((xt[k,1]-1)^2*evS.Qsi[N+1,1]                 for k=1:horzLen+1) +
    					sum(sum((u[(k-1)*N+n,1])^2*evS.Ri[n,1]           for n=1:N) for k=1:horzLen+1)
    	fGap= objFun(dLogadmm.Sn[:,p+1],dLogadmm.Xt[:,p+1],dLogadmm.Un[:,p+1])-cSolnl.objVal
    	#fGap= objFun(Sn[:,p],Xt[:,p],Un[:,p])-fStar
    	snGap=norm((dLogadmm.Sn[:,p+1]-cSolnl.sn),2)
        unGap=norm((dLogadmm.Un[:,p+1]-cSolnl.un),2)
    	constGap=norm(dLogadmm.couplConst[:,p+1],2)
    	itGap = norm(dLogadmm.Lam[:,p+1]-dLogadmm.Lam[:,p],2)
    	convGap = norm(dLogadmm.Lam[:,p+1]-cSolnl.lamCoupl,2)
    	dCMadmm.objVal[p,1]=abs(fGap)
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
    		@printf("fGap     %e after %g iterations\n\n",fGap,p)

    	end
    end

    return (dLogadmm,dCMadmm,convIt)
end

#aladin
function nlEValad(N::Int,S::Int,horzLen::Int,maxIt::Int,evS::scenarioStruct,cSolnl::centralSolutionStruct, relaxed=false)
    #initialize with current states
    sn0=evS.s0
    xt0=evS.t0

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
    σS=ones(N,1)/10
    σI=1/N
    σT=1/10000
    Hu=2*evS.Ri#*(1+rand())
    Hs=2*evS.Qsi#*(1+rand())
    Hi=1e-6
    Ht=1e-6
    ρALAD=1
    ρRate=1.15
    muALAD=10^8

    #.1/.1 looks good expect for const
    #.01/2 increasing ρ helps const
    #.1/1.1 best so far???

    dCMalad=convMetricsStruct()
    dLogalad=itLogNL()
    convCheck=zeros(maxIt+1,1)

    lambda0=1000*ones(horzLen+1,1)
    vt0=ones(horzLen+1,1)
    vi0=ones((horzLen+1),1)
    vu0=.01*ones(N*(horzLen+1),1)
    vs0=.5*ones(N*(horzLen+1),1)
    # lambda0=lamCurrStarNL
    # vt0=xtStarNL
    # vi0=itotalStarNL
    # vu0=uStarNL
    # vs0=snStarNL
    dLogalad.Vu[:,1]=vu0 #initial guess goes in column 1
    dLogalad.Vs[:,1]=vs0 #initial guess goes in column 1
    dLogalad.Vi[:,1]=vi0
    dLogalad.Vt[:,1]=vt0
    dLogalad.Lam[:,1]=lambda0

    ρALADp=ρALAD*ones(1,maxIt+1)
    ΔY=zeros(1,maxIt+1)

    for p=1:maxIt
        @printf "Starting iteration %g \n" p

        #solve decoupled
        @sync @parallel for evInd=1:N
            ind=[evInd]
            for k=1:horzLen
                append!(ind,k*N+evInd)
            end
            lambda=dLogalad.Lam[:,p]
            evVu=dLogalad.Vu[ind,p]
            evVs=dLogalad.Vs[ind,p]
            target=zeros((horzLen+1),1)
            target[(evS.Kn[evInd,1]-(stepI-1)):1:length(target),1]=evS.Snmin[evInd,1]
            #evM = Model(solver = GurobiSolver(NumericFocus=3))
            evM = Model(solver = IpoptSolver())
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
            @assert statusEVM==:Optimal "ALAD EV NLP NL optimization not solved to optimality"

    		# kappaMax=-getdual(curKappaMax)
    		# kappaMin=-getdual(curKappaMin)
            # socMax=-getdual(socKappaMax)
            # socMin=-getdual(socKappaMin)
            uVal=getvalue(u)
            snVal=getvalue(sn)

            cValMax=abs.(uVal-evS.imax[evInd,1]).<tolU
            cValMin=abs.(uVal-evS.imin[evInd,1]).<tolU
            dLogalad.Cu[ind,p+1]=1cValMax-1cValMin
            # cVal=zeros(length(ind),1)
            # cVal[kappaMax.>0]=1
            # cVal[kappaMin.<0]=-1
            # Cu[ind,p+1]=cVal

            cValMax=abs.(snVal-1).<tolS
            cValMin=abs.(snVal-target).<tolS
            dLogalad.Cs[ind,p+1]=1cValMax-1cValMin
            # cVal=zeros(length(ind),1)
            # cVal[socMax.>0]=1
            # cVal[socMin.<0]=-1
            # Cs[ind,p+1]=cVal

            dLogalad.Sn[ind,p+1]=snVal
    		dLogalad.Un[ind,p+1]=uVal
            dLogalad.Gu[ind,p+1]=2*evS.Ri[evInd,1]*uVal
            dLogalad.Gs[ind,p+1]=2*evS.Qsi[evInd,1]*snVal-2*evS.Qsi[evInd,1]
            #Gu[collect(evInd:N:length(Gu[:,p+1])),p+1]=σU[evInd,1]*(evVu-uVal)+lambda
            #Gs[collect(evInd:N:length(Gs[:,p+1])),p+1]=σN[evInd,1]*(evVs-snVal)-lambda
        end

        #N+1 decoupled problem aka transformer current
        tM = Model(solver = IpoptSolver())
        @variable(tM,itotal[1:(horzLen+1)])
        @variable(tM,xt[1:(horzLen+1)])
        @objective(tM,Min, sum(-dLogalad.Lam[k,p]*itotal[k]+
                  ρALADp[1,p]/2*σI*(itotal[k]-dLogalad.Vi[k,p])^2+
                  ρALADp[1,p]/2*σT*(xt[k]-dLogalad.Vt[k,p])^2  for k=1:(horzLen+1)))
        if relaxed
            @constraint(tM,tempCon1,xt[1]-evS.τP*xt0-evS.γP*(itotal[1])^2-evS.ρP*evS.w[stepI*2,1]>=0)
            @constraint(tM,tempCon2[k=1:horzLen],xt[k+1]-evS.τP*xt[k]-evS.γP*(itotal[k+1])^2-evS.ρP*evS.w[stepI*2+k*2,1]>=0)
        else
            @NLconstraint(tM,tempCon1,xt[1]-evS.τP*xt0-evS.γP*(itotal[1])^2-evS.ρP*evS.w[stepI*2,1]==0)
            @NLconstraint(tM,tempCon2[k=1:horzLen],xt[k+1]-evS.τP*xt[k]-evS.γP*(itotal[k+1])^2-evS.ρP*evS.w[stepI*2+k*2,1]==0)
        end
        if noTlimit==false
        	@constraint(tM,upperTCon,xt.<=evS.Tmax)
        end
        @constraint(tM,lowerTCon,xt.>=0)
        @constraint(tM,KappaMin,itotal.>=0)
        @constraint(tM,KappaMax,itotal.<=evS.ItotalMax)
        TT = STDOUT # save original STDOUT stream
        redirect_stdout()
        statusTM = solve(tM)
        redirect_stdout(TT)
        @assert statusTM==:Optimal "ALAD XFRM NL optimization not solved to optimality"

        # kappaMax=-getdual(KappaMin)
        # kappaMin=-getdual(KappaMax)
        # tMax=-getdual(upperTCon)
        # tMin=-getdual(lowerTCon)
        lambdaTemp=[-getdual(tempCon1);-getdual(tempCon2)]
        iVal=getvalue(itotal)
        xtVal=getvalue(xt)

        cValMax=abs.(iVal-evS.ItotalMax).<tolI
        cValMin=abs.(iVal-0).<tolI
        dLogalad.Ci[:,p+1]=1cValMax-1cValMin
        # cVal=zeros(length(ind),1)
        # cVal[kappaMax.>0]=1
        # cVal[kappaMin.<0]=-1
        # Ci[:,p+1]=cVal


        cValMax=abs.(xtVal-evS.Tmax).<tolT
        cValMin=abs.(xtVal-0).<tolT
        dLogalad.Ct[:,p+1]=1cValMax-1cValMin
        # cVal=zeros(length(ind),1)
        # cVal[tMax.>0]=1
        # cVal[tMin.<0]=-1
        # Ct[:,p+1]=cVal

        dLogalad.Xt[:,p+1]=xtVal
        dLogalad.Itotal[:,p+1]=iVal
        dLogalad.Gi[:,p+1]=0
        #Gz[:,p+1]=σZ*(Vz[:,p]-zVal)-repeat(-Lam[:,p],inner=S)

        for k=1:horzLen+1
            dLogalad.uSum[k,p+1]=sum(dLogalad.Un[(k-1)*N+n,p+1] for n=1:N)
            dLogalad.couplConst[k,p+1]=dLogalad.uSum[k,p+1] + evS.w[(k-1)*2+(stepI*2-1),1] - dLogalad.Itotal[k,p+1]
        end

        #check for convergence
        constGap=norm(dLogalad.couplConst[:,p+1],1)
        cc=norm(vcat((dLogalad.Vu[:,p]-dLogalad.Un[:,p+1]),(dLogalad.Vi[:,p]-dLogalad.Itotal[:,p+1])),1)
        #convCheck=ρALADp*norm(vcat(repmat(σU,horzLen+1,1).*(Vu[:,p]-Un[:,p+1]),σZ*(Vz[:,p]-Z[:,p+1])),1)
        objFun(sn,xt,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,1]     for n=1:N) for k=1:horzLen+1) +
                        sum((xt[k,1]-1)^2*evS.Qsi[N+1,1]                 for k=1:horzLen+1) +
                        sum(sum((u[(k-1)*N+n,1])^2*evS.Ri[n,1]           for n=1:N) for k=1:horzLen+1)
        fGap= abs(objFun(dLogalad.Sn[:,p+1],dLogalad.Xt[:,p+1],dLogalad.Un[:,p+1])-cSolnl.objVal)
        fGap2= abs((objFun(dLogalad.Sn[:,p+1],dLogalad.Xt[:,p+1],dLogalad.Un[:,p+1])-cSolnl.objVal)/cSolnl.objVal)
        snGap=norm((dLogalad.Sn[:,p+1]-cSolnl.sn),2)
        unGap=norm((dLogalad.Un[:,p+1]-cSolnl.un),2)
        itGap = norm(dLogalad.Lam[:,p]-dLogalad.Lam[:,max(p-1,1)],2)
        convGap = norm(dLogalad.Lam[:,p]-cSolnl.lamCoupl,2)
        convGap2 = norm(round.(dLogalad.Lam[:,p]-cSolnl.lamCoupl,10)./cSolnl.lamCoupl,2) #need some rounding here if both are 1e-8

        dCMalad.objVal[p,1]=fGap
        dCMalad.sn[p,1]=snGap
        dCMalad.un[p,1]=unGap
        dCMalad.lamIt[p,1]=itGap
        dCMalad.couplConst[p,1]=constGap
        dCMalad.lam[p,1]=convGap
        convCheck[p,1]=cc
        if  constGap<=epsilon && convCheck<=epsilon
            @printf "Converged after %g iterations\n" p
            convIt=p+1
            break
        else
            @printf "lastGap    %e after %g iterations\n" itGap p
            @printf "convLamGap %e after %g iterations\n" convGap p
            @printf "convCheck  %e after %g iterations\n" cc p
            @printf "constGap   %e after %g iterations\n" constGap p
            #@printf "snGap      %e after %g iterations\n" snGap p
            @printf("fGap       %e after %g iterations\n\n",fGap,p)
        end

        #coupled QP
        cM = Model(solver = IpoptSolver())
        @variable(cM,dUn[1:(N)*(horzLen+1)])
        @variable(cM,dSn[1:(N)*(horzLen+1)])
        @variable(cM,dI[1:(horzLen+1)])
        @variable(cM,dXt[1:(horzLen+1)])
        #@variable(cM,relaxS[1:(horzLen+1)])
        #objExp=objExp+Lam[:,p]'*relaxS+muALAD/2*sum(relaxS[k,1]^2 for k=1:horzLen+1)
    	@objective(cM,Min, sum(sum(0.5*dUn[(k-1)*N+i,1]^2*Hu[i,1]+dLogalad.Gu[(k-1)*N+i,p+1]*dUn[(k-1)*N+i,1]+
                                   0.5*dSn[(k-1)*N+i,1]^2*Hs[i,1]+dLogalad.Gs[(k-1)*N+i,p+1]*dSn[(k-1)*N+i,1] for i=1:N)+
                               0.5*dI[k,1]^2*(Hi-lambdaTemp[k,1]*2*evS.γP)+
                               0.5*dXt[k,1]^2*Ht   for k=1:(horzLen+1)))
        @constraint(cM,currCon[k=1:horzLen+1],sum(dLogalad.Un[(k-1)*(N)+n,p+1]+dUn[(k-1)*(N)+n,1] for n=1:N)-
                                                    (dLogalad.Itotal[k,p+1]+dI[k])==-evS.w[(k-1)*2+1])#+relaxS[k,1])
        @constraint(cM,stateCon1[n=1:N],dSn[n,1]==evS.ηP[n,1]*dUn[n,1])
        @constraint(cM,stateCon2[k=1:horzLen,n=1:N],dSn[n+(k)*(N),1]==dSn[n+(k-1)*(N),1]+evS.ηP[n,1]*dUn[n+(k)*(N),1])
        @constraint(cM,tempCon1,dXt[1,1]==2*evS.γP*dLogalad.Itotal[1,p+1]*dI[1])
        @constraint(cM,tempCon2[k=1:horzLen],dXt[k+1,1]==evS.τP*dXt[k,1]+2*evS.γP*dLogalad.Itotal[k+1,p+1]*dI[k+1,1])
        @constraint(cM,dLogalad.Ci[:,p+1].*dI.<=0)
        @constraint(cM,dLogalad.Cu[:,p+1].*dUn.<=0)
        @constraint(cM,dLogalad.Cs[:,p+1].*dSn.<=0)
        @constraint(cM,dLogalad.Ct[:,p+1].*dXt.<=0)

    	TT = STDOUT # save original STDOUT stream
        redirect_stdout()
        statusM = solve(cM)
        redirect_stdout(TT)
        @assert statusM==:Optimal "ALAD Central optimization not solved to optimality"

        #update step
        # Lam[:,p+1]=-getdual(currCon)
        α1=1
        α2=1
        α3=1
        #α1=α1/ceil(p/2)
        #Lam[:,p+1]=Lam[:,p]+alpha3*(-getdual(currCon)-Lam[:,p])
        dLogalad.Lam[:,p+1]=max.(dLogalad.Lam[:,p]+α3*(-getdual(currCon)-dLogalad.Lam[:,p]),0)

        dLogalad.Vu[:,p+1]=dLogalad.Vu[:,p]+α1*(dLogalad.Un[:,p+1]-dLogalad.Vu[:,p])+α2*getvalue(dUn)
        dLogalad.Vi[:,p+1]=dLogalad.Vi[:,p]+α1*(dLogalad.Itotal[:,p+1]-dLogalad.Vi[:,p])+α2*getvalue(dI)
        dLogalad.Vs[:,p+1]=dLogalad.Vs[:,p]+α1*(dLogalad.Sn[:,p+1]-dLogalad.Vs[:,p])+α2*getvalue(dSn)
        dLogalad.Vt[:,p+1]=dLogalad.Vt[:,p]+α1*(dLogalad.Xt[:,p+1]-dLogalad.Vt[:,p])+α2*getvalue(dXt)

        ρALADp[1,p+1]=min(ρALADp[1,p]*ρRate,1e6) #increase ρ every iteration
        ΔY[1,p+1]=norm(vcat(getvalue(dUn),getvalue(dI),getvalue(dSn),getvalue(dXt)),Inf)
    end

    return dLogalad,dCMalad,convIt,ΔY,convCheck
end
