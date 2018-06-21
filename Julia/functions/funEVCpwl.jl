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

    cSol=centralSolution(xt=xtRaw,un=uRaw,sn=snRaw,z=zRaw,
                        itotal=itotal,uSum=uSum,zSum=zSum,
                        objVal=getobjectivevalue(centralModel),
                        lamTemp=lambdaUpperT,lamCoupl=lambdaCurr,
                        Tactual=Tactual)
    return cSol
end

#dual
function pwlEVdual(N::Int,S::Int,horzLen::Int,maxIt::Int,updateMethod,evS::scenarioStruct,cSol::centralSolution)

    #initialize
    sn0=evS.s0
    xt0=evS.t0
    dCM=convMetrics()
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
    alphaP=alpha*ones(maxIt,1)

    #initialize with guess
    lambda0=1000*ones(horzLen+1,1)
    #lambda0=lamCurrStar
    dLog.Lam[:,1]=lambda0

    #iterate at each time step until convergence
    for p=1:maxIt-1
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



#ALADIN
