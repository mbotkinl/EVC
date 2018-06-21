#functions to run NL EVC problems


#central
function nlEVcentral(N,S,horzLen,evS::scenarioStruct)

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
    @NLconstraint(cModel,tempCon1,xt[1,1]==evS.τP*xt0+evS.γP*(itotal[1])^2+evS.ρP*evS.w[stepI*2,1])
    @NLconstraint(cModel,tempCon2[k=1:horzLen],xt[k+1,1]==evS.τP*xt[k,1]+evS.γP*(itotal[k+1])^2+evS.ρP*evS.w[stepI*2+k*2,1]) #check id index???
    @constraint(cModel,currCon[k=1:horzLen+1],0==-sum(u[(k-1)*(N)+n,1] for n=1:N)-evS.w[(k-1)*2+1]+itotal[k]) #fix for MPC
    @constraint(cModel,sn.<=1)
    @constraint(cModel,sn.>=target)
    if noTlimit==0
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

    toc()

    uRaw=getvalue(u)
    snRaw=getvalue(sn)
    xtRaw=getvalue(xt)
    itotalRaw=getvalue(itotal)

    if noTlimit==0
    	kappaUpperT=-getdual(upperTCon)
    else
    	kappaUpperT=zeros(horzLen+1,1)
    end
    lambdaCurr=-getdual(currCon)
    lambdaTemp=[-getdual(tempCon1);-getdual(tempCon2)]

    uSum=zeros(horzLen+1,1)
    for k=1:horzLen+1
        uSum[k,1]=sum(uRaw[(k-1)*N+n,1] for n=1:N)
    end

    cSol=centralSolution(xt=xtRaw,un=uRaw,sn=snRaw,
                        itotal=itotalRaw,uSum=uSum,
                        objVal=getobjectivevalue(cModel),
                        lamTemp=lambdaTemp,lamCoupl=lambdaCurr)
    return cSol
end
