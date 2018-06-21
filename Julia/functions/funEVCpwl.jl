#functions to run PWL EVC problems


#central
function pwlEVcentral(N,S,horzLen,evS::scenarioStruct)
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
