#EVC with ALADIN for PWL convex relaxation
#Micah Botkin-Levy
#5/22/18

#formulation try 2
#u, sn, xt, and z are all in "y", v is the auxilliary variable corresponding to x in literature
#current constraint is coupling

tic()
#initialize with current states
sn0=s0
xt0=T0

stepI = 1;
horzLen=K1
epsilon = 1e-8
tolU=1e-6
tolS=1e-8
tolT=1e-4
tolI=1e-6
maxIt=50
convIt=maxIt
ConvALAD=zeros(maxIt,1)
constConvALAD=zeros(maxIt,1)
itConvALAD=zeros(maxIt,1)
fConvALAD=zeros(maxIt,1)
snConvALAD=zeros(maxIt,1)
avgN=zeros(maxIt,1)

#ALADIN tuning and initial guess
sigmaU=1*ones(N,1)
sigmaS=ones(N,1)/10
sigmaI=1/N
sigmaT=1/10000
Hz=1e-6
Ht=1e-6
rhoALAD=1e-1
rhoRate=1.1
muALAD=10^6
# lambda0=5*rand(Truncated(Normal(0), 0, 1), horzLen+1)
# vt0=Tmax*rand(Truncated(Normal(0), 0, 1), (horzLen+1))
# vz0=ItotalMax/1000*rand(Truncated(Normal(0), 0, 1), S*(horzLen+1))
# vu0=imax[1,1]*0.8*rand(Truncated(Normal(0), 0, 1), N*(horzLen+1))
# vs0=rand(Truncated(Normal(0), 0, 1), N*(horzLen+1))

lambda0=ones(horzLen+1,1)
vt0=ones(horzLen+1,1)
vz0=ones(S*(horzLen+1),1)
vu0=.01*ones(N*(horzLen+1),1)
vs0=.5*ones(N*(horzLen+1),1)

# lambda0=lamCurrStar
# vt0=xtStar
# vz0=zStar
# vu0=uStar
# vs0=snStar

#save matrices
Un=SharedArray{Float64}(N*(horzLen+1),maxIt) #row are time (N states for k=1, them N states for k=2),  columns are iteration
Sn=SharedArray{Float64}(N*(horzLen+1),maxIt)  #row are time,  columns are iteration
Xt=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration
Z=zeros(S*(horzLen+1),maxIt)  #row are time,  columns are iteration
Lam=zeros((horzLen+1),maxIt) #(rows are time, columns are iteration)
#Lam[:,1]=max.(lambda0,0)
Lam[:,1]=lambda0

rhoALADp=rhoALAD*ones(1,maxIt)

#auxillary variables
Vu=zeros((N)*(horzLen+1),maxIt) #row are time,  columns are iteration
Vu[:,1]=vu0 #initial guess goes in column 1
Vs=zeros((N)*(horzLen+1),maxIt) #row are time,  columns are iteration
Vs[:,1]=vs0 #initial guess goes in column 1
Vz=zeros(S*(horzLen+1),maxIt) #row are time,  columns are iteration
Vz[:,1]=vz0
Vt=zeros((horzLen+1),maxIt) #row are time,  columns are iteration
Vt[:,1]=vt0

#Gradian Vectors
Gu=SharedArray{Float64}(N*(horzLen+1),maxIt) #row are time (N states for k=1, them N states for k=2),  columns are iteration
Gs=SharedArray{Float64}(N*(horzLen+1),maxIt) #row are time (N states for k=1, them N states for k=2),  columns are iteration
Gz=zeros(S*(horzLen+1),maxIt) #row are time (N states for k=1, them N states for k=2),  columns are iteration

#Jacobian C Vectors
Cs=SharedArray{Float64}(N*(horzLen+1),maxIt)  #row are time,  columns are iteration
Cu=SharedArray{Float64}(N*(horzLen+1),maxIt)  #row are time,  columns are iteration
Cz=zeros(S*(horzLen+1),maxIt)  #row are time,  columns are iteration
Ct=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration

uSum=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration
zSum=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration
currConst=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration

for p=1:maxIt-1

    #solve decoupled
    @sync @parallel for evInd=1:N

        ind=[evInd]
        for k=1:horzLen
            append!(ind,k*N+evInd)
        end

        lambda=Lam[:,p]
        evVu=Vu[ind,p]
        evVs=Vs[ind,p]
        #evV=zeros(horzLen+1,1)
        #target=zeros((horzLen+1),1)
        #target[(Kn[evInd,1]-(stepI-1)):1:length(target),1]=Snmin[evInd,1]
        evM = Model(solver = GurobiSolver(NumericFocus=3))
        @variable(evM,sn[1:(horzLen+1)])
        @variable(evM,u[1:(horzLen+1)])
        #@objective(evM,Min, sum(objFun(sn,u)+sum(lambda[k,1]*(u[k,1]-evV[k,1]) for k=1:horzLen+1)+rho_p/2*sum((u[k,1]-evV[k,1])^2 for k=1:horzLen+1)))
        @objective(evM,Min,sum((sn[k,1]-1)^2*Qsi[evInd,1]+(u[k,1])^2*Ri[evInd,1]+
                                lambda[k,1]*(u[k,1])+
                                rhoALADp[1,p]/2*(u[k,1]-evVu[k,1])*sigmaU[evInd,1]*(u[k,1]-evVu[k,1])+
                                rhoALADp[1,p]/2*(sn[k,1]-evVs[k,1])*sigmaS[evInd,1]*(sn[k,1]-evVs[k,1]) for k=1:horzLen+1))
        @constraint(evM,sn[1,1]==sn0[evInd,1]+etaP[evInd,1]*u[1,1])
        @constraint(evM,[k=1:horzLen],sn[k+1,1]==sn[k,1]+etaP[evInd,1]*u[k+1,1])
        @constraint(evM,socKappaMax,sn.<=1)
        @constraint(evM,socKappaMin,sn.>=target[ind])
        @constraint(evM,curKappaMax,u.<=imax[evInd,1])
        @constraint(evM,curKappaMin,u.>=imin[evInd,1])


        TT = STDOUT # save original STDOUT stream
        redirect_stdout()
        status = solve(evM)
        redirect_stdout(TT)
        if status!=:Optimal
            println("solver issues with EV NLP")
            return
        else
			kappaMax=-getdual(curKappaMax)
			kappaMin=-getdual(curKappaMin)
            socMax=-getdual(socKappaMax)
            socMin=-getdual(socKappaMin)
            uVal=getvalue(u)
            snVal=getvalue(sn)

            cValMax=abs.(uVal-imax[evInd,1]).<tolU
            cValMin=abs.(uVal-imin[evInd,1]).<tolU
            # cVal=kappaMax
            # cVal[cVal.>0]=1
            # cVal=kappaMin
            # cVal[cVal.<0]=-1
            Cu[ind,p+1]=1cValMax-1cValMin


            cValMax=abs.(snVal-1).<tolS
            cValMin=abs.(snVal-target[ind]).<tolS
            # cVal=socMax
            # cVal[cVal.>0]=1
            # cVal=socMin
            # cVal[cVal.<0]=-1
            Cs[ind,p+1]=1cValMax-1cValMin

            Sn[ind,p+1]=snVal
    		Un[ind,p+1]=uVal
            Gu[ind,p+1]=2*Ri[evInd,1]*uVal
            Gs[ind,p+1]=2*Qsi[evInd,1]*snVal-2*Qsi[evInd,1]
            #Gu[collect(evInd:N:length(Gu[:,p+1])),p+1]=sigmaU[evInd,1]*(evVu-uVal)+lambda
            #Gs[collect(evInd:N:length(Gs[:,p+1])),p+1]=sigmaN[evInd,1]*(evVs-snVal)-lambda
        end
    end

    #N+1 decoupled problem aka transformer current
    tM = Model(solver = GurobiSolver(NumericFocus=3))
    #tM = Model(solver = IpoptSolver())
    @variable(tM,z[1:(S)*(horzLen+1)])
    @variable(tM,xt[1:(horzLen+1)])
    tMobj=sum(-Lam[k,p]*sum(z[(k-1)*(S)+s] for s=1:S)+
              rhoALADp[1,p]/2*sigmaZ*sum((z[(k-1)*(S)+s]-Vz[(k-1)*(S)+s,p])^2 for s=1:S)+
              rhoALADp[1,p]/2*sigmaT*(xt[k]-Vt[k,p])^2  for k=1:(horzLen+1))
    @objective(tM,Min, tMobj)
    @constraint(tM,tempCon1,xt[1]==tauP*xt0+gammaP*deltaI*sum((2*s-1)*z[s] for s=1:S)+rhoP*w[stepI*2,1])
    @constraint(tM,tempCon2[k=1:horzLen],xt[k+1]==tauP*xt[k]+gammaP*deltaI*sum((2*s-1)*z[(k)*(S)+s] for s=1:S)+rhoP*w[stepI*2+k*2,1])
    # @constraint(tM,tempCon1,xt[1]/gammaP==tauP*xt0/gammaP+deltaI*sum((2*s-1)*z[s] for s=1:S)+rhoP*w[stepI*2,1]/gammaP)
    # @constraint(tM,tempCon2[k=1:horzLen],xt[k+1]/gammaP==tauP*xt[k]/gammaP+deltaI*sum((2*s-1)*z[(k)*(S)+s] for s=1:S)+rhoP*w[stepI*2+k*2,1]/gammaP)
    if noTlimit==0
    	@constraint(tM,upperTCon,xt.<=Tmax)
    end
    @constraint(tM,lowerTCon,xt.>=0)
    @constraint(tM,pwlKappaMin,z.>=0)
    @constraint(tM,pwlKappaMax,z.<=deltaI)
    TT = STDOUT # save original STDOUT stream
    redirect_stdout()
    statusTM = solve(tM)
    redirect_stdout(TT)
    if statusTM!=:Optimal
        println("solver issues with XFRM NLP")
        return
    else
		#kappaMax=-getdual(pwlKappaMax)
		#kappaMin=-getdual(pwlKappaMin)
        #tMax=-getdual(upperTCon)
        #tMin=-getdual(lowerTCon)
        zVal=getvalue(z)
        xtVal=getvalue(xt)

        cValMax=abs.(zVal-deltaI).<tolZ
        cValMin=abs.(zVal-0).<tolZ
        # cVal=kappaMax
        # cVal=kappaMin
        # cVal[cVal.<0]=-1
        # cVal[cVal.>0]=1
        Cz[:,p+1]=1cValMax-1cValMin

        cValMax=abs.(xtVal-Tmax).<tolT
        cValMin=abs.(xtVal-0).<tolT
        # cVal=tMin
        # cVal[cVal.<0]=-1
        # cVal=tMax
        # cVal[cVal.>0]=1
        Ct[:,p+1]=1cValMax-1cValMin

        Xt[:,p+1]=xtVal
        Z[:,p+1]=zVal
        Gz[:,p+1]=0
        #Gz[:,p+1]=sigmaZ*(Vz[:,p]-zVal)-repeat(-Lam[:,p],inner=S)
    end

    for k=1:horzLen+1
        uSum[k,p+1]=sum(Un[(k-1)*N+n,p+1] for n=1:N)
        zSum[k,p+1]=sum(Z[(k-1)*(S)+s,p+1] for s=1:S)
        currConst[k,p+1]=uSum[k,p+1] + w[(k-1)*2+(stepI*2-1),1] - zSum[k,p+1]
    end

    #check for convergence
    constGap=norm(currConst[:,p+1],1)
    convCheck=norm(vcat((Vu[:,p]-Un[:,p+1]),(Vz[:,p]-Z[:,p+1])),1)
    #convCheck=rhoALAD*norm(vcat(repmat(sigmaU,horzLen+1,1).*(Vu[:,p]-Un[:,p+1]),sigmaZ*(Vz[:,p]-Z[:,p+1])),1)
    avgN[p,1]=mean(convCheck)
    objFun(sn,xt,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*Qsi[n,1]     for n=1:N) for k=1:horzLen+1) +
                    sum((xt[k,1]-1)^2*Qsi[N+1,1]                 for k=1:horzLen+1) +
                    sum(sum((u[(k-1)*N+n,1])^2*Ri[n,1]           for n=1:N) for k=1:horzLen+1)
    fGap= abs(objFun(Sn[:,p+1],Xt[:,p+1],Un[:,p+1])-fStar)
    snGap=norm((Sn[:,p+1]-snStar),2)
    itGap = norm(Lam[:,p]-Lam[:,max(p-1,1)],2)
    convGap = norm(Lam[:,p]-lamCurrStar,2)
    fConvALAD[p,1]=fGap
    snConvALAD[p,1]=snGap
    itConvALAD[p,1]=itGap
    constConvALAD[p,1]=constGap
    ConvALAD[p,1]=convGap
    if  constGap<=epsilon && convCheck<=epsilon
        @printf "Converged after %g iterations\n" p
        convIt=p+1
        break
    else
        @printf "lastGap    %e after %g iterations\n" itGap p
        @printf "convLamGap %e after %g iterations\n" convGap p
        @printf "convCheck  %e after %g iterations\n" norm(convCheck,1) p
        @printf "constGap   %e after %g iterations\n" constGap p
        @printf "snGap      %e after %g iterations\n" snGap p
        @printf("fGap       %e after %g iterations\n\n",fGap,p)
    end


    #coupled QP
    cM = Model(solver = GurobiSolver(Presolve=0,NumericFocus=3))
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
    objExp=sum(sum(0.5*dUn[(k-1)*N+i,1]^2*2*Ri[i,1]+Gu[(k-1)*N+i,p+1]*dUn[(k-1)*N+i,1]+
                   0.5*dSn[(k-1)*N+i,1]^2*2*Qsi[i,1]+Gs[(k-1)*N+i,p+1]*dSn[(k-1)*N+i,1] for i=1:N) for k=1:(horzLen+1))
    objExp=objExp+sum(0.5*dZ[k,1]^2*Hz+
                      0.5*dXt[k,1]^2*Ht   for k=1:(horzLen+1))
    #objExp=objExp+Lam[:,p]'*relaxS+muALAD/2*sum(relaxS[k,1]^2 for k=1:horzLen+1)
	@objective(cM,Min, objExp)
    # @constraint(cM,currCon[k=1:horzLen+1],w[(k-1)*2+1]+relaxS[k,1]==-sum(Un[(k-1)*(N)+n,p+1]+dUn[(k-1)*(N)+n,1] for n=1:N)+
    #                                          sum(Z[(k-1)*(S)+s,p+1]+dZ[(k-1)*(S)+s,1] for s=1:S))
    @constraint(cM,currCon[k=1:horzLen+1],sum(Un[(k-1)*(N)+n,p+1]+dUn[(k-1)*(N)+n,1] for n=1:N)-
                                             sum(Z[(k-1)*(S)+s,p+1]+dZ[(k-1)*(S)+s,1] for s=1:S)==-w[(k-1)*2+1])#+relaxS[k,1])
    #local equality constraints C*(X+deltaX)=0 is same as C*deltaX=0 since we already know CX=0
    @constraint(cM,stateCon1[n=1:N],dSn[n,1]==etaP[n,1]*dUn[n,1])
    @constraint(cM,stateCon2[k=1:horzLen,n=1:N],dSn[n+(k)*(N),1]==dSn[n+(k-1)*(N),1]+etaP[n,1]*dUn[n+(k)*(N),1])
    @constraint(cM,tempCon1,dXt[1,1]==gammaP*deltaI*sum((2*s-1)*dZ[s,1] for s=1:S))
    @constraint(cM,tempCon2[k=1:horzLen],dXt[k+1,1]==tauP*dXt[k,1]+gammaP*deltaI*sum((2*s-1)*dZ[k*S+s,1] for s=1:S))

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

    @constraint(cM,Cz[:,p+1].*dZ.<=0)
    @constraint(cM,Cu[:,p+1].*dUn.<=0)
    @constraint(cM,Cs[:,p+1].*dSn.<=0)
    @constraint(cM,Ct[:,p+1].*dXt.<=0)


	TT = STDOUT # save original STDOUT stream
    redirect_stdout()
    statusM = solve(cM)
    redirect_stdout(TT)
    if statusM!=:Optimal
        println("solver issues with Central QP")
        break
    end

    #update step
    # Lam[:,p+1]=-getdual(currCon)
    alpha1=1
    alpha2=1
    alpha3=1
    #alpha3=alpha3/ceil(p/2)
    #Lam[:,p+1]=Lam[:,p]+alpha3*(-getdual(currCon)-Lam[:,p])
    Lam[:,p+1]=max.(Lam[:,p]+alpha3*(-getdual(currCon)-Lam[:,p]),0)

    Vu[:,p+1]=Vu[:,p]+alpha1*(Un[:,p+1]-Vu[:,p])+alpha2*getvalue(dUn)
    Vz[:,p+1]=Vz[:,p]+alpha1*(Z[:,p+1]-Vz[:,p])+alpha2*getvalue(dZ)
    Vs[:,p+1]=Vs[:,p]+alpha1*(Sn[:,p+1]-Vs[:,p])+alpha2*getvalue(dSn)
    Vt[:,p+1]=Vt[:,p]+alpha1*(Xt[:,p+1]-Vt[:,p])+alpha2*getvalue(dXt)

    # Vu[:,p+1]=Un[:,p+1]+getvalue(dUn)
    # Vz[:,p+1]=Z[:,p+1]+getvalue(dZ)
    # Vs[:,p+1]=Sn[:,p+1]+getvalue(dSn)
    # Vt[:,p+1]=Xt[:,p+1]+getvalue(dXt)

    rhoALADp[1,p+1]=rhoALADp[1,p]*rhoRate #increase rho every iteration
end

println("plotting....")
xPlot=zeros(horzLen+1,N)
for ii= 1:N
	xPlot[:,ii]=Sn[collect(ii:N:length(Sn[:,convIt])),convIt]
end

#plot(x=1:horzLen+1,y=xPlot2[:,ii])
# plot(x=1:Kn[ii,1],y=xPlot2[1:Kn[ii,1],ii])

pd1alad=plot(xPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV SOC"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_decentral_ALADIN_SOC.png", 24inch, 12inch), pd1alad) end


uPlot=zeros(horzLen+1,N)
for ii= 1:N
	uPlot[:,ii]=Un[collect(ii:N:length(Un[:,convIt])),convIt]
end

pd2alad=plot(uPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV Current (A)"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_decentral_ALADIN_Curr.png", 24inch, 12inch), pd2alad) end

pd3alad=plot(layer(x=1:horzLen+1,y=Xt[:,convIt],Geom.line,Theme(default_color=colorant"blue")),
		#layer(x=1:horzLen+1,y=Tactual,Geom.line,Theme(default_color=colorant"green")),
		yintercept=[Tmax],Geom.hline(color=["red"],style=:dot),
		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",key_position = :top,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
		#Guide.manual_color_key("", ["PWL Temp", "Actual Temp"], ["blue", "green"]))
if drawFig==1 draw(PNG(path*"J_decentral_ALADIN_Temp.png", 24inch, 12inch), pd3alad) end

pd4alad=plot(layer(x=1:horzLen+1,y=Lam[:,convIt],Geom.line),
		layer(x=1:horzLen+1,y=lamCurrStar,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
		Guide.xlabel("Time"), Guide.ylabel(raw"Lambda ($/A)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_decentral_ALADIN_Lam.png", 24inch, 12inch), pd4alad) end


lamPlotalad=plot(Lam[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=lamCurrStar,Geom.line,Theme(default_color=colorant"black",line_width=4pt)),
			Guide.xlabel("Time"), Guide.ylabel("Lambda"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
uSumPlotalad=plot(uSum[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=uSumStar,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("U sum"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
zSumPlotalad=plot(zSum[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=zSumStar,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("Z sum"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
constPlotalad2=plot(currConst[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel("curr constraint diff"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_ALADIN_LamConv.png", 36inch, 12inch), lamPlotalad) end

convItPlotalad=plot(x=1:convIt,y=itConvALAD[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
convPlotalad=plot(x=1:convIt,y=ConvALAD[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("central lambda gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
constPlotalad=plot(x=1:convIt,y=constConvALAD[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("central lambda gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
fPlotalad=plot(x=1:convIt-1,y=fConvALAD[1:convIt-1,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("obj function gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_ALADIN_Conv.png", 36inch, 12inch), vstack(convItPlotalad,convPlotalad,fPlotalad)) end
