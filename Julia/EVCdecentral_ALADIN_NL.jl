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
epsilon = 1e-12
tolU=1e-4
tolS=1e-8
tolT=1e-4
tolI=1e-6
# tolU=1e-4
# tolS=1e-8
# tolT=1e-4
# tolI=1e-4
maxIt=50
convIt=maxIt
ConvALAD=zeros(maxIt,1)
constConvALAD=zeros(maxIt,1)
itConvALAD=zeros(maxIt,1)
fConvALAD=zeros(maxIt,1)
snConvALAD=zeros(maxIt,1)
convCheck=zeros(maxIt,1)

#ALADIN tuning and initial guess
##sigmas are tuned to i_n [.010kA]
sigmaU=1*ones(N,1)
sigmaS=ones(N,1)/10
sigmaI=1/N
sigmaT=1/10000

##sigmas are tuned to i_n [10A]
# sigmaU=1*ones(N,1)
# sigmaS=ones(N,1)*100
# sigmaI=1/N
# sigmaT=1/10

Hi=1e-6
Ht=1e-6
rhoALAD=1
rhoRate=1.15
muALAD=10^8

# convItPlotaladNL
# convPlotaladNL
# constPlotaladNL
# fPlotaladNL

#.1/.1 looks good expect for const
#.01/2 increasing rho helps const
#.1/1.1 best so far???

lambda0=ones(horzLen+1,1)
vt0=ones(horzLen+1,1)
vi0=ones((horzLen+1),1)
vu0=.01*ones(N*(horzLen+1),1)
vs0=.5*ones(N*(horzLen+1),1)

# lambda0=lamCurrStarNL
# vt0=xtStarNL
# vi0=itotalStarNL
# vu0=uStarNL
# vs0=snStarNL

#save matrices
Un=SharedArray{Float64}(N*(horzLen+1),maxIt) #row are time (N states for k=1, them N states for k=2),  columns are iteration
Sn=SharedArray{Float64}(N*(horzLen+1),maxIt)  #row are time,  columns are iteration
Xt=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration
Itotal=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration
Lam=zeros((horzLen+1),maxIt) #(rows are time, columns are iteration)
#Lam[:,1]=max.(lambda0,0)
Lam[:,1]=lambda0

rhoALADp=rhoALAD*ones(1,maxIt)

#auxillary variables
Vu=zeros((N)*(horzLen+1),maxIt) #row are time,  columns are iteration
Vu[:,1]=vu0 #initial guess goes in column 1
Vs=zeros((N)*(horzLen+1),maxIt) #row are time,  columns are iteration
Vs[:,1]=vs0 #initial guess goes in column 1
Vi=zeros((horzLen+1),maxIt) #row are time,  columns are iteration
Vi[:,1]=vi0
Vt=zeros((horzLen+1),maxIt) #row are time,  columns are iteration
Vt[:,1]=vt0

#Gradian Vectors
Gu=SharedArray{Float64}(N*(horzLen+1),maxIt) #row are time (N states for k=1, them N states for k=2),  columns are iteration
Gs=SharedArray{Float64}(N*(horzLen+1),maxIt) #row are time (N states for k=1, them N states for k=2),  columns are iteration
Gi=zeros((horzLen+1),maxIt) #row are time (N states for k=1, them N states for k=2),  columns are iteration

#Jacobian C Vectors
Cs=SharedArray{Float64}(N*(horzLen+1),maxIt)  #row are time,  columns are iteration
Cu=SharedArray{Float64}(N*(horzLen+1),maxIt)  #row are time,  columns are iteration
Ci=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration
Ct=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration

uSum=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration
currConst=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration
deltaY=zeros(1,maxIt)

for p=1:maxIt-1
    @printf "Starting iteration %g \n" p

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
            Cu[ind,p+1]=1cValMax-1cValMin
            # cVal=zeros(length(ind),1)
            # cVal[kappaMax.>0]=1
            # cVal[kappaMin.<0]=-1
            # Cu[ind,p+1]=cVal

            cValMax=abs.(snVal-1).<tolS
            cValMin=abs.(snVal-target[ind]).<tolS
            Cs[ind,p+1]=1cValMax-1cValMin
            # cVal=zeros(length(ind),1)
            # cVal[socMax.>0]=1
            # cVal[socMin.<0]=-1
            # Cs[ind,p+1]=cVal

            Sn[ind,p+1]=snVal
    		Un[ind,p+1]=uVal
            Gu[ind,p+1]=2*Ri[evInd,1]*uVal
            Gs[ind,p+1]=2*Qsi[evInd,1]*snVal-2*Qsi[evInd,1]
            #Gu[collect(evInd:N:length(Gu[:,p+1])),p+1]=sigmaU[evInd,1]*(evVu-uVal)+lambda
            #Gs[collect(evInd:N:length(Gs[:,p+1])),p+1]=sigmaN[evInd,1]*(evVs-snVal)-lambda
        end
    end

    #N+1 decoupled problem aka transformer current
    tM = Model(solver = IpoptSolver())
    @variable(tM,itotal[1:(horzLen+1)])
    @variable(tM,xt[1:(horzLen+1)])
    tMobj=sum(-Lam[k,p]*itotal[k]+
              rhoALADp[1,p]/2*sigmaI*(itotal[k]-Vi[k,p])^2+
              rhoALADp[1,p]/2*sigmaT*(xt[k]-Vt[k,p])^2  for k=1:(horzLen+1))
    @objective(tM,Min, tMobj)
    @NLconstraint(tM,tempCon1,xt[1]-tauP*xt0-gammaP*(itotal[1])^2-rhoP*w[stepI*2,1]==0)
    @NLconstraint(tM,tempCon2[k=1:horzLen],xt[k+1]-tauP*xt[k]-gammaP*(itotal[k+1])^2-rhoP*w[stepI*2+k*2,1]==0)
    if noTlimit==0
    	@constraint(tM,upperTCon,xt.<=Tmax)
    end
    @constraint(tM,lowerTCon,xt.>=0)
    @constraint(tM,KappaMin,itotal.>=0)
    @constraint(tM,KappaMax,itotal.<=ItotalMax)
    TT = STDOUT # save original STDOUT stream
    redirect_stdout()
    statusTM = solve(tM)
    redirect_stdout(TT)
    if statusTM!=:Optimal
        println("solver issues with XFRM NLP")
        return
    else
        kappaMax=-getdual(KappaMin)
        kappaMin=-getdual(KappaMax)
        tMax=-getdual(upperTCon)
        tMin=-getdual(lowerTCon)
        lambdaTemp=[-getdual(tempCon1);-getdual(tempCon2)]
        iVal=getvalue(itotal)
        xtVal=getvalue(xt)

        cValMax=abs.(iVal-ItotalMax).<tolI
        cValMin=abs.(iVal-0).<tolI
        Ci[:,p+1]=1cValMax-1cValMin
        # cVal=zeros(length(ind),1)
        # cVal[kappaMax.>0]=1
        # cVal[kappaMin.<0]=-1
        # Ci[:,p+1]=cVal


        cValMax=abs.(xtVal-Tmax).<tolT
        cValMin=abs.(xtVal-0).<tolT
        Ct[:,p+1]=1cValMax-1cValMin
        # cVal=zeros(length(ind),1)
        # cVal[tMax.>0]=1
        # cVal[tMin.<0]=-1
        # Ct[:,p+1]=cVal

        Xt[:,p+1]=xtVal
        Itotal[:,p+1]=iVal
        Gi[:,p+1]=0
        #Gz[:,p+1]=sigmaZ*(Vz[:,p]-zVal)-repeat(-Lam[:,p],inner=S)
    end

    for k=1:horzLen+1
        uSum[k,p+1]=sum(Un[(k-1)*N+n,p+1] for n=1:N)
        currConst[k,p+1]=uSum[k,p+1] + w[(k-1)*2+(stepI*2-1),1] - Itotal[k,p+1]
    end

    #check for convergence
    constGap=norm(currConst[:,p+1],1)
    cc=norm(vcat((Vu[:,p]-Un[:,p+1]),(Vi[:,p]-Itotal[:,p+1])),1)
    #convCheck=rhoALAD*norm(vcat(repmat(sigmaU,horzLen+1,1).*(Vu[:,p]-Un[:,p+1]),sigmaZ*(Vz[:,p]-Z[:,p+1])),1)
    objFun(sn,xt,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*Qsi[n,1]     for n=1:N) for k=1:horzLen+1) +
                    sum((xt[k,1]-1)^2*Qsi[N+1,1]                 for k=1:horzLen+1) +
                    sum(sum((u[(k-1)*N+n,1])^2*Ri[n,1]           for n=1:N) for k=1:horzLen+1)
    fGap= abs(objFun(Sn[:,p+1],Xt[:,p+1],Un[:,p+1])-fStarNL)
    snGap=norm((Sn[:,p+1]-snStarNL),2)
    itGap = norm(Lam[:,p]-Lam[:,max(p-1,1)],2)
    convGap = norm(Lam[:,p]-lamCurrStarNL,2)
    fConvALAD[p,1]=fGap
    snConvALAD[p,1]=snGap
    itConvALAD[p,1]=itGap
    constConvALAD[p,1]=constGap
    ConvALAD[p,1]=convGap
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
    objExp=sum(sum(0.5*dUn[(k-1)*N+i,1]^2*2*Ri[i,1]+Gu[(k-1)*N+i,p+1]*dUn[(k-1)*N+i,1]+
                   0.5*dSn[(k-1)*N+i,1]^2*2*Qsi[i,1]+Gs[(k-1)*N+i,p+1]*dSn[(k-1)*N+i,1] for i=1:N) for k=1:(horzLen+1))
    objExp=objExp+sum(0.5*dI[k,1]^2*(Hi-lambdaTemp[k,1]*2*gammaP)+
                      0.5*dXt[k,1]^2*Ht   for k=1:(horzLen+1))
    #objExp=objExp+Lam[:,p]'*relaxS+muALAD/2*sum(relaxS[k,1]^2 for k=1:horzLen+1)
	@objective(cM,Min, objExp)
    @constraint(cM,currCon[k=1:horzLen+1],sum(Un[(k-1)*(N)+n,p+1]+dUn[(k-1)*(N)+n,1] for n=1:N)-(Itotal[k,p+1]+dI[k])==-w[(k-1)*2+1])#+relaxS[k,1])
    @constraint(cM,stateCon1[n=1:N],dSn[n,1]==etaP[n,1]*dUn[n,1])
    @constraint(cM,stateCon2[k=1:horzLen,n=1:N],dSn[n+(k)*(N),1]==dSn[n+(k-1)*(N),1]+etaP[n,1]*dUn[n+(k)*(N),1])
    @constraint(cM,tempCon1,dXt[1,1]==2*gammaP*Itotal[1,p+1]*dI[1])
    @constraint(cM,tempCon2[k=1:horzLen],dXt[k+1,1]==tauP*dXt[k,1]+2*gammaP*Itotal[k+1,p+1]*dI[k+1,1])
    @constraint(cM,Ci[:,p+1].*dI.<=0)
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
    Vi[:,p+1]=Vi[:,p]+alpha1*(Itotal[:,p+1]-Vi[:,p])+alpha2*getvalue(dI)
    Vs[:,p+1]=Vs[:,p]+alpha1*(Sn[:,p+1]-Vs[:,p])+alpha2*getvalue(dSn)
    Vt[:,p+1]=Vt[:,p]+alpha1*(Xt[:,p+1]-Vt[:,p])+alpha2*getvalue(dXt)

    rhoALADp[1,p+1]=rhoALADp[1,p]*rhoRate #increase rho every iteration

    deltaY[1,p+1]=norm(vcat(getvalue(dUn),getvalue(dI),getvalue(dSn),getvalue(dXt)),Inf)
end

println("plotting....")
xPlot=zeros(horzLen+1,N)
uPlot=zeros(horzLen+1,N)
for ii= 1:N
	xPlot[:,ii]=Sn[collect(ii:N:length(Sn[:,convIt])),convIt]
    uPlot[:,ii]=Un[collect(ii:N:length(Un[:,convIt])),convIt]
end

#plot(x=1:horzLen+1,y=xPlot2[:,ii])
# plot(x=1:Kn[ii,1],y=xPlot2[1:Kn[ii,1],ii])

pd1aladNL=plot(xPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV SOC"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_decentral_ALADIN_SOC.png", 24inch, 12inch), pd1alad) end

pd2aladNL=plot(uPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV Current (A)"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_decentral_ALADIN_Curr.png", 24inch, 12inch), pd2alad) end

pd3aladNL=plot(layer(x=1:horzLen+1,y=Xt[:,convIt],Geom.line,Theme(default_color=colorant"blue")),
		#layer(x=1:horzLen+1,y=Tactual,Geom.line,Theme(default_color=colorant"green")),
		yintercept=[Tmax],Geom.hline(color=["red"],style=:dot),
		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",key_position = :top,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
		#Guide.manual_color_key("", ["PWL Temp", "Actual Temp"], ["blue", "green"]))
if drawFig==1 draw(PNG(path*"J_decentral_ALADIN_Temp.png", 24inch, 12inch), pd3alad) end

pd4aladNL=plot(layer(x=1:horzLen+1,y=Lam[:,convIt],Geom.line),
		layer(x=1:horzLen+1,y=lamCurrStarNL,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
		Guide.xlabel("Time"), Guide.ylabel(raw"Lambda ($/A)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_decentral_ALADIN_Lam.png", 24inch, 12inch), pd4alad) end


lamPlotaladNL=plot(Lam[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=lamCurrStarNL,Geom.line,Theme(default_color=colorant"black",line_width=4pt)),
			Guide.xlabel("Time"), Guide.ylabel("Lambda"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
uSumPlotaladNL=plot(uSum[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=uSumStarNL,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("U sum"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
iPlotaladNL=plot(Itotal[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=itotalStarNL,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("Z sum"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
xtPlotaladNL=plot(Xt[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=xtStarNL,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("Z sum"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
constPlotaladNL2=plot(currConst,x=Row.index,y=Col.value,color=Col.index,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel("curr constraint diff"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_ALADIN_LamConv.png", 36inch, 12inch), lamPlotalad) end

activeSet=zeros(convIt,1)
setChanges=zeros(convIt,1)
for ii=2:convIt
    activeSet[ii,1]=sum(abs.(Cs[:,ii]))+sum(abs.(Ct[:,ii]))+
              sum(abs.(Cu[:,ii]))+sum(abs.(Ci[:,ii]))
    setChanges[ii,1]=sum(abs.(Cs[:,ii]-Cs[:,ii-1]))+sum(abs.(Ct[:,ii]-Ct[:,ii-1]))+
                     sum(abs.(Cu[:,ii]-Cu[:,ii-1]))+sum(abs.(Ci[:,ii]-Ci[:,ii-1]))
end
activeSetPlot=plot(x=2:convIt,y=activeSet[2:convIt],Geom.line,
                   Guide.xlabel("Iteration"), Guide.ylabel("Total Active inequality constraints",orientation=:vertical),
                   Coord.Cartesian(xmin=2,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
       			   minor_label_font_size=26pt,key_label_font_size=26pt))
setChangesPlot=plot(x=3:convIt,y=setChanges[3:convIt],Geom.line,
                    Guide.xlabel("Iteration"), Guide.ylabel("Changes in Active inequality constraints",orientation=:vertical),
                    Coord.Cartesian(xmin=3,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
        			minor_label_font_size=26pt,key_label_font_size=26pt))
solChangesplot=plot(layer(x=2:convIt,y=deltaY[2:convIt],Geom.line),
                    layer(x=2:convIt,y=convCheck[2:convIt],Geom.line),Scale.y_log)

convItPlotaladNL=plot(x=1:convIt,y=itConvALAD[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
convPlotaladNL=plot(x=1:convIt,y=ConvALAD[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("central lambda gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
constPlotaladNL=plot(x=1:convIt,y=constConvALAD[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("consensus gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
fPlotaladNL=plot(x=1:convIt-1,y=fConvALAD[1:convIt-1,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("obj function gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_ALADIN_Conv.png", 36inch, 12inch), vstack(convItPlotalad,convPlotalad,fPlotalad)) end
