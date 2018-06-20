#EVC with ALADIN for PWL convex relaxation
#Micah Botkin-Levy
#5/22/18

#formulation try 2
#u, sn, xt, and z are all in "y", v is the auxilliary variable corresponding to x in literature
#current constraint is coupling

tic()
#pull out a few key variables
N=evS.N
S=evS.S

#initialize with current states
sn0=evS.s0
xt0=evS.t0

stepI = 1
epsilon = 1e-12
tolU=1e-4
tolS=1e-8
tolT=1e-4
tolI=1e-6
maxIt=30
convIt=maxIt

#ALADIN tuning and initial guess
##σs are tuned to i_n [.010kA]
σU=1*ones(N,1)
σS=ones(N,1)/10
σI=1/N
σT=1/10000
Hu=2*evS.Ri*(1+rand())
Hs=2*evS.Qsi*(1+rand())
Hi=1e-6
Ht=1e-6
ρALAD=1
ρRate=1.15
muALAD=10^8

#.1/.1 looks good expect for const
#.01/2 increasing ρ helps const
#.1/1.1 best so far???

include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCsetup.jl")
dCMalad=convMetrics()
dLogalad=itLogNL()
convCheck=zeros(maxIt,1)

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
dLogalad.Vu[:,1]=vu0 #initial guess goes in column 1
dLogalad.Vs[:,1]=vs0 #initial guess goes in column 1
dLogalad.Vi[:,1]=vi0
dLogalad.Vt[:,1]=vt0
dLogalad.Lam[:,1]=lambda0

ρALADp=ρALAD*ones(1,maxIt)
ΔY=zeros(1,maxIt)

for p=1:maxIt-1
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
        evM = Model(solver = GurobiSolver(NumericFocus=3))
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
    @NLconstraint(tM,tempCon1,xt[1]-evS.τP*xt0-evS.γP*(itotal[1])^2-evS.ρP*evS.w[stepI*2,1]==0)
    @NLconstraint(tM,tempCon2[k=1:horzLen],xt[k+1]-evS.τP*xt[k]-evS.γP*(itotal[k+1])^2-evS.ρP*evS.w[stepI*2+k*2,1]==0)
    if noTlimit==0
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
    fGap= abs(objFun(dLogalad.Sn[:,p+1],dLogalad.Xt[:,p+1],dLogalad.Un[:,p+1])-fStarNL)
    fGap2= abs((objFun(dLogalad.Sn[:,p+1],dLogalad.Xt[:,p+1],dLogalad.Un[:,p+1])-fStarNL)/fStarNL)
    snGap=norm((dLogalad.Sn[:,p+1]-snStarNL),2)
    itGap = norm(dLogalad.Lam[:,p]-dLogalad.Lam[:,max(p-1,1)],2)
    convGap = norm(dLogalad.Lam[:,p]-lamCurrStarNL,2)
    convGap2 = norm(round.(dLogalad.Lam[:,p]-lamCurrStarNL,10)./lamCurrStarNL,2) #need some rounding here if both are 1e-8

    dCMalad.objVal[p,1]=fGap
    dCMalad.sn[p,1]=snGap
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

println("plotting....")
xPlot=zeros(horzLen+1,N)
uPlot=zeros(horzLen+1,N)
for ii= 1:N
	xPlot[:,ii]=dLogalad.Sn[collect(ii:N:length(dLogalad.Sn[:,convIt])),convIt]
    uPlot[:,ii]=dLogalad.Un[collect(ii:N:length(dLogalad.Un[:,convIt])),convIt]
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
		Guide.xlabel("Time"), Guide.ylabel("PEV Current (kA)"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_decentral_ALADIN_Curr.png", 24inch, 12inch), pd2alad) end

pd3aladNL=plot(layer(x=1:horzLen+1,y=dLogalad.Xt[:,convIt],Geom.line,Theme(default_color=colorant"blue")),
		#layer(x=1:horzLen+1,y=Tactual,Geom.line,Theme(default_color=colorant"green")),
		yintercept=[evS.Tmax],Geom.hline(color=["red"],style=:dot),
		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",key_position = :top,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
		#Guide.manual_color_key("", ["PWL Temp", "Actual Temp"], ["blue", "green"]))
if drawFig==1 draw(PNG(path*"J_decentral_ALADIN_Temp.png", 24inch, 12inch), pd3alad) end

pd4aladNL=plot(layer(x=1:horzLen+1,y=dLogalad.Lam[:,convIt],Geom.line),
		layer(x=1:horzLen+1,y=lamCurrStarNL,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
		Guide.xlabel("Time"), Guide.ylabel(raw"Lambda ($/kA)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_decentral_ALADIN_Lam.png", 24inch, 12inch), pd4alad) end


lamPlotaladNL=plot(dLogalad.Lam[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=lamCurrStarNL,Geom.line,Theme(default_color=colorant"black",line_width=4pt)),
			Guide.xlabel("Time"), Guide.ylabel("Lambda"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
uSumPlotaladNL=plot(dLogalad.uSum[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=uSumStarNL,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("U sum"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
iPlotaladNL=plot(dLogalad.Itotal[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=itotalStarNL,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("Z sum"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
xtPlotaladNL=plot(dLogalad.Xt[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=xtStarNL,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("Z sum"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
constPlotaladNL2=plot(dLogalad.couplConst,x=Row.index,y=Col.value,color=Col.index,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel("curr constraint diff"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_ALADIN_LamConv.png", 36inch, 12inch), lamPlotalad) end

activeSet=zeros(convIt,1)
setChanges=zeros(convIt,1)
for ii=2:convIt
    activeSet[ii,1]=sum(abs.(dLogalad.Cs[:,ii]))+sum(abs.(dLogalad.Ct[:,ii]))+
              sum(abs.(dLogalad.Cu[:,ii]))+sum(abs.(dLogalad.Ci[:,ii]))
    setChanges[ii,1]=sum(abs.(dLogalad.Cs[:,ii]-dLogalad.Cs[:,ii-1]))+sum(abs.(dLogalad.Ct[:,ii]-dLogalad.Ct[:,ii-1]))+
                     sum(abs.(dLogalad.Cu[:,ii]-dLogalad.Cu[:,ii-1]))+sum(abs.(dLogalad.Ci[:,ii]-dLogalad.Ci[:,ii-1]))
end
activeSetPlot=plot(x=2:convIt,y=activeSet[2:convIt],Geom.line,
                   Guide.xlabel("Iteration"), Guide.ylabel("Total Active inequality constraints",orientation=:vertical),
                   Coord.Cartesian(xmin=2,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
       			   minor_label_font_size=26pt,key_label_font_size=26pt))
setChangesPlot=plot(x=3:convIt,y=setChanges[3:convIt],Geom.line,
                    Guide.xlabel("Iteration"), Guide.ylabel("Changes in Active inequality constraints",orientation=:vertical),
                    Coord.Cartesian(xmin=3,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
        			minor_label_font_size=26pt,key_label_font_size=26pt))
solChangesplot=plot(layer(x=2:convIt,y=ΔY[2:convIt],Geom.line),
                    layer(x=2:convIt,y=convCheck[2:convIt],Geom.line),Scale.y_log)

convItPlotaladNL=plot(x=1:convIt,y=dCMalad.lamIt[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
convPlotaladNL=plot(x=1:convIt,y=dCMalad.lam[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("central lambda gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
constPlotaladNL=plot(x=1:convIt,y=dCMalad.couplConst[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("consensus gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
fPlotaladNL=plot(x=1:convIt-1,y=dCMalad.objVal[1:convIt-1,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("obj function gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_ALADIN_Conv.png", 36inch, 12inch), vstack(convItPlotalad,convPlotalad,fPlotalad)) end
