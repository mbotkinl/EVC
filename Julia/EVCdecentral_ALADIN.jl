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

#initialize
sn0=evS.s0
xt0=evS.t0

stepI = 1
epsilon = 1e-8
tolU=1e-4
tolS=1e-8
tolT=1e-4
tolZ=1e-6
maxIt=50

convIt=maxIt

#ALADIN tuning and initial guess
σU=1*ones(N,1)
σS=ones(N,1)/10 #for kA
#σS=ones(N,1)*100 #for A
σZ=1/N
σT=1/10000 #for kA
#σT=1/10  #for A

# Hu=2*evS.Ri
# Hs=2*evS.Qsi
Hu=2*evS.Ri *((1.5-2.5)*rand()+2.5)
Hs=2*evS.Qsi *((1.5-2.5)*rand()+2.5)
Hz=1e-6
Ht=1e-6

ρALAD=1
ρRate=1.15
muALAD=10^8


include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCsetup.jl")
dCMalad=convMetrics()
dLogalad=itLogPWL()
convCheck=zeros(maxIt,1)

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

deltaY=zeros(1,maxIt)
ρALADp=ρALAD*ones(1,maxIt)

for p=1:maxIt-1

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
    fGap= abs(objFun(dLogalad.Sn[:,p+1],dLogalad.Xt[:,p+1],dLogalad.Un[:,p+1])-fStar)
    snGap=norm((dLogalad.Sn[:,p+1]-snStar),2)
    unGap=norm((dLogalad.Un[:,p+1]-uStar),2)
    itGap = norm(dLogalad.Lam[:,p]-dLogalad.Lam[:,max(p-1,1)],2)
    convGap = norm(dLogalad.Lam[:,p]-lamCurrStar,2)
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
    deltaY[1,p+1]=norm(vcat(getvalue(dUn),getvalue(dZ),getvalue(dSn),getvalue(dXt)),Inf)
end

println("plotting....")
xPlot=zeros(horzLen+1,N)
uPlot=zeros(horzLen+1,N)
for ii= 1:N
	xPlot[:,ii]=dLogalad.Sn[collect(ii:N:length(dLogalad.Sn[:,convIt])),convIt]
    uPlot[:,ii]=dLogalad.Un[collect(ii:N:length(dLogalad.Un[:,convIt])),convIt]
end

pd1alad=plot(xPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV SOC"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_decentral_ALADIN_SOC.png", 24inch, 12inch), pd1alad) end

pd2alad=plot(uPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV Current (kA)"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_decentral_ALADIN_Curr.png", 24inch, 12inch), pd2alad) end

pd3alad=plot(layer(x=1:horzLen+1,y=dLogalad.Xt[:,convIt],Geom.line,Theme(default_color=colorant"blue")),
		#layer(x=1:horzLen+1,y=Tactual,Geom.line,Theme(default_color=colorant"green")),
		yintercept=[evS.Tmax],Geom.hline(color=["red"],style=:dot),
		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",key_position = :top,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
		#Guide.manual_color_key("", ["PWL Temp", "Actual Temp"], ["blue", "green"]))
if drawFig==1 draw(PNG(path*"J_decentral_ALADIN_Temp.png", 24inch, 12inch), pd3alad) end

pd4alad=plot(layer(x=1:horzLen+1,y=dLogalad.Lam[:,convIt],Geom.line),
		layer(x=1:horzLen+1,y=lamCurrStar,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
		Guide.xlabel("Time"), Guide.ylabel(raw"Lambda ($/kA)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_decentral_ALADIN_Lam.png", 24inch, 12inch), pd4alad) end


lamPlotalad=plot(dLogalad.Lam[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=lamCurrStar,Geom.line,Theme(default_color=colorant"black",line_width=4pt)),
			Guide.xlabel("Time"), Guide.ylabel("Lambda"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
uSumPlotalad=plot(dLogalad.uSum[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=uSumStar,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("U sum"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
zSumPlotalad=plot(dLogalad.zSum[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=zSumStar,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("Z sum"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
constPlotalad2=plot(dLogalad.couplConst[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel("curr constraint diff"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_ALADIN_LamConv.png", 36inch, 12inch), lamPlotalad) end

activeSet=zeros(convIt,1)
setChanges=zeros(convIt,1)
for ii=2:convIt
    activeSet[ii,1]=sum(abs.(dLogalad.Cs[:,ii]))+sum(abs.(dLogalad.Ct[:,ii]))+
              sum(abs.(dLogalad.Cu[:,ii]))+sum(abs.(dLogalad.Cz[:,ii]))
    setChanges[ii,1]=sum(abs.(dLogalad.Cs[:,ii]-dLogalad.Cs[:,ii-1]))+sum(abs.(dLogalad.Ct[:,ii]-dLogalad.Ct[:,ii-1]))+
                     sum(abs.(dLogalad.Cu[:,ii]-dLogalad.Cu[:,ii-1]))+sum(abs.(dLogalad.Cz[:,ii]-dLogalad.Cz[:,ii-1]))
end
activeSetPlot=plot(x=2:convIt,y=activeSet[2:convIt],Geom.line,
                   Guide.xlabel("Iteration"), Guide.ylabel("Total Active inequality constraints",orientation=:vertical),
                   Coord.Cartesian(xmin=2,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
       			   minor_label_font_size=26pt,key_label_font_size=26pt))
setChangesPlot=plot(x=3:convIt,y=setChanges[3:convIt],Geom.line,
                    Guide.xlabel("Iteration"), Guide.ylabel("Changes in Active inequality constraints",orientation=:vertical),
                    Coord.Cartesian(xmin=3,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
        			minor_label_font_size=26pt,key_label_font_size=26pt))
solChangesplot=plot(layer(x=2:convIt,y=deltaY[2:convIt],Geom.line,Theme(default_color=colorant"green")),
                    layer(x=2:convIt,y=convCheck[2:convIt],Geom.line,Theme(default_color=colorant"red")),
                    Scale.y_log,Guide.manual_color_key("", ["ΔY","y-x"], ["green","red"]))

convItPlotalad=plot(x=1:convIt,y=dCMalad.lamIt[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
convPlotalad=plot(x=1:convIt,y=dCMalad.lam[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("central lambda gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
constPlotalad=plot(x=1:convIt,y=dCMalad.couplConst[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("const gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
fPlotalad=plot(x=1:convIt-1,y=dCMalad.objVal[1:convIt-1,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("obj function gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_ALADIN_Conv.png", 36inch, 12inch), vstack(convItPlotalad,convPlotalad,fPlotalad)) end
