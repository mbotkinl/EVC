#EVC with ADMM for PWL convex relaxation
#Micah Botkin-Levy
#5/22/18

#formulation try 2
#u, sn, xt, and z are all in "x" v is the auxilliary variable corresponding to z in literature
#current constraint is coupling

tic()

#pull out a few key variables
N=evS.N
S=evS.S

#initialize
sn0=evS.s0
xt0=evS.t0

stepI = 1;
convChk = 1e-16
maxIt=50
convIt=maxIt

#admm  initial parameters and guesses
#ρADMM=10.0^(0)
ρADMM=10^6 #for kA
#ρADMM=1    #for A


#u w and z are one index ahead of x. i.e the x[k+1]=x[k]+η*u[k+1]
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCsetup.jl")
dCMadmm=convMetrics()
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

for p in 1:maxIt-1
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
	fGap= abs(objFun(dLogadmm.Sn[:,p+1],dLogadmm.Xt[:,p+1],dLogadmm.Un[:,p+1])-fStar)
	#fGap= objFun(Sn[:,p],Xt[:,p],Un[:,p])-fStar
	snGap=norm((dLogadmm.Sn[:,p+1]-snStar),2)
	unGap=norm((dLogadmm.Un[:,p+1]-uStar),2)
	constGap=norm(dLogadmm.couplConst[:,p+1],2)
	itGap = norm(dLogadmm.Lam[:,p+1]-dLogadmm.Lam[:,p],2)
	convGap = norm(dLogadmm.Lam[:,p+1]-lamCurrStar,2)
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
toc()


println("plotting....")
xPlot=zeros(horzLen+1,N)
uPlot=zeros(horzLen+1,N)
for ii= 1:N
	xPlot[:,ii]=dLogadmm.Sn[collect(ii:N:length(dLogadmm.Sn[:,convIt])),convIt]
	uPlot[:,ii]=dLogadmm.Un[collect(ii:N:length(dLogadmm.Un[:,convIt])),convIt]
end


pd1admm=plot(xPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV SOC"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_decentral_ADMM_SOC.png", 24inch, 12inch), pd1admm) end

pd2admm=plot(uPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV Current (A)"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_decentral_ADMM_Curr.png", 24inch, 12inch), pd2admm) end

pd3admm=plot(layer(x=1:horzLen+1,y=dLogadmm.Xt[:,convIt],Geom.line,Theme(default_color=colorant"blue")),
		layer(x=1:horzLen+1,y=dLogadmm.Tactual[:,convIt],Geom.line,Theme(default_color=colorant"green")),
		yintercept=[evS.Tmax],Geom.hline(color=["red"],style=:dot),
		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",key_position = :top,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt),
		Guide.manual_color_key("", ["PWL Temp", "Actual Temp"], ["blue", "green"]))
if drawFig==1 draw(PNG(path*"J_decentral_ADMM_Temp.png", 24inch, 12inch), pd3admm) end

pd4admm=plot(x=1:horzLen+1,y=dLogadmm.Lam[:,convIt],Geom.line,
		layer(x=1:horzLen+1,y=lamCurrStar,Geom.line,Theme(default_color=colorant"black")),
		Guide.xlabel("Time"), Guide.ylabel(raw"Lambda ($/A)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_decentral_ADMM_Lam.png", 24inch, 12inch), pd4admm) end


lamPlotadmm=plot(dLogadmm.Lam[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=lamCurrStar,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("Lambda"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
constPlotadmm2=plot(dLogadmm.couplConst[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel("curr constraint diff"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
zSumPlotadmm=plot(dLogadmm.zSum[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=zSumStar,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("Z sum"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
uSumPlotadmm=plot(dLogadmm.uSum[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=uSumStar,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("U sum"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_ADMM_LamConv.png", 36inch, 12inch), lamPlotadmm) end

convItPlotadmm=plot(x=1:convIt,y=dCMadmm.lamIt[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
convPlotadmm=plot(x=1:convIt,y=dCMadmm.lam[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("central lambda gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
fPlotadmm=plot(x=1:convIt-1,y=dCMadmm.objVal[1:convIt-1,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("obj function gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
constPlotadmm=plot(x=1:convIt,y=dCMadmm.couplConst[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("curr constraint Gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_ADMM_Conv.png", 36inch, 12inch), vstack(convItPlotadmm,convPlotadmm,fPlotadmm,constPlotadmm)) end
