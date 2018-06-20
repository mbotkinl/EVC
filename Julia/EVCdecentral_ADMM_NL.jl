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

#initialize with current states
sn0=evS.s0
xt0=evS.t0

stepI = 1
convChk = 1e-8
maxIt=50
convIt=maxIt

#admm  initial parameters and guesses
#ρADMM=10.0^(0)
ρADMM=5e5



#u w and z are one index ahead of x. i.e the x[k+1]=x[k]+eta*u[k+1]
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCsetup.jl")
dCMadmm=convMetrics()
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

for p in 1:maxIt-1
	#ρADMMp = ρADMM*ceil(p/2)
    ρADMMp = ρADMM

    #x minimization eq 7.66 in Bertsekas
    @sync @parallel for evInd=1:N
        evV=dLogadmm.Vu[collect(evInd:N:length(dLogadmm.Vu[:,p])),p]
        target=zeros((horzLen+1),1)
		target[(evS.Kn[evInd,1]-(stepI-1)):1:length(target),1]=evS.Snmin[evInd,1]
    	evM = Model(solver = GurobiSolver())
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
    @NLconstraint(tM,tempCon1,xt[1,1]==evS.τP*xt0+evS.γP*(itotal[1,1])^2+evS.ρP*evS.w[stepI*2,1])
    @NLconstraint(tM,tempCon2[k=1:horzLen],xt[k+1,1]==evS.τP*xt[k,1]+evS.γP*(itotal[k+1,1])^2+evS.ρP*evS.w[stepI*2+k*2,1])
    if noTlimit==0
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
		#Lam[k,p+1]=max.(Lam[k,p]+ρADMMp/(S*(N))*(currConst[k,1]),0)
		dLogadmm.Lam[k,p+1]=dLogadmm.Lam[k,p]+ρADMMp/(S*(N))*(dLogadmm.couplConst[k,p+1])
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
	fGap= objFun(dLogadmm.Sn[:,p+1],dLogadmm.Xt[:,p+1],dLogadmm.Un[:,p+1])-fStarNL
	#fGap= objFun(Sn[:,p],Xt[:,p],Un[:,p])-fStar
	snGap=norm((dLogadmm.Sn[:,p+1]-snStarNL),2)
	constGap=norm(dLogadmm.couplConst[:,p+1],2)
	itGap = norm(dLogadmm.Lam[:,p+1]-dLogadmm.Lam[:,p],2)
	convGap = norm(dLogadmm.Lam[:,p+1]-lamCurrStarNL,2)
	dCMadmm.objVal[p,1]=abs(fGap)
	dCMadmm.sn[p,1]=snGap
	dCMadmm.couplConst[p,1]=constGap
	dCMadmm.lamIt[p,1]=itGap
	dCMadmm.lam[p,1]=convGap
	if(itGap <= convChk )
		@printf "Converged after %g iterations\n" p
		convIt=p+1
		break
	else
		@printf "lastGap %e after %g iterations\n" itGap p
		@printf "convGap %e after %g iterations\n" convGap p
		@printf "constGap %e after %g iterations\n" constGap p
        @printf "snGap %e after %g iterations\n" snGap p
		@printf("fGap %e after %g iterations\n\n",fGap,p)

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

pd1NLadmm=plot(xPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV SOC"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_centralNL_ADMM_SOC.png", 24inch, 12inch), pd1NLadmm) end

pd2NLadmm=plot(uPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV Current (A)"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_centralNL_ADMM_Curr.png", 24inch, 12inch), pd2NLadmm) end

pd3NLadmm=plot(layer(x=1:horzLen+1,y=dLogadmm.Xt[:,convIt],Geom.line,Theme(default_color=colorant"blue")),
		#layer(x=1:horzLen+1,y=Tactual,Geom.line,Theme(default_color=colorant"green")),
		yintercept=[evS.Tmax],Geom.hline(color=["red"],style=:dot),
		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",key_position = :top,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
		#Guide.manual_color_key("", ["PWL Temp", "Actual Temp"], ["blue", "green"]))
if drawFig==1 draw(PNG(path*"J_centralNL_ADMM_Temp.png", 24inch, 12inch), pd3NLadmm) end

pd4NLadmm=plot(x=1:horzLen+1,y=dLogadmm.Lam[:,convIt],Geom.line,
		Guide.xlabel("Time"), Guide.ylabel(raw"Lambda ($/A)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_centralNL_ADMM_Lam.png", 24inch, 12inch), pd4NLadmm) end

fName="J_Central.png"


lamPlotNLadmm=plot(dLogadmm.Lam[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
            layer(x=1:horzLen+1,y=lamCurrStarNL,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("Lambda",orientation=:vertical),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
constPlotNLadmm2=plot(dLogadmm.couplConst[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel("curr constraint diff",orientation=:vertical),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_ADMM_LamConv.png", 36inch, 12inch), lamPlotadmm) end

convItPlotNLadmm=plot(x=1:convIt,y=dCMadmm.lamIt[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap",orientation=:vertical),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
convPlotNLadmm=plot(x=1:convIt,y=dCMadmm.lam[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap",orientation=:vertical),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
fPlotadmm=plot(x=1:convIt-1,y=dCMadmm.objVal[1:convIt-1,1],Geom.line,#Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap",orientation=:vertical),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
constPlotadmm=plot(x=1:convIt,y=dCMadmm.couplConst[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap",orientation=:vertical),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_ADMM_Conv.png", 36inch, 12inch), convPlotadmm) end
