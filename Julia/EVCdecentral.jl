#Micah Botkin-Levy
#4/10/18

tic()
#pull out a few key variables
N=evS.N
S=evS.S

#initialize
sn0=evS.s0
xt0=evS.t0

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
maxIt=50
convIt=maxIt

#u w and z are one index ahead of x. i.e the x[k+1]=x[k]+eta*u[k+1]
alphaP=alpha*ones(maxIt,1)

include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//structEVC.jl")
dCM=convMetrics()
dLog=itLogPWL()

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
	fGap= objFun(dLog.Sn[:,p+1],dLog.Un[:,p+1])-fStar
	snGap=norm((dLog.Sn[:,p+1]-snStar),2)
	unGap=norm((dLog.Un[:,p+1]-uStar),2)
	itGap = norm(dLog.Lam[:,p+1]-dLog.Lam[:,p],2)
	if updateMethod=="fastAscent"
		convGap = norm(dLog.Lam[:,p+1]-lamTempStar,2)
	else
		convGap = norm(dLog.Lam[:,p+1]-lamCurrStar,2)
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

s=Symbol(@sprintf("dCM_%s",updateMethod))
v=Symbol(@sprintf("dCM"))
@eval(($s)=($v))

println("plotting....")
xPlot=zeros(horzLen+1,N)
uPlotd=zeros(horzLen+1,N)
for ii= 1:N
	xPlot[:,ii]=dLog.Sn[collect(ii:N:length(dLog.Sn[:,convIt])),convIt]
	uPlotd[:,ii]=dLog.Un[collect(ii:N:length(dLog.Un[:,convIt])),convIt]
end

pd1=plot(xPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV SOC"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=24pt,line_width=3pt,
		minor_label_font_size=20pt,key_label_font_size=20pt))
if drawFig==1 draw(PNG(path*"J_"*updateMethod*"_SOC.png", 24inch, 12inch), pd1) end

pd2=plot(uPlotd,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV Current (kA)"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=24pt,line_width=3pt,
		minor_label_font_size=20pt,key_label_font_size=20pt))
if drawFig==1 draw(PNG(path*"J_"*updateMethod*"_Curr.png", 24inch, 12inch), pd2) end

pd3=plot(layer(x=1:horzLen+1,y=dLog.Tactual[:,convIt],Geom.line,Theme(default_color=colorant"green",line_width=3pt)),
		layer(yintercept=[evS.Tmax],Geom.hline(color=["red"],style=:dot),Theme(line_width=3pt)),
		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=24pt,
		minor_label_font_size=20pt,key_label_font_size=20pt))
if updateMethod=="dualAscent"
	push!(pd3,layer(x=1:horzLen+1,y=dLog.Xt[:,convIt],Geom.line,Theme(default_color=colorant"blue",line_width=3pt)))
	push!(pd3,Guide.manual_color_key("", ["PWL Temp", "Actual Temp"], ["blue", "green"]))
	push!(pd3,Theme(key_position = :top,background_color=colorant"white",major_label_font_size=24pt,line_width=3pt,
	minor_label_font_size=20pt,key_label_font_size=20pt))
end
if drawFig==1 draw(PNG(path*"J_"*updateMethod*"_Temp.png", 24inch, 12inch), pd3) end

if updateMethod=="fastAscent"
	lamLabel=raw"Lambda ($/K)"
else
	lamLabel=raw"Lambda ($/kA)"
end
pd4=plot(layer(x=1:horzLen+1,y=dLog.Lam[:,convIt],Geom.line,Theme(line_width=3pt)),
		Guide.xlabel("Time"), Guide.ylabel(lamLabel),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=24pt,
		minor_label_font_size=20pt,key_label_font_size=20pt))
if drawFig==1 draw(PNG(path*"J_"*updateMethod*"_Lam.png", 24inch, 12inch), pd4) end


#fName="J_Decentral_notfast.png"
#fName="J_Decentral_fast.png"
#draw(PNG(path*fName, 13inch, 14inch), vstack(pd1,pd2,pd3,pd4))


uSumPlotalad=plot(dLog.uSum[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=uSumStar,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("U sum"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
zSumPlotalad=plot(dLog.zSum[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=zSumStar,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("Z sum"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
lamPlot=plot(dLog.Lam[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=lamCurrStar,Geom.line,Theme(default_color=colorant"black",line_width=4pt)),
			Guide.xlabel("Time"), Guide.ylabel("Lambda"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_"*updateMethod*"_LamConv.png", 36inch, 12inch), lamPlot) end

fPlot=plot(x=1:convIt-1,y=dCM.objVal[1:convIt-1,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("obj function gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
convItPlot=plot(x=1:convIt,y=dCM.lamIt[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
convPlot=plot(x=1:convIt,y=dCM.lam[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("Lambda Star Gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_"*updateMethod*"_Conv.png", 36inch, 12inch), convPlot) end



#compare central and decentral current agg
aggU=plot(layer(x=1:horzLen+1,y=sum(uPlot[:,i] for i=1:N),Geom.line,Theme(default_color=colorant"blue")),
		layer(x=1:horzLen+1,y=sum(uPlotd[:,i] for i=1:N),Geom.line,Theme(default_color=colorant"green")),
		Guide.xlabel("Time"), Guide.ylabel("PEV Current (kA)"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white"),
		Guide.manual_color_key("", ["Central", "Decentral"], ["blue", "green"]))
#draw(PNG(path*"aggPlot_fast.png", 13inch, 8inch), aggU)
