#Micah Botkin-Levy
#4/10/18

tic()
#pull out a few key variables
N=evS.N
S=evS.S

#initialize with current states
sn0=evS.s0
xt0=evS.t0

lambda0=1000*ones(horzLen+1,1)
#lambda0=lamCurrStarNL
#lambda0=max.(lamCurrStarNL,0)

if updateMethod=="fastAscent"
	#alpha = 0.1 #for A
	alpha = 5e3 #for kA
	alphaDivRate=2
	minAlpha=1e-6
	#alphaRate=.99
else
	#alpha = .01 #for A
	alpha = 5e3 #for kA
	alphaDivRate=2
	minAlpha=1e-6
	#alphaRate=.99
end

stepI = 1
convChk = 1e-8
maxIt=50
convIt=maxIt

alphaP=alpha*ones(maxIt,1)
ConvDual=zeros(maxIt,1)
itConvDual=zeros(maxIt,1)
constConvDual=zeros(maxIt,1)
fConvDual=zeros(maxIt,1)
snConvDual=zeros(maxIt,1)
unConvDual=zeros(maxIt,1)

#u w and z are one index ahead of x. i.e the x[k+1]=x[k]+eta*u[k+1]
Lam=zeros((horzLen+1),maxIt) #(rows are time, columns are iteration)
Lam[:,1]=lambda0
Xt=zeros((horzLen+1),maxIt) #rows are time
Itotal=zeros((horzLen+1),maxIt)
Sn=SharedArray{Float64}(N*(horzLen+1),maxIt)
Un=SharedArray{Float64}(N*(horzLen+1),maxIt)
uSum=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration

#iterate at each time step until convergence
for p=1:maxIt-1
    #solve subproblem for each EV
	@sync @parallel for evInd=1:N
		ind=[evInd]
		for k=1:horzLen
			append!(ind,k*N+evInd)
		end
        target=zeros((horzLen+1),1)
		target[(evS.Kn[evInd,1]-(stepI-1)):1:length(target),1]=evS.Snmin[evInd,1]
        evM=Model(solver = IpoptSolver())
        @variable(evM,un[1:horzLen+1])
        @variable(evM,sn[1:horzLen+1])
		@objective(evM,Min,sum((sn[k,1]-1)^2*evS.Qsi[evInd,1]+(un[k,1])^2*evS.Ri[evInd,1]+Lam[k,p]*un[k,1] for k=1:horzLen+1))
		@constraint(evM,sn[1,1]==sn0[evInd,1]+evS.ηP[evInd,1]*un[1,1]) #fix for MPC loop
		@constraint(evM,[k=1:horzLen],sn[k+1,1]==sn[k,1]+evS.ηP[evInd,1]*un[k+1,1]) #check K+1
        @constraint(evM,sn.<=1)
        @constraint(evM,sn.>=target)
        @constraint(evM,un.<=evS.imax[evInd,1])
        @constraint(evM,un.>=evS.imin[evInd,1])

		TT = STDOUT # save original STDOUT stream
		redirect_stdout()
        status = solve(evM)
		redirect_stdout(TT)

		if status!=:Optimal
            break
        else
            Sn[ind,p+1]=getvalue(sn) #solved state goes in next time slot
            Un[ind,p+1]=getvalue(un) #current go
        end
    end

	if updateMethod=="dualAscent"
	    #solve coordinator problem
		coorM=Model(solver = IpoptSolver())
	    @variable(coorM,itotal[1:(horzLen+1)])
	    @variable(coorM,xt[1:(horzLen+1)])
	    @objective(coorM,Min,-sum(Lam[k,p]*itotal[k,1] for k=1:(horzLen+1)))
		@NLconstraint(coorM,xt[1,1]==evS.τP*xt0+evS.γP*(itotal[1])^2+evS.ρP*evS.w[2,1]) #fix for MPC loop
		@NLconstraint(coorM,[k=1:horzLen],xt[k+1,1]==evS.τP*xt[k,1]+evS.γP*(itotal[k+1,1])^2+evS.ρP*evS.w[k*2+2,1])
		if noTlimit==0
			@constraint(coorM,upperTCon,xt.<=evS.Tmax)
		end
	    @constraint(coorM,xt.>=0)
	    @constraint(coorM,itotal.<=evS.ItotalMax)
	    @constraint(coorM,itotal.>=0)
		TT = STDOUT # save original STDOUT strea
		redirect_stdout()
	    statusC = solve(coorM)
		redirect_stdout(TT)

	    if statusC!=:Optimal
	        break
		else
			 Xt[:,p+1]=getvalue(xt)
			 Itotal[:,p+1]=getvalue(itotal)
		end

	    #grad of lagragian
		gradL=zeros(horzLen+1,1)
		for k=1:horzLen+1
			uSum[k,p+1]=sum(Un[(k-1)*N+n,p+1] for n=1:N)
			gradL[k,1]=uSum[k,p+1] + evS.w[(k-1)*2+(stepI*2-1),1] - Itotal[k,p+1]
		end
		constConvDual[p,1]=norm(gradL,2)
	end

	if updateMethod=="fastAscent"
		#fast ascent
		if noTlimit==0
			gradL=Tactual[:,p+1]-evS.Tmax*ones(horzLen+1,1)
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

	#Lam[:,p+1]=Lam[:,p]+alpha_p*gradL
    Lam[:,p+1]=max.(Lam[:,p]+alphaP[p+1,1]*gradL,0)

	#check convergence
	objFun(sn,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,1]     for n=1:N) for k=1:horzLen+1) +
					sum(sum((u[(k-1)*N+n,1])^2*evS.Ri[n,1]           for n=1:N) for k=1:horzLen+1)
	fGap= objFun(Sn[:,p+1],Un[:,p+1])-fStarNL
	snGap=norm((Sn[:,p+1]-snStarNL),2)
	unGap=norm((Un[:,p+1]-uStarNL),2)
	itGap = norm(Lam[:,p+1]-Lam[:,p],2)
	if updateMethod=="fastAscent"
		convGap = norm(Lam[:,p+1]-lamTempStarNL,2)
	else
		convGap = norm(Lam[:,p+1]-lamCurrStarNL,2)
	end
	fConvDual[p,1]=abs(fGap)
	snConvDual[p,1]=snGap
	unConvDual[p,1]=unGap
	itConvDual[p,1]=itGap
	ConvDual[p,1]=convGap
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

for name in ["f","sn","un","it",""]
	s=Symbol(@sprintf("%sConvDual_%s",name,updateMethod))
	v=Symbol(@sprintf("%sConvDual",name))
	@eval(($s)=($v))
end

println("plotting....")
xPlot=zeros(horzLen+1,N)
uPlotd=zeros(horzLen+1,N)
for ii= 1:N
	xPlot[:,ii]=Sn[collect(ii:N:length(Sn[:,convIt])),convIt]
	uPlotd[:,ii]=Un[collect(ii:N:length(Un[:,convIt])),convIt]
end

pd1nl=plot(xPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV SOC"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=24pt,line_width=3pt,
		minor_label_font_size=20pt,key_label_font_size=20pt))
if drawFig==1 draw(PNG(path*"J_"*updateMethod*"_SOC.png", 24inch, 12inch), pd1) end

pd2nl=plot(uPlotd,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV Current (kA)"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1,ymin=0),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=24pt,line_width=3pt,
		minor_label_font_size=20pt,key_label_font_size=20pt))
if drawFig==1 draw(PNG(path*"J_"*updateMethod*"_Curr.png", 24inch, 12inch), pd2) end

pd3nl=plot(layer(yintercept=[evS.Tmax],Geom.hline(color=["red"],style=:dot),Theme(line_width=3pt)),
		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=24pt,
		minor_label_font_size=20pt,key_label_font_size=20pt))
if updateMethod=="dualAscent"
	push!(pd3nl,layer(x=1:horzLen+1,y=Xt[:,convIt],Geom.line,Theme(default_color=colorant"blue",line_width=3pt)))
	push!(pd3nl,Guide.manual_color_key("", ["PWL Temp", "Actual Temp"], ["blue", "green"]))
	push!(pd3nl,Theme(key_position = :top,background_color=colorant"white",major_label_font_size=24pt,line_width=3pt,
	minor_label_font_size=20pt,key_label_font_size=20pt))
end
if drawFig==1 draw(PNG(path*"J_"*updateMethod*"_Temp.png", 24inch, 12inch), pd3) end

if updateMethod=="fastAscent"
	lamLabel=raw"Lambda ($/K)"
else
	lamLabel=raw"Lambda ($/kA)"
end
pd4nl=plot(layer(x=1:horzLen+1,y=Lam[:,convIt],Geom.line,Theme(line_width=3pt)),
		Guide.xlabel("Time"), Guide.ylabel(lamLabel),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=24pt,
		minor_label_font_size=20pt,key_label_font_size=20pt))
if drawFig==1 draw(PNG(path*"J_"*updateMethod*"_Lam.png", 24inch, 12inch), pd4) end

#fName="J_Decentral_notfast.png"
#fName="J_Decentral_fast.png"
#draw(PNG(path*fName, 13inch, 14inch), vstack(pd1,pd2,pd3,pd4))

uSumPlotNL=plot(uSum[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=uSumStarNL,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("U sum"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
iPlotNL=plot(Itotal[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=itotalStarNL,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("Z sum"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
lamPlotNL=plot(Lam[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=lamCurrStarNL,Geom.line,Theme(default_color=colorant"black",line_width=4pt)),
			Guide.xlabel("Time"), Guide.ylabel("Lambda"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_"*updateMethod*"_LamConv.png", 36inch, 12inch), lamPlot) end

fPlotNL=plot(x=1:convIt-1,y=fConvDual[1:convIt-1,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("obj function gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
convItPlotNL=plot(x=1:convIt,y=itConvDual[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap",orientation=:vertical),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
convPlotNL=plot(x=1:convIt,y=ConvDual[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("Lambda Star Gap",orientation=:vertical),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
constPlotNL=plot(x=1:convIt,y=constConvDual[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("Const Gap",orientation=:vertical),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_"*updateMethod*"_Conv.png", 36inch, 12inch), convPlot) end
