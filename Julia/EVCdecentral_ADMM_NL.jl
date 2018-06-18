#EVC with ADMM for PWL convex relaxation
#Micah Botkin-Levy
#5/22/18

#formulation try 2
#u, sn, xt, and z are all in "x" v is the auxilliary variable corresponding to z in literature
#current constraint is coupling

tic()
#initialize with current states
sn0=s0
xt0=T0

stepI = 1;
horzLen=K1
convChk = 1e-8
numIteration=100
convIt=numIteration
Conv=zeros(numIteration,1)
itConv=zeros(numIteration,1)
constConv=zeros(numIteration,1)
fConv=zeros(numIteration,1)
snConv=zeros(numIteration,1)

#admm  initial parameters and guesses
#ρADMM=10.0^(0)
ρADMM=1e6

# lambda0=lamCurrStarNL
# vi0=-itotalStarNL
# vu0=uStarNL

lambda0=1000*ones(horzLen+1,1)
vi0=-ones(horzLen+1,1)
vu0=.01*ones(N*(horzLen+1),1)

#u w and z are one index ahead of x. i.e the x[k+1]=x[k]+eta*u[k+1]
Lam=zeros((horzLen+1),numIteration) #(rows are time, columns are iteration)

Un=SharedArray{Float64}(N*(horzLen+1),numIteration) #row are time,  columns are iteration
Sn=SharedArray{Float64}(N*(horzLen+1),numIteration)  #row are time,  columns are iteration
Xt=zeros((horzLen+1),numIteration)  #row are time,  columns are iteration
Itotal=zeros((horzLen+1),numIteration)  #row are time,  columns are iteration
Vu=zeros((N)*(horzLen+1),numIteration) #row are time,  columns are iteration
Vi=zeros((horzLen+1),numIteration)

Vi[:,1]=vi0
Vu[:,1]=vu0
Lam[:,1]=lambda0

#for debugging
CC=zeros((horzLen+1),numIteration)  #row are time,  columns are iteration

for p in 1:numIteration-1
	#ρ_p = ρADMM/ceil(p/2)
    ρ_p = ρADMM

    #x minimization eq 7.66 in Bertsekas
    @sync @parallel for evInd=1:N
        evV=Vu[collect(evInd:N:length(Vu[:,p])),p]
        target=zeros((horzLen+1),1)
		target[(Kn[evInd,1]-(stepI-1)):1:length(target),1]=Snmin[evInd,1]
    	evM = Model(solver = GurobiSolver())
    	@variable(evM,sn[1:(horzLen+1)])
    	@variable(evM,u[1:(horzLen+1)])
		@objective(evM,Min, sum((sn[k,1]-1)^2*Qsi[evInd,1]+(u[k,1])^2*Ri[evInd,1]+
                                Lam[k,p]*(u[k,1]-evV[k,1])+
                                ρ_p/2*(u[k,1]-evV[k,1])^2 for k=1:horzLen+1))
        @constraint(evM,sn[1,1]==sn0[evInd,1]+ηP[evInd,1]*u[1,1])
        @constraint(evM,[k=1:horzLen],sn[k+1,1]==sn[k,1]+ηP[evInd,1]*u[k+1,1])
    	@constraint(evM,sn.<=1)
    	@constraint(evM,sn.>=target)
        @constraint(evM,u.<=imax[evInd,1])
        @constraint(evM,u.>=imin[evInd,1])
    	TT = STDOUT # save original STDOUT stream
    	redirect_stdout()
    	status = solve(evM)
    	redirect_stdout(TT)
    	if status!=:Optimal
    	    return
    	else
    		Sn[collect(evInd:N:length(Sn[:,p+1])),p+1]=getvalue(sn)
    		Un[collect(evInd:N:length(Un[:,p+1])),p+1]=getvalue(u)
    	end
    end

    #N+1 decoupled problem aka transformer current
    tM = Model(solver = IpoptSolver())
    #@variable(tM,z[1:(S)*(horzLen+1)])
    @variable(tM,xt[1:(horzLen+1)])
    @variable(tM,itotal[1:(horzLen+1)])
	# constFun1(u,v)=sum(Lam[k,p]*(u[k,1]-v[k,1])  for k=1:(horzLen+1))
	# constFun2(u,v)=ρADMM/2*sum((u[k,1]-v[k,1])^2  for k=1:(horzLen+1))
    # @objective(tM,Min, constFun1(-itotal,Vi[:,p])+constFun2(-itotal,Vi[:,p]))
    @objective(tM,Min,sum(Lam[k,p]*(-itotal[k,1]-Vi[k,p])+
                          ρADMM/2*(-itotal[k,1]-Vi[k,p])^2  for k=1:(horzLen+1)))
    @NLconstraint(tM,tempCon1,xt[1,1]==τP*xt0+γP*(itotal[1,1])^2+ρP*w[stepI*2,1])
    @NLconstraint(tM,tempCon2[k=1:horzLen],xt[k+1,1]==τP*xt[k,1]+γP*(itotal[k+1,1])^2+ρP*w[stepI*2+k*2,1])
    if noTlimit==0
    	@constraint(tM,upperTCon,xt.<=Tmax)
    end
    @constraint(tM,xt.>=0)
    @constraint(tM,itotal.>=0)
    @constraint(tM,itotal.<=ItotalMax)

    TT = STDOUT # save original STDOUT stream
    redirect_stdout()
    status = solve(tM)
    redirect_stdout(TT)
    if status!=:Optimal
        return
    else
        Xt[:,p+1]=getvalue(xt)
        Itotal[:,p+1]=getvalue(itotal)
    end

    #lambda update eq 7.68
    currConst=zeros(horzLen+1,1)
	for k=1:horzLen+1
		uSum[k,p+1]=sum(Un[(k-1)*N+n,p+1] for n=1:N)
		currConst[k,1]=uSum[k,p+1] + w[(k-1)*2+(stepI*2-1),1] - Itotal[k,p+1]
		#Lam[k,p+1]=max.(Lam[k,p]+ρ_p/(S*(N))*(currConst[k,1]),0)
		Lam[k,p+1]=Lam[k,p]+ρ_p/(S*(N))*(currConst[k,1])
	end
	CC[:,p+1]=currConst

    #v upate eq 7.67
    for k=1:horzLen+1
        Vu[(k-1)*N+collect(1:N),p+1]=min.(max.(Un[(k-1)*N+collect(1:N),p+1]+(Lam[k,p]-Lam[k,p+1])/ρ_p,imin),imax)
        #Vi[(k-1)*(S)+collect(1:S),p+1]=-Z[(k-1)*(S)+collect(1:S),p+1]+(Lam[k,p]-Lam[k,p+1])/ρ_p
		Vi[k,p+1]=max.(min.(-Itotal[k,p+1]+(Lam[k,p]-Lam[k,p+1])/ρ_p,0),-ItotalMax)
    end


    #check convergence
	objFun(sn,xt,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*Qsi[n,1]     for n=1:N) for k=1:horzLen+1) +
					sum((xt[k,1]-1)^2*Qsi[N+1,1]                 for k=1:horzLen+1) +
					sum(sum((u[(k-1)*N+n,1])^2*Ri[n,1]           for n=1:N) for k=1:horzLen+1)
	fGap= objFun(Sn[:,p+1],Xt[:,p+1],Un[:,p+1])-fStarNL
	#fGap= objFun(Sn[:,p],Xt[:,p],Un[:,p])-fStar
	snGap=norm((Sn[:,p+1]-snStarNL),2)
	constGap=norm(currConst,2)
	itGap = norm(Lam[:,p+1]-Lam[:,p],2)
	convGap = norm(Lam[:,p+1]-lamCurrStarNL,2)
	fConv[p,1]=fGap
	snConv[p,1]=snGap
	constConv[p,1]=constGap
	itConv[p,1]=itGap
	Conv[p,1]=convGap
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
for ii= 1:N
	xPlot[:,ii]=Sn[collect(ii:N:length(Sn[:,convIt])),convIt]
end

#plot(x=1:horzLen+1,y=xPlot2[:,ii])
# plot(x=1:Kn[ii,1],y=xPlot2[1:Kn[ii,1],ii])

pd1NLadmm=plot(xPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV SOC"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_centralNL_ADMM_SOC.png", 24inch, 12inch), pd1NLadmm) end


uPlot=zeros(horzLen+1,N)
for ii= 1:N
	uPlot[:,ii]=Un[collect(ii:N:length(Un[:,convIt])),convIt]
end

pd2NLadmm=plot(uPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV Current (A)"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_centralNL_ADMM_Curr.png", 24inch, 12inch), pd2NLadmm) end

pd3NLadmm=plot(layer(x=1:horzLen+1,y=Xt[:,convIt],Geom.line,Theme(default_color=colorant"blue")),
		#layer(x=1:horzLen+1,y=Tactual,Geom.line,Theme(default_color=colorant"green")),
		yintercept=[Tmax],Geom.hline(color=["red"],style=:dot),
		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",key_position = :top,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
		#Guide.manual_color_key("", ["PWL Temp", "Actual Temp"], ["blue", "green"]))
if drawFig==1 draw(PNG(path*"J_centralNL_ADMM_Temp.png", 24inch, 12inch), pd3NLadmm) end

pd4NLadmm=plot(x=1:horzLen+1,y=Lam[:,convIt],Geom.line,
		Guide.xlabel("Time"), Guide.ylabel(raw"Lambda ($/A)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_centralNL_ADMM_Lam.png", 24inch, 12inch), pd4NLadmm) end

fName="J_Central.png"


lamPlotNLadmm=plot(Lam[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
            layer(x=1:horzLen+1,y=lamCurrStarNL,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("Lambda"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
constPlotNLadmm2=plot(CC[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel("curr constraint diff"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_ADMM_LamConv.png", 36inch, 12inch), lamPlotadmm) end

convItPlotNLadmm=plot(x=1:convIt,y=itConv[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
convPlotNLadmm=plot(x=1:convIt,y=Conv[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
fPlotadmm=plot(x=1:convIt-1,y=fConv[1:convIt-1,1],Geom.line,#Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
constPlotadmm=plot(x=1:convIt,y=constConv[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_ADMM_Conv.png", 36inch, 12inch), convPlotadmm) end
