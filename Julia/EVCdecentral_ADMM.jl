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
ConvADMM=zeros(maxIt,1)
itConvADMM=zeros(maxIt,1)
constConvADMM=zeros(maxIt,1)
fConvADMM=zeros(maxIt,1)
snConvADMM=zeros(maxIt,1)
unConvADMM=zeros(maxIt,1)

#admm  initial parameters and guesses
#ρADMM=10.0^(0)
ρADMM=10^6 #for kA
#ρADMM=1    #for A


lambda0=1000*ones(horzLen+1,1)
vz0=-ones(S*(horzLen+1),1)
vu0=.01*ones(N*(horzLen+1),1)
#vz0=-zStar
#vu0=uStar
#lambda0=lamCurrStar

#u w and z are one index ahead of x. i.e the x[k+1]=x[k]+η*u[k+1]
Un=SharedArray{Float64}(N*(horzLen+1),maxIt) #row are time,  columns are iteration
#Un[:,1]=u0
Lam=zeros((horzLen+1),maxIt) #(rows are time, columns are iteration)
Sn=SharedArray{Float64}(N*(horzLen+1),maxIt)  #row are time,  columns are iteration
Xt=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration
Z=zeros(S*(horzLen+1),maxIt)  #row are time,  columns are iteration
Tactual=zeros((horzLen+1),maxIt) #rows are time
Vu=zeros((N)*(horzLen+1),maxIt) #row are time,  columns are iteration
Vz=zeros(S*(horzLen+1),maxIt)

Lam[:,1]=lambda0
Vz[:,1]=vz0
Vu[:,1]=vu0

#for debugging
CC=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration
ZS=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration
US=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration


for p in 1:maxIt-1
	#try
		#ρ_p = ρADMM/ceil(p/2)
		ρI = ρADMM
	    #x minimization eq 7.66 in Bertsekas
	    @sync @parallel for evInd=1:N
			lambda=Lam[:,p]
	        evV=Vu[collect(evInd:N:length(Vu[:,p])),p]
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
	    tM = Model(solver = GurobiSolver())
	    @variable(tM,z[1:(S)*(horzLen+1)])
	    @variable(tM,xt[1:(horzLen+1)])
	    # constFun1(u,v)=sum(Lam[k,p]*sum(u[(k-1)*(S)+s,1]-v[(k-1)*(S)+s,1] for s=1:S)  for k=1:(horzLen+1))
	    # constFun2(u,v)=ρI/2*sum(sum((u[(k-1)*(S)+s,1]-v[(k-1)*(S)+s,1])*(u[(k-1)*(S)+s,1]-v[(k-1)*(S)+s,1]) for s=1:S)  for k=1:(horzLen+1))
	    # @objective(tM,Min, constFun1(-z,Vz[:,p])+constFun2(-z,Vz[:,p]))
		@objective(tM,Min,sum(Lam[k,p]*(sum(-z[(k-1)*(S)+s,1] for s=1:S)-sum(Vz[(k-1)*(S)+s,p] for s=1:S)) +
							ρI/2*(sum(-z[(k-1)*(S)+s,1] for s=1:S)-sum(Vz[(k-1)*(S)+s,p] for s=1:S))^2  for k=1:(horzLen+1)))
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
	    status = solve(tM)
	    redirect_stdout(TT)
	    if status!=:Optimal
	        return
	    else
	        Xt[:,p+1]=getvalue(xt)
	        Z[:,p+1]=getvalue(z)
	    end

	    #lambda update eq 7.68
	    currConst=zeros(horzLen+1,1)
		zSum=zeros(horzLen+1,1)
		uSum=zeros(horzLen+1,1)
		#zSum=getvalue(zSum)
		for k=1:horzLen+1
			uSum[k,1]=sum(Un[(k-1)*N+n,p+1] for n=1:N)
			zSum[k,1]=sum(Z[(k-1)*(S)+s,p+1] for s=1:S)
			currConst[k,1]= uSum[k,1] + evS.w[(k-1)*2+(stepI*2-1),1] - zSum[k,1]
			#Lam[k,p+1]=max.(Lam[k,p]+ρADMM/(horzLen+1)*(currConst[k,1]),0)
			Lam[k,p+1]=Lam[k,p]+ρI/(S*(N))*(currConst[k,1])
		end

		#calculate actual temperature from nonlinear model of XFRM
		Tactual[1,p+1]=evS.τP*xt0+evS.γP*zSum[1,1]^2+evS.ρP*evS.w[2,1] #fix for mpc
		for k=1:horzLen
			Tactual[k+1,p+1]=evS.τP*Tactual[k,p+1]+evS.γP*zSum[k+1,1]^2+evS.ρP*evS.w[k*2+2,1]  #fix for mpc
		end

		CC[:,p+1]=currConst
		ZS[:,p+1]=zSum
		US[:,p+1]=uSum

	    #v upate eq 7.67
	    for k=1:horzLen+1
	        Vu[(k-1)*N+collect(1:N),p+1]=min.(max.(Un[(k-1)*N+collect(1:N),p+1]+(Lam[k,p]-Lam[k,p+1])/ρI,evS.imin),evS.imax)
	        Vz[(k-1)*(S)+collect(1:S),p+1]=max.(min.(-Z[(k-1)*(S)+collect(1:S),p+1]+(Lam[k,p]-Lam[k,p+1])/ρI,0),-evS.deltaI)
			#Vz[k,p+1]=-zSum[k,1]+(Lam[k,p]-Lam[k,p+1])/ρADMM
	    end

	    #check convergence
		objFun(sn,xt,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,1]     for n=1:N) for k=1:horzLen+1) +
						sum((xt[k,1]-1)^2*evS.Qsi[N+1,1]                 for k=1:horzLen+1) +
						sum(sum((u[(k-1)*N+n,1])^2*evS.Ri[n,1]           for n=1:N) for k=1:horzLen+1)
		fGap= abs(objFun(Sn[:,p+1],Xt[:,p+1],Un[:,p+1])-fStar)
		#fGap= objFun(Sn[:,p],Xt[:,p],Un[:,p])-fStar
		snGap=norm((Sn[:,p+1]-snStar),2)
		unGap=norm((Un[:,p+1]-uStar),2)
		constGap=norm(currConst,2)
		itGap = norm(Lam[:,p+1]-Lam[:,p],2)
		convGap = norm(Lam[:,p+1]-lamCurrStar,2)
		fConvADMM[p,1]=fGap
		snConvADMM[p,1]=snGap
		unConvADMM[p,1]=unGap
		constConvADMM[p,1]=constGap
		itConvADMM[p,1]=itGap
		ConvADMM[p,1]=convGap
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
	# catch e
	# 	@printf "error %s after %g iterations\n" e p
	# end
end
toc()




println("plotting....")
xPlot=zeros(horzLen+1,N)
uPlot=zeros(horzLen+1,N)
for ii= 1:N
	xPlot[:,ii]=Sn[collect(ii:N:length(Sn[:,convIt])),convIt]
	uPlot[:,ii]=Un[collect(ii:N:length(Un[:,convIt])),convIt]
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

pd3admm=plot(layer(x=1:horzLen+1,y=Xt[:,convIt],Geom.line,Theme(default_color=colorant"blue")),
		layer(x=1:horzLen+1,y=Tactual[:,convIt],Geom.line,Theme(default_color=colorant"green")),
		yintercept=[evS.Tmax],Geom.hline(color=["red"],style=:dot),
		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",key_position = :top,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt),
		Guide.manual_color_key("", ["PWL Temp", "Actual Temp"], ["blue", "green"]))
if drawFig==1 draw(PNG(path*"J_decentral_ADMM_Temp.png", 24inch, 12inch), pd3admm) end

pd4admm=plot(x=1:horzLen+1,y=Lam[:,convIt],Geom.line,
		layer(x=1:horzLen+1,y=lamCurrStar,Geom.line,Theme(default_color=colorant"black")),
		Guide.xlabel("Time"), Guide.ylabel(raw"Lambda ($/A)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_decentral_ADMM_Lam.png", 24inch, 12inch), pd4admm) end


lamPlotadmm=plot(Lam[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=lamCurrStar,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("Lambda"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
constPlotadmm2=plot(CC[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel("curr constraint diff"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
zSumPlotadmm=plot(ZS[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=zSumStar,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("Z sum"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
uSumPlotadmm=plot(US[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=uSumStar,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("U sum"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_ADMM_LamConv.png", 36inch, 12inch), lamPlotadmm) end

convItPlotadmm=plot(x=1:convIt,y=itConvADMM[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
convPlotadmm=plot(x=1:convIt,y=ConvADMM[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("central lambda gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
fPlotadmm=plot(x=1:convIt-1,y=fConvADMM[1:convIt-1,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("obj function gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
constPlotadmm=plot(x=1:convIt,y=constConvADMM[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("curr constraint Gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_ADMM_Conv.png", 36inch, 12inch), vstack(convItPlotadmm,convPlotadmm,fPlotadmm,constPlotadmm)) end
