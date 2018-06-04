#Micah Botkin-Levy
#4/10/18

tic()
#initialize
#initialize with current states
xn0=s0
xt0=T0

d = Truncated(Normal(0), 0, 50)
lambda0=rand(d, K1+1)
#lambdaGuess=2
#lambda0=ones(K1+1,1)*lambdaGuess
#lambda0=[0;linspace(7,0,K)]
#lambda0=rand(K1+1)/2
lambda=lambda0
if updateMethod=="fastAscent"
	alpha = 0.1
	#alpha=0.001
else
	alpha = 0.01
	#alpha= 0.001
end

stepI = 1;
horzLen=K1
convChk = 1e-16
numIteration=1000
convIt=numIteration
ConvDual=zeros(numIteration,1)
itConvDual=zeros(numIteration,1)
constConvDual=zeros(numIteration,1)
fConvDual=zeros(numIteration,1)
xnConvDual=zeros(numIteration,1)
unConvDual=zeros(numIteration,1)

#u w and z are one index ahead of x. i.e the x[k+1]=x[k]+eta*u[k+1]
Lam=zeros((horzLen+1),numIteration) #(rows are time, columns are iteration)
Lam[:,1]=lambda0
Xt=zeros((horzLen+1),numIteration) #rows are time
Tactual=zeros((horzLen+1),numIteration) #rows are time
#Xt[1,1]=T0
#Xn=zeros(N*(horzLen+1),numIteration) #row are time,columns are iteration
#Xn=zeros((horzLen+1),N) #row are time,columns are iteration
#Xn[:,1]=s0
#Un=zeros(N*(horzLen+1),numIteration) #row are time, columns are iteration
#Un=zeros((horzLen+1),N) #row are time, columns are iteration


# Xn=zeros(N*(horzLen+1),numIteration)
# Un=zeros(N*(horzLen+1),numIteration)


Xn=SharedArray{Float64}(N*(horzLen+1),numIteration)
Un=SharedArray{Float64}(N*(horzLen+1),numIteration)

#iterate at each time step until convergence
for p=1:numIteration-1
    #solve subproblem for each EV
	@sync @parallel for evInd=1:N
	#for evInd=1:N
        target=zeros((horzLen+1),1)
		target[(Kn[evInd,1]-(stepI-1)):1:length(target),1]=Sn[evInd,1]
        evM=Model(solver = GurobiSolver(OutputFlag=0))
        @variable(evM,un[1:horzLen+1])
        @variable(evM,xn[1:horzLen+1])
        objFun(x,u)=sum((x[k,1]-1)^2*Qsi[evInd,1] for k=1:horzLen+1) +
        			sum((u[k,1])^2*Ri[evInd,1]    for k=1:horzLen+1) +
                    sum(lambda[k,1]*u[k,1]        for k=1:horzLen+1)
        @objective(evM,Min, objFun(xn,un))
		@constraint(evM,xn[1,1]==xn0[evInd,1]+eta[evInd,1]*un[1,1]) #fix for MPC loop
		@constraint(evM,[k=1:horzLen],xn[k+1,1]==xn[k,1]+eta[evInd,1]*un[k+1,1]) #check K+1
        @constraint(evM,xn.<=1)
        @constraint(evM,xn.>=target)
        @constraint(evM,un.<=imax[evInd,1])
        @constraint(evM,un.>=imin[evInd,1])

		TT = STDOUT # save original STDOUT stream
		redirect_stdout()
        status = solve(evM)
		redirect_stdout(TT)

		if status!=:Optimal
            break
        else
            Xn[collect(evInd:N:N*(horzLen+1)),p+1]=getvalue(xn) #solved state goes in next time slot
            Un[collect(evInd:N:N*(horzLen+1)),p+1]=getvalue(un) #current go
        end
    end


	if updateMethod=="dualAscent"
	    #solve coordinator problem
	    coorM=Model(solver = GurobiSolver(OutputFlag=0))
	    @variable(coorM,z[1:S*(horzLen+1)])
	    @variable(coorM,xt[1:horzLen+1])
	    @objective(coorM,Min,-sum(lambda[k,1]*sum(z[(k-1)*S+ii,1] for ii=1:S) for k=1:horzLen+1))
		@constraint(coorM,xt[1,1]==tau*xt0+gamma*deltaI*sum((2*m+1)*z[m+1,1] for m=0:S-1)+rho*w[2,1]) #fix for MPC loop
		@constraint(coorM,[k=1:horzLen],xt[k+1,1]==tau*xt[k,1]+gamma*deltaI*sum((2*m+1)*z[k*S+(m+1),1] for m=0:S-1)+rho*w[k*2+2,1])
		if noTlimit==0
			@constraint(coorM,upperTCon,xt.<=Tmax)
		end
	    @constraint(coorM,xt.>=0)
	    @constraint(coorM,z.<=deltaI)
	    @constraint(coorM,z.>=0)
		TT = STDOUT # save original STDOUT stream
		redirect_stdout()
	    status = solve(coorM)
		redirect_stdout(TT)
	    if status!=:Optimal
	        break
		else
			 Xt[:,p+1]=getvalue(xt);
		end

	    #grad of lagragian
		zSum=zeros(horzLen+1,1)
		gradL=zeros(horzLen+1,1)
		for k=1:horzLen+1
			zSum[k,1]=sum(getvalue(z)[(k-1)*(S)+(1:S)])
			gradL[k,1]=sum(Un[(k-1)*N+n,p+1] for n=1:N) + w[(k-1)*2+(stepI*2-1),1] - zSum[k,1]
		end
		constConvDual[p,1]=norm(gradL,2)
	end

	#calculate actual temperature from nonlinear model of XFRM
	ztotal=zeros(horzLen+1,1)
	for k=1:horzLen+1
		ztotal[k,1]=sum(Un[(k-1)*N+n,p+1]    for n=1:N) + w[(k-1)*2+(stepI*2-1),1]
	end
	Tactual[1,p+1]=tau*xt0+gamma*ztotal[1,1]^2+rho*w[2,1] #fix for mpc
	for k=1:horzLen
		Tactual[k+1,p+1]=tau*Tactual[k,p+1]+gamma*ztotal[k+1,1]^2+rho*w[k*2+2,1]  #fix for mpc
	end


	if updateMethod=="fastAscent"
		#fast ascent
		if noTlimit==0
			gradL=Tactual[:,p+1]-Tmax*ones(horzLen+1,1)
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
	if updateMethod=="fastAscent"
		alpha_p = alpha/ceil(p/2)
		#alpha_p = alpha/(p*5)
	else
		alpha_p = alpha/ceil(p/2)
		#alpha_p = alpha/(p*5)
	end

    lambda_new=max.(lambda+alpha_p*gradL,0)
	#lambda_new=lambda+alpha_p*gradL
    Lam[:,p+1]=lambda_new
    lambda=lambda_new

	#check convergence
	objFun(xn,u)=sum(sum((xn[(k-1)*(N)+n,1]-1)^2*Qsi[n,1]     for n=1:N) for k=1:horzLen+1) +
					sum(sum((u[(k-1)*N+n,1])^2*Ri[n,1]           for n=1:N) for k=1:horzLen+1)
	fGap= objFun(Xn[:,p+1],Xt[:,p+1],Un[:,p+1])-fStar
	fGap= objFun(Xn[:,p+1],Un[:,p+1])-fStar
	xnGap=norm((Xn[:,p+1]-xnStar),2)
	unGap=norm((Un[:,p+1]-uStar),2)
	itGap = norm(Lam[:,p+1]-Lam[:,p],2)
	convGap = norm(Lam[:,p+1]-lamTempStar,2)
	fConvDual[p,1]=fGap
	xnConvDual[p,1]=xnGap
	unConvDual[p,1]=unGap
	itConvDual[p,1]=itGap
	ConvDual[p,1]=convGap
	if(convGap <= convChk )
		@printf "Converged after %g iterations\n" p
		convIt=p
		break
	else
		@printf "lastGap %e after %g iterations\n" itGap p
		@printf "convGap %e after %g iterations\n" convGap p
        @printf "xnGap %e after %g iterations\n" xnGap p
		@printf("fGap %e after %g iterations\n\n",fGap,p)

	end
end
toc()

println("plotting....")
xPlot=zeros(horzLen+1,N)
for ii= 1:N
	xPlot[:,ii]=Xn[collect(ii:N:length(Xn[:,convIt])),convIt]
end

pd1=plot(xPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV SOC"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=24pt,line_width=3pt,
		minor_label_font_size=20pt,key_label_font_size=20pt))
if drawFig==1 draw(PNG(path*"J_"*updateMethod*"_SOC.png", 24inch, 12inch), pd1) end


uPlot=zeros(horzLen+1,N)
for ii= 1:N
	uPlot[:,ii]=Un[collect(ii:N:length(Un[:,convIt])),convIt]
end

pd2=plot(uPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV Current (A)"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=24pt,line_width=3pt,
		minor_label_font_size=20pt,key_label_font_size=20pt))
if drawFig==1 draw(PNG(path*"J_"*updateMethod*"_Curr.png", 24inch, 12inch), pd2) end

pd3=plot(layer(x=1:horzLen+1,y=Tactual[:,convIt],Geom.line,Theme(default_color=colorant"green",line_width=3pt)),
		layer(yintercept=[Tmax],Geom.hline(color=["red"],style=:dot),Theme(line_width=3pt)),
		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=24pt,
		minor_label_font_size=20pt,key_label_font_size=20pt))
if updateMethod=="dualAscent"
	push!(pd3,layer(x=1:horzLen+1,y=Xt[:,convIt],Geom.line,Theme(default_color=colorant"blue",line_width=3pt)))
	push!(pd3,Guide.manual_color_key("", ["PWL Temp", "Actual Temp"], ["blue", "green"]))
	push!(pd3,Theme(key_position = :top,background_color=colorant"white",major_label_font_size=24pt,line_width=3pt,
	minor_label_font_size=20pt,key_label_font_size=20pt))
end
if drawFig==1 draw(PNG(path*"J_"*updateMethod*"_Temp.png", 24inch, 12inch), pd3) end

if updateMethod=="fastAscent"
	lamLabel=raw"Lambda ($/K)"
else
	lamLabel=raw"Lambda ($/A)"
end
pd4=plot(layer(x=1:horzLen+1,y=Lam[:,convIt],Geom.line,Theme(line_width=3pt)),
		Guide.xlabel("Time"), Guide.ylabel(lamLabel),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=24pt,
		minor_label_font_size=20pt,key_label_font_size=20pt))
if drawFig==1 draw(PNG(path*"J_"*updateMethod*"_Lam.png", 24inch, 12inch), pd4) end


#fName="J_Decentral_notfast.png"
#fName="J_Decentral_fast.png"
#draw(PNG(path*fName, 13inch, 14inch), vstack(pd1,pd2,pd3,pd4))



lamPlot=plot(Lam[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel("Lambda"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_"*updateMethod*"_LamConv.png", 36inch, 12inch), lamPlot) end

convItPlot=plot(x=1:convIt,y=itConvDual[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
convPlot=plot(x=1:convIt,y=ConvDual[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_"*updateMethod*"_Conv.png", 36inch, 12inch), convPlot) end



#compare central and decentral current agg
# aggU=plot(layer(x=1:horzLen+1,y=sum(uPlot[:,i] for i=1:N),Geom.line,Theme(default_color=colorant"blue")),
# 		layer(x=1:horzLen+1,y=sum(Un[:,i] for i=1:N),Geom.line,Theme(default_color=colorant"green")),
# 		Guide.xlabel("Time"), Guide.ylabel("PEV Current (A)"),
# 		Coord.Cartesian(xmin=0,xmax=horzLen+1),
# 		Theme(background_color=colorant"white"),
# 		Guide.manual_color_key("", ["Central", "Decentral (fast)"], ["blue", "green"]))
#draw(PNG(path*"aggPlot_fast.png", 13inch, 8inch), aggU)


#compare dual ascent and fast ascent convergence

# dualConv=Conv

# fastConv=Conv
# convPlotcomp=plot(layer(x=1:convIt,y=dualConv[1:convIt,1],Geom.line,Theme(default_color=colorant"blue")),
# 			  layer(x=1:convIt,y=fastConv[1:convIt,1],Geom.line,Theme(default_color=colorant"green")),
# 			  Guide.xlabel("Iteration"), Guide.ylabel("Convergance Gap"),
# 			  Coord.Cartesian(xmin=0,xmax=10),Theme(background_color=colorant"white"),
# 		  	  Guide.manual_color_key("", ["Dual Ascent", "Fast Ascent"], ["blue", "green"]))
# draw(PNG(path*"convComp_fast.png", 13inch, 8inch), aggU)
