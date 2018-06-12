#Micah Botkin-Levy
#4/10/18

tic()

lambda0=1
#lambda0=lamCurrStar1

if updateMethod=="fastAscent"
	alpha = 0.1
else
	alpha = 500
end

stepI = 1;
convChk = 1e-8
maxIt=1000
convIt=maxIt


ConvDual=zeros(maxIt,1)
itConvDual=zeros(maxIt,1)
constConvDual=zeros(maxIt,1)
fConvDual=zeros(maxIt,1)
snConvDual=zeros(maxIt,1)
unConvDual=zeros(maxIt,1)

#u w and z are one index ahead of x. i.e the x[k+1]=x[k]+eta*u[k+1]
Lam=zeros(1,maxIt) #(rows are time, columns are iteration)
Lam[1,1]=lambda0
Xt=zeros(1,maxIt) #rows are time
Tactual=zeros(1,maxIt) #rows are time

Sn=SharedArray{Float64}(N,maxIt)
Un=SharedArray{Float64}(N,maxIt)

#iterate at each time step until convergence
for p=1:maxIt-1
    #solve subproblem for each EV
	@sync @parallel for evInd=1:N
        evM=Model(solver = GurobiSolver(OutputFlag=0))
        @variable(evM,un)
        @variable(evM,sn)
        objFun(x,u)=sum((x-1)^2*Qsi[evInd,1] +
        			(u)^2*Ri[evInd,1]   +
                    Lam[1,p]*u)
        @objective(evM,Min, objFun(sn,un))
		@constraint(evM,sn==sn0[evInd,1]+etaP[evInd,1]*un)
        @constraint(evM,sn<=1)
        @constraint(evM,sn>=target)
        @constraint(evM,un<=imax[evInd,1])
        @constraint(evM,un>=imin[evInd,1])

		TT = STDOUT # save original STDOUT stream
		redirect_stdout()
        status = solve(evM)
		redirect_stdout(TT)

		if status!=:Optimal
            break
        else
            Sn[evInd,p+1]=getvalue(sn) #solved state goes in next time slot
            Un[evInd,p+1]=getvalue(un) #current go
        end
    end


	if updateMethod=="dualAscent"
	    #solve coordinator problem
	    #coorM=Model(solver = GurobiSolver(OutputFlag=0))
		coorM=Model(solver = IpoptSolver())
	    @variable(coorM,z[1:S])
	    @variable(coorM,xt)
	    @objective(coorM,Min,sum(-Lam[1,p]*sum(z[s,1] for s=1:S)))
		@constraint(coorM,xt==tauP*xt0+gammaP*deltaI*sum((2*m+1)*z[m+1,1] for m=0:S-1)+rhoP*w[2,1])
		if noTlimit==0
			@constraint(coorM,upperTCon,xt<=Tmax)
		end
	    @constraint(coorM,xt>=0)
	    @constraint(coorM,z.<=deltaI)
	    @constraint(coorM,z.>=0)
		TT = STDOUT # save original STDOUT stream
		redirect_stdout()
	    status = solve(coorM)
		redirect_stdout(TT)
	    if status!=:Optimal
	        return
		else
			 Xt[1,p+1]=getvalue(xt)
		end

	    #grad of lagragian
		zSum=sum(getvalue(z)[s] for s=1:S)
		gradL=sum(Un[n,p+1] for n=1:N) + w[1,1] - zSum
		constConvDual[p,1]=norm(gradL,2)
	end

	#calculate actual temperature from nonlinear model of XFRM
	ztotal=sum(Un[n,p+1]  for n=1:N) + w[1,1]
	Tactual[1,p+1]=tauP*xt0+gammaP*ztotal^2+rhoP*w[2,1]

	if updateMethod=="fastAscent"
		gradL=Tactual[1,p+1]-Tmax
	end


    #update lambda
	if updateMethod=="fastAscent"
		alpha_p = alpha/ceil(p/2)
		#alpha_p = alpha/(p*5)
	else
		alpha_p = alpha/ceil(p/2)
		#alpha_p = alpha/(p*5)
	end

	#lambda_new=lambda+alpha_p*gradL
    Lam[1,p+1]=max.(Lam[1,p]+alpha_p*gradL,0)

	#check convergence
	objFun(sn,u)=sum((sn[n,1]-1)^2*Qsi[n,1]   for n=1:N)  +
					sum((u[n,1])^2*Ri[n,1]   for n=1:N)
	fGap= objFun(Sn[:,p+1],Un[:,p+1])-fStar1
	snGap=norm((Sn[:,p+1]-snStar1),2)
	unGap=norm((Un[:,p+1]-uStar1),2)
	itGap = norm(Lam[:,p+1]-Lam[:,p],2)
	if updateMethod=="fastAscent"
		convGap = norm(Lam[:,p+1]-lamTempStar1,2)
	else
		convGap = norm(Lam[:,p+1]-lamCurrStar1,2)
	end
	fConvDual[p,1]=norm(fGap,2)
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



lamPlot=plot(x=1:convIt,y=Lam[1,1:convIt])

fPlot=plot(x=1:convIt-1,y=fConvDual[1:convIt-1,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("obj function gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
convItPlot=plot(x=1:convIt,y=itConvDual[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
convPlot=plot(x=1:convIt,y=ConvDual[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("Lambda Star Gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
