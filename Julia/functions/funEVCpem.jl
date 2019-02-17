# run PEM

function pemEVClocal(n,stepI,horzLen,packLen,evS,pemSol)
	epsilon=1e-3
	#desiredSOC=1 # for now
	desiredSOC=evS.Snmin[n] #this should be a ratio?
	prevSOC= if stepI>1 pemSol.Sn[stepI-1,n] else evS.s0[n] end

	# clean up this logic***
	if stepI==1
		existPack=false
	elseif ((prevSOC-desiredSOC)>=-epsilon)
		existPack=false
	elseif (pemSol.Un[stepI-1,n]>0)
		if stepI<=packLen
			prevChar=stepI-1
			existPack=true
		elseif (pemSol.Un[stepI-(packLen),n]==0)
			prevChar=sum(pemSol.Un[(stepI-1):(stepI-(packLen)),n])
			existPack=true
		else
			existPack=false
		end
	else
		existPack=false
	end

	if existPack==true # still using a packet
		Req=-packLen+prevChar
	else
		if (prevSOC-1)>=-epsilon #100% SOC reached
			ratio=0
			Req=0
		elseif (prevSOC-desiredSOC)>=-epsilon #desired SOC reached
			ratio=0
			Req=2 #could use a probability here
		else
			# ratio=(desiredSOC-prevSOC)/(evS.ηP[n]*evS.imax[n]*(evS.Kn[n]-(stepI-1)))
			ratio=(desiredSOC-prevSOC)/(evS.ηP[n]*evS.imax[n]*(evS.Kn[n]-(stepI)))
			if ratio>=1 # opt out (need to charge for rest of time)
				ratio=1
				Req=-packLen
			else
				#mu[stepI,n] = 1/mttr*((desiredSOC-ratio[stepI,n])/(ratio[stepI,n]-0))*((setSOC-0)/(desiredSOC-setSOC))
				mu = 1/mttr*((ratio-0)/(1-ratio))*((1-setSOC)/(setSOC-0))
				P = min(max(1-exp(-mu*evS.Ts),0),1)
				t=rand()
				Req=if (t>(1-P)) 1 else  0  end
			end
		end
	end
	#return ratio,mu,P,Req
	return Req
end

function pemEVCcoord(stepI,horzLen,packLen,evS,pemSol,Req)
	poolSol=Int(min(round(N*packLen),100))
	slackWeight=1e8
	extraWeight=round(1/(2*N*(horzLen+1)),digits=6)

	#receive requests and forecast for the next packLen intervals
	prevT = if stepI>1 pemSol.Tactual[stepI-1] else evS.t0 end
    requiredCh=Int.(Req[stepI,:].<0) # all negative numbers
	requiredInd=findall(x->x==1,requiredCh)
	extraInd=findall(x->x==2,Req[stepI,:])
	optOffInd=findall(x->x==0,Req[stepI,:]) #did not request

	m = Model(solver = GurobiSolver(PoolSearchMode=1,PoolSolutions=poolSol,PoolGap=0,TimeLimit=9/10*evS.Ts))
	@variable(m,u[1:horzLen+1,1:N],Bin) #binary charge variable
	@variable(m,Test[1:horzLen+1])
	@variable(m,slackT)
	objExp=slackWeight*slackT
	if !isempty(extraInd)
		objExp=objExp-extraWeight*sum(sum(u[:,n] for n in extraInd))
	end
	if !isempty(setdiff(1:N,extraInd))
		objExp=objExp-sum(sum(u[:,n] for n in setdiff(1:N,extraInd)))
	end
	@objective(m,Min,objExp)
	@constraint(m,tempCon1,Test[1]>=evS.τP*prevT+evS.γP*(sum(u[1,n]*evS.imax[n] for n=1:N)+evS.iD_pred[stepI])^2+evS.ρP*evS.Tamb[stepI])
	@constraint(m,tempCon2[kk=1:horzLen],Test[kk+1]>=evS.τP*Test[kk]+evS.γP*(sum(u[kk+1,n]*evS.imax[n] for n=1:N)+evS.iD_pred[stepI+kk])^2+evS.ρP*evS.Tamb[stepI+kk])
	@constraint(m,Test.<=evS.Tmax+slackT)
	@constraint(m,slackT>=0)
	@constraint(m,optOnC[nn=1:length(requiredInd)],sum(u[kk,requiredInd[nn]] for kk=1:Int(min(abs(Req[stepI,requiredInd[nn]]),horzLen+1)))==min(abs(Req[stepI,requiredInd[nn]]),horzLen+1))
	@constraint(m,optOffC[kk=1:horzLen+1],sum(u[kk,n] for n in optOffInd)==0)


	if solverSilent
        @suppress_out begin
			statusC = solve(m)
        end
    else
		statusC = solve(m)
    end


	if statusC==:Optimal
		#take random solution if multiple
		solCount=Gurobi.get_sol_count(getrawsolver(m))
		if solCount>1
			println(solCount," multiple solutions found")
			solNum=rand(1:solCount-1)
			setparam!(getrawsolver(m),"SolutionNumber",solNum)
			getparam(getrawsolver(m),"SolutionNumber")
			sol=Gurobi.get_dblattrarray(getrawsolver(m),"Xn",1,N)
		elseif solCount==1
			sol=getvalue(u)[1,:]
		else
			sol=requiredCh
		end
		Rec=sol
	else
		#return required EVs to charges
		Rec=requiredCh
	end

	#R>=1 and packRec get charged
	for n=1:N
		pemSol.Un[stepI,n]= if Rec[n]==1 evS.imax[n] else  0  end
		prevSOC= if stepI>1 pemSol.Sn[stepI-1,n] else evS.s0[n] end
		pemSol.Sn[stepI,n]=prevSOC+evS.ηP[n]*pemSol.Un[stepI,n]
	end

	return nothing
end

function pemEVCstep(stepI,evS,pemSol,silent)


	horzLen=min(packLen,evS.K-stepI)
 	Req=zeros(evS.K,N)

	#println(stepI)
	#send requests
    for n=1:N
		Req[stepI,n]=pemEVClocal(n,stepI,horzLen,packLen,evS,pemSol)
	end

	pemEVCcoord(stepI,horzLen,packLen,evS,pemSol,Req)

	# actual dynamics
	prevT = if stepI>1 pemSol.Tactual[stepI-1] else evS.t0 end
	pemSol.uSum[stepI] = round.(sum(pemSol.Un[stepI,n] for n=1:N),digits=8)
	pemSol.Itotal[stepI] = round.(pemSol.uSum[stepI] + evS.iD_actual[stepI],digits=8)
	pemSol.Tactual[stepI] = round.(evS.τP*prevT+evS.γP*pemSol.Itotal[stepI]^2+evS.ρP*evS.Tamb[stepI],digits=6)
	return nothing
end

function pemEVC(evS::scenarioStruct,slack::Bool,silent::Bool)
	pemSol=solutionStruct(K=evS.K,N=evS.N,S=evS.S)

	Juno.progress() do id
		for stepI=1:evS.K
			@info "$(Dates.format(Dates.now(),"HH:MM:SS")): $(stepI) of $(evS.K)....\n" progress=stepI/evS.K _id=id
			@printf "%s: time step %g of %g....\n" Dates.format(Dates.now(),"HH:MM:SS") stepI evS.K
			try
				pemEVCstep(stepI,evS,pemSol,silent)
			catch e
				@printf "error: %s" e
				break
			end
		end
	end

	pemSol.objVal[1,1]=sum(sum((pemSol.Sn[stepI,n]-1)^2*evS.Qsi[n,1]+(pemSol.Un[stepI,n])^2*evS.Ri[n,1] for n=1:N) for stepI=1:(evS.K))
	return pemSol
end
