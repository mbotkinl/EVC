# run PEM

function pemEVClocal(n,stepI,horzLen,packLen,evS,pemSol, startPacket)
	epsilon=1e-3
	#desiredSOC=1 # for now
	desiredSOC=evS.Snmin[n] #this should be a ratio?
	prevSOC= if stepI>1 pemSol.Sn[stepI-1,n] else evS.s0[n] end

	if stepI==1
		existPack=false
	elseif ((prevSOC-desiredSOC)>=-epsilon)
		existPack=false
	elseif sum(startPacket[max(stepI-(packLen-1),1):stepI-1,n])>=1
		existPack=true
		prevChar = packLen-1 # this is not dynamic need to count how many timesteps back the "1" is in startpacket
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
	backgroundForecast = min(20,evS.K-stepI)
	#receive requests and forecast for the next packLen intervals
	prevT = if stepI>1 pemSol.Tactual[stepI-1] else evS.t0 end
    requiredCh=Int.(Req.<0) # all negative numbers
	requiredInd=findall(x->x==1,requiredCh)
	extraInd=findall(x->x==2,Req)
	optOffInd=findall(x->x==0,Req) #did not request

	if stepI==280
		Rec=requiredCh
	else

		m = Model(solver = GurobiSolver(PoolSearchMode=1,PoolSolutions=poolSol,PoolGap=0,TimeLimit=9/10*evS.Ts))
		@variable(m,u[1:horzLen,1:N],Bin) #binary charge variable
		@variable(m,Test[1:backgroundForecast])
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
		@constraint(m,tempCon2[kk=1:horzLen-1],Test[kk+1]>=evS.τP*Test[kk]+evS.γP*(sum(u[kk+1,n]*evS.imax[n] for n=1:N)+evS.iD_pred[stepI+kk])^2+evS.ρP*evS.Tamb[stepI+kk])
		@constraint(m,tempCon3[kk=horzLen:backgroundForecast-1],Test[kk+1]>=evS.τP*Test[kk]+evS.ρP*evS.Tamb[stepI+kk])
		@constraint(m,Test.<=evS.Tmax+slackT)
		@constraint(m,slackT>=0)
		@constraint(m,optOnC[nn=1:length(requiredInd)],sum(u[kk,requiredInd[nn]] for kk=1:Int(min(abs(Req[requiredInd[nn]]),horzLen)))==min(abs(Req[requiredInd[nn]]),horzLen))
		@constraint(m,optOffC[kk=1:horzLen],sum(u[kk,n] for n in optOffInd)==0)

		if solverSilent
	        @suppress_out begin
				statusC = solve(m)
	        end
	    else
			statusC = solve(m)
	    end

		pemSol.timeSolve[stepI,1]=getsolvetime(m)

		if statusC==:Optimal
			# temperature slack active
			if getvalue(slackT)>1e-2
				println("PEM coordinator: temp slack active")
			end

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
			println("Infeasible coordinator")
			#return required EVs to charges
			Rec=requiredCh
		end
	end


	newPackets = zeros(N)
	#R>=1 and packRec get charged
	for n=1:N
		pemSol.Un[stepI,n]= if Rec[n]==1 evS.imax[n] else  0  end
		prevSOC= if stepI>1 pemSol.Sn[stepI-1,n] else evS.s0[n] end
		pemSol.Sn[stepI,n]=prevSOC+evS.ηP[n]*pemSol.Un[stepI,n]

		newPackets[n] = if ((Req[n]>=0) && (Rec[n]==1)) 1 else 0 end
	end

	return newPackets
end

function pemEVCstep(stepI,evS,pemSol,startPacket,silent)
	#function pemEVCstep(stepI,evS,pemSol,silent)

	Req=zeros(N)
	horzLen=min(packLen,evS.K-stepI)

	#println(stepI)
	#send requests
    for n=1:N
		Req[n]=pemEVClocal(n,stepI,horzLen,packLen,evS,pemSol, startPacket)
	end

	newPackets = pemEVCcoord(stepI,horzLen,packLen,evS,pemSol,Req)

	# actual dynamics
	prevT = if stepI>1 pemSol.Tactual[stepI-1] else evS.t0 end
	pemSol.uSum[stepI] = round.(sum(pemSol.Un[stepI,n] for n=1:N),digits=8)
	pemSol.Itotal[stepI] = round.(pemSol.uSum[stepI] + evS.iD_actual[stepI],digits=8)
	pemSol.Tactual[stepI] = round.(evS.τP*prevT+evS.γP*pemSol.Itotal[stepI]^2+evS.ρP*evS.Tamb[stepI],digits=6)
	return Req, newPackets
end

function pemEVC(evS::scenarioStruct,slack::Bool,silent::Bool)
	pemSol=solutionStruct(K=evS.K,N=evS.N,S=evS.S)
	Req=zeros(evS.K,N)
	startPacket=zeros(evS.K,N)

	Juno.progress() do id
		for stepI=1:evS.K
			@info "$(Dates.format(Dates.now(),"HH:MM:SS")): $(stepI) of $(evS.K)....\n" progress=stepI/evS.K _id=id
			@printf "%s: time step %g of %g....\n" Dates.format(Dates.now(),"HH:MM:SS") stepI evS.K
			try
				Req[stepI,:], startPacket[stepI,:]=pemEVCstep(stepI,evS,pemSol,startPacket,silent);
			catch e
				@printf "error: %s" e
				break
			end
		end
	end

	pemSol.objVal[1,1]=sum(sum((pemSol.Sn[stepI,n]-1)^2*evS.Qsi[n,1]+(pemSol.Un[stepI,n])^2*evS.Ri[n,1] for n=1:N) for stepI=1:(evS.K))
	return pemSol, Req
end

#HUB PEM
function pemHublocal(h,stepI,horzLen,packLen,hubS,pemSol)
	epsilon=1e-3
	#desiredSOC=1 # for now
	# desiredSOC, desiredKn=findmax(hubS.eDepart_min[stepI:hubS.K,h])
	# desiredSOC, desiredKn=findmax(hubS.eMax[stepI:hubS.K,h])
	ii=findfirst(hubS.eDepart_min[stepI:hubS.K,h].>0)
	if isnothing(ii)
		desiredKn=stepI+horzLen
	else
		desiredKn=stepI+ii-1
	end
	desiredSOC=hubS.eMax[stepI+horzLen,h]

	prevSOC= if stepI>1 pemSol.E[stepI-1,h] else hubS.e0[h] end

	# clean up this logic***
	if stepI==1
		existPack=false
	elseif ((prevSOC-desiredSOC)>=-epsilon)
		existPack=false
	elseif (pemSol.U[stepI-1,h]>0)
		if stepI<=packLen
			prevChar=stepI-1
			existPack=true
		elseif (pemSol.U[stepI-(packLen),h]==0)  # this is not correct need to fix
			prevChar=sum(pemSol.U[(stepI-1):(stepI-(packLen)),h])
			existPack=true
		else
			existPack=false
		end
	else
		existPack=false
	end

	if existPack==true # still using a packet
		ReqH=-packLen+prevChar
	else
		if (prevSOC-hubS.eMax[stepI,h])>=-epsilon #100% SOC reached
			ratio=0
			ReqH=0
		elseif (prevSOC-desiredSOC)>=-epsilon #desired SOC reached
			ratio=0
			ReqH=2 #could use a probability here
		else
			# ratio=(desiredSOC-prevSOC)/(hubS.ηP[n]*hubS.imax[n]*(hubS.Kn[n]-(stepI-1)))
			ratio=(desiredSOC-prevSOC)/(hubS.ηP[h]*hubS.uMax[stepI,h]*(desiredKn-(stepI)))
			if ratio>=1 # opt out (need to charge for rest of time)
				ratio=1
				ReqH=-packLen
			else
				#mu[stepI,n] = 1/mttr*((desiredSOC-ratio[stepI,n])/(ratio[stepI,n]-0))*((setSOC-0)/(desiredSOC-setSOC))
				mu = 1/mttr*((ratio-0)/(1-ratio))*((1-setSOC)/(setSOC-0))
				P = min(max(1-exp(-mu*hubS.Ts),0),1)
				t=rand()
				ReqH=if (t>(1-P)) 1 else  0  end
			end
		end
	end
	#return ratio,mu,P,Req
	return ReqH
end

function pemHubcoord(stepI,horzLen,packLen,hubS,pemSol,ReqI)
	poolSol=Int(min(round(H*packLen),100))
	slackWeight=1e8
	extraWeight=round(1/(2*H*(horzLen+1)),digits=6)

	#receive requests and forecast for the next packLen intervals
	prevT = if stepI>1 pemSol.Tactual[stepI-1] else hubS.t0 end
    requiredCh=Int.(ReqI.<0) # all negative numbers
	requiredInd=findall(x->x==1,requiredCh)
	extraInd=findall(x->x==2,ReqI)
	optOffInd=findall(x->x==0,ReqI) #did not request

	m = Model(solver = GurobiSolver(PoolSearchMode=1,PoolSolutions=poolSol,PoolGap=0,TimeLimit=9/10*hubS.Ts))
	@variable(m,u[1:horzLen+1,1:H],Bin) #binary charge variable
	@variable(m,Test[1:horzLen+1])
	@variable(m,slackT)
	objExp=slackWeight*slackT
	if !isempty(extraInd)
		objExp=objExp-extraWeight*sum(sum(u[:,h] for h in extraInd))
	end
	if !isempty(setdiff(1:H,extraInd))
		objExp=objExp-sum(sum(u[:,h] for h in setdiff(1:H,extraInd)))
	end
	@objective(m,Min,objExp)

	#need to index uMax
	@constraint(m,tempCon1,Test[1]>=hubS.τP*prevT+hubS.γP*(sum(u[1,h]*hubS.uMax[stepI,h] for h=1:H)+hubS.iD_pred[stepI])^2+hubS.ρP*hubS.Tamb[stepI])
	@constraint(m,tempCon2[kk=1:horzLen],Test[kk+1]>=hubS.τP*Test[kk]+hubS.γP*(sum(u[kk+1,h]*hubS.uMax[stepI+kk,h] for h=1:H)+hubS.iD_pred[stepI+kk])^2+hubS.ρP*hubS.Tamb[stepI+kk])
	@constraint(m,Test.<=hubS.Tmax+slackT)
	@constraint(m,slackT>=0)
	#@constraint(m,currMax[h=1:H,kk=1:horzLen],sum(u[kk,h]*hubS.uMax[stepI+kk-1,h] for h=1:H)+hubS.iD_pred[stepI+kk-1]<=hubS.ItotalMax)
	@constraint(m,optOnC[nn=1:length(requiredInd)],sum(u[kk,requiredInd[nn]] for kk=1:Int(min(abs(ReqI[requiredInd[nn]]),horzLen+1)))==min(abs(ReqI[requiredInd[nn]]),horzLen+1))
	@constraint(m,optOffC[kk=1:horzLen+1],sum(u[kk,h] for h in optOffInd)==0)


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
			sol=Gurobi.get_dblattrarray(getrawsolver(m),"Xn",1,H)
		elseif solCount==1
			sol=getvalue(u)[1,:]
		else
			sol=requiredCh
		end
		Rec=sol
	else
		println(" Infeasible solution found")
		#return required EVs to charges
		Rec=requiredCh
	end

	#R>=1 and packRec get charged
	eDepart_min=hubS.eDepart_min[stepI,:]
	eArrive_actual=hubS.eArrive_actual[stepI,:]
	for h=1:H
		pemSol.U[stepI,h]= if Rec[h]==1 hubS.uMax[stepI,h] else  0  end
		prevSOC= if stepI>1 pemSol.E[stepI-1,h] else hubS.e0[h] end
		#pemSol.E[stepI,h]=prevSOC+hubS.ηP[h]*pemSol.U[stepI,h]-(eDepart_min[h]+eΔ[h])+eArrive_actual[h])
		pemSol.E[stepI,h]=prevSOC+hubS.ηP[h]*pemSol.U[stepI,h]-eDepart_min[h]+eArrive_actual[h]
	end

	return nothing
end

function pemHubstep(stepI,hubS,pemSol,silent)

	ReqI=zeros(H)
	horzLen=min(packLen,hubS.K-stepI)

	#println(stepI)
	#send requests
    for h=1:H
		ReqI[h]=pemHublocal(h,stepI,horzLen,packLen,hubS,pemSol)
	end

	pemHubcoord(stepI,horzLen,packLen,hubS,pemSol,ReqI)

	# actual dynamics
	prevT = if stepI>1 pemSol.Tactual[stepI-1] else hubS.t0 end
	pemSol.uSum[stepI] = round.(sum(pemSol.U[stepI,h] for h=1:H),sigdigits=roundSigFigs)
	pemSol.Itotal[stepI] = round.(pemSol.uSum[stepI] + hubS.iD_actual[stepI],sigdigits=roundSigFigs)
	pemSol.Tactual[stepI] = round.(hubS.τP*prevT+hubS.γP*pemSol.Itotal[stepI]^2+hubS.ρP*hubS.Tamb[stepI],sigdigits=roundSigFigs)
	return ReqI
end

function pemHub(hubS,silent::Bool)
	pemSol=hubSolutionStruct(K=hubS.K,H=H)
	Req=zeros(hubS.K,H)

	Juno.progress() do id
		for stepI=1:hubS.K
			@info "$(Dates.format(Dates.now(),"HH:MM:SS")): $(stepI) of $(hubS.K)....\n" progress=stepI/hubS.K _id=id
			@printf "%s: time step %g of %g....\n" Dates.format(Dates.now(),"HH:MM:SS") stepI hubS.K
			try
				Req[stepI,:]=pemHubstep(stepI,hubS,pemSol,silent)
			catch e
				@printf "error: %s" e
				break
			end
		end
	end

	#pemSol.objVal[1,1]=sum(sum((pemSol.Sn[stepI,n]-1)^2*hubS.Qsi[n,1]+(pemSol.Un[stepI,n])^2*hubS.Ri[n,1] for n=1:N) for stepI=1:(hubS.K))
	return pemSol, Req
end
