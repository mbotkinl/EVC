using Distributions
using Gurobi

#pull out a few key variables
N=evS.N
S=evS.S
horzLen=evS.K1

#check indexes***
sn=zeros(horzLen+1,N)
un=zeros(horzLen+1,N)
uSum=zeros(horzLen+1,1)
Itotal=zeros(horzLen+1,1)
ratio=zeros(horzLen+1,N)
mu=zeros(horzLen+1,N)
P=zeros(horzLen+1,N)
Req=zeros(horzLen+1,N)
Rec=zeros(horzLen+1,N)
T=zeros(horzLen+1,1)

# sn[1:N]=evS.s0
# T[1]=evS.t0
epsilon=1e-3
packLen=3 #number of time steps
mttr=300

tic()
for k =1:horzLen+1
	#println(k)
	#send requests

    for n=1:N
		# clean up this logic***
		if k==1
			existPack=false
		elseif (un[k-1,n]>0)
			if k<=packLen
				prevChar=k-1
				existPack=true
			elseif (un[k-(packLen),n]==0)
				prevChar=sum(un[(k-1):(k-(packLen)),n])
				existPack=true
			else
				existPack=false
			end
		else
			existPack=false
		end

		if existPack==true # still using a packet
			ratio[k,n]=1
			Req[k,n]=-packLen+prevChar
		else
			desiredSOC=1 # for now
			setSOC=0.5 # for now
			# desiredSOC=evS.Snmin[n] #this should be a ratio?
			prevSOC= if k>1 sn[k-1,n] else evS.s0[n] end
	        if (prevSOC-desiredSOC)>=-epsilon #desiredSOC satisfied
				ratio[k,n]=0
				Req[k,n]=0
				#elseif k>=evS.Kn[n] # opt out (not satisified)
			elseif ratio[k,n]>=1 # opt out (need to charge for rest of time)
				Req[k,n]=-packLen
	        else
				ratio[k,n]=(desiredSOC-prevSOC)/(evS.ηP[n]*evS.imax[n]*(evS.Kn[n]-(k-1)))
				#mu[k,n] = 1/mttr*((desiredSOC-ratio[k,n])/(ratio[k,n]-0))*((setSOC-0)/(desiredSOC-setSOC))
				mu[k,n] = 1/mttr*((ratio[k,n]-0)/(desiredSOC-ratio[k,n]))*((desiredSOC-setSOC)/(setSOC-0))
				#mu[k,n] = ratio[k,n]
				P[k,n] = min(max(1-exp(-mu[k,n]*evS.Ts),0),1)
				t=rand()
			  	Req[k,n]=if (t>(1-P[k,n])) 1 else  0  end
	        end
		end
	end

	prevT = if k>1 T[k-1] else evS.t0 end
    requiredCh=Int.(Req[k,:].<0)
	#requiredCh=Int.(ratio[k,:].>=1)


	#receive requests and forecast for the next packLen intervals

	# Test=zeros(packLen,1)
	# Iest=sum(abs(Req[k,n])*evS.imax[n] for n=1:N)
	# Test[1] = evS.τP*prevT+evS.γP*(Iest+evS.iD[k])^2+evS.ρP*evS.Tamb[k]
	# for kk=2:packLen
	# 	Test[kk] = evS.τP*Test[kk-1]+evS.γP*(Iest+evS.iD[kk])^2+evS.ρP*evS.Tamb[k+kk]
	# end
	#for now
	# if all(Test.<=evS.Tmax)
	# 	Rec[k,:]=abs.(Req[k,:])
	# else # accept the only the opt outs aka -1
	# 	Rec[k,:]=Int.(Req[k,:].<0)
	# end

	requiredInd=find(x->x==1,requiredCh)
	optOff=find(x->x==0,Req[k,:])

	m = Model(solver = GurobiSolver())
	@variable(m,u[1:packLen,1:N],Bin)
	@variable(m,Test[1:packLen])
	@objective(m,Min,-sum(u))
	@constraint(m,tempCon1,Test[1]>=evS.τP*prevT+evS.γP*(sum(u[1,n]*evS.imax[n] for n=1:N)+evS.iD[k])^2+evS.ρP*evS.Tamb[k])
	@constraint(m,tempCon2[kk=2:packLen],Test[kk]>=evS.τP*Test[kk-1]+evS.γP*(sum(u[kk,n]*evS.imax[n] for n=1:N)+evS.iD[k+kk])^2+evS.ρP*evS.Tamb[k+kk])
	@constraint(m,Test.<=evS.Tmax)
	# @constraint(m,optOn[kk=1:packLen],sum(u[n,kk] for n in requiredInd[kk])==length(requiredInd[kk]))
	@constraint(m,optOn[nn=1:length(requiredInd)],sum(u[kk,requiredInd[nn]] for kk=1:Int(abs(Req[k,requiredInd[nn]])))==abs(Req[k,requiredInd[nn]]))
	@constraint(m,optOff[kk=1:packLen],sum(u[kk,n] for n in optOff)==0)

	TT = STDOUT
	redirect_stdout()
	statusC = solve(m)
	redirect_stdout(TT)

	# m = Model(solver = GurobiSolver())
	# @variable(m,u[1:5],Bin)
	# @objective(m,Min,-sum(u))
	# @constraint(m,sum(u)==3)
	# statusC = solve(m)
    # getvalue(u)
	#


	if statusC==:Optimal
		Rec[k,:]=getvalue(u)[1,:]
	else
		#return
		Rec[k,:]=requiredCh
	end

	#R>=1 and packRec get charged
	for n=1:N
		un[k,n]= if Rec[k,n]==1 evS.imax[n] else  0  end
		prevSOC= if k>1 sn[k-1,n] else evS.s0[n] end
		sn[k,n]=prevSOC+evS.ηP[n]*un[k,n]
	end

	# actual dynamics
	uSum[k] = sum(un[k,n] for n=1:N)
	Itotal[k] = uSum[k] + evS.iD[k]
	T[k] = evS.τP*prevT+evS.γP*Itotal[k]^2+evS.ρP*evS.Tamb[k]
end
timeT=toc()
timeT=timeT/horzLen

objVal=sum(sum((sn[k,n]-1)^2*evS.Qsi[n,1]+(un[k,n])^2*evS.Ri[n,1] for n=1:N) for k=1:(horzLen+1))

pemSol=centralSolutionStruct(Xt=T,Un=un,Sn=sn, Itotal=Itotal,uSum=uSum,objVal=objVal)

filename = "d_PEM_N$(N)"
# save
if saveResults==1 saveRun(path,filename,timeT, evS,pemSol) end
# load
if loadResults==1
	loadF=JLD.load(path*filename*".jld")
	evS=loadF["scenario"]
	pemSol=loadF["solution"]
end

p1pem=plot(pemSol.Sn,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV SOC"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1,ymax=1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=24pt,line_width=3pt,
		minor_label_font_size=20pt,key_label_font_size=20pt))

p2pem=plot(pemSol.Un,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV Current (kA)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=24pt,
		minor_label_font_size=20pt,line_width=3pt,key_label_font_size=20pt))

p3pem=plot(x=1:horzLen+1,y=pemSol.Xt,Geom.line,Theme(default_color=colorant"green"),
		yintercept=[evS.Tmax],Geom.hline(color=["red"],style=:dot),
		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",key_position = :top,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))

pRpem=plot(ratio,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("R"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=24pt,
		minor_label_font_size=20pt,line_width=3pt,key_label_font_size=20pt))


checkDesiredStates(pemSol.Sn,evS.Kn,evS.Snmin)
checkDesiredStates(pemSol.Sn,evS.Kn,ones(N))
