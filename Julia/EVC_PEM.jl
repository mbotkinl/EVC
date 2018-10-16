using Distributions
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

	#send requests
    for n=1:N

		# clean up this logic***
		if k==1
			existPack=false
		elseif (un[k-1,n]>0)
			if k<=packLen
				existPack=true
			elseif (un[k-(packLen),n]==0)
				existPack=true
			else
				existPack=false
			end
		else
			existPack=false
		end
		if existPack==true # still using a packet
			Req[k,n]=-1 # check if this should by pass the estimation step??***
			# un[k,n]=evS.imax[n]
		else
			desiredSOC=1 # for now
			setSOC=0.5 # for now
			# desiredSOC=evS.Snmin[n]
			prevSOC= if k>1 sn[k-1,n] else evS.s0[n] end
	        if (prevSOC-desiredSOC)>=-epsilon #desiredSOC satisfied
				ratio[k,n]=0
				Req[k,n]=0
			elseif k>=evS.Kn[n] # opt out (not satisified)
				ratio[k,n]=1
				Req[k,n]=-1 # check if this should by pass the estimation step??***
	        else
				ratio[k,n]=(desiredSOC-prevSOC)/(evS.ηP[n]*evS.imax[n]*(evS.Kn[n]-(k-1)))
				mu[k,n] = 1/mttr*((desiredSOC-ratio[k,n])/(ratio[k,n]-0))*((setSOC-0)/(desiredSOC-setSOC))
				P[k,n] = min(max(1-exp(-mu[k,n]*evS.Ts),0),1)
				t=rand()
			  	Req[k,n]=if (t>(1-P[k,n])) 1 else  0  end
	        end

			#decided by Temp and probability
			#un[k,n]=if (T[k]<evS.Tmax) & (t>(1-R[k,n])) evS.imax[n] else  0  end
		end
        #sn[k,n]=sn[k-1,n]+evS.ηP[n]*un[k,n]
	end

	prevT = if k>1 T[k-1] else evS.t0 end

	#receive requests and forecast for the next packLen intervals
	Test=zeros(packLen,1)
	Iest=sum(abs(Req[k,n])*evS.imax[n] for n=1:N)
	Test[1] = evS.τP*prevT+evS.γP*(Iest+evS.iD[k])^2+evS.ρP*evS.Tamb[k]
	for kk=2:packLen
		Test[kk] = evS.τP*Test[kk-1]+evS.γP*(Iest+evS.iD[kk])^2+evS.ρP*evS.Tamb[k+kk]
	end

	#for now
	if all(Test.<=evS.Tmax)
		Rec[k,:]=abs.(Req[k,:])
	else # accept the only the opt outs aka -1
		Rec[k,:]=Int.(Req[k,:].<0)
	end

	#R>=1 and packRec get charged
	for n=1:N
		un[k,n]=if Rec[k,n]==1 evS.imax[n] else  0  end
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

p1=plot(pemSol.Sn,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV SOC"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1,ymax=1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=24pt,line_width=3pt,
		minor_label_font_size=20pt,key_label_font_size=20pt))

p2=plot(pemSol.Un,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV Current (kA)"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=24pt,
		minor_label_font_size=20pt,line_width=3pt,key_label_font_size=20pt))

p3=plot(x=1:horzLen+1,y=pemSol.Xt,Geom.line,Theme(default_color=colorant"green"),
		yintercept=[evS.Tmax],Geom.hline(color=["red"],style=:dot),
		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",key_position = :top,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))

pR=plot(ratio,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("R"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=24pt,
		minor_label_font_size=20pt,line_width=3pt,key_label_font_size=20pt))


checkDesiredStates(pemSol.Sn,evS.Kn,evS.Snmin)
checkDesiredStates(pemSol.Sn,evS.Kn,ones(N))
