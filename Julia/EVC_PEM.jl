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
R=zeros(horzLen,N)
T=zeros(horzLen+1,1)

sn[1:N]=evS.s0
T[1]=evS.t0
epsilon=1e-3

tic()
for k =1:horzLen

    for n =1:N
		desiredSOC=1 # for now
		# desiredSOC=evS.Snmin[n]
        if (sn[k,n]-desiredSOC)>=-epsilon #desiredSOC satisfied
			R[k,n]=0 # not sure if this is what we want***
		elseif k>=evS.Kn[n] # opt out (not satisified)
			R[k,n]=1
        else
			R[k,n]=(desiredSOC-sn[k,n])/(evS.ηP[n]*evS.imax[n]*(evS.Kn[n]-k))
        end
        t=rand()
		#decided by Temp and probability
		un[k+1,n]=if (T[k]<evS.Tmax) & (t>(1-R[k,n])) evS.imax[n] else  0  end
        sn[k+1,n]=sn[k,n]+evS.ηP[n]*un[k+1,n]
	end

	uSum[k+1] = sum(un[k+1,n] for n=1:N)
	Itotal[k+1] = uSum[k+1] + evS.w[(k-1)*2+1] 	#add background current***
	T[k+1] = evS.τP*T[k]+evS.γP*Itotal[k+1]^2+evS.ρP*evS.w[k*2+2,1]
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

pR=plot(R,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("R"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=24pt,
		minor_label_font_size=20pt,line_width=3pt,key_label_font_size=20pt))


checkDesiredStates(pemSol.Sn,evS.Kn,evS.Snmin)
checkDesiredStates(pemSol.Sn,evS.Kn,ones(N))
