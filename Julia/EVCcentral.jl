#Micah Botkin-Levy
#4/8/18

tic()
#initialize with current states
sn0=s0
xt0=T0

stepI=1
horzLen=K1

#desired SOC
target=zeros(N*(horzLen+1),1);
for ii=1:N
   cur=Kn[ii]-(stepI-1)
   ind=max(0,(cur-1)*N)+ii:N:length(target)
   target[ind]=Snmin[ii,1]
end


println("setting up model")
centralModel = Model(solver = GurobiSolver(Presolve=0,BarHomogeneous=1,NumericFocus=3))
#m = Model(solver = IpoptSolver())

#u w and z are one index ahead of x. i.e the x[k+1]=x[k]+eta*u[k+1]
@variable(centralModel,u[1:N*(horzLen+1)])
@variable(centralModel,sn[1:(N)*(horzLen+1)])
@variable(centralModel,xt[1:(horzLen+1)])
@variable(centralModel,z[1:S*(horzLen+1)])

println("obj")
objFun(sn,xt,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*Qsi[n,1]     for n=1:N) for k=1:horzLen+1) +
				sum((xt[k,1]-1)^2*Qsi[N+1,1]                 for k=1:horzLen+1) +
				sum(sum((u[(k-1)*N+n,1])^2*Ri[n,1]           for n=1:N) for k=1:horzLen+1)
@objective(centralModel,Min, objFun(sn,xt,u))

println("constraints")
@constraint(centralModel,stateCon1,sn[1:N,1].==sn0[1:N,1]+ηP[:,1].*u[1:N,1])
@constraint(centralModel,stateCon2[k=1:horzLen,n=1:N],sn[n+(k)*(N),1]==sn[n+(k-1)*(N),1]+ηP[n,1]*u[n+(k)*(N),1])
@constraint(centralModel,tempCon1,xt[1,1]==τP*xt0+γP*deltaI*sum((2*s-1)*z[s,1] for s=1:S)+ρP*w[stepI*2,1])
@constraint(centralModel,tempCon2[k=1:horzLen],xt[k+1,1]==τP*xt[k,1]+γP*deltaI*sum((2*s-1)*z[k*S+s,1] for s=1:S)+ρP*w[stepI*2+k*2,1])
@constraint(centralModel,currCon[k=1:horzLen+1],0==-sum(u[(k-1)*(N)+n] for n=1:N)-w[(k-1)*2+1]+sum(z[(k-1)*(S)+s] for s=1:S))
@constraint(centralModel,sn.<=1)
@constraint(centralModel,sn.>=target)
if noTlimit==0
	@constraint(centralModel,upperTCon,xt.<=Tmax)
end
@constraint(centralModel,xt.>=0)
@constraint(centralModel,upperCCon,u.<=repmat(imax,horzLen+1,1))
@constraint(centralModel,u.>=repmat(imin,horzLen+1,1))
@constraint(centralModel,z.>=0)
@constraint(centralModel,z.<=deltaI)

println("solving....")
statusM = solve(centralModel)
if statusM!=:Optimal
	@printf "Failed %s \n" statusM
    return
else
	toc()

	uRaw=getvalue(u)
	snRaw=getvalue(sn)
	xtRaw=getvalue(xt)
	zRaw=getvalue(z)
	f=objFun(snRaw,xtRaw,uRaw)


	#calculate actual temp
	Tactual=zeros(horzLen+1,1)
	ztotal=zeros(horzLen+1,1)
	for k=1:horzLen+1
		ztotal[k,1]=sum((uRaw[(k-1)*N+n,1]) for n=1:N) + w[(k-1)*2+1,1]
	end
	Tactual[1,1]=τP*T0+γP*ztotal[1,1]^2+ρP*w[2,1]
	for k=1:horzLen
		Tactual[k+1,1]=τP*Tactual[k,1]+γP*ztotal[k+1,1]^2+ρP*w[k*2+2,1]
	end

	#getobjectivevalue(m)
	lambdaCurr=-getdual(currCon)
	#lambdaTemp=[getdual(tempCon1);getdual(tempCon2)]
	#lambdaState=[getdual(stateCon1)';getdual(stateCon2)]
	if noTlimit==0
		lambdaUpperT=-getdual(upperTCon)
	else
		lambdaUpperT=zeros(horzLen+1,1)
	end

	#lambdaUpperC=-getdual(upperCCon)

	xtStar=xtRaw
    snStar=snRaw
	uStar=uRaw
	zStar=zRaw
    fStar=getobjectivevalue(centralModel)
    lamTempStar=lambdaUpperT
	lamCurrStar=lambdaCurr

    uSumStar=zeros(horzLen+1,1)
	zSumStar=zeros(horzLen+1,1)
	for k=1:horzLen+1
		zSumStar[k,1]=sum(zStar[(k-1)*(S)+s,1] for s=1:S)
        uSumStar[k,1]=sum(uStar[(k-1)*N+n,1] for n=1:N)
	end

	println("plotting....")
	xPlot=zeros(horzLen+1,N)
	xPlot2=zeros(horzLen+1,N)
	for ii= 1:N
		xPlot[:,ii]=snRaw[collect(ii:N:length(snRaw))]
		xPlot2[:,ii]=(Snmin[ii,1]-xPlot[:,ii])./(Kn[ii,1]-(1:1:length(xPlot[:,ii])))
	end

	#plot(x=1:horzLen+1,y=xPlot2[:,ii])
	# plot(x=1:Kn[ii,1],y=xPlot2[1:Kn[ii,1],ii])

	p1=plot(xPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel("PEV SOC"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),
			Theme(background_color=colorant"white",key_position = :none,major_label_font_size=24pt,line_width=3pt,
			minor_label_font_size=20pt,key_label_font_size=20pt))
	p1b=plot(xPlot2,x=Row.index,y=Col.value,color=Col.index,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel("PEV SOC"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),
			Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
			minor_label_font_size=16pt,key_label_font_size=16pt))
	if drawFig==1 draw(PNG(path*"J_central_SOC.png", 24inch, 12inch), p1) end


	uPlot=zeros(horzLen+1,N)
	for ii= 1:N
		uPlot[:,ii]=uRaw[collect(ii:N:length(uRaw))]
	end

	p2=plot(uPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel("PEV Current (kA)"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),
			Theme(background_color=colorant"white",key_position = :none,major_label_font_size=24pt,
			minor_label_font_size=20pt,line_width=3pt,key_label_font_size=20pt))
	if drawFig==1 draw(PNG(path*"J_central_Curr.png", 24inch, 12inch), p2) end

	p3=plot(layer(x=1:horzLen+1,y=xtRaw,Geom.line,Theme(default_color=colorant"blue")),
			layer(x=1:horzLen+1,y=Tactual,Geom.line,Theme(default_color=colorant"green")),
			yintercept=[Tmax],Geom.hline(color=["red"],style=:dot),
			Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)",orientation=:vertical),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",key_position = :top,major_label_font_size=18pt,
			minor_label_font_size=16pt,key_label_font_size=16pt),
			Guide.manual_color_key("", ["PWL Temp", "Actual Temp"], ["blue", "green"]))
	if drawFig==1 draw(PNG(path*"J_central_Temp.png", 24inch, 12inch), p3) end

	p4b=plot(x=1:horzLen+1,y=lambdaUpperT,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel(raw"Lambda ($/K)",orientation=:vertical),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=18pt,
			minor_label_font_size=16pt,key_label_font_size=16pt))
	p4=plot(x=1:horzLen+1,y=lambdaCurr	,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel(raw"Lambda ($/kA)",orientation=:vertical),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=18pt,
			minor_label_font_size=16pt,key_label_font_size=16pt))
	if drawFig==1 draw(PNG(path*"J_central_Lam.png", 24inch, 12inch), p4) end

	fName="J_Central.png"

	#draw(PNG(path*fName, 13inch, 14inch), vstack(p1,p2,p3,p4))



	#uPlotNoLim=uPlot
	#do this more elegantly
	# aggU=plot(layer(x=1:horzLen+1,y=sum(uPlot[:,i] for i=1:N),Geom.line,Theme(default_color=colorant"blue")),
	# 		layer(x=1:horzLen+1,y=sum(uPlotNoLim[:,i] for i=1:N),Geom.line,Theme(default_color=colorant"red")),
	# 		Guide.xlabel("Time"), Guide.ylabel("PEV Current (kA)"),
	# 		Coord.Cartesian(xmin=0,xmax=horzLen+1),
	# 		Theme(background_color=colorant"white"),
	# 		Guide.manual_color_key("", ["Aggregate Current", "Unconstrained Aggregate Current"], ["blue","red"]))
	# fName="aggCentral_noLim.png"
	#  draw(PNG(path*fName, 13inch, 14inch), aggU)


	# p3nolimit=plot(layer(x=1:horzLen+1,y=xtRaw,Geom.line,Theme(default_color=colorant"blue")),
	# 		layer(x=1:horzLen+1,y=Tactual,Geom.line,Theme(default_color=colorant"green")),
	# 		layer(x=1:horzLen+1,y=noLimitt,Geom.line,Theme(default_color=colorant"red")),
	# 		yintercept=[Tmax],Geom.hline(color=["red"],style=:dot),
	# 		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)",orientation=:vertical),
	# 		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",key_position = :top),
	# 		Guide.manual_color_key("", ["PWL Temp", "Actual Temp","Unconstrained Temp"], ["blue", "green","red"]))
	# fName="Temp_w_nolimit.png"
	# draw(PNG(path*fName, 13inch, 14inch), p3nolimit)

end
