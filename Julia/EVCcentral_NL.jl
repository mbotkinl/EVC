#Micah Botkin-Levy
#4/8/18

tic()
#initialize with current states
sn0=s0
xt0=T0

#add mpc loop here ??
stepI=1
horzLen=K1

println("setting up model")
m = Model(solver = IpoptSolver())
#m = Model(solver = IpoptSolver(linear_solver = "ma86"))

#u w and z are one index ahead of x. i.e the x[k+1]=x[k]+eta*u[k+1]
@variable(m,u[1:N*(horzLen+1)])
@variable(m,sn[1:(N)*(horzLen+1)])
@variable(m,xt[1:(horzLen+1)])
@variable(m,itotal[1:(horzLen+1)])

#desired SOC
target=zeros(N*(horzLen+1),1);
for ii=1:N
   cur=Kn[ii]-(stepI-1)
   ind=max(0,(cur-1)*N)+ii:N:length(target)
   target[ind]=Snmin[ii,1]
end

println("obj")
objFun(sn,xt,u)=sum(sum((sn[(k-1)*(N)+n,1]-1)^2*Qsi[n,1]     for n=1:N) for k=1:horzLen+1) +
				sum((xt[k,1]-1)^2*Qsi[N+1,1]                 for k=1:horzLen+1) +
				sum(sum((u[(k-1)*N+n,1])^2*Ri[n,1]           for n=1:N) for k=1:horzLen+1)
@objective(m,Min, objFun(sn,xt,u))

println("constraints")
@constraint(m,stateCon1,sn[1:N,1].==sn0[1:N,1]+etaP[:,1].*u[1:N,1])
@constraint(m,stateCon2[k=1:horzLen,n=1:N],sn[n+(k)*(N),1]==sn[n+(k-1)*(N),1]+etaP[n,1]*u[n+(k)*(N),1])
@constraint(m,tempCon1,xt[1,1]==tauP*xt0+gammaP*(itotal[1])^2+rhoP*w[stepI*2,1])
@constraint(m,tempCon2[k=1:horzLen],xt[k+1,1]==tauP*xt[k,1]+gammaP*(itotal[k+1])^2+rhoP*w[stepI*2+k*2,1]) #check id index???
@constraint(m,currCon[k=1:horzLen+1],0==-sum(u[(k-1)*(N)+n,1] for n=1:N)-w[(k-1)*2+1]+itotal[k]) #fix for MPC
@constraint(m,sn.<=1)
@constraint(m,sn.>=target)
if noTlimit==0
	@constraint(m,upperTCon,xt.<=Tmax)
end
@constraint(m,xt.>=0)
@constraint(m,upperCCon,u.<=repmat(imax,horzLen+1,1))
@constraint(m,u.>=repmat(imin,horzLen+1,1))
@constraint(m,itotal.<=ItotalMax)
@constraint(m,itotal.>=0)

println("solving....")
status = solve(m)
if status!=:Optimal
	@printf "Failed %s \n" status
    return
else
	toc()

	uRaw=getvalue(u)
	snRaw=getvalue(sn)
	xtRaw=getvalue(xt)
    itotalRaw=getvalue(itotal)

	if noTlimit==0
		lambdaUpperT=-getdual(upperTCon)
	else
		lambdaUpperT=zeros(horzLen+1,1)
	end
    lambdaCurr=-getdual(currCon)

    xtStarNL=xtRaw
    snStarNL=snRaw
    uStarNL=uRaw
    itotalStarNL=itotalRaw
    fStarNL=getobjectivevalue(m)
    lamTempStarNL=lambdaUpperT
    lamCurrStarNL=lambdaCurr

    uSumStarNL=zeros(horzLen+1,1)
    for k=1:horzLen+1
        uSumStarNL[k,1]=sum(uStarNL[(k-1)*N+n,1] for n=1:N)
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

	p1nl=plot(xPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel("PEV SOC"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),
			Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
			minor_label_font_size=16pt,key_label_font_size=16pt))
	p1bnl=plot(xPlot2,x=Row.index,y=Col.value,color=Col.index,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel("PEV SOC"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),
			Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
			minor_label_font_size=16pt,key_label_font_size=16pt))
	if drawFig==1 draw(PNG(path*"J_centralNL_SOC.png", 24inch, 12inch), p1nl) end


	uPlot=zeros(horzLen+1,N)
	for ii= 1:N
		uPlot[:,ii]=uRaw[collect(ii:N:length(uRaw))]
	end

	p2nl=plot(uPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel("PEV Current (A)"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),
			Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
			minor_label_font_size=16pt,key_label_font_size=16pt))
	if drawFig==1 draw(PNG(path*"J_centralNL_Curr.png", 24inch, 12inch), p2nl) end

	p3nl=plot(layer(x=1:horzLen+1,y=xtRaw,Geom.line,Theme(default_color=colorant"blue")),
			yintercept=[Tmax],Geom.hline(color=["red"],style=:dot),
			Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)",orientation=:vertical),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",key_position = :top,major_label_font_size=18pt,
			minor_label_font_size=16pt,key_label_font_size=16pt))
	if drawFig==1 draw(PNG(path*"J_centralNL_Temp.png", 24inch, 12inch), p3nl) end

	p4nlb=plot(x=1:horzLen+1,y=lambdaUpperT,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel(raw"Lambda ($/K)",orientation=:vertical),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=18pt,
			minor_label_font_size=16pt,key_label_font_size=16pt))
    p4nl=plot(x=1:horzLen+1,y=lambdaCurr,Geom.line,
            Guide.xlabel("Time"), Guide.ylabel(raw"Lambda ($/A)",orientation=:vertical),
            Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=18pt,
            minor_label_font_size=16pt,key_label_font_size=16pt))
	if drawFig==1 draw(PNG(path*"J_centralNL_Lam.png", 24inch, 12inch), p4nl) end

	fName="J_Central.png"

	#draw(PNG(path*fName, 13inch, 14inch), vstack(p1,p2,p3,p4))



	#uPlotNoLim=uPlot
	#do this more elegantly
	# aggU=plot(layer(x=1:horzLen+1,y=sum(uPlot[:,i] for i=1:N),Geom.line,Theme(default_color=colorant"blue")),
	# 		layer(x=1:horzLen+1,y=sum(uPlotNoLim[:,i] for i=1:N),Geom.line,Theme(default_color=colorant"red")),
	# 		Guide.xlabel("Time"), Guide.ylabel("PEV Current (A)"),
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
