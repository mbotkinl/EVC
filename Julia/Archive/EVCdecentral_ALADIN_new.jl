#EVC with ALADIN for PWL convex relaxation
#Micah Botkin-Levy
#5/22/18

#formulation try 2
#u, sn, xt, and z are all in "y", v is the auxilliary variable corresponding to x in literature
#current constraint is coupling

tic()

#compile functions
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCdecentral_ALADIN.jl")

#initialize with current states
sn0=s0
xt0=T0

stepI = 1;
horzLen=K1
epsilon = 1e-8
tolU=1e-4
tolS=1e-8
tolT=1e-4
tolZ=1e-6
maxIt=30
convIt=maxIt
ConvALAD=zeros(maxIt,1)
constConvALAD=zeros(maxIt,1)
itConvALAD=zeros(maxIt,1)
fConvALAD=zeros(maxIt,1)
snConvALAD=zeros(maxIt,1)
unConvALAD=zeros(maxIt,1)
convCheck=zeros(maxIt,1)

#ALADIN tuning and initial guess
σU=1*ones(N,1)
σS=ones(N,1)/10 #for kA
#σS=ones(N,1)*100 #for A
σZ=1/N
σT=1/10000 #for kA
#σT=1/10  #for A

Hz=1e-6
Ht=1e-6
ρALAD=1
ρRate=1.15
muALAD=10^8
# lambda0=5*rand(Truncated(Normal(0), 0, 1), horzLen+1)
# vt0=Tmax*rand(Truncated(Normal(0), 0, 1), (horzLen+1))
# vz0=ItotalMax/1000*rand(Truncated(Normal(0), 0, 1), S*(horzLen+1))
# vu0=imax[1,1]*0.8*rand(Truncated(Normal(0), 0, 1), N*(horzLen+1))
# vs0=rand(Truncated(Normal(0), 0, 1), N*(horzLen+1))

lambda0=1000*ones(horzLen+1,1)
vt0=ones(horzLen+1,1)
vz0=ones(S*(horzLen+1),1)
vu0=.01*ones(N*(horzLen+1),1)
vs0=.5*ones(N*(horzLen+1),1)

# lambda0=lamCurrStar+
# vt0=xtStar
# vz0=zStar
# vu0=uStar
# vs0=snStar

#save matrices
Un=SharedArray{Float64}(N*(horzLen+1),maxIt) #row are time (N states for k=1, them N states for k=2),  columns are iteration
Sn=SharedArray{Float64}(N*(horzLen+1),maxIt)  #row are time,  columns are iteration
Xt=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration
Z=zeros(S*(horzLen+1),maxIt)  #row are time,  columns are iteration
Lam=zeros((horzLen+1),maxIt) #(rows are time, columns are iteration)
#Lam[:,1]=max.(lambda0,0)
Lam[:,1]=lambda0

ρALADp=ρALAD*ones(1,maxIt)

#auxillary variables
Vu=zeros((N)*(horzLen+1),maxIt) #row are time,  columns are iteration
Vu[:,1]=vu0 #initial guess goes in column 1
Vs=zeros((N)*(horzLen+1),maxIt) #row are time,  columns are iteration
Vs[:,1]=vs0 #initial guess goes in column 1
Vz=zeros(S*(horzLen+1),maxIt) #row are time,  columns are iteration
Vz[:,1]=vz0
Vt=zeros((horzLen+1),maxIt) #row are time,  columns are iteration
Vt[:,1]=vt0

#Gradian Vectors
Gu=SharedArray{Float64}(N*(horzLen+1),maxIt) #row are time (N states for k=1, them N states for k=2),  columns are iteration
Gs=SharedArray{Float64}(N*(horzLen+1),maxIt) #row are time (N states for k=1, them N states for k=2),  columns are iteration
Gz=zeros(S*(horzLen+1),maxIt) #row are time (N states for k=1, them N states for k=2),  columns are iteration
Gt=zeros(S*(horzLen+1),maxIt) #row are time (N states for k=1, them N states for k=2),  columns are iteration

#Jacobian C Vectors
Cs=SharedArray{Float64}(N*(horzLen+1),maxIt)  #row are time,  columns are iteration
Cu=SharedArray{Float64}(N*(horzLen+1),maxIt)  #row are time,  columns are iteration
Cz=zeros(S*(horzLen+1),maxIt)  #row are time,  columns are iteration
Ct=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration

#central QP solutions
save_dUn=zeros(N*(horzLen+1),maxIt)  #row are time,  columns are iteration
save_dXt=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration
save_dSn=zeros(N*(horzLen+1),maxIt)  #row are time,  columns are iteration
save_dZ=zeros(S*(horzLen+1),maxIt)  #row are time,  columns are iteration
lamQP=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration

uSum=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration
zSum=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration
currConst=zeros((horzLen+1),maxIt)  #row are time,  columns are iteration
deltaY=zeros(1,maxIt)



#initialize optimization problems
println("Setting up optimization problems")
#EV nlp
evM = Model(solver = IpoptSolver())
@variable(evM,sn[1:(horzLen+1)])
@variable(evM,u[1:(horzLen+1)])

# @variable(evM,Qn==1)
# @variable(evM,Rn==1)
# @variable(evM,σUn==1)
# @variable(evM,σSn==1)
# @variable(evM,ηPn==1)
# @variable(evM,sn0n==0)
# @variable(evM,imaxn==1)
# @variable(evM,iminn==0)
# @variable(evM,evVu[1:(horzLen+1)]==1)
# @variable(evM,evVs[1:(horzLen+1)]==1)
# @variable(evM,lambdaEV[1:(horzLen+1)]==1)
# @variable(evM,targetEV[1:(horzLen+1)]==1)
# @variable(evM,ρALADev==1)

@NLparameter(evM,Qn==1)
@NLparameter(evM,Rn==1)
@NLparameter(evM,σUn==1)
@NLparameter(evM,σSn==1)
@NLparameter(evM,ηPn==1)
@NLparameter(evM,sn0n==0)
@NLparameter(evM,imaxn==1)
@NLparameter(evM,iminn==0)
@NLparameter(evM,evVu[1:(horzLen+1)]==1)
@NLparameter(evM,evVs[1:(horzLen+1)]==1)
@NLparameter(evM,lambdaEV[1:(horzLen+1)]==1)
@NLparameter(evM,targetEV[1:(horzLen+1)]==1)
@NLparameter(evM,ρALADev==1)

@NLobjective(evM,Min,sum((sn[k,1]-1)^2*Qn+(u[k,1])^2*Rn+
                        lambdaEV[k,1]*(u[k,1])+
                        ρALADev/2*(u[k,1]-evVu[k,1])*σUn*(u[k,1]-evVu[k,1])+
                        ρALADev/2*(sn[k,1]-evVs[k,1])*σSn*(sn[k,1]-evVs[k,1]) for k=1:horzLen+1))
@NLconstraint(evM,sn[1,1]==sn0n+ηPn*u[1,1])
@NLconstraint(evM,[k=1:horzLen],sn[k+1,1]==sn[k,1]+ηPn*u[k+1,1])
@constraint(evM,socKappaMax,sn.<=1)
@NLconstraint(evM,socKappaMin[k=1:horzLen],sn[k]>=targetEV[k])
@NLconstraint(evM,curKappaMax[k=1:horzLen],u[k]<=imaxn)
@NLconstraint(evM,curKappaMin[k=1:horzLen],u[k]>=iminn)
JuMP.build(evM)
# TT = STDOUT # save original STDOUT stream
# redirect_stdout()
# statusEV = solve(evM)
# redirect_stdout(TT)


#XFRM nlp
tM = Model(solver = IpoptSolver())
#optimization variables
@variable(tM,z[1:(S)*(horzLen+1)])
@variable(tM,xt[1:(horzLen+1)])
#input variables
#@variable(tM,lambdaXFRM[1:(horzLen+1)]==1)
#@variable(tM,tempVz[1:S*(horzLen+1)]==1)
#@variable(tM,tempVt[1:(horzLen+1)]==1)
#@variable(tM,ρALADXFRM==1)   #fix this????
@NLparameter(tM,lambdaXFRM[1:(horzLen+1)]==1)
@NLparameter(tM,tempVz[1:S*(horzLen+1)]==1)
@NLparameter(tM,tempVt[1:(horzLen+1)]==1)
@NLparameter(tM,ρALADXFRM==1)

@NLobjective(tM,Min, sum(-lambdaXFRM[k,1]*sum(z[(k-1)*(S)+s] for s=1:S)+
                    ρALADXFRM/2*σZ*(sum(z[(k-1)*(S)+s] for s=1:S)-sum(tempVz[(k-1)*(S)+s] for s=1:S))^2+
                    ρALADXFRM/2*σT*(xt[k]-tempVt[k,1])^2  for k=1:(horzLen+1)))
@constraint(tM,tempCon1,xt[1]-τP*xt0-γP*deltaI*sum((2*s-1)*z[s] for s=1:S)-ρP*w[stepI*2,1]==0) #need to add w to variables for MPC???
@constraint(tM,tempCon2[k=1:horzLen],xt[k+1]-τP*xt[k]-γP*deltaI*sum((2*s-1)*z[(k)*(S)+s] for s=1:S)-ρP*w[stepI*2+k*2,1]==0)
if noTlimit==0
    @constraint(tM,upperTCon,xt.<=Tmax)
end
@constraint(tM,lowerTCon,xt.>=0)
@constraint(tM,pwlKappaMin,z.>=0)
@constraint(tM,pwlKappaMax,z.<=deltaI)
TT = STDOUT # save original STDOUT stream
redirect_stdout()
statusTM = solve(tM)
redirect_stdout(TT)

#xfrmNLPmodel,lambda,tempVz,tempVt,ρALADp=xfrmModel_init()

#Central QP
cM = Model(solver = IpoptSolver())
@variable(cM,dUn[1:(N)*(horzLen+1)])
@variable(cM,dSn[1:(N)*(horzLen+1)])
@variable(cM,dZ[1:(S)*(horzLen+1)])
@variable(cM,dXt[1:(horzLen+1)])

# @variable(cM,tempGu[1:N*(horzLen+1)]==1)
# @variable(cM,tempGs[1:N*(horzLen+1)]==1)
# @variable(cM,tempZ[1:S*(horzLen+1)]==1)
# @variable(cM,tempUn[1:N*(horzLen+1)]==1)
# @variable(cM,tempCz[1:S*(horzLen+1)]==1)
# @variable(cM,tempCs[1:N*(horzLen+1)]==1)
# @variable(cM,tempCu[1:N*(horzLen+1)]==1)
# @variable(cM,tempCt[1:(horzLen+1)]==1)

@NLparameter(cM,tempGu[1:N*(horzLen+1)]==1)
@NLparameter(cM,tempGs[1:N*(horzLen+1)]==1)
@NLparameter(cM,tempZ[1:S*(horzLen+1)]==1)
@NLparameter(cM,tempUn[1:N*(horzLen+1)]==1)
@NLparameter(cM,tempCz[1:S*(horzLen+1)]==1)
@NLparameter(cM,tempCs[1:N*(horzLen+1)]==1)
@NLparameter(cM,tempCu[1:N*(horzLen+1)]==1)
@NLparameter(cM,tempCt[1:(horzLen+1)]==1)

# objExp=sum(sum(0.5*dUn[(k-1)*N+i,1]^2*2*Ri[i,1]+tempGu[(k-1)*N+i,1]*dUn[(k-1)*N+i,1]+
#                0.5*dSn[(k-1)*N+i,1]^2*2*Qsi[i,1]+tempGs[(k-1)*N+i,1]*dSn[(k-1)*N+i,1] for i=1:N) for k=1:(horzLen+1))
# objExp=objExp+sum(0.5*dZ[k,1]^2*Hz+
#                   0.5*dXt[k,1]^2*Ht   for k=1:(horzLen+1))
# @NLobjective(cM,Min, objExp)
@NLobjective(cM,Min, sum(sum(0.5*dUn[(k-1)*N+i,1]^2*2*Ri[i,1]+tempGu[(k-1)*N+i,1]*dUn[(k-1)*N+i,1]+
               0.5*dSn[(k-1)*N+i,1]^2*2*Qsi[i,1]+tempGs[(k-1)*N+i,1]*dSn[(k-1)*N+i,1] for i=1:N) +
               0.5*dZ[k,1]^2*Hz+ 0.5*dXt[k,1]^2*Ht   for k=1:(horzLen+1)))

@NLconstraint(cM,currCon[k=1:horzLen+1],sum(tempUn[(k-1)*(N)+n,1]+dUn[(k-1)*(N)+n,1] for n=1:N)-
                                         sum(tempZ[(k-1)*(S)+s,1]+dZ[(k-1)*(S)+s,1] for s=1:S)==-w[(k-1)*2+1])#+relaxS[k,1])
@constraint(cM,stateCon1[n=1:N],dSn[n,1]==ηP[n,1]*dUn[n,1])
@constraint(cM,stateCon2[k=1:horzLen,n=1:N],dSn[n+(k)*(N),1]==dSn[n+(k-1)*(N),1]+ηP[n,1]*dUn[n+(k)*(N),1])
@constraint(cM,tempCon1,dXt[1,1]==γP*deltaI*sum((2*s-1)*dZ[s,1] for s=1:S))
@constraint(cM,tempCon2[k=1:horzLen],dXt[k+1,1]==τP*dXt[k,1]+γP*deltaI*sum((2*s-1)*dZ[k*S+s,1] for s=1:S))
@NLconstraint(cM,[ii=1:length(tempCz)],tempCz[ii,1]*dZ[ii]<=0)
@NLconstraint(cM,[ii=1:length(tempCu)],tempCu[ii,1]*dUn[ii]<=0)
@NLconstraint(cM,[ii=1:length(tempCs)],tempCs[ii,1]*dSn[ii]<=0)
@NLconstraint(cM,[k=1:horzLen+1],tempCt[k,1]*dXt[k]<=0)

TT = STDOUT # save original STDOUT stream
redirect_stdout()
statusCM = solve(cM)
redirect_stdout(TT)

for p=1:maxIt-1
    runALADit(p)
end

println("plotting....")
xPlot=zeros(horzLen+1,N)
uPlot=zeros(horzLen+1,N)
for ii= 1:N
	xPlot[:,ii]=Sn[collect(ii:N:length(Sn[:,convIt])),convIt]
    uPlot[:,ii]=Un[collect(ii:N:length(Un[:,convIt])),convIt]
end

pd1alad=plot(xPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV SOC"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_decentral_ALADIN_SOC.png", 24inch, 12inch), pd1alad) end

pd2alad=plot(uPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV Current (kA)"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_decentral_ALADIN_Curr.png", 24inch, 12inch), pd2alad) end

pd3alad=plot(layer(x=1:horzLen+1,y=Xt[:,convIt],Geom.line,Theme(default_color=colorant"blue")),
		#layer(x=1:horzLen+1,y=Tactual,Geom.line,Theme(default_color=colorant"green")),
		yintercept=[Tmax],Geom.hline(color=["red"],style=:dot),
		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",key_position = :top,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
		#Guide.manual_color_key("", ["PWL Temp", "Actual Temp"], ["blue", "green"]))
if drawFig==1 draw(PNG(path*"J_decentral_ALADIN_Temp.png", 24inch, 12inch), pd3alad) end

pd4alad=plot(layer(x=1:horzLen+1,y=Lam[:,convIt],Geom.line),
		layer(x=1:horzLen+1,y=lamCurrStar,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
		Guide.xlabel("Time"), Guide.ylabel(raw"Lambda ($/kA)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_decentral_ALADIN_Lam.png", 24inch, 12inch), pd4alad) end


lamPlotalad=plot(Lam[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=lamCurrStar,Geom.line,Theme(default_color=colorant"black",line_width=4pt)),
			Guide.xlabel("Time"), Guide.ylabel("Lambda"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
uSumPlotalad=plot(uSum[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=uSumStar,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("U sum"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
zSumPlotalad=plot(zSum[:,2:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=zSumStar,Geom.line,Theme(default_color=colorant"black",line_width=3pt)),
			Guide.xlabel("Time"), Guide.ylabel("Z sum"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
constPlotalad2=plot(currConst[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel("curr constraint diff"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_ALADIN_LamConv.png", 36inch, 12inch), lamPlotalad) end

activeSet=zeros(convIt,1)
setChanges=zeros(convIt,1)
for ii=2:convIt
    activeSet[ii,1]=sum(abs.(Cs[:,ii]))+sum(abs.(Ct[:,ii]))+
              sum(abs.(Cu[:,ii]))+sum(abs.(Cz[:,ii]))
    setChanges[ii,1]=sum(abs.(Cs[:,ii]-Cs[:,ii-1]))+sum(abs.(Ct[:,ii]-Ct[:,ii-1]))+
                     sum(abs.(Cu[:,ii]-Cu[:,ii-1]))+sum(abs.(Cz[:,ii]-Cz[:,ii-1]))
end
activeSetPlot=plot(x=2:convIt,y=activeSet[2:convIt],Geom.line,
                   Guide.xlabel("Iteration"), Guide.ylabel("Total Active inequality constraints",orientation=:vertical),
                   Coord.Cartesian(xmin=2,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
       			   minor_label_font_size=26pt,key_label_font_size=26pt))
setChangesPlot=plot(x=3:convIt,y=setChanges[3:convIt],Geom.line,
                    Guide.xlabel("Iteration"), Guide.ylabel("Changes in Active inequality constraints",orientation=:vertical),
                    Coord.Cartesian(xmin=3,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
        			minor_label_font_size=26pt,key_label_font_size=26pt))
solChangesplot=plot(layer(x=2:convIt,y=deltaY[2:convIt],Geom.line),
                    layer(x=2:convIt,y=convCheck[2:convIt],Geom.line),Scale.y_log)

convItPlotalad=plot(x=1:convIt,y=itConvALAD[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
convPlotalad=plot(x=1:convIt,y=ConvALAD[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("central lambda gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
constPlotalad=plot(x=1:convIt,y=constConvALAD[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("const gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
fPlotalad=plot(x=1:convIt-1,y=fConvALAD[1:convIt-1,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("obj function gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_ALADIN_Conv.png", 36inch, 12inch), vstack(convItPlotalad,convPlotalad,fPlotalad)) end
