#EVC with ALADIN for PWL convex relaxation
#Micah Botkin-Levy
#5/22/18

#formulation try 2
#u, xn, xt, and z are all in "y", v is the auxilliary variable corresponding to x in literature
#current constraint is coupling


datafile="jld" #"mat" #"jld" #"n"
drawFig=0
noTlimit=0
if datafile in ["mat" "jld"]; N=30 end

println("Loading Packages...")

using Gadfly
using JuMP
using Gurobi
using Cairo #for png output
using Fontconfig
using Distributions
#using ProximalOperators

if datafile=="mat"
	using MAT #to read in scenarios from matlab
elseif datafile=="jld"
	using JLD
end

if datafile in ["mat" "jld"]
	println("Reading in Data...")

	function string_as_varname(s::String,v::Any)
		 s=Symbol(s)
		 @eval (($s) = ($v))
	end

	#read in mat scenario
	path="C:\\Users\\micah\\Documents\\uvm\\Research\\EVC code\\N$(N)\\"
	file="EVCscenarioN$(N)."*datafile
	if datafile=="mat"
		vars = matread(path*file)
	elseif datafile=="jld"
		vars=load(path*file)
	end
	varnames=keys(vars)
	varNum=length(varnames)
	varKeys=collect(varnames)
	varValues=collect(values(vars))

	for i =1:varNum
		n=varKeys[i]
		v=varValues[i]
		if datafile=="mat"
			if n in ["N" "K" "S"]
				v=convert(Int, v)
			end
		end
		string_as_varname(n,v)
	end
	println("done reading in")

	if datafile=="mat"
		Kn=convert(Array{Int,2},Kn)
	end
end

tic()
#initialize with current states
xn0=s0
xt0=T0

stepI = 1;
horzLen=K1
epsilon = 1e-4
numIteration=20
convIt=numIteration
Conv=zeros(numIteration,1)
itConv=zeros(numIteration,1)
constConv=zeros(numIteration,1)
fConv=zeros(numIteration,1)
xnConv=zeros(numIteration,1)
avgN=zeros(numIteration,1)

#ALADIN tuning and initial guess
#H=Qsi
H=vcat(.001*Qsi[1:N,1],1)
#H=vcat(10*ones(N,1),1)


#vn0=rand(Truncated(Normal(0), 0, 1), N*(horzLen+1))
#vn0=uStar
vn0=max.(uStar + rand(Truncated(Normal(0), -0.1, 0.1), N*(horzLen+1)),0)
#vz0=ItotalMax*rand(Truncated(Normal(0), 0, 1), S*(horzLen+1))
#vz0=zStar
vz0=max.(zStar-rand(Truncated(Normal(0), 0, 5), S*(horzLen+1)),0)
lambda0=rand(Truncated(Normal(0), 0, 5), horzLen+1)
#lambda0=lamCurrStar
#lambda0=zeros(horzLen+1,1)

#save matrices
Gn=SharedArray{Float64}(N*(horzLen+1),numIteration) #row are time (N states for k=1, them N states for k=2),  columns are iteration
Gz=zeros(S*(horzLen+1),numIteration) #row are time (N states for k=1, them N states for k=2),  columns are iteration
Un=SharedArray{Float64}(N*(horzLen+1),numIteration) #row are time (N states for k=1, them N states for k=2),  columns are iteration
Xn=SharedArray{Float64}(N*(horzLen+1),numIteration)  #row are time,  columns are iteration
Xt=zeros((horzLen+1),numIteration)  #row are time,  columns are iteration
Z=zeros(S*(horzLen+1),numIteration)  #row are time,  columns are iteration
Vn=zeros((N)*(horzLen+1),numIteration) #row are time,  columns are iteration
Vn[:,1]=vn0 #initial guess goes in column 1
Vz=zeros(S*(horzLen+1),numIteration) #row are time,  columns are iteration
Vz[:,1]=vz0
Lam=zeros((horzLen+1),numIteration) #(rows are time, columns are iteration)
Lam[:,1]=lambda0

ZS=zeros((horzLen+1),numIteration)  #row are time,  columns are iteration

for p=1:numIteration-1

    #solve decoupled
    @sync @parallel for evInd=1:N

        lambda=Lam[:,p]
        evV=Vn[collect(evInd:N:length(Vn[:,p])),p]
        #evV=zeros(horzLen+1,1)
        target=zeros((horzLen+1),1)
        target[(Kn[evInd,1]-(stepI-1)):1:length(target),1]=Sn[evInd,1]
        evM = Model(solver = GurobiSolver())
        @variable(evM,xn[1:(horzLen+1)])
        @variable(evM,u[1:(horzLen+1)])
        objFun(xn,u)=sum((xn[k,1]-1)^2*Qsi[evInd,1] for k=1:horzLen+1) +
                        sum((u[k,1])^2*Ri[evInd,1]     for k=1:horzLen+1)
        #@objective(evM,Min, sum(objFun(xn,u)+sum(lambda[k,1]*(u[k,1]-evV[k,1]) for k=1:horzLen+1)+rho_p/2*sum((u[k,1]-evV[k,1])^2 for k=1:horzLen+1)))
        @objective(evM,Min, sum(objFun(xn,u)+sum(lambda[k,1]*(u[k,1]) for k=1:horzLen+1)+1/2*sum((u[k,1]-evV[k,1])*H[evInd,1]*(u[k,1]-evV[k,1]) for k=1:horzLen+1)))
        @constraint(evM,xn[1,1]==xn0[evInd,1]+eta[evInd,1]*u[1,1])
        @constraint(evM,[k=1:horzLen],xn[k+1,1]==xn[k,1]+eta[evInd,1]*u[k+1,1])
        @constraint(evM,xn.<=1)
        @constraint(evM,xn.>=target)
        @constraint(evM,u.<=imax[evInd,1])
        @constraint(evM,u.>=imin[evInd,1])


        TT = STDOUT # save original STDOUT stream
        redirect_stdout()
        status = solve(evM)
        redirect_stdout(TT)
        if status!=:Optimal
            return
        else
            Xn[collect(evInd:N:length(Xn[:,p+1])),p+1]=getvalue(xn)
    		Un[collect(evInd:N:length(Un[:,p+1])),p+1]=getvalue(u)
            Gn[collect(evInd:N:length(Gn[:,p+1])),p+1]=H[evInd,1]*(evV-getvalue(u))-lambda
            #update H here???
        end
    end


    #N+1 decoupled problem aka transformer current
    tM = Model(solver = GurobiSolver())
    @variable(tM,z[1:(S)*(horzLen+1)])
    @variable(tM,xt[1:(horzLen+1)])
    constFun1(u,v)=sum(Lam[k,p]*sum(-u[(k-1)*(S)+s,1] for s=1:S)  for k=1:(horzLen+1))
    constFun2(u,v)=1/2*sum(sum(u[(k-1)*(S)+s,1]-v[(k-1)*(S)+s,1] for s=1:S)*H[N+1,1]*sum(u[(k-1)*(S)+s,1]-v[(k-1)*(S)+s,1] for s=1:S) for k=1:(horzLen+1))
    @objective(tM,Min, constFun1(z,Vz[:,p])+constFun2(z,Vz[:,p]))
    @constraint(tM,tempCon1,xt[1,1]==tau*xt0+gamma*deltaI*sum((2*m+1)*z[m+1,1] for m=0:S-1)+rho*w[stepI*2,1])
    @constraint(tM,tempCon2[k=1:horzLen],xt[k+1,1]==tau*xt[k,1]+gamma*deltaI*sum((2*m+1)*z[k*S+(m+1),1] for m=0:S-1)+rho*w[stepI*2+k*2,1])
    if noTlimit==0
    	@constraint(tM,upperTCon,xt.<=Tmax)
    end
    @constraint(tM,xt.>=0)
    @constraint(tM,z.>=0)
    @constraint(tM,z.<=deltaI)
    TT = STDOUT # save original STDOUT stream
    redirect_stdout()
    status = solve(tM)
    redirect_stdout(TT)
    if status!=:Optimal
        return
    else
        Xt[:,p+1]=getvalue(xt)
        Z[:,p+1]=getvalue(z)
        Gz[:,p+1]=H[N+1,1]*(Vz[:,p]-getvalue(z))-repeat(-Lam[:,p],inner=S)
    end

	for k=1:horzLen+1
		ZS[k,p+1]=sum(Z[(k-1)*(S)+s,p+1] for s=1:S)
	end

    #check for convergence
    convCheck=norm.(vcat(Vn[:,p]-Un[:,p+1],Vz[:,p]-Z[:,p+1]))
    avgN[p,1]=mean(convCheck)
    objFun(xn,xt,u)=sum(sum((xn[(k-1)*(N)+n,1]-1)^2*Qsi[n,1]     for n=1:N) for k=1:horzLen+1) +
                    sum((xt[k,1]-1)^2*Qsi[N+1,1]                 for k=1:horzLen+1) +
                    sum(sum((u[(k-1)*N+n,1])^2*Ri[n,1]           for n=1:N) for k=1:horzLen+1)
    fGap= abs(objFun(Xn[:,p+1],Xt[:,p+1],Un[:,p+1])-fStar)
    xnGap=norm((Xn[:,p+1]-xnStar),2)
    itGap = norm(Lam[:,p]-Lam[:,max(p-1,1)],2)
    convGap = norm(Lam[:,p]-lamCurrStar,2)
    fConv[p,1]=fGap
    xnConv[p,1]=xnGap
    itConv[p,1]=itGap
    Conv[p,1]=convGap
    if all(convCheck.<=epsilon)
        @printf "Converged after %g iterations\n" p
        convIt=p+1
        break
    else
        @printf "lastGap %e after %g iterations\n" itGap p
        @printf "convGap %e after %g iterations\n" convGap p
        @printf "xnGap %e after %g iterations\n" xnGap p
        @printf("fGap %e after %g iterations\n\n",fGap,p)

    end

    #coupled QP
    cM = Model(solver = GurobiSolver())
    @variable(cM,dUn[1:(N)*(horzLen+1)])
    @variable(cM,dZ[1:(S)*(horzLen+1)])
    coupledObj(deltaY,Hi,gi)=1/2*deltaY'*Hi*deltaY+gi'*deltaY
    @objective(cM,Min, sum(sum(coupledObj(dUn[collect(n:N:length(dUn[:,1])),1],H[n,1],Gn[collect(n:N:length(Gn[:,p+1])),p+1]) for n=1:N)+
                            coupledObj(dZ,H[N+1,1],Gz[:,p+1])))
    @constraint(cM,currCon[k=1:horzLen+1],0==-sum(Un[(k-1)*(N)+n,p+1]+dUn[(k-1)*(N)+n,1] for n=1:N)-w[(k-1)*2+1]+
                                             sum(Z[(k-1)*(S)+s,p+1]+dZ[(k-1)*(S)+s,1] for s=1:S))
 	@constraint(cM,(Z[:,p+1]+dZ).>=0)
    @constraint(cM,(Z[:,p+1]+dZ).<=deltaI)
	@constraint(cM,(Un[:,p+1]+dUn).<=repmat(imax,horzLen+1,1))
	@constraint(cM,(Un[:,p+1]+dUn).>=repmat(imin,horzLen+1,1))
    TT = STDOUT # save original STDOUT stream
    redirect_stdout()
    status = solve(cM)
    redirect_stdout(TT)
    if status!=:Optimal
        return
    end

    #update step
    Lam[:,p+1]=-getdual(currCon)
    Vn[:,p+1]=Un[:,p+1]+getvalue(dUn)
    Vz[:,p+1]=Z[:,p+1]+getvalue(dZ)

end


println("plotting....")
xPlot=zeros(horzLen+1,N)
for ii= 1:N
	xPlot[:,ii]=Xn[collect(ii:N:length(Xn[:,convIt])),convIt]
end

#plot(x=1:horzLen+1,y=xPlot2[:,ii])
# plot(x=1:Kn[ii,1],y=xPlot2[1:Kn[ii,1],ii])

pd1alad=plot(xPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV SOC"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_central_ALADIN_SOC.png", 24inch, 12inch), pd1alad) end


uPlot=zeros(horzLen+1,N)
for ii= 1:N
	uPlot[:,ii]=Un[collect(ii:N:length(Un[:,convIt])),convIt]
end

pd2alad=plot(uPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV Current (A)"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_central_ALADIN_Curr.png", 24inch, 12inch), pd2alad) end

pd3alad=plot(layer(x=1:horzLen+1,y=Xt[:,convIt],Geom.line,Theme(default_color=colorant"blue")),
		#layer(x=1:horzLen+1,y=Tactual,Geom.line,Theme(default_color=colorant"green")),
		yintercept=[Tmax],Geom.hline(color=["red"],style=:dot),
		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",key_position = :top,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
		#Guide.manual_color_key("", ["PWL Temp", "Actual Temp"], ["blue", "green"]))
if drawFig==1 draw(PNG(path*"J_central_ALADIN_Temp.png", 24inch, 12inch), pd3alad) end

pd4alad=plot(layer(x=1:horzLen+1,y=Lam[:,convIt],Geom.line),
		layer(x=1:horzLen+1,y=lamCurrStar,Geom.line,Theme(default_color=colorant"black")),
		Guide.xlabel("Time"), Guide.ylabel(raw"Lambda ($/A)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_central_ALADIN_Lam.png", 24inch, 12inch), pd4alad) end


lamPlotalad=plot(Lam[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			layer(x=1:horzLen+1,y=lamCurrStar,Geom.line,Theme(default_color=colorant"black")),
			Guide.xlabel("Time"), Guide.ylabel("Lambda"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
zSumPlotalad=plot(ZS[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
				layer(x=1:horzLen+1,y=zSumStar,Geom.line,Theme(default_color=colorant"black")),
			Guide.xlabel("Time"), Guide.ylabel("Z sum"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_ALADIN_LamConv.png", 36inch, 12inch), lamPlotalad) end

convItPlotalad=plot(x=1:convIt,y=itConv[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
convPlotalad=plot(x=1:convIt,y=Conv[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("central lambda gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
fPlotalad=plot(x=1:convIt-1,y=fConv[1:convIt-1,1],Geom.line,#Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("obj function gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_ALADIN_Conv.png", 36inch, 12inch), convPlotalad) end
