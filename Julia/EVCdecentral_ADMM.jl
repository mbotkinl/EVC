#EVC with ADMM for PWL convex relaxation
#Micah Botkin-Levy
#5/22/18

#formulation try 2
#u, xn, xt, and z are all in "x" v is the auxilliary variable corresponding to z in literature
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
convChk = 1e-6
numIteration=1000
convIt=numIteration
Conv=zeros(numIteration,1)
itConv=zeros(numIteration,1)
constConv=zeros(numIteration,1)
fConv=zeros(numIteration,1)
xnConv=zeros(numIteration,1)

#admm  initial parameters and guesses
rhoADMM=10^0
d = Truncated(Normal(0), 0, 5)
lambda0=rand(d, horzLen+1)
#u0=repmat(imax,horzLen+1,1) #guess max charging
#u0=rand(d, N*(horzLen+1)) #make better guess here (and has to be within limit)

#u w and z are one index ahead of x. i.e the x[k+1]=x[k]+eta*u[k+1]
Un=zeros(N*(horzLen+1),numIteration) #row are time,  columns are iteration
#Un[:,1]=u0
Lam=zeros((horzLen+1),numIteration) #(rows are time, columns are iteration)
Lam[:,1]=lambda0
Xn=zeros(N*(horzLen+1),numIteration)  #row are time,  columns are iteration
Xt=zeros((horzLen+1),numIteration)  #row are time,  columns are iteration
Z=zeros(S*(horzLen+1),numIteration)  #row are time,  columns are iteration
Vn=zeros((N)*(horzLen+1),numIteration) #row are time,  columns are iteration
vn0=uStar
Vn[:,1]=vn0
Vz=zeros((horzLen+1),numIteration)
Vz0=zeros((horzLen+1),1)
for k=1:horzLen+1
    Vz0[k,1]=sum(zStar[(k-1)*(S)+(1:S)])
end
Vz[:,1]=Vz0


for p=1:numIteration-1

    #x minimization eq 7.66 in Bertsekas
    @parallel for evInd=1:N
        evV=Vn[collect(evInd:N:length(Vn[:,p])),p]
        target=zeros((horzLen+1),1)
        target[Kn[evInd,1]:1:length(target),1]=s0[evInd,1]
    	evM = Model(solver = GurobiSolver())
    	@variable(evM,xn[1:(horzLen+1)])
    	@variable(evM,u[1:(horzLen+1)])
    	objFun(xn,xt,u)=sum((xn[k,1]-1)^2*Qsi[evInd,1] for k=1:horzLen+1) +
        			    sum((u[k,1])^2*Ri[evInd,1]    for k=1:horzLen+1)
        constFun1(u,v)=sum(Lam[k,p]*(u[k,1]-v[k,1])  for k=1:horzLen+1)
        constFun2(u,v)=rhoADMM/2*sum((u[k,1]-v[k,1])^2  for k=1:horzLen+1)
    	@objective(evM,Min, objFun(xn,xt,u)+constFun1(u,evV)+constFun2(u,evV))
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
    		Xn[collect(evInd:N:length(Xn[:,p])),p+1]=getvalue(xn)
    		Un[collect(evInd:N:length(Un[:,p])),p+1]=getvalue(u)
    	end
    end

    #N+1 decoupled problem aka transformer current
    tM = Model(solver = GurobiSolver())
    @variable(tM,z[1:(S)*(horzLen+1)])
    @variable(tM,xt[1:(horzLen+1)])
    @variable(tM,zSum[1:(horzLen+1)])
    constFun1(u,v)=sum(Lam[k,p]*(u[k,1]-v[k,1])  for k=1:horzLen+1)
    constFun2(u,v)=rhoADMM/2*sum((u[k,1]-v[k,1])^2  for k=1:horzLen+1)
    @objective(tM,Min, constFun1(zSum,Vz)+constFun2(zSum,Vz))
    @constraint(tM,tempCon1,xt[1,1]==tau*xt0+gamma*deltaI*sum((2*m+1)*z[m+1,1] for m=0:S-1)+rho*w[stepI*2,1])
    @constraint(tM,tempCon2[k=1:horzLen],xt[k+1,1]==tau*xt[k,1]+gamma*deltaI*sum((2*m+1)*z[k*S+(m+1),1] for m=0:S-1)+rho*w[stepI*2+k*2,1])
    if noTlimit==0
    	@constraint(tM,upperTCon,xt.<=Tmax)
    end
    @constraint(tM,xt.>=0)
    @constraint(tM,z.>=0)
    @constraint(tM,z.<=deltaI)
    @constraint(tM,zC[k=1:horzLen+1],zSum[k,1]==sum(z[(k-1)*(S)+(1:S)]))

    TT = STDOUT # save original STDOUT stream
    redirect_stdout()
    status = solve(tM)
    redirect_stdout(TT)
    if status!=:Optimal
        return
    else
        Xt[:,p+1]=getvalue(xt)
        Z[:,p+1]=getvalue(z)
    end

    #lambda update eq 7.68
    currConst=zeros(horzLen+1,1)
	for k=1:horzLen+1
		currConst[k,1]=sum(Un[(k-1)*N+n,p+1] for n=1:N) + w[(k-1)*2+(stepI*2-1),1] - sum(Z[(k-1)*(S)+(1:S),p+1])
		Lam[k,p+1]=Lam[k,p]+rhoADMM/(horzLen+1)*(currConst[k,1])
	end


    #v upate eq 7.67
    for k=1:horzLen+1
        Vn[(k-1)*N+1:N,p+1]=Un[(k-1)*N+1:N,p+1]+(Lam[k,p]-Lam[k,p+1])/rhoADMM
        Vz[k,p+1]=sum(Z[(k-1)*(S)+(1:S),p+1])+(Lam[k,p]-Lam[k,p+1])/rhoADMM
    end


    #check convergence
	fGap= objFun(Xn[:,p+1],Xt[:,p+1],Un[:,p+1])-fStar
	xnGap=norm((Xn[:,p+1]-xnStar),2)
	constGap=norm(currConst,2)
	itGap = norm(Lam[:,p+1]-Lam[:,p],2)
	convGap = norm(Lam[:,p+1]-lamCurrStar,2)
	fConv[p,1]=fGap
	xnConv[p,1]=xnGap
	constConv[p,1]=constGap
	itConv[p,1]=itGap
	Conv[p,1]=convGap
	if(itGap <= convChk )
		@printf "Converged after %g iterations\n" p
		convIt=p
		break
	else
		@printf "lastGap %e after %g iterations\n" itGap p
		@printf "convGap %e after %g iterations\n" convGap p
        @printf "xnGap %e after %g iterations\n\n" xnGap p

	end
end
toc()




println("plotting....")
xPlot=zeros(horzLen+1,N)
for ii= 1:N
	xPlot[:,ii]=Xn[collect(ii:N:length(Xn[:,convIt])),convIt]
end

#plot(x=1:horzLen+1,y=xPlot2[:,ii])
# plot(x=1:Kn[ii,1],y=xPlot2[1:Kn[ii,1],ii])

pd1admm=plot(xPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV SOC"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_central_ADMM_SOC.png", 24inch, 12inch), pd1admm) end


uPlot=zeros(horzLen+1,N)
for ii= 1:N
	uPlot[:,ii]=Un[collect(ii:N:length(Un[:,convIt])),convIt]
end

pd2admm=plot(uPlot,x=Row.index,y=Col.value,color=Col.index,Geom.line,
		Guide.xlabel("Time"), Guide.ylabel("PEV Current (A)"),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),
		Theme(background_color=colorant"white",key_position = :none,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_central_ADMM_Curr.png", 24inch, 12inch), pd2admm) end

pd3admm=plot(layer(x=1:horzLen+1,y=Xt[:,convIt],Geom.line,Theme(default_color=colorant"blue")),
		#layer(x=1:horzLen+1,y=Tactual,Geom.line,Theme(default_color=colorant"green")),
		yintercept=[Tmax],Geom.hline(color=["red"],style=:dot),
		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",key_position = :top,major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
		#Guide.manual_color_key("", ["PWL Temp", "Actual Temp"], ["blue", "green"]))
if drawFig==1 draw(PNG(path*"J_central_ADMM_Temp.png", 24inch, 12inch), pd3admm) end

pd4admm=plot(x=1:horzLen+1,y=Lam[:,convIt],Geom.line,
		Guide.xlabel("Time"), Guide.ylabel(raw"Lambda ($/A)",orientation=:vertical),
		Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=18pt,
		minor_label_font_size=16pt,key_label_font_size=16pt))
if drawFig==1 draw(PNG(path*"J_central_ADMM_Lam.png", 24inch, 12inch), pd4admm) end

fName="J_Central.png"


lamPlotadmm=plot(Lam[:,1:convIt],x=Row.index,y=Col.value,color=Col.index,Geom.line,
			Guide.xlabel("Time"), Guide.ylabel("Lambda"),Guide.ColorKey(title="Iteration"),
			Coord.Cartesian(xmin=0,xmax=horzLen+1),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_ADMM_LamConv.png", 36inch, 12inch), lamPlotadmm) end

convItPlotadmm=plot(x=1:convIt,y=itConv[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
convPlotadmm=plot(x=1:convIt,y=Conv[1:convIt,1],Geom.line,Scale.y_log10,
			Guide.xlabel("Iteration"), Guide.ylabel("2-Norm Lambda Gap"),
			Coord.Cartesian(xmin=0,xmax=convIt),Theme(background_color=colorant"white",major_label_font_size=30pt,line_width=2pt,
			minor_label_font_size=26pt,key_label_font_size=26pt))
if drawFig==1 draw(PNG(path*"J_ADMM_Conv.png", 36inch, 12inch), convPlotadmm) end
