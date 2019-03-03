
function test()
    for i = 1:10
       global z
       z = i
    end
end
test()



function innerloop()
    cc=cc+1
end

function wrapper()
    # cc=1
    for i=1:10
        global cc
        cc=2
        #innerloop()
    end
end



wrapper()

using JuMP
using Gurobi
using Dates
using Suppressor
using Printf
silent=false
solverSilent=false

Juno.progress() do id
      for i = 1:10
          @printf "%s: time step %g of %g....\n" Dates.format(Dates.now(),"HH:MM:SS") i 10
          @info "$(Dates.format(Dates.now(),"HH:MM:SS")): $(i) of $(10)....\n" progress=i/10 _id=id
          sleep(0.5)

          m=Model(solver=GurobiSolver())
          @variable(m,x)
          @objective(m,Min,sum(x))
          @constraint(m,x>=0)

          if solverSilent
              @suppress_out begin
                  status=solve(m)
              end
          else
              status=solve(m)
          end

      end
  end


using JuMP
using Gurobi

m=Model(solver=GurobiSolver())
@variable(m,q)
@variable(m,e[1:2])
@variable(m,p)
@variable(m,h,Bin)
@variable(m,l,Bin)
@objective(m,Min,-q)
@constraint(m,q<=(h*0.9+l*.7)*p)
@constraint(m,h+l==1)
@constraint(m,q>=0)
@constraint(m,e[2]==e[1]+q*.5)
@constraint(m,e[1]==5)
@constraint(m,p>=0)
@constraint(m,p<=100)

status=solve(m)




#checking Pauls XFRM equation continous
using QuadGK
using Plots; pyplot()
b=0.0149
t0=0
L=25e3/3
Ta=40
T0=40
alpha = 0.178
a=0.000939e-6;
# t=1:100
t_range=1:13*60
T=zeros(length(t_range))
for i=1:length(t_range)
    t=t_range[i]
    f(tau) = exp(-b*(t-tau))*(a*L^2+b*Ta+alpha)
    # f(tau) = exp(-b*(t-tau))*(0.000939*L^2+b*Ta+0.178)
    (T_1,E)=quadgk(f,0,t, rtol=1e-9)
    T[i]=exp(-b*t)*T0+T_1
end
plot(t_range,T)




#testing new model in Celsius disctete
using Plots;pyplot()
xfrmR  = 10e3/3                          # single phase transformer rating VA
Vac = 240                              # PEV battery rms voltage --- V [used in PEV kW -> kA conversion]
Vtf = 8320                             # distr-level transformer rms voltage --- V [used in inelastic kW -> kA conv]
Ntf   = Vtf/Vac                        # pole-top transformer turns ratio
Tamb=60
T0=Tamb
power_weight = 0.000939/(1e3/3/2.5)^2/60
beta = 0.0149/60
alpha = 0.178/60
# Ts=3 #minutes
# K=round(Int,13*60/Ts)
Ts=3*60 #minutes
K=round(Int,13*60*60/Ts)
τP = exp(- Ts*beta)
ρP = (1 - τP)
curr_weight = power_weight*Vac^2
γP = 1/beta*ρP*curr_weight
ItotalMax =  (xfrmR/Vtf)*Ntf#kA
T=zeros(K,1)
u=ItotalMax*ones(K,1)
#p=xfrmR*ones(K,1)
T[1]=T0
for k=1:K-1
    T[k+1,1]=τP*T[k,1]+γP*u[k]^2+ρP*(Tamb+alpha/beta)
    #T[k+1,1]=τP*T[k,1]+γP*p[k]^2+ρP*(Tamb+1/beta*alpha)  #works
end
plot(collect(Ts*(1:K-1)),T)
#plot(T)








#while debugging scenario
using Gurobi
#initialize
t0=evS.t0
s0=evS.s0

horzLen=evS.K1
K=evS.K
N=evS.N
S=evS.S
#S=8
cSol=solutionStruct(K=K,N=N,S=S)
cSave=centralLogStruct(logLength=length(saveLogInd),horzLen=horzLen,N=N,S=S)

stepI=1
deltaI=evS.ItotalMax/S
horzLen=min(evS.K1,K-stepI)

#desired SOC
target=zeros(N*(horzLen+1),1);
for ii=1:N
   cur=evS.Kn[ii]-(stepI-1)
   ind=(max(0,(cur-1)*N)+ii):N:length(target)
   target[ind].=evS.Snmin[ii,1]
end

centralModel = Model(solver = GurobiSolver(NumericFocus=3))
@variable(centralModel,sn[1:(N)*(horzLen+1)])
@variable(centralModel,u[1:(N)*(horzLen+1)])
@variable(centralModel,t[1:(horzLen+1)])
@variable(centralModel,z[1:S*(horzLen+1)])
objExp =sum((sn[n,1]-1)^2*evS.Qsi[n,1]+(u[n,1])^2*evS.Ri[n,1] for n=1:N)
for k=2:horzLen+1
	append!(objExp,sum((sn[(k-1)*(N)+n,1]-1)^2*evS.Qsi[n,1]+(u[(k-1)*N+n,1])^2*evS.Ri[n,1]  for n=1:N))
end
@objective(centralModel,Min,objExp)
@constraint(centralModel,stateCon1,sn[1:N,1].==s0[1:N,1]+evS.ηP[:,1].*u[1:N,1])
@constraint(centralModel,stateCon2[k=1:horzLen,n=1:N],sn[n+(k)*(N),1]==sn[n+(k-1)*(N),1]+evS.ηP[n,1]*u[n+(k)*(N),1])
@constraint(centralModel,tempCon1,t[1,1]==evS.τP*t0+evS.γP*deltaI*sum((2*s-1)*z[s,1] for s=1:S)+evS.ρP*evS.Tamb[stepI,1])
@constraint(centralModel,tempCon2[k=1:horzLen],t[k+1,1]==evS.τP*t[k,1]+evS.γP*deltaI*sum((2*s-1)*z[k*S+s,1] for s=1:S)+evS.ρP*evS.Tamb[stepI+k,1])
@constraint(centralModel,currCon[k=1:horzLen+1],0==-sum(u[(k-1)*(N)+n] for n=1:N)-evS.iD_pred[stepI+(k-1)]+sum(z[(k-1)*(S)+s] for s=1:S))
@constraint(centralModel,sn.<=1)
@constraint(centralModel,sn.>=target)
@constraint(centralModel,upperTCon,t.<=evS.Tmax)
@constraint(centralModel,t.>=0)
@constraint(centralModel,upperCCon,u.<=repeat(evS.imax,horzLen+1,1))
@constraint(centralModel,u.>=repeat(evS.imin,horzLen+1,1))
@constraint(centralModel,z.>=0)
@constraint(centralModel,z.<=deltaI)

status = solve(centralModel)
@assert status==:Optimal "Central optimization not solved to optimality"

uRaw=getvalue(u)
snRaw=getvalue(sn)
tRaw=getvalue(t)
zRaw=getvalue(z)

#calculate actual temp
Tactual=zeros(horzLen+1,1)
itotal=zeros(horzLen+1,1)
for k=1:horzLen+1
	itotal[k,1]=sum((uRaw[(k-1)*N+n,1]) for n=1:N) + evS.iD_actual[stepI+(k-1),1]
	#itotal[k,1]= evS.iD_actual[stepI+(k-1),1]
end
Tactual[1,1]=evS.τP*t0+evS.γP*itotal[1,1]^2+evS.ρP*evS.Tamb[stepI,1]
for k=1:horzLen
	Tactual[k+1,1]=evS.τP*Tactual[k,1]+evS.γP*itotal[k+1,1]^2+evS.ρP*evS.Tamb[stepI+k,1]
end
lambdaCurr=-getdual(currCon)
if noTlimit==false
	lambdaUpperT=-getdual(upperTCon)
else
	lambdaUpperT=zeros(horzLen+1,1)
end

uSum=zeros(horzLen+1,1)
zSum=zeros(horzLen+1,1)
for k=1:horzLen+1
	zSum[k,1]=sum(zRaw[(k-1)*(S)+s,1] for s=1:S)
	uSum[k,1]=sum(uRaw[(k-1)*N+n,1] for n=1:N)
end

xPlot=zeros(horzLen+1,N)
uPlot=zeros(horzLen+1,N)
for ii= 1:N
	xPlot[:,ii]=snRaw[collect(ii:N:length(snRaw))]
	uPlot[:,ii]=uRaw[collect(ii:N:length(uRaw))]
end

p1=plot(xPlot,xlabel="Time",ylabel="PEV SOC",xlims=(1,horzLen+1))

p2=plot(uPlot,xlabel="Time",ylabel="PEV Current (kA)",xlims=(1,horzLen+1))

p3=plot(hcat(tRaw,Tactual,evS.Tamb_raw[2:horzLen+2]),label=["Pred Temp" "Actual Temp" "Ambient Temp"],xlims=(1,horzLen+1),xlabel="Time",ylabel="Temp (K)")
plot!(p3,1:horzLen+1,evS.Tmax*ones(horzLen+1),label="XFRM Limit",line=(:dash,:red))

aggU=plot(hcat(uSum,itotal,evS.iD_actual[1:horzLen+1]),label=["Central" "Total" "iD"],xlims=(0,horzLen+1),xlabel="Time",ylabel="Current (kA)")
plot!(aggU,1:horzLen+1,evS.ItotalMax*ones(horzLen+1),label="XFRM Limit",line=(:dash,:red))
plot(uSum)
p4b=plot(lambdaUpperT,xlabel="Time",ylabel=raw"Lambda ($/K)",xlims=(1,evS.K),legend=false)
p4=plot(lambdaCurr,xlabel="Time",ylabel=raw"Lambda ($/kA)",xlims=(1,evS.K),legend=false)
