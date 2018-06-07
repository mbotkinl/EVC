using JuMP
using Gurobi

N=1
horzLen=94
S=50
xt0=[370]
sn0=[0.753]
stepI=1
Kn=[84]
Snmin=[.834405]
Qsi=[90;0]
Ri=[0.00013]
tauP=0.888
gammaP=4.88e-6
rhoP=0.111
etaP=0.000457
deltaI=80

println("setting up model")
m = Model(solver = GurobiSolver(Presolve=0,BarHomogeneous=1))

#u w and z are one index ahead of x. i.e the x[k+1]=x[k]+eta*u[k+1]
@variable(m,u[1:N*(horzLen+1)])
@variable(m,sn[1:(N)*(horzLen+1)])
@variable(m,xt[1:(horzLen+1)])
@variable(m,z[1:S*(horzLen+1)])

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
@constraint(m,stateCon1,sn[1:N,1].==sn0[1:N,1]+etaP.*u[1:N,1])
@constraint(m,stateCon2[k=1:horzLen,n=1:N],sn[n+(k)*(N),1]==sn[n+(k-1)*(N),1]+etaP[n,1]*u[n+(k)*(N),1])
@constraint(m,tempCon1,xt[1,1].==tauP*xt0+gammaP*deltaI*sum((2*m+1)*z[m+1,1] for m=0:S-1)+rhoP*370)
@constraint(m,tempCon2[k=1:horzLen],xt[k+1,1]==tauP*xt[k,1]+gammaP*deltaI*sum((2*m+1)*z[k*S+(m+1),1] for m=0:S-1)+rhoP*370)
@constraint(m,currCon[k=1:horzLen+1],0==-sum(u[(k-1)*(N)+n] for n=1:N)+sum(z[(k-1)*(S)+s] for s=1:S))
@constraint(m,sn.<=1)
@constraint(m,sn.>=target)
@constraint(m,upperTCon,xt.<=370.01 )
@constraint(m,xt.>=0)
@constraint(m,upperCCon,u.<=22)
@constraint(m,u.>=0)
@constraint(m,z.>=0)
@constraint(m,z.<=deltaI)

println("solving....")
statusM = solve(m)

compV=getvalue(z)[1,1]


# for i=1:100000000
#     if !(compV==compV)
#         println("not equal!!")
#         break
#     end
# end
#
# tS=@sprintf("val: %f",compV)
# for i=1:100000000
#     if !(@sprintf("val: %f",compV)==tS)
#         println("not equal!!")
#         break
#     end
# end

#@printf "val: %f" compV

# for ii=1:2000
#     println(compV)
#     @printf "val: %f \n" compV
#     #compV
# end
