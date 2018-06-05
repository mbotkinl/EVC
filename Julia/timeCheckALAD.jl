

using JuMP
using Gurobi
using Distributions


# function timeCheckModel2(opt)
#
#     S=50
#     horzLen=283
#     xt0=368
#     tauP=0.888
#     gammaP=4.88e-6
#     rhoP=0.111
#     Vz=4000*rand(Truncated(Normal(0), 0, 1), S*(horzLen+1))
#     tM = Model(solver = GurobiSolver())
#     @variable(tM,z[1:(S)*(horzLen+1)])
#     @variable(tM,xt[1:(horzLen+1)])
#     if opt==1
#         objExp=0*z[1,1]
#         for k=1:(horzLen+1)
#             objExp=objExp+sum(-z[(k-1)*(S)+s,1] for s=1:S)
#             zDiff=sum(z[(k-1)*(S)+s,1]-Vz[(k-1)*(S)+s,1] for s=1:S)
#             objExp=objExp+1/2*zDiff^2
#         end
#         @objective(tM,Min,objExp)
#     else
#         constFun1(u,v)=sum(sum(u[(k-1)*(S)+s,1] for s=1:S)  for k=1:(horzLen+1))
#         constFun2(u,v)=sum(sum((u[(k-1)*(S)+s,1]-v[(k-1)*(S)+s,1])*(u[(k-1)*(S)+s,1]-v[(k-1)*(S)+s,1]) for s=1:S) for k=1:(horzLen+1))
#         @objective(tM,Min, constFun1(-z,Vz[:,1])+constFun2(z,Vz[:,1]))
#     end
#     @constraint(tM,tempCon1,xt[1,1]==tauP*xt0+gammaP*80*sum((2*m+1)*z[m+1,1] for m=0:S-1)+rhoP*363)
#     @constraint(tM,tempCon2[k=1:horzLen],xt[k+1,1]==tauP*xt[k,1]+gammaP*80*sum((2*m+1)*z[k*S+(m+1),1] for m=0:S-1)+rhoP*363)
#     @constraint(tM,upperTCon,xt.<=393)
#     @constraint(tM,xt.>=0)
#     @constraint(tM,z.>=0)
#     @constraint(tM,z.<=80)
#     status = solve(tM)
# end


function timeCheckModel(opt)

    S=50
    horzLen=283
    N=30
    Un=20*0.8*rand(Truncated(Normal(0), 0, 1), N*(horzLen+1))
    Z=4000*rand(Truncated(Normal(0), 0, 1), S*(horzLen+1))
    Gn=-5*rand(Truncated(Normal(0), 0, 1), 8520)
    Gz=1300*rand(Truncated(Normal(0), 0, 1), 14200)

    cM = Model(solver = GurobiSolver())
    @variable(cM,dUn[1:(N)*(horzLen+1)])
    @variable(cM,dZ[1:(S)*(horzLen+1)])
    coupledObj1(deltaY,gi)=1/2*deltaY'*deltaY+gi'*deltaY
    function coupledObj2(deltaY,gi)
        l1=length(deltaY)
        expp=0
        for ind=1:l1
            expp=expp+deltaY[ind,1]*deltaY[ind,1]+gi[ind,1]*deltaY[ind,1]
        end
        return expp
    end

    coupledObj3a(deltaY)=1/2*deltaY'*deltaY
    coupledObj3b(deltaY,gi)=gi'*deltaY


    if opt==1
        @objective(cM,Min, append!(sum(coupledObj1(dUn[collect(n:N:length(dUn[:,1])),1],Gn[collect(n:N:length(Gn[:,1])),1]) for n=1:N),
                                coupledObj1(dZ,Gz[:,1])))
    elseif opt==2

        @objective(cM,Min, append!(sum(coupledObj2(dUn[collect(n:N:length(dUn[:,1])),1],Gn[collect(n:N:length(Gn[:,1])),1]) for n=1:N),
                                coupledObj2(dZ,Gz[:,1])))

    elseif opt==3
        sumObj=sum(coupledObj3a(dUn[collect(n:N:length(dUn[:,1])),1]) for n=1:N)+
        sum(coupledObj3b(dUn[collect(n:N:length(dUn[:,1])),1],Gn[collect(n:N:length(Gn[:,1])),1]) for n=1:N)+
        coupledObj3a(dZ)+coupledObj3b(dZ,Gz[:,1])
        @objective(cM,Min, sumObj)
    end

    @constraint(cM,currCon[k=1:horzLen+1],0==-sum(Un[(k-1)*(N)+n,1]+dUn[(k-1)*(N)+n,1] for n=1:N)-770+
                                             sum(Z[(k-1)*(S)+s,1]+dZ[(k-1)*(S)+s,1] for s=1:S))
 	@constraint(cM,(Z[:,1]+dZ).>=0)
    @constraint(cM,(Z[:,1]+dZ).<=80)
	@constraint(cM,(Un[:,1]+dUn).<=20)
	@constraint(cM,(Un[:,1]+dUn).>=0)
    status = solve(cM)
end

Profile.clear()
timeCheckModel(3)
@profile timeCheckModel(3)
Juno.profiler()
