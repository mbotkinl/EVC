using Mosek
using JuMP

A=[3;2]
b=5
c=[12;7]
m1 = Model(solver = MosekSolver())
@variable(m1, x[1:2] )
@variable(m1, t)
@objective(m1, Min, t)
@constraint(m1, soc, norm( x[i] for i=1:2 ) <= t)
@constraint(m1,eqCon,A'*x==b)
status = solve(m1)
x1=getvalue(x)
s1=getdual(soc)
e1=-getdual(eqCon)
p1=getobjectivevalue(m1)

m2 = Model(solver = MosekSolver())
@variable(m2, x[1:2])
@objective(m2, Min, sum(c[i]*x[i]^2 for i=1:2 ))
@constraint(m2,eqCon,A'*x==b)
status = solve(m2)
x2=getvalue(x)
e2=-getdual(eqCon)
p2=getobjectivevalue(m2)


all(round.(abs.(x1-x2),digits=6).==0)
round(abs(p1^2-p2),digits=6)==0
round(abs(e1-e2),digits=6)==0



#dual of SOCP
tM=Model(solver = MosekSolver())
@variable(tM, v)
@objective(tM,Max,5v)
@constraint(tM,13v^2<=1)
#@constraint(tM,soc,norm(Av)<=1)
status = solve(tM)
getvalue(v)
getobjectivevalue(tM)





m3 = Model(solver = MosekSolver())
@variable(m3, x[1:2] )
@variable(m3, t)
@objective(m3, Min, t)
objCon=[2*sqrt(c[1])x[1];2*sqrt(c[2])*x[2];t-1]
@constraint(m3, soc, norm( objCon) <= t+1)
@constraint(m3,eqCon,A'*x==b)
status = solve(m3)
x3=getvalue(x)
s3=getdual(soc)
e3=-getdual(eqCon)
p3=getobjectivevalue(m3)
