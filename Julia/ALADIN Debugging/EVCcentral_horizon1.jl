
sn0=s0
xt0=T0

#set up Problem
m = Model(solver = IpoptSolver())

#u w and z are one index ahead of x. i.e the x[k+1]=x[k]+eta*u[k+1]
@variable(m,u[1:N])
@variable(m,sn[1:N])
@variable(m,xt)
@variable(m,z[1:S])
objFun(sn,u)=sum((sn[n,1]-1)^2*Qsi[n,1]     for n=1:N) +
			 sum((u[n,1])^2*Ri[n,1]           for n=1:N)
@objective(m,Min, objFun(sn,u))
#@constraint(m,stateCon1[n=1:N],sn[n,1]/etaP[n,1]==sn0[n,1]/etaP[n,1]+u[n,1])
@constraint(m,stateCon1,sn[1:N,1].==sn0[1:N,1]+etaP[:,1].*u[1:N,1])
#@constraint(m,tempCon1,xt/gammaP==tauP*xt0/gammaP+deltaI*sum((2*m+1)*z[m+1,1] for m=0:S-1)+rhoP*w[2,1]/gammaP)
@constraint(m,tempCon1,xt==tauP*xt0+gammaP*deltaI*sum((2*m+1)*z[m+1,1] for m=0:S-1)+rhoP*w[2,1])
@constraint(m,currCon,sum(u[n] for n=1:N)-sum(z[s] for s=1:S)==-w[1])
@constraint(m,sn.<=1)
@constraint(m,sn.>=target)
@constraint(m,upperTCon,xt<=Tmax)
@constraint(m,xt>=0)
@constraint(m,upperCCon,u.<=imax)
@constraint(m,u.>=imin)
@constraint(m,z.>=0)
@constraint(m,z.<=deltaI)

println("solving....")
statusM = solve(m)
if statusM!=:Optimal
    return
else
	xtStar1=getvalue(xt)
	snStar1=getvalue(sn)
	uStar1=getvalue(u)
	zStar1=getvalue(z)
	fStar1=getobjectivevalue(m)
	lamTempStar1=-getdual(upperTCon)
	lamCurrStar1=-getdual(currCon)
end
