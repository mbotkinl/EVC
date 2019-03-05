#test discretization
T1=evS.K1*3600
Ts=evS.Ts
T0=evS.t0
t0=0
Tamb=evS.Tamb
#Tamb=370
iD=evS.iD_actual
#iD=0
# u=cSolnl.Un
horzLen=evS.K1


#
# m   = 2000                             # transformer mass in kg
# C   = 450*m                            # heat cap. thermal mass J/K ----- spec. heat cap. of C = {carbon steel, iron, veg. oil} = {490, 450, 1670} J/(kg*K)
# Rh   = 1070e-4/(35*5*(m/7870)^(2/3))   # heat outflow resistance K/W : R = 0.1073 (K*m^2/W)/(A_s), rule of thumb calculation
# Rw  = 1
# Vac = 240                              # PEV battery rms voltage --- V [used in PEV kW -> kA conversion]
# Vtf = 8320                             # distr-level transformer rms voltage --- V [used in inelastic kW -> kA conv]
# Ntf   = Vtf/Vac                        # pole-top transformer turns ratio
#
# beta=1/(C*Rh)

#t=range(t0,stop=T1,step=1)
timeSteps=1:horzLen+1


# using DifferentialEquations
# ff(u,p,t) = 1/C*((Rw/Ntf)*(1)^2-1/Rh*(u-370))
# #ff(u,p,t) =-1/Rh*(u-370)
#
# u0=370.0
# tspan = (1,(horzLen+1)*Ts)
# prob = ODEProblem(ff,u0,tspan)
# sol = DifferentialEquations.solve(prob,Tsit5(),reltol=1e-12,abstol=1e-12,saveat=Ts)
# actualT=sol.u
#
# plot(actualT)
#



#continous time
using QuadGK
intR=zeros(length(timeSteps),1)
for i=1:length(timeSteps)
    ind=timeSteps[i]
    t=ind*Ts
    f=x -> exp(-beta*(t-x))*(Rw/(C*Ntf)*(1000)^2+beta*370)
    rr,tol=quadgk(f, t0, t, rtol=1e-8)
    intR[i]=rr
end
actualT = exp.(-beta*(timeSteps*Ts.-t0))*T0.+intR

stepI=1
iTotal=zeros(horzLen+1,1)
for k=1:horzLen+1
    # iTotal[k,1]=sum(u[(k-1)*N+n,1] for n=1:N) + iD[stepI+(k-1),1]
    #iTotal[k,1]=2*cos(k*Ts)^2
    iTotal[k,1]=1
end


#old discretization
τP1 = 1 - Ts/(Rh*C)
ρP1 = 1 - τP1            # no units, ambient-to-temp param: 1/RC
γP1 = Ts*Rw/(C*Ntf)*1000^2     # K/kW, ohmic losses-to-temp parameter
discT1=zeros(length(timeSteps),1)
# discT1[1,1]=τP1*T0+γP1*iTotal[1,1]^2+ρP1*Tamb[stepI,1] #fix for mpc
discT1[1,1]=τP1*T0+γP1*iTotal[1,1]^2+ρP1*370
for k=1:horzLen
    #discT1[k+1,1]=τP1*discT1[k,1]+γP1*iTotal[k+1,1]^2+ρP1*Tamb[stepI+k,1]  #fix for mpc
    discT1[k+1,1]=τP1*discT1[k,1]+γP1*iTotal[k+1,1]^2+ρP1*370  #fix for mpc
end

#new discretization
τP2 = exp(- Ts/(Rh*C))
ρP2 = 1 - τP2            # no units, ambient-to-temp param: 1/RC
γP2 = Rh*Rw/(Ntf)*ρP2*1000^2    # K/kW, ohmic losses-to-temp parameter
discT2=zeros(length(timeSteps),1)
# discT2[1,1]=τP2*T0+γP2*iTotal[1,1]^2+ρP2*Tamb[stepI,1] #fix for mpc
discT2[1,1]=τP2*T0+γP2*iTotal[1,1]^2+ρP2*370 #fix for mpc
for k=1:horzLen
    #discT2[k+1,1]=τP2*discT2[k,1]+γP2*iTotal[k+1,1]^2+ρP2*Tamb[stepI+k,1]  #fix for mpc
    discT2[k+1,1]=τP2*discT2[k,1]+γP2*iTotal[k+1,1]^2+ρP2*370  #fix for mpc
end

p1=plot(hcat(discT1,discT2))


p1=plot(hcat(actualT,discT1,discT2),labels=["CT" "D1" "D2"])


function rms(A)
   s = 0.0
   for a in A
      s += a*a
   end
   return sqrt(s / length(A))
end

 rms(actualT-discT1)

 rms(actualT-discT2)
