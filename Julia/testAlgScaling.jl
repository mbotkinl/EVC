#Test how algorithms scale
# run EVC_init before hand

using Gurobi
include("C://Users//micah//Documents//uvm//Research//EVC code//Julia//functions//funEVCpwl.jl")

fname = "central_N$(N)"
println("Reading in Central Sim")
loadF=load(path*fname*".jld2")
cSol=loadF["solution"]
cSave=loadF["centralLog"]

Ns=[90 100 110]
#Ns=[100]

stepI=1
t0=evS.t0
s0=evS.s0
K=evS.K
horzLen=min(evS.K1,K-stepI)
S=evS.S
roundSigFigs=30
maxIt=5000
silent=true
saveLogInd=[]
p=1

compareSpeedTable = DataFrame(name=Int[],dualTime=Int[],dualConv=Int[],admmTime=Int[],
                              admmConv=Int[],aladTime=Int[],aladConv=Int[])

for i=1:length(Ns)
    @printf "%s: calc %g of %g....\n" Dates.format(Dates.now(),"HH:MM:SS") i length(Ns)
    global N=Ns[i]
    global evS

    Qsi=repeat(evS.Qsi,1,evS.K)
    Ri=repeat(evS.Ri,1,evS.K)
    evS=scenarioStruct(N,evS.Ts,evS.K1,evS.K2,evS.K,evS.S,evS.ItotalMax,evS.deltaI,evS.Tmax,evS.imin,evS.imax,
                        evS.a,evS.b_kWh,evS.ηP,evS.τP,evS.ρP,evS.γP,evS.s0,evS.t0,evS.Snmin,evS.Kn,evS.iD_pred,evS.iD_actual,
                        evS.Tamb,evS.Tamb_raw,Qsi,Ri,evS.β)


    dSol=solutionStruct(K=K,N=N,S=S)
    dCM=convMetricsStruct(maxIt=maxIt,logLength=length(saveLogInd))

    # decentral
    global prevLam=5e5*ones(evS.K1+1,1)
    global alphaDivRate=2
    global alpha0 = 1e4 #for kA
    dLog=itLogPWL(horzLen=horzLen,N=N,S=S)
    p=1
    timeStart=now()
    while (p<=maxIt && round(now()-timeStart,Second)<=Dates.Second(9/10*evS.Ts))
        #global p
        if p==1
            itLam=prevLam
            #global alpha0=max(min(maximum(prevLam)/5,1e6),1e-3)
        else
            itLam=round.(dLog.Lam[:,(p-1)],sigdigits=roundSigFigs)
        end
        cFlag=runEVDualIt(p,stepI,evS,itLam,dLog,dCM,dSol,cSave,roundSigFigs,silent)
        global convIt=p
        if cFlag
            break
        end
        p+=1
    end
    dualTime = round(now()-timeStart,Second)
    dualConv= convIt
    @printf "Dual Tlim hit: %s \n" any(abs.(dLog.Tpred[:,convIt] .- evS.Tmax) .<1e-3)

    #ADMM
    global prevLam=2e4*ones(evS.K1+1,1)
    global prevVz=-evS.deltaI/2*ones(evS.S*(evS.K1+1),1)
    global prevVu=.02*ones(evS.N*(evS.K1+1),1)
    global ρADMMp = 2e5
    global ρDivRate=1.02
    dLogadmm=itLogPWL(horzLen=horzLen,N=N,S=S)
    p=1
    timeStart=now()
    while (p<=maxIt && round(now()-timeStart,Second)<=Dates.Second(9/10*evS.Ts))
        #global p
        if p==1
            itLam=prevLam
            itVu=prevVu
            itVz=prevVz
            itρ=ρADMMp
        else
            itLam=round.(dLogadmm.Lam[:,(p-1)],sigdigits=roundSigFigs)
            itVu=round.(dLogadmm.Vu[:,(p-1)],sigdigits=roundSigFigs)
            itVz=round.(dLogadmm.Vz[:,(p-1)],sigdigits=roundSigFigs)
            itρ=round.(dLogadmm.itUpdate[1,(p-1)],sigdigits=roundSigFigs)
        end
        cFlag=runEVADMMIt(p,stepI,evS,itLam,itVu,itVz,itρ,dLogadmm,dCM,dSol,cSave,roundSigFigs,silent)
        global convIt=p
        if cFlag
            break
        end
        p+=1
    end
    admmTime = round(now()-timeStart,Second)
    admmConv= convIt
    @printf "ADMM Tlim hit: %s \n" any(abs.(dLogadmm.Tpred[:,convIt] .- evS.Tmax) .<1e-3)

    #ALADIN
    global prevLam=1e5*ones(evS.K1+1,1)
    global prevVt=evS.Tmax*ones(evS.K1+1,1)
    global prevVz=.01*ones(evS.S*(evS.K1+1),1)
    global prevVu=.01*ones(evS.N*(evS.K1+1),1)
    global prevVs=.5*ones(evS.N*(evS.K1+1),1)
    if eqForm
        global ρALADp = 1e6
        global ρRate=1.15
    else
        global ρALADp = 1
        global ρRate=1.1
    end
    global reg_weight=1e-3
    global reg=false
    dLogalad=itLogPWL(horzLen=horzLen,N=N,S=S)
    p=1
    timeStart=now()
    while (p<=maxIt && round(now()-timeStart,Second)<=Dates.Second(9/10*evS.Ts))
        #global p
        if p==1
            itLam=prevLam
            itVu=prevVu
            itVz=prevVz
            itVs=prevVs
            itVt=prevVt
            itρ=ρALADp
        else
            itLam=round.(dLogalad.Lam[:,(p-1)],sigdigits=roundSigFigs)
            itVu=round.(dLogalad.Vu[:,(p-1)],sigdigits=roundSigFigs)
            itVz=round.(dLogalad.Vz[:,(p-1)],sigdigits=roundSigFigs)
            itVs=round.(dLogalad.Vs[:,(p-1)],sigdigits=roundSigFigs)
            itVt=round.(dLogalad.Vt[:,(p-1)],sigdigits=roundSigFigs)
            itρ=round.(dLogalad.itUpdate[1,(p-1)],sigdigits=roundSigFigs)
        end
        cFlag=runEVALADIt(p,stepI,evS,itLam,itVu,itVz,itVs,itVt,itρ,dLogalad,dCM,dSol,cSave,eqForm,roundSigFigs,silent)
        global convIt=p
        if cFlag
            break
        end
        p+=1
    end
    aladTime = round(now()-timeStart,Second)
    aladConv= convIt
    @printf "ALAD Tlim hit: %s \n" any(abs.(dLogalad.Tpred[:,convIt] .- evS.Tmax) .<1e-3)

    stats = [N dualTime.value dualConv admmTime.value admmConv aladTime.value aladConv]
    push!(compareSpeedTable,stats)
end
