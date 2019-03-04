# helper functions for EVC code
using DataFrames
using FileIO
using Dates
using LaTeXStrings
# using Gadfly
# using Cairo #for png output
# using Fontconfig

function rebuildEVS(evS)

    a=zeros(N,1)
    b_kWh=zeros(N,1)
    Tamb_raw=evS.Tamb
    newS=scenarioStruct(evS.N,evS.Ts,evS.K1,evS.K2,evS.K,evS.S,evS.ItotalMax,evS.deltaI,evS.Tmax,evS.imin,evS.imax,
                        a,b_kWh,evS.ηP,evS.τP,evS.ρP,evS.γP,evS.s0,evS.t0,evS.Snmin,evS.Kn,evS.iD_pred,evS.iD_actual,
                        evS.Tamb,Tamb_raw,evS.Qsi,evS.Ri,evS.β)
    #
    # newS=scenarioStruct(N=evS.N,Ts=evS.Ts,K1=evS.K1,K2=evS.K2,K=evS.K,S=evS.S,ItotalMax=evS.ItotalMax,deltaI=evS.deltaI,
    #                     Tmax=evS.Tmax,imin=evS.imin,imax=evS.imax,a=evS.a,b_kWh=evS.b_kWh,ηP=evS.ηP,τP=evS.τP,ρP=evS.ρP,
    #                     γP=evS.γP,s0=evS.s0,t0=evS.t0,Snmin=evS.Snmin,Kn=evS.Kn,iD_pred=evS.iD_pred,iD_actual=evS.iD_actual,
    #                     Tamb=evS.Tamb,Tamb_raw=evS.Tamb_raw,Qsi=evS.Qsi,Ri=evS.Ri,β=evS.β)

    if newS==true
        save(path*"EVCscenarioN$(N).jld2","evScenario",newS)
    end

    return newS
end

function saveRun(path::String, filename::String, time::Float64, solution; cSave=centralLogStruct(logLength=1,horzLen=1,N=1,S=1),
    convMetrics=convMetricsStruct(maxIt=1,logLength=1))

    saveFile=path*filename*".jld2"
    save(saveFile,"runTime", time, "solution", solution,"centralLog",cSave,"convMetrics", convMetrics)
end

function clips()
    t=clipboard()
    return t[2:length(t)-1]
end

function clipr()
    t=clipboard()
    loadF=load(t[2:length(t)-1])

    return loadF
end

function pubPlot(p;upscale=8,thickscale=2,dpi=300,sizeWH=(800,600))
    plot!(p,size=(sizeWH[1]*upscale,sizeWH[2]*upscale),thickness_scaling=thickscale*upscale,dpi=dpi)
end

function renameLabel(label)
    if occursin("PEM",label)
        new="PEM"
    elseif occursin("ADMM",label)
        new="ADMM"
    elseif occursin("ALADIN",label)
        new="ALADIN"
    elseif occursin("dual",label)
        new="dual"
    elseif occursin("central_NL",label)
        new="NL"
    elseif occursin("central",label)
        new="PWL"
    end
    return new
end

function checkDesiredStates(Sn,Kn,Snmin)
    epsilon=1e-2
    flag=true
    N=length(Kn)

    #if size(Sn)[2]>1 #columns are for each vehicle
    if length(size(Sn))>1
        for n=1:N
            if (Sn[Kn[n],n]-Snmin[n])<-epsilon
                flag=false
                #println(n)
            end
        end
    else #all vehicles in one columns
        for n=1:N
            if (Sn[(Kn[n]-1)*N+n]-Snmin[n])<-epsilon
                flag=false
            end
        end
    end

    return flag
end
