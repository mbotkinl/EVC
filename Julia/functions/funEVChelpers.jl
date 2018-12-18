# helper functions for EVC code
using DataFrames
using FileIO
using Dates

# using Gadfly
# using Cairo #for png output
# using Fontconfig

function saveRun(path::String, filename::String, time::Float64, scenario, solution, convMetrics=convMetricsStruct(), convIt=1)
    save(path*filename*".jld2","runTime", time, "scenario", scenario, "solution", solution,
    "convMetrics", convMetrics, "convIt", convIt)
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

function readRuns(path)
    files = filter(x->occursin(".jld2",x), readdir(path))
    files = filter(x->occursin("_",x), files) # avoid evScenario
    noLimFile = filter(x->occursin("noLim",x), files)
    cFile = filter(x->occursin("central",x), files)
    dFiles= setdiff(files,cFile)

    if length(cFile)>0
        cFile=setdiff(cFile,noLimFile)
        cRun=load(path*cFile[1])
    else
        cRun=Nothing
    end

    if length(noLimFile)>0
        noLim=load(path*noLimFile[1])
    else
        noLim=Nothing
    end

    runs=Dict{String,Any}()
    for ii=1:length(dFiles)
        runs[dFiles[ii]]=load(path*dFiles[ii])
    end
    return cRun, runs, noLim
end

function compareRunsGraph(runs, cRun, noLim, saveF::Bool, lowRes::Bool)
    runNames=collect(keys(runs))
    cSol=cRun["solution"]
    numIt=size(runs[runNames[1]]["convMetrics"].lam)[1]
    Klen=size(cSol.Tactual)[1]
    P=length(runNames)
    evS=cRun["scenario"]
    N=evS.N

    objPerc = zeros(numIt,P)
    objConv = zeros(numIt,P)
    lamConv = zeros(numIt,P)
    lamRMSE = zeros(numIt,P)
    lamInfNorm = zeros(numIt,P)
    Lam = zeros(Klen,P)
    T = zeros(Klen,P)
    Sn = zeros(Klen*N,P)
    uSum = zeros(Klen,P)
    snSum = zeros(Klen,P)
    snAvg = zeros(Klen,P)

    for i in 1:length(runNames)
        println(runNames[i])
        runI=runs[runNames[i]]
        cIt=runI["convIt"]
        if typeof(runI["solution"]) == solutionStruct #PEM
            lamRMSE[:,i].=NaN
            lamInfNorm[:,i].=NaN
            #objPerc[1:numIt,i].=abs.(cSol.objVal.-runI["solution"].objVal)/cSol.objVal*100
            Lam[:,i]=runI["solution"].lamCoupl
            T[:,i]=runI["solution"].Tactual[:,cIt]
        else
            lamRMSE[:,i]=[sqrt(1/Klen*sum((runI["solution"].Lam[k,it]-cSol.lamCoupl[k])^2/abs(cSol.lamCoupl[k]) for k=1:numIt)) for it=1:numIt]
            lamInfNorm[:,i]=[maximum(abs.(runI["solution"].Lam[:,it]-cSol.lamCoupl)) for it=1:numIt]
            objPerc[1:numIt,i]=abs.(cSol.objVal.-runI["solution"].objVal[1:numIt]')/cSol.objVal*100
            Lam[:,i]=runI["solution"].Lam[:,cIt]

            if typeof(runI["solution"])==itLogPWL
                T[:,i]=runI["solution"].Tactual[:,cIt] #for PWL
            else
                T[:,i]=runI["solution"].Xt[:,cIt] # for NL
            end
            objConv[:,i]=runI["convMetrics"].obj
            lamConv[:,i]=runI["convMetrics"].lam
        end
        uSum[:,i]=runI["solution"].uSum[:,cIt]
        if size(runI["solution"].Sn)[1]>Klen+1
            Sn[:,i]=runI["solution"].Sn[:,cIt]
            snSum[:,i]=[sum(Sn[N*(k-1)+n,i] for n=1:N) for k in 1:Klen]# sum up across N
            snAvg[:,i]=[mean(Sn[N*(k-1)+n,i] for n=1:N) for k in 1:Klen]# mean across N
        else
            Sn[:,i]=runI["solution"].Sn'[:]
            snSum[:,i]=[sum(runI["solution"].Sn[k,:]) for k in 1:Klen]# sum up across N
            snAvg[:,i]=[mean(runI["solution"].Sn[k,:]) for k in 1:Klen]# mean across N
        end

        # fill in NaN for greater than convIt
        objPerc[cIt:numIt,i].=NaN
        lamRMSE[cIt:numIt,i].=NaN
        lamInfNorm[cIt:numIt,i].=NaN
        objConv[cIt:numIt,i].=NaN
        lamConv[cIt:numIt,i].=NaN
    end

    plotLabels=permutedims(runNames)
    allColors=get_color_palette(:auto, plot_color(:white), P+1)
    plotColors=allColors[1:P]'

    stT1=Time(20,0)
    endT1=Time(23,59)
    stT2=Time(0,0)
    endT2=Time(10,0)
    Xlabels=vcat(collect(stT1:Dates.Second(round(evS.Ts)):endT1),collect(stT2:Dates.Second(round(evS.Ts)):endT2))
    #Xlabels=vcat(collect(stT1:Dates.Minute(3):endT1),collect(stT2:Dates.Minute(3):endT2))
    xticks=(1:40:Klen,Dates.format.(Xlabels[1:40:Klen],"HH:MM"))

    #Iteration plots
    lamConvPlot=plot(lamConv,xlabel="",ylabel="2-norm Lambda gap",labels=plotLabels,yscale=:log10,
                    seriescolor=plotColors,legend=false,xlims=(0,numIt))
    lamRMSEPlot=plot(lamRMSE,xlabel="",ylabel="Relative RMSE Lambda Gap",labels=plotLabels,yscale=:log10,
                    seriescolor=plotColors,xlims=(0,numIt))
    lamInfNormPlot=plot(lamInfNorm,xlabel="",ylabel="Max Lambda Gap",labels=plotLabels,yscale=:log10,
                        seriescolor=plotColors,legend=false,xlims=(0,numIt))
    objNormPlot=plot(objConv,xlabel="",ylabel="Objective Value Magintude Gap",labels=plotLabels,yscale=:log10,
                    seriescolor=plotColors',legend=false,xlims=(0,numIt))
    objPercPlot=plot(objPerc,xlabel="Iteration",ylabel="Objective Value Percentage Gap",labels=plotLabels,yscale=:log10,
                    seriescolor=plotColors,legend=false,xlims=(0,numIt))

    #Time plots
    tempPlot=plot(1:Klen,cSol.Tactual*1000,label="",seriescolor=:black,linewidth=4,linealpha=0.25,xlims=(0,Klen),
                    xlabel="",ylabel="Temp (K)",xticks=xticks)
    plot!(tempPlot,1:Klen,evS.Tmax*ones(Klen)*1000,label="XFRM Limit",line=(:dash,:red))
    plot!(tempPlot,T*1000,labels="",seriescolor=plotColors)
    if noLim !=nothing
        plot!(tempPlot,1:Klen,noLim["solution"].Tactual*1000,label="",seriescolor=allColors[P+1])
    end

    uSumPlot=plot(1:Klen,cSol.uSum,xlabel="",ylabel="Current Sum (kA)",xlims=(0,Klen),labels="Central",
                  seriescolor=:black,linewidth=4,linealpha=0.25,xticks=xticks)
    plot!(uSumPlot,uSum,labels=plotLabels,seriescolor=plotColors,legend=false)
    if noLim !=nothing
        plot!(uSumPlot,1:Klen,noLim["solution"].uSum,label="Uncoordinated",seriescolor=allColors[P+1])
    end

    iD=evS.iD_actual[1:Klen].+10
    loadPlot=plot(1:Klen,cSol.uSum+iD,xlabel="",ylabel="Total Load (kA)",xlims=(0,Klen),labels="",
                  seriescolor=:black,linewidth=4,linealpha=0.25,xticks=xticks)
    plot!(loadPlot,1:Klen,iD,label="Background Demand",line=(:dash))
    plot!(loadPlot,uSum.+iD,labels="",seriescolor=plotColors)
    if noLim !=nothing
        plot!(loadPlot,1:Klen,noLim["solution"].uSum+iD,label="",seriescolor=allColors[P+1])
    end

    target=zeros(Klen)
    for k in 1:Klen
        ind = evS.Kn.<=k
        target[k]=sum(evS.Snmin[ind])
    end

    snSumCentral=[sum(cSol.Sn[N*(k-1)+n,1] for n=1:N) for k in 1:Klen]# sum up across N
    snSumPlot=plot(1:Klen,snSumCentral,xlabel="",ylabel="SOC Sum",xlims=(0,Klen),labels="",
                   seriescolor=:black,linewidth=4,linealpha=0.25,xticks=xticks)
    plot!(snSumPlot,1:Klen,target,label="SOC Target",line=(:dash))
    plot!(snSumPlot,snSum,labels="",seriescolor=plotColors)
    if noLim !=nothing
        snSumNoLim=[sum(noLim["solution"].Sn[N*(k-1)+n,1] for n=1:N) for k in 1:Klen]# sum up across N
        plot!(snSumPlot,1:Klen,snSumNoLim,label="",seriescolor=allColors[P+1])
    end


    snAvgCentral=[mean(cSol.Sn[N*(k-1)+n,1] for n=1:N) for k in 1:Klen]# mean across N
    snAvgPlot=plot(1:Klen,snAvgCentral,xlabel="",ylabel="Avg.SOC",xlims=(0,Klen),labels=plotLabels,
                   seriescolor=:black,linewidth=4,linealpha=0.25,xticks=xticks)
    plot!(snAvgPlot,snAvg,labels=plotLabels,seriescolor=plotColors)
    if noLim !=nothing
        snAvgNoLim=[mean(noLim["solution"].Sn[N*(k-1)+n,1] for n=1:N) for k in 1:Klen]# sum up across N
        plot!(snAvgPlot,1:Klen,snAvgNoLim,label="Uncoordinated",seriescolor=allColors[P+1])
    end

    lamPlot=plot(1:Klen,cSol.lamCoupl,xlabel="Time",ylabel=raw"Lambda ($/kA)",xlims=(0,Klen),labels="Central",
                   seriescolor=:black,linewidth=4,linealpha=0.25,xticks=xticks)
    plot!(lamPlot,Lam,labels=plotLabels,seriescolor=plotColors)
    if noLim !=nothing
        plot!(lamPlot,1:Klen,noLim["solution"].lamCoupl,label="Uncoordinated",seriescolor=allColors[P+1])
    end
    # R plot***
    Rmax=zeros(Klen,P)
    Ravg=zeros(Klen,P)
    for p=1:P
        for k=1:Klen
            R=[max((evS.Snmin[n,1]-Sn[N*(k-1)+n,p]),0)./(evS.Kn[n,1]-k) for n=1:N]
            Rmax[k,p]=maximum(R)
            Ravg[k,p]=mean(R)
        end
    end


    Rmax=plot(Rmax,xlabel="Time",ylabel="R Max",xlims=(0,Klen),labels=plotLabels,seriescolor=plotColors,legend=false)


    resPlot=plot(tempPlot,loadPlot,snSumPlot,lamPlot,layout=(4,1))
    if lowRes
        pubPlot(resPlot,thickscale=0.4,sizeWH=(400,300),dpi=40)
    else
        pubPlot(resPlot,thickscale=0.8,sizeWH=(800,600),dpi=100)
    end

    convPlot=plot(lamRMSEPlot,lamInfNormPlot,objPercPlot,layout=(3,1))
    if lowRes
        pubPlot(convPlot,thickscale=0.4,sizeWH=(400,300),dpi=40)
    else
        pubPlot(convPlot,thickscale=0.8,sizeWH=(800,600),dpi=100)
    end

    # p=plot(lamRMSEPlot, tempPlot,lamInfNormPlot, loadPlot, objPercPlot,lamPlot,layout=(3,2))
    # if lowRes
    #     pubPlot(p,thickscale=0.4,sizeWH=(500,300),dpi=40)
    # else
    #     pubPlot(p,thickscale=0.8,sizeWH=(1000,600),dpi=100)
    # end

    if saveF savefig(resPlot,path*"resPlot.png") end
    if saveF savefig(convPlot,path*"convPlot.png") end

    return resPlot, convPlot
end

function compareRunsTable(runs)
    # compareTable = DataFrame(name=String[],time=Float64[],cLamDiff=Float64[],lamDiff=Float64[],
    # cObjDiff=Float64[],objDiff=Float64[])
    compareTable = DataFrame(name=String[],timeTotal=Float64[],timePerIt=Float64[],convIt=Int64[])
    for keyI in keys(runs)
        println(keyI)
        local loadF=runs[keyI]
        local timeT=loadF["runTime"]
        #cm=loadF["convMetrics"]
        local convIt=loadF["convIt"]
        #ind= if convIt>0 convIt-1 else length(cm.lam) end
        #stats = [key timeT minimum(cm.lam[1:ind]) cm.lam[ind-1]-cm.lam[ind-2] minimum(cm.obj[1:ind]) cm.obj[ind-1]-cm.obj[ind-2]]
        local stats = [keyI timeT timeT/convIt convIt]
        push!(compareTable,stats)
    end

    return compareTable
end


function checkDesiredStates(Sn,Kn,Snmin)
    epsilon=1e-3
    flag=true
    N=length(Kn)

    #if size(Sn)[2]>1 #columns are for each vehicle
    if length(size(Sn))>1
        for n=1:N
            if (Sn[Kn[n],n]-Snmin[n])<-epsilon
                flag=false
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
