# helper functions for EVC code
using DataFrames
using FileIO
using Dates
using LaTeXStrings
# using Gadfly
# using Cairo #for png output
# using Fontconfig

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

function readRuns(path)
    files = filter(x->occursin(".jld2",x), readdir(path))
    evFile = filter(x->occursin("scenario",x), files)
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

    if length(evFile)>0
        evS=load(path*evFile[1])["evScenario"]
    else
        evS=Nothing
    end

    runs=Dict{String,Any}()
    for ii=1:length(dFiles)
        runs[dFiles[ii]]=load(path*dFiles[ii])
    end
    return cRun, runs, noLim, evS
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
    end
    return new
end

function compareRunsGraph(runs, cRun, noLim, saveF::Bool, lowRes::Bool)
    runNames=collect(keys(runs))
    cSol=cRun["solution"]
    numIt=size(runs[runNames[1]]["convMetrics"].lam)[1]
    Klen=size(cSol.Tactual)[1]
    P=length(runNames)
    N=evS.N

    Lam = zeros(Klen,P)
    T = zeros(Klen,P)
    Sn = zeros(Klen*N,P)
    uSum = zeros(Klen,P)
    snSum = zeros(Klen,P)
    snAvg = zeros(Klen,P)
    convIt = zeros(Klen,P)

    for i in 1:length(runNames)
        println(runNames[i])
        runI=runs[runNames[i]]
        convIt[:,i]=runI["solution"].convIt
        Lam[:,i]=runI["solution"].lamCoupl
        T[:,i]=runI["solution"].Tactual[:,1]
        uSum[:,i]=runI["solution"].uSum[:,1]
        if size(runI["solution"].Sn)[1]>Klen+1
            Sn[:,i]=runI["solution"].Sn[:,cIt]
            snSum[:,i]=[sum(Sn[N*(k-1)+n,i] for n=1:N) for k in 1:Klen]# sum up across N
            snAvg[:,i]=[mean(Sn[N*(k-1)+n,i] for n=1:N) for k in 1:Klen]# mean across N
        else
            Sn[:,i]=runI["solution"].Sn'[:]
            snSum[:,i]=[sum(runI["solution"].Sn[k,:]) for k in 1:Klen]# sum up across N
            snAvg[:,i]=[mean(runI["solution"].Sn[k,:]) for k in 1:Klen]# mean across N
        end
    end

    plotLabels=copy(permutedims(runNames))

    for i=1:length(plotLabels)
        plotLabels[i]=renameLabel(plotLabels[i])
    end

    allColors=get_color_palette(:auto, plot_color(:white), P+1)
    plotColors=allColors[1:P]'

    stT1=Time(20,0)
    endT1=Time(23,59)
    stT2=Time(0,0)
    endT2=Time(10,0)
    Xlabels=vcat(collect(stT1:Dates.Second(round(evS.Ts)):endT1),collect(stT2:Dates.Second(round(evS.Ts)):endT2))
    #Xlabels=vcat(collect(stT1:Dates.Minute(3):endT1),collect(stT2:Dates.Minute(3):endT2))
    xticks=(1:40:Klen,Dates.format.(Xlabels[1:40:Klen],"HH:MM"))

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

    snSumCentral=sum(cSol.Sn,dims=2)# sum up across N
    snSumPlot=plot(snSumCentral,xlabel="",ylabel="SOC Sum",xlims=(0,Klen),labels="",
                   seriescolor=:black,linewidth=4,linealpha=0.25,xticks=xticks)
    plot!(snSumPlot,target,label="SOC Target",line=(:dash))
    plot!(snSumPlot,snSum,labels="",seriescolor=plotColors)
    if noLim !=nothing
        plot!(snSumPlot,sum(noLim["solution"].Sn,dims=2),label="",seriescolor=allColors[P+1])
    end


    # snAvgCentral=[mean(cSol.Sn[N*(k-1)+n,1] for n=1:N) for k in 1:Klen]# mean across N
    # snAvgPlot=plot(1:Klen,snAvgCentral,xlabel="",ylabel="Avg.SOC",xlims=(0,Klen),labels=plotLabels,
    #                seriescolor=:black,linewidth=4,linealpha=0.25,xticks=xticks)
    # plot!(snAvgPlot,snAvg,labels=plotLabels,seriescolor=plotColors)
    # if noLim !=nothing
    #     snAvgNoLim=[mean(noLim["solution"].Sn[N*(k-1)+n,1] for n=1:N) for k in 1:Klen]# sum up across N
    #     plot!(snAvgPlot,1:Klen,snAvgNoLim,label="Uncoordinated",seriescolor=allColors[P+1])
    # end

    convItPlot=plot(convIt,xlabel="Time",ylabel=raw"Lambda ($/kA)",xlims=(0,Klen),xticks=xticks,labels=plotLabels)


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


    #resPlot=plot(tempPlot,loadPlot,snSumPlot,lamPlot,layout=(4,1))
    resPlot=plot(tempPlot,loadPlot,lamPlot,layout=(3,1))
    if lowRes
        pubPlot(resPlot,thickscale=0.4,sizeWH=(400,300),dpi=40)
    else
        pubPlot(resPlot,thickscale=0.8,sizeWH=(800,600),dpi=100)
    end

    if saveF savefig(resPlot,path*"resPlot.png") end

    return resPlot
end

function compareConvGraph_New(evS;ind=1)
    runNames=collect(keys(runs))
    K=evS.K
    N=evS.N

    #ignore PEM
    pemRun = filter(x->occursin("PEM",x), runNames)
    runNames=setdiff(runNames,pemRun)
    P=length(runNames)
    numIt=size(runs[runNames[1]]["convMetrics"].coupl1Norm)[1]

    convNames=["coupl1Norm","lamIt2Norm","objAbs","objPerc",
                "lam1Norm","lam2Norm","lamInfNorm",
                "un1Norm","un2Norm","unInfNorm",
                "t1Norm","t2Norm","tInfNorm",
                "z1Norm","z2Norm","zInfNorm",]
    convDict=Dict()
    for ii=1:length(convNames)
        convDict[convNames[ii]]=zeros(numIt,P)
        for i=1:length(runNames)
            cIt=Int(runs[runNames[i]]["convMetrics"].convIt[ind])
            temp=getfield(runs[runNames[i]]["convMetrics"],Symbol(convNames[ii]))
            convDict[convNames[ii]][:,i]=temp[:,ind]
            # fill in NaN for greater than convIt
            convDict[convNames[ii]][cIt:numIt,i].=NaN
        end
    end

    plotLabels=copy(permutedims(runNames))
    for i=1:length(plotLabels)
        plotLabels[i]=renameLabel(plotLabels[i])
    end

    #internal metrics
    couplPlot=plot(convDict["coupl1Norm"],legend=false,yscale=:log10,ylabel= L"||i_d[k+1]+\sum_{n=1}^N i_n[k+1]-\sum_{m=1}^M i_m^{PW}[k+1]||_1")
    dualPlot=plot(convDict["lamIt2Norm"],labels=plotLabels,yscale=:log10,xlabel="Iteration",ylabel=L"||\lambda^{(p)}-\lambda^{(p-1)}||_2")

    #external
    un1Plot=plot(convDict["un1Norm"],labels=plotLabels,yscale=:log10,xlabel="Iteration",ylabel=L"||u_n^{(p)}-u_n^{*}||_1")
    unInfPlot=plot(convDict["unInfNorm"],labels=plotLabels,yscale=:log10,xlabel="Iteration",ylabel=L"||u_n^{(p)}-u_n^{*}||_{\infty}")
end

function compareConvGraph()
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

    for i in 1:length(runNames)
        println(runNames[i])
        runI=runs[runNames[i]]

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

        # fill in NaN for greater than convIt
        objPerc[cIt:numIt,i].=NaN
        lamRMSE[cIt:numIt,i].=NaN
        lamInfNorm[cIt:numIt,i].=NaN
        objConv[cIt:numIt,i].=NaN
        lamConv[cIt:numIt,i].=NaN
    end

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

    convPlot=plot(lamRMSEPlot,lamInfNormPlot,objPercPlot,layout=(3,1))
    if lowRes
        pubPlot(convPlot,thickscale=0.4,sizeWH=(400,300),dpi=40)
    else
        pubPlot(convPlot,thickscale=0.8,sizeWH=(800,600),dpi=100)
    end

    if saveF savefig(convPlot,path*"convPlot.png") end
end

function compareRunsTable(runs,evS)
    # compareTable = DataFrame(name=String[],time=Float64[],cLamDiff=Float64[],lamDiff=Float64[],
    # cObjDiff=Float64[],objDiff=Float64[])
    compareTable = DataFrame(name=String[],timeTotal=Float64[],avgTimePerTs=Float64[],avgtimePerIt=Float64[],avgconvIt=Float64[],maxConvIt=Float64[])
    for keyI in keys(runs)
        println(keyI)
        loadF=runs[keyI]
        Klen=evS.K
        timeT=loadF["runTime"]
        convIt=loadF["solution"].convIt
        #cm=loadF["convMetrics"]
        #local convIt=loadF["convIt"]
        #ind= if convIt>0 convIt-1 else length(cm.lam) end
        #stats = [key timeT minimum(cm.lam[1:ind]) cm.lam[ind-1]-cm.lam[ind-2] minimum(cm.obj[1:ind]) cm.obj[ind-1]-cm.obj[ind-2]]
        avgTimeTs=timeT/Klen
        avgTime=timeT/max(sum(convIt),Klen)
        avgConv=max(sum(convIt),Klen)/Klen
        maxConv=maximum(convIt)
        stats = [keyI timeT avgTimeTs avgTime avgConv maxConv]
        push!(compareTable,stats)
    end

    return compareTable
end

function compareHubsGraph(runs, cRun, noLim, saveF::Bool, lowRes::Bool)
    runNames=collect(keys(runs))
    cSol=cRun["solution"]
    Klen=size(cSol.Tactual)[1]
    P=length(runNames)
    hubS=cRun["scenario"]
    H=hubS.H

    Lam = zeros(Klen,P)
    T = zeros(Klen,P)
    uSum = zeros(Klen,P)
    eSum = zeros(Klen,P)
    eAvg = zeros(Klen,P)
    convIt = zeros(Klen,P)

    for i in 1:length(runNames)
        println(runNames[i])
        runI=runs[runNames[i]]
        convIt[:,i]=runI["solution"].convIt
        Lam[:,i]=runI["solution"].Lam
        T[:,i]=runI["solution"].Tactual[:,1]
        uSum[:,i]=runI["solution"].uSum[:,1]
        eSum[:,i]=[sum(runI["solution"].E[k,:]) for k in 1:Klen]# sum up across N
        eAvg[:,i]=[mean(runI["solution"].E[k,:]) for k in 1:Klen]# mean across N
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

    #Time plots
    tempPlot=plot(1:Klen,cSol.Tactual*1000,label="",seriescolor=:black,linewidth=4,linealpha=0.25,xlims=(0,Klen),
                    xlabel="",ylabel="Temp (K)",xticks=xticks)
    plot!(tempPlot,1:Klen,hubS.Tmax*ones(Klen)*1000,label="XFRM Limit",line=(:dash,:red))
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

    convItPlot=plot(convIt,xlabel="Time",ylabel=raw"Lambda ($/kA)",xlims=(0,Klen),xticks=xticks,labels=plotLabels)


    lamPlot=plot(1:Klen,cSol.Lam,xlabel="Time",ylabel=raw"Lambda ($/kA)",xlims=(0,Klen),labels="Central",
                   seriescolor=:black,linewidth=4,linealpha=0.25,xticks=xticks)
    plot!(lamPlot,Lam,labels=plotLabels,seriescolor=plotColors)
    if noLim !=nothing
        plot!(lamPlot,1:Klen,noLim["solution"].Lam,label="Uncoordinated",seriescolor=allColors[P+1])
    end

    resPlot=plot(tempPlot,loadPlot,lamPlot,layout=(3,1))
    if lowRes
        pubPlot(resPlot,thickscale=0.4,sizeWH=(400,300),dpi=40)
    else
        pubPlot(resPlot,thickscale=0.8,sizeWH=(800,600),dpi=100)
    end

    if saveF savefig(resPlot,path*"resPlot.png") end

    return resPlot, convPlot
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

function calcPrivacy(N,K,S)
    aladMaxIt=10
    admmMaxIt=100
    dualMaxIt=1000

    bpf= 64 #bits per double
    bpi= 16 #bits per int
    bpb= 1#bits per bool

    pem=N*bpi # 1 Req Int per EV
    pem+= 1*bpf # 1 float for temperature

    dual=N*K*bpf #i_n
    dual+=K*bpf #sum iPWL
    dual+=K*bpf # price to transformer
    dual+=K*N*bpf # price to EVs
    dual*=dualMaxIt

    admm=N*K*bpf #i_n
    admm+=K*S*bpf # iPWL
    admm+=K*bpf # price to transformer
    admm+=K*N*bpf # price to EVs
    admm+=K*S*bpf # aux variable to transformer
    admm+=K*N*bpf # aux variable to EVs
    admm*=admmMaxIt

    alad=N*K*bpf #i_n
    alad+=N*K*bpf #g_i
    alad+=N*K*bpf #g_s
    alad+=4*N*K*bpb # all Cs
    alad+=K*S*bpf # iPWL
    # alad+=S*K*bpf #g_z
    # alad+=K*bpf #g_t
    alad+=2*S*K*bpb+2*K*bpb # all Cs
    alad+=K*bpf # price to transformer
    alad+=K*N*bpf # price to EVs
    alad+=K*S*bpf+K*bpf # aux variable to transformer
    alad+=2*K*N*bpf # aux variables to EVs
    alad*=aladMaxIt

    pem/1e3 #kilobit
    dual/1e9 #gigabit
    admm/1e6 #megabit
    alad/1e6 #megabit
end

function othergraphs()
    c1=plot(hcat(dCM.couplConst,dCMadmm.couplConst,dCMalad.couplConst),labels=["Dual Ascent" "ADMM" "ALADIN"],
        xlabel="",ylabel="2-Norm Coupling Gap",yscale=:log10,legend=false)
    c2=plot(hcat(dCM.lamIt,dCMadmm.lamIt,dCMalad.lamIt),labels=["Dual Ascent" "ADMM" "ALADIN"],
        xlabel="Iterations",ylabel="2-Norm Lambda Gap",yscale=:log10)
    convPlot=plot(c1,c2,layout=(2,1))
    if lowRes
        pubPlot(convPlot,thickscale=0.4,sizeWH=(400,300),dpi=40)
    else
        pubPlot(convPlot,thickscale=0.8,sizeWH=(800,600),dpi=100)
    end

    if saveF savefig(convPlot,path*"dCMPlot.png") end
end
