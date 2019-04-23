
using Plots.PlotMeasures
#run comparison
#path = clips()
path=path*"ForVis\\"

path=path*"PWL\\"

cRun,runs, noLim=readRuns(path);
lowRes=true
saveF=false

#resPlot=compareRunsGraph(runs, cRun, noLim, saveResults,lowRes)
#cTable=compareRunsTable(runs,evS)

# for NL/PWL compare
cRun,runs, noLim=readRuns(path,true);
nl_pwlCompare(evS, runs,savefig,lowRes)

# @time runALADit(1)
#@time testALAD(1)

# Profile.clear()
# runALADit(1)
# @profile runALADit(1)
# Juno.profiler()


# Profile.clear()
# testDual(1)
# @profile testDual(1)
# Juno.profiler()


function readRuns(path,multiCentral=false)
    files = filter(x->occursin(".jld2",x), readdir(path))
    #evFile = filter(x->occursin("scenario",x), files)
    files = filter(x->occursin("_",x), files) # avoid evScenario

    # pemRun = filter(x->occursin("PEM",x), files)
    # files=setdiff(files,pemRun)
    noLimFile = filter(x->occursin("noLim",x), files)
    if multiCentral
        cFile=[]
        dFiles=files
    else
        cFile = filter(x->occursin("central",x), files)
        dFiles= setdiff(files,cFile)
    end

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

    # if length(evFile)>0
    #     evS=load(path*evFile[1])["evScenario"]
    # else
    #     evS=Nothing
    # end

    runs=Dict{String,Any}()
    for ii=1:length(dFiles)
        runs[dFiles[ii]]=load(path*dFiles[ii])
    end
    return cRun, runs, noLim#, evS
end

function compareRunsGraph(runs, cRun, noLim, saveF::Bool, lowRes::Bool)
    runNames=collect(keys(runs))
    cSol=cRun["solution"]
    numIt=size(runs[runNames[1]]["convMetrics"].coupl1Norm)[1]
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

    Plots.scalefontsizes(1.2)



    cAlpha=0.25
    lWidth=1.5
    #Time plots
    tempPlot=plot(1:Klen,cSol.Tactual,label="",seriescolor=:black,linealpha=cAlpha,
                    xlabel="",ylabel="Temp (C)",xticks=xticks,linewidth=8, title="a")
    plot!(tempPlot,1:Klen,evS.Tmax*ones(Klen),label="XFRM Limit",line=(:dash,:red),linewidth=2)
    plot!(tempPlot,T,labels="",seriescolor=plotColors,linewidth=lWidth)
    if noLim !=nothing
        plot!(tempPlot,1:Klen,noLim["solution"].Tactual,label="",seriescolor=allColors[P+1],linewidth=lWidth)
    end

    uSumPlot=plot(1:Klen,cSol.uSum,xlabel="",ylabel="Current Sum (kA)",xlims=(0,Klen),labels="Central",
                  seriescolor=:black,linewidth=6,linealpha=cAlpha,xticks=xticks)
    plot!(uSumPlot,uSum,labels=plotLabels,seriescolor=plotColors,legend=false,linewidth=lWidth)
    if noLim !=nothing
        plot!(uSumPlot,1:Klen,noLim["solution"].uSum,label="Uncoordinated",seriescolor=allColors[P+1],linewidth=lWidth)
    end

    iD=evS.iD_actual[1:Klen]
    loadPlot=plot(1:Klen,(cSol.uSum+iD)/Ntf*1000,xlabel="",ylabel="Current (A)",labels="",
                  seriescolor=:black,linewidth=6,linealpha=cAlpha,xticks=xticks, title="b")
    plot!(loadPlot,1:Klen,iD/Ntf*1000,label="Background Demand",line=(:dash),linewidth=2)
    #plot!(loadPlot,1:Klen,evS.ItotalMax/Ntf*1000*ones(Klen),label="Current Limit",line=(:dash,:red),linewidth=3)
    plot!(loadPlot,(uSum.+iD)/Ntf*1000,labels="",seriescolor=plotColors,linewidth=lWidth)
    if noLim !=nothing
        plot!(loadPlot,1:Klen,(noLim["solution"].uSum+iD)/Ntf*1000,label="",seriescolor=allColors[P+1],linewidth=lWidth)
    end

    target=zeros(Klen)
    for k in 1:Klen
        ind = evS.Kn.<=k
        target[k]=sum(evS.Snmin[ind])
    end

    snSumCentral=sum(cSol.Sn,dims=2)# sum up across N
    snSumPlot=plot(snSumCentral,xlabel="",ylabel="SoC Sum",xlims=(0,Klen),labels="",
                   seriescolor=:black,linewidth=4,linealpha=cAlpha,xticks=xticks)
    plot!(snSumPlot,target,label="SoC Target",line=(:dash))
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

    lamPlot=plot(1:Klen,cSol.lamCoupl/1000,xlabel="Time",ylabel=raw"Lambda ($/A)",labels="Central",
                   seriescolor=:black,linewidth=6,linealpha=cAlpha,xticks=xticks, title="c")
    plot!(lamPlot,Lam/1000,labels=plotLabels,seriescolor=plotColors,linewidth=lWidth)
    if noLim !=nothing
        plot!(lamPlot,1:Klen,noLim["solution"].lamCoupl/1000,label="Uncoordinated",seriescolor=allColors[P+1],linewidth=lWidth)
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
        pubPlot(resPlot,thickscale=1,sizeWH=(1000,600),dpi=40)
    end

    if saveF savefig(resPlot,path*"resPlot.png") end


    convItPlot=plot(convIt,xlabel="Time",ylabel="Number of Iterations",xlims=(0,Klen),xticks=xticks,labels=plotLabels,linewidth=2)
    savefig(convItPlot,path*"itPlot.png")

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
            convDict[convNames[ii]][1:cIt,i]=temp[1:cIt,ind]
            # fill in NaN for greater than convIt
            convDict[convNames[ii]][cIt+1:numIt,i].=NaN
        end
    end

    plotLabels=copy(permutedims(runNames))

    #
    # plotLabels=["NL" "PWL"]
    #plotLabels=["Equality ALADIN" "Inequlity ALADIN"]

    for i=1:length(plotLabels)
        plotLabels[i]=renameLabel(plotLabels[i])
    end


    Plots.scalefontsizes(1.2)

    maxIt=150
    #internal metrics
    couplLabel= L"||i_d[k+1]+\sum_{n=1}^N i_n[k+1]-\sum_{m=1}^M i_m^{PW}[k+1]||_1"
    couplLabel="1-Norm Coupling Gap"
    couplPlot=plot(convDict["coupl1Norm"],legend=false,xlims=(0,maxIt),yscale=:log10,ylabel=couplLabel,linewidth=2,title="a")
    dualLabel=L"||\lambda^{(p)}-\lambda^{(p-1)}||_2"
    dualLabel="2-Norm Lambda It Gap"
    dualPlot=plot(convDict["lamIt2Norm"],labels=plotLabels,xlims=(0,maxIt),yscale=:log10,xlabel="Iteration",
                    ylabel=dualLabel,linewidth=2,title="b")
    intConvPlot=plot(couplPlot,dualPlot,layout=(2,1))
    pubPlot(intConvPlot,thickscale=1,sizeWH=(1000,800),dpi=100)
    pname="convPlot.png"
    pname="nl_pwl_i$(ind).png"
    pname="eq_ineq_i$(ind).png"
    if saveF savefig(intConvPlot,path*pname) end

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


    savefPerc[savefPerc.==0] .=NaN
    convPlot= plot(savefPerc,xlabel="Iteration",ylabel="|o^p-o^*|/o^*",labels=["Dual" "ADMM" "ALADIN"],yscale=:log10, xlims=(0,200))
    savefig(convPlot,path*"savefPerc"*".png")

    savefItPerc[savefItPerc.==0] .=NaN
    convPlot= plot(savefItPerc,xlabel="Iteration",ylabel="|o^p-o^(p-1)|/o^(p-1)",labels=["Dual" "ADMM" "ALADIN"],yscale=:log10, xlims=(0,200))
    savefig(convPlot,path*"savefItPerc"*".png")


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
    Klen=evS.K

    compareSpeedTable = DataFrame(name=String[],timeTotal=Float64[],avgTimePerTs=Float64[],avgtimePerIt=Float64[],avgconvIt=Float64[],maxConvIt=Float64[])
    for keyI in keys(runs)
        println(keyI)
        loadF=runs[keyI]
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
        stats = [renameLabel(keyI) timeT avgTimeTs avgTime avgConv maxConv]
        push!(compareSpeedTable,stats)
    end
    #add central
    avgTimeTs=cRun["runTime"]/Klen
    stats = ["Central" cRun["runTime"] avgTimeTs avgTimeTs 1 1]
    push!(compareSpeedTable,stats)

    #EVC
    cUn=cRun["solution"].Un
    cLam = cRun["solution"].lamCoupl
    comparePerfTable = DataFrame(name=String[],curr2Norm=Float64[],currRMSE=Float64[],lam2Norm=Float64[],lamRMSE=Float64[])
    for keyI in keys(runs)
        println(keyI)
        loadF=runs[keyI]
        curr2Norm=norm(cUn*1000-loadF["solution"].Un*1000,2)
        # norm(cUn*1000,2)
        # norm(cUn-loadF["solution"].Un,2)
        # norm(cUn-loadF["solution"].Un,1)
        # norm((cUn-loadF["solution"].Un)./cUn,2)
        # norm((cUn-loadF["solution"].Un)./cUn,1)
        currRMSE=sqrt(sum((cUn*1000 .-loadF["solution"].Un*1000).^2)/length(cUn))

        lam2Norm=norm(cLam/1000-loadF["solution"].lamCoupl/1000,2)
        #norm((cLam/1000-loadF["solution"].lamCoupl/1000)./cLam/1000,2)
        lamRMSE=sqrt(sum((cLam/1000 .-loadF["solution"].lamCoupl/1000).^2)/length(cLam))

        stats = [renameLabel(keyI) curr2Norm currRMSE lam2Norm lamRMSE]
        push!(comparePerfTable,stats)
    end

    #Hub
    cU=cRun["solution"].U
    cLam = cRun["solution"].Lam
    comparePerfTable = DataFrame(name=String[],curr2Norm=Float64[],currRMSE=Float64[],lam2Norm=Float64[],lamRMSE=Float64[])
    for keyI in keys(runs)
        println(keyI)
        loadF=runs[keyI]
        curr2Norm=norm(cU*1000-loadF["solution"].U*1000,2)
        currRMSE=sqrt(sum((cU*1000 .-loadF["solution"].U*1000).^2)/length(cU))

        lam2Norm=norm(cLam/1000-loadF["solution"].Lam/1000,2)
        lamRMSE=sqrt(sum((cLam/1000 .-loadF["solution"].Lam/1000).^2)/length(cLam))

        stats = [renameLabel(keyI) curr2Norm currRMSE lam2Norm lamRMSE]
        push!(comparePerfTable,stats)
    end


    return compareSpeedTable, comparePerfTable
end

function compareHubsGraph(runs, cRun, noLim, saveF::Bool, lowRes::Bool)
    runNames=collect(keys(runs))
    cSol=cRun["solution"]
    Klen=size(cSol.Tactual)[1]
    P=length(runNames)
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
        uSum[:,i]=sum(runI["solution"].U,dims=2)
        # uSum[:,i]=runI["solution"].uSum[:,1]
        eSum[:,i]=[sum(runI["solution"].E[k,:]) for k in 1:Klen]# sum up across N
        eAvg[:,i]=[mean(runI["solution"].E[k,:]) for k in 1:Klen]# mean across N
    end

    plotLabels=permutedims(runNames)
    for i=1:length(plotLabels)
        plotLabels[i]=renameLabel(plotLabels[i])
    end
    allColors=get_color_palette(:auto, plot_color(:white), P+1)
    plotColors=allColors[1:P]'


    Plots.scalefontsizes(1.2)



    cAlpha=0.25
    lWidth=1.5
    #Time plots
    tempPlot=plot(1:Klen,cSol.Tactual,label="",seriescolor=:black,linewidth=8,linealpha=cAlpha,
                    xlabel="",ylabel="Temp (C)",xticks=xticks,title="a")
    plot!(tempPlot,1:Klen,hubS.Tmax*ones(Klen),label="XFRM Limit",line=(:dash,:red))
    plot!(tempPlot,T,labels="",seriescolor=plotColors,linewidth=lWidth)
    if noLim !=nothing
        plot!(tempPlot,1:Klen,noLim["solution"].Tactual,label="",seriescolor=allColors[P+1],linewidth=lWidth)
    end

    iD=hubS.iD_actual[1:Klen]
    loadPlot=plot(1:Klen,(cSol.uSum+iD)/Ntf,xlabel="",ylabel="Current (kA)",labels="",
                  seriescolor=:black,linewidth=8,linealpha=cAlpha,xticks=xticks,title="b")
    plot!(loadPlot,1:Klen,iD/Ntf,label="Background Demand",line=(:dash))
    plot!(loadPlot,(uSum.+iD)/Ntf,labels="",seriescolor=plotColors,linewidth=lWidth)
    if noLim !=nothing
        plot!(loadPlot,1:Klen,(noLim["solution"].uSum+iD)/Ntf,label="",seriescolor=allColors[P+1])
    end

    lamPlot=plot(1:Klen,cSol.Lam/1000,xlabel="Time",ylabel=raw"Lambda ($/A)",labels="Central",
                   seriescolor=:black,linewidth=8,linealpha=cAlpha,xticks=xticks,title="c")
    plot!(lamPlot,Lam/1000,labels=plotLabels,seriescolor=plotColors,linewidth=lWidth)
    if noLim !=nothing
        plot!(lamPlot,1:Klen,noLim["solution"].Lam,label="Uncoordinated",seriescolor=allColors[P+1])
    end

    resPlot=plot(tempPlot,loadPlot,lamPlot,layout=(3,1))
    if lowRes
        pubPlot(resPlot,thickscale=0.4,sizeWH=(400,300),dpi=40)
    else
        pubPlot(resPlot,thickscale=1,sizeWH=(1000,600),dpi=40)
    end

    if saveF savefig(resPlot,path*"hubPlot.png") end


    convItPlot=plot(convIt,xlabel="Time",ylabel="Number of Iterations",xticks=xticks,labels=plotLabels,linewidth=2,yscale=:log10)
    if lowRes
        pubPlot(convItPlot,thickscale=0.4,sizeWH=(400,300),dpi=40)
    else
        pubPlot(convItPlot,thickscale=1.5,sizeWH=(1000,600),dpi=40)
    end

    if saveF savefig(convItPlot,path*"hubConvPlot.png") end

    return resPlot, convPlot
end

function calcPrivacy(N,K1,S)
    aladMaxIt=10
    admmMaxIt=100
    dualMaxIt=1000

    bpf= 64 #bits per double
    bpi= 16 #bits per int
    bpb= 1#bits per bool

    pemIt=N*bpi # 1 Req Int per EV
    pemIt+=N*bpi # 1 Rec Int per EV
    pemEV=pemIt/N  # data per EV
    pemIt+= 1*bpf # 1 float for temperature
    pem=pemIt

    dualIt=N*K1*bpf #i_n
    dualIt+=K1*N*bpf # price to EVs
    dualEV=dualIt/N*dualMaxIt # data per EV
    dualIt+=K1*bpf #sum iPWL
    dualIt+=K1*bpf # price to transformer
    dual=dualIt*dualMaxIt

    admmIt=N*K1*bpf #i_n
    admmIt+=K1*N*bpf # price to EVs
    admmIt+=K1*N*bpf # aux variable to EVs
    admmEV=admmIt/N*admmMaxIt # data per EV
    admmIt+=K1*S*bpf # iPWL
    admmIt+=K1*bpf # price to transformer
    admmIt+=K1*S*bpf # aux variable to transformer
    admm=admmIt*admmMaxIt

    aladIt=N*K1*bpf #i_n
    aladIt+=N*K1*bpf #g_i
    aladIt+=N*K1*bpf #g_s
    aladIt+=4*N*K1*bpb # all Cs
    aladIt+=K1*N*bpf # price to EVs
    aladIt+=2*K1*N*bpf # aux variables to EVs
    aladEV=aladIt/N*aladMaxIt # data per EV
    aladIt+=K1*S*bpf # iPWL
    # alad+=S*K*bpf #g_z
    # alad+=K*bpf #g_t
    aladIt+=2*S*K1*bpb+2*K1*bpb # all Cs
    aladIt+=K1*bpf # price to transformer
    aladIt+=K1*S*bpf+K1*bpf # aux variable to transformer
    alad=aladIt*aladMaxIt


    # per iteration total data
    aladIt/1e6 #megabit
    dualIt/1e9 #gigabit
    admmIt/1e6 #megabit
    pemIt/1e3 #kilobit


    # total time step per EV
    aladEV/1e3 #kilobit
    dualEV/1e6 #megabit
    admmEV/1e6 #megabit
    pemEV #bit

    # percent of average US internet
    bandwidth=18.7e6*180
    aladEV/bandwidth*100
    admmEV/bandwidth*100
    dualEV/bandwidth*100
    pemEV/bandwidth*100

    # total time step total data
    alad/1e6 #megabit
    dual/1e9 #gigabit
    admm/1e6 #megabit
    pem/1e3 #kilobit
end

function pwlSegPlot()
    #plot central PWL solutions

    I=evS.ItotalMax
    S=evS.S
    deltaI=I/S

    Z=cSol.Z
    Tactual=cSol.Tactual
    uSum=cSol.uSum
    K=evS.K
    Tpred=cSol.Tpred

    # Z=zeros(horzLen+1,S)
    # for ii= 1:S
    # 	Z[:,ii]=zRaw[collect(ii:S:length(zRaw))]
    # end
    # Tactual=Tactual
    # uSum=uSum
    # K=horzLen+1
    # Tpred=tRaw

    solX=sum(Z,dims=2)
    solY=zeros(K,1)
    for i=1:S
        global solY[:,1]=solY[:,1].+deltaI*(2*i-1)*Z[:,i]
    end
    tightT=abs.(Tpred.-evS.Tmax).<1e-2
    indT=findlast(tightT.==true)[1]

    dX=100
    x=range(0+I/(S*dX),I,length=S*dX)
    #segP=plot(x,x.^2,label=L"i^2")
    segs=zeros(length(x))
    prevY=0
    for i=1:S
        tempY=(i-1)*dX+1:i*dX
        alpha=(2*i-1)*deltaI
        segs[tempY]=x[1:dX]*alpha.+prevY
        global prevY=segs[i*dX]
    end

    Plots.scalefontsizes(1.2)
    mSize1=6
    mSize2=4

    segP=plot(x,x.^2,label=L"i^2",legend=:topleft)
    plot!(segP,x,segs,label="PWL Approx",xlabel="Current (kA)")
    scatter!(segP,solX[1:indT],solY[1:indT],label="Before Overload",markersize=mSize1,markerstrokewidth=0,
                seriescolor="red",xlims=(minimum(solX)-deltaI,maximum(solX)+deltaI))
    scatter!(segP,solX[indT+1:K],solY[indT+1:K],label="After Overload",markersize=mSize2,markerstrokewidth=0,seriescolor="green")
    ylabel!("e[k]")
    xlims!(segP,(0,25))


    segP2=plot(x,x.^2,label=L"i^2",linewidth=2,widen=true,legend=:topleft)
    plot!(segP2,x,segs,label="PWL Approx",xlabel="Current (kA)",linewidth=2)
    scatter!(segP2,solX[1:indT],solY[1:indT],label="Before Overload",markersize=mSize1,markerstrokewidth=0,
                seriescolor="red",xlims=(minimum(solX)-deltaI,maximum(solX)+deltaI))
    ylims!(segP2,(380,500))
    xlims!(segP2,(19.7,22))


    segPlots=plot(segP,segP2,layout=(1,2),size=(1200,800))
    savefig(segPlots,path*"pwlSegPlot.png")



    #check one spot
    i=3
    pwlI=sum(Z[i,:])
    nlI=uSum[i]+evS.iD_actual[i]
    nlI2=(nlI)^2
    pwlI2=solY[i]
    full_segs = Int(floor(pwlI/deltaI))
    prevY=(full_segs*deltaI)^2
    pwlSeg=prevY+deltaI*(2*(full_segs+1)-1)*(pwlI-full_segs*deltaI)

    nlI2<pwlSeg<pwlI2
end

function nl_pwlCompare(evS, runs, saveF::Bool, lowRes::Bool)
    runNames=collect(keys(runs))
    Klen=evS.K
    P=length(runNames)
    N=evS.N

    Lam = zeros(Klen,P)
    T = zeros(Klen,P)
    Sn = zeros(Klen*N,P)
    uSum = zeros(Klen,P)
    snSum = zeros(Klen,P)
    snAvg = zeros(Klen,P)

    for i in 1:length(runNames)
        print(runNames[i])
        runI=runs[runNames[i]]
        print(": ")
        println(runI["runTime"]/60)
    end

    for i in 1:length(runNames)
        println(runNames[i])
        runI=runs[runNames[i]]
        Lam[:,i]=runI["solution"].lamCoupl
        T[:,i]=runI["solution"].Tactual[:,1]
        uSum[:,i]=runI["solution"].uSum[:,1]
        Sn[:,i]=runI["solution"].Sn'[:]
        snSum[:,i]=[sum(runI["solution"].Sn[k,:]) for k in 1:Klen]# sum up across N
        snAvg[:,i]=[mean(runI["solution"].Sn[k,:]) for k in 1:Klen]# mean across N
    end

    plotLabels=copy(permutedims(runNames))
    plotLabels=["Uncoordinated" "Centrally Coordinated"]
    plotLabels=["NL" "PWL"]

    allColors=get_color_palette(:auto, plot_color(:white), P+1)
    plotColors=allColors[1:P]'

    Plots.scalefontsizes(1.2)


    #Time plots
    tempPlot=plot(T,label="",seriescolor=plotColors,
                    xlabel="",ylabel="Temp (C)",xticks=xticks)
    plot!(tempPlot,evS.Tmax*ones(Klen),label="XFRM Limit",line=(:dash,:red))

    iD=evS.iD_actual
    loadPlot=plot(1:Klen,(uSum.+iD)/Ntf*1000,xlabel="",ylabel="Total Load (A)",labels="",
                  seriescolor=plotColors,xticks=xticks)
    plot!(loadPlot,1:Klen,(iD)/Ntf*1000,label="Background Demand",line=(:dash))


    lamPlot=plot(Lam/1000,xlabel="Time",ylabel=raw"Lambda ($/A)",labels=plotLabels,
                   seriescolor=plotColors,xticks=xticks)

    #resPlot=plot(tempPlot,loadPlot,snSumPlot,lamPlot,layout=(4,1))
    resPlot=plot(tempPlot,loadPlot,lamPlot,layout=(3,1))
    if lowRes
        pubPlot(resPlot,thickscale=0.4,sizeWH=(400,300),dpi=40)
    else
        pubPlot(resPlot,thickscale=1,sizeWH=(1000,600),dpi=40)
    end

    pName="nl_pwl.png"
    pName="central_uncoord.png"
    if saveF savefig(resPlot,path*pName) end

    return resPlot
end

function forecastError()

    #plotLabels
    plotLabels=["Perfect Forecast" "Forecast Error"]
    iD=uSum
    iD[:,1]=(iD[:,1] .+ evS.iD_pred)/Ntf*1000
    iD[:,2]=(iD[:,2] .+ iD_actual)/Ntf*1000

    loadPlot=plot(1:Klen,iD,xlabel="",ylabel="Total Load (A)",labels="",
                  seriescolor=plotColors,xticks=xticks)
    plot!(loadPlot,1:Klen,hcat(evS.iD_pred/Ntf*1000,iD_actual/Ntf*1000),label=["Predicted Background Demand" "Actual Background Demand"],line=(:dash))
    savefig(resPlot,path*"iDError.png")
end

function scenarioPlots(evS)

    #back calculate
    # b_options=[40 60 75 100]
    # ratio = evS.ηP/(evS.Ts*240*1000)*3.6e6
    # b_close=round.(0.85 ./ ratio)
    # ind=1
    # b_kWh=zeros(Int,N)
    # for i in 1:evS.N
    #     #println(i)
    #     b_kWh[i]=Int(b_options[abs.(b_close[i].-b_options).<=7][1])
    # end
    b_kWh=evS.b_kWh

    # b_hist=histogram(b_kWh,nbins=20,legend=false,xlabel="EV Battery Size (kWh)",ylabel= "Number of EVs")
    # imax_hist=histogram(evS.imax,nbins=12,legend=false,xlabel="EV Max Charging Power (kW)",ylabel= "")
    # battParamsPlot=plot(b_hist,imax_hist,layout=(1,2))
    tickS=42.5:10:102.5
    xticksHist=(tickS,string.(collect(40:10:100)))
    battParamsPlot=histogram2d(b_kWh,evS.imax*1000,nbins=20,xlabel="EV Battery Size (kWh)",ylabel="EV Max Charging Power (A)",
        colorbar_title="Number of EVs",thickness_scaling=1.8,dpi=100,size=(1000,500),bar_edges=false,xticks=xticksHist)
    savefig(battParamsPlot,path*"Scenario\\battParamsPlot.png")

    s0Plot=histogram(evS.s0,nbins=40,legend=false,xlabel="EV Initial SoC (%)",ylabel= "Number of EVs",title="(a)")

    departTime=evS.Kn*evS.Ts/(60*60).+(20-24)
    stT=Time(6,30)
    endT=Time(10,0)
    XlabelsAM=vcat(collect(stT:Dates.Second(round(evS.Ts)):endT))
    xticksAM=(6.5:0.5:10,Dates.format.(XlabelsAM[1:10:71],"HH:MM"))


    minSoCPlot=histogram2d(departTime,evS.Snmin,nbins=12,xlabel="Time",ylabel="Minimum Departure SoC",
        colorbar_title="Number of EVs",title="(b)",xticks=xticksAM)
    evPlot=plot(s0Plot,minSoCPlot,layout=(1,2),thickness_scaling=1.5,dpi=100,size=(1000,400))
    savefig(evPlot,path*"Scenario\\evParamsPlot.png")

    num_homes =1000
    idPlot=plot(evS.iD_actual*240/num_homes,xticks=xticks,widen=false,legend=false,ylabel="Demand (kW)",title="(a)")
    #mean(iD_pred*240/num_homes)
    #maximum(iD_pred*240/num_homes)
    tambPlot=plot(evS.Tamb_raw,xticks=xticks,legend=false,xlabel="Time",ylabel="Temperature (C)",widen=false,title="(b)")
    backgroundPlot=plot(idPlot,tambPlot,layout=(2,1))
    pubPlot(backgroundPlot,thickscale=1.4,sizeWH=(800,400),dpi=100)
    savefig(backgroundPlot,path*"Scenario\\backgroundParamsPlot.png")
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


    #comparing dual variables
    dualComp=plot(hcat(cSol.lamCoupl/1000,cSol.lamTemp/1000),labels=["Coupling Constraint Dual" "Temperature Limit Dual"],
                    xlabel="Time",xticks=xticks,xlims=(0,evS.K))
    pubPlot(dualComp,thickscale=1.5,sizeWH=(800,400),dpi=100)
    savefig(dualComp,path*"dualCompPlot.png")
end
