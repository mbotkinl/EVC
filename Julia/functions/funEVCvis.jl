
using Plots.PlotMeasures
#run comparison
#path = clips()
path=path*"ForVis\\"
cRun,runs, noLim=readRuns(path);
lowRes=true
savefig=false

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

    # Plots.scalefontsizes(1.2)

    #Time plots
    tempPlot=plot(1:Klen,cSol.Tactual,label="",seriescolor=:black,linealpha=0.4,xlims=(0,Klen),
                    xlabel="",ylabel="Temp (C)",xticks=xticks,linewidth=8)
    plot!(tempPlot,1:Klen,evS.Tmax*ones(Klen),label="XFRM Limit",line=(:dash,:red),linewidth=4)
    plot!(tempPlot,T,labels="",seriescolor=plotColors,linewidth=2)
    if noLim !=nothing
        plot!(tempPlot,1:Klen,noLim["solution"].Tactual,label="",seriescolor=allColors[P+1],linewidth=2)
    end

    uSumPlot=plot(1:Klen,cSol.uSum,xlabel="",ylabel="Current Sum (kA)",xlims=(0,Klen),labels="Central",
                  seriescolor=:black,linewidth=6,linealpha=0.4,xticks=xticks)
    plot!(uSumPlot,uSum,labels=plotLabels,seriescolor=plotColors,legend=false,linewidth=2)
    if noLim !=nothing
        plot!(uSumPlot,1:Klen,noLim["solution"].uSum,label="Uncoordinated",seriescolor=allColors[P+1],linewidth=2)
    end

    iD=evS.iD_actual[1:Klen]
    loadPlot=plot(1:Klen,cSol.uSum+iD,xlabel="",ylabel="Total Load (kA)",xlims=(0,Klen),labels="",
                  seriescolor=:black,linewidth=6,linealpha=0.4,xticks=xticks)
    plot!(loadPlot,1:Klen,iD,label="Background Demand",line=(:dash),linewidth=3)
    plot!(loadPlot,1:Klen,evS.ItotalMax*ones(Klen),label="Current Limit",line=(:dash,:red),linewidth=3)
    plot!(loadPlot,uSum.+iD,labels="",seriescolor=plotColors,linewidth=2)
    if noLim !=nothing
        plot!(loadPlot,1:Klen,noLim["solution"].uSum+iD,label="",seriescolor=allColors[P+1],linewidth=2)
    end

    target=zeros(Klen)
    for k in 1:Klen
        ind = evS.Kn.<=k
        target[k]=sum(evS.Snmin[ind])
    end

    snSumCentral=sum(cSol.Sn,dims=2)# sum up across N
    snSumPlot=plot(snSumCentral,xlabel="",ylabel="SoC Sum",xlims=(0,Klen),labels="",
                   seriescolor=:black,linewidth=4,linealpha=0.25,xticks=xticks)
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

    lamPlot=plot(1:Klen,cSol.lamCoupl/1000,xlabel="Time",ylabel=raw"Lambda ($/A)",xlims=(0,Klen),labels="Central",
                   seriescolor=:black,linewidth=6,linealpha=0.4,xticks=xticks)
    plot!(lamPlot,Lam/1000,labels=plotLabels,seriescolor=plotColors,linewidth=2)
    if noLim !=nothing
        plot!(lamPlot,1:Klen,noLim["solution"].lamCoupl/1000,label="Uncoordinated",seriescolor=allColors[P+1],linewidth=2)
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
        pubPlot(resPlot,thickscale=1,sizeWH=(1000,600),dpi=100)
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
    # plotLabels=["Equality ALADIN" "Inequlity ALADIN"]

    for i=1:length(plotLabels)
        plotLabels[i]=renameLabel(plotLabels[i])
    end


    #Plots.scalefontsizes(1.2)

    maxIt=150
    #internal metrics
    couplLabel= L"||i_d[k+1]+\sum_{n=1}^N i_n[k+1]-\sum_{m=1}^M i_m^{PW}[k+1]||_1"
    couplLabel="1-Norm Coupling Gap"
    couplPlot=plot(convDict["coupl1Norm"],legend=false,xlims=(0,maxIt),yscale=:log10,ylabel=couplLabel,linewidth=2)
    dualLabel=L"||\lambda^{(p)}-\lambda^{(p-1)}||_2"
    dualLabel="2-Norm Lambda It Gap"
    dualPlot=plot(convDict["lamIt2Norm"],labels=plotLabels,xlims=(0,maxIt),yscale=:log10,xlabel="Iteration",ylabel=dualLabel,linewidth=2)
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
        stats = [renameLabel(keyI) timeT avgTimeTs avgTime avgConv maxConv]
        push!(compareTable,stats)
    end

    #add central
    avgTimeTs=cRun["runTime"]/Klen
    stats = ["Central" cRun["runTime"] avgTimeTs avgTimeTs 1 1]
    push!(compareTable,stats)

    return compareTable
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

    #Time plots
    tempPlot=plot(1:Klen,cSol.Tactual,label="",seriescolor=:black,linewidth=6,linealpha=0.25,xlims=(0,Klen),
                    xlabel="",ylabel="Temp (C)",xticks=xticks)
    plot!(tempPlot,1:Klen,hubS.Tmax*ones(Klen),label="XFRM Limit",line=(:dash,:red))
    plot!(tempPlot,T,labels="",seriescolor=plotColors,linewidth=2)
    if noLim !=nothing
        plot!(tempPlot,1:Klen,noLim["solution"].Tactual,label="",seriescolor=allColors[P+1])
    end

    iD=hubS.iD_actual[1:Klen]
    loadPlot=plot(1:Klen,cSol.uSum+iD,xlabel="",ylabel="Total Load (kA)",xlims=(0,Klen),labels="",
                  seriescolor=:black,linewidth=6,linealpha=0.25,xticks=xticks)
    plot!(loadPlot,1:Klen,iD,label="Background Demand",line=(:dash))
    plot!(loadPlot,uSum.+iD,labels="",seriescolor=plotColors,linewidth=2)
    if noLim !=nothing
        plot!(loadPlot,1:Klen,noLim["solution"].uSum+iD,label="",seriescolor=allColors[P+1])
    end

    lamPlot=plot(1:Klen,cSol.Lam/1000,xlabel="Time",ylabel=raw"Lambda ($/A)",xlims=(0,Klen),labels="Central",
                   seriescolor=:black,linewidth=6,linealpha=0.25,xticks=xticks)
    plot!(lamPlot,Lam/1000,labels=plotLabels,seriescolor=plotColors,linewidth=2)
    if noLim !=nothing
        plot!(lamPlot,1:Klen,noLim["solution"].Lam,label="Uncoordinated",seriescolor=allColors[P+1])
    end

    resPlot=plot(tempPlot,loadPlot,lamPlot,layout=(3,1))
    if lowRes
        pubPlot(resPlot,thickscale=0.4,sizeWH=(400,300),dpi=40)
    else
        pubPlot(resPlot,thickscale=1,sizeWH=(1000,600),dpi=100)
    end

    if saveF savefig(resPlot,path*"hubPlot.png") end


    convItPlot=plot(convIt,xlabel="Time",ylabel="Number of Iterations",xlims=(0,Klen),xticks=xticks,labels=plotLabels,linewidth=2)
    if lowRes
        pubPlot(convItPlot,thickscale=0.4,sizeWH=(400,300),dpi=40)
    else
        pubPlot(convItPlot,thickscale=1.5,sizeWH=(1000,600),dpi=100)
    end

    if saveF savefig(convItPlot,path*"hubConvPlot.png") end

    return resPlot, convPlot
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


    segP=plot(x,x.^2,label=L"i^2",legend=:topleft)
    plot!(segP,x,segs,label="PWL Approx",xlabel="Current (kA)")
    scatter!(segP,solX[1:indT],solY[1:indT],label="Before Overload",markersize=10,markerstrokewidth=0,
                seriescolor="red",xlims=(minimum(solX)-deltaI,maximum(solX)+deltaI))
    scatter!(segP,solX[indT+1:K],solY[indT+1:K],label="After Overload",markersize=5,markerstrokewidth=0,seriescolor="green")
    ylabel!("e[k]")
    xlims!(segP,(0,25))


    segP2=plot(x,x.^2,label=L"i^2",linewidth=2,widen=true,legend=:topleft)
    plot!(segP2,x,segs,label="PWL Approx",xlabel="Current (kA)",linewidth=2)
    scatter!(segP2,solX[1:indT],solY[1:indT],label="Before Overload",markersize=15,markerstrokewidth=0,
                seriescolor="red",xlims=(minimum(solX)-deltaI,maximum(solX)+deltaI))
    ylims!(segP2,(380,500))
    xlims!(segP2,(19.5,22))


    segPlots=plot(segP,segP2,layout=(1,2))
    pubPlot(segPlots,thickscale=1.5,sizeWH=(1200,400),dpi=100)
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
    tempPlot=plot(1:Klen,T,label="",seriescolor=plotColors,xlims=(0,Klen),
                    xlabel="",ylabel="Temp (K)",xticks=xticks)
    plot!(tempPlot,1:Klen,evS.Tmax*ones(Klen),label="XFRM Limit",line=(:dash,:red))

    iD=evS.iD_actual
    loadPlot=plot(1:Klen,uSum.+iD,xlabel="",ylabel="Total Load (kA)",xlims=(0,Klen),labels="",
                  seriescolor=plotColors,xticks=xticks)
    plot!(loadPlot,1:Klen,iD,label="Background Demand",line=(:dash))


    lamPlot=plot(1:Klen,Lam/1000,xlabel="Time",ylabel=raw"Lambda ($/A)",xlims=(0,Klen),labels=plotLabels,
                   seriescolor=plotColors,xticks=xticks)

    #resPlot=plot(tempPlot,loadPlot,snSumPlot,lamPlot,layout=(4,1))
    resPlot=plot(tempPlot,loadPlot,lamPlot,layout=(3,1))
    if lowRes
        pubPlot(resPlot,thickscale=0.4,sizeWH=(400,300),dpi=40)
    else
        pubPlot(resPlot,thickscale=1.4,sizeWH=(1000,600),dpi=100)
    end

    if saveF savefig(resPlot,path*"nl_pwl.png") end

    return resPlot
end

function forecastError()

    #plotLabels
    plotLabels=["Perfect Forecast" "Forecast Error"]
    iD=uSum
    iD[:,1]=iD[:,1] .+ evS.iD_pred
    iD[:,2]=iD[:,2] .+ iD_actual

    loadPlot=plot(1:Klen,iD,xlabel="",ylabel="Total Load (kA)",xlims=(0,Klen),labels="",
                  seriescolor=plotColors,xticks=xticks)
    plot!(loadPlot,1:Klen,hcat(evS.iD_pred,iD_actual),label=["Predicted Background Demand" "Actual Background Demand"],line=(:dash))
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
    battParamsPlot=histogram2d(b_kWh,evS.imax*1000,nbins=20,xlabel="EV Battery Size (kWh)",ylabel="EV Max Charging Power (A)",
        colorbar_title="Number of EVs",thickness_scaling=1.8,dpi=100,size=(1000,500),bar_edges=false)
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
