# helper functions for EVC code
using JLD2
using DataFrames
using FileIO
# using Gadfly
# using Cairo #for png output
# using Fontconfig

function saveRun(path::String, filename::String, time::Float64, scenario::scenarioStruct, solution, convMetrics=convMetricsStruct(), convIt=1)
    save(path*filename*".jld2","runTime", time, "scenario", scenario, "solution", solution,
    "convMetrics", convMetrics, "convIt", convIt)
end

function clips()
    t=clipboard()
    return t[2:length(t)-1]
end

function clipr()
    t=clipboard()
    return JLD.load(t[2:length(t)-1])
end

function readRuns(path)
    files = filter(x->contains(x,".jld"), readdir(path))
    files = filter(x->contains(x,"_"), files) # avoid evScenario
    cFile = filter(x->contains(x,"central"), files)
    cRun=JLD.load(path*cFile[1])

    dFiles= setdiff(files,cFile)
    runs=Dict()
    for file in dFiles
        runs[file]=JLD.load(path*file)
    end
end

# until converted to Gadfly***
# function compareRunsGraph(runs, cRun)
#     names=collect(keys(runs))
#     cSol=cRun["solution"]
#     numIt=size(runs[names[1]]["convMetrics"].lam)[1]
#     Klen=size(cSol.Xt)[1]
#     P=length(names)
#     evS=cRun["scenario"]
#     N=evS.N
#
#     objPerc = zeros(numIt+1,P)
#     objConv = zeros(numIt,P)
#     lamConv = zeros(numIt,P)
#     lamRMSE = zeros(numIt,P)
#     lamInfNorm = zeros(numIt,P)
#     T = zeros(Klen,P)
#     Sn = zeros(Klen*N,P)
#     uSum = zeros(Klen,P)
#     snSum = zeros(Klen,P)
#
#     for i in 1:length(names)
#         println(names[i])
#         run=runs[names[i]]
#         cIt=run["convIt"]
#         objPerc[:,i]=abs.(cSol.objVal-run["solution"].objVal')/cSol.objVal*100
#         objConv[:,i]=run["convMetrics"].obj
#         lamConv[:,i]=run["convMetrics"].lam
#         if typeof(run["solution"]) == centralSolutionStruct #PEM
#             lamRMSE[:,i]=zeros(numIt)
#             lamInfNorm[:,i]=zeros(numIt)
#         else
#             lamRMSE[:,i]=[sqrt(1/Klen*sum((run["solution"].Lam[k,it]-cSol.lamCoupl[k])^2/abs(cSol.lamCoupl[k]) for k=1:Klen)) for it=1:numIt]
#             lamInfNorm[:,i]=[maximum(abs.(run["solution"].Lam[:,it]-cSol.lamCoupl)) for it=1:numIt]
#         end
#         T[:,i]=run["solution"].Xt[:,cIt]
#         uSum[:,i]=run["solution"].uSum[:,cIt]
#         if size(run["solution"].Sn)[1]>Klen+1
#             Sn[:,i]=run["solution"].Sn[:,cIt]
#             snSum[:,i]=[sum(Sn[N*(k-1)+n,i] for n=1:N) for k in 1:Klen]# sum up across N
#         else
#             Sn[:,i]=run["solution"].Sn'[:]
#             snSum[:,i]=[sum(run["solution"].Sn[k,:]) for k in 1:Klen]# sum up across N
#         end
#     end
#
#     # ty=[match(r"central|d",names[i]).match for i in 1:length(names)]
#     # form = [match(r"NL|N|relax",names[i]).match for i in 1:length(names)]
#     #lamConv=DataFrame(hcat(names,transpose(lamConv)))
#     lamConv=DataFrame(lamConv)
#     names!(lamConv.colindex, map(Symbol, names))
#     lamRMSE=DataFrame(lamRMSE)
#     names!(lamRMSE.colindex, map(Symbol, names))
#     lamInfNorm=DataFrame(lamInfNorm)
#     names!(lamInfNorm.colindex, map(Symbol, names))
#     objConv=DataFrame(objConv)
#     names!(objConv.colindex, map(Symbol, names))
#     objPerc=DataFrame(objPerc)
#     names!(objPerc.colindex, map(Symbol, names))
#     T=DataFrame(T)
#     names!(T.colindex, map(Symbol, names))
#     uSum=DataFrame(uSum)
#     names!(uSum.colindex, map(Symbol, names))
#     snSum=DataFrame(snSum)
#     names!(snSum.colindex, map(Symbol, names))
#
#     lamPlot=plot(lamConv,x=Row.index,y=Col.value,Geom.point,Geom.line, color=Col.index,Scale.y_log10,
#     			Guide.xlabel("Iteration"), Guide.ylabel("2-norm Lambda gap",orientation=:vertical),
#     			Theme(background_color=colorant"white",major_label_font_size=15pt,line_width=2pt,
#     			minor_label_font_size=20pt,key_label_font_size=20pt,key_position=:top))
#     lamRMSEPlot=plot(lamRMSE,x=Row.index,y=Col.value,Geom.point,Geom.line, color=Col.index,Scale.y_log10,
#     			Guide.xlabel("Iteration"), Guide.ylabel("Relative RMSE Lambda Gap",orientation=:vertical),
#     			Theme(background_color=colorant"white",major_label_font_size=15pt,line_width=2pt,
#     			minor_label_font_size=20pt,key_label_font_size=20pt,key_position=:top))
#     lamInfNormPlot=plot(lamInfNorm,x=Row.index,y=Col.value,Geom.point,Geom.line, color=Col.index,Scale.y_log10,
#     			Guide.xlabel("Iteration"), Guide.ylabel("Max Lambda Gap",orientation=:vertical),
#     			Theme(background_color=colorant"white",major_label_font_size=15pt,line_width=2pt,
#     			minor_label_font_size=20pt,key_label_font_size=20pt,key_position=:top))
#
#     objNormPlot=plot(objConv,x=Row.index,y=Col.value,Geom.point,Geom.line, color=Col.index,Scale.y_log10,
#     			Guide.xlabel("Iteration"), Guide.ylabel("Objective Value Magintude Gap",orientation=:vertical),
#     			Theme(background_color=colorant"white",major_label_font_size=15pt,line_width=2pt,
#     			minor_label_font_size=20pt,key_label_font_size=12pt,key_position=:none))
#     objPercPlot=plot(objPerc,x=Row.index,y=Col.value,Geom.point,Geom.line, color=Col.index,Scale.y_log10,Coord.Cartesian(xmin=0,xmax=numIt+1),
#                 Guide.xlabel("Iteration"), Guide.ylabel("Objective Value Percentage Gap",orientation=:vertical),
#                 Theme(background_color=colorant"white",major_label_font_size=15pt,line_width=2pt,
#                 minor_label_font_size=20pt,key_label_font_size=12pt,key_position=:none))
#
#     tempPlot=plot(layer(T,x=Row.index,y=Col.value,Geom.point,Geom.line,color=Col.index),
#     		layer(yintercept=[evS.Tmax],Geom.hline(color=["red"],style=:dot),Theme(line_width=3pt)),
#     		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)",orientation=:vertical),
#     		Coord.Cartesian(xmin=0,xmax=Klen),Theme(background_color=colorant"white",major_label_font_size=24pt,
#     		minor_label_font_size=20pt,key_label_font_size=12pt,key_position=:none))
#
#     uSumPlot=plot(layer(uSum,x=Row.index,y=Col.value,Geom.point,Geom.line,color=Col.index),
#     		# layer(yintercept=[evS.Tmax],Geom.hline(color=["red"],style=:dot),Theme(line_width=3pt)),
#     		Guide.xlabel("Time"), Guide.ylabel("Current Sum (kA)",orientation=:vertical),
#             Guide.colorkey(title="",labels=names),
#     		Coord.Cartesian(xmin=0,xmax=Klen),Theme(background_color=colorant"white",major_label_font_size=24pt,
#     		minor_label_font_size=20pt,key_label_font_size=12pt,key_position=:none))
#
#     target=zeros(Klen)
#     for k in 1:Klen
#         ind = evS.Kn.<=k
#         target[k]=sum(evS.Snmin[ind])
#     end
#
#     snSumPlot=plot(layer(snSum,x=Row.index,y=Col.value,Geom.point,Geom.line,color=Col.index),
#     		layer(x=1:Klen,y=target,Geom.line,Theme(line_width=3pt,default_color=colorant"red",line_style=:dot)),
#     		Guide.xlabel("Time"), Guide.ylabel("SOC Sum",orientation=:vertical),
#             Guide.colorkey(title="",labels=names),
#     		Coord.Cartesian(xmin=0,xmax=Klen),Theme(background_color=colorant"white",major_label_font_size=24pt,
#     		minor_label_font_size=20pt,key_label_font_size=12pt,key_position=:none))
#
#     # R plot***
#     Rmax=zeros(Klen,P)
#     Ravg=zeros(Klen,P)
#     for p=1:P
#         for k=1:Klen
#             R=[max((evS.Snmin[n,1]-Sn[N*(k-1)+n,p]),0)./(evS.Kn[n,1]-k) for n=1:N]
#             Rmax[k,p]=maximum(R)
#             Ravg[k,p]=mean(R)
#         end
#     end
#
#     Rmax=DataFrame(Rmax)
#     Ravg=DataFrame(Ravg)
#     Rplot=plot(layer(Ravg,x=Row.index,y=Col.value,Geom.point,Geom.line,color=Col.index),
#             Guide.xlabel("Time"), Guide.ylabel("R Max",orientation=:vertical),
#             Guide.colorkey(title="",labels=names),
#             Coord.Cartesian(xmin=0,xmax=Klen),Theme(background_color=colorant"white",major_label_font_size=24pt,
#             minor_label_font_size=20pt,key_label_font_size=12pt,key_position=:none))
#
#     p=vstack(lamRMSEPlot, objPercPlot, tempPlot, uSumPlot,snSumPlot)
#     fName="compPlot.png"
#     draw(PNG(path*fName, 28inch, 20inch),p)
# end

function compareRunsTable(runs)
    # compareTable = DataFrame(name=String[],time=Float64[],cLamDiff=Float64[],lamDiff=Float64[],
    # cObjDiff=Float64[],objDiff=Float64[])
    compareTable = DataFrame(name=String[],timeTotal=Float64[],timePerIt=Float64[],convIt=Int64[])
    for key in keys(runs)
        println(key)
        loadF=runs[key]
        timeT=loadF["time"]
        #cm=loadF["convMetrics"]
        convIt=loadF["convIt"]
        #ind= if convIt>0 convIt-1 else length(cm.lam) end
        #stats = [key timeT minimum(cm.lam[1:ind]) cm.lam[ind-1]-cm.lam[ind-2] minimum(cm.obj[1:ind]) cm.obj[ind-1]-cm.obj[ind-2]]
        stats = [key timeT timeT/convIt convIt]
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
