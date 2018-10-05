# helper functions for EVC code
using JLD
using DataFrames
using Gadfly
using Cairo #for png output
using Fontconfig

function saveRun(path::String, filename::String, time::Float64, scenario::scenarioStruct, solution, convMetrics=convMetricsStruct(), convIt=1)
    JLD.save(path*filename*".jld","time", time, "scenario", scenario, "solution", solution,
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
    files = filter(x->contains(x,"_"), readdir(path))
    runs=Dict()
    for file in files
        runs[file]=JLD.load(path*file)
    end
end


function compareRunsGraph(runs)
    names=collect(keys(runs))
    numIt=size(runs[names[1]]["convMetrics"].lam)[1]
    Klen=size(runs[names[1]]["solution"].Xt)[1]
    P=length(names)
    evS=runs[names[1]]["scenario"]
    N=evS.N

    objConv = zeros(numIt,P)
    lamConv = zeros(numIt,P)
    T = zeros(Klen,P)
    Sn = zeros(Klen*N,P)
    uSum = zeros(Klen,P)
    snSum = zeros(Klen,P)

    for i in 1:length(names)
        run=runs[names[i]]
        objConv[:,i]=run["convMetrics"].objVal
        lamConv[:,i]=run["convMetrics"].lam
        T[:,i]=run["solution"].Xt[:,run["convIt"]]
        uSum[:,i]=run["solution"].uSum[:,run["convIt"]]
        Sn[:,i]=run["solution"].Sn[:,run["convIt"]]
        snSum[:,i]=[sum(Sn[N*(k-1)+n,i] for n=1:N) for k in 1:Klen]# sum up across N
    end

    # ty=[match(r"central|d",names[i]).match for i in 1:length(names)]
    # form = [match(r"NL|N|relax",names[i]).match for i in 1:length(names)]
    #lamConv=DataFrame(hcat(names,transpose(lamConv)))
    lamConv=DataFrame(lamConv)
    names!(lamConv.colindex, map(Symbol, names))
    objConv=DataFrame(objConv)
    names!(objConv.colindex, map(Symbol, names))
    T=DataFrame(T)
    names!(T.colindex, map(Symbol, names))
    uSum=DataFrame(uSum)
    names!(uSum.colindex, map(Symbol, names))
    snSum=DataFrame(snSum)
    names!(snSum.colindex, map(Symbol, names))

    lamPlot=plot(lamConv,x=Row.index,y=Col.value,Geom.point,Geom.line, color=Col.index,Scale.y_log10,
    			Guide.xlabel("Iteration"), Guide.ylabel("2-norm Lambda gap",orientation=:vertical),
    			Theme(background_color=colorant"white",major_label_font_size=24pt,line_width=2pt,
    			minor_label_font_size=20pt,key_label_font_size=20pt,key_position=:top))

    objPlot=plot(objConv,x=Row.index,y=Col.value,Geom.point,Geom.line, color=Col.index,Scale.y_log10,
    			Guide.xlabel("Iteration"), Guide.ylabel("2-norm Objective Value Gap",orientation=:vertical),
    			Theme(background_color=colorant"white",major_label_font_size=24pt,line_width=2pt,
    			minor_label_font_size=20pt,key_label_font_size=12pt,key_position=:none))

    tempPlot=plot(layer(T,x=Row.index,y=Col.value,Geom.point,Geom.line,color=Col.index),
    		layer(yintercept=[evS.Tmax],Geom.hline(color=["red"],style=:dot),Theme(line_width=3pt)),
    		Guide.xlabel("Time"), Guide.ylabel("Xfrm Temp (K)",orientation=:vertical),
    		Coord.Cartesian(xmin=0,xmax=Klen),Theme(background_color=colorant"white",major_label_font_size=24pt,
    		minor_label_font_size=20pt,key_label_font_size=12pt,key_position=:none))

    uSumPlot=plot(layer(uSum,x=Row.index,y=Col.value,Geom.point,Geom.line,color=Col.index),
    		# layer(yintercept=[evS.Tmax],Geom.hline(color=["red"],style=:dot),Theme(line_width=3pt)),
    		Guide.xlabel("Time"), Guide.ylabel("Current Sum",orientation=:vertical),
            Guide.colorkey(title="",labels=names),
    		Coord.Cartesian(xmin=0,xmax=Klen),Theme(background_color=colorant"white",major_label_font_size=24pt,
    		minor_label_font_size=20pt,key_label_font_size=12pt,key_position=:none))

    target=zeros(Klen)
    for k in 1:Klen
        ind = evS.Kn.<=k
        target[k]=sum(evS.Snmin[ind])
    end

    snSumPlot=plot(layer(snSum,x=Row.index,y=Col.value,Geom.point,Geom.line,color=Col.index),
    		layer(x=1:Klen,y=target,Geom.line,Theme(line_width=3pt,default_color=colorant"red",line_style=:dot)),
    		Guide.xlabel("Time"), Guide.ylabel("SOC Sum",orientation=:vertical),
            Guide.colorkey(title="",labels=names),
    		Coord.Cartesian(xmin=0,xmax=Klen),Theme(background_color=colorant"white",major_label_font_size=24pt,
    		minor_label_font_size=20pt,key_label_font_size=12pt,key_position=:none))

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

    Rmax=DataFrame(Rmax)
    Ravg=DataFrame(Ravg)
    Rplot=plot(layer(Rmax,x=Row.index,y=Col.value,Geom.point,Geom.line,color=Col.index),
            Guide.xlabel("Time"), Guide.ylabel("R Max",orientation=:vertical),
            Guide.colorkey(title="",labels=names),
            Coord.Cartesian(xmin=0,xmax=Klen),Theme(background_color=colorant"white",major_label_font_size=24pt,
            minor_label_font_size=20pt,key_label_font_size=12pt,key_position=:none))

    p=vstack(lamPlot, objPlot, tempPlot, uSumPlot,snSumPlot)
    fName="compPlot.png"
    draw(PNG(path*fName, 28inch, 20inch),p)
end

function compareRunsTable(runs)
    compareTable = DataFrame(name=String[],time=Float64[],cLamDiff=Float64[],lamDiff=Float64[],
    cObjDiff=Float64[],objDiff=Float64[])
    for key in keys(runs)
        loadF=runs[key]
        timeT=loadF["time"]
        cm=loadF["convMetrics"]
        convIt=loadF["convIt"]
        ind= if convIt>0 convIt-1 else length(cm.lam) end
        stats = [key timeT minimum(cm.lam[1:ind]) cm.lam[ind-1]-cm.lam[ind-2] minimum(cm.objVal[1:ind]) cm.objVal[ind-1]-cm.objVal[ind-2]]
        push!(compareTable,stats)
    end

    return compareTable
end
