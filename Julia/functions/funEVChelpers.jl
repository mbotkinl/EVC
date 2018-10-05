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
    K=size(runs[names[1]]["solution"].Xt)[1]
    P=length(names)
    evS=runs[names[1]]["scenario"]

    objConv = zeros(numIt,P)
    lamConv = zeros(numIt,P)
    T = zeros(K,P)
    uSum = zeros(K,P)

    for i in 1:length(names)
        run=runs[names[i]]
        objConv[:,i]=run["convMetrics"].objVal
        lamConv[:,i]=run["convMetrics"].lam
        T[:,i]=run["solution"].Xt[:,run["convIt"]]
        uSum[:,i]=run["solution"].uSum[:,run["convIt"]]
    end

    # ty=[match(r"central|d",names[i]).match for i in 1:length(names)]
    # form = [match(r"NL|N|relax",names[i]).match for i in 1:length(names)]
    #lamConv=DataFrame(hcat(names,transpose(lamConv)))
    lamConv=DataFrame(lamConv)
    objConv=DataFrame(objConv)
    T=DataFrame(T)
    uSum=DataFrame(uSum)

    lamPlot=plot(lamConv,x=Row.index,y=Col.value,Geom.point,Geom.line, color=Col.index,
    			Guide.xlabel("Iteration"), Guide.ylabel("2-norm Lambda gap",orientation=:vertical),
                Guide.colorkey(title="",labels=names),Scale.y_log10,
    			Theme(background_color=colorant"white",major_label_font_size=24pt,line_width=2pt,
    			minor_label_font_size=20pt,key_label_font_size=20pt,key_position=:top))

    objPlot=plot(objConv,x=Row.index,y=Col.value,Geom.point,Geom.line, color=Col.index,
    			Guide.xlabel("Iteration"), Guide.ylabel("2-norm Objective Value Gap",orientation=:vertical),
                Guide.colorkey(title="",labels=names),Scale.y_log10,
    			Theme(background_color=colorant"white",major_label_font_size=24pt,line_width=2pt,
    			minor_label_font_size=20pt,key_label_font_size=12pt,key_position=:none))

    tempPlot=plot(layer(T,x=Row.index,y=Col.value,Geom.point,Geom.line,color=Col.index),
    		layer(yintercept=[evS.Tmax],Geom.hline(color=["red"],style=:dot),Theme(line_width=3pt)),
    		Guide.xlabel("Iteration"), Guide.ylabel("Xfrm Temp (K)",orientation=:vertical),
            Guide.colorkey(title="",labels=names),
    		Coord.Cartesian(xmin=0,xmax=K),Theme(background_color=colorant"white",major_label_font_size=24pt,
    		minor_label_font_size=20pt,key_label_font_size=12pt,key_position=:none))

    uSumPlot=plot(layer(uSum,x=Row.index,y=Col.value,Geom.point,Geom.line,color=Col.index),
    		# layer(yintercept=[evS.Tmax],Geom.hline(color=["red"],style=:dot),Theme(line_width=3pt)),
    		Guide.xlabel("Iteration"), Guide.ylabel("Current Sum",orientation=:vertical),
            Guide.colorkey(title="",labels=names),
    		Coord.Cartesian(xmin=0,xmax=K),Theme(background_color=colorant"white",major_label_font_size=24pt,
    		minor_label_font_size=20pt,key_label_font_size=12pt,key_position=:none))


    p=vstack(lamPlot, objPlot, tempPlot, uSumPlot)
    fName="compPlot.png"
    draw(PNG(path*fName, 28inch, 14inch),p)
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
