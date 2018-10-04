# helper functions for EVC code
using JLD2
using DataFrames

function saveRun(path::String, filename::String, time::Float64, scenario::scenarioStruct, solution, convMetrics=convMetricsStruct(), convIt=0)
    @save(path*filename*".jld","time", time, "scenario", scenario, "solution", solution,
    "convMetrics", convMetrics, "convIt", convIt)
end

function clips()
    t=clipboard()
    return t[2:length(t)-1]
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
    K=length(runs[names[1]]["convMetrics"].lam)
    P=length(names)
    lamConv = DataFrame()
    for name in names
        run=runs[name]
        lamConv[:name]=run["convMetrics"].lam
    end
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
        stats = [file timeT minimum(cm.lam[1:ind]) cm.lam[ind-1]-cm.lam[ind-2] minimum(cm.objVal[1:ind]) cm.objVal[ind-1]-cm.objVal[ind-2]]
        push!(compareTable,stats)
    end

    return compareTable
end
