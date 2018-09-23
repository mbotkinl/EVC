# helper functions for EVC code

function saveRun(path::String, filename::String, time::Float64, scenario::scenarioStruct, solution, convMetrics=convMetricsStruct(), convIt=0)
    JLD.save(path*filename*".jld","time", time, "scenario", scenario, "solution", solution,
    "convMetrics", convMetrics, "convIt", convIt)
end


function clips()
end
