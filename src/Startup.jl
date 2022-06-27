using Distributed
using Test
using BenchmarkTools
const numCores = 6
if (nworkers() != numCores)
    rmprocs(procs()[2:end])
    addprocs(numCores)
end
@everywhere begin
    include("../src/SimulateDynamics.jl")
    include("../src/GenerateGraphs.jl")
    include("../src/Partition.jl")
    include("../src/ReduceNetwork.jl")
    include("../src/EvaluateError.jl")
    include("../src/visualize_landscape.jl")
    include("../src/LocalSearch.jl")
end

using Plots
using GraphRecipes
using Dates
