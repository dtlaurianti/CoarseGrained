using Distributed
using Test
using BenchmarkTools
const numCores = 104
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
    include("../src/LocalSearch.jl")
end
include("../src/visualize_landscape.jl")

using Plots
using GraphRecipes
using Dates
