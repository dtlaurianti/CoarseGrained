using Distributed
rmprocs(procs()[2:end])
const numCores = 6
addprocs(numCores)
@everywhere begin
    using Test
    using BenchmarkTools
    include("../src/SimulateDynamics.jl")
    include("../src/GenerateGraphs.jl")
    include("../src/Partition.jl")
    include("../src/ReduceNetwork.jl")
    include("../src/EvaluateError.jl")
    include("../src/visualize_landscape.jl")
end
