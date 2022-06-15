using Distributed
rmprocs(procs()[2:end])
addprocs(6)
@everywhere begin
    using Test
    using BenchmarkTools
    include("../src/SimulateDynamics.jl")
    include("../src/GenerateGraphs.jl")
    include("../src/ReduceNetwork.jl")
    include("../src/Partition.jl")
    include("../src/EvaluateError.jl")
    include("../src/visualize_landscape.jl")
end
#include("SimulateDynamicsTest.jl")
#include("GenerateGraphsTest.jl")
#include("ReduceNetworkTest.jl")
#include("EvaluateErrorTest.jl")
#include("VisualizeLandscapeTest.jl")
include("PartitionTest.jl")
