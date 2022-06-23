include("../src/Startup.jl")
println("Beginning Testing...")
#include("SimulateDynamicsTest.jl")
#include("GenerateGraphsTest.jl")
#include("ReduceNetworkTest.jl")
#include("EvaluateErrorTest.jl")
#include("VisualizeLandscapeTest.jl")
#include("PartitionTest.jl")
#include("RunMultiModelsTest.jl")

println("Testing Complete.")
