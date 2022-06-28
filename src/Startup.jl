using Distributed
using Test
using BenchmarkTools
using SharedArrays
const numCores = 10
if (nworkers() != numCores)
    rmprocs(procs()[2:end])
    addprocs(numCores)
end

@everywhere using MatrixNetworks
@everywhere using DifferentialEquations
@everywhere using LinearAlgebra
@everywhere using SparseArrays

include("../src/SimulateDynamics.jl")
include("../src/GenerateGraphs.jl")
include("../src/Partition.jl")
include("../src/ReduceNetwork.jl")
include("../src/EvaluateError.jl")
include("../src/LocalSearch.jl")
include("../src/visualize_landscape.jl")
include("../src/NoPlotLandscape.jl")

using Plots
using GraphRecipes
using Dates
