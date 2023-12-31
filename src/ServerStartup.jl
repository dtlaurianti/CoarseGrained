using Distributed
using Test
using BenchmarkTools
using SharedArrays
const numCores = 104
if (nworkers() != numCores)
    rmprocs(procs()[2:end])
    addprocs(numCores)
end

@everywhere using MatrixNetworks
@everywhere using DifferentialEquations
@everywhere using LinearAlgebra
@everywhere using SparseArrays
@everywhere using StatsBase
@everywhere using CSV
@everywhere using DataFrames

using Test
using BenchmarkTools
using Plots
using StatsPlots
using GraphRecipes
using Dates
using ProgressBars
using CSV
using DataFrames
using NetworkLayout

include("../src/SimulateDynamics.jl")
include("../src/GenerateGraphs.jl")
include("../src/Partition.jl")
include("../src/ReduceNetwork.jl")
include("../src/EvaluateError.jl")
include("../src/LocalSearch.jl")
include("../src/VisualizeLandscape.jl")
