using Distributed

const numCores = 10
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

include("../src/SimulateDynamics.jl")
include("../src/GenerateGraphs.jl")
include("../src/Partition.jl")
include("../src/ReduceNetwork.jl")
include("../src/EvaluateError.jl")
include("../src/LocalSearch.jl")
include("../src/VisualizeLandscape.jl")
include("../src/DataVis.jl")
include("../src/PyPlot.jl")
include("../src/AnalyzeLandscape.jl")


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
