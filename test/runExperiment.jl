#include("../src/ServerStartup.jl")
include("../src/Startup.jl")

#=
This code is inteded to be a simple way to run an entire experiment from one file.
Simply fill in the parameters below to match whatever kind of experiment you'd like 
to run, then run this file from the REPL with the following invocation: 
include("test/runExperiment.jl")
=#

#Options: linear_model, SI_model, SIS_model, kuramoto_model, LotkaVolterra_model, linear_opinions, nonlinear_opinions
model = SIS_model

# Enter as a string. Options: Line, Cycle, Grid, GNP, SBM, CM 
NetworkType = "SBM"

numOriginalNodes = 20
numReducedNodes = 10
numPartitions = 100

#If using GNP (ignore if not using GNP):
GNPProbability = 0.1

#If using SBM (ignore if not using SBM):
#cint must be less than half of the number of nodes and cext must be less than cint
cint = 8
cext = 4

#if using configuration model (CM) (ignore if not using CM):
degreeArray = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5]

#true if you'd like to plot, false if not:
plot=true

#Radius-based search parameters:
xRadius = 1.0
yRadius = 0.5
#leave as empty string if you want a random starting partition to be used
startingPartition = ""

#That's it! Don't modify anything below the line
#============================================================================================#

Part = generateRandomPartitions(numOriginalNodes, numReducedNodes, numPartitions)
SBMsize = trunc(Int, numOriginalNodes/2)

if NetworkType == "Line"
    Network = line_graph(numOriginalNodes)
elseif NetworkType == "Cycle"
    Network = cycle_graph(numOriginalNodes)
elseif NetworkType == "Grid"
    Network = grid_graph(numOriginalNodes)
elseif NetworkType == "GNP"
    Network = gnp_graph(numOriginalNodes;p=GNPProbability)
elseif NetworkType == "SBM"
    Network = MatrixNetwork(sparse(stochastic_block_model(cint, cext, [SBMsize, SBMsize])))
elseif NetworkType == "CM"
    Network = cm_graph(numOriginalNodes, degreeArray)
else
    return "Unsupported Network Type. Check for Typos."
end 

dt = now()
DT = Dates.format(dt, "mm-dd_HH-MM-SS")
timeString = "test" * DT
dataString = "PART" * timeString
localMinPath = "data/visualization_data/" * dataString * ".csv"
surfaceplots(Part, Network, numOriginalNodes, modelType=model, save_to_string=timeString, plotting=plot)
if startingPartition != ""
    findLocalMinimum(localMinPath, xRadius, yRadius, startingPartition=startingPartition)
else
    findLocalMinimum(localMinPath, xRadius, yRadius)
end