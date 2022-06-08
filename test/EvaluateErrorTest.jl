using Plots
using GraphRecipes
include("../src/ReduceNetwork.jl")
@testset "aggregateTimeSeries_test" begin
    MatNet = line_graph(10)
    x = [1, 2, 3, 4, 5, 5, 4, 3, 2, 1]
    Partition1 = spectralClustering(MatNet, 7)
    sol = simulateODEonGraph(MatNet, x)
    cMatNet = compressAdjacencyMatrix(MatNet, Partition1)
    cx = compressInitialCondition(x, Partition1)
    csol = simulateODEonGraph(cMatNet, cx)


    plt1 = plot(10,xlim=(0,1000),ylim=(0,5),title="line_graph linear_model", legend=false)
    plt2 = plot(7,xlim=(0,1000),ylim=(0,5),title="line_graph linear_model compressed",legend=false)
    agg = aggregateTimeSeries(sol, Partition1)
    plt3 = plot(7,xlim=(0,1000),ylim=(0,5),title="aggregate error",legend=false)
    display(agg)
    anim1 = @animate for i=1:1000
        push!(plt1,sol[i])
        push!(plt2,csol[i])
        push!(plt3,agg[:,i])
        plt = plot(plt1, plt2, plt3, layout = (3, 1))
    end every 10
    display(gif(anim1))
end
