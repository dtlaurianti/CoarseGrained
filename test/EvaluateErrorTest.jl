using Plots
using GraphRecipes
include("../src/ReduceNetwork.jl")
@testset "aggregateTimeSeries_test" begin
    Random.seed!(trunc(Int, time() * 1000000))
    n = 50
    cn = trunc(Int64, n/2)
    MatNet = gnp_graph(n)
    x = rand(n)
    Partition1 = agglomerationReduction(MatNet, cn)
    display(Partition1)
    sol = simulateODEonGraph(MatNet, x)
    cMatNet = compressAdjacencyMatrix(MatNet, Partition1)
    cx = compressInitialCondition(x, Partition1)
    csol = simulateODEonGraph(cMatNet, cx)

    method=:stress
    scale = 1/log(n)
    colors = collect(values(sort(Partition1)))
    gplt1 = graphplot(sparse(MatNet), title="uncompressed initial", method=method,
        node_weights=x, markercolor = colors, thickness_scaling=scale, palette=distinguishable_colors(n))
    gplt2 = graphplot(sparse(MatNet), title="uncompressed final", method=method,
        node_weights=sol[1000], markercolor = colors, thickness_scaling=scale, palette=distinguishable_colors(n))
    gplt1c = graphplot(sparse(cMatNet), title="compressed initial", method=method,
        node_weights=cx, markercolor = Vector(1:cn), thickness_scaling=scale, palette=distinguishable_colors(cn))
    gplt2c = graphplot(sparse(cMatNet), title="compressed final", method=method,
        node_weights=csol[1000], markercolor = Vector(1:cn), thickness_scaling=scale, palette=distinguishable_colors(cn))
    display(gplt1)
    display(gplt2)
    display(gplt1c)
    display(gplt2c)

    plt1 = plot(n,xlim=(0,1000),ylim=(0,25),title="line_graph linear_model", legend=false)
    plt2 = plot(cn,xlim=(0,1000),ylim=(0,25),title="line_graph linear_model compressed",legend=false)
    agg = aggregateTimeSeries(sol, Partition1)
    plt3 = plot(cn,xlim=(0,1000),ylim=(0,25),title="aggregate error",legend=false)
    anim2 = @animate for i=1:1000
        push!(plt1,sol[i])
        push!(plt2,csol[i])
        push!(plt3,agg[:,i])
        plt = plot(plt1, plt2, plt3, layout = (3, 1))
    end every 10

    #display(gif(anim1))
    display(gif(anim2))
end
