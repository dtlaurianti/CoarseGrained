G = gnp_graph(20)
P = agglomerationReduction(G, 10)
x = rand(20)

@testset "getLoss_test" begin
    getLoss(G, P, x, linear_model, 10, .1; ω=rand(20), γ=0.1)
end
#=
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

    plt1 = Plots.plot(n,xlim=(0,1000),ylim=(0,25),title="line_graph linear_model", legend=false)
    plt2 = Plots.plot(cn,xlim=(0,1000),ylim=(0,25),title="line_graph linear_model compressed",legend=false)
    agg = aggregateTimeSeries(sol, Partition1)
    plt3 = Plots.plot(cn,xlim=(0,1000),ylim=(0,25),title="aggregate error",legend=false)
    anim2 = @animate for i=1:1000
        push!(plt1,sol[i])
        push!(plt2,csol[i])
        push!(plt3,agg[:,i])
        plt = Plots.plot(plt1, plt2, plt3, layout = (3, 1))
    end every 10

    #display(gif(anim1))
    display(gif(anim2))
end
=#
#=
Part = generateRandomPartitions(10, 7, 10)
Partition1 = Part[1]
LG = line_graph(10)

@testset "Efficiency Testing" begin
    display(@benchmark getLoss(LG, Partition1, ones(10), linear_model, 10, 0.01, ϵ=-0.3))
    @profile getLoss(LG, Partition1, ones(10), linear_model, 10, 0.01, ϵ=-0.3)
    pprof(;webport=58697)
end

=#
