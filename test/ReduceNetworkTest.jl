
#=
@testset "getSupernodeSizes_test" begin
    Part = generateRandomPartitions(10, 7, 10)
    display(Part[1])
    display(Part[6])
    display(getSupernodeSizes(Part[1]))
    display(getSupernodeSizes(Part[6]))
end

@testset "compressInitialCondition_test" begin
    initialConditions = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    #initialConditions = rand(10)
    Part = generateRandomPartitions(10, 7, 10)
    Partition1 = Part[1]
    display(Partition1)
    display(getSupernodeSizes(Partition1))
    display(compressInitialCondition(initialConditions, Partition1))
end

@testset "compressAdjacencyMatrix_test" begin
    Part = generateRandomPartitions(10, 7, 10)
    Partition1 = Part[1]
    MatNet = line_graph(10)
    display(sparse(MatNet))
    gplt1 = graphplot(sparse(MatNet))
    display(gplt1)
    MatNet2 = compressAdjacencyMatrix(MatNet, Partition1)
    display(sparse(MatNet2))
    gplt2 = graphplot(sparse(MatNet2))
    display(gplt2)
end
=#
@testset "presentation_visuals" begin
    n = 10000
    cn = 10
    G = sbm_graph(n)
    x = rand(n)
    cx = rand(cn)
    P = agglomerationReduction(G, cn)
    CG = compressAdjacencyMatrix(G, P)


    colors = collect(values(sort(P)))

    gplt = graphplot(sparse(G), thickness_scaling=scale, markercolor = colors, node_weights=x, palette=distinguishable_colors(n))
    @time sol = simulateODEonGraph(G, x, ϵ=0.5, dynamical_function=linear_model)
    plt = Plots.plot(sol,title="sbm_graph linear_model", palette=distinguishable_colors(n), seriescolor = transpose(colors))
    display(Plots.plot(plt, gplt))

    gplt2 = graphplot(sparse(CG), thickness_scaling=scale, markercolor = Vector(1:cn), node_weights=cx, palette=distinguishable_colors(cn))
    @time sol2 = simulateODEonGraph(CG, cx, ϵ=0.5, dynamical_function=linear_model)
    plt2 = Plots.plot(sol2,title="sbm_graph linear_model compressed", palette=distinguishable_colors(cn))
    display(Plots.plot(plt2, gplt2))

end
