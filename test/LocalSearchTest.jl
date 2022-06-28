@testset "neighbors_tests" begin
    p = generateRandomPartitions(10,5,1)[1]
    display(p)
    sample = getAdjacentSample(p, 10, 5, 10)
    display(sample)

    p2 = generateRandomPartitions(5,3,1)[1]
    display(p2)
    neighbors = getAdjacentPartitions(p2, 5, 3)
    display(neighbors)
    G = gnp_graph(5, p=.5)
    x = rand(5)
    p3 = iterativeImprovement(G, p2, 1, x, linear_model, 10, 0.1)
    display(p3)
end
