@testset "neighbors_tests" begin
    p = generateRandomPartitions(10,5,1)[1]
    display(p)
    sample = getNeighborSample(p, 10, 5, 10)
    display(sample)

    p2 = generateRandomPartitions(5,3,1)[1]
    display(p2)
    neighbors = getNeighbors(p2, 5, 3)
    display(neighbors)
end
