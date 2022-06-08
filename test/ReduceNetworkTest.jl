using Plots
using GraphRecipes

@userplot NetworkPlot
@recipe function f(np::NetworkPlot)
end

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