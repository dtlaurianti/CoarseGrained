@testset "aggregateTimeSeries_test" begin
    Part = generateRandomPartitions(10, 7, 10)
    Partition1 = Part[1]
    MatNet = line_graph(5)
    x = [1, 2, 3, 4, 5]
    sol = simulateODEonGraph(MatNet, x)
    display(aggregateTimeSeries(sol, Partition1))
end