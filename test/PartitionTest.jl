#=
@testset "Partition Functions Test" begin
    print("generateRandomPartitions(10, 5, 10)): ")
    display(generateRandomPartitions(10, 5, 10))
    print("spectralClustering(gnp_graph(10), 5): ")
    display(spectralClustering(gnp_graph(10), 5))
    print("agglomerationReduction(gnp_graph(10), 5): ")
    display(agglomerationReduction(gnp_graph(10), 5))
    print("exhaustivePartition(5): ")
    display(exhaustivePartition(5))
    print("kpartition(5,4): ")
    display(kPartition(5, 4))
end
=#
@testset "Efficiency Testing" begin
    display(@benchmark generateRandomPartitions(10, 5, 100))
    println(nprocs())
    display(@benchmark generateRandomPartitionsFast(10, 5, 100))
end
