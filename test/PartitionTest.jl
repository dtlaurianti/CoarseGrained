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
    #display(@benchmark generateRandomPartitions(10, 5, 1000))
    #display(@benchmark agglomerationReduction(gnp, 50) setup=(gnp=gnp_graph(100)) seconds=10)
    display(@benchmark exhaustivePartition(3))
    @profiler begin
        for i=1:100
            exhaustivePartition(3)
        end
    end
end