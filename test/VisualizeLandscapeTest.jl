@testset "dict_to_array_test" begin
    Part = generateRandomPartitions(10, 7, 10)
    Partition1 = Part[1]
    #display(Partition1)
    a = dict_to_array(Part)
    b = partition_dict_to_array(Partition1)
    display(a)
    display(b)
    @test a[1] == b
end

@testset "variation_of_information_test" begin
    Part = generateRandomPartitions(10, 7, 10)
    Arr = dict_to_array(Part)
    display(variation_of_information(Arr[1], Arr[2]))
end

@testset "surfaceplots" begin
    Part = generateRandomPartitions(10, 7, 1000)
    LG = line_graph(10)
    GNP = gnp_graph(10;p=0.5)
    surfaceplots(Part, LG)
end