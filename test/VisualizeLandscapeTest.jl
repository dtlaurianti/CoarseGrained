@testset "dict_to_array_test" begin
    Part = generateRandomPartitions(10, 7, 10)
    Partition1 = Part[1]
    display(Partition1)
    display(dict_to_array(Part))
    display(partition_dict_to_array(Partition1))
    
end