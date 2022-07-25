#=
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
=#
#=
@testset "variation_of_information_test" begin
    Part = generateRandomPartitions(10, 7, 10)
    Arr = dict_to_array(Part)
    display(variation_of_information(Arr[1], Arr[2]))
end
=#

@testset "surfaceplots" begin
    
    numOriginalNodes = 1000
    Part = generateRandomPartitions(numOriginalNodes, 650, 7000)
    #LG = line_graph(numOriginalNodes)
    #GNP = gnp_graph(numOriginalNodes;p=0.1)
    #CM = cm_graph(numOriginalNodes, [1, 1, 1, 1, 2, 2, 2, 3, 3, 4])
    SBM = MatrixNetwork(sparse(stochastic_block_model(300, 150, [500, 500])))
    dt = now()

    DT = Dates.format(dt, "mm-dd_HH-MM-SS")
    timeString = "test" * DT
    #Uncomment one of the following depending on if you want the results to be saved
    #to a CSV file or not
    surfaceplots(Part, SBM, numOriginalNodes, modelType=linear_model, save_to_string=timeString, plotting=true)
    
    #=In normal terminal, (with R installed)
    call Rscript --vanilla ~/Documents/GitHub/CoarseGrained/src/make_smoothdata.R string"
    to make smooth data. Then plot it with the function invokation below (justreplace the name)=#
    #=
    plt = plot_smoothed_surface("./data/visualization_data/test07-19_15-47-54.csv")
    
    df = DataFrame(CSV.File("./data/visualization_data/PARTtest07-14_10-20-00.csv"))
    x = df.x
    y = df.y
    z = df.z
    parts = df.partition
    println(x, y, z, parts)
    plot!(plt, )
    =#
    
end


#=
@testset "optimization plots" begin
    numOriginalNodes = 8
    numReducedNodes = 4
    partitions = kPartition(numOriginalNodes, numReducedNodes)
    SBM = sbm_graph(numOriginalNodes)
    x,y,z = surfaceplots(partitions, SBM, numOriginalNodes, save_to_string="8_4_full_space", plotting=true)
end
=#
#=
@testset "Efficiency Testing" begin
    #display(@benchmark generateRandomPartitions(15, 10, 500))
    #display(@benchmark line_graph(15))
    #display(@benchmark surfaceplots(generateRandomPartitions(15, 10, 500), line_graph(15), 15))
end
=#
