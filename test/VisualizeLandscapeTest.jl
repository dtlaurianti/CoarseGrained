
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
    Part = generateRandomPartitions(25, 10, 100)
    #LG = line_graph(10)
    GNP = gnp_graph(25;p=0.5)
    #CM = cm_graph(10, [1, 1, 1, 1, 2, 2, 2, 3, 3, 4])
    dt = now()
    DT = Dates.format(dt, "mm-dd_HH-MM-SS")
    string = "test" * DT
    #Uncomment one of the following depending on if you want the results to be saved
    #to a CSV file or not
    #surfaceplots(Part, GNP, 25, save_to_string=string)
    surfaceplots(Part, GNP, 25, modelType=SIS_model)
    println(string)

    #=In normal terminal, (with R installed)
    call Rscript --vanilla ~/Documents/GitHub/CoarseGrained/src/make_smoothdata.R string"
    to make smooth data. Then plot it with the function invokation below (justreplace the name)=#
    #plot_smoothed_surface("../data/visualization_data/test06-16_12-01-35_smooth.csv")
end


@testset "Efficiency Testing" begin
    #display(@benchmark generateRandomPartitions(15, 10, 500))
    #display(@benchmark line_graph(15))
    #display(@benchmark surfaceplots(generateRandomPartitions(15, 10, 500), line_graph(15), 15))
end
