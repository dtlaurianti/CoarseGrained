@testset "GetXYZ" begin
    numOriginalNodes = 300
    Part = generateRandomPartitions(numOriginalNodes, 100, 1000)
    #LG = line_graph(numOriginalNodes)
    GNP = gnp_graph(numOriginalNodes;p=0.5)
    #CM = cm_graph(numOriginalNodes, [1, 1, 1, 1, 2, 2, 2, 3, 3, 4])
    dt = now()

    DT = Dates.format(dt, "mm-dd_HH-MM-SS")
    timeString = "test" * DT
    #Uncomment one of the following depending on if you want the results to be saved
    #to a CSV file or nott
    #surfaceplots(Part, GNP, 25, save_to_string=string)
    @time GetXYZ(Part, GNP, numOriginalNodes, modelType=linear_model, save_to_string=timeString)
    println(string)

    #=In normal terminal, (with R installed)
    call Rscript --vanilla ~/Documents/GitHub/CoarseGrained/src/make_smoothdata.R string"
    to make smooth data. Then plot it with the function invokation below (justreplace the name)=#
    #plot_smoothed_surface("./data/visualization_data/test06-23_11-54-30_smooth.csv")
end