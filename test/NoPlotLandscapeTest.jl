@testset "GetXYZ" begin
    Part = generateRandomPartitions(100, 50, 100)
    #LG = line_graph(10)
    GNP = gnp_graph(100;p=0.5)
    #CM = cm_graph(10, [1, 1, 1, 1, 2, 2, 2, 3, 3, 4])
    dt = now()

    DT = Dates.format(dt, "mm-dd_HH-MM-SS")
    timeString = "test" * DT
    #Uncomment one of the following depending on if you want the results to be saved
    #to a CSV file or not
    #surfaceplots(Part, GNP, 25, save_to_string=string)
    GetXYZ(Part, GNP, 100, modelType=kuramoto_model, save_to_string=timeString)
    println(string)

    #=In normal terminal, (with R installed)
    call Rscript --vanilla ~/Documents/GitHub/CoarseGrained/src/make_smoothdata.R string"
    to make smooth data. Then plot it with the function invokation below (justreplace the name)=#
    #plot_smoothed_surface("./data/visualization_data/test06-23_11-54-30_smooth.csv")
end