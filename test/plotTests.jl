x = zeros(1)
y = []
times = []
Part1 = generateRandomPartitions(15, 10, 10);
GNP1 = gnp_graph(15;p=0.5)
trash = @elapsed(GetXYZ(Part1, GNP1, 15, modelType=linear_model))

for numOriginalNodes in 100:50:1000
    push!(x, numOriginalNodes)
    Part = generateRandomPartitions(numOriginalNodes, 10, 10);
    #LG = line_graph(numOriginalNodes)
    GNP = gnp_graph(numOriginalNodes;p=0.5)
    #CM = cm_graph(numOriginalNodes, [1, 1, 1, 1, 2, 2, 2, 3, 3, 4])
    push!(times, @elapsed(GetXYZ(Part, GNP, numOriginalNodes, modelType=linear_model)))
    push!(times, @elapsed(GetXYZ(Part, GNP, numOriginalNodes, modelType=linear_model)))
    push!(times, @elapsed(GetXYZ(Part, GNP, numOriginalNodes, modelType=linear_model)))
    push!(y, median(times))
end

pushfirst!(y, 0)
figure = Plots.plot(x, y, title="Time Complexity With Respect To numOriginalNodes", xlabel="numOriginalNodes",
ylabel = "time (seconds)")
png(figure, "OriginalNodesPlot.png")

x = zeros(1)
y = []
times = []

for numReducedNodes in 100:50:1000
    push!(x, numReducedNodes)
    Part = generateRandomPartitions(1000, numReducedNodes, 10);
    #LG = line_graph(numOriginalNodes)
    GNP = gnp_graph(1000;p=0.5)
    #CM = cm_graph(numOriginalNodes, [1, 1, 1, 1, 2, 2, 2, 3, 3, 4])
    push!(times, @elapsed(GetXYZ(Part, GNP, 1000, modelType=linear_model)))
    push!(times, @elapsed(GetXYZ(Part, GNP, 1000, modelType=linear_model)))
    push!(times, @elapsed(GetXYZ(Part, GNP, 1000, modelType=linear_model)))
    push!(y, median(times))
end

pushfirst!(y, 0)
figure = Plots.plot(x, y, title="Time Complexity With Respect To numReducedNodes", 
xlabel="numReducedNodes", ylabel = "time (seconds)")
png(figure, "ReducedNodesPlot.png")

x = zeros(1)
y = []
times = []

for numPartitions in 100:50:1000
    push!(x, numPartitions)
    Part = generateRandomPartitions(20, 10, numPartitions);
    #LG = line_graph(numOriginalNodes)
    GNP = gnp_graph(20;p=0.5)
    #CM = cm_graph(numOriginalNodes, [1, 1, 1, 1, 2, 2, 2, 3, 3, 4])
    push!(times, @elapsed(GetXYZ(Part, GNP, 20, modelType=linear_model)))
    push!(times, @elapsed(GetXYZ(Part, GNP, 20, modelType=linear_model)))
    push!(times, @elapsed(GetXYZ(Part, GNP, 20, modelType=linear_model)))
    push!(y, median(times))
end

pushfirst!(y, 0)
figure = Plots.plot(x, y, title="Time Complexity With Respect To numPartitions", 
xlabel="numPartitions", ylabel = "time (seconds)")
png(figure, "numPartitionsPlot.png")