x = []
y = []
times = []
Part1 = generateRandomPartitions(15, 10, 10);
GNP1 = gnp_graph(15;p=0.5)
for i in 1:5
    trash0 = @elapsed(GetXYZ(Part1, GNP1, 15, modelType=linear_model))
end

for numOriginalNodes in 100:100:15000
    push!(x, numOriginalNodes)
    Part = generateRandomPartitions(numOriginalNodes, 10, 10);
    #LG = line_graph(numOriginalNodes)
    GNP = gnp_graph(numOriginalNodes;p=0.5)
    #CM = cm_graph(numOriginalNodes, [1, 1, 1, 1, 2, 2, 2, 3, 3, 4])
    for i in 1:4
        push!(times, @elapsed(GetXYZ(Part, GNP, numOriginalNodes, modelType=linear_model)))
    end
    push!(y, median(times))
end

figure = Plots.plot(x, y, title="Time Complexity With Respect To numOriginalNodes", xlabel="numOriginalNodes",
ylabel = "time (seconds)", margin=9Plots.mm, ylims = (0, Inf))
png(figure, "OriginalNodesPlot.png")

df = DataFrame(["x" => x, "y" => y])
loc = "./data/time_complexity_data/original_nodes.csv"
CSV.write(loc, df)

x = []
y = []
times = []

for i in 1:5
    trash2 = @elapsed(GetXYZ(Part1, GNP1, 15, modelType=linear_model))
end

for numReducedNodes in 100:100:15000
    push!(x, numReducedNodes)
    Part = generateRandomPartitions(15100, numReducedNodes, 10);
    #LG = line_graph(numOriginalNodes)
    GNP = gnp_graph(15100;p=0.5)
    #CM = cm_graph(numOriginalNodes, [1, 1, 1, 1, 2, 2, 2, 3, 3, 4])
    for i in 1:4
        push!(times, @elapsed(GetXYZ(Part, GNP, 15100, modelType=linear_model)))
    end
    push!(y, median(times))
end

figure = Plots.plot(x, y, title="Time Complexity With Respect To numReducedNodes", 
xlabel="numReducedNodes", ylabel = "time (seconds)", margin=9Plots.mm, ylims = (0, Inf))
png(figure, "ReducedNodesPlot.png")

df2 = DataFrame(["x" => x, "y" => y])
loc = "./data/time_complexity_data/reduced_nodes.csv"
CSV.write(loc, df2)

x = []
y = []
times = []

for i in 1:5
    trash3 = @elapsed(GetXYZ(Part1, GNP1, 15, modelType=linear_model))
end

for numPartitions in 100:100:15100
    push!(x, numPartitions)
    Part = generateRandomPartitions(20, 10, numPartitions);
    #LG = line_graph(numOriginalNodes)
    GNP = gnp_graph(20;p=0.5)
    #CM = cm_graph(numOriginalNodes, [1, 1, 1, 1, 2, 2, 2, 3, 3, 4])
    for i in 1:4
        push!(times, @elapsed(GetXYZ(Part, GNP, 20, modelType=linear_model)))
    end
    push!(y, median(times))
end

figure = Plots.plot(x, y, title="Time Complexity With Respect To numPartitions", 
xlabel="numPartitions", ylabel = "time (seconds)", margin=9Plots.mm, ylims = (0, Inf))
png(figure, "numPartitionsPlot.png")

df3 = DataFrame(["x" => x, "y" => y])
loc = "./data/time_complexity_data/numPartitions.csv"
CSV.write(loc, df3)