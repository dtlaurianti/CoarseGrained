using ProgressBars
x = []
y = []
times = []
Part1 = generateRandomPartitions(100, 50, 500);
GNP1 = gnp_graph(100;p=0.5)
for i in 1:5
    trash0 = @elapsed(GetXYZ(Part1, GNP1, 100, modelType=linear_model))
end

println("Reduced Nodes Progress:")
R = 200:200:10000
O = 400:200:10200
for numReducedNodes, numOriginal in R, O
    numOriginal = 400
    GNP = gnp_graph(numOriginal; p=0.5)
    Part = generateRandomPartitions(numOriginal, numReducedNodes, 10);
    push!(x, numReducedNodes)
    for i in 1:4
        push!(times, @elapsed(getXYZ(Part, GNP, numOriginal, modelType=linear_model)))
    end
    push!(y, median(times))
end

figure = Plots.plot(x, y, title="Time Complexity With Respect To numReducedNodes",
xlabel="numReducedNodes", ylabel = "time (seconds)", margin=9Plots.mm, ylims = (0, Inf))
png(figure, "ReducedNodesPlot.png")

df2 = DataFrame(["x" => x, "y" => y])
loc = "./data/time_complexity_data/reduced_nodes.csv"
CSV.write(loc, df2)
println("Complete!")

x = []
y = []
times = []

for i in 1:5
    trash2 = @elapsed(GetXYZ(Part1, GNP1, 100, modelType=linear_model))
end

println("Original Nodes Progress:")
for numOriginalNodes in ProgressBar(200:200:20000)
    push!(x, numOriginalNodes)
    Part = generateRandomPartitions(numOriginalNodes, 10, 10);
    #LG = line_graph(numOriginalNodes)
    GNP = gnp_graph(numOriginalNodes;p=0.5)
    #CM = cm_graph(numOriginalNodes, [1, 1, 1, 1, 2, 2, 2, 3, 3, 4])
    for i in 1:4
        push!(times, @elapsed(getXYZ(Part, GNP, numOriginalNodes, modelType=linear_model)))
    end
    push!(y, median(times))
end

figure = Plots.plot(x, y, title="Time Complexity With Respect To numOriginalNodes", xlabel="numOriginalNodes",
ylabel = "time (seconds)", margin=9Plots.mm, ylims = (0, Inf))
png(figure, "OriginalNodesPlot.png")

df = DataFrame(["x" => x, "y" => y])
loc = "./data/time_complexity_data/original_nodes.csv"
CSV.write(loc, df)
println("Complete!")

x = []
y = []
times = []

for i in 1:5
    trash3 = @elapsed(GetXYZ(Part1, GNP1, 100, modelType=linear_model))
end

println("Partitions Progress:")
for numPartitions in ProgressBar(200:200:20000)
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
println("Complete!")
