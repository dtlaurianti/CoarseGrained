#=
@testset "Partition_tests" begin
    c = 5
    G = gnp_graph(10, p=0.5, directed=false)
    Ps = generateRandomPartitions(10, 5, c)
    u = rand(10)

    gplts, cgplts = plotPartitions(G, Ps, u)
    plt = Plots.plot(gplts..., cgplts..., layout=(2,c))
    display(plt)
end
=#
@testset "Dynamics_tests" begin
    #c = 4
    G = gnp_graph(10, p=0.25, directed=true)
    #G = cycle_graph(10)
    #Ps = generateRandomPartitions(10, 5, c)
    P = generateRandomPartitions(10, 5, 1)[1]
    u = rand(10)
    #scale = 1/(2*log(c, 10))

    anim = animatePartitionDynamics(G, P, u)
    display(anim)
    #=
    gplts, cgplts, dplts = plotPartitionsDynamics(G, Ps, u)
    plt = Plots.plot(gplts..., cgplts..., dplts..., layout=(3, c), thickness_scaling=scale)
    display(plt)
    =#
end