#=
@testset "neighbors_tests" begin
    p = generateRandomPartitions(10,5,1)[1]
    display(p)
    sample = getAdjacentSample(p, 10, 5, 10)
    display(sample)

    p2 = generateRandomPartitions(5,3,1)[1]
    display(p2)
    neighbors = getAdjacentPartitions(p2, 5, 3)
    display(neighbors)
    G = gnp_graph(5, p=.5)
    x = rand(5)
    p3 = iterativeImprovement(G, p2, 1, x, linear_model, 10, 0.1)
    display(p3)
    loss = getLoss(G, p3, x, linear_model, 10, 0.1)
    println("Iterative Loss: ", loss)

end
=#
#=
@testset "conversion_tests" begin
    p = generateRandomPartitions(10,5,1)[1]
    display(p)
    P = dict_to_matrix(p, 10, 5)
    display(P)
    p2 = matrix_to_dict(P, 10, 5)
    display(p2)
end
=#

@testset "genetic_tests" begin
    #=
    p1 = generateRandomPartitions(10,5,1)[1]
    p2 = generateRandomPartitions(10,5,1)[1]
    pc = supernodeBucketCross(p1, p2, 10, 5)
    println("Parent 1: ", p1)
    println("Parent 2: ", p2)
    println("Child: ", pc)

    pcm = randomWalkMutate(pc, 10, 5, 0.75)
    println("Mutant Child: ", pcm)
    =#
    #=
    gen1 = generateRandomPartitions(100, 50, 512)
    x = rand(100)
    GNP = gnp_graph(100)
    f = 10
    agg1 = agglomerationReduction(GNP, 50)
    gen1[1] = agg1

    @time genf = geneticImprovement(GNP, gen1, f, 0.9, x, linear_model, 10, 0.1)
    gen1_loss = 0
    for i = 1:512
        gen1_loss += getLoss(GNP, gen1[i], x, linear_model, 10, 0.1)
    end
    genf_loss = 0
    for i = 1:512
        genf_loss += getLoss(GNP, genf[i], x, linear_model, 10, 0.1)
    end
    println("Avg Loss Gen1: ", gen1_loss/512)
    println("Avg Loss Gen$f: ", genf_loss/512)
    =#
    c = 100
    n = 5
    k = 3
    f0 = 10
    f1 = 100
    f2 = 500
    f3 = 1000
    gen1 = generateRandomPartitions(n, k, c)
    x = rand(n)
    G = gnp_graph(n, p=.5)
    genf0 = geneticImprovement(G, gen1, f0, 0.1, x, linear_model, 10, 0.1)
    genf1 = geneticImprovement(G, genf0, f1-f0, 0.1, x, linear_model, 10, 0.1)
    genf2 = geneticImprovement(G, genf1, f2-f1, 0.1, x, linear_model, 10, 0.1)
    genf3 = geneticImprovement(G, genf2, f3-f2, 0.1, x, linear_model, 10, 0.1)
    gen1_losses = getLossBatch(G, gen1, x, linear_model, 10, 0.1)
    gen1_avg = sum(gen1_losses/c)
    genf0_losses = getLossBatch(G, genf0, x, linear_model, 10, 0.1)
    genf0_avg = sum(genf0_losses/c)
    genf1_losses = getLossBatch(G, genf1, x, linear_model, 10, 0.1)
    genf1_avg = sum(genf1_losses/c)
    genf2_losses = getLossBatch(G, genf2, x, linear_model, 10, 0.1)
    genf2_avg = sum(genf2_losses/c)
    genf3_losses = getLossBatch(G, genf3, x, linear_model, 10, 0.1)
    genf3_avg = sum(genf3_losses/c)

    plt = violin(1:100, gen1_losses, label = "Gen 1, Avg = $gen1_avg")
    violin!(plt, 1:100, genf0_losses, label = "Gen $f0, Avg = $genf0_avg")
    violin!(plt, 1:100, genf1_losses, label = "Gen $f1, Avg = $genf1_avg")
    violin!(plt, 1:100, genf2_losses, label = "Gen $f2, Avg = $genf2_avg")
    violin!(plt, 1:100, genf3_losses, label = "Gen $f3, Avg = $genf3_avg")
    #=
    iter_partitions = pmap(part->iterativeImprovement(G, part, 1, x, linear_model, 10, 0.1), gen1)
    iter_losses = pmap(part->getLoss(G, part, x, linear_model, 10, 0.1), iter_partitions)
    iter_avg = sum(iter_losses)

    violin!(plt, 1:100, iter_losses, label = "Greedy, Avg = $iter_avg")
    =#
    display(plt)
end

#=
@testset "iterative_tests" begin
    part = generateRandomPartitions(6,4,1)[1]
    G = gnp_graph(6, p=0.5)
    x = rand(6)
    part2 = iterativeImprovementDynamic(G, part, 1, x, linear_model, 10, 0.1)
end
=#
