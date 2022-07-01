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
end

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
    p1 = generateRandomPartitions(10,5,1)[1]
    p2 = generateRandomPartitions(10,5,1)[1]
    pc = supernodeBucketCross(p1, p2, 10, 5)
    println("Parent 1: ", p1)
    println("Parent 2: ", p2)
    println("Child: ", pc)

    pcm = randomWalkMutate(pc, 10, 5, 0.75)
    println("Mutant Child: ", pcm)

    gen1 = generateRandomPartitions(100, 50, 512)
    x = rand(100)
    GNP = gnp_graph(100)

    @time gen3 = geneticImprovement(GNP, gen1, 16, 0.75, x, linear_model, 10, 0.1)
    @time gen3fast = geneticImprovementFast(GNP, gen1, 16, 0.75, x, linear_model, 10, 0.1)
    #=
    for part in gen1
        println("Gen 1: ", part)
    end
    for part in gen3
        println("Gen 3: ", part)
    end
    =#
end
