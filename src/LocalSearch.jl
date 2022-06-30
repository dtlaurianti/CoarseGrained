
#Function: geneticImprovement
#Parameters: A, MatrixNetwork to partition and run dynamics on
#            dynamical_function, the dynamic model to run on the network
#            partitions, the partitions we will evolve from
#            generations, the number of reproduction cycles to run
#            crossover_prob, the chance of a partition to randomly changing into a adjacent partition
#            mutation_prob, the chance of a partition to randomly changing into a adjacent partition
#            initial_condition, the initial variable values
#            dynamical_function, the dynamical model we will run getLoss using
#            tmax, when to end the model
#            dt, the length of the timesteps in our simulation
#Purpose: To find the partition with local minimum loss using a basic iterative approach
#Return value: returns a partition that is an approximate local minimum of the partition space
function geneticImprovement(A::MatrixNetwork, partitions::Array{Dict{Integer, Integer}}, generations::Integer, crossover_prob::Float64, mutation_prob::Float64, initial_condition::Vector, dynamical_function::Function, tmax::Number, dt::Number; function_args...)
    c = length(partitions)
    n = length(partitions[1])
    k = length(unique(values(partitions[1])))

    crossover_prob = pweights(1-crossover_prob, crossover_prob)
    mutation_prob = pweights(1-mutation_prob, mutation_prob)

    function_args = Dict(function_args)
    individuals = partitions
    for _=1:generations
        # reproduction phase
        loss_log = []
        # store the magnitude of loss for each partition
        for individual in individuals
            append!(loss, log(getLoss(A, individual, initial_condition, dynamical_function, tmax, dt, function_args...)))
        end
        loss_sum = sum(loss)
        prob = []
        # each partition reproduces with probability scaled to the magnitude of its loss
        for i = 1:c
            append!(prob, loss_log[i]/loss_sum)
        end
        weights = pweights(prob)
        children = []
        for i = 1:c
            append!(children, individuals[StatsBase.sample(1:c, weights)])
        end

        # crossing phase
        for i = 1:c
            swap = StatsBase.sample([0,1], crossover_prob, )
        end

        # mutation phase
        for i = 1:c
            if StatsBase.sample([0,1], mutation_prob) == 1
                children[c] =
        end
    end
    finalPartitions = []
    for individual in individuals
        push!(finalPartitions, individual)
    end
    return individuals

end

#=
#Function: dict_to_matrix
#Parameters: p, the partition to convert to matrix form
#            n, the number of nodes in the partition
#            k, the number of supernodes in the partition
#Purpose: To find the partition with local minimum loss using a basic iterative approach
#Return value: returns a partition that is an approximate local minimum of the partition space
function dict_to_matrix(p::Dict{Integer, Integer}, n::Integer, k::Integer)
    P = zeros(n, k)
    for i = 1:n
        P[i, p[i]] = 1
    end
    return P
end

#Function: matrix_to_dict
#Parameters: P, the partition to convert to ditionary form
#            n, the number of nodes in the partition
#            k, the number of supernodes in the partition
#Purpose: To find the partition with local minimum loss using a basic iterative approach
#Return value: returns a partition that is an approximate local minimum of the partition space
function matrix_to_dict(P::Matrix, n::Integer, k::Integer)
    p = Dict{Integer, Integer}()
    for i = 1:n
        for j = 1:k
            if P[i, j] == 1
                p[i] = j
            end
        end
    end
    return p
end
=#



function findLocalMinimum(A::MatrixNetwork, p::Dict{Integer, Integer}, depth::Integer, initial_condition::Vector, dynamical_function::Function, tmax::Number, dt::Number; function_args...)
    n = length(p)
    k = length(unique(values(p)))
    function_args = Dict(function_args)
    minloss = getLoss(A, p, initial_condition, dynamical_function, tmax, dt, function_args...)
    display(minloss)
    p2 = p
    # iteratively and greedily choose the least adjacent partition until we find a local minimum
    while true
        neighborhood = getNeighborhood(p, n, k, depth) #replace with radius-based search
        # find the adjacent partition with the least loss
        for neighbor in neighborhood
            loss = getLoss(A, neighbor, initial_condition, dynamical_function, tmax, dt, function_args...)#just use z values
            if loss < minloss
                p2 = neighbor
                minloss = loss
            end
        end
        # if we haven't found a better partition in the surrounding partitions we return out current partition
        if p2 == p
            break
        end
        # update the point we are at
        p = p2
        display(minloss)
    end
    return p
end

#Function: iterativeImprovement
#Parameters: A, MatrixNetwork to partition and run dynamics on
#            dynamical_function, the dynamic model to run on the network
#            p, the partition to start from
#            depth, how many adjacent partitions away our neighborhoods should look
#            initial_condition, the initial variable values
#            dynamical_function, the dynamical model we will run getLoss using
#            tmax, when to end the model
#            dt, the length of the timesteps in our simulation
#Purpose: To find the partition with local minimum loss using a basic iterative approach
#Return value: returns a partition that is a local minimum of the partition space
function iterativeImprovement(A::MatrixNetwork, p::Dict{Integer, Integer}, depth::Integer, initial_condition::Vector, dynamical_function::Function, tmax::Number, dt::Number; function_args...)
    n = length(p)
    k = length(unique(values(p)))
    function_args = Dict(function_args)
    minloss = getLoss(A, p, initial_condition, dynamical_function, tmax, dt, function_args...)
    display(minloss)
    p2 = p
    # iteratively and greedily choose the least adjacent partition until we find a local minimum
    while true
        neighborhood = getNeighborhood(p, n, k, depth)
        # find the adjacent partition with the least loss
        for neighbor in neighborhood
            loss = getLoss(A, neighbor, initial_condition, dynamical_function, tmax, dt, function_args...)
            if loss < minloss
                p2 = neighbor
                minloss = loss
            end
        end
        # if we haven't found a better partition in the surrounding partitions we return out current partition
        if p2 == p
            break
        end
        # update the point we are at
        p = p2
        display(minloss)
    end
    return p
end

#Function: getAdjacentPartitions
#Parameters: p, the partition to start from
#            n, the number of nodes
#            k, the number of supernodes
#Purpose: To return a list of the adjacent partitions
#Return value: returns a list of partitions adjacent to the input partition
function getAdjacentPartitions(p::Dict{Integer, Integer}, n::Integer, k::Integer)
    adjacents = []
    # apply every possible single node change and record the resulting partitions
    for i=1:n
        for j=1:k
            if p[i] != j
                append!(adjacents, [copy(p)])
                adjacents[end][i] = j
                if length(unique(values(adjacents[end]))) != k
                    pop!(adjacents)
                end
            end
        end
    end
    return adjacents
end

#Function: getNeighborhood
#Parameters: p, the partition to start from
#            n, the number of nodes
#            k, the number of supernodes
#            depth, the number of partitions away to look
#Purpose: To return a list of the neighboring partitions
#Return value: returns a list of partitions up to depth away from the input partition
function getNeighborhood(p::Dict{Integer, Integer}, n::Integer, k::Integer, depth::Integer)
    neighbors = Set([p])
    for _=1:depth
        for partition in neighbors
            union!(neighbors, getAdjacentPartitions(partition, n, k))
        end
    end
    pop!(neighbors, p)
    return neighbors
end

#Function: getAdjacentSample
#Parameters: p, the partition to start from
#            n, the number of nodes
#            k, the number of supernodes
#            s, the number of sample partitions
#Purpose: To return a sample list of the adjacent partitions
#Return value: returns a sample list of partitions adjacent to the input partition
function getAdjacentSample(p::Dict{Integer, Integer}, n::Integer, k::Integer, s::Integer=1)
    sample = []
    # change the supernode we partition one node into to create s new partitions adjacent to this partition
    for _=1:s
        append!(sample, [copy(p)])
        while true
            ns = rand(1:n)
            ks = rand(1:k)
            if p[ns] != ks
                sample[end][ns] = ks
                if length(unique(values(sample[end]))) != k
                    pop!(sample)
                    append!(sample, [copy(p)])
                else
                    break
                end
            end
        end
    end
    return sample
end
