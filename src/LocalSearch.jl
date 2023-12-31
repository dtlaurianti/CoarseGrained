@everywhere using DataStructures


#Function: supernodeBucketCross
#Parameters: parent1, the partition to cross
#            parent2, the partition to cross
#            n, the number of nodes in the partitions
#            k, the number of supernodes in the partitions
#Purpose: To cross-breed two partitions to produce a partition with some qualities of both
#Return value: returns a partition that is a combination of parent1 and parent2
@everywhere function supernodeBucketCross(parent1::Dict{Integer, Integer}, parent2::Dict{Integer, Integer}, n::Integer, k::Integer)
    # store the parent nodes in reversed vector notation
    buckets1 = Vector{Pair}()
    for (node,supernode) in parent1
        push!(buckets1, supernode=>node)
    end
    buckets2 = Vector{Pair}()
    for (node,supernode) in parent2
        push!(buckets2, supernode=>node)
    end

    # identically partitioned nodes are automatically placed in the child partition
    child = intersect(buckets1, buckets2)
    # differently partitioned nodes will have the nodes picked randomly from either to be placed in the child partition
    buckets = symdiff(buckets1, buckets2)
    # priority queue to control the order that supernode buckets are pulled from
    pq = PriorityQueue{Integer, Integer}()
    # initialize the priority queue with the supernodes and the number of nodes in their buckets
    nodecounts = [0 for _=1:k]
    samecounts = [0 for _=1:k]
    for pair in buckets
        nodecounts[pair.first] += 10
    end
    # and count the number of nodes in each child supernode
    for pair in child
        nodecounts[pair.first] += 10
        samecounts[pair.first] += 10
        # give tiebreaker to supernodes that haven't recieved a node yet
        if nodecounts[pair.first] % 10 == 0
            nodecounts[pair.first] += 1
        end
    end
    # if there are nodes in the the buckets, put those buckets in the priority queue
    for i = 1:k
        if nodecounts[i] > 1 + samecounts[i]
            enqueue!(pq, i=>nodecounts[i])
        end
    end

    # until we have assigned all nodes to a supernode in the child,
    # use the priority queue to select a bucket to pull from and add that node to the child partition
    while !isempty(pq)
        # choose which supernode bucket to select from using the priority queue
        supernode = peek(pq)
        # get the bucket of nodes corresponding to that supernode
        nodes = filter(pair->pair.first==supernode.first, buckets)
        # choose a node at random from the bucket
        node = rand(nodes)
        # add that node/supernode pair to the child partition
        push!(child, node)
        # get the supernode buckets containing that node
        supernodes_containing = filter(pair->pair.second==node.second, buckets)
        # give tiebreaker to supernodes that haven't recieved a node yet
        if nodecounts[node.first] % 10 == 0
            nodecounts[node.first] += 1
        end
        # delete the priority queue entry for each supernode containing that node and enqueue a new entry with one less node counted for
        for (sn, _) in supernodes_containing
            nodecounts[sn] -= 10
            delete!(pq, sn)
            if nodecounts[sn] > 1 + samecounts[sn]
                enqueue!(pq, sn=>nodecounts[sn])
            end
        end
        # delete all instances containing the chosen node from the buckets
        filter!(pair->pair.second!=node.second, buckets)
    end
    # restore the child partition to standard Dictionary notation
    childDict = Dict{Integer, Integer}()
    for pair in child
        childDict[pair.second] = pair.first
    end

    return cleanPartition(childDict)
end

#Function: randomWalkMutate
#Parameters: individual, the partition to mutate
#            n, the number of nodes in the partitions
#            k, the number of supernodes in the partitions
#            mutation_prob, the probability of taking another step on the random walk
#Purpose: To mutate a partition
#Return value: returns a partition that is slightly different than the input partition
@everywhere function randomWalkMutate(individual::Dict{Integer,Integer}, n::Integer, k::Integer, mutation_prob::Float64)
    weights = pweights([1-mutation_prob, mutation_prob])
    # each time we flip heads with chance mutation_prob walk to an adjacent partition
    while true
        if StatsBase.sample([0,1], weights) == 1
            individual = randomAdjacentPartition(individual, n, k)
        else
            break
        end
    end
    return individual
end

#Function: geneticImprovement
#Parameters: A, MatrixNetwork to partition and run dynamics on
#            dynamical_function, the dynamic model to run on the network
#            partitions, the partitions we will evolve from, should be an even number
#            generations, the number of reproduction cycles to run
#            mutation_prob, the chance of a partition to randomly changing into a adjacent partition
#            initial_condition, the initial variable values
#            dynamical_function, the dynamical model we will run getLoss using
#            tmax, when to end the model
#            dt, the length of the timesteps in our simulation
#Purpose: To find the partition with local minimum loss using a genetic evolution algorithm
#Return value: returns a set of partitions that evolved from the original set of partitions
@everywhere function geneticImprovement(A::MatrixNetwork, partitions::Array{Dict{Integer, Integer}}, generations::Integer, mutation_prob::Float64, initial_condition::Vector, dynamical_function::Function, tmax::Number, dt::Number; function_args...)
    c = length(partitions)
    n = length(partitions[1])
    k = length(unique(values(partitions[1])))
    originalTimeSeries = simulateODEonGraph(A, initial_condition; dynamical_function=dynamical_function, tmax=tmax, dt=dt, function_args...)

    function_args = Dict(function_args)
    individuals = partitions
    for _=1:generations
        # reproduction phase

        # store the magnitude of loss for each partition
        losses = pmap(individual->getLossPartition(originalTimeSeries, A, individual, initial_condition, dynamical_function, tmax, dt, function_args...), individuals)
        loss_sum = sum(losses)
        prob = Vector{Float64}()
        # each partition reproduces with probability scaled to the magnitude of its loss
        for i = 1:c
            append!(prob, loss_sum/(losses[i]*c))
        end
        weights = aweights(prob)
        # display(weights)
        children = []
        for i = 1:c
            push!(children, individuals[StatsBase.sample(1:c, weights)])
        end
        individuals = children

        # crossing phase
        for i = 1:2:c
            child1 = supernodeBucketCross(individuals[i], individuals[i+1], n, k)
            child2 = supernodeBucketCross(individuals[i+1], individuals[i], n, k)
            individuals[i] = child1
            individuals[i+1] = child2
        end

        # mutation phase
        individuals = pmap(individual->randomWalkMutate(individual, n, k, mutation_prob), individuals)
    end
    return collect(individuals)
end


function countcsvlines(file)
    n = 0
    for row in CSV.Rows(file)
        n += 1
    end
    return n
end

function findLocalMinimum(xyzpData::String, xradius::Number, yradius::Number; startingPartition::String="")
    Random.seed!(trunc(Int, time() * 1000000))
    save = []
    save2 = []
    lossDifference = []
    z1 = 0
    z2 = 0
    minSoFar = ""
    csv_reader = CSV.File(xyzpData)
    if !isempty(startingPartition)
        currentPart = 0
        for row in csv_reader
            currentPart += 1
            if cmp(row.partition, startingPartition) == 0
                break
            end
        end
    else
        numRows = countcsvlines(xyzpData)
        currentPart = rand(2:numRows)
    end
    i, x, y, z = 0, 0, 0, 0
    part = ""
    for row in csv_reader
        i += 1
        if i == currentPart
            x = row.x; y = row.y; z = row.z; part = row.partition; partArray = row.partitionArray
            push!(save, part)
            push!(save2, partArray)
            z1 = z
            z2 = z
            push!(lossDifference, (z1 - z2))
            break
        end
    end
    @label start
    xmin = x - xradius; xmax = x + xradius; ymin = y - yradius; ymax = y + yradius; zcurrent = z
    j = 0
    minIndex = 0
    for row in csv_reader
        j += 1
        if row.x >= xmin && row.x <= xmax && row.y >= ymin && row.y <= ymax
            if row.z < zcurrent
                zcurrent = row.z
                minIndex = j
            end
        end
    end
    k = 0
    for row in csv_reader
        k += 1
        if k == minIndex
            if row.partition != minSoFar
                x = row.x; y = row.y; z = row.z; part = row.partition; partArray = row.partitionArray
                minSoFar = part
                z2 = z
                push!(save, part)
                push!(save2, partArray)
                push!(lossDifference, (z1 - z2))
                @goto start
            else
                z2 = z
                push!(lossDifference, (z1 - z2))
                dt = now()
                DT = Dates.format(dt, "mm-dd_HH-MM-SS")
                df = DataFrame(["partition" => save, "partArray" => save2, "lossDifference" => lossDifference])
                loc = "./data/localMin_data/" * DT * ".csv"
                CSV.write(loc, df)
                return row.partition
                break
            end
        end
    end
    dt = now()
    DT = Dates.format(dt, "mm-dd_HH-MM-SS")
    df = DataFrame(["partition" => save, "partArray" => save2, "lossDifference" => lossDifference])
    loc = "./data/localMin_data/" * DT * ".csv"
    CSV.write(loc, df)
    return part
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
@everywhere function iterativeImprovement(A::MatrixNetwork, p::Dict{Integer, Integer}, depth::Integer, initial_condition::Vector, dynamical_function::Function, tmax::Number, dt::Number; function_args...)
    n = length(p)
    k = length(unique(values(p)))
    function_args = Dict(function_args)
    originalTimeSeries = simulateODEonGraph(A, initial_condition; dynamical_function=dynamical_function, tmax=tmax, dt=dt, function_args...)

    minloss = getLossPartition(originalTimeSeries, A, p, initial_condition, dynamical_function, tmax, dt, function_args...)

    p2 = p
    # iteratively and greedily choose the least adjacent partition until we find a local minimum
    while true
        neighborhood = getNeighborhood(p, n, k, depth)
        # find the adjacent partition with the least loss
        for neighbor in neighborhood
            loss = getLossPartition(originalTimeSeries, A, neighbor, initial_condition, dynamical_function, tmax, dt, function_args...)
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
    end
    return p
end

#Function: getAdjacentPartitions
#Parameters: p, the partition to start from
#            n, the number of nodes
#            k, the number of supernodes
#Purpose: To return a list of the adjacent partitions
#Return value: returns a list of partitions adjacent to the input partition
@everywhere function getAdjacentPartitions(p::Dict{Integer, Integer}, n::Integer, k::Integer)
    adjacents = []
    # apply every possible single node change and record the resulting partitions
    for i=1:n
        for j=1:k
            if p[i] != j
                append!(adjacents, [copy(p)])
                adjacents[end][i] = j
                if length(unique(values(adjacents[end]))) != k
                    pop!(adjacents)
                else
                    adjacents[end] = cleanPartition(adjacents[end])
                end
            end
        end
    end
    return unique(adjacents)
end

#Function: getNeighborhood
#Parameters: p, the partition to start from
#            n, the number of nodes
#            k, the number of supernodes
#            depth, the number of partitions away to look
#Purpose: To return a list of the neighboring partitions
#Return value: returns a list of partitions up to depth away from the input partition
@everywhere function getNeighborhood(p::Dict{Integer, Integer}, n::Integer, k::Integer, depth::Integer)
    # create a vector of set for partitions a certain number of changes away
    neighbors = [Set() for _=1:depth+1]
    # the set of 0 changes away is just p
    union!(neighbors[1], [p])
    # make one change to all partitions in the set of i-1 changes to get the set of partitions with i changes
    for i=1:depth
        for partition in neighbors[i]
            union!(neighbors[i+1], getAdjacentPartitions(partition, n, k))
        end
    end
    # consolidate all layers into one set
    neighbors = union(neighbors...)
    # remove p
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
@everywhere function getAdjacentSample(p::Dict{Integer, Integer}, n::Integer, k::Integer, s::Integer=1)
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
                    sample[end] = cleanPartition(sample[end])
                    break
                end
            end
        end
    end
    return sample
end

#Function: randomAdjacentPartition
#Parameters: p, the partition to start from
#            n, the number of nodes
#            k, the number of supernodes
#Purpose: To return a random adjacent partition
#Return value: returns a partition adjacent to the input partition
@everywhere function randomAdjacentPartition(p::Dict{Integer, Integer}, n::Integer, k::Integer)
    sample = copy(p)
    # change the supernode we partition one node into to create a new partition adjacent to the input partition
    while true
        ns = rand(1:n)
        ks = rand(1:k)
        if p[ns] != ks
            sample[ns] = ks
            if length(unique(values(sample))) != k
                sample[ns] = p[ns]
            else
                break
            end
        end
    end
    return cleanPartition(sample)
end
