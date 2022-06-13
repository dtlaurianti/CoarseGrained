using SparseArrays
using Clustering
using ResumableFunctions
using Random
using StatsBase

# generate numPartitions random partitions mapping a network of size originalSize
# to a network of size reducedSize
function generateRandomPartitions(originalSize::Integer, reducedSize::Integer, numPartitions::Integer)
    partitionList = Array{Dict{Integer,Integer}, 1}(undef, numPartitions)
    Random.seed!(trunc(Int, time() * 1000000))
    for index in 1:numPartitions
        partitionAccepted = false
        while !partitionAccepted
            partition = Dict{Integer, Integer}()
            # assigns each node in the original nodes to a new node in the reduced nodes
            labels = StatsBase.sample(1:reducedSize, originalSize)
            for node in 1:originalSize
                partition[node] = labels[node]
            end
            # check that each partition is nonempty
            if length(Set(labels)) == reducedSize
                partitionAccepted = true
                partitionList[index]=partition
            end
        end
    end
    return partitionList
end

function spectralClustering(A::MatrixNetwork, reducedSize::Integer)
    # get the eigen decomposition of the diagonal of A
    X = eigen(Diagonal(diag(sparse(A))))
    # compute clusters with k centroids
    labels = kmeans(X.vectors[:,1:reducedSize]', reducedSize, init=:kmpp).assignments
    # return the partition
    return Dict{Integer, Integer}(i => labels[i] for i in 1:size(A,1))
end

function mEEP(A::MatrixNetwork, reducedSize::Integer)
    return print("under dev")
end

# works only for undirected networks as of now.
function agglomerationReduction(A::MatrixNetwork, reducedSize::Integer)
    A = sparse(A)
    partition = Dict{Integer, Integer}()
    n = size(A, 1)
    nodeIds = 1:n
    # populate the partition
    for nodeId in nodeIds
        partition[nodeId] = nodeId
    end

    k = sum(A, dims=1)
    m = sum(k)
    Q = -sum(k.^2)/m^2

    while length(Set(values(partition))) > reducedSize
        partition, Q = greedyMerge(A, partition, Q, k)
    end

    # make the supernode IDs number 1:n
    cleanedPartition = Dict{Integer, Integer}()
    for (node, group) in partition
        try
            partition[node] = cleanedPartition[group]
        catch
            try
                cleanedPartition[group] = maximum(collect(values(cleanedPartition))) + 1
            catch
                cleanedPartition[group] = 1
            end
            partition[node] = cleanedPartition[group]
        end
    end
    return partition
end

function greedyMerge(A::SparseMatrixCSC, partition::Dict, Q::Number, k::Matrix)
    m = sum(k)
    # iterate over unigue groups
    groupIds = collect(Set(values(partition)))
    numGroups = length(groupIds)
    maxQ = -Inf
    # loop over all supernodes and find the two most closely connected ?
    maxGroupId1 = 1
    maxGroupId2 = 1
    for groupIndex1 in 1:numGroups
        for groupIndex2 in 1:(groupIndex1-1)
            groupId1 = groupIds[groupIndex1]
            groupId2 = groupIds[groupIndex2]

            group1 = [nodeId for (nodeId, group) in partition if group == groupId1]
            group2 = [nodeId for (nodeId, group) in partition if group == groupId2]

            eUV = 0
            aU = 0
            aV = 0
            for i in group1
                for j in group2
                    try
                        # the edge between the two supernodes
                        eUV += A[i,j]/m # because directed edge list
                    catch
                    end
                end
            end
            # calculate the magnitude of supernode U
            for i in group1
                aU += k[i]/m
            end

            # calculate the magnitude of supernode V
            for i in group2
                aV += k[i]/m
            end

            # calculate the change in the graph
            newQ = Q + 2*(eUV - aU*aV)
            if newQ > maxQ
                maxQ = newQ
                maxGroupId1 = groupId1
                maxGroupId2 = groupId2
            end
        end
    end
    newPartition = copy(partition)
    for node in collect(keys(newPartition))
        if newPartition[node] == maxGroupId2
            newPartition[node] = maxGroupId1
        end
    end
    return newPartition, maxQ
end

# WARNING: not super efficient at the moment!
function exhaustivePartition(n::Integer)
    #=
    input: number of nodes n
    output: dictionary of all possible partitions
    =#
    allPartitions = Dict{Integer,Dict{Integer,Integer}}()
    nodeIds = 1:n

    # compute all integer partitions
    println("start")
    display(partitionNodes(nodeIds))
    println("helloout")
    for (index, part) in enumerate(partitionNodes(nodeIds))
        println("iter")
        partition = Dict{Integer,Integer}()
        for (supernodeId, subnodeIds) in enumerate(part)
            for nodeId in subnodeIds
                partition[nodeId] = supernodeId
            end
        end
        allPartitions[index] = partition
    end
    return allPartitions
end

@resumable function partitionNodes(nodeIds)
    println("call")
    #=
    input:  list of nodeIds
    output: list of all possible partitions
    =#
    if length(nodeIds) == 1
        println("exit")
        @yield nodeIds
        return
    end

    first = nodeIds[1]
    for smaller in partitionNodes(nodeIds[2:end])
        # insert first in each of the subpartition's subsets
        for (k, subset) in enumerate(smaller)
            println("hello")
            display([smaller[1:k] [[ first ] subset] smaller[k+1:end]])
            @yield [smaller[1:k] [[ first ] subset] smaller[k+1:end]]
        end
        # put first in its own subset
        display([[[first]] smaller])
        @yield [[[first]] smaller]
    end
    println("done")
end

@resumable function kPartitionNodesAll(nodeIds, k::Integer)
    #=
    generate all possible partitions of nodeIds into k supernodes
    includes empty supernodes that need to be removed
    =#
    n = length(nodeIds)
    if k==1
        @yield [nodeIds]
    elseif n == k
        @yield [[node] for node in nodeIds]
    else
        first = nodeIds[1]
        rest = nodeIds[2:end]
        for subset in kPartitionNodesAll(rest, k-1)
            @yield [[first] subset]
        end
        for subset in kPartitionNodesAll(rest, k)
            for i in 1:length(subset)
                @yield [[[first] subset[i]] subset[1:i] subset[i+1:end]]
            end
        end
    end
end

@resumable function kPartitionNodes(nodeIds, k::Integer)
    #=
    generate all possible partitions of nodeIds into k supernodes
    using kPartitionNodesAll, then remove empy supernodes
    =#
    # remove empty supernodes
    for part in kPartitionNodesAll(nodeIds, k)
        if all(supernode for supernode in part)    # check if all supernodes are populated
            @yield part
        end
    end
end

function kPartition(n::Integer, k::Integer)
    #=
    input: number of nodes n
    output: dictionary of all possible partitions with a specified number of supernodes (k)
    =#

    # initialize
    kPartitions = []
    nodeIds = [1:n]
    # list all ways to split n nodes into k supernodes, then format as a partition
    for (index, part) in enumerate(kPartitionNodes(nodeIds, k))
        partition = Dict{Integer,Integer}() # initialize partition
        for (supernodeId, subnodeIds) in enumerate(part) # for each list of partitioned nodes
            for nodeId in subnodeIds                     # for each node
                partition[nodeId] = supernodeId          # assign id of the supernode it belongs to
            end
        end
        append!(kPartitions, [partition])                  # store partition in the list of partitions
        # kPartitions[index] = partition                 # store partition in the dictionary of partitions
    end
    return kPartitions
end
