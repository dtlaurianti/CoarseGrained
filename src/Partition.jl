using SparseArrays
using Clustering
using Random
using StatsBase
using Distributed
using SharedArrays
using Base.Threads

# generate numPartitions random partitions mapping a network of size originalSize
# to a network of size reducedSize

#Function: generateRandomPartitions
#Parameters: originalSize, integer number of nodes in original network
#            reducedSize, integer number of nodes to map the original nodes to. Should be smaller than originalSize
#            numPartitions, integer number of different partitions to make from originalSize nodes to reducedSize nodes
#Purpose: To generate numPartitions many mappings from originalSize nodes to reducedSize nodes.
#Return value: An array of dictionaries where there are originalSize many keys associated with reducedSize
#              many values such that multiple keys will be associated with the same value(s) representing
#              multiple nodes in the original network being mapped to the same "supernode" in the reduced network.
#              There are numPartitions many dictionaries in the array, each representing a different mapping of
#              originalSize many nodes to reducedSize many nodes.
function generateRandomPartitions(originalSize::Integer, reducedSize::Integer, numPartitions::Integer)
    partitionList = Array{Dict{Integer,Integer}, 1}(undef, numPartitions)
    Random.seed!(trunc(Int, time() * 1000000))
    @threads for index in 1:numPartitions
        partitionAccepted = false
        while !partitionAccepted
            partition = Dict{Integer, Integer}()
            # assigns each node in the original nodes to a new node in the reduced nodes
            labels = StatsBase.sample(1:reducedSize, originalSize)
            for node in 1:originalSize
                partition[node] = labels[node]
            end
            # check that each partition is nonempty
            if length(unique(labels)) == reducedSize
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
    # the supernodes of the partition
    groupIds = collect(unique(values(partition)))
    numGroups = length(groupIds)
    orderedPartition = [[] for i=1:numGroups]
    for (nodeId, groupId) in partition
        append!(orderedPartition[groupId], nodeId)
    end
    # reduce the number of supernodes by one per each iteration of greedyMerge
    while length(orderedPartition) > reducedSize
        orderedPartitionpartition, Q = greedyMerge(A, orderedPartition, Q, k)
    end

    aggPartition = Dict{Integer, Integer}()
    for supernode in 1:length(orderedPartition)
        for node in orderedPartition[supernode]
            aggPartition[node] = supernode
        end
    end
    return aggPartition
end

function greedyMerge(A::SparseMatrixCSC, orderedPartition::Vector{Vector{Any}}, Q::Number, k::Matrix)
    m = sum(k)
    numGroups = length(orderedPartition)
    lk = ReentrantLock()
    maxQ = Dict(i=>-Inf for i=1:nthreads())
    # loop over all supernodes and find the two most closely connected ?
    maxGroupId1 = Dict(i=>1 for i=1:nthreads())
    maxGroupId2 = Dict(i=>1 for i=1:nthreads())

    @threads for groupIndex1 in 1:numGroups
        for groupIndex2 in 1:(groupIndex1-1)
            eUV = 0
            aU = 0
            aV = 0
            for i in orderedPartition[groupIndex1]
                for j in orderedPartition[groupIndex2]
                    # the edge between the two supernodes
                    eUV += A[i,j]/m # because directed edge list
                end
            end
            # calculate the magnitude of supernode U
            for i in orderedPartition[groupIndex1]
                aU += k[i]/m
            end

            # calculate the magnitude of supernode V
            for i in orderedPartition[groupIndex2]
                aV += k[i]/m
            end

            # calculate the change in the graph
            newQ = Q + 2*(eUV - aU*aV)
            id = threadid()
            if newQ > maxQ[id]
                maxQ[id] = newQ
                maxGroupId1[id] = groupIndex1
                maxGroupId2[id] = groupIndex2
            end
        end
    end
    maxQInd = maximum(maxQ)[1]
    # merge the two nodes into one supernode
    node = orderedPartition[maxGroupId1[maxQInd]]
    deleteat!(orderedPartition, maxGroupId1[maxQInd])
    append!(orderedPartition[maxGroupId2[maxQInd]], node)
    return orderedPartition, maxQ[maxQInd]
end

# Takes the number of nodes n and returns a Dictionary of all possible partitions
# The Dict should be of length Σ S(n)  (Stirling Number of the Second Kind)
function exhaustivePartition(n::Integer)
    allPartitions = Dict{Integer,Dict{Integer, Integer}}()
    c = 1
    # storing the partition in a way that lets us find the next partition in O(1) time
    curPartition = Dict{Integer, Integer}(i=>1 for i=1:n)
    allPartitions[c] = copy(curPartition)
    # storing the max supernode a node can be incremented to
    m = [0 for _=1:n]
    while curPartition[n] < n
        for i = 2:n
            m[i] = max(m[i-1], curPartition[i-1])
        end
        iter = false
        ind = n
        # get the rightmost node that we can increment
        while !iter
            # if this node can be moved to the next supernode or a new supernode can be created
            if curPartition[ind] <= m[ind]
                iter = true
                # move the node into the next supernode
                curPartition[ind] += 1
                # move all nodes after this node back to the first supernode
                for i=(ind+1):n
                    curPartition[i] = 1
                end
            else
                ind -= 1
            end
        end
        c += 1

        # this line is the slow part of this function, but no speedup was to be had by simply working inside allPartitions from the start
        # not sure if there is a better way to do this
        allPartitions[c] = copy(curPartition)
    end
    return allPartitions
end

# Gives a Dict of all partitions of n nodes with k supernodes
# the number of partitions should be equal to S(n, k) (Stirling Number of the Second Kind)
function kPartition(n::Integer, k::Integer)
    allPartitions = Dict{Integer,Dict{Integer, Integer}}()
    c = 0
    # storing the partition in a way that lets us find the next partition in O(1) time
    curPartition = Dict{Integer, Integer}(i=>1 for i=1:n)
    # storing the max supernode a node can be incremented to
    m = [0 for _=1:n]
    while curPartition[n] < n
        for i = 2:n
            m[i] = max(m[i-1], curPartition[i-1])
        end
        iter = false
        ind = n
        # get the rightmost node that we can increment
        while !iter
            # if this node can be moved to the next supernode or a new supernode can be created
            if curPartition[ind] <= m[ind]
                iter = true
                # move the node into the next supernode
                curPartition[ind] += 1
                # move all nodes after this node back to the first supernode
                for i=(ind+1):n
                    curPartition[i] = 1
                end
                #curPartition[ind+1:end] .= 1
            else
                ind -= 1
            end
        end
        if maximum(values(curPartition)) == k
            c += 1
            allPartitions[c] = copy(curPartition)
        end
    end
    return allPartitions
end
