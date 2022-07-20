using Clustering
using Random
using Distributed
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
        partition = Dict{Integer, Integer}()
        # assigns each node in the original nodes to a new node in the reduced nodes
        labels = StatsBase.sample(1:reducedSize, originalSize)
        # specifically assigns each supernode one random node so they are nonempty
        assignindices = StatsBase.sample(1:originalSize, reducedSize, replace=false)
        for i=1:reducedSize
            labels[assignindices[i]] = i
        end
        partitionList[index] = cleanPartition(Dict{Integer,Integer}(zip(1:originalSize, labels)))
    end
    return partitionList
end

#Function: spectralClustering
#Parameters: A, MatrixNetwork containing the network to be compressed
#            reducedSize, integer number of nodes to map the original nodes to. Should be smaller than originalSize
#Purpose: To reduce a network to the specified number of nodes using the spectral clustering method
#Return value: returns a partition Dictionary specifying the supernodes of the reduced network
function spectralClustering(A::MatrixNetwork, reducedSize::Integer)
    # get the eigen decomposition of the diagonal of A
    X = eigen(Diagonal(diag(sparse(A))))
    # compute clusters with k centroids
    labels = kmeans(X.vectors[:,1:reducedSize]', reducedSize, init=:kmpp).assignments
    # return the partition
    return cleanPartition(Dict{Integer, Integer}(i => labels[i] for i in 1:size(A,1)))
end

function mEEP(A::MatrixNetwork, reducedSize::Integer)
    return print("under dev")
end

#Function: agglomerationReduction
#Parameters: A, MatrixNetwork containing the network to be compressed
#            reducedSize, integer number of nodes to map the original nodes to. Should be smaller than originalSize
#Purpose: To reduce a network to the specified number of nodes using the agglomeration reduction method
#Return value: returns a partition Dictionary specifying the supernodes of the reduced network
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

    # convert the array-format partition into a dictionary-format partition
    aggPartition = Dict{Integer, Integer}()
    for supernode in 1:length(orderedPartition)
        for node in orderedPartition[supernode]
            aggPartition[node] = supernode
        end
    end
    aggPartition = cleanPartition(aggPartition)
    return aggPartition
end

#Function: greedyMerge
#Parameters: A, SparseMatrixCSC containing the network to be compressed
#            orderedPartition, a partition in Array format storing the iterative reductions we've made so far to the network
#            Q, a property of A calculated in the last iteration which we use to measure the change in A
#            k, the degree of the nodes in the original network
#Purpose: To reduce a network by combining two nodes/supernodes into one supernode in a way that minimizes the change in the network,
# agglomerationReduction repeatedly calls this method until the network is reduced to the required size
#Return value: returns a partition Dictionary specifying the supernodes of the reduced network
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

#Function: exhaustivePartition
#Parameters: n, the number of nodes in the resultant partition
#Return value: returns a Dictionary of partition Dictionaries containing all possible partitions of n nodes
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

#= Stirling Numbers of the Second Kind
n\k	0 	1 	2 	3 	4 	5 	6 	7 	8 	9 	10
0 	1
1 	0 	1
2 	0 	1 	1
3 	0 	1 	3 	1
4 	0 	1 	7 	6 	1
5 	0 	1 	15 	25 	10 	1
6 	0 	1 	31 	90 	65 	15 	1
7 	0 	1 	63 	301 	350 	140 	21 	1
8 	0 	1 	127 	966 	1701 	1050 	266 	28 	1
9 	0 	1 	255 	3025 	7770 	6951 	2646 	462 	36 	1
10 	0 	1 	511 	9330 	34105 	42525 	22827 	5880 	750 	45 	1
=#

#Function: kPartition
#Parameters: n, the number of nodes in the original network
#            k, the number of supernodes in the compressed network
#Purpose: To reduce an arbitrary network of n nodes to k nodes
#Return value: returns a Dictionary of Partition Dictionaries containing all possible reductions from n to k nodes
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
            else
                # otherwise move on to the node to the left
                ind -= 1
            end
        end
        # if we have created exactly k supernodes accept the partition
        if maximum(values(curPartition)) == k
            c += 1
            allPartitions[c] = copy(curPartition)
        end
    end
    return allPartitions
end

#Function: cleanPartition
#Parameters: partition, the partition to be cleaned
#Purpose: To convert a partition into standard form
#Return value: Returns a partition in standard format
#=
@everywhere function cleanPartition(partition::Dict{Integer,Integer})
    cleanedPartition = Dict()
    for (node, supernode) in partition
        try # if the supernode has already been reassigned use that
            partition[node] = cleanedPartition[supernode]
        catch # if it has not been reassign it
            try # to the next value availible
                cleanedPartition[supernode] = maximum(values(cleanedPartition)) + 1
            catch # or 1 if it is the first supernode we've cleaned
                cleanedPartition[supernode] = 1
            end # then assign that node to the new label
            partition[node] = cleanedPartition[supernode]
        end
    end
    return partition
end
=#
@everywhere function cleanPartition(partition::Dict{Integer,Integer})
    cleanedPartition = Dict()
    maxv = 1
    for (node, supernode) in sort(collect(partition))
        if haskey(cleanedPartition, supernode)
            partition[node] = cleanedPartition[supernode]
        else
            cleanedPartition[supernode] = maxv
            maxv += 1
            partition[node] = cleanedPartition[supernode]
        end
    end
    return partition
end
