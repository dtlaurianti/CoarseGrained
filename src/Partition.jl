using SparseArrays
using Clustering
using ResumableFunctions
using Parameters

# generate numPartitions random partitions mapping a network of size originalSize
# to a network of size reducedSize
function generateRandomPartitions(originalSize::Integer, reducedSize::Integer, numPartitions::Integer)
    partitionList = []
    for index in 1:numPartitions
        partitionAccepted = false
        while !partitionAccepted
            partition = Dict{Integer, Integer}
            # assigns each node in the original nodes to a new node in the reduced nodes
            labels = [1:reducedSize for i=1:originalSize]
            for node in 1:originalSize
                partition[node] = labels[node]
            end
            # check that each partition is nonempty
            if length(Set(labels)) == reducedSize
                partitionAccepted = True
            end
        end
        append!(partitionList, partition)
    end
    return partitionList
end

function spectralClustering(A::MatrixNetwork, reducedSize::Integer)
    # get the eigen decomposition of the diagonal of A
    X = eigen(spdiagm(diag(sparse(A))))
    # compute clusters with k centroids
    labels = kmeans(X.vectors[:,1:reducedSize], reducedSize, init=:kmpp)[1]
    # return the partition
    return Dict{Integer, Integer}[i => labels[i] for i in 1:size(A,1)]
end

function mEEP(A::MatrixNetwork, reducedSize::Integer)
    return print("under dev")
end

# works only for undirected networks as of now.
function agglomerationReduction(A::MatrixNetwork, reducedSize::Integer)
    partition = Dict{Integer, Integer}()
    n = size(A, 1)
    nodeIds = 1:n
    # populate the partition
    for nodeId in nodeIds
        partition[nodeId] = nodeId
    end

    while length(Set(values(partition))) > reducedSize
        partition, Q = greedyMerge(A, partition)
    end
    cleanedPartition = Dict{Integer, Integer}()
    for (node, group) in partition
        try
            partition[node] = cleanedPartition[group]
        catch
            try
                cleanedPartition[group] = max(cleanedPartition.values()) + 1
            catch
                cleanedPartition[group] = 0
            partition[node] = cleanedPartition[group]
            end
        end
    end
    return partition
end

function greedyMerge(A::MatrixNetwork, partition::Dict)
    k = sum(A, dims=1)
    Q = -sum(k.^2)/m^2
    m = sum(k)
    # iterate over unigue groups
    groupIds = collect(Set(values(partition)))
    numGroups = length(groupIds)
    maxQ = -Inf
    for groupIndex1 in 1:numGroups
        for groupIndex2 in 1:groupIndex1
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
                        eUV += sparse(A)[i,j]/m # because directed edge list
                    catch
                    end
                end
            end
            for i in group1
                aU += k[i]/m
            end

            for i in group2
                aV += k[i]/m
            end

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
    nodeIds = [1:n]

    # compute all integer partitions
    for (index, part) in enumerate(partitionNodes(nodeIds))
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

@resumable function partitionNodes(nodeIds::Array{Integer})
    #=
    input:  list of nodeIds
    output: list of all possible partitions
    =#
    if length(nodeIds) == 1
        @yield nodeIds
        return
    end

    first = nodeIds[1]
    for smaller in partitionNodes(nodeIds[2:end])
        # insert first in each of the subpartition's subsets
        for (k, subset) in enumerate(smaller)
            @yield [smaller[1:k] [[ first ] subset] smaller[k+1:end]]
        end
        # put first in its own subset
        @yield [[[first]] smaller]
    end
end

@resumable function kPartition(n::Integer, k::Integer)
    #=
    input: number of nodes n
    output: dictionary of all possible partitions with a specified number of supernodes (k)
    =#
    @resumable function kPartitionNodesAll(nodeIds::Array{Integer}, k::Integer)
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

    @resumable function kPartitionNodes(nodeIds::Array, k::Integer)
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
        append!(kPartitions, partition)                  # store partition in the list of partitions
        # kPartitions[index] = partition                 # store partition in the dictionary of partitions
    end                
    return kPartitions
end
