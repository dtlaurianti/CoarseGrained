using MatrixNetworks
using LinearAlgebra

function getSupernodeSizes(partition)
  supernodeSizes = dict()
  for supernode in partition.values()
      try
          supernodeSizes[supernode] += 1
      catch
          supernodeSizes[supernode] = 1
      end
    end
  return supernodeSizes
end

# Reduce Initial Conditions - REQUIRES SEQUENTIAL (no missing) SUPERNODE LABELS
function compressInitialCondition(initialConditions, partition)
    supernodeSizes = getSupernodeSizes(partition)
    compressedInitialConditions = zeros(length(supernodeSizes))
    for i in 1:length(initialConditions)
        compressedInitialConditions[partition[i]] += initialConditions[i]./supernodeSizes[partition[i]]
    end
    return compressedInitialConditions
end

# Reduce Networks with a partition
function compressAdjacencyMatrix(A, partition)
    #Reduce matrix using spectral method from Gfeller et al. (2008).
    numGroups = length(set(partition.values()))
    groupSizes = getSupernodeSizes(partition)
    K = zeros((length(groupSizes),length(A)))
    R = zeros((length(A),length(groupSizes)))
    for i in 1:length(A)
        R[i,partition[i]] = 1
        K[partition[i],i] = 1/groupSizes[partition[i]]
    end
    return K⋅A⋅R
end

######## Partitioning ########
function generateRandomPartitions(originalSize, reducedSize, numPartitions)
    partitionList = []
    for index in 1:numPartitions
        partitionAccepted = false
        while !partitionAccepted
            partition = dict()
            #labels = random.choices(range(reducedSize), k=originalSize)
            labels = []
            for i in 1:originalSize
                push!(labels, sample(1:reducedSize))
            end
            for node in 1:originalSize
                partition[node] = labels[node]
            end
            if length(set(labels)) == reducedSize
                partitionAccepted = true
            end
        end
        push!(partitionList, partition)
    end 
    return partitionList
end

function spectralClustering(A, reducedSize)
    l, X = eigvals(diag(A))
    #Figure out how to do this
    labels = k_means(X[:,:reducedSize], reducedSize, init='k-means++')[1]
    return {i:labels[i] for i in range(len(A))}
end

function mEEP(A, reducedSize)
    return print("under dev")
end

# works only for undirected networks as of now.
function agglomerationReduction(A, reducedSize)
    partition = dict()
    n = size(A, 1)
    for nodeId in 1:n
        partition[nodeId] = nodeId
    end
    k = sum(A, 1)
    m = sum(k)
    Q = -sum(k.^2)./m.^2
    while length(set(values(partition))) > reducedSize
        partition, Q = greedyMerge(A, Q, partition, k)
        println(Q)
    end
    cleanedPartition = dict()
    items = []
    for i in partition
        push!(items, i)
    end
    for node, group in items
        try
            partition[node] = cleanedPartition[group]
        catch
            try
                cleanedPartition[group] = max(values(cleanedPartition)) + 1
            catch
                cleanedPartition[group] = 0
            end
            partition[node] = cleanedPartition[group]
        end
    end
    return partition
 end

#not finished
function greedyMerge(A, Q, partition, k)
    # iterate over unique groups
    m = sum(k)
    groupIds = list(set(values(partition)))
    numGroups = length(groupIds)
    maxQ = -Inf
    for groupIndex1 in 1:numGroups
        for groupIndex2 in 1:groupIndex1
            groupId1 = groupIds[groupIndex1]
            groupId2 = groupIds[groupIndex2]

            group1 = [nodeId for nodeId, group in partition.items() if group == groupId1]
            group2 = [nodeId for nodeId, group in partition.items() if group == groupId2]

            eUV = 0
            aU = 0
            aV = 0
            for i in group1:
                for j in group2:
                    try
                        eUV += A[i,j]/m # because directed edge list
                    catch
                        pass
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
    newPartition = dict(partition)
    for node in list(newPartition.keys()):
        if newPartition[node] == maxGroupId2
            newPartition[node] = maxGroupId1
        end
    end
    return newPartition, maxQ
end

#not finished
# WARNING: not super efficient at the moment!
function exhaustivePartition(n)
    #=
    input: number of nodes n
    output: dictionary of all possible partitions
    =#
    allPartitions = {}
    nodeIds = list(range(n))

    # compute all integer partitions
    for index, part in enumerate(partitionNodes(nodeIds),1):
        # print(index, part)

        partition = dict()
        for supernodeId, subnodeIds in enumerate(part):
            for nodeId in subnodeIds:
                partition[nodeId] = supernodeId

        allPartitions[index] = partition

    return allPartitions
end

function partitionNodes(nodeIds)
    #=
    input:  list of nodeIds
    output: list of all possible partitions
    =#
    if len(nodeIds) == 1:
        #TODO: determine how to get python-like yield functionality in Julia 
        yield [ nodeIds ]
        return
    end
    first = nodeIds[0]
    for smaller in partitionNodes(nodeIds[1:]):
        # insert first in each of the subpartition's subsets
        for k, subset in enumerate(smaller)
            yield smaller[:k] + [[ first ] + subset] + smaller[k+1:]
        # put first in its own subset
        yield [[ first ]] + smaller
end