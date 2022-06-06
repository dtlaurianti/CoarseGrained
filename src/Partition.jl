# generate numPartitions random partitions mapping a network of size originalSize
# to a network of size reducedSize
function generateRandomPartitions(originalSize, reducedSize, numPartitions::Integer):
    partitionList = []
    for index in 1:numPartitions:
        partitionAccepted = false
        while !partitionAccepted:
            partition = Dict{Integer, Integer}
            # assigns each node in the original nodes to a new node in the reduced nodes
            labels = [1:reducedSize for i=1:originalSize]
            for node in 1:originalSize:
                partition[node] = labels[node]
            end
            # check that each partition is nonempty
            if size(Set(labels),1) == reducedSize:
                partitionAccepted = True
            end
        end
        append!(partitionList, partition)
    end
    return partitionList
end

function spectralClustering(A::MatrixNetwork, reducedSize::Integer):
    l, X = np.linalg.eig(np.diag(A))
    labels = k_means(X[:,:reducedSize], reducedSize, init='k-means++')[1]
    return {i:labels[i] for i in range(len(A))}

def mEEP(A, reducedSize):
    return print("under dev")

# works only for undirected networks as of now.
def agglomerationReduction(A, reducedSize):
    partition = dict()
    n = np.size(A, axis=0)
    nodeIds = range(n)
    for nodeId in nodeIds:
        partition[nodeId] = nodeId

    k = np.sum(A, axis=0)
    m = np.sum(k)
    Q = -np.sum(np.power(k, 2))/m**2
    while len(set(partition.values())) > reducedSize:
        partition, Q = greedyMerge(A, Q, partition, k)
        print(Q)
    cleanedPartition = dict()
    for node, group in partition.items():
        try:
            partition[node] = cleanedPartition[group]
        except:
            try:
                cleanedPartition[group] = max(cleanedPartition.values()) + 1
            except:
                cleanedPartition[group] = 0
            partition[node] = cleanedPartition[group]
    return partition

def greedyMerge(A, Q, partition, k):
    # iterate over unigue groups
    m = np.sum(k)
    groupIds = list(set(partition.values()))
    numGroups = len(groupIds)
    maxQ = -math.inf
    for groupIndex1 in range(numGroups):
        for groupIndex2 in range(groupIndex1):
            groupId1 = groupIds[groupIndex1]
            groupId2 = groupIds[groupIndex2]

            group1 = [nodeId for nodeId, group in partition.items() if group == groupId1]
            group2 = [nodeId for nodeId, group in partition.items() if group == groupId2]

            eUV = 0
            aU = 0
            aV = 0
            for i in group1:
                for j in group2:
                    try:
                        eUV += A[i,j]/m # because directed edge list
                    except:
                        pass
            for i in group1:
                aU += k[i]/m

            for i in group2:
                aV += k[i]/m

            newQ = Q + 2*(eUV - aU*aV)

            if newQ > maxQ:
                maxQ = newQ
                maxGroupId1 = groupId1
                maxGroupId2 = groupId2
    newPartition = dict(partition)
    for node in list(newPartition.keys()):
        if newPartition[node] == maxGroupId2:
            newPartition[node] = maxGroupId1
    return newPartition, maxQ

# WARNING: not super efficient at the moment!
def exhaustivePartition(n):
    '''
    input: number of nodes n
    output: dictionary of all possible partitions
    '''
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

def partitionNodes(nodeIds):
    '''
    input:  list of nodeIds
    output: list of all possible partitions
    '''
    if len(nodeIds) == 1:
        yield [ nodeIds ]
        return

    first = nodeIds[0]
    for smaller in partitionNodes(nodeIds[1:]):
        # insert first in each of the subpartition's subsets
        for k, subset in enumerate(smaller):
            yield smaller[:k] + [[ first ] + subset] + smaller[k+1:]
        # put first in its own subset
        yield [[ first ]] + smaller

def kPartition(n, k):
    '''
    input: number of nodes n
    output: dictionary of all possible partitions with a specified number of supernodes (k)
    '''
    def kPartitionNodesAll(nodeIds, k):
        '''
        generate all possible partitions of nodeIds into k supernodes
        includes empty supernodes that need to be removed
        '''
        n = len(nodeIds)
        if k ==1:
            yield [nodeIds]
        elif n == k:
            yield [[node] for node in nodeIds]
        else:
            first, *rest = nodeIds
            for subset in kPartitionNodesAll(rest, k-1):
                yield [[first], *subset]
            for subset in kPartitionNodesAll(rest, k):
                for i in range(len(subset)):
                    yield [[first] + subset[i]] + subset[:i] + subset[i+1:]

    def kPartitionNodes(nodeIds, k):
        '''
        generate all possible partitions of nodeIds into k supernodes
        using kPartitionNodesAll, then remove empy supernodes
        '''
        # remove empty supernodes
        for part in kPartitionNodesAll(nodeIds, k):
            if all(supernode for supernode in part):      # check if all supernodes are populated
                yield part

    # initialize
    kPartitions = []
    nodeIds = list(range(n))
    # list all ways to split n nodes into k supernodes, then format as a partition
    for index, part in enumerate(kPartitionNodes(nodeIds, k),1):
        partition = dict() # initialize partition
        for supernodeId, subnodeIds in enumerate(part): # for each list of partitioned nodes
            for nodeId in subnodeIds:                   # for each node
                partition[nodeId] = supernodeId         # assign id of the supernode it belongs to
        kPartitions.append(partition)                   # store partition in the list of partitions
        # kPartitions[index] = partition                # store partition in the dictionary of partitions
    return kPartitions
