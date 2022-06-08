using MatrixNetworks
using LinearAlgebra
include("Partition.jl")

# counts the number of nodes in each supernode of a partition
function getSupernodeSizes(partition::Dict{Integer,Integer})
  supernodeSizes = Dict{Integer,Integer}()
  for supernode in values(partition)
      supernodeSizes[supernode] = get!(supernodeSizes, supernode, 0)+1
    end
  return supernodeSizes
end

# Reduce Initial Conditions - REQUIRES SEQUENTIAL (no missing) SUPERNODE LABELS
function compressInitialCondition(initialConditions::Vector, partition::Dict{Integer,Integer})
    supernodeSizes = getSupernodeSizes(partition)
    compressedInitialConditions = zeros(length(supernodeSizes))
    for i in 1:length(initialConditions)
        # scales the initial variables of each node by the number of nodes incorporated into the supernode
        # and then sets the initial variables of the supernode as a combination of all nodes in its partition
        compressedInitialConditions[partition[i]] += initialConditions[i]/supernodeSizes[partition[i]]
    end
    return compressedInitialConditions
end

# Reduce Networks with a partition
function compressAdjacencyMatrix(A::MatrixNetwork, partition::Dict{Integer,Integer})
    A = sparse(A)
    #Reduce matrix using spectral method from Gfeller et al. (2008).
    numGroups = length(Set(values(partition)))
    groupSizes = getSupernodeSizes(partition)
    K = zeros((length(groupSizes),size(A,1)))
    R = zeros((size(A,1),length(groupSizes)))
    for i in 1:size(A,1)
        R[i,partition[i]] = 1
        K[partition[i],i] = 1/groupSizes[partition[i]]
    end
    return MatrixNetwork(sparse(K*A*R))
end
