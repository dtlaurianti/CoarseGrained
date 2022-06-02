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
    R = zeros((len(A),len(groupSizes)))
    for i in 1:length(A)
        R[i,partition[i]] = 1
        K[partition[i],i] = 1/groupSizes[partition[i]]
    end
    return K⋅A⋅R
end
