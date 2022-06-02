using MatrixNetworks
using LinearAlgebra

function getSupernodeSizes(partition)
  supernodeSizes = dict()
  for supernode in list(partition.values())
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
    compressedInitialConditions = np.zeros(len(supernodeSizes))
    for i in range(len(initialConditions))
        compressedInitialConditions[partition[i]] += initialConditions[i]/supernodeSizes[partition[i]]
    end
    return compressedInitialConditions
end
