#Function: getSupernodeSizes
#Parameters: partition, the partition that specifies the supernodes
#Purpose: Counts the number of nodes in each supernode of a partition
#Return value: A Vector of the number of nodes in each supernode
@everywhere function getSupernodeSizes(partition::Dict{Integer,Integer})
  supernodeSizes = Dict{Integer,Integer}()
  for supernode in values(partition)
      supernodeSizes[supernode] = get!(supernodeSizes, supernode, 0)+1
    end
  return supernodeSizes
end

#Function: compressInitialCondition
#Parameters: initialConditions, a Vector of the inital value of the node variables
#            partition, a Dict that specifies the new coarse-grained network
#Purpose: creates a new initial value for each new supernode based on the initial values of the nodes it contains
#Return value: A Vector of the new initial conditions with length equal to the number of supernodes
# Reduce Initial Conditions - REQUIRES SEQUENTIAL (no missing) SUPERNODE LABELS
@everywhere function compressInitialCondition(initialConditions::Vector, partition::Dict{Integer,Integer})
    supernodeSizes = getSupernodeSizes(partition)
    compressedInitialConditions = zeros(length(supernodeSizes))
    for i in 1:length(initialConditions)
        # scales the initial variables of each node by the number of nodes incorporated into the supernode
        # and then sets the initial variables of the supernode as a combination of all nodes in its partition
        compressedInitialConditions[partition[i]] += initialConditions[i]/supernodeSizes[partition[i]]
    end
    return compressedInitialConditions
end

#Function: compressArguments
#Parameters: partition, a Dict that specifies the new coarse-grained network
#            function_args, the var-kwargs that was passed into getLoss
#Purpose: creates a new initial set of arguments for the dynamical model with appropriately sized Vectors for the partitioned network
#Return value: A Dict of the modified function_args
@everywhere function compressArguments(partition::Dict{Integer,Integer}, function_args...)
    function_args = Dict(function_args)
    # If function_args is specifying ω which is a vector argument we need to compress it
    if haskey(function_args, :ω)
        supernodeSizes = getSupernodeSizes(partition)
        ω = [0.0 for _=1:length(supernodeSizes)]
        # compute the ω for supernodes as the average of their contained nodes' ω
        for i in 1:length(function_args[:ω])
            ω[partition[i]] +=  function_args[:ω][i]/supernodeSizes[partition[i]]
        end
        function_args[:ω] = ω
    end
    return function_args
end


#Function: compressAdjacencyMatrix
#Parameters: A, a MatrixNetwork that we wish to compress according to the partition
#            partition, the partition that specifies the supernodes in the created network
#Purpose: To compress a network according to a partition
#Return value: A new, compressed MatrixNetwork
@everywhere function compressAdjacencyMatrix(A::MatrixNetwork, partition::Dict{Integer,Integer})
    A = sparse(A)
    #Reduce matrix using spectral method from Gfeller et al. (2008).
    numGroups = length(unique(values(partition)))
    groupSizes = getSupernodeSizes(partition)
    K = zeros((length(groupSizes),size(A,1)))
    R = zeros((size(A,1),length(groupSizes)))
    for i in 1:size(A,1)
        R[i,partition[i]] = 1
        K[partition[i],i] = 1/groupSizes[partition[i]]
    end
    return MatrixNetwork(sparse(K*A*R))
end

#Function: compressNodeCoordinates
#Parameters: x, y, coordinates
#            partition, the partition that specifies the supernodes in the created network
#Purpose: To compress the coordinates of nodes according to a partition
#Return value: A new, compressed coordinate system
@everywhere function compressNodeCoordinates(x::Vector,y::Vector, partition::Dict{Integer,Integer})
    supernodeSizes = getSupernodeSizes(partition)
    n = length(partition)
    c = length(supernodeSizes)
    cx = Vector{Float64}()
    cy = Vector{Float64}()
    # each supernodes new x, y values are just an average of its nodes
    for supernode=1:c
        append!(cx, mean([x[n] for (n,sn) in partition if sn==supernode]))
        append!(cy, mean([y[n] for (n,sn) in partition if sn==supernode]))
    end
    return cx, cy
end
