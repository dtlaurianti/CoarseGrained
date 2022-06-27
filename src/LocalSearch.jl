#Function: iterativeImprovement
#Parameters: A, MatrixNetwork to partition and run dynamics on
#            dynamical_function, the dynamic model to run on the network
#            p, the partition to start from
#Purpose: To find the partition with local minimum loss using a basic iterative approach
#Return value: returns a partition that is a local minimum of the partition space
function iterativeImprovement(A::MatrixNetwork, dynamical_function::Function, p::Dict{Integer, Integer})
    n = length(p)
    k = length(unique(values(p)))


end

#Function: getNeighbors
#Parameters: p, the partition to start from
#Purpose: To return a list of the neighboring partitions
#Return value: returns a list of partitions adjacent to the input partition
function getNeighbors(p::Dict{Integer, Integer}, n::Integer, k::Integer)
    neighbors = []
    for i=1:n
        for j=1:k
            if p[i] != j
                append!(neighbors, [copy(p)])
                neighbors[end][i] = j
            end
        end
    end
    return neighbors
end

#Function: getNeighbors
#Parameters: p, the partition to start from
#Purpose: To return a list of the neighboring partitions
#Return value: returns a sample list of partitions adjacent to the input partition
function getNeighborSample(p::Dict{Integer, Integer}, n::Integer, k::Integer, s::Integer)
    sample = []
    for _=1:s
        append!(sample, [copy(p)])
        while true
            ns = rand(1:n)
            ks = rand(1:k)
            if p[ns] != ks
                sample[end][ns] = ks
                break
            end
        end
    end
    return sample
end
