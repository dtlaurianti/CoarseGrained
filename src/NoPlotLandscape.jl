using Clustering
using MultivariateStats
using DataFrames
using CSV
using ScikitLearn
using StatsBase

#Function: dict_to_array
#Parameters: partitions, an array of dictionaries each representing a partition from some number of nodes to some
#            smaller number of nodes.
#Purpose: to create an array of arrays each containing smaller arrays representing nodes in the reduced network.
#         These smaller arrays contain all of the original nodes that were mapped to that node in the reduced network.
#Return value: An array of arrays representing the reduced network for each partition in the original partitions array.
function dict_to_array(partitions::Array{Dict{Integer,Integer}})
    Arr = []
    for i in 1:length(partitions)
        partition = partition_dict_to_array(partitions[i])
        append!(Arr, [partition])
    end
    return Arr
end

#Function: partition_dict_to_array
#Parameters: partition, a dictionary representing a partition from some number of nodes to some
#            smaller number of nodes.
#Purpose: to create an array of smaller arrays representing nodes in the reduced network.
#         These smaller arrays contain all of the original nodes that were mapped to that node in the reduced network.
#         This function performs the same task as dict_to_array, but for just one partition dictionary as opposed to
#         an array of them.
#Return value: An array representing the nodes of the reduced network defined in the partition dictionary
function partition_dict_to_array(partition::Dict{Integer,Integer})
    # the quantity of nodes in each supernode
    supernodeSizes = getSupernodeSizes(partition)
    # create an array of supernodes where each supernode is represented as an array of the nodes in its partition
    partitionArray = [[] for i in 1:length(supernodeSizes)]
    # add each node to its supernodes array
    for (node, supernode) in partition
        push!(partitionArray[supernode], node)
    end
    return partitionArray
end

#Function variation_of_information
#
#Input: Two partitions in nested array format
#Output: The distance between the two partitions
#
# Meila, M. (2007). Comparing clusterings-an information
#   based distance. Journal of Multivariate Analysis, 98,
#   873-895. doi:10.1016/j.jmva.2006.11.013
# https://gist.github.com/jwcarr/626cbc80e0006b526688

#Function: variation_of_information
#Parameters: X, an array of arrays of integers
#            Y, an array of arrays of integers
#Purpose: To provide a measure of the difference between two partitions
#Return value: a positive, floating point value that is higher the more different
#              X and Y are.
@everywhere function variation_of_information(X::Vector{Vector{Any}},Y::Vector{Vector{Any}})
  n = float(sum([size(x,1) for x in X]))
  σ = 0.0
  for x in X
    # p = the ratio of nodes in the supernode x to nodes in the graph
    p = size(x,1) / n
    for y in Y
      # q = the ratio of nodes in the supernode y to nodes in the graph
      q = size(y,1) / n
      # r = the ratio of nodes in both x & y to nodes in the graph
      r = (length(unique(vcat(x,y))) / n)
      # if x & y share at least one node they are seen as comparable supernodes
      if r > 0.0
        # add to the distance between partitions
        σ += r * (log(2, r / p) + log(2, r / q))
      end
    end
  end
  return abs(σ)
end


function GetXYZ(partitions::Vector{Dict{Integer, Integer}}, A, NumOriginalNodes; save_to_string="", modelType::Function=linear_model)
    listModelArgs = Dict(:ϵ=>-3/NumOriginalNodes, :β=>0.5, :γ=>0.5, :ω=>rand(NumOriginalNodes), :K=>0.5, :d=>0.5, :c=>0.5, :b=>0.5)
    #convert dictionary to an array
    Arr = dict_to_array(partitions)

    #calculate distance matrix
    num_par = length(partitions)
    D = SharedArray{Float64}((num_par,num_par))
    @sync @distributed for i in 1:num_par
        for j in 1:num_par
            # the dissimilarity matrix
            D[i, j] = variation_of_information(Arr[i],Arr[j])
        end
    end

    D = Array(D)
    #calculate MDS on disimilarity matrix
    embedding = StatsBase.fit(MDS, D, distances=true, maxoutdim=2)
    X_transformed = StatsBase.predict(embedding)
    #Format data
    x = X_transformed[1,:]
    y = X_transformed[2,:]

    #Calculate z dimension
    z = SharedArray{Float64}((num_par))
    @sync @distributed for i in 1:num_par
        # using hard-coded model and parameters, possibly want to make the outer function accept those parameters?
      loss = getLoss(A, partitions[i], rand(NumOriginalNodes), modelType, NumOriginalNodes, 0.01; listModelArgs...)
      z[i] = loss
      println(i)
    end
    z = Array(z)

    #Save vector data if we want to smooth it later
    if !isempty(save_to_string)
      df = DataFrame(["x" => x, "y" => y, "z" => z])
      loc = "./data/visualization_data/" * save_to_string * ".csv"
      CSV.write(loc, df)
      df = DataFrame(["x" => x, "y" => y, "z" => z, "partition" => partitions])
      loc = "./data/visualization_data/PART" * save_to_string * ".csv"
      CSV.write(loc, df)
    end
end
