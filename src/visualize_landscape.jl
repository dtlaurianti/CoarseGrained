#January 2021
#Author(s): Nicole Eikmeier, Alice Schwarze, Nicholas Landry, Phil Chodrow, Mari Kawakatsu, Benji Zusman

#This method is inspired/taken by/from authors Leto Peel, Daniel B. Larremore, and Aaron Clauset in the paper 'The ground truth about metadata and community detection in networks'

# STEPS:
# 1) Partition Sampling
#    Sample a number of partitions and save the log-likelihood
# 2) Data Projection
#    Project the K^N dimensional partition data down to two dimensions using Multi-dimensional Scaling (MDS)
#    https://scikit-learn.org/stable/modules/generated/sklearn.manifold.MDS.html
#    This gives a two-dimensional representation of the partition space that preserves the variation of information between partitions
# 3) Surface Interpolation
#    Fit an Interpolated surface to the data
#    May need to smooth using Gaussian noise

using Clustering
using MultivariateStats
using DataFrames
using CSV
using ScikitLearn
using StatsBase
using Plots
using PyPlot
using SharedArrays


# convert a partition from our dictionary supernode format to a nested array format
# functionality changed from the dict_to_array function because our dissimilarity formula
# is impartial to the ordering of the nodes.
#= Old:
#Function dict_to_array
#
#Input: A list of dictionary partitions
#Output: An array of partitions
function dict_to_array(dictionary::Array{Dict,1})
    Arr = []
    for i in range(size(dictionary,1))
        s = findmax(dictionary[i])[1]
        x = [[] for i in 1:(s+1)]
        for (key, value) in dictionary
            v = dictionary[i][key]
            append!(x[v], key)
        end
        append!(Arr, x)
    end
  return Arr
end
=#

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
function variation_of_information(X::Vector{Vector{Any}},Y::Vector{Vector{Any}})
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



#Function: surfaceplots
#Parameters: partitions, an array of dictionaries representing network partitions
#            A -- MatrixNetwork representation of a network
#            NumOriginalNodes -- the number of nodes in A, which should be the same as the number
#            of original nodes in all of the partitions.
#            save_to_string -- (optional) name of the string that the plotted data should be saved to
#            in CSV format. If left empty, the data will not be saved to a file.
#            modelType -- (optional) denotes which model will be used to calculate loss.
#            Default model is linear_model.
#Purpose: To plot a 3d surface representing the loss landscape of a range of partitions
#Return value: none. Plots a graph and saves the (x, y, z) data in a CSV file if save_to_string
#              is provided a value.
function surfaceplots(partitions::Vector{Dict{Integer, Integer}}, A, NumOriginalNodes; save_to_string="", modelType::Function=linear_model)
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
    z = zeros(num_par)
    for i in 1:num_par
        # using hard-coded model and parameters, possibly want to make the outer function accept those parameters?
      loss = getLoss(A, partitions[i], rand(NumOriginalNodes), modelType, NumOriginalNodes, 0.01; listModelArgs...)
      z[i] = loss
    end

    #Save vector data if we want to smooth it later
    if !isempty(save_to_string)
      df = DataFrame(["x" => x, "y" => y, "z" => z])
      loc = "./data/visualization_data/" * save_to_string * ".csv"
      CSV.write(loc, df)
      df = DataFrame(["x" => x, "y" => y, "z" => z, "partition" => partitions])
      loc = "./data/visualization_data/PART" * save_to_string * ".csv"
      CSV.write(loc, df)
    end

    # plot surface
    pyplot()
    pygui(true)
    plt = Plots.plot(x, y, z, st=:surface, extra_kwargs=Dict(:subplot=>Dict("3d_colorbar_axis" => [0.85, 0.05, 0.05, 0.9])))
    display(plt)

end

#function plot_smoothed_surface
#
#Input: A csv file which has been smoothed by the R function
#Output: No output, just a plot
function plot_smoothed_surface(data)
  #Grab data from file
  df = DataFrame(CSV.File(data))
  x = df.x
  y = df.y
  z = df.z

  # plot surface
  pyplot()
  pygui(true)
  plt = Plots.plot(x, y, z, st=:surface, extra_kwargs=Dict(:subplot=>Dict("3d_colorbar_axis" => [0.85, 0.05, 0.05, 0.9])))
  display(plt)
end


# function subplot_smoothed_surface
# a subplot version of plot_smoothed_surface
# Input: A csv file which has been smoothed by the R function, (optional) ax
# Output: A subplot
function subplot_smoothed_surface(data, fig, ax=None)
  #Grab data from file
  data = pd.read_csv(data)
  x = data['x']
  y = data['y']
  z = data['z']

  #Plot surface
  triang = mtri.Triangulation(x,y)
  if ax != None
      ax = fig.add_subplot(1,2,ax,projection="3d")
  else
      ax = fig.add_subplot(1,2,1,projection="3d")
  end
  surf = ax.plot_trisurf(triang,z,cmap=plt.cm.CMRmap,antialiased=True)
  return surf
end

#Test Examples
#G = nx.fast_gnp_random_graph(10, 0.3)
#A = nx.to_numpy_array(G)
#P = RN.generateRandomPartitions(10,7,1000)

#surfaceplots(P,A,"n10_p7_s1000")

# plot_smoothed_surface("data/visualization_data/test.csv")
