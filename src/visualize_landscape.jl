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

include("EvaluateError.jl")
include("ReduceNetwork.jl")
include("GenerateGraphs.jl")
include("SimulateDynamics.jl")



#Function variation_of_information
#
#Input: Two partitions
#Output: The distance between the two partitions
#
# Meila, M. (2007). Comparing clusterings-an information
#   based distance. Journal of Multivariate Analysis, 98,
#   873-895. doi:10.1016/j.jmva.2006.11.013
# https://gist.github.com/jwcarr/626cbc80e0006b526688
function variation_of_information(X, Y)
  n = float(sum([size(x,1) for x in X]))
  σ = 0.0
  for x in X
    p = size(x,1) / n
    for y in Y
      q = size(y,1) / n
      r = size(Set(x) & Set(y)) / n
      if r > 0.0
        σ += r * (log(r / p, 2) + log(r / q, 2))
      end
    end
  end
  return abs(σ)
end

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

#Function surfaceplots
#
#Inputs:
#  partitions is a dictionary of partitions
#  A is an adjacency matrix for a graph
#
#Outputs:
#  Plots a choppy 3D Surface
#
#Example of how to run:
#  G = nx.fast_gnp_random_graph(10, 0.3)
#  A = nx.to_numpy_array(G)
#  P = RN.generateRandomPartitions(10,7,100)
#  surfaceplots(P,A)
#
<<<<<<< HEAD
function surfaceplots(partitions, A, save_to_string=None):
=======
function surfaceplots(partitions, A, save_to_string=None)
>>>>>>> cbef35f036f5f01f8db8602e7fafe4a28826e941
    #convert dictionary to an array
    Arr = dict_to_array(partitions)

    #calculate distance matrix
    num_par = len(partitions)
    D = [[0 for i in range(num_par)] for j in range(num_par)]
    for i in range(num_par)
        for j in range(num_par)
          D[i][j] = variation_of_information(Arr[i],Arr[j])
        end
    end

    #calculate MDS on disimilarity matrix
    embedding = MDS(n_components=2, dissimilarity="precomputed")
    X_transformed = embedding.fit_transform(D)

    #Format data
    x = X_transformed[:,0]
    y = X_transformed[:,1]

    #Calculate z dimension
    z = zeros(num_par)
    for i in range(num_par)
      loss = EvaluateError.getLoss(A, partitions[i], ones(10), SimulateDynamics.linear_model, 10, 0.01, {"epsilon":-0.3})
      z[i] = loss
    end

    #Save vector data if we want to smooth it later
    if save_to_string != None
      df = pd.DataFrame({'x' : x, 'y' : y, 'z' : z})
      loc = "data/visualization_data/" + save_to_string + ".csv"
      df.to_csv(loc, index = False)
    end

    #Plot surface
    triang = mtri.Triangulation(x,y)
    fig = plt.figure()
    ax = fig.add_subplot(1,2,2,projection="3d")
    surf = ax.plot_trisurf(triang,z,cmap=plt.cm.CMRmap,antialiased=True)

    plt.show()
end

#Function surfaceplots2
#
#Inputs:
#  partitions is a dictionary of partitions
#  A is an adjacency matrix for a graph
#  parameters is a dictionary of simulation parameters (see PlotLandscape.py for example)
#Outputs:
#  Plots a choppy 3D Surface
#
#Example of how to run:
#  G = nx.fast_gnp_random_graph(10, 0.3)
#  A = nx.to_numpy_array(G)
#  P = RN.generateRandomPartitions(10,7,100)
#  surfaceplots(P,A)
#
function surfaceplots2(parameters, numSamples, save_to_string=None)
    #Generate partitions
    partitions = ReduceNetwork.generateRandomPartitions( parameters["originalSize"], parameters["reducedSize"], numSamples)

    #convert dictionary to an array
    Arr = dict_to_array(partitions)

    #calculate distance matrix
    num_par = length(partitions)
    D = [[0 for i in range(num_par)] for j in range(num_par)]
    for i in range(num_par)
        for j in range(num_par)
          D[i][j] = variation_of_information(Arr[i],Arr[j])
        end
    end

    #calculate MDS on disimilarity matrix
    embedding = MDS(n_components=2, dissimilarity="precomputed")
    X_transformed = embedding.fit_transform(D)

    #Format data
    x = X_transformed[:,0]
    y = X_transformed[:,1]

    #Generate graph (must be connected)
    isAccepted = false
    while !isAccepted
        if parameters["graphType"] == "random"
            G = nx.fast_gnp_random_graph(parameters["originalSize"], parameters["graphArgs"]["p"])
        elseif parameters["graphType"] == "cycle"
            G = GenerateGraphs.cycle_graph(parameters["originalSize"], directed=False)
        end
        if nx.is_connected(G)
            isAccepted = True
        end
    end

    #Convert to adj matrix
    A = nx.to_numpy_array(G)

    # select model
    if parameters["modelType"] == "linear_model"
        modelFunc = SimulateDynamics.linear_model
    elseif parameters["modelType"] == "SIS"
        modelFunc = SimulateDynamics.SIS_model
    elseif parameters["modelType"] == "SI"
        modelFunc = SimulateDynamics.SI_model
    elseif parameters["modelType"] == "kuramoto_model"
        modelFunc = SimulateDynamics.kuramoto_model
    end

    tmax, tinc = parameters["tmax"], parameters["tinc"]
    modelArgs  = parameters["modelArgs"]
    initial_condition = ones(10)

    #Calculate z dimension
    z = zeros(num_par)
    for i in range(num_par)
      loss = EE.getLoss(A, partitions[i], initial_condition, modelFunc, tmax, tinc, modelArgs)
      z[i] = loss
    end

    #Save vector data if we want to smooth it later
    if save_to_string != None
      df = pd.DataFrame({'x' : x, 'y' : y, 'z' : z})
      loc = "data/visualization_data/" + save_to_string + ".csv"
      df.to_csv(loc, index = False)
    end

    #Plot surface
    triang = mtri.Triangulation(x,y)
    fig = plt.figure()
    ax = fig.add_subplot(1,2,2,projection="3d")
    surf = ax.plot_trisurf(triang,z,cmap=plt.cm.CMRmap,antialiased=True)

    plt.show()
end

#function plot_smoothed_surface
#
#Input: A csv file which has been smoothed by the R function
#Output: No output, just a plot
function plot_smoothed_surface(data)
  #Grab data from file
  data = pd.read_csv(data)
  x = data['x']
  y = data['y']
  z = data['z']

  #Plot surface
  triang = mtri.Triangulation(x,y)
  fig = plt.figure()
  ax = fig.add_subplot(1,2,2,projection="3d")
  surf = ax.plot_trisurf(triang,z,cmap=plt.cm.CMRmap,antialiased=True)

  plt.show()
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

<<<<<<< HEAD
#plot_smoothed_surface("data/visualization_data/test.csv")
=======
# plot_smoothed_surface("data/visualization_data/test.csv")
>>>>>>> 8cf4b000d3b0996eb627efbf83d9d2447b83f887
