using LinearAlgebra
using MatrixNetworks
using LoopVectorization
using BenchmarkTools
using SparseArrays
using Random
using Graphs

# returns true if A has paths between every node pair
function isConnected(A::MatrixNetwork)
  A = sparse(A)
  n = size(A,1)
  # sum the number of k-walks between each node pair
  Asum = sum(A^k for k=1:(n-1))
  # check if any of the nodes have no path between them
  return !(Asum .≈ 0)
end

#Function: line_graph
#Parameters: n, integer number of nodes
#            edge_weight, (optional) integer/floating point weight of edges in graph
#            directed, (optional) bloolean (true/false) denoting whether the graph should be directed
#Purpose: To generate a line graph with n nodes connected in a line.
#Return value: MatrixNetwork representation of a line graph
# example graphs (non random)
function line_graph(n::Int; edge_weight::Number=1.0, directed::Bool=true)
  if n < 0
    throw(DomainError("n must be >= 0."))
  elseif n == 0
    return MatrixNetwork(sparse(Array{Float64}(undef, 0, 0)))
  else
    G = zeros(n, n)
    if directed
      # Create a directed path from the first node to the last node
      for i=1:(n-1)
        G[i,i+1] = edge_weight
      end
    else
      # Create an undirected path from the first node to the last node
      for i=1:(n-1)
        G[i,i+1] = edge_weight
        G[i+1,i] = edge_weight
      end
    end
    return MatrixNetwork(sparse(G))
  end
end

#Function: cycle_graph
#Parameters: n, integer number of nodes
#            edge_weight, (optional) integer/floating point weight of edges in graph
#            directed, (optional) bloolean (true/false) denoting whether the graph should be directed
#Purpose: To generate a cycle graph with n nodes connected in a cycle that includes all nodes. 
#         All nodes have degree 2.
#Return value: MatrixNetwork representation of a cycle graph
function cycle_graph(n::Int; edge_weight::Number=1.0, directed::Bool=true)
  if n < 0
    throw(DomainError("n must be >= 0."))
  elseif n == 0
    return MatrixNetwork(sparse(Array{Float64}(undef, 0, 0)))
  else
    G = zeros(n, n)
    if directed
      # Create a directed cycle from the first node to the last node
      for i=1:(n-1)
        G[i, i+1] = edge_weight
      end
      G[n, 1] = edge_weight
    else
      # Create an undirected cycle from the first node to the last node
      for i=1:(n-1)
        G[i,i+1] = edge_weight
        G[i+1,i] = edge_weight
      end
      G[n, 1] = edge_weight
      G[1, n] = edge_weight
    end
    return MatrixNetwork(sparse(G))
  end
end

#Function: grid_graph
#Parameters: n, integer number of nodes
#            edge_weight, (optional) integer/floating point weight of edges in graph
#            directed, (optional) bloolean (true/false) denoting whether the graph should be directed
#Purpose: To generate a grid-shaped graph that is as close to a perfect 2-d lattice as possible with n nodes
#Return value: MatrixNetwork representation of a grid graph
function grid_graph(n::Int; edge_weight::Number=1.0, directed::Bool=false)
  #=A square (or close to square) 2d lattice with a number of nodes close to `n`. The
  directed version of the grid graph has a source and a sink node (e.g., everything
  flows from top left to bottom right on the grid).=#
  if n < 0
    throw(DomainError("n must be >= 0."))
  elseif n == 0
    return MatrixNetwork(sparse(Array{Float64}(undef, 0, 0)))
  else
    # Choose rectangular dimensions to have a number of nodes close to n
    m = sqrt(n)
    mint = trunc(Int, m)
    if m == mint
      m1, m2 = mint, mint
    else
      m1, m2 = mint+1, mint
    end
    n = m1*m2
    G = zeros(n, n)
    for i=1:n
      if i+m2 <= n
        G[i, i+m2] = edge_weight
      end
      if i+1 <=n && (i+1)%m2 != 1
        G[i, i+1] = edge_weight
      end
    end
    if !directed
      for i=1:n, j=1:n
        G[j, i] = G[i, j]
      end
    end
    return MatrixNetwork(sparse(G))
  end
end

# random graphs

#Function: gnp_graph
#Parameters: n, integer number of nodes
#            p, floating point number (maximum 1.0) that represents the probability that any node is connected to any other node
#            directed, (optional) bloolean (true/false) denoting whether the graph should be directed
#            edge_weight, (optional) integer/floating point weight of edges in graph
#Purpose: To generate a random graph using a user-specified probability that nodes will be connected to each other.
#         This is basically the Erdos-Renyi model.
#Return value: MatrixNetwork representation of a gnp graph
function gnp_graph(n::Int; p::AbstractFloat=0.1, directed::Bool=true, edge_weight::Number=1.0)
  if directed
    return MatrixNetwork(sprand(n,n,p,bitrand).*edge_weight)
  else
    G = sprand(n,n,p,bitrand)
    # Mirroring Upper and Lower triangles to make the network undirected
    return max.(G, G').*edge_weight
  end
end

#Function: sbm_graph
#Parameters: n, integer number of nodes
#            communities, integer number of distinct communities within the stochastic block model
#            p_within, floating point number (maximum 1.0) that represents the probability that any node is connected
#            to any other node within the same community
#            p_between, floating point number (maximum 1.0) that represents the probability that any node is connected
#            to any other node outside of its community
#            directed, (optional) bloolean (true/false) denoting whether the graph should be directed
#            edge_weight, (optional) integer/floating point weight of edges in graph
#Purpose: To generate a random stochastic block model graph with n nodes and communities communities
#         using the user-specified probabilities of within and between community connections.
#Return value: MatrixNetwork representation of an sbm graph
function sbm_graph(n; communities=4, p_within=0.2, p_between=0.05, edge_weight=1.0, directed=true)
  if directed
    if communities == 1
      G = gnp_graph(n, p=p_within)
      return G
    end
    G = sparse(gnp_graph(n÷communities, p=p_within))
    i = 1
    while i != communities
      if (i + 1) == communities
        G = blockdiag(G, sparse(gnp_graph(n÷communities + n%communities, p=p_within)))
        i += 1
      else
        G = blockdiag(G, sparse(gnp_graph(n÷communities, p=p_within)))
        i += 1
      end
    end
  else
    if communities == 1
      G = gnp_graph(n, p=p_within, directed=false)
      return G
    end
    G = sparse(gnp_graph(n÷communities, p=p_within, directed=false))
    i = 1
    while i != communities
      if (i + 1) == communities
        G = blockdiag(G, sparse(gnp_graph(n÷communities + n%communities, p=p_within, directed=false)))
        i += 1
      else
        G = blockdiag(G, sparse(gnp_graph(n÷communities, p=p_within, directed=false)))
        i += 1
      end
    end
  end
  i = 1
  mask = sparse(ones(n÷communities, n÷communities))
  while i != communities
    if (i + 1) == communities
      mask = blockdiag(mask, sparse(ones(n÷communities + n%communities, n÷communities + n%communities)))
      i+= 1
    else
      mask = blockdiag(mask, sparse(ones(n÷communities, n÷communities)))
      i+= 1
    end
  end
  G = Matrix(G)
  mask = Matrix(mask)
  num_rows = n
  num_columns = n
  for i in 1:num_rows
    for j in 1:num_columns
      if G[i, j] == 0 || G[i, j] == 0.0 && mask[i, j] == 0 || mask[i, j] == 0.0
        portion = p_between * 10000
        if rand(big.(1:10000)) <= portion
          G[i, j] = 1
        end
      end
    end
  end
  return MatrixNetwork(sparse(G))
end

#Function: cm_graph
#Parameters: n, integer number of nodes
#            degreeArray, an array of integers specifying the degree of each node in the graph.
#            degreeArray must not contain any number greater than or equal to n, degreeArray
#            must be of length n, and the sum of all elements in degreeArray must be an even number
#Purpose: To generate a random configuration model graph with n nodes each having degree specified in degreeArray
#Return value: MatrixNetwork representation of a configuration model network
function cm_graph(n, degreeArray)
  G = random_configuration_model(n, degreeArray)
  return MatrixNetwork(sparse((G)))
end

# TODO: add something like a feed-forward neural-network graph?
# or something else with hierarchical or layered structure?
