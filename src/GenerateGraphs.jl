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
  # sum the number of k-walks between each node pair, add I to not require a path from a node to itself
  Asum = sum(A^k for k=1:(n-1)) + I
  # check if any of the nodes have no path between them
  return sum(Asum .≈ 0) == 0
end

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

function gnp_graph(n::Int; p::AbstractFloat=0.1, directed::Bool=true, edge_weight::Number=1.0)
  if directed
    return MatrixNetwork(sprand(n,n,p,bitrand).*edge_weight)
  else
    G = sprand(n,n,p,bitrand)
    # Mirroring Upper and Lower triangles to make the network undirected
    return max.(G, G').*edge_weight
  end
end

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

function cm_graph(n; max_degree=5)
  nodesPerDegree = n ÷ max_degree
  degreeArray = Array{Int64, 1}(undef, n)
  degree = 1
  iterator = 1
  for i in 1:n
    degreeArray[i] = degree
    #push!(degreeArray, degree)
    if iterator >= nodesPerDegree
      if i == (n-1)
        degree -= 1
      end
      iterator = 1
      degree += 1
      @goto skip
    end
    iterator += 1
    @label skip
  end
  G = random_configuration_model(n, degreeArray)
  return MatrixNetwork(sparse((G)))
end

# TODO: add something like a feed-forward neural-network graph?
# or something else with hierarchical or layered structure?
