using LinearAlgebra
using MatrixNetworks

# example graphs (non random) 
function line_graph(n, edge_weight=1.0, directed=True):
  if directed:
    G = nx.DiGraph()
    nx.add_path(G, range(n), weight=edge_weight)
  else:
    G = nx.Graph()
    nx.add_path(G, range(n), weight=edge_weight)
  return G


def cycle_graph(n, edge_weight=1.0, directed=True):
  if directed:
    G = nx.DiGraph()
    nx.add_cycle(G, range(n), weight=edge_weight)
  else:
    G = nx.Graph()
    nx.add_cycle(G, range(n), weight=edge_weight)
  return G


def grid_graph(n, edge_weight=1.0, directed=False):
  '''A square (or close to square) 2d lattice with a number of nodes close to `n`. The
  directed version of the grid graph has a source and a sink node (e.g., everything
  flows from top left to bottom right on the grid).'''

  m = np.sqrt(n)
  if m == np.round(m):
    m1, m2 = int(m), int(m)
  else:
    m1, m2 = int(m)+1, int(m)
  G = nx.grid_2d_graph(m1, m2, periodic=False)

  if directed:
    H = nx.DiGraph()
    for e in G.edges():
      H.add_edge(min(e), max(e), weight=edge_weight)
    return H
  else:
    nx.set_edge_attributes(G, edge_weight, "weight")
    return G


# random graphs

def gnp_graph(n, p=0.1, directed=True, seed=None, edge_weight=1.0):
  G = nx.gnp_random_graph(n, p, directed=directed, seed=seed)
  nx.set_edge_attributes(G, edge_weight, "weight")
  return G


def sbm_graph(n, communities=4, p_within=0.2, p_between=0.05, seed=None, edge_weight=1.0, directed=True):

  # make array of block probabilities
  p = (p_within*np.eye(communities)
       +p_between*(np.ones(communities)-np.eye(communities)))

  # make list of community sizes
  size = int(n/communities)
  sizes = size*np.ones(communities, dtype=int)
  for s in sizes:
    if sum(sizes) < n:
      s += 1
    else:
      break

  # make graph
  G = nx.stochastic_block_model(sizes, p, nodelist=None, seed=seed,
                                directed=directed, selfloops=False, sparse=False)
  nx.set_edge_attributes(G, edge_weight, "weight")

  return G


def cm_graph(n, max_degree=5, directed=True, seed=None):
  '''(Directed) configuration model where there is roughly the same number of
  nodes with (in-)degree 1 as with (in-)degree 2 ... as with (in-)degree
  `max_degree`. For directed networks, the out-degree distribution is the same
  as the in-degree distribution but the in-degree and out-degree of a node are u
  ncorrelated.'''

  # make list of community sizes
  size = int(n/max_degree)
  sizes = size*np.ones(max_degree, dtype=int)
  for s in sizes:
    if sum(sizes) < n:
      s += 1
    else:
      break

  # make degree sequence
  indegrees = np.concatenate([(i+1)*np.ones(s, dtype=int) for i, s in enumerate(sizes)])
  if directed:
    outdegrees = np.copy(indegrees)
    np.random.shuffle(outdegrees)

    # make graph
    G = nx.directed_configuration_model(indegrees, outdegrees, create_using=nx.DiGraph(), seed=seed)

  else:
    # make graph
    G = nx.configuration_model(indegrees, seed=seed)

  return G

# TODO: add something like a feed-forward neural-network graph?
# or something else with hierarchical or layered structure?
