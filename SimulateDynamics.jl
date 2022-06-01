using MatrixNetworks
using LinearAlgebra

function linear_model(t, x, A, ϵ)
  (ϵ*A-I)⋅x
end

function SIS_model(t, x, A, γ, β)
  nrow = size(A, 1)
  -γ*x + β*[(ones(n)-x).*(A⋅x)]
end

function SI_model(t, x, A, beta)
  n = size(A, 1)
  beta*((ones(n)-x).*(A⋅x))
end

def kuramoto_model(t, x, A, omega=None, K=1):
    n = np.size(A,axis=0)
    dxdt = np.ones(n)*omega #omega.copy()
    for i in range(n):
        S = np.sin(np.subtract(x[i]*np.ones(n),x))
        dxdt += K*np.multiply(A[i,:],S)
    return dxdt

def LotkaVolterra_model(t, x, A, omega):
  n = np.size(A, axis=0)
  return omega*x + np.multiply(x, A.dot(x))

def linear_opinions(t, x, A, c=1):
  '''
  Simplest model of linear opinion dynamics; c = constant, for now
  '''
  Dout = np.diag(np.sum(A, axis=0))            # compute outdegree matrix
  Din  = np.diag(np.sum(A, axis=1))            # compute indegree matrix
  L =  Dout + Din - A - np.matrix.transpose(A)  # compute graph laplacian
  return -c*L.dot(x)

def nonlinear_opinions(t, x, A, d=0.1, u=1, b=0):
  '''
  An instance of the model in https://arxiv.org/abs/2009.04332
  Parameters:
      d = resistance,
      u = attention to social influence
      b = bias
  For now, d, u, b are constants, but we can eventually vary by individual
  '''
  n = np.size(A, axis=0)
  return - d*x + u*np.tanh(A.dot(x)) + b*np.ones(n)