using MatrixNetworks
using LinearAlgebra

function linear_model(t, x, A, ϵ)
  (ϵ*A-I)⋅x
end

function SIS_model(t, x, A, γ, β)
  nrow = size(A, 1)
  -γ*x + β*[(ones(n)-x).*(A⋅x)]
end

function SI_model(t, x, A, β)
  n = size(A, 1)
  β*((ones(n)-x).*(A⋅x))
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

function linear_opinions(t, x, A, c=1):
  #=
  Simplest model of linear opinion dynamics; c = constant, for now
  =#
  Dout = diag(sum(A, 1))            # compute outdegree matrix
  Din  = diag(sum(A, 2))            # compute indegree matrix
  L =  Dout + Din - A - transpose(A)  # compute graph laplacian
  return -c*L⋅x

function nonlinear_opinions(t, x, A, d=0.1, u=1, b=0):
  #=
  An instance of the model in https://arxiv.org/abs/2009.04332
  Parameters:
      d = resistance,
      u = attention to social influence
      b = bias
  For now, d, u, b are constants, but we can eventually vary by individual
  =#
  n = size(A, 1)
  return - d*x + u*tanh(A⋅x)) + b*ones(n)
