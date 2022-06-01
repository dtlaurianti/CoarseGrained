using MatrixNetworks
using LinearAlgebra

function linear_model(t, x, A, epsilon)
    (epsilon*A-np.eye(len(A))).dot(x)
end

function SIS_model(t, x, A, gamma, beta)
  n = np.size(A, axis=0)
  -gamma*x + beta*np.multiply((np.ones(n)-x), A.dot(x))
end

function SI_model(t, x, A, beta)
  n = np.size(A, axis=0)
  beta*np.multiply((np.ones(n)-x), A.dot(x))
end
