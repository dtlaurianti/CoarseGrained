using MatrixNetworks
using LinearAlgebra

function linear_model(t, x, A, epsilon)
    (epsilon*A-I)⋅x
end

function SIS_model(t, x, A, gamma, beta)
  n = np.size(A, axis=0)
  -gamma*x + beta*np.multiply((np.ones(n)-x), A.dot(x))
end

function SI_model(t, x, A, beta)
  n = size(A, 1)
  beta*multiply((ones(n)-x), A⋅x)
end
