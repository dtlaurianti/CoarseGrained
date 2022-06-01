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
  beta*multiply((ones(n)-x), A⋅x)
end
