using MatrixNetworks
using LinearAlgebra
using DifferentialEquations
using SparseArrays

# p = A
function linear_model(du::Vector, u::Vector, p::SparseMatrixCSC, t::Number)
  du .= (p-I)*u
end
#=
function SIS_model(t::Number, x::Vector, A::MatrixNetwork, γ::Number, β::Number)
  nrow = size(A, 1)
  return -γ.*x + β.*[(ones(n).-x).*(A⋅x)]
end
=#
# p = (A, γ, β)
function SIS_model(du::Vector, u::Vector, p::Tuple{SparseMatrixCSC, Number, Number}, t::Number)
  du .= -p[2].*u + p[3].*[(ones(n)-u).*(p[1]*u)]
end

function SI_model(t::Number, x::Vector, A::MatrixNetwork, β::Number)
  n = size(A, 1)
  return β.*((ones(n).-x).*(A⋅x))
end

function kuramoto_model(t::Number, x::Vector, A; ω::Number=nothing, K::Number=1)
    n = size(A, 1)
    dxdt = ones(n).*ω
    for i in 1:n
        S = sin(x[i].*ones(n).-(x))
        dxdt .+= K.*(A[i,:].*(S))
    end
    return dxdt
end

function LotkaVolterra_model(t::Number, x::Vector, A::MatrixNetwork, ω::Number)
  n = size(A, 1)
  return ω.*x + (x.*(A⋅x))
end

function linear_opinions(t::Number, x::Vector, A::MatrixNetwork, c::Number=1)
  #=
  Simplest model of linear opinion dynamics; c = constant, for now
  =#
  Dout = diag(sum(A, 1))            # compute outdegree matrix
  Din  = diag(sum(A, 2))            # compute indegree matrix
  L =  Dout + Din - A - transpose(A)  # compute graph laplacian
  return -c.*L⋅x
end

function nonlinear_opinions(t::Number, x::Vector, A::MatrixNetwork; d::Number=0.1, u::Number=1, b::Number=0)
  #=
  An instance of the model in https://arxiv.org/abs/2009.04332
  Parameters:
      d = resistance,
      u = attention to social influence
      b = bias
  For now, d, u, b are constants, but we can eventually vary by individual
  =#
  n = size(A, 1)
  return - d.*x + u.*tanh.(A⋅x) + b.*ones(n)
end

function simulateODEonGraph(A::MatrixNetwork, initial_condition::Vector; dynamical_function::Function=linear_model, tmax::Number=10, dt::Number=0.01, function_args...)
  tspan = (0, tmax)
  prob = ODEProblem(dynamical_function, initial_condition, tspan, sparse(A))
  sol = solve(prob, saveat=dt)
  # time_series = solve_ivp(t, x -> dynamical_function(t, x, A, function_args), (0, tmax), initial_condition, t_eval=t)
  return sol
end

function foldername(dynSys_string, dynSys_args, graphModel_string, graphModel_args)
  # function to generate
  dynSys_args_string = join([key*string(get(dynSys_args,key)) for key in sort(keys(dynSys_args))], "_")
  graphModel_args_string = join([key*string(get(graphModel_args, key)) for key in sort(keys(graphModel_args))], "_")

  string1 = dynSys_string*"_"*dynSys_args_string
  string2 = graphModel_string*"_"*graphModel_args_string

  return "data/"*string1*"/"*string2
end
