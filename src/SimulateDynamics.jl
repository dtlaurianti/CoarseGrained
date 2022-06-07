using MatrixNetworks
using LinearAlgebra
using DifferentialEquations
using SparseArrays

#=
# create a common Type for passing parameters
Base.@kwdef struct Model_Parameters
  A::SparseMatrixCSC
  β::Number=0
  γ::Number=0
  ω::Number=0
  K::Number=1
  d::Number=0.1
  c::Number=1
  b::Number=0
end
=#

function getModelParameters(A::SparseMatrixCSC;  β::Number=0, γ::Number=0, ω::Number=0, K::Number=1, d::Number=0.1, c::Number=1, b::Number=0)
   Model_Parameters = [A, β, γ, ω, K, d, c, b]
   return Model_Parameters
end

function linear_model(du::Vector, u::Vector, p::Array, t::Number)
  du .= (p[1]-I)*u
end

function SIS_model(du::Vector, u::Vector, p::Array, t::Number)
  du .= -p[3].*u + p[2].*(ones(length(u))-u).*(p[1]*u)
end

function SI_model(du::Vector, u::Vector, p::Array, t::Number)
  du .= p[2].*(ones(length(u))-u).*(p[1]*u)
end

function kuramoto_model(du::Vector, u::Vector, p::Array, t::Number)
    n = size(p[1], 1)
    dxdt = ones(n).*p[4]
    for i in 1:n
        S = sin.(u[i].*ones(n).-(u))
        dxdt .+= p[5].*(p[1][i,:].*(S))
    end
    return dxdt
end

function LotkaVolterra_model(du::Vector, u::Vector, p::Array, t::Number)
  n = size(p[1], 1)
  return p[4].*u + (u.*(p[1]⋅u))
end

function linear_opinions(du::Vector, u::Vector, p::Array, t::Number)
  #=
  Simplest model of linear opinion dynamics; c = constant, for now
  =#
  Dout = diag(sum(p[1], 1))            # compute outdegree matrix
  Din  = diag(sum(p[1], 2))            # compute indegree matrix
  L =  Dout + Din - p[1] - transpose(p[1])  # compute graph laplacian
  return -c.*L⋅x
end

function nonlinear_opinions(du::Vector, u::Vector, p::Array, t::Number)
  #=
  An instance of the model in https://arxiv.org/abs/2009.04332
  Parameters:
      d = resistance,
      c = attention to social influence
      b = bias
  For now, d, c, b are constants, but we can eventually vary by individual
  =#
  n = size(p[1], 1)
  return - p[6].*u + p[7].*tanh.(p[1]⋅u) + p[8].*ones(n)
end

function simulateODEonGraph(A::MatrixNetwork, initial_condition::Vector; dynamical_function::Function=linear_model, tmax::Number=10, dt::Number=0.01, function_args...)
  tspan = (0, tmax)
  # splat the MatrixNetwork and additional parameters into one parameters Tuple
  p = getModelParameters(sparse(A); function_args...)
  prob = ODEProblem(dynamical_function, initial_condition, tspan, p)
  sol = solve(prob, saveat=dt)
  # time_series = solve_ivp(t, x -> dynamical_function(t, x, A, function_args), (0, tmax), initial_condition, t_eval=t)
  return sol
end

function foldername(dynSys_string, dynSys_args, graphModel_string, graphModel_args)
  # function to generate
  dynSys_args_string = join([key*string(get(dynSys_args, key, "Error")) for key in sort(keys(dynSys_args))], "_")
  graphModel_args_string = join([key*string(get(graphModel_args, key, "Error")) for key in sort(keys(graphModel_args))], "_")

  string1 = dynSys_string*"_"*dynSys_args_string
  string2 = graphModel_string*"_"*graphModel_args_string

  return "data/"*string1*"/"*string2
end
