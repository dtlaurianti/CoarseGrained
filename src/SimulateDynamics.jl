using MatrixNetworks
using LinearAlgebra
using DifferentialEquations
using SparseArrays
using Parameters

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

function linear_model(du::Vector, u::Vector, p::Model_Parameters, t::Number)
  du .= (p.A-I)*u
end

function SIS_model(du::Vector, u::Vector, p::Model_Parameters, t::Number)
  du .= -p.γ.*u + p.β.*(ones(length(u))-u).*(p.A*u)
end

function SI_model(du::Vector, u::Vector, p::Model_Parameters, t::Number)
  du .= p.β.*(ones(length(u))-u).*(p.A*u)
end

function kuramoto_model(du::Vector, u::Vector, p::Model_Parameters, t::Number)
    n = size(p.A, 1)
    dxdt = ones(n).*p.ω
    for i in 1:n
        S = sin.(u[i].*ones(n).-(u))
        dxdt .+= p.K.*(p.A[i,:].*(S))
    end
    return dxdt
end

function LotkaVolterra_model(du::Vector, u::Vector, p::Model_Parameters, t::Number)
  n = size(p.A, 1)
  return p.ω.*u + (u.*(p.A⋅u))
end

function linear_opinions(du::Vector, u::Vector, p::Model_Parameters, t::Number)
  #=
  Simplest model of linear opinion dynamics; c = constant, for now
  =#
  Dout = diag(sum(p.A, 1))            # compute outdegree matrix
  Din  = diag(sum(p.A, 2))            # compute indegree matrix
  L =  Dout + Din - p.A - transpose(p.A)  # compute graph laplacian
  return -c.*L⋅x
end

function nonlinear_opinions(du::Vector, u::Vector, p::Model_Parameters, t::Number)
  #=
  An instance of the model in https://arxiv.org/abs/2009.04332
  Parameters:
      d = resistance,
      c = attention to social influence
      b = bias
  For now, d, c, b are constants, but we can eventually vary by individual
  =#
  n = size(p.A, 1)
  return - p.d.*u + p.c.*tanh.(p.A⋅u) + p.b.*ones(n)
end

function simulateODEonGraph(A::MatrixNetwork, initial_condition::Vector; dynamical_function::Function=linear_model, tmax::Number=10, dt::Number=0.01, function_args...)
  tspan = (0, tmax)
  # splat the MatrixNetwork and additional parameters into one parameters Tuple
  p = Model_Parameters(A=sparse(A); function_args...)
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