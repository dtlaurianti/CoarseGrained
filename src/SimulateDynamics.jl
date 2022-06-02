using MatrixNetworks
using LinearAlgebra

function linear_model(t::Number, x::Vector, A::MatrixNetwork, ϵ::Number)
  return (ϵ.*A-I)⋅x
end

function SIS_model(t, x, A, γ, β)
  nrow = size(A, 1)
  return -γ*x + β*[(ones(n)-x).*(A⋅x)]
end

function SI_model(t, x, A, β)
  n = size(A, 1)
  return β*((ones(n)-x).*(A⋅x))
end

function kuramoto_model(t, x, A; ω=nothing, K=1)
    n = size(A, 1)
    dxdt = ones(n)*ω
    for i in 1:n
        S = sin(x[i]*ones(n).-(x))
        dxdt += K*(A[i,:].*(S))
    end
    return dxdt
end

function LotkaVolterra_model(t, x, A, ω)
  n = size(A, 1)
  return ω*x + (x.*(A⋅(x)))
end

function linear_opinions(t, x, A, c=1)
  #=
  Simplest model of linear opinion dynamics; c = constant, for now
  =#
  Dout = diag(sum(A, 1))            # compute outdegree matrix
  Din  = diag(sum(A, 2))            # compute indegree matrix
  L =  Dout + Din - A - transpose(A)  # compute graph laplacian
  return -c*L⋅x
end

function nonlinear_opinions(t, x, A; d=0.1, u=1, b=0)
  #=
  An instance of the model in https://arxiv.org/abs/2009.04332
  Parameters:
      d = resistance,
      u = attention to social influence
      b = bias
  For now, d, u, b are constants, but we can eventually vary by individual
  =#
  n = size(A, 1)
  return - d*x + u*tanh(A⋅x) + b*ones(n)
end

function simulateODEonGraph(A, initial_condition; dynamical_function=linear_model, tmax=10, dt=0.01, function_args...)
  t = range(0, tmax, dt)
  time_series = solve_ivp(t, x -> dynamical_function(t, x, A, function_args), (0, tmax), initial_condition, t_eval=t)
  return time_series
end

function foldername(dynSys_string, dynSys_args, graphModel_string, graphModel_args)
  # function to generate
  dynSys_args_string = join([key*string(get(dynSys_args,key)) for key in sort(keys(dynSys_args))], "_")
  graphModel_args_string = join([key*string(get(graphModel_args, key)) for key in sort(keys(graphModel_args))], "_")

  string1 = dynSys_string*"_"*dynSys_args_string
  string2 = graphModel_string*"_"*graphModel_args_string

  return "data/"*string1*"/"*string2
end
