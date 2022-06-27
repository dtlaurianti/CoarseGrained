# create a common Type for passing parameters
# A is the matrix that represents our network
# all other parameters are optional,
# and the model will only use the parameters it needs
@everywhere Base.@kwdef struct Model_Parameters
  A::SparseMatrixCSC
  ϵ::Number=1
  β::Number=0
  γ::Number=0
  ω::Vector=[missing]
  K::Number=1
  d::Number=0.1
  c::Number=1
  b::Number=0
end

# write functions to modify du element-wise,
# assigning du to another vector will not modify it in-place but reassign its pointer (?)
# this causes unexpected behavior in the simulateODEonGraph function



#Function: linear_model
#Parameters: du, passed by the ODESolver, modified in place
#            u, passed by the ODESolver, our initial_conditions on first pass, then whatever they become
#            p, the struct which we will use A and ϵ from
#            t, the time variable ODESolver will use
#Purpose: to simulate a linear model with the given inputs
#Return value: Modifies in place, return value not used
@everywhere function linear_model(du::Vector, u::Vector, p::Model_Parameters, t::Number)
  du .= (p.ϵ.*p.A-I)*u
end

#=
function SIS_model(t, x, A, gamma, beta)
  n = np.size(A, axis=0)
  return -gamma*x + beta*np.multiply((np.ones(n)-x), A.dot(x))
  =#


  #Function: SIS_model
  #Parameters: du, passed by the ODESolver, modified in place
  #            u, passed by the ODESolver, our initial_conditions on first pass, then whatever they become
  #            p, the struct which we will use A, γ, and β from
  #            t, the time variable ODESolver will use
  #Purpose: to simulate an SIS model with the given inputs
  #Return value: Modifies in place, return value not used
@everywhere function SIS_model(du::Vector, u::Vector, p::Model_Parameters, t::Number)
  du .= -p.γ.*u + p.β.*((ones(length(u))-u).*(p.A*u))
end
#=
function SI_model(t, x, A, beta)
  n = np.size(A, axis=0)
  return beta*np.multiply((np.ones(n)-x), A.dot(x))
=#


#Function: SI_model
#Parameters: du, passed by the ODESolver, modified in place
#            u, passed by the ODESolver, our initial_conditions on first pass, then whatever they become
#            p, the struct which we will use A and β from
#            t, the time variable ODESolver will use
#Purpose: to simulate an SI model with the given inputs
#Return value: Modifies in place, return value not used
@everywhere function SI_model(du::Vector, u::Vector, p::Model_Parameters, t::Number)
  du .= p.β.*((ones(length(u))-u).*(p.A*u))
end


#Function: kuramoto_model
#Parameters: du, passed by the ODESolver, modified in place
#            u, passed by the ODESolver, our initial_conditions on first pass, then whatever they become
#            p, the struct which we will use A, K, and ω from
#            t, the time variable ODESolver will use
#Purpose: to simulate an Kuramoto model with the given inputs
#Return value: Modifies in place, return value not used
@everywhere function kuramoto_model(du::Vector, u::Vector, p::Model_Parameters, t::Number)
    n = size(p.A, 1)
    du .= ones(n).*p.ω
    for i in 1:n
        S = sin.(u[i].*ones(n).-(u))
        du .+= p.K.*(Vector(p.A[i,:]).*(S))
    end
end


#Function: LotkaVolterra_model
#Parameters: du, passed by the ODESolver, modified in place
#            u, passed by the ODESolver, our initial_conditions on first pass, then whatever they become
#            p, the struct which we will use A and ω from
#            t, the time variable ODESolver will use
#Purpose: to simulate an Lotka-Volterra model with the given inputs
#Return value: Modifies in place, return value not used
@everywhere function LotkaVolterra_model(du::Vector, u::Vector, p::Model_Parameters, t::Number)
  n = size(p.A, 1)
  du .= p.ω.*u + u.*(p.A*u)
end


#Function: linear_opinions
#Parameters: du, passed by the ODESolver, modified in place
#            u, passed by the ODESolver, our initial_conditions on first pass, then whatever they become
#            p, the struct which we will use A and c from
#            t, the time variable ODESolver will use
#Purpose: to simulate an linear opinions model with the given inputs
#Return value: Modifies in place, return value not used
@everywhere function linear_opinions(du::Vector, u::Vector, p::Model_Parameters, t::Number; c=1)
  #=
  Simplest model of linear opinion dynamics; c = constant, for now
  =#
  Dout = diag(sum(p.A, dims=1))            # compute outdegree matrix
  Din  = diag(sum(p.A, dims=2))            # compute indegree matrix
  L =  Dout .+ Din .- p.A .- transpose(p.A)  # compute graph laplacian
  du .= -p.c.*sum((L*u), dims=2)
end


#Function: nonlinear_opinions
#Parameters: du, passed by the ODESolver, modified in place
#            u, passed by the ODESolver, our initial_conditions on first pass, then whatever they become
#            p, the struct which we will use A, d, c, and b from
#            t, the time variable ODESolver will use
#Purpose: to simulate an nonlinear opinion model with the given inputs
#Return value: Modifies in place, return value not used
@everywhere function nonlinear_opinions(du::Vector, u::Vector, p::Model_Parameters, t::Number)
  #=
  An instance of the model in https://arxiv.org/abs/2009.04332
  Parameters:
      d = resistance,
      c = attention to social influence
      b = bias
  For now, d, c, b are constants, but we can eventually vary by individual
  =#
  n = size(p.A, 1)
  du .= -p.d.*u + p.c.*tanh.(p.A*u) + p.b.*ones(n)
end


#Function: SIS_model
#Parameters: A, the network to simulate dynamics on
#            initial_condition, the vector of initial variables for each node
#            dynamical_function, the model which we will use to calculate the dynamics
#            tmax, when the simulation should end
#            dt, the timesteps to numerically compute for the simulation
#            function_args, the extra inputs to the dynamical model (ϵ, γ, β, K, ω, d, c, b)
#Purpose: to simulate the dynamics on a network using the chosen dynamical model and conditions
#Return value: returns an ODESolution or ODECompositeSolution type, accessible through an array interface
@everywhere function simulateODEonGraph(A::MatrixNetwork, initial_condition::Vector; dynamical_function::Function=linear_model, tmax::Number=10, dt::Number=0.01, function_args...)
  tspan = (0, tmax)
  # splat the MatrixNetwork and additional parameters into one parameters Tuple
  p = Model_Parameters(A=sparse(A); function_args...)
  prob = ODEProblem(dynamical_function, initial_condition, tspan, p)
  sol = solve(prob, saveat=dt)
  # time_series = solve_ivp(t, x -> dynamical_function(t, x, A, function_args), (0, tmax), initial_condition, t_eval=t)
  return sol
end
