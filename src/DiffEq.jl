using DifferentialEquations
using Plots
include("SimulateDynamics.jl")

A = [0.0 1.0 0.0 0.0
     1.0 0.0 1.0 0.0
     0.0 1.0 0.0 1.0
     0.0 0.0 1.0 0.0]

Am = sparse(A)

u0 = rand(4,1)

tspan = (0.0,1.0)

p=getModelParameters(Am)

prob = ODEProblem(linear_model, u0, p, tmax=1)
sol = solve(prob)
#plot(sol)
