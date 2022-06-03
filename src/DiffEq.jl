using DifferentialEquations
using Plots
using ParameterizedFunctions

A = [0.0 1.0 0.0 0.0
     1.0 0.0 1.0 0.0
     0.0 1.0 0.0 1.0
     0.0 0.0 1.0 0.0]

Am = MatrixNetwork(sparse(A))

u0 = rand(4,1)

tspan = (0.0,1.0)

function linear_model(du, u, p, t)
  du .= (p-I)*u
end

prob = ODEProblem(linear_model, u0, tspan, Am)
sol = solve(prob)
plot(sol)
