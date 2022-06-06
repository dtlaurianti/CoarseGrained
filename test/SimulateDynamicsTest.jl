using Plots
using GraphRecipes

@userplot NetworkPlot
@recipe function f(np::NetworkPlot)
end

LG = line_graph(5)
CG = cycle_graph(5)
GG = grid_graph(4; directed=true)

x = [1, 1, 1, 1, 0]
sol1 = simulateODEonGraph(CG, x)
plt1 = plot(5,xlim=(0,1000),ylim=(0,1),title="line_graph linear_model")
anim1 = @animate for i=1:1000
    push!(plt1,sol1[i])
end every 10
display(gif(anim1))
gplt1 = graphplot(sparse(LG))
display(gplt1)
sol2 = simulateODEonGraph(LG, x, γ=0.1, β=0.1, dynamical_function=SIS_model)
sol3 = simulateODEonGraph(LG, x, β=0.1, dynamical_function=SI_model)
sol4 = simulateODEonGraph(LG, x, dynamical_function=kuramoto_model)
