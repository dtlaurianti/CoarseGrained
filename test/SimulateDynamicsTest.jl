using Plots
using GraphRecipes

n = 16

LG = line_graph(n)
CG = cycle_graph(n)
GG = grid_graph(n; directed=true)
GNPG = gnp_graph(n, p=0.5)

x = rand(n)


scale = 1/log(n)

sol1 = simulateODEonGraph(GG, x)
plt1 = Plots.plot(sol1,xlim=(0,10),ylim=(0,5),title="grid_graph linear_model")
gplt1 = graphplot(sparse(GG), thickness_scaling=scale)
display(Plots.plot(plt1, gplt1))

#=
anim1 = @animate for i=1:1000
    push!(plt1,sol1[i])
end every 10
display(gif(anim1))

gplt1 = graphplot(sparse(LG))
display(gplt1)
=#
sol2 = simulateODEonGraph(GNPG, x, γ=0.1, β=0.1, dynamical_function=SIS_model)
plt2 = Plots.plot(sol2,xlim=(0,10),ylim=(0,1),title="gnp_graph SIS_model")
display(Plots.plot(plt2, gplt1))
sol3 = simulateODEonGraph(GNPG, x, β=0.1, dynamical_function=SI_model)
plt3 = Plots.plot(sol3,xlim=(0,10),ylim=(0,1),title="gnp_graph SI_model")
display(Plots.plot(plt3, gplt1))
sol4 = simulateODEonGraph(GNPG, x, K=0.1, dynamical_function=kuramoto_model)
plt4 = Plots.plot(sol4,xlim=(0,10),ylim=(0,1),title="gnp_graph Kuramoto_model")
display(Plots.plot(plt4, gplt1))
