using Plots
using GraphRecipes

n = 16

LG = line_graph(n)
CG = cycle_graph(n)
GG = grid_graph(n; directed=true)
GNPG = gnp_graph(n, p=0.5)

x = rand(n)


scale = 1/log(n)

sol1 = simulateODEonGraph(CG, x)
plt1 = Plots.plot(sol1,title="cycle_graph linear_model", palette=distinguishable_colors(n))
gplt1 = graphplot(sparse(CG), thickness_scaling=scale, markercolor = Vector(1:n), node_weights=x, palette=distinguishable_colors(n), names = 1:n)
display(Plots.plot(plt1, gplt1))

#=
anim1 = @animate for i=1:1000
    push!(plt1,sol1[i])
end every 10
display(gif(anim1))

gplt1 = graphplot(sparse(LG))
display(gplt1)
=#
gplt2 = graphplot(sparse(GNPG), thickness_scaling=scale, markercolor = Vector(1:n), node_weights=x, palette=distinguishable_colors(n))
sol2 = simulateODEonGraph(GNPG, x, γ=0.1, β=0.1, dynamical_function=SIS_model)
plt2 = Plots.plot(sol2,title="gnp_graph SIS_model", palette=distinguishable_colors(n))
display(Plots.plot(plt2, gplt2))

sol3 = simulateODEonGraph(GNPG, x, β=0.1, dynamical_function=SI_model)
plt3 = Plots.plot(sol3,title="gnp_graph SI_model", palette=distinguishable_colors(n))
display(Plots.plot(plt3, gplt2))

sol4 = simulateODEonGraph(CG, x, K=0.1, ω=rand(n), dynamical_function=kuramoto_model)
plt4 = Plots.plot(sol4,title="gnp_graph Kuramoto_model", palette=distinguishable_colors(n))
display(Plots.plot(plt4, gplt1))
