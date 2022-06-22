
n = 16

LG = line_graph(n)
CG = cycle_graph(n)
GG = grid_graph(n; directed=true)
GNPG = gnp_graph(n, p=0.5)
SBMG = sbm_graph(n)

x = rand(n)


scale = 1/log(n)

sol1 = simulateODEonGraph(LG, x)
plt1 = Plots.plot(sol1,title="line_graph linear_model", palette=distinguishable_colors(n))
gplt1 = graphplot(sparse(LG), thickness_scaling=scale, markercolor = Vector(1:n), node_weights=x, palette=distinguishable_colors(n), names = 1:n)
display(Plots.plot(plt1, gplt1))

#=
anim1 = @animate for i=1:1000
    push!(plt1,sol1[i])
end every 10
display(gif(anim1))

gplt1 = graphplot(sparse(LG))
display(gplt1)
=#
gplt2 = graphplot(sparse(CG), thickness_scaling=scale, markercolor = Vector(1:n), node_weights=x, palette=distinguishable_colors(n))
sol2 = simulateODEonGraph(CG, x, γ=0.1, β=0.1, dynamical_function=SIS_model)
plt2 = Plots.plot(sol2,title="cycle_graph SIS_model", palette=distinguishable_colors(n))
display(Plots.plot(plt2, gplt2))

gplt3 = graphplot(sparse(GG), thickness_scaling=scale, markercolor = Vector(1:n), node_weights=x, palette=distinguishable_colors(n))
sol3 = simulateODEonGraph(GG, x, β=0.1, dynamical_function=SI_model)
plt3 = Plots.plot(sol3,title="grid_graph SI_model", palette=distinguishable_colors(n))
display(Plots.plot(plt3, gplt3))

gplt4 = graphplot(sparse(GNPG), thickness_scaling=scale, markercolor = Vector(1:n), node_weights=x, palette=distinguishable_colors(n))
sol4 = simulateODEonGraph(GNPG, x, K=0.1, ω=rand(n), dynamical_function=kuramoto_model)
plt4 = Plots.plot(sol4,title="gnp_graph Kuramoto_model", palette=distinguishable_colors(n))
display(Plots.plot(plt4, gplt4))

gplt5 = graphplot(sparse(SBMG), thickness_scaling=scale, markercolor = Vector(1:n), node_weights=x, palette=distinguishable_colors(n))
sol5 = simulateODEonGraph(SBMG, x, K=0.1, ω=rand(n), dynamical_function=LotkaVolterra_model)
plt5 = Plots.plot(sol5,title="sbm_graph LotkaVolterra_model", palette=distinguishable_colors(n))
display(Plots.plot(pl54, gplt5))
