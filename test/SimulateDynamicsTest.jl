#=
n = 10

LG = line_graph(n)
CG = cycle_graph(n)
GG = grid_graph(n; directed=true)
GNPG = gnp_graph(n, p=0.5)
SBMG = sbm_graph(n)

x = rand(n)


scale = 1/log(n)

sol1 = simulateODEonGraph(CG, x)
plt1 = Plots.plot(sol1,title="cycle_graph linear_model", palette=distinguishable_colors(n))
gplt1 = graphplot(sparse(CG), thickness_scaling=scale, markercolor = Vector(1:n), node_weights=x, palette=distinguishable_colors(n), names = 1:n)
display(Plots.plot(plt1, gplt1))

display(plt1)

plta = Plots.plot(n, xlim=(0,1000), ylim=(0,1), legend=false, title="cycle_graph linear_model", palette=distinguishable_colors(n))
anim1 = @animate for i=1:1000
    push!(plta,sol1[i])
end every 10
display(gif(anim1))
=#
#=
gplt1 = graphplot(sparse(LG))
display(gplt1)


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
sol5 = simulateODEonGraph(SBMG, x, ω=rand(n), dynamical_function=LotkaVolterra_model)
plt5 = Plots.plot(sol5,title="sbm_graph LotkaVolterra_model", palette=distinguishable_colors(n))
display(Plots.plot(plt5, gplt5))

sol6 = simulateODEonGraph(GNPG, x, c=.1, dynamical_function=linear_opinions)
plt6 = Plots.plot(sol6,title="gnp_graph linear_opinions_model", palette=distinguishable_colors(n))
display(Plots.plot(plt6, gplt4))

sol7 = simulateODEonGraph(GNPG, x, d=.1, c=.1, b=.1, dynamical_function=nonlinear_opinions)
plt7 = Plots.plot(sol7,title="gnp_graph nonlinear_opinions_model", palette=distinguishable_colors(n))
display(Plots.plot(plt7, gplt4))
=#

@testset "efficiency_testing" begin
    n=5
    G = gnp_graph(n)
    x = rand(n)
    @time simulateODEonGraph(G, x)
    @time simulateODEonGraphMap(G, x)
    #@time simulateODEonGraphStatic(G, x)
end
