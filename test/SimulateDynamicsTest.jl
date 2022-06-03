using Plots

LG = line_graph(5)
CG = cycle_graph(5)
GG = grid_graph(4; directed=true)

x = [1, 1, 1, 1, 1]
sol1 = simulateODEonGraph(LG, x)
sol2 = simulateODEonGraph(LG, x, γ=0.1, β=0.1, dynamical_function=SIS_model)
sol3 = simulateODEonGraph(LG, x, β=0.1, dynamical_function=SI_model)
sol4 = simulateODEonGraph(LG, x, dynamical_function=kuramoto_model)
