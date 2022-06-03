using Plots

LG = line_graph(5)
CG = cycle_graph(5)
GG = grid_graph(4; directed=true)

x = [1, 1, 1, 1, 1]
plot(simulateODEonGraph(LG, x))
@testset "linear_model_tests" begin
    @test plot(simulateODEonGraph(LG, x) == true
end
