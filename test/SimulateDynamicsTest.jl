using Plots

LG = line_graph(5)
CG = cycle_graph(5)
GG = grid_graph(4; directed=true)

x = [1, 0, 0, 0, 0]

@testset "linear_model_tests" begin
    @test dump(simulateODEonGraph(LG, x)) == true
end
