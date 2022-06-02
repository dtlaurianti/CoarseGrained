LG = line_graph(5)
CG = cycle_graph(5)
GG = grid_graph(4; directed=true)

x = [1, 1, 1, 1, 1]

linear_model(5, x, LG, 5)

@testset "linear_model_tests" begin
    @test true
