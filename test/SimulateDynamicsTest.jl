LG = line_graph(5)
CG = cycle_graph(5)
GG = grid_graph(4; directed=true)

x = [1, 2, 3, 4, 5]

linear_model(5, x, LG, 5)

@testset "linear_model_tests" begin
    @test true
