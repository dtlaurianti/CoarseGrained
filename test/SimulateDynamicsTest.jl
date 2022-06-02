LG = line_graph(5)
CG = cycle_graph(5)
GG = grid_graph(4; directed=true)

linear_model(5, 5, LG, 5)

@testset "linear_model_tests" begin
    @test true
