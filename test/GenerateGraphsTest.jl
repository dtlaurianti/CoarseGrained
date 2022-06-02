@testset "line_graph_tests" begin
    @test line_graph(5) == MatrixNetwork{Float64}(5, [1, 2, 3, 4, 5, 5], [2, 3, 4, 5], [1.0, 1.0, 1.0, 1.0])
    @test line_graph(4; edge_weight=3) == true
    @test line_graph(3; directed=false) == true
    @test line_graph(0; directed=false) == true
    @test_throws DomainError line_graph(-2; directed=false) == true
end

@testset "cycle_graph_tests" begin
    @test cycle_graph(5) == true
    @test cycle_graph(4; edge_weight=3) == true
    @test cycle_graph(3; directed=false) == true
    @test cycle_graph(0; directed=false) == true
    @test_throws DomainError cycle_graph(-2; directed=false) == true
end

@testset "grid_graph_tests" begin
    @test grid_graph(4; directed=true) == true
    @test grid_graph(4; edge_weight=3) == true
    @test grid_graph(3) == true
    @test grid_graph(0) == true
    @test_throws DomainError grid_graph(-2; directed=false) == true
end
