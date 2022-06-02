@testset "line_graph_tests" begin
    @test broadcast(convert, Int, line_graph(5)) == [0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; 0 0 0 0 0]
    @test broadcast(convert, Int, line_graph(4; edge_weight=3)) == [0 3 0 0; 0 0 3 0; 0 0 0 3; 0 0 0 0]
    @test broadcast(convert, Int, line_graph(3; directed=false)) == [0 1 0; 1 0 1; 0 1 0]
    @test broadcast(convert, Int, line_graph(0; directed=false)) == []
    @test_throws DomainError broadcast(convert, Int, line_graph(-2; directed=false))
end

@testset "cycle_graph_tests" begin
    @test broadcast(convert, Int, cycle_graph(5)) == [0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; 1 0 0 0 0]
    @test broadcast(convert, Int, cycle_graph(4; edge_weight=3)) == [0 3 0 0; 0 0 3 0; 0 0 0 3; 3 0 0 0]
    @test broadcast(convert, Int, cycle_graph(3; directed=false)) == [0 1 1; 1 0 1; 1 1 0]
    @test broadcast(convert, Int, cycle_graph(0; directed=false)) == []
    @test_throws DomainError broadcast(convert, Int, cycle_graph(-2; directed=false))
end

@testset "grid_graph_tests" begin
    @test broadcast(convert, Int, grid_graph(4; directed=true)) == [0 1 1 0; 0 0 0 1; 0 0 0 1; 0 0 0 0]
    @test broadcast(convert, Int, grid_graph(4; edge_weight=3)) == [0 3 3 0; 3 0 0 3; 3 0 0 3; 0 3 3 0]
    @test broadcast(convert, Int, grid_graph(3)) == [0 1; 1 0]
    @test broadcast(convert, Int, grid_graph(0)) == []
    @test_throws DomainError broadcast(convert, Int, grid_graph(-2; directed=false))
end
