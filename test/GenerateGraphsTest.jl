@testset "line_graph_tests" begin
    @test broadcast(convert, Int, line_graph(5)) == [0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; 0 0 0 0 0]
    @test broadcast(convert, Int, line_graph(4; edge_weight=3)) == [0 3 0 0; 0 0 3 0; 0 0 0 3; 0 0 0 0]
    @test broadcast(convert, Int, line_graph(3; directed=false)) ==[0 1 0; 1 0 1; 0 1 0]
end
