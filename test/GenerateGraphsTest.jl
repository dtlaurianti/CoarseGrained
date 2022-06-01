@testset "line_graph_tests" begin
    @test line_graph(5) == [[0,1,0,0,0],[0,0,1,0,0][0,0,0,1,0][0,0,0,0,1][0,0,0,0,0]]
    @test line_graph(4; edge_weight=3) == [[0,3,0,0][0,0,3,0][0,0,0,3][0,0,0,0]]
    @test line_graph(3; directed=False) ==[[0,1,0],[1,0,1],[0,1,0]]
end
