@testset "line_graph_tests" begin
    @test norm(sparse(line_graph(5))) - norm(sparse(MatrixNetwork{Float64}(5, [1, 2, 3, 4, 5, 5], [2, 3, 4, 5], [1.0, 1.0, 1.0, 1.0]))) ≈ 0
    @test norm(sparse(line_graph(4; edge_weight=3))) - norm(sparse(MatrixNetwork{Float64}(4, [1, 2, 3, 4, 4], [2, 3, 4], [3.0, 3.0, 3.0]))) ≈ 0
    @test norm(sparse(line_graph(3; directed=false))) - norm(sparse(MatrixNetwork{Float64}(3, [1, 2, 4, 5], [2, 1, 3, 2], [1.0, 1.0, 1.0, 1.0]))) ≈ 0
    @test norm(sparse(line_graph(0; directed=false))) - norm(sparse(MatrixNetwork{Float64}(0, [1], Int64[], Float64[]))) ≈ 0
    @test_throws DomainError line_graph(-2; directed=false) == true
end

@testset "cycle_graph_tests" begin
    @test norm(sparse(cycle_graph(5))) - norm(sparse(MatrixNetwork{Float64}(5, [1, 2, 3, 4, 5, 6], [2, 3, 4, 5, 1], [1.0, 1.0, 1.0, 1.0, 1.0]))) ≈ 0
    @test norm(sparse(cycle_graph(4; edge_weight=3))) - norm(sparse(MatrixNetwork{Float64}(4, [1, 2, 3, 4, 5], [2, 3, 4, 1], [3.0, 3.0, 3.0, 3.0]))) ≈ 0
    @test norm(sparse(cycle_graph(3; directed=false))) - norm(sparse(MatrixNetwork{Float64}(3, [1, 3, 5, 7], [2, 3, 1, 3, 1, 2], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]))) ≈ 0
    @test norm(sparse(cycle_graph(0; directed=false))) - norm(sparse(MatrixNetwork{Float64}(0, [1], Int64[], Float64[]))) ≈ 0
    @test_throws DomainError cycle_graph(-2; directed=false)
end

@testset "grid_graph_tests" begin
    @test norm(sparse(grid_graph(5))) - norm(sparse(MatrixNetwork{Float64}(6, [1, 3, 5, 8, 11, 13, 15], [2, 3, 1, 4, 1, 4, 5, 2, 3, 6, 3, 6, 4, 5], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]))) ≈ 0
    @test norm(sparse(grid_graph(4; edge_weight=3))) - norm(sparse(MatrixNetwork{Float64}(4, [1, 3, 5, 7, 9], [2, 3, 1, 4, 1, 4, 2, 3], [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0]))) ≈ 0
    @test norm(sparse(grid_graph(3; directed=false))) - norm(sparse(MatrixNetwork{Float64}(2, [1, 2, 3], [2, 1], [1.0, 1.0]))) ≈ 0
    @test norm(sparse(grid_graph(0; directed=false))) - norm(sparse(MatrixNetwork{Float64}(0, [1], Int64[], Float64[]))) ≈ 0
    @test_throws DomainError grid_graph(-2; directed=false) == true
end

@testset "gnp_graph_tests" begin
    display(sparse(gnp_graph(5; edge_weight=3)))
    display(sparse(gnp_graph(5;p=0.5)))
    display(sparse(gnp_graph(5;directed=false)))
end

@testset "sbm_graph_tests" begin
    display(Matrix(sparse(sbm_graph(64))))
    display(Matrix(sparse(sbm_graph(64, p_within = 1.0, p_between = 0.5))))
    display(Matrix(sparse(sbm_graph(27, communities = 2, p_within = 1.0))))
    display(Matrix(sparse(sbm_graph(27, communities = 1))))
end

@testset "cm_graph_tests" begin
    display((cm_graph(10, 5)))
    display((cm_graph(20, 4)))
    display((cm_graph(28, 7)))
end