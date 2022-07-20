using Profile, PProf, Plots
#=
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
=#

@testset "other_sbm_tests" begin
    display(Matrix(sparse(stochastic_block_model(150, 50, [500, 500]))))
    display(Matrix(sparse(stochastic_block_model(100, 30, [200, 200]))))
end
#=
@testset "cm_graph_tests" begin
    display((cm_graph(10, [1, 1, 2, 2, 3, 3, 4, 4, 5, 5])))
    display((cm_graph(10, [1, 1, 1, 1, 2, 2, 2, 3, 3, 4])))
    display((cm_graph(6, [1, 1, 1, 2, 2, 3])))
end
=#
#=
@testset "Efficiency Testing" begin
    #display(@benchmark line_graph(15))
    #@profile line_graph(15)
    #pprof(;webport=58699)
    #display(@benchmark cycle_graph(10000))
    #@profile cycle_graph(10000)
    #pprof(;webport=58699)
    #display(@benchmark grid_graph(5000))
    #@profile grid_graph(5000)
    #pprof(;webport=58699)
    #display(@benchmark gnp_graph(5000))
    #@profile gnp_graph(5000)
    #pprof(;webport=58699)
    #display(@benchmark sbm_graph(64, p_within = 1.0, p_between = 0.5))
    #@profile sbm_graph(64, p_within = 1.0, p_between = 0.5)
    #pprof(;webport=58699)
    #display(@benchmark cm_graph(25, [20, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10]))
    #@profile cm_graph(25, [20, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10])
    #pprof(;webport=58698)
end
=#