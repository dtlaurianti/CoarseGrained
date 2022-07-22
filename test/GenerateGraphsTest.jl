#using Profile, PProf, Plots
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
#=
@testset "other_sbm_tests" begin
    display(Matrix(sparse(stochastic_block_model(150, 50, [500, 500]))))
    display(Matrix(sparse(stochastic_block_model(100, 30, [200, 200]))))
end
=#
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

@testset "layout_tests" begin
    n = 10
    c = trunc(Int64, n/2)
    #scale = 1/log(n, 10)
    G = gnp_graph(n, p=0.25, directed=false)
    #display(sparse(G))
    u = rand(n)

    avgd = sum(sparse(G))/n
    nodeds = [sum(sparse(G)[i,:]) for i=1:n]
    list_central = filter(x->nodeds[x] > avgd, 1:n)
    layout = NetworkLayout.shell(sparse(G), nlist=[list_central, ])
    #display(layout)

    x,y=Vector{Float64}(),Vector{Float64}()
    max_dist = 0
    for point in layout
        append!(x, point[1])
        append!(y, point[2])
        dist = sqrt(point[1]^2 + point[2]^2)
        if max_dist < dist
            max_dist = dist
        end
    end
    nodesize = max_dist/2


    gplt = graphplot(sparse(G), x=x, y=y, markercolor = Vector(1:n), nodesize=nodesize, node_weights=u, palette=distinguishable_colors(n), title="Original Network")

    Ps = generateRandomPartitions(n,c,2)

    G1 = compressAdjacencyMatrix(G, Ps[1])
    colors1 = collect(values(sort(Ps[1])))
    gplt1 = graphplot(sparse(G), x=x, y=y, markercolor = colors1, nodesize=nodesize, node_weights=u, palette=distinguishable_colors(n), title="Random Compression")

    G2 = compressAdjacencyMatrix(G, Ps[2])
    colors2 = collect(values(sort(Ps[2])))
    gplt2 = graphplot(sparse(G), x=x, y=y, markercolor = colors2, nodesize=nodesize, node_weights=u, palette=distinguishable_colors(n), title="Random Compression")

    Pa = agglomerationReduction(G, c)

    G3 = compressAdjacencyMatrix(G, Pa)
    colors3 = collect(values(sort(Pa)))
    gplt3 = graphplot(sparse(G), x=x, y=y, markercolor = colors3, nodesize=nodesize, node_weights=u, palette=distinguishable_colors(n), title = "Agglomeration Reduction")

    plt = Plots.plot(gplt, gplt3, gplt2, gplt1)
    display(plt)
end

#=
@testset "animation_tests" begin
    n = 16
    cn = trunc(Int64, n/2)
    #scale = 1/log(n, 10)
    G = gnp_graph(n, p=0.25, directed=false)
    #G = grid_graph(n, directed=false)
    #display(sparse(G))
    u = rand(n)

    avgd = sum(sparse(G))/n
    nodeds = [sum(sparse(G)[i,:]) for i=1:n]
    list_central = filter(x->nodeds[x] > avgd, 1:n)
    layout = NetworkLayout.shell(sparse(G), nlist=[list_central, ])
    #display(layout)

    x,y=Vector{Float64}(),Vector{Float64}()
    max_dist = 0
    for point in layout
        append!(x, point[1])
        append!(y, point[2])
        dist = sqrt(point[1]^2 + point[2]^2)
        if max_dist < dist
            max_dist = dist
        end
    end
    nodesize = max_dist/2

    sol = simulateODEonGraph(G, u)
    #display(sol)
    Pa = agglomerationReduction(G, cn)
    CG = compressAdjacencyMatrix(G, Pa)
    cu = compressInitialCondition(u, Pa)
    csol = simulateODEonGraph(CG, cu)
    colors = collect(values(sort(Pa)))


    gplt = graphplot(sparse(G), x=x, y=y, markercolor = colors, nodesize=nodesize, node_weights=u, palette=distinguishable_colors(n), title="Original Network")


    cx,cy = compressNodeCoordinates(x,y,Pa)

    cgplt = graphplot(sparse(CG), x=cx, y=cy, markercolor = Vector(1:cn), nodesize=nodesize, node_weights=cu, palette=distinguishable_colors(cn), title="Agglomeration Reduction")

    anim = @animate for i=1:250
        gplt = graphplot(sparse(G), x=x, y=y, markercolor = colors, nodesize=nodesize, node_weights=sol[i], palette=distinguishable_colors(n), title="Original Network")
        cgplt = graphplot(sparse(CG), x=cx, y=cy, markercolor = Vector(1:cn), nodesize=nodesize, node_weights=csol[i], palette=distinguishable_colors(cn), title="Agglomeration Reduction")
        plt = Plots.plot(gplt, cgplt)
    end every 5
    display(gif(anim))
end
=#
#idea to prevent scaling of nodes all increasing together, something like this
#log(mean(sol[i]./u))
