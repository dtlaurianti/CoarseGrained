#Function: getLayout
#Parameters: A, the network to generate a layout of
#            layout_func, the function to choose the location of the nodes
#Purpose: To create a layout for the given network
#Return value: two vectors x and y of the node coordinates
function getLayout(A::MatrixNetwork, layout_func=NetworkLayout.shell)
    if layout_func==NetworkLayout.shell
        avgd = sum(sparse(A))/n
        nodeds = [sum(sparse(A)[i,:]) for i=1:n]
        list_central = filter(x->nodeds[x] > avgd, 1:n)
        layout = NetworkLayout.shell(sparse(A), nlist=[list_central, ])
    else
        layout = layout_func(sparse(A))
    end
    x,y=Vector{Float64}(),Vector{Float64}()
    for point in layout
        append!(x, point[1])
        append!(y, point[2])
    end
    return x, y
end
#Function: getNodeSize
#Parameters: x, y, the coordinates of the nodes
#            scale, a scaling factor for bigger or smaller nodes
#Purpose: To create a layout for the given network
#Return value: two vectors x and y of the node coordinates
function getNodeSize(x, y; scale=0.5)
    for i=1:length(x)
        dist = sqrt(x[i]^2 + y[i]^2)
        if max_dist < dist
            max_dist = dist
        end
    end
    return max_dist*scale
end

#Function: plotNetwork
#Parameters: A, the network to plot
#            u, initial condition vector to scale the nodes
#            title, the plot title
#            layout_func, the function to choose the location of the nodes
#            colors, the colors of the nodes
#Purpose: To create a basic graph plot of the input network
#Return value: a graphplot
function plotNetwork(A::MatrixNetwork; u::Vector{Number}=nothing, title::String="", layout_func=NetworkLayout.shell, colors::Vector=nothing)
    if layout_func==NetworkLayout.shell
        avgd = sum(sparse(A))/n
        nodeds = [sum(sparse(A)[i,:]) for i=1:n]
        list_central = filter(x->nodeds[x] > avgd, 1:n)
        layout = NetworkLayout.shell(sparse(A), nlist=[list_central, ])
    else
        layout = layout_func(sparse(A))
    end

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

    if colors==nothing
        colors = Vector(1:n)
    end

    gplt = graphplot(
        sparse(A),
        x=x, y=y,
        markercolor=colors,
        nodesize=nodesize, node_weights=u,
        palette=distinguishable_colors(n),
        title=title)
    return gplt
end

#Parameters: A, the network to plot
#            x, y, coordinates of the nodes
#            max_dist, the distance of the furthest node from the center
#            u, initial condition vector to scale the nodes
#            title, the plot title
#            layout_func, the function to choose the location of the nodes
#            colors, the colors of the nodes
function plotNetwork(A::MatrixNetwork, x::Vector, y::Vector; nodesize::Number=0, u=nothing, title::String="", colors::Vector=nothing)
    if nodesize==0
        nodesize=getNodeSize(x, y)
    end

    if colors==nothing
        colors = Vector(1:n)
    end
    gplt = graphplot(
        sparse(A),
        x=x, y=y,
        markercolor = Vector(1:n),
        nodesize=nodesize, node_weights=u,
        palette=distinguishable_colors(n),
        title=title)
    return gplt
end

#Function: plotPartition
#Parameters: A, the network to plot
#            P, the partition to compress the network with
#            u, initial condition vector to scale the nodes
#            x, y, the coordinates of the nodes
#            nodesize, the size of nodes in the precompression graphplot
#            cnodesize, the size of nodes in the postcompression graphplot
#            colors, the colors of the nodes
#            title, the plot title
#Purpose: To create a plot of a network both before and after partitioning, showing the nodes that are combined into supernodes
#Return value: Two graphplots
function plotPartition(A::MatrixNetwork, P::Dict{Integer,Integer}, u::Vector, x::Vector, y::Vector; nodesize::Number=0, colors::Vector=nothing, title::String="")
    n = length(u)
    # create the post-compression network, variables, and layout
    CA = compressAdjacencyMatrix(A, P)
    cu = compressInitialCondition(u, P)
    cx, cy = compressNodeCoordinates(x, y, P)
    cn = length(cu)
    # compute the colors of the nodes in the original network to match the colors of nodes in the same supernode
    if colors==nothing
        colors = Vector(1:n)
    end
    colorInd = collect(values(sort(OrderedDict(P))))
    ccolors = [colors[colorInd[i]] for i=1:n]
    # compute the scale of the nodes in the network
    if nodesize==0
        nodesize=getNodeSize(x, y)
    end
    cnodesize=getNodeSize(x, y)
    # generate the precompression network graphplot
    gplt = graphplot(
        sparse(A),
        x=x, y=y,
        markercolor=ccolors,
        nodesize=nodesize, node_weights=u,
        palette=distinguishable_colors(n),
        title=title+" Precompression")
    # generate the postcompression network graphplot
    cgplt = graphplot(
        sparse(CA),
        x=cx, y=cy,
        markercolor=colors,
        nodesize=cnodesize, node_weights=cu,
        palette=distinguishable_colors(cn),
        title=title+" Postcompression")

    return gplt, cgplt
end



function plotPartitions(A::MatrixNetwork, partitions::Vector{Dict{Integer, Integer}}; titles::Vector{String}=nothing, layout_func::Function=NetworkLayout.shell)
    if titles==nothing
        titles = ["" for _=1:length(partitions)]
    end
    layout = getLayout(A, layout_func)
    return
end

#Function: plotDynamics
#Parameters: n, the number of variables
#            sol, the ODE solution
#            title, the plot title
#Purpose: To create a plot of the input network
#Return value: a graphplot
function plotDynamics(n, sol; title="")
    return Plots.plot(sol, title=title, palette = distinguishable_colors(n))
end

#Parameters: A, MatrixNetwork that represents the original network
#            partitions, a Vector of Dicts that specify the supernodes in the compressed networks
#            initial_condition, a Vector of the node variables in the original network
#            title, the title for the plot
#            dynamical_function, the method that we will use to calculate the dynamics on the network and it's compressed version
#            tmax, the final t value to compute up to
#            dt, the length of the time steps
#            function_args, a var-kwargs of the inputs to the model
function plotDynamics(A::MatrixNetwork, initial_condition::Vector; title="", dynamical_function::Function=linear_model, tmax::Number=10, dt::Number=0.01, function_args...)
    sol = simulateODEonGraph(A, initial_condition; dynamical_function=dynamical_function, tmax=tmax, dt=dt, function_args...)
    return Plots.plot(sol, title=title, palette = distinguishable_colors(n))
end
