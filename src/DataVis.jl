#Function: getLayout
#Parameters: A, the network to generate a layout of
#            layout_func, the function to choose the location of the nodes
#Purpose: To create a layout for the given network
#Return value: two vectors x and y of the node coordinates
function getLayout(A::MatrixNetwork, layout_func=NetworkLayout.shell)
    n = size(sparse(A), 1)
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
    max_dist = 0
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
function plotNetwork(A::MatrixNetwork; u=nothing, title::String="", layout_func=NetworkLayout.shell, colors::Vector=Vector())
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

    if colors==Vector()
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
function plotNetwork(A::MatrixNetwork, x::Vector, y::Vector; nodesize::Number=0, u=nothing, title::String="", colors::Vector=Vector())
    if nodesize==0
        nodesize=getNodeSize(x, y)
    end

    if colors==Vector()
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
function plotPartition(A::MatrixNetwork, P::Dict{Integer,Integer}, u::Vector, x::Vector, y::Vector; nodesize::Number=0, colors::Vector=Vector(), title::String="")
    n = length(u)
    scale = 1/log(n, 10)
    # create the post-compression network, variables, and layout
    CA = compressAdjacencyMatrix(A, P)
    cu = compressInitialCondition(u, P)
    cx, cy = compressNodeCoordinates(x, y, P)
    cn = length(cu)
    cscale = 1/log(cn, 10)
    # compute the colors of the nodes in the original network to match the colors of nodes in the same supernode
    if colors==Vector()
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
        thickness_scaling=scale,
        title=title*" Pre")
    # generate the postcompression network graphplot
    cgplt = graphplot(
        sparse(CA),
        x=cx, y=cy,
        markercolor=colors,
        nodesize=cnodesize, node_weights=cu,
        palette=distinguishable_colors(cn),
        thickness_scaling=scale,
        title=title*" Post")

    return gplt, cgplt
end


#Function: plotPartition
#Parameters: A, the network to plot
#            partitions, the partitions to compress the network with
#            u, initial condition vector to scale the nodes
#            titles, the plot titles
#            layout_func, the function to compute the layout with
#Purpose: To create a plot of a network both before and after partitioning, showing the nodes that are combined into supernodes
#Return value: a list of graphplots of before and after partitioning for each partition
function plotPartitions(A::MatrixNetwork, partitions::Vector{Dict{Integer, Integer}}, u::Vector; titles::Vector{String}=Vector{String}(), layout_func::Function=NetworkLayout.shell)
    if titles==Vector{String}()
        titles = ["" for _=1:length(partitions)]
    end
    x,y = getLayout(A, layout_func)
    if titles==Vector{String}()
        for i = 1:length(partitions)
            push!(titles, "P$i")
        end
    end
    gplts = []
    cgplts = []
    titles = copy(titles)
    for partition in partitions
        gplt, cgplt = plotPartition(A, partition, u, x, y, title=popfirst!(titles))
        push!(gplts, gplt)
        push!(cgplts, cgplt)
    end
    return gplts, cgplts
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
    function_args = Dict(function_args)
    n = size(A, 1)
    sol = simulateODEonGraph(A, initial_condition; dynamical_function=dynamical_function, tmax=tmax, dt=dt, function_args...)
    return Plots.plot(sol, title=title, palette = distinguishable_colors(n))
end

#Function: plotPartitionsDynamics
#Parameters: A, the network to plot
#            partitions, the partitions to compress the network with
#            u, the variable vector to scale the nodes
#            titles, the plot titles
#            layout_func, the function to compute the layout with
#            dynamical_function, the method that we will use to calculate the dynamics on the network and it's compressed version
#            tmax, the final t value to compute up to
#            dt, the length of the time steps
#            function_args, a var-kwargs of the inputs to the model
#Purpose: To create a plot of a network both before and after partitioning, showing the nodes that are combined into supernodes, as well as the dynamics of the network
#Return value: a list of graphplots of before and after partitioning for each partition
function plotPartitionsDynamics(A::MatrixNetwork, partitions::Vector{Dict{Integer, Integer}}, u::Vector; titles::Vector{String}=Vector{String}(), layout_func::Function=NetworkLayout.shell, dynamical_function::Function=linear_model, tmax::Number=10, dt::Number=0.01, function_args...)
    if titles==Vector{String}()
        titles = ["P$i" for i=1:length(partitions)]
    end
    titles = copy(titles)
    gplts, cgplts = plotPartitions(A, partitions, u, titles=titles, layout_func=layout_func)
    dplts = []
    for partition in partitions
        CA = compressAdjacencyMatrix(A, partition)
        cu = compressInitialCondition(u, partition)
        function_args = compressArguments(partition, function_args...)
        push!(dplts, plotDynamics(CA, cu, title=popfirst!(titles)*" Dyn", dynamical_function=dynamical_function, tmax=tmax, dt=dt, function_args...))
    end
    return gplts, cgplts, dplts
end

#Function: animatePartitionsDynamics
#Parameters: A, the network to plot
#            partition, the partition to compress the network with
#            u, the variable vector to scale the nodes
#            titles, the plot titles
#            layout_func, the function to compute the layout with
#            dynamical_function, the method that we will use to calculate the dynamics on the network and it's compressed version
#            tmax, the final t value to compute up to
#            dt, the length of the time steps
#            function_args, a var-kwargs of the inputs to the model
#Purpose: To create an animation of a network both before and after partitioning, showing the nodes that are combined into supernodes, as well as the dynamics of the network
#Return value: a gif showing the three plots of the precompression network, postcompression network, and dynamics
function animatePartitionDynamics(A::MatrixNetwork, partition::Dict{Integer, Integer}, u::Vector; title::String="", layout_func::Function=NetworkLayout.shell, dynamical_function::Function=linear_model, tmax::Number=10, dt::Number=0.01, function_args...)
    x, y = getLayout(A, layout_func)
    n = size(A, 1)
    CA = compressAdjacencyMatrix(A, partition)
    cn = size(CA,1)
    cu = compressInitialCondition(u, partition)
    cfunction_args = compressArguments(partition, function_args...)
    sol = simulateODEonGraph(A, u; dynamical_function=dynamical_function, tmax=tmax, dt=dt, function_args...)
    csol = simulateODEonGraph(CA, cu; dynamical_function=dynamical_function, tmax=tmax, dt=dt, cfunction_args...)
    tsteps = trunc(Int64, tmax/dt)
    ymax = maximum(maximum(sol.u))
    cymax = maximum(maximum(csol.u))
    scale = 1/(2*log(cn, 10))
    palette = distinguishable_colors(cn)
    colorInd = collect(values(sort(OrderedDict(partition))))
    cpalette = [palette[colorInd[i]] for i=1:n]
    anim = @animate for i=1:tsteps
        gplt, cgplt = plotPartition(A, partition, sol[i], x, y, title=title)
        subplt = Plots. plot(gplt, cgplt, layout=(2,1))
        dplt = Plots.plot(sol[1:i], xlim=(0, tmax), ylim=(0, ymax), title=title*" Pre Dyn", palette=cpalette, leg=false)
        cdplt = Plots.plot(csol[1:i], xlim=(0, tmax), ylim=(0, ymax), title=title*" Post Dyn", palette=palette, leg=false)
        plt = Plots.plot(gplt, dplt, cgplt, cdplt, layout=(2,2), thickness_scaling=scale)
    end every 10
    return gif(anim)
end

function animatePartitionsDynamics(A::MatrixNetwork, partitions::Vector{Dict{Integer, Integer}}, u::Vector; title::String="", layout_func::Function=NetworkLayout.shell, dynamical_function::Function=linear_model, tmax::Number=10, dt::Number=0.01, function_args...)
    anims = []
    for partition in partitions
        push!(anims, animatePartitionDynamics(A, partition, u, title=title, layout_func=layout_func, dynamical_function=dynamical_function, tmax=tmax, dt=dt, function_args...))
    end
    return anims
end
