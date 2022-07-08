using JLD2

# network parameters
originalSize = 10
reducedSize  = 9
#=
options for graphType, graphArgs
    'random', {'p': p} --> generates G = GenerateGraphs.gnp_graph(originalSize, p)
    "cycle", {}:       --> generates G = GenerateGraphs.cycle_graph(originalSize, directed=False)
    "line", {}:        --> generates G = GenerateGraphs.line_graph(originalSize, directed=False)
=#
# graphType, graphArgs = "random", (p=0.4)
# graphType, graphArgs = "cycle", ()
graphType, graphArgs = "line", Dict()

# model parameters
#=
options for modelType, modelArgs
    "linear_model", {"epsilon": epsilon}        --> runs SimulateDynamics.linear_model
    "SIS_model", {"gamma": gamma, "beta": beta} --> runs SimulateDynamics.SIS
    "SI_model", {beta': beta}                   --> runs SimulateDynamics.SI
    "kuramoto_model", {"omega": omega, 'K': K}  --> runs SimulateDynamics.kuramoto_model
=#
# all models
# listModelType   = ["linear_model", "kuramoto_model", "SIS_model", "SI_model", "linear_opinions","nonlinear_opinions"]
# listModelAbbrev = ["Lin", "Kur", "SIS", "SI", "LinOp","NLOp"]
# listModelArgs   = [{"epsilon":-3/originalSize}, {"omega":1, "K":1}, {"beta":1, "gamma":1}, {"beta":1}, {"c":1}, {"d":0.1, "u":1, "b":0}]

# select models
listModelType   = ["linear_model"]
listModelAbbrev = ["Lin"]
listModelArgs   = [Dict(:ϵ=>-3/originalSize)]

# simulation parameters
tmax    = 10
tinc    = 0.1
numRuns = 1 # number of initial conditions to simulate

###### NO NEED TO TOUCH BELOW ######
# Generate a list of partitions to test
listOfPartitions  = kPartition(originalSize, reducedSize)

# Generate a list of initial conditions
listOfICs = [rand(originalSize) for run in 1:numRuns]

for id in 1:length(listModelType)
    # specify model parameters
    modelType = listModelType[id]
    modelAbbrev, modelArgs = listModelAbbrev[id], listModelArgs[id]

    # store parameters
    parameters = (originalSize=originalSize, reducedSize=reducedSize,
        graphType=graphType,graphArgs=graphArgs, modelType=modelType,
        modelArgs=modelArgs, modelAbbrev=modelAbbrev, tmax=tmax, tinc=tinc)

    # select model
    if modelType == "linear_model"
        modelFunc = linear_model
    elseif modelType == "SIS_model"
        modelFunc = SIS_model
    elseif modelType == "SI_model"
        modelFunc = SI_model
    elseif modelType == "kuramoto_model"
        modelFunc = kuramoto_model
    elseif modelType == "linear_opinions"
        modelFunc = linear_opinions
    elseif modelType == "nonlinear_opinions"
        modelFunc = nonlinear_opinions
    end


    # construct network
    isAccepted = false
    # declare A outside loop so it stays in scope
    A = nothing
    # generate a graph with the specified parameters, making sure it is connected
    while !isAccepted
        if graphType == "random"
            A = gnp_graph(originalSize, graphArgs...)
        elseif graphType == "cycle"
            A = cycle_graph(originalSize, directed=false)
        elseif graphType == "line"
            A = line_graph(originalSize, directed=false)
        end
        if isConnected(A)
            isAccepted = true
        end
    end
    # Run simulations numRuns times with different initial conditions
    for run in 1:numRuns
        # initial_condition = np.random.rand(originalSize)
        # initial_condition = np.concatenate((np.ones(1),np.zeros(originalSize-1)))
        initial_condition = listOfICs[run]
        # for each partition, create a task to solve for the loss, and pass the Future to the task to the losses RemoteChannel
        losses = getLossBatch(A, listOfPartitions, initial_condition, modelFunc, tmax, tinc; modelArgs...)

        # specify folder and check if it already exists
        modelArgsString = ""
        for (key,value) in modelArgs
            modelArgsString *= "_"*string(key)*"="*string(value)
        end
        graphArgsString = ""
        for (key,value) in graphArgs
            graphArgsString *= "_"*string(key)*"="*string(value)
        end

        # generate a place to store the file of the data
        foldername = "data/model_data/"*modelType*modelArgsString*"/"*graphType*graphArgsString
        if !isdir(foldername)
             mkpath(foldername) # make folder if it doesn't exist
        end
        # prep filename and data
        filename = foldername * "/" * string(originalSize) * "_" * string(reducedSize) * "_run" * string(run)

        # save all the data in a .jld2 file which we can later import
        # to recreate the data in their original Julia Object forms
        jldsave("$filename.jld2";
            A,
            partitions = listOfPartitions,
            losses = losses,
            parameters,
            initial_condition
        )

        #=
        data = Dict()
        data["A"] = A
        data["partitions"] = listOfPartitions
        # get all the loss data from the RemoteChannel for this run
        data["losses"] = [fetch(losses[i]) for i=1:length(listOfPartitions)]
        data["parameters"] = parameters
        data["initial_condition"] = initial_condition
        store("./data/model_data/testid$(id)run$(run).pkl", dump(data))

        #open(filename, "wb") do file
        #    pickle.dump(data, file)
        #end
        =#
    end
end
