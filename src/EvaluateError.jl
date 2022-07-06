
#Function: aggregateTimeSeries
#Parameters: originalTimeSeries, the ODESolution for the uncompressed network
#            partition, a Dict specifying the compression of the network
#            type, a keyword choosing the method to compute the aggregate time series
#Purpose: To take a time series from an uncompressed network and compress that series to the size of the supernodes for comparison with a compressed network
#Return value: Time series with the compressed version of originalTimeSeries
@everywhere function aggregateTimeSeries(originalTimeSeries, partition::Dict{Integer,Integer}, type="average")
    if type == "average"
        supernodeSizes = getSupernodeSizes(partition)
        reducedSize = length(supernodeSizes)
        originalSize = size(originalTimeSeries,1)
        numTimeSteps = size(originalTimeSeries,2)
        aggregatedTimeSeries = zeros((reducedSize, numTimeSteps))
        for node in 1:originalSize
          aggregatedTimeSeries[partition[node], :] .+= originalTimeSeries[node, :]./supernodeSizes[partition[node]]
        end
    end
    return aggregatedTimeSeries
end

#Function: computeDynamicalError
#Parameters: originalTimeSeries, the ODESolution for the uncompressed network
#            reducedTimeSeries, the ODESolution for the compressed network
#            partition, a Dict specifying the supernodes of the compressed network
#Purpose: To compute the difference in the time series for a network and it's compressed version
#Return value: Float of loss value
@everywhere function computeDynamicalError(originalTimeSeries, reducedTimeSeries, partition::Dict{Integer,Integer})
  numTimeSteps = size(reducedTimeSeries, 2)
  aggregatedTimeSeries = aggregateTimeSeries(originalTimeSeries, partition)
  return sum(sum((reducedTimeSeries - aggregatedTimeSeries).^2))/numTimeSteps
end

#Function: computeIndividualError
#Parameters: originalTimeSeries, the ODESolution for the uncompressed network
#            reducedTimeSeries, the ODESolution for the compressed network
#            partition, a Dict specifying the supernodes of the compressed network
#Purpose: To compute the difference in the time series for a network and it's compressed version
#Return value: Float of loss value
@everywhere function computeIndividualError(originalTimeSeries, reducedTimeSeries, partition::Dict{Integer,Integer})
  originalSize = size(reducedTimeSeries, 1)
  loss = 0
  for node in 1:originalSize
    loss += lossFunction(originalTimeSeries[node, :], reducedTimeSeries[partition[node], :])
  end
  return loss
end

#Function: getLoss
#Parameters: timeseries1, an ODESolution
#            timeseries2, an ODESolution
#            type, a keyword for the way that we will calculate the loss
#Purpose: To compute the difference in two time series using the specified method
#Return value: Float of loss value
@everywhere function lossFunction(timeseries1, timeseries2, type="L2")
  if type == "L2"
    return sum((timeseries1 .- timeseries2).^2)/size(timeseries1, 1)
  end
end


#Function: getLoss
#Parameters: A, MatrixNetwork that represents the original network
#            partition, a Dict that specifies the supernodes in the compressed network
#            initial_condition, a Vector of the node variables in the original network
#            dynamical_function, the method that we will use to calculate the dynamics on the network and it's compressed version
#            tmax, the final t value to compute up to
#            dt, the length of the time steps
#            function_args, a var-kwargs of the inputs to the model
#Purpose: To compute the loss caused by coarse-graining a network according to a partition
#Return value: Float of loss value
@everywhere function getLoss(A::MatrixNetwork, partition::Dict{Integer, Integer}, initial_condition::Vector, dynamical_function::Function, tmax::Number, dt::Number; function_args...)
  function_args = Dict(function_args)
  # create the initial conditions for the compressed version of A
  compressed_initial_condition = compressInitialCondition(initial_condition, partition)
  # create the model arguments for the compressed version of A
  compressed_args = compressArguments(partition, function_args...)
  # create the compressed version of A
  reducedA = compressAdjacencyMatrix(A, partition)
  # compute the time series
  originalTimeSeries = simulateODEonGraph(A, initial_condition; dynamical_function=dynamical_function, tmax=tmax, dt=dt, function_args...)
  reducedTimeSeries = simulateODEonGraph(reducedA, compressed_initial_condition; dynamical_function=dynamical_function, tmax=tmax, dt=dt, compressed_args...)

  # calculate the difference between the time series
  #loss = computeIndividualError(originalTimeSeries, reducedTimeSeries, partition)
  loss = computeDynamicalError(originalTimeSeries, reducedTimeSeries, partition)
  return loss
end

#Function: getLossPartition
#Parameters: originalTimeSeries, the time series of the uncompressed network that we precompute
#            A, MatrixNetwork that represents the original network
#            partition, a Dict that specifies the supernodes in the compressed network
#            initial_condition, a Vector of the node variables in the original network
#            dynamical_function, the method that we will use to calculate the dynamics on the network and it's compressed version
#            tmax, the final t value to compute up to
#            dt, the length of the time steps
#            function_args, a var-kwargs of the inputs to the model
#Purpose: To compute the loss caused by coarse-graining a network according to a partition
#Return value: Float of loss value
@everywhere function getLossPartition(originalTimeSeries, A::MatrixNetwork, partition::Dict{Integer, Integer}, initial_condition::Vector, dynamical_function::Function, tmax::Number, dt::Number; loss_function::Function=computeDynamicalError, function_args...)
  function_args = Dict(function_args)
  # create the initial conditions for the compressed version of A
  compressed_initial_condition = compressInitialCondition(initial_condition, partition)
  # create the model arguments for the compressed version of A
  compressed_args = compressArguments(partition, function_args...)
  # create the compressed version of A
  reducedA = compressAdjacencyMatrix(A, partition)
  # compute the time series
  reducedTimeSeries = simulateODEonGraph(reducedA, compressed_initial_condition; dynamical_function=dynamical_function, tmax=tmax, dt=dt, compressed_args...)

  # calculate the difference between the time series
  #loss = computeIndividualError(originalTimeSeries, reducedTimeSeries, partition)
  loss = loss_function(originalTimeSeries, reducedTimeSeries, partition)
  return loss
end

#Function: getLossBatch
#Parameters: A, MatrixNetwork that represents the original network
#            partitions, a Vector of Dicts that specify the supernodes in the compressed networks
#            initial_condition, a Vector of the node variables in the original network
#            dynamical_function, the method that we will use to calculate the dynamics on the network and it's compressed version
#            tmax, the final t value to compute up to
#            dt, the length of the time steps
#            function_args, a var-kwargs of the inputs to the model
#Purpose: To compute the loss caused by coarse-graining a network according to a partition
#Return value: Vector of Float of loss value
@everywhere function getLossBatch(A::MatrixNetwork, partitions::Vector{Dict{Integer, Integer}}, initial_condition::Vector, dynamical_function::Function, tmax::Number, dt::Number; loss_function::Function=computeDynamicalError, function_args...)
  originalTimeSeries = simulateODEonGraph(A, initial_condition; dynamical_function=dynamical_function, tmax=tmax, dt=dt, function_args...)
  return pmap(part->getLossPartition(originalTimeSeries, A, part, initial_condition, dynamical_function, tmax, dt; loss_function=loss_function, function_args...), partitions)
end
