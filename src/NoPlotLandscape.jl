using Clustering
using MultivariateStats
using DataFrames
using CSV
using ScikitLearn
using StatsBase
using SharedArrays

function GetXYZ(partitions::Vector{Dict{Integer, Integer}}, A, NumOriginalNodes; save_to_string="", modelType::Function=linear_model)
    listModelArgs = Dict(:ϵ=>-3/NumOriginalNodes, :β=>0.5, :γ=>0.5, :ω=>rand(NumOriginalNodes), :K=>0.5, :d=>0.5, :c=>0.5, :b=>0.5)
    #convert dictionary to an array
    Arr = dict_to_array(partitions)

    #calculate distance matrix
    num_par = length(partitions)
    D = SharedArray{Float64}((num_par,num_par))
    @sync @distributed for i in 1:num_par
        for j in 1:num_par
            # the dissimilarity matrix
            D[i, j] = variation_of_information(Arr[i],Arr[j])
        end
    end

    D = Array(D)
    #calculate MDS on disimilarity matrix
    embedding = StatsBase.fit(MDS, D, distances=true, maxoutdim=2)
    X_transformed = StatsBase.predict(embedding)
    #Format data
    x = X_transformed[1,:]
    y = X_transformed[2,:]

    #Calculate z dimension
    z = zeros(num_par)
    for i in 1:num_par
        # using hard-coded model and parameters, possibly want to make the outer function accept those parameters?
      loss = getLoss(A, partitions[i], rand(NumOriginalNodes), modelType, NumOriginalNodes, 0.01; listModelArgs...)
      z[i] = loss
    end

    #Save vector data if we want to smooth it later
    if !isempty(save_to_string)
      df = DataFrame(["x" => x, "y" => y, "z" => z])
      loc = "./data/visualization_data/" * save_to_string * ".csv"
      CSV.write(loc, df)
      df = DataFrame(["x" => x, "y" => y, "z" => z, "partition" => partitions])
      loc = "./data/visualization_data/PART" * save_to_string * ".csv"
      CSV.write(loc, df)
    end
end
