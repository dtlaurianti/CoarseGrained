function getDictSet(partitionData::String, partitionNumber::Integer)
    csv_reader = CSV.File(partitionData)
    i = 1
    str = ""
    for row in csv_reader
        if i == partitionNumber
            str = row.partition
        end
        i+= 1
    end
    cleanStr = replace(str, "Dict{Integer, Integer}("=>"")
    arr = split(cleanStr, ", ")
    result = Set(arr)
    return result
end

function getPartArraySet(partitionData::String, partitionNumber::Integer)
    csv_reader = CSV.File(partitionData)
    i = 1
    str = ""
    for row in csv_reader
        if i == partitionNumber
            str = row.partArray
        end
        i+= 1
    end
    cleanStr = replace(str, "Vector{Any}[["=>"")
    arr = split(cleanStr, "], [")
    result = Set(arr)
    return result
end

function getMaxDictPartition(partitionData::String, PartitionNumber::Integer)
    
end

function batchIntersect(partitionData::String; startingPartition=1, endingPartition=0)
    if endingPartition == 0
        max = countcsvlines(partitionData) - 1
    else
        max = endingPartition
    end
    DictSets = Vector()
    PartArraySets = Vector()
    for i in startingPartition:max
        DictSet = replace(getDictSet(partitionData, i), "Set(SubString{String}["=>"")
        PartArraySet = replace(getPartArraySet(partitionData, i), "Set(SubString{String}["=>"")
        push!(DictSets, DictSet)
        push!(PartArraySets, PartArraySet)
    end

    result1 = intersect(DictSets)
    result2 = intersect(PartArraySets)
    return result1, result2
end