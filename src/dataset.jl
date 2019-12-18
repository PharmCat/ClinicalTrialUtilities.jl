struct DataSet{T <: AbstractData}
    data::Vector{T}
end
function Base.iterate(iter::DataSet)
    return Base.iterate(iter.data)
end
function Base.iterate(iter::DataSet, i::Int)
    return Base.iterate(iter.data, i)
end
