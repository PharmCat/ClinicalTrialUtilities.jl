struct DataSet{T <: AbstractData}
    data::Vector{T}
end
#length
function Base.length(a::DataSet)::Int
    return length(a.data)
end
#Iterate
function Base.iterate(iter::DataSet)
    return Base.iterate(iter.data)
end
function Base.iterate(iter::DataSet, i::Int)
    return Base.iterate(iter.data, i)
end
#Eltype
function Base.eltype(obj::DataSet)
    return eltype(obj.data)
end
#Getindex
function Base.getindex(a::DataSet{T}, i::Int)::T where T
    return a.data[i]
end
function Base.getindex(a::DataSet{T}, d::Pair)::T where T
    for i = 1:length(a)
        if d ∈ a[i].sort return a[i] end
    end
    throw(ArgumentError("Element not found!"))
end
function Base.getindex(a::DataSet{T}, d::Dict)::T where T
    for i = 1:length(a)
        if d ∈ a[i].sort return a[i] end
    end
    throw(ArgumentError("Element not found!"))
end
function Base.getindex(a::DataSet{T}, d::Tuple{Vararg{Pair}})::T where T
    return a[Dict(d)]
end


function Base.getindex(a::DataSet{T}, i::Int, s::Symbol)::Real where T
    return a.data[i].result[s]
end
#
