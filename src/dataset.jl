struct DataSet{T <: AbstractData}
    data::Vector{T}
end

#!!!
function in(a::Dict, b::Dict)
    k = collect(keys(a))
    if any(x -> x  ∉  collect(keys(b)), k) return false end
    for i = 1:length(a)
        if a[k[i]] != b[k[i]] return false end
    end
    return true
end
function in(a::Pair, b::Dict)
    if a[1]  ∉  collect(keys(b)) return false end
    if a[2] != b[a[1]] return false end
    return true
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
function Base.findall(a::DataSet{T}, sort::Dict) where T
    inds = Vector{Int}(undef, 0)
    for i = 1:length(a)
        if sort ∈ a.data[i].sort push!(inds, i) end
    end
    inds
end
#
function Base.deleteat!(a::DataSet{T}, i::Int) where T
    deleteat!(a.data, i)
    return a
end

function Base.deleteat!(a::DataSet{T}, inds::AbstractArray) where T
    deleteat!(a.data, inds)
    return a
end

function Base.deleteat!(a::DataSet{T}, inds::Dict) where T
    deleteat!(a.data, findall(a, inds))
    return a
end
