abstract type AbstractData end
abstract type AbstractDataSet end

"""
    Descriptive statistics type
"""
struct Descriptive <: AbstractData
    var::Union{Symbol, Nothing}
    varname::Union{String, Nothing}
    sort::Dict
    #data::NamedTuple{}
    data::NamedTuple{S,T} where S where T <: Tuple{Vararg{Number}}
end
function Base.show(io::IO, obj::Descriptive)
    println(io, "Descriptive statistics")
    if obj.var !== nothing print(io, "Variable: ", obj.var) end
    if obj.varname !== nothing println(io, "(", obj.varname,")") else println(io, "") end
    if obj.sortval !== nothing println(io, "Group: ", obj.sortval) end
    maxcol = 0
    for k in keys(obj.data)
        if length(string(k)) > maxcol maxcol = length(string(k)) end
    end
    for k in keys(obj.data)
        spv = Vector{Char}(undef, maxcol-length(string(k)))
        spv .= ' '
        val = obj.data[k]
        if typeof(val) <: AbstractFloat val = round(val; digits = 8) end
        val =  string(val)
        println(io, string(k), String(spv), " | ", val)
        spv = Vector{Char}(undef, maxcol)
        spv .= '_'
        println(io, String(spv), "___", "____________")
    end
end
function Base.getindex(a::Descriptive, s::Symbol)::Real
    return a.data[s]
end
struct DataSet{T <: AbstractData}
    data::Vector{T}
end
function Base.show(io::IO, obj::DataSet{Descriptive})
    println(io, "Descriptive statistics set")
    for i = 1:length(obj.data)
        print(io, i, " | ", string(obj.data[i].var))
        if obj.data[i].sortval !== nothing println(io, " | ", string(obj.data[i].sortval)) end
        println(io, "____________________")
    end
end
function Base.getindex(a::DataSet{Descriptive}, i::Int64)::Descriptive
    return a.data[i]
end
function Base.getindex(a::DataSet{Descriptive}, i::Int64, s::Symbol)::Real
    return a.data[i].data[s]
end
function Base.length(data::DataSet{T}) where T <: AbstractData
    return length(data.data)
end
"""
    Descriptive statistics
"""
function descriptive(data::DataFrame;
    sort::Union{Symbol, Array{T,1}} = Array{Symbol,1}(undef,0),
    vars::Union{Symbol, Array{T,1}},
    stats::Union{Symbol, Array{T,1}, Tuple{Vararg{Symbol}}} = :default)::DataSet{Descriptive} where T <: Union{Symbol, String}

    stats = checkstats(stats)

    if isa(vars, Symbol) vars = [vars] end
    if eltype(vars) <: String vars = Symbol.(vars) end

    if isa(sort, Symbol) sort = [sort] end
    if eltype(sort) <: String sort = Symbol.(sort) end

    d = Array{Descriptive, 1}(undef, 0) # Temp Descriptive array

    if isa(sort, Array) && length(sort) == 0
        #pushvardescriptive!(d, vars, Matrix(data[:,vars]), Dict(), stats)
        mx = Matrix(data[:,vars])
        for v  = 1:length(vars)  #For each variable in list
            push!(d, Descriptive(vars[v], nothing, Dict(), NamedTuple{stats}(Tuple(descriptive_(mx[:, v], stats)))))
        end
        return DataSet(d)
    end

    sortlist = unique(data[:, sort])
    #sort!(sortlist, sort)
    sortlist = Matrix(sort!(sortlist, sort)) #more allocs less time (or not :( )

    for i = 1:size(sortlist, 1)
        sortval = Tuple(sortlist[i,:])
        mx = getsortedmatrix(data; datacol=vars, sortcol=sort, sortval=sortval)
        for v  = 1:length(vars)  #For each variable in list
            push!(d, Descriptive(vars[v], nothing, Dict(sort .=> sortval), NamedTuple{stats}(Tuple(descriptive_(mx[:, v], stats)))))
        end
        #pushvardescriptive!(d, vars, mx, sortval, stats)  #push variable descriprives for mx
    end
    return DataSet(d)
end
function descriptive(data::Array{T, 1}; stats::Union{Symbol, Vector, Tuple} = :default, var = nothing, varname = nothing, sort = Dict())::Descriptive where T <: Real
    stats = checkstats(stats)
    return Descriptive(var, varname, sort, NamedTuple{stats}(Tuple(descriptive_(data, stats))))
end
#=
"""
"""
@inline function usort(data::DataFrame, sort)
    if size(data, 1) == 0 throw(ArgumentError("DataFrame size equal 0")) end
    sortlist = Array{Array{Any,1},1}(undef, 0)
    push!(sortlist, getrow(data, sort,1))
    if size(data, 1) == 1 return sortlist end
    for i = 2:size(data, 1)
        v = getrow(data, sort, i)
        if v ∉ sortlist push!(sortlist, v) end
    end
    return sortlist
end
"""
"""
@inline function getrow(df::DataFrame, sort::Array{T, 1}, i::Int)::Array{Any, 1} where T
    v = Array{Any, 1}(undef,0)
    for c = 1:length(sort)
        push!(v, df[i, sort[c]])
    end
    return v
end
=#
"""
    Check if all statistics in allstat list. return stats tuple
"""
@inline function checkstats(stats::Union{Symbol, Array{T,1}, Tuple{Vararg{Symbol}}})::Tuple{Vararg{Symbol}} where T <: Union{Symbol, String}
    allstat = (:n, :min, :max, :range, :mean, :var, :sd, :sem, :cv, :harmmean, :geomean, :geovar, :geosd, :geocv, :skew, :ses, :kurt, :sek, :uq, :median, :lq, :iqr, :mode)
    if isa(stats, Symbol)
        if stats == :default stats = (:n, :mean, :sd, :sem, :uq, :median, :lq)
        elseif stats == :all stats = allstat
        else stats = Tuple(stats) end
    end
    stats = Tuple(Symbol.(stats))
    if any(x -> x  ∉  allstat, stats) throw(ArgumentError("stats element not in allstat list")) end
    return stats
end
"""
    Push in d Descriptive obj in mx vardata
"""
@inline function pushvardescriptive!(d::Array{Descriptive, 1}, vars::Array{Symbol, 1}, mx::Union{DataFrame, Matrix{T}}, sortval::Union{Tuple{Vararg{Any}}, Nothing}, stats::Tuple{Vararg{Symbol}}) where T<: Real
    for v  = 1:length(vars)  #For each variable in list
        push!(d, Descriptive(vars[v], nothing, sortval, NamedTuple{stats}(Tuple(descriptive_(mx[:, v], stats)))))
    end
end
"""
    Check if data row sortcol equal sortval
"""
@inline function checksort(data::DataFrame, row::Int, sortcol::Array{Symbol, 1}, sortval::Tuple{Vararg{Any}})::Bool
    for i = 1:length(sortcol)
        if data[row, sortcol[i]] != sortval[i] return false end
    end
    return true
end
"""
    Return matrix of filtered data (datacol) by sortcol with sortval
"""
@inline function getsortedmatrix(data::DataFrame; datacol::Array{Symbol,1}, sortcol::Array{Symbol,1}, sortval::Tuple{Vararg{Any}})::Matrix{Real}
    result  = Array{Real, 1}(undef, 0)
    for c = 1:size(data, 1) #For each line in data
        if checksort(data, c, sortcol, sortval)
            for i = 1:length(datacol)
                @inbounds push!(result, data[c, datacol[i]])
            end
        end
    end
    return Matrix(reshape(result, length(datacol), :)')
end

#Statistics calculation
@inline function descriptive_(data::Array{T,1}, stats::Union{Tuple{Vararg{Symbol}}, Array{Symbol,1}})::Array{Real,1} where T <: Real
    deleteat!(data, findall(x->x === NaN || x === nothing || x === missing, data))

    dn         = nothing
    #if (@isdefined dn) end
    dmin       = nothing
    dmax       = nothing
    drange     = nothing
    dmean      = nothing
    dvar       = nothing
    dsd        = nothing
    dsem       = nothing
    dcv        = nothing
    dgeomean   = nothing
    dgeovar    = nothing
    dgeosd     = nothing
    dgeocv     = nothing
    dharmmean  = nothing
    duq        = nothing
    dlq        = nothing
    sesv       = nothing
    sekv       = nothing
    #dirq       = nothing
    #dmode      = nothing
    if length(data) > 0
        sarray = Array{Real,1}(undef, 0)
    elseif length(data) <= 2
        sesv = NaN
        sekv = NaN
    elseif length(data) <= 3
        sekv = NaN
    elseif length(data) == 0
        sarray = Array{Real,1}(undef, length(stats))
        return sarray .= NaN
    end

    if any(x -> (x == :geomean || x == :geocv || x == :geosd || x == :geovar), stats)
        if any(x -> (x <= 0), data)
            dgeomean   = NaN
            dgeovar    = NaN
            dgeosd     = NaN
            dgeocv     = NaN
        else
            logsdata = log.(data)
        end
    end
    for s in stats
        if s == :n
            if dn === nothing dn = length(data) end
            push!(sarray, dn)
        elseif s == :sum
            push!(sarray, sum(data))
        elseif s == :min
            if dmin === nothing dmin = minimum(data) end
            push!(sarray, dmin)
        elseif s == :max
            if dmax === nothing dmax = maximum(data) end
            push!(sarray, dmax)
        elseif s == :range
            if dmax === nothing dmax = max(data) end
            if dmin === nothing dmin = max(data) end
            push!(sarray, abs(dmax-dmin))
        elseif s == :mean
            if dmean === nothing dmean = mean(data) end
            push!(sarray, dmean)
        elseif s == :var
            if dvar === nothing dvar = var(data) end
            push!(sarray, dvar)
        elseif s == :sd
            if dvar === nothing dvar = var(data) end
            if dsd  === nothing dsd  = sqrt(dvar) end
            push!(sarray, dsd)
        elseif s == :sem
            if dvar === nothing dvar = var(data) end
            if dn === nothing dn = length(data) end
            push!(sarray, sqrt(dvar/dn))
        elseif s == :cv
            if dmean === nothing dmean = mean(data) end
            if dvar === nothing dvar = var(data) end
            if dsd  === nothing dsd  = sqrt(dvar) end
            push!(sarray, dsd/dmean*100)
        elseif s == :harmmean
            if dharmmean === nothing
                if any(x -> (x == 0), data)
                    dharmmean = NaN
                else
                    dharmmean = harmmean(data)
                end
            end
            push!(sarray, dharmmean)
        elseif s == :geomean
            if dgeomean === nothing dgeomean = exp(mean(logsdata)) end
            push!(sarray, dgeomean)
        elseif s == :geovar
            if dgeovar === nothing dgeovar = var(logsdata) end
            push!(sarray, dgeovar)
        elseif s == :geosd
            if dgeovar === nothing dgeovar = var(logsdata) end
            if dgeosd === nothing dgeosd = sqrt(dgeovar) end
            push!(sarray, dgeosd)
        elseif s == :geocv
            if dgeovar === nothing dgeovar = var(logsdata) end
            if dgeocv === nothing dgeocv = sqrt(exp(dgeovar)-1)*100 end
            push!(sarray, dgeocv)
        elseif s == :skew
            push!(sarray, skewness(data))
        elseif s == :ses
            if dn === nothing dn = length(data) end
            sesv = ses(dn)
            push!(sarray, sesv)
        elseif s == :kurt
            push!(sarray, kurtosis(data))
        elseif s == :sek
            if dn === nothing dn = length(data) end
            if sesv === nothing sesv = ses(dn) end
            sekv = sek(dn; ses = sesv)
            push!(sarray, sekv)
        elseif s == :uq
            if duq  === nothing duq  = percentile(data, 75) end
            push!(sarray, duq)
        elseif s == :median
            push!(sarray,  median(data))
        elseif s == :lq
            if dlq  === nothing dlq  = percentile(data, 25) end
            push!(sarray, dlq)
        elseif s == :iqr
            if duq  === nothing duq  = percentile(data, 75) end
            if dlq  === nothing dlq  = percentile(data, 25) end
            push!(sarray, abs(duq-dlq))
        elseif s == :mode
            push!(sarray, mode(data))
        end
    end
    return sarray
end

@inline function ses(data::AbstractVector)::Float64
    n = length(data)
    ses(n)
end
@inline function ses(n::Int)::Float64
    return sqrt(6 * n *(n - 1) / ((n - 2) * (n + 1) * (n + 3)))
end

function sek(data::AbstractVector; ses::T = ses(data))::Float64 where T <: Real
    n = length(data)
    sek(n; ses = ses)
end
function sek(n::Int; ses::T = ses(n))::Float64 where T <: Real
    return 2 * ses * sqrt((n * n - 1)/((n - 3) * (n + 5)))
end
