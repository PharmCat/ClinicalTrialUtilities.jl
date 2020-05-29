"""
    Descriptive statistics type
"""
struct Descriptive <: AbstractData
    var::Union{Symbol, Nothing}
    varname::Union{String, Nothing}
    sort::Dict
    result::Dict
end
function Base.show(io::IO, obj::Descriptive)
    println(io, "Descriptive statistics")
    if obj.var     !== nothing print(io, "Variable: ", obj.var) end
    if obj.varname !== nothing println(io, "(", obj.varname,")") else println(io, "") end
    if obj.sort    !== nothing println(io, "Group: ", obj.sort) end
    maxcol = 0
    for k in keys(obj.result)
        if length(string(k)) > maxcol maxcol = length(string(k)) end
    end
    for k in keys(obj.result)
        spv = Vector{Char}(undef, maxcol-length(string(k)))
        spv .= ' '
        val = obj.result[k]
        if typeof(val) <: AbstractFloat val = round(val; digits = 8) end
        val =  string(val)
        println(io, string(k), String(spv), " | ", val)
        spv = Vector{Char}(undef, maxcol)
        spv .= '_'
        println(io, String(spv), "___", "____________")
    end
end
function Base.getindex(a::Descriptive, s::Symbol)::Real
    return a.result[s]
end
function Base.show(io::IO, obj::DataSet{Descriptive})
    println(io, "Descriptive statistics set")
    for i = 1:length(obj.data)
        print(io, i, " | ", string(obj.data[i].var))
        if obj.data[i].sort !== nothing println(io, " | ", string(obj.data[i].sort)) end
        println(io, "____________________")
    end
end
"""
    descriptive(data::DataFrame;
        sort::Union{Symbol, Array{T,1}} = Array{Symbol,1}(undef,0),
        vars = [],
        stats::Union{Symbol, Array{T,1}, Tuple{Vararg{Symbol}}} = :default)::DataSet{Descriptive} where T <: Union{Symbol, String}

Descriptive statistics.

- ``sort`` sorting columns
- ``vars`` variabels
- ``stats`` statistics
"""
function descriptive(data::DataFrame;
    sort::Union{Symbol, Array{T,1}} = Array{Symbol,1}(undef,0),
    vars = [],
    stats::Union{Symbol, Array{T,1}, Tuple{Vararg{Symbol}}} = :default)::DataSet{Descriptive} where T <: Union{Symbol, String}

    stats = checkstats(stats)

    if isa(vars, Array) && length(vars) == 0
        vars = filter(x -> x ∉ sort, names(data))
        del  = []
        for i = 1:length(vars)
            if !(eltype(data[!, vars[i]]) <: Union{Missing, Real}) push!(del, i) end
        end
        if length(del) > 0
            deleteat!(vars, del)
        end
    end
    if isa(vars, Symbol)
        vars = [vars]
    elseif !(isa(vars, Array))
        vars = [Symbol(vars)]
    end
    if isa(vars, Array) && !(eltype(vars) <: Symbol) vars = Symbol.(vars) end

    if isa(sort, Symbol) sort = [sort] end
    if eltype(sort) <: String sort = Symbol.(sort) end

    d = Array{Descriptive, 1}(undef, 0) # Temp Descriptive array

    if isa(sort, Array) && length(sort) == 0
        #pushvardescriptive!(d, vars, Matrix(data[:,vars]), Dict(), stats)
        mx = Matrix(data[:,vars])
        for v  = 1:length(vars)  #For each variable in list
            push!(d, Descriptive(vars[v], nothing, Dict(:Variable => v), descriptive_(mx[:, v], stats)))
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
            dict = Dict{Symbol, Any}(sort .=> sortval)
            dict[:Variable]   =  vars[v]
            push!(d, Descriptive(vars[v], nothing, dict, descriptive_(mx[:, v], stats)))
        end
        #pushvardescriptive!(d, vars, mx, sortval, stats)  #push variable descriprives for mx
    end
    return DataSet(d)
end
function descriptive(data::Array{T, 1}; stats::Union{Symbol, Vector, Tuple} = :default, var = nothing, varname = nothing, sort = Dict())::Descriptive where T <: Real
    stats = checkstats(stats)
    return Descriptive(var, varname, sort, descriptive_(data, stats))

end

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
        push!(d, Descriptive(vars[v], nothing, sortval, descriptive_(mx[:, v], stats)))
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

function notnan(x)
    return !(x === NaN || x === nothing || x === missing)
end

@inline function descriptive_(data::Vector{T}, stats::Union{Tuple{Vararg{Symbol}}, Array{Symbol,1}}) where T <: Real

    #=
    dlist = findall(x -> x === NaN || x === nothing || x === missing, data)
    if length(dlist) > 0
        data = copy(data)
        deleteat!(data, dlist)
    end
    =#

    data       = data[notnan.(data)]

    dn         = nothing
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
    #dgeocv     = nothing
    #dharmmean  = nothing
    duq        = nothing
    dlq        = nothing
    sesv       = nothing





    dict = Dict{Symbol, Real}()


    if length(data) == 0
        for i in stats
            dict[i] = NaN
        end
        return dict
    end

    if :sesv in stats
        if length(data) <= 2
            sesv = NaN
            dict[:sesv] = sesv
        end
    end
    if :sekv in stats
        if length(data) <= 3
            sekv = NaN
            dict[:sekv] = sekv
        end

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
            dict[s] = dn
        elseif s == :sum
            dict[s] = sum(data)
        elseif s == :min
            if dmin === nothing dmin = minimum(data) end
            dict[s] = dmin
        elseif s == :max
            if dmax === nothing dmax = maximum(data) end
            dict[s] = dmax
        elseif s == :range
            if dmax === nothing dmax = minimum(data) end
            if dmin === nothing dmin = maximum(data) end
            dict[s] =  abs(dmax-dmin)
        elseif s == :mean
            if dmean === nothing dmean = mean(data) end
            dict[s] = dmean
        elseif s == :var
            if dvar === nothing dvar = var(data) end
            dict[s] = dvar
        elseif s == :sd
            if dvar === nothing dvar = var(data) end
            if dsd  === nothing dsd  = sqrt(dvar) end
            dict[s] = dsd
        elseif s == :sem
            if dvar === nothing dvar = var(data) end
            if dn === nothing dn = length(data) end
            dict[s] = sqrt(dvar/dn)
        elseif s == :cv
            if dmean === nothing dmean = mean(data) end
            if dvar === nothing dvar = var(data) end
            if dsd  === nothing dsd  = sqrt(dvar) end
            dict[s] = dsd/dmean*100
        elseif s == :harmmean
            if any(x -> (x == 0), data)
                dict[s] = NaN
            else
                dict[s] = harmmean(data)
            end
        elseif s == :geomean
            if dgeomean === nothing dgeomean = exp(mean(logsdata)) end
            dict[s] =  dgeomean
        elseif s == :geovar
            if dgeovar === nothing dgeovar = var(logsdata) end
            dict[s] = dgeovar
        elseif s == :geosd
            if dgeovar === nothing dgeovar = var(logsdata) end
            if dgeosd  === nothing dgeosd = sqrt(dgeovar) end
            dict[s] = dgeosd
        elseif s == :geocv
            if dgeovar === nothing dgeovar = var(logsdata) end
            dict[s] = sqrt(exp(dgeovar)-1)*100
        elseif s == :skew
            dict[s] = skewness(data)
        elseif s == :ses
            if length(data) <= 2
                dict[s] = NaN
            else
                if dn   === nothing dn   = length(data) end
                if sesv === nothing sesv = ses(dn) end
                dict[s] = sesv
            end
        elseif s == :kurt
            dict[s] = kurtosis(data)
        elseif s == :sek
            if length(data) <= 3
                dict[s] = NaN
            else
                if dn   === nothing dn = length(data) end
                if sesv === nothing sesv = ses(dn) end
                dict[s] = sek(dn; ses = sesv)
            end
        elseif s == :uq
            if duq  === nothing duq  = percentile(data, 75) end
            dict[s] = duq
        elseif s == :median
            dict[s] = median(data)
        elseif s == :lq
            if dlq  === nothing dlq  = percentile(data, 25) end
            dict[s] = dlq
        elseif s == :iqr
            if duq  === nothing duq  = percentile(data, 75) end
            if dlq  === nothing dlq  = percentile(data, 25) end
            dict[s] = abs(duq-dlq)
        elseif s == :mode
            dict[s] = mode(data)
        end
    end
    return dict
end

@inline function ses(data::AbstractVector)
    n = length(data)
    ses(n)
end
@inline function ses(n::Int)
    return sqrt(6 * n *(n - 1) / ((n - 2) * (n + 1) * (n + 3)))
end

@inline function sek(data::AbstractVector; ses::T = ses(data)) where T <: Real
    n = length(data)
    sek(n; ses = ses)
end
@inline  function sek(n::Int; ses::T = ses(n)) where T <: Real
    return 2 * ses * sqrt((n * n - 1)/((n - 3) * (n + 5)))
end
