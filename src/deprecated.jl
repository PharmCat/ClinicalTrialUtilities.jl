function designProp(type::Symbol)
    if type == :parallel
        return Parallel()
    elseif type == :d2x2
        return Crossover(x -> x - 2, 0.5, 2)
    elseif type == :d2x2x3
        return Crossover(x -> 2 * x - 3, 0.375, 2)
    elseif type == :d2x2x4
        return Crossover(x -> 3 * x - 4, 0.25, 2)
    elseif type == :d2x4x4
        return  Crossover(x -> 3 * x - 4, 0.0625, 4)
    elseif type == :d2x3x3
        return Crossover(x -> 2 * x - 3, 1/6, 3)
    elseif type == :d2x4x2
        return Crossover(x -> x - 2, 0.5, 4)
    elseif type == :d3x3
        return Crossover(x -> 2 * x - 4, 2/9, 3)
    elseif type == :d3x6x3
        return Crossover(x -> 2 * x - 4, 1/18, 6)
    else throw(ArgumentError("Design type not known!")) end
end

#Exceptions
struct CTUException <: Exception
    n::Int
    var::String
end

Base.showerror(io::IO, e::CTUException) = print("CTU Exception code: ", e.n, " Message: ", e.var);

#Descriptives
function descriptive_deprecated(data::DataFrame;
    sort::Union{Symbol, Array{T, 1}, Nothing} = nothing,
    vars::Union{Symbol, Array{T, 1},  Nothing} = nothing,
    stats::Union{Symbol, Array{T, 1}} = [:n, :mean, :sd, :sem, :uq, :median, :lq])::DataFrame where T <: Union{Symbol, String}
    #Filtering
    dfnames = names(data) # Col names of dataframe
    #Filter sort
    if isa(sort, Array)
        sort = Symbol.(sort)
        filter!(x->x in dfnames, sort)
        if length(sort) == 0 sort = nothing end
    elseif isa(sort, Symbol)
        if any(x -> x == sort, dfnames)
            sort = [sort]
        else
            sort =  nothing
        end
    else
        sort =  nothing
    end
    #Filters vars
    if isa(vars, Symbol) vars = [vars] end
    if isa(vars, Array)
        vars = Symbol.(vars)
        filter!(x->x in dfnames, vars)
        filter!(x-> !(x in sort), vars)
        if length(vars) > 0
            for i = 1:length(vars)
                if !(eltype(data[:, vars[1]]) <: Real) deleteat!(vars, i) end
            end
        end
        if length(vars) == 0 error("No variables of type Real[] in dataset found! Check vars array!") end
    else
        vars = Array{Symbol,1}(undef, 0)
        for i in dfnames
            if eltype(data[:, i]) <: Real push!(vars, i) end
        end
        if sort !== nothing filter!(x-> !(x in sort), vars) end
        if length(vars) == 0 error("Real[] columns not found!") end
    end
    #Filter statistics array
    if stats == :all
        stats = [:n, :min, :max, :range, :mean, :var, :sd, :sem, :cv, :harmmean, :geomean, :geovar, :geosd, :geocv, :skew, :kurt, :uq, :median, :lq, :iqr, :mode]
    elseif isa(stats, Array)
        stats = Symbol.(stats)
        filter!(x->x in [:n, :min, :max, :range, :mean, :var, :sd, :sem, :cv, :harmmean, :geomean, :geovar, :geosd, :geocv, :skew, :ses, :kurt, :sek, :uq, :median, :lq, :iqr, :mode], stats)
        if length(stats) == 0
            stats = [:n, :mean, :sd, :sem, :uq, :median, :lq]
            @warn "Error in stats, default used."
        end
    else
        stats = [:n, :mean, :sd, :sem, :uq, :median, :lq]
        @warn "Unknown stats, default used."
    end
    #End filtering

    #construct datasets
    dfs   = DataFrame(vars = Symbol[]) #Init DataFrames
    sdata = DataFrame()                #Temp dataset for sorting

    if isa(sort, Array)
        for i in sort                  #if sort - make sort columns in dataset
            dfs[:, i] = eltype(data[:, i])[]
        end
    end
    for i in stats                     #make columns for statistics
        dfs[:, i] = Real[]
    end
    for i in vars                      #var columns for temp dataset
        sdata[:, i] = Real[]
    end
    if !isa(sort, Array)               #if no sort - for vars without sorting
        for v in vars
            sarray = Array{Any,1}(undef, 0)
            deleteat!(data[:, v], findall(x->x === NaN || x === nothing, data[:, v]))
            if length(data[:, v]) > 1
                push!(sarray, v)
                append!(sarray, descriptive_(data[:, v], stats))
                push!(dfs, sarray)
            end
        end
        return dfs
    end

    #If sort exist
    sortlist = unique(data[:, sort])
    sort!(sortlist, tuple(sort))
    for i = 1:size(sortlist, 1) #For each line in sortlist
            if size(sdata, 1) > 0 deleterows!(sdata, 1:size(sdata, 1)) end
            for c = 1:size(data, 1) #For each line in data
                if data[c, sort] == sortlist[i,:]
                    push!(sdata, data[c, vars])  #tmp ds constructing
                end
            end
        for v in vars #For each variable in list
            sarray = Array{Any,1}(undef, 0)
            push!(sarray, v)
            for c in sort
                push!(sarray, sortlist[i, c])
            end
            append!(sarray, descriptive_(sdata[:, v], stats))
            push!(dfs, sarray)
        end
    end
    return dfs
end
