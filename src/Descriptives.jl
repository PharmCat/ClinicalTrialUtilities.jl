
function descriptives(data::DataFrame; sort = NaN, vars = NaN, stats = [:n, :mean, :sd, :sem, :uq, :median, :lq])::DataFrame
    #Filtering
    if isa(sort, Array)
        sort = Symbol.(sort)
        filter!(x->x in names(data), sort)
        if length(sort) == 0 sort = NaN end
    end
    if isa(vars, Array)
        vars = Symbol.(vars)
        filter!(x->x in names(data), vars)
        filter!(x-> !(x in sort), vars)
        if length(vars) > 0
            for i = 1:length(vars)
                if !(eltype(data[vars[1]]) <: Real) deleteat!(vars, i) end
            end
        end
        if length(vars) == 0 error("No variables of type Real[] in dataset found! Check vars array!") end
    else
        vars = Array{Symbol,1}(undef, 0)
        for i in names(data)
            if eltype(data[i]) <: Real push!(vars, i) end
        end
        filter!(x-> !(x in sort), vars)
        if length(vars) == 0 error("Real[] columns not found!") end
    end
    #End filtering
    dfs   = DataFrame(vars = Symbol[]) #Init DataFrames
    sdata = DataFrame()                #Temp dataset for sorting
    if isa(sort, Array)
        for i in sort                  #if sort - make sort columns in dataset
            dfs[i] = eltype(data[i])[]
        end
    end
    for i in stats                     #make columns for statistics
        dfs[i] = Real[]
    end
    for i in vars                      #var columns for temp dataset
        sdata[i] = Real[]
    end
    if !isa(sort, Array)               #if no sort
        for v in vars
            sarray = Array{Any,1}(undef, 0)
            deleteat!(data[v], findall(x->x === NaN || x === nothing, data[v]))
            if length(data[v]) > 1
                push!(sarray, v)
                descriptives_!(sarray, data[v], stats)
                push!(dfs, sarray)
            end
        end
        return dfs
    end
    #If sort...
    sortlist = unique(data[sort])
    for v in vars #For each variable in list
        for i = 1:size(sortlist, 1) #For each line in sortlist
            if size(sdata, 1) > 0 deleterows!(sdata, 1:size(sdata, 1)) end
            for c = 1:size(data, 1) #For each line in data
                if data[c, sort] == sortlist[i,:]
                    push!(sdata, data[c, vars])
                end
            end
        sarray = Array{Any,1}(undef, 0)
        push!(sarray, v)
            for c in sort
                push!(sarray, sortlist[i, c])
            end
            deleteat!(sdata[v], findall(x->x === NaN || x === nothing, sdata[v]))
            if length(sdata[v]) > 1
                descriptives_!(sarray, sdata[v], stats)
                push!(dfs, sarray)
            end
        end
    end
    return dfs
end

#Statistics calculation
function descriptives_!(sarray, data, stats)

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
    dgeocv     = nothing
    dharmmean  = nothing
    duq        = nothing
    dlq        = nothing
    #dirq       = nothing
    #dmode      = nothing

    if any(x -> (x == :geomean || x == :geomeancv || x == :geomeansd), stats)
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
        elseif s == :min
            if dmin === nothing dmin = min(data) end
            push!(sarray, dmin)
        elseif s == :max
            if dmax === nothing dmax = max(data) end
            push!(sarray, dmax)
        elseif s == :range
            if dmax === nothing dmax = max(data) end
            if dmim === nothing dmin = max(data) end
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
        elseif s == :kurt
            push!(sarray, kurtosis(data))
        elseif s == :uq
            if duq  === nothing duq  = percentile(data, 75) end
            push!(sarray, duq)
        elseif s == :median
            push!(sarray,  median(data))
        elseif s == :lq
            if dlq  === nothing dlq  = percentile(data, 25) end
            push!(sarray, dlq)
        elseif s == :irq
            if duq  === nothing duq  = percentile(data, 75) end
            if dlq  === nothing dlq  = percentile(data, 25) end
            push!(sarray, abs(duq-dlq))
        elseif s == :mode
            push!(sarray, mode(data))
        end
    end
end
