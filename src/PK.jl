module PK

import ..NCA

export nca

using DataFrames

function nca(data; conc=:Concentration, time=:Time, sort = NaN, stack = false)
    if isa(conc, String)  conc = Symbol(conc) end
    if isa(time, String)  time = Symbol(time) end

    if sort === NaN
        nca, rsqn, elim = nca_(data, conc, time)
    else
        sortlist = unique(data[sort])
        nca  = DataFrame(AUClast = Float64[], Cmax = Float64[], Tmax = Float64[], AUMClast = Float64[], MRTlast = Float64[], Kel = Float64[], HL = Float64[], Rsq = Float64[], AUCinf = Float64[], AUCpct = Float64[])
        elim = DataFrame(Start = Int[], End = Int[], b = Float64[], a = Float64[], Rsq = Float64[])
        for i = 1:size(sortlist, 1) #For each line in sortlist
            datai = DataFrame(Concentration = Float64[], Time = Float64[])
            for c = 1:size(data, 1) #For each line in data
                if data[c, sort] == sortlist[i,:]
                    push!(datai, [data[c, conc], data[c, time]])
                end
            end
            ncares = nca_(datai, :Concentration, :Time)
            append!(nca, ncares[1])
            append!(elim, DataFrame(ncares[3][ncares[2],:]))
        end
        nca  = hcat(sortlist, nca)
        elim = hcat(sortlist, elim)
    end
    return NCA(nca, elim, DataFrame(), "", "")
end

function nca_(data, conc, time)
    if length(unique(data[:,time])) != length(data[:,time])
        return  DataFrame(AUClast = [NaN], Cmax = [NaN], Tmax = [NaN], AUMClast = [NaN], MRTlast = [NaN], Kel = [NaN], HL = [NaN], Rsq = [NaN], AUCinf = [NaN], AUCpct = [NaN])
    end
    pklog = "NCA analysis: \n"
    pklog *= "Concentration column: "*string(conc)*"\n"
    pklog *= "Time column: "*string(time)*"\n"
    sort!(data, [time])
    pklog *= "Sorting by Time... \n"
    n::Int     = nrow(data)
    cmax = maximum(data[conc])               #Cmax
    clast =  data[n, conc]
    pklog *= "Cmax = "*string(cmax)*"\n"
    auc   = 0                                  #AUClast
    aumc  = 0                                 #AUMClast
    mrt   = 0
    tmax  = NaN                               #Tmax
    tmaxn::Int = 0                             #Time point Tmax
    kel   = NaN
    hl    = NaN
    rsq   = NaN
    aucinf = NaN
    aucinfpct = NaN
    mrtlast = NaN

    pklog *= "AUC calculation:\n"

    for i = 2:nrow(data)
        auc_ = (data[i - 1, conc] + data[i, conc]) * (data[i, time] - data[i - 1, time])/2
        aumc_ = (data[i - 1, conc]*data[i - 1, time] + data[i, conc]*data[i, time]) * (data[i, time] - data[i - 1, time])/2
        auc += auc_
        aumc += aumc_
        pklog *= "AUC("*string(i-1)*") = "*string(auc_)*"; Sum = "*string(auc)*"\n"
    end
    mrt = aumc / auc
    pklog *= "AUClast = "*string(auc)*"\n"
    for i = 1:nrow(data)
        if data[i, conc] == cmax
            tmax = data[i, time]
            tmaxn = i
            pklog *= "Tmax = "*string(tmax)*"\n"
            break
        end
    end
    if n - tmaxn >= 3
        keldata = DataFrame(Start = Int[], End = Int[], b = Float64[], a = Float64[], Rsq = Float64[])
        logconc = log.(data[conc])
        for i::Int = tmax+1:n-2
            sl = slope(data[i:n,time], logconc[i:n])
            push!(keldata, [i, n, sl[1], sl[2], sl[3]])
        end

        rsq, rsqn = findmax(keldata[:Rsq])
        kel = abs(keldata[rsqn,:a])
        hl  = 0.6931471805599453/kel
        aucinf = auc + clast/kel
        aucinfpct = (aucinf - auc) / aucinf * 100.0
    end

    return DataFrame(AUClast = [auc], Cmax = [cmax], Tmax = [tmax], AUMClast = [aumc], MRTlast = [mrt], Kel = [kel], HL = [hl], Rsq = [rsq], AUCinf = [aucinf], AUCpct = [aucinfpct]), rsqn, keldata
end

function slope(x,y)
    n   = length(x)
    Σxy = sum(x .* y)
    Σx  = sum(x)
    Σy  = sum(y)
    Σx2 = sum( x .* x)
    Σy2 = sum( y .* y)

    a   = (Σy * Σx2 - Σx * Σxy)/(n * Σx2 - Σx^2)
    b   = (n * Σxy - Σx * Σy)/(n * Σx2 - Σx^2)
    r2  = (n * Σxy - Σx * Σy)^2/((n * Σx2 - Σx^2)*(n * Σy2 - Σy^2))
    return a, b, r2
end #end slope
end #end module
