module PK

using DataFrames

import ..NCA, ..LOG2, ..AbstractData, ..DataSet, ..DataSort, Base.push!, Base.getindex

abstract type AbstractSubject <: AbstractData end

export nca

struct KelData <: AbstractData
    s::Array{Int, 1}
    e::Array{Int, 1}
    a::Array{Float64, 1}
    b::Array{Float64, 1}
    r::Array{Float64, 1}
    function KelData(s, e, a, b, r)::KelData
        new(s, e, a, b, r)::KelData
    end
    function KelData()::KelData
        new(Array{Int, 1}(undef, 0), Array{Int, 1}(undef, 0), Array{Float64, 1}(undef, 0), Array{Float64, 1}(undef, 0), Array{Float64, 1}(undef, 0))::KelData
    end
end
function Base.push!(keldata::KelData, s, e, a, b, r)
    push!(keldata.s, s)
    push!(keldata.e, e)
    push!(keldata.a, a)
    push!(keldata.b, b)
    push!(keldata.r, r)
end

mutable struct ElimRange
    kelstart::Int
    kelend::Int
    kelexcl::Array{Int, 1}
    function ElimRange(kelstart, kelend, kelexcl)
        new(kelstart, kelend, kelexcl)::ElimRange
    end
    function ElimRange()
        new(0, 0, Array{Int, 1}(undef, 0))::ElimRange
    end
end

mutable struct DoseTime
    dose::Real
    time::Real
    tau::Real
    function DoseTime(dose, time, tau)
        new(dose, time, tau)::DoseTime
    end
    function DoseTime()
        new(NaN, NaN, NaN)::DoseTime
    end
end


mutable struct PKSubject <: AbstractSubject
    time::Array
    conc::Array
    kelauto::Bool
    kelrange::ElimRange
    dosetime::DoseTime
    keldata::KelData
    sort::DataSort
    function PKSubject(time::Vector, conc::Vector, kelauto::Bool, kelrange::ElimRange, dosetime::DoseTime, keldata::KelData; sort = DataSort())
        new(time, conc, kelauto, kelrange, dosetime, keldata, sort)::PKSubject
    end
    function PKSubject(time::Vector, conc::Vector, kelauto::Bool, kelrange, dosetime; sort = DataSort())
        new(time, conc, kelauto, kelrange, dosetime, KelData(), sort)::PKSubject
    end
    function PKSubject(time::Vector, conc::Vector; sort = DataSort())
        new(time, conc, true, ElimRange(), DoseTime(), KelData(), sort)::PKSubject
    end
    function PKSubject(data::DataFrame; time::Symbol, conc::Symbol, timesort::Bool = false, sort = DataSort())
        if timesort sort!(data, time) end
        #Check double time
        new(copy(data[!,time]), copy(data[!,conc]), true, ElimRange(), DoseTime(), KelData(), sort)::PKSubject
    end
end

mutable struct PDSubject <: AbstractSubject
    time::Array
    resp::Array
    bl::Real
    th::Real
    sort::DataSort
    function PDSubject(time, resp, bl, th; sort = DataSort())
        new(time, resp, bl, th, sort)::PDSubject
    end
    function PDSubject(time, resp; sort = DataSort())
        new(time, resp, 0, NaN, sort)::PDSubject
    end
end

struct PKPDProfile{T <: AbstractSubject} <: AbstractData
    subject::T
    method
    result::Dict
    function PKPDProfile(subject::T, result; method = nothing) where T <: AbstractSubject
        new{T}(subject, method, result)
    end
end

struct LimitRule
    lloq::Float64
    btmax::Float64
    atmax::Float64
    function LimitRule(lloq, btmax, atmax)
        new(lloq, btmax, atmax)::LimitRule
    end
    function LimitRule()
        new(NaN, NaN, NaN)::LimitRule
    end
end

function Base.getindex(a::PKPDProfile{T}, s::Symbol)::Real where T <: AbstractSubject
    return a.result[s]
end
function Base.getindex(a::DataSet{PKPDProfile{T}}, i::Int64) where T <: AbstractSubject
    return a.data[i]
end
function Base.getindex(a::DataSet{PKPDProfile{T}}, i::Int64, s::Symbol)::Real where T <: AbstractSubject
    return a.data[i].result[s]
end
function Base.getindex(a::DataSet{T}, i::Int64) where T <: AbstractSubject
    return a.data[i]
end
using DataFrames
    #Makoid C, Vuchetich J, Banakar V. 1996-1999. Basic Pharmacokinetics.
    function nca(data::DataFrame; conc = NaN, effect = NaN, time=:Time, sort = NaN, calcm = :lint, bl::Float64 = NaN, th::Float64 = NaN)::NCA
        columns  = DataFrames.names(data)
        errorlog = ""
        errorn   = 0
        errors   = Int[]
        if typeof(findfirst(x ->  x  == conc, columns)) == Nothing && typeof(findfirst(x ->  x  == effect, columns)) == Nothing
            errorn+=1
            errorlog *= string(errorn)*": Concentration column not found\n"
            push!(errors, 1)
        end
        if typeof(findfirst(x ->  x  == time, columns)) == Nothing
            errorn+=1
            errorlog *= string(errorn)*": Time column not found\n"
            push!(errors, 2)
        end
        if  !(calcm == :lint || calcm == :logt || calcm == :luld)
            errorn+=1
            errorlog *= string(errorn)*": Calculation method not found\n"
            push!(errors, 3)
        end
        if sort !== NaN
            if isa(sort, Array)
                for i = 1:length(sort)
                    if typeof(findfirst(x ->  x  == sort[i], columns)) == Nothing
                        errorn+=1
                        errorlog *= string(errorn)*": Sort column not found\n"
                        push!(errors, 4)
                    end
                end
            else
                errorn+=1
                errorlog *= string(errorn)*": Sort is not array\n"
                push!(errors, 5)
            end
        end
        if errorn > 0 return NCA(DataFrame(), DataFrame(), DataFrame(), "", errorlog, errors) end
        if isa(conc, String)  conc = Symbol(conc) end
        if isa(time, String)  time = Symbol(time) end

        res  = DataFrame()
        elim = DataFrame()

        if sort === NaN
            if effect !== NaN
                res   = ncapd(data; conc = effect, time = time, calcm = calcm, bl = bl, th = th)
            elseif conc !== NaN
                res, rsqn, elim = nca_(data, conc, time, calcm)
            end
        else
            #PK NCA
            if conc !== NaN
                sortlist = unique(data[:, sort])
                res  = DataFrame(AUClast = Float64[], Cmax = Float64[], Tmax = Float64[], AUMClast = Float64[], MRTlast = Float64[], Kel = Float64[], HL = Float64[], Rsq = Float64[], AUCinf = Float64[], AUCpct = Float64[])
                elim = DataFrame(Start = Int[], End = Int[], b = Float64[], a = Float64[], Rsq = Float64[])
                for i = 1:size(sortlist, 1) #For each line in sortlist
                    datai = DataFrame(Concentration = Float64[], Time = Float64[])
                    for c = 1:size(data, 1) #For each line in data
                        if data[c, sort] == sortlist[i,:]
                            push!(datai, [data[c, conc], data[c, time]])
                        end
                    end
                    ncares = nca_(datai, :Concentration, :Time, calcm)
                    append!(res, ncares[1])
                    if ncares[2] !== NaN
                        append!(elim, DataFrame(ncares[3][ncares[2], :]))
                    else
                        append!(elim, DataFrame(Start = [0], End = [0], b = [NaN], a = [NaN], Rsq = [NaN]))
                    end
                end
                elim = hcat(sortlist, elim)
            #PD NCA
            elseif effect !== NaN && bl !== NaN
                sortlist = unique(data[:, sort])
                res      = DataFrame(RMAX = Float64[], TH = Float64[], BL = Float64[], AUCABL = Float64[], AUCBBL = Float64[], AUCATH = Float64[], AUCBTH = Float64[], AUCBLNET = Float64[], AUCTHNET = Float64[], AUCDBLTH = Float64[], TABL = Float64[], TBBL = Float64[], TATH = Float64[], TBTH = Float64[])
                for i = 1:size(sortlist, 1) #For each line in sortlist
                    datai = DataFrame(Concentration = Float64[], Time = Float64[])
                    for c = 1:size(data, 1) #For each line in data
                        if data[c, sort] == sortlist[i,:]
                            #println(conc, time)
                            push!(datai, [data[c, effect], data[c, time]])
                        end
                    end
                    pdres  = ncapd(datai; conc = :Concentration, time = :Time, calcm = calcm, bl = bl, th = th)
                    append!(res, pdres)
                end
            end
            res  = hcat(sortlist, res)
        end
        return NCA(res, elim, DataFrame(), "", "", [0])
    end

    function nca_(data::DataFrame, conc::Symbol, time::Symbol, calcm = :lint)
        if length(unique(data[:,time])) != length(data[:,time])
            return  DataFrame(AUClast = [NaN], Cmax = [NaN], Tmax = [NaN], AUMClast = [NaN], MRTlast = [NaN], Kel = [NaN], HL = [NaN], Rsq = [NaN], AUCinf = [NaN], AUCpct = [NaN])
        end
        ncares = Dict(:Nums => 0, :Cmax => NaN, :Tmax => NaN, :Tmaxn => 0, :Tlast => NaN, :Clast => NaN, :AUClast => NaN, :AUMClast => NaN, :MRTlast => NaN, :Kel => NaN, :HL => NaN, :Rsq => NaN, :AUCinf => NaN, :AUCpct => NaN)

        pklog = "NCA analysis: \n"
        pklog *= "Concentration column: "*string(conc)*"\n"
        pklog *= "Time column: "*string(time)*"\n"
        #sort!(data, [time])
        ncarule!(data, conc, time, LimitRule(eps(), 0, NaN))

        n::Int     = nrow(data)
        cmax = maximum(data[:, conc])               #Cmax
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
        rsqn  = NaN
        aucinf = NaN
        aucinfpct = NaN
        mrtlast = NaN
    # --- Tmax ---
        cmax, tmax, tmaxn = ctmax(data, conc, time)
        # --- AUC AUMC ---
        for i = 2:nrow(data)
            if calcm == :lint
                auc_  = linauc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])    # (data[i - 1, conc] + data[i, conc]) * (data[i, time] - data[i - 1, time])/2
                aumc_ = linaumc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])  # (data[i - 1, conc]*data[i - 1, time] + data[i, conc]*data[i, time]) * (data[i, time] - data[i - 1, time])/2
            elseif calcm == :logt
                if i > tmaxn
                    if data[i, conc] < data[i - 1, conc] &&  data[i, conc] > 0 &&  data[i - 1, conc] > 0
                        auc_  = logauc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                        aumc_ = logaumc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                    else
                        auc_  = linauc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                        aumc_ = linaumc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                    end
                else
                    auc_ = linauc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                    aumc_ = linaumc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                end
            elseif calcm == :luld
                if data[i, conc] < data[i - 1, conc] &&  data[i, conc] > 0 &&  data[i - 1, conc] > 0
                    auc_  = logauc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                    aumc_ = logaumc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                else
                    auc_  = linauc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                    aumc_ = linaumc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                end
            end
            auc   += auc_
            aumc  += aumc_
            pklog *= "AUC("*string(i-1)*") = "*string(auc_)*"; Sum = "*string(auc)*"\n"
        end
        # --- MRT ---
        mrt    = aumc / auc
        pklog *= "AUClast = "*string(auc)*"\n"

        # --- Kel, HL, ets
        if n - tmaxn >= 3
            keldata = DataFrame(Start = Int[], End = Int[], b = Float64[], a = Float64[], Rsq = Float64[])
            logconc = log.(data[:, conc])
            for i::Int = tmaxn+1:n-2
                sl = slope(data[i:n,time], logconc[i:n])
                if sl[2] < 0
                    push!(keldata, [i, n, sl[1], sl[2], sl[3]])
                end
            end
            if nrow(keldata) > 0
                rsq, rsqn = findmax(keldata[:, :Rsq])
                kel = abs(keldata[rsqn,:a])
                hl  = LOG2/kel
                aucinf = auc + clast/kel
                aucinfpct = (aucinf - auc) / aucinf * 100.0
            end
        else
            keldata = nothing
        end
        # arsq = 1 - (1-rsq)*(rsqn-1)/(rsqn-2)
        return DataFrame(AUClast = [auc], Cmax = [cmax], Tmax = [tmax], AUMClast = [aumc], MRTlast = [mrt], Kel = [kel], HL = [hl], Rsq = [rsq], AUCinf = [aucinf], AUCpct = [aucinfpct]), rsqn, keldata
    end

    function ncarule!(data::DataFrame, conc::Symbol, time::Symbol, rule::LimitRule)
        sort!(data, [time])

        cmax, tmax, tmaxn = ctmax(data, conc, time)
        filterv = Array{Int, 1}(undef, 0)
        for i = 1:nrow(data)
            if data[i, conc] < rule.lloq
                if i <= tmaxn
                    data[i, conc] = rule.btmax
                else
                    data[i, conc] = rule.atmax
                end
            end
        end
        for i = 1:nrow(data)
            if data[i, conc] === NaN push!(filterv, i) end
        end
        if length(filterv) > 0
            deleterows!(data, filterv)
        end
    end

    """
        cmax
        tmax
        tmaxn
    """
    function ctmax(data::DataFrame, conc::Symbol, time::Symbol)
        cmax  = maximum(data[!, conc])
        tmax  = NaN
        tmaxn = 0
        for i = 1:nrow(data)
            if data[i, conc] == cmax
                tmax  = data[i, time]
                tmaxn = i
                break
            end
        end
        return cmax, tmax, tmaxn
    end
    function ctmax(data::PKSubject)
        cmax  = maximum(data.conc)
        tmax  = NaN
        tmaxn = 0
        for i = 1:obsnum(data)
            if data.conc[i] == cmax
                tmax  = data.time[i]
                tmaxn = i
                break
            end
        end
        return cmax, tmax, tmaxn
    end

    function ncapd(data::DataFrame; conc::Symbol, time::Symbol, calcm::Symbol = :lint, bl::Real, th::Real)

        aucabl = 0
        aucbbl = 0
        aucath = 0
        aucbth = 0
        aucdblth = 0
        tabl     = 0
        tbbl     = 0
        tath     = 0
        tbth     = 0
        for i = 2:nrow(data)
            #BASELINE
            if data[i - 1, conc] <= bl && data[i, conc] <= bl
                tbbl   += data[i, time] - data[i - 1, time]
                aucbbl += linauc(data[i - 1, time], data[i, time], bl - data[i - 1, conc], bl - data[i, conc])
            elseif data[i - 1, conc] <= bl &&  data[i, conc] > bl
                tx      = linpredict(data[i - 1, conc], data[i, conc], bl, data[i - 1, time], data[i, time])
                tbbl   += tx - data[i - 1, time]
                tabl   += data[i, time] - tx
                aucbbl += (tx - data[i - 1, time]) * (bl - data[i - 1, conc]) / 2
                aucabl += (data[i, time] - tx) * (data[i, conc] - bl) / 2 #Ok
            elseif data[i - 1, conc] > bl &&  data[i, conc] <= bl
                tx      = linpredict(data[i - 1, conc], data[i, conc], bl, data[i - 1, time], data[i, time])
                tbbl   += data[i, time] - tx
                tabl   += tx - data[i - 1, time]
                aucbbl += (data[i, time] - tx) * (bl - data[i, conc]) / 2
                aucabl += (tx - data[i - 1, time]) * (data[i - 1, conc] - bl) / 2
            elseif data[i - 1, conc] > bl &&  data[i, conc] > bl
                tabl   += data[i, time] - data[i - 1, time]
                aucabl += linauc(data[i - 1, time], data[i, time], data[i - 1, conc] - bl, data[i, conc] - bl)
            end
            #THRESHOLD
            if th !== NaN
                if data[i - 1, conc] <= th && data[i, conc] <= th
                    tbth   += data[i, time] - data[i - 1, time]
                    aucbth += linauc(data[i - 1, time], data[i, time], th - data[i - 1, conc], th - data[i, conc])
                elseif data[i - 1, conc] <= th &&  data[i, conc] > th
                    tx      = linpredict(data[i - 1, conc], data[i, conc], th, data[i - 1, time], data[i, time])
                    tbth   += tx - data[i - 1, time]
                    tath   += data[i, time] - tx
                    aucbth += (tx - data[i - 1, time]) * (th - data[i - 1, conc]) / 2
                    aucath += (data[i, time] - tx) * (data[i, conc] - th) / 2 #Ok
                elseif data[i - 1, conc] > th &&  data[i, conc] <= th
                    tx      = linpredict(data[i - 1, conc], data[i, conc], th, data[i - 1, time], data[i, time])
                    tbth   += data[i, time] - tx
                    tath   += tx - data[i - 1, time]
                    aucbth += (data[i, time] - tx) * (th - data[i, conc]) / 2
                    aucath += (tx - data[i - 1, time]) * (data[i - 1, conc] - th) / 2
                elseif data[i - 1, conc] > th &&  data[i, conc] > th
                    tath   += data[i, time] - data[i - 1, time]
                    aucath += linauc(data[i - 1, time], data[i, time], data[i - 1, conc] - th, data[i, conc] - th)
                end
            end
        end
        if bl > th
            aucdblth = aucath - aucabl
        elseif bl < th
            aucdblth = aucabl - aucath
        end
        return DataFrame(RMAX = [maximum(data[!, conc])], TH = [th], BL = [bl], AUCABL = [aucabl], AUCBBL = [aucbbl], AUCATH = [aucath], AUCBTH = [aucbth], AUCBLNET = [aucabl-aucbbl], AUCTHNET = [aucath-aucbth], AUCDBLTH = [aucdblth], TABL = [tabl], TBBL = [tbbl], TATH = [tath], TBTH = [tbth])
    end

    function slope(x::Vector, y::Vector)::Tuple{Float64, Float64, Float64}
        if length(x) != length(y) throw(ArgumentError("Unequal vector length!")) end
        n   = length(x)
        Σxy = sum(x .* y)
        Σx  = sum(x)
        Σy  = sum(y)
        Σx2 = sum(x .* x)
        Σy2 = sum(y .* y)

        a   = (Σy * Σx2 - Σx * Σxy)/(n * Σx2 - Σx^2)
        b   = (n * Σxy - Σx * Σy)/(n * Σx2 - Σx^2)
        r2  = (n * Σxy - Σx * Σy)^2/((n * Σx2 - Σx^2)*(n * Σy2 - Σy^2))
        return a, b, r2
    end #end slope
        #Linear trapezoidal auc
    function linauc(t₁, t₂, c₁, c₂)::Float64
        return (t₂-t₁)*(c₁+c₂)/2
    end
        #Linear trapezoidal aumc
    function linaumc(t₁, t₂, c₁, c₂)::Float64
        return (t₂-t₁)*(t₁*c₁+t₂*c₂)/2
    end
        #Log trapezoidal auc
    function logauc(t₁, t₂, c₁, c₂)::Float64
        return  (t₂-t₁)*(c₂-c₁)/log(c₂/c₁)
    end
        #Log trapezoidal aumc
    function logaumc(t1, t2, c1, c2)::Float64
        return (t2-t1) * (t2*c2-t1*c1) / log(c2/c1) - (t2-t1)^2 * (c2-c1) / log(c2/c1)^2
    end

    function linpredict(a1::Real, a2::Real, ax::Real, b1::Real, b2::Real)::Float64
        return abs((ax - a1) / (a2 - a1))*(b2 - b1) + b1
    end

    function logtpredict(c1::Real, c2::Real, cx::Real, t1::Real, t2::Real)::Float64
        return log(cx/c1)/log(c2/c1)*(t2-t1)+t1
    end

    function logcpredict(t1::Real, t2::Real, tx::Real, c1::Real, c2::Real)::Float64
        return exp(log(c1) + abs((tx-t1)/(t2-t1))*(log(c2) - log(c1)))
    end


    #---------------------------------------------------------------------------
    function obsnum(data::T) where T <:  AbstractSubject
        return length(data.time)
    end
    function  obsnum(keldata::KelData)
        return length(keldata.a)
    end

    function nca!(data::PKSubject; calcm = :lint)
        result           = Dict(:Obsnum => 0, :Cmax => 0, :Tmax => 0, :Tmaxn => 0, :Tlast => 0, :Clast => 0, :AUClast => 0, :AUMClast => 0, :MRTlast => 0, :Kel => NaN, :HL => NaN, :Rsq => NaN, :Rsqn => NaN, :AUCinf => NaN, :AUCpct => NaN)
        result[:Obsnum]  = obsnum(data)
        result[:Cmax]    = maximum(data.conc)
        result[:Tlast]   = data.time[result[:Obsnum]]
        result[:Clast]   = data.conc[result[:Obsnum]]
        result[:Cmax], result[:Tmax], result[:Tmaxn] = ctmax(data)
        for i = 2:result[:Obsnum]
            if calcm == :lint
                auc   =  linauc(data.time[i - 1], data.time[i], data.conc[i - 1], data.conc[i])    # (data[i - 1, conc] + data[i, conc]) * (data[i, time] - data[i - 1, time])/2
                aumc  = linaumc(data.time[i - 1], data.time[i], data.conc[i - 1], data.conc[i])    # (data[i - 1, conc]*data[i - 1, time] + data[i, conc]*data[i, time]) * (data[i, time] - data[i - 1, time])/2
            elseif calcm == :logt
                if i > result[:Tmaxn]
                    if data.conc[i] < data.conc[i - 1] &&  data.conc[i] > 0 &&  data.conc[i - 1] > 0
                        auc   =  logauc(data.time[i - 1], data.time[i], data.conc[i - 1], data.conc[i])
                        aumc  = logaumc(data.time[i - 1], data.time[i], data.conc[i - 1], data.conc[i])
                    else
                        auc   =  linauc(data.time[i - 1], data.time[i], data.conc[i - 1], data.conc[i])
                        aumc  = linaumc(data.time[i - 1], data.time[i], data.conc[i - 1], data.conc[i])
                    end
                else
                    auc   =  linauc(data.time[i - 1], data.time[i], data.conc[i - 1], data.conc[i])
                    aumc  = linaumc(data.time[i - 1], data.time[i], data.conc[i - 1], data.conc[i])
                end
            elseif calcm == :luld
                if data.conc[i] < data.conc[i - 1] &&  data.conc[i] > 0 &&  data.conc[i - 1] > 0
                    auc   =  logauc(data.time[i - 1], data.time[i], data.conc[i - 1], data.conc[i])
                    aumc  = logaumc(data.time[i - 1], data.time[i], data.conc[i - 1], data.conc[i])
                else
                    auc   =  linauc(data.time[i - 1], data.time[i], data.conc[i - 1], data.conc[i])
                    aumc  = linaumc(data.time[i - 1], data.time[i], data.conc[i - 1], data.conc[i])
                end
            end
            result[:AUClast]  += auc
            result[:AUMClast] += aumc
        end
        result[:MRTlast]       = result[:AUMClast] / result[:AUClast]
        keldata                = KelData(Int[], Int[], Float64[], Float64[], Float64[])
        if data.kelauto
            if result[:Obsnum] - result[:Tmaxn] >= 3
                logconc = log.(data.conc)
                cmask   = Array{Bool, 1}(undef, result[:Obsnum])
                for i  = result[:Tmaxn] + 1:result[:Obsnum] - 2
                    cmask                    .= false
                    cmask[i:result[:Obsnum]] .= true
                    sl = slope(data.time[cmask], logconc[cmask])
                    if sl[2] < 0
                        push!(keldata, i, result[:Obsnum], sl[2], sl[1], sl[3])
                    end
                end
            end
        else
            logconc = log.(data.conc)
            cmask   = Array{Bool, 1}(undef, result[:Obsnum])
            cmask  .= false
            cmask[data.kelrange.kelstart:data.kelrange.kelend] .= true
            cmask[data.kelrange.kelexcl] .= false
            sl = slope(data.time[cmask], logconc[cmask])
            push!(keldata, data.kelrange.kelstart, data.kelrange.kelend, sl[2], sl[1], sl[3])
        end
        if  obsnum(keldata) > 0
            result[:Rsq], result[:Rsqn] = findmax(keldata.r)
            data.kelrange.kelstart   = keldata.s[result[:Rsqn]]
            data.kelrange.kelend     = keldata.e[result[:Rsqn]]
            data.keldata    = keldata
            result[:Kel]    = abs(keldata.a[result[:Rsqn]])
            result[:HL]     = LOG2 / result[:Kel]
            result[:AUCinf] = result[:AUClast] + result[:Clast] / result[:Kel]
            result[:AUCpct] = (result[:AUCinf] - result[:AUClast]) / result[:AUCinf] * 100.0
        end
        return PKPDProfile(data, result; method = calcm)
    end
    function nca!(data::DataSet{PKSubject}; calcm = :lint)
        results = Array{PKPDProfile, 1}(undef, 0)
        for i = 1:length(data.data)
            push!(results, nca!(data.data[i]; calcm = calcm))
        end
        return DataSet(Tuple(results))
    end
    #---------------------------------------------------------------------------
    function nca!(data::PDSubject)::PKPDProfile{PDSubject}
        result = Dict(:Obsnum => 0, :RMAX => 0, :TH => 0, :BL => 0, :AUCABL => 0, :AUCBBL => 0, :AUCATH => NaN, :AUCBTH => NaN, :AUCBLNET => 0, :AUCTHNET => 0, :AUCDBLTH => NaN, :TABL => 0, :TBBL => 0, :TATH => NaN, :TBTH => NaN)
        result[:Obsnum] = obsnum(data)
        result[:TH] = data.th
        result[:BL] = data.bl
        for i = 2:obsnum(data) #BASELINE
            if data.resp[i - 1] <= result[:BL] && data.resp[i] <= result[:BL]
                result[:TBBL]   += data.time[i,] - data.time[i - 1]
                result[:AUCBBL] += linauc(data.time[i - 1], data.time[i], result[:BL] - data.resp[i - 1], result[:BL] - data.resp[i])
            elseif data.resp[i - 1] <= result[:BL] &&  data.resp[i] > result[:BL]
                tx      = linpredict(data.resp[i - 1], data.resp[i], result[:BL], data.time[i - 1], data.time[i])
                result[:TBBL]   += tx - data.time[i - 1]
                result[:TABL]   += data.time[i] - tx
                result[:AUCBBL] += (tx - data.time[i - 1]) * (result[:BL] - data.resp[i - 1]) / 2
                result[:AUCABL] += (data.time[i] - tx) * (data.resp[i] - result[:BL]) / 2 #Ok
            elseif data.resp[i - 1] > result[:BL] &&  data.resp[i] <= result[:BL]
                tx      = linpredict(data.resp[i - 1], data.resp[i], result[:BL], data.time[i - 1], data.time[i])
                result[:TBBL]   += data.time[i] - tx
                result[:TABL]   += tx - data.time[i - 1]
                result[:AUCBBL] += (data.time[i] - tx) * (result[:BL] - data.resp[i]) / 2
                result[:AUCABL] += (tx - data.time[i - 1]) * (data.resp[i - 1] - result[:BL]) / 2
            elseif data.resp[i - 1] > result[:BL] &&  data.resp[i] > result[:BL]
                result[:TABL]   += data.time[i] - data.time[i - 1]
                result[:AUCABL]     += linauc(data.time[i - 1], data.time[i], data.resp[i - 1] - result[:BL], data.resp[i] - result[:BL])
            end
        end #BASELINE

        #THRESHOLD
        if result[:TH] !== NaN
            result[:AUCATH]   = 0
            result[:AUCBTH]   = 0
            result[:TATH]     = 0
            result[:TBTH]     = 0
            result[:AUCTHNET] = 0
            result[:AUCDBLTH] = 0
            for i = 2:obsnum(data)
                if data.resp[i - 1] <= result[:TH] && data.resp[i] <= result[:TH]
                    result[:TBTH]   += data.time[i] - data.time[i - 1]
                    result[:AUCBTH] += linauc(data.time[i - 1], data.time[i], result[:TH] - data.resp[i - 1], result[:TH] - data.resp[i])
                elseif data.resp[i - 1] <= result[:TH] &&  data.resp[i] > result[:TH]
                    tx      = linpredict(data.resp[i - 1], data.resp[i], result[:TH], data.time[i - 1], data.time[i])
                    result[:TBTH]   += tx - data.time[i - 1]
                    result[:TATH]   += data.time[i] - tx
                    result[:AUCBTH] += (tx - data.time[i - 1]) * (result[:TH] - data.resp[i - 1]) / 2
                    result[:AUCATH] += (data.time[i] - tx) * (data.resp[i] - result[:TH]) / 2 #Ok
                elseif data.resp[i - 1] > result[:TH] &&  data.resp[i] <= result[:TH]
                    tx      = linpredict(data.resp[i - 1], data.resp[i], result[:TH], data.time[i - 1], data.time[i])
                    result[:TBTH]   += data.time[i] - tx
                    result[:TATH]   += tx - data.time[i - 1]
                    result[:AUCBTH] += (data.time[i] - tx) * (result[:TH] - data.resp[i]) / 2
                    result[:AUCATH] += (tx - data.time[i - 1]) * (data.resp[i - 1] - result[:TH]) / 2
                elseif data.resp[i - 1] > result[:TH] &&  data.resp[i] > result[:TH]
                    result[:TATH]   += data.time[i] - data.time[i - 1]
                    result[:AUCATH] += linauc(data.time[i - 1], data.time[i], data.resp[i - 1] - result[:TH], data.resp[i] - result[:TH])
                end
            end
            if result[:BL] > result[:TH]
                result[:AUCDBLTH] = result[:AUCATH] - result[:AUCABL]
            else
                result[:AUCDBLTH] = result[:AUCABL] -result[:AUCATH]
            end
            result[:AUCTHNET] = result[:AUCATH] - result[:AUCBTH]
        end
        result[:AUCBLNET] = result[:AUCABL] - result[:AUCBBL]
        return PKPDProfile(data, result)
    end

    function nca!(data::DataSet{PDSubject})
        results = Array{PKPDProfile, 1}(undef, 0)
        for i = 1:length(data.data)
            push!(results, nca!(data.data[i]))
        end
        return DataSet(Tuple(results))
    end
    #---------------------------------------------------------------------------
    function pkimport(data::DataFrame, sort::Array, rule::LimitRule; conc::Symbol, time::Symbol)
        sortlist = unique(data[:, sort])
        results  = Array{PKSubject, 1}(undef, 0)
        for i = 1:size(sortlist, 1) #For each line in sortlist
            datai = DataFrame(Time = Float64[], Conc = Float64[])
            for c = 1:size(data, 1) #For each line in data
                if data[c, sort] == sortlist[i,:]
                    push!(datai, [data[c, time], data[c, conc]])
                end
            end
            if rule.lloq !== NaN
                ncarule!(datai, :Conc, :Time, rule)
            else
                sort!(datai, :Time)
            end
            push!(results, PKSubject(datai[!, :Time], datai[!, :Conc], sort = DataSort(sort, collect(sortlist[i,:]))))
        end
        return DataSet(Tuple(results))
    end

    function pkimport(data::DataFrame, sort::Array; time::Symbol, conc::Symbol)
        rule = LimitRule()
        return pkimport(data, sort, rule; conc = conc, time = time)
    end

    function pkimport(data::DataFrame, rule::LimitRule; time::Symbol, conc::Symbol)
        rule = LimitRule()
        datai = ncarule!(copy(data[!,[time, conc]]), conc, time, rule)
        return DataSet(tuple(PKSubject(datai[!, time], datai[!, conc])))
    end
    function pkimport(data::DataFrame; time::Symbol, conc::Symbol)
        datai = sort(data[!,[time, conc]], time)
        return DataSet(tuple(PKSubject(datai[!, time], datai[!, conc])))
    end
    #---------------------------------------------------------------------------
    function pdimport(data::DataFrame, sort::Array; resp::Symbol, time::Symbol, bl = 0, th = NaN)
        sortlist = unique(data[:, sort])
        results  = Array{PDSubject, 1}(undef, 0)
        for i = 1:size(sortlist, 1) #For each line in sortlist
            datai = DataFrame(Time = Float64[], Resp = Float64[])
            for c = 1:size(data, 1) #For each line in data
                if data[c, sort] == sortlist[i,:]
                    push!(datai, [data[c, time], data[c, resp]])
                end
            end
            sort!(datai, :Time)
            push!(results, PDSubject(datai[!, :Time], datai[!, :Resp], bl, th, sort = DataSort(sort, collect(sortlist[i,:]))))
        end
        return DataSet(Tuple(results))
    end
    function pdimport(data::DataFrame; resp::Symbol, time::Symbol, bl = 0, th = NaN)
        datai = sort(data[!,[time, resp]], time)
        return DataSet(tuple(PDSubject(datai[!, time], datai[!, resp], bl, th)))
    end
    #---------------------------------------------------------------------------
    function ncarule!(data::PKSubject, rule::LimitRule)
    end
    function ncarule!(data::DataSet{PKSubject}, rule::LimitRule)
    end
    function applyncarule!(data::PKPDProfile{PKSubject}, rule::LimitRule)
    end
    function applyncarule!(data::DataSet{PKPDProfile{PKSubject}}, rule::LimitRule)
    end
    function setelimrange!(data::PKSubject, range::ElimRange)
    end
    function setelimrange!(data::DataSet{PKSubject}, range::ElimRange)
    end
    function setelimrange!(data::DataSet{PKSubject}, subj::Array{Int,1}, range::ElimRange)
    end
    function setelimrange!(data::DataSet{PKSubject}, subj::Int, range::ElimRange)
    end
    function applyelimrange!(data::PKPDProfile{PKSubject}, range::ElimRange)
    end
    function applyelimrange!(data::DataSet{PKPDProfile{PKSubject}}, range::ElimRange)
    end
    function applyelimrange!(data::DataSet{PKPDProfile{PKSubject}}, subj::Array{Int,1}, range::ElimRange)
    end
    function applyelimrange!(data::DataSet{PKPDProfile{PKSubject}}, subj::Int, range::ElimRange)
    end
    function setdosetime!(data::PKSubject, dosetime::DoseTime)
    end
    function setdosetime!(data::DataSet{PKSubject}, dosetime::DoseTime)
    end

    function setth!(data::DataSet{PDSubject}, th)
        for i = 1:length(data)
            data[1].th = th
        end
        data
    end
    function setbl!(data::DataSet{PDSubject}, bl)
        for i = 1:length(data)
            data[1].bl = bl
        end
        data
    end

end #end module
