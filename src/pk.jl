#Pharmacokinetics
#Makoid C, Vuchetich J, Banakar V. 1996-1999. Basic Pharmacokinetics.

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
#!!!
struct KelData <: AbstractData
    s::Array{Int, 1}
    e::Array{Int, 1}
    a::Array{Number, 1}
    b::Array{Number, 1}
    r::Array{Number, 1}
    ar::Array{Number, 1}
    function KelData(s, e, a, b, r, ar)::KelData
        new(s, e, a, b, r, ar)::KelData
    end
    function KelData()::KelData
        new(Int[], Int[], Real[], Real[], Real[], Real[])::KelData
    end
end
function Base.push!(keldata::KelData, s, e, a, b, r, ar)
    push!(keldata.s, s)
    push!(keldata.e, e)
    push!(keldata.a, a)
    push!(keldata.b, b)
    push!(keldata.r, r)
    push!(keldata.ar, ar)
end

mutable struct ElimRange
    kelstart::Int
    kelend::Int
    kelexcl::Vector{Int}
    function ElimRange(kelstart::Int, kelend::Int, kelexcl::Vector{Int})::ElimRange
        if kelstart < 0 throw(ArgumentError("Kel start point < 0")) end
        if kelend   < 0 throw(ArgumentError("Kel endpoint < 0")) end
        if any(x -> x < 0, kelexcl) throw(ArgumentError("Exclude point < 0")) end
        new(kelstart, kelend, kelexcl)::ElimRange
    end
    function ElimRange(kelstart::Int, kelend::Int)
        ElimRange(kelstart, kelend, Vector{Int}(undef, 0))
    end
    function ElimRange()
        ElimRange(0, 0, Vector{Int}(undef, 0))
    end
end

struct DoseTime
    dose::Number
    time::Number
    tau::Number
    function DoseTime(;dose::Number = NaN, time::Number = 0, tau::Number = NaN)
        new(dose, time, tau)::DoseTime
    end
    function DoseTime(dose)
        new(dose, 0, NaN)::DoseTime
    end
    function DoseTime(dose, time)
        new(dose, time, NaN)::DoseTime
    end
    function DoseTime(dose, time, tau)
        new(dose, time, tau)::DoseTime
    end
end


mutable struct PKSubject <: AbstractSubject
    time::Array
    obs::Array
    kelauto::Bool
    kelrange::ElimRange
    dosetime::DoseTime
    keldata::KelData
    sort::Dict
    function PKSubject(time::Vector, conc::Vector, kelauto::Bool, kelrange::ElimRange, dosetime::DoseTime, keldata::KelData; sort = Dict())
        new(time, conc, kelauto, kelrange, dosetime, keldata, sort)::PKSubject
    end
    function PKSubject(time::Vector, conc::Vector, kelauto::Bool, kelrange, dosetime; sort = Dict())
        new(time, conc, kelauto, kelrange, dosetime, KelData(), sort)::PKSubject
    end
    function PKSubject(time::Vector, conc::Vector; sort = Dict())
        new(time, conc, true, ElimRange(), DoseTime(NaN, 0 * time[1]), KelData(), sort)::PKSubject
    end
    function PKSubject(data::DataFrame; time::Symbol, conc::Symbol, timesort::Bool = false, sort = Dict())
        if timesort sort!(data, time) end
        #Check double time
        new(copy(data[!,time]), copy(data[!,conc]), true, ElimRange(), DoseTime(), KelData(), sort)::PKSubject
    end
end

mutable struct PDSubject <: AbstractSubject
    time::Array
    obs::Array
    bl::Real
    th::Real
    sort::Dict
    function PDSubject(time, resp, bl, th; sort = Dict())
        new(time, resp, bl, th, sort)::PDSubject
    end
    function PDSubject(time, resp; sort = Dict())
        new(time, resp, 0, NaN, sort)::PDSubject
    end
end

mutable struct PKPDProfile{T <: AbstractSubject} <: AbstractData
    subject::T
    method
    result::Dict
    sort::Dict
    function PKPDProfile(subject::T, result; method = nothing) where T <: AbstractSubject
        new{T}(subject, method, result, subject.sort)
    end
end

struct LimitRule
    lloq::Real
    btmax::Real
    atmax::Real
    nan::Real
    rm::Bool
    function LimitRule(lloq, btmax, atmax, nan, rm)
        new(lloq, btmax, atmax, nan, rm)::LimitRule
    end
    function LimitRule(lloq, btmax, atmax, nan)
        new(lloq, btmax, atmax, nan, false)::LimitRule
    end
    function LimitRule(lloq, btmax, atmax)
        new(lloq, btmax, atmax, NaN, true)::LimitRule
    end
    function LimitRule()
        new(NaN, NaN, NaN, NaN, false)::LimitRule
    end
end

#=
function Base.getindex(a::PKPDProfile{T}, s::Symbol)::Real where T <: AbstractSubject
    return a.result[s]
end
function Base.getindex(a::DataSet{PKPDProfile}, i::Int)
    return a.data[i]
end
function Base.getindex(a::DataSet{PKPDProfile}, i::Int, s::Symbol)::Real
    return a.data[i].result[s]
end
function Base.getindex(a::DataSet{PKPDProfile}, d::Pair)
    for i = 1:length(a)
        if d ∈ a[i].subject.sort return a[i] end
    end
end
function Base.getindex(a::DataSet{T}, i::Int) where T <: AbstractSubject
    return a.data[i]
end
=#

#-------------------------------------------------------------------------------
function obsnum(data::T) where T <:  AbstractSubject
    return length(data.time)
end
function obsnum(keldata::KelData)
    return length(keldata.a)
end
function length(data::T) where T <:  AbstractSubject
    return length(data.time)
end
function  length(keldata::KelData)
    return length(keldata.a)
end
#-------------------------------------------------------------------------------
function Base.show(io::IO, obj::DataSet{PKSubject})
    println(io, "DataSet: Pharmacokinetic subject")
    println(io, "Length: $(length(obj))")
    for i = 1:length(obj)
        println(io, "Subject $(i): ", obj[i].sort)
    end
end
function Base.show(io::IO, obj::DataSet{PDSubject})
    println(io, "DataSet: Pharmacodynamic subject")
    println(io, "Length: $(length(obj))")
    for i = 1:length(obj)
        println(io, "Subject $(i): ", obj[i].sort)
    end
end
function Base.show(io::IO, obj::DataSet{PKPDProfile})
    println(io, "DataSet: Pharmacokinetic/Pharmacodynamic profile")
    println(io, "Length: $(length(obj))")
    for i = 1:length(obj)
        println(io, "Subject $(i): ", obj[i].subject.sort)
    end
end
function Base.show(io::IO, obj::PKSubject)
    println(io, "Pharmacokinetic subject")
    println(io, "Observations: $(length(obj))")
    println(io,  obj.kelrange)
    println(io, "Time   Concentration")
    for i = 1:length(obj)
        println(io, obj.time[i], " => ", obj.obs[i])
    end
end
function Base.show(io::IO, obj::PDSubject)
    println(io, "Pharmacodynamic subject")
    println(io, "Observations: $(length(obj))")
    println(io, "Time   Responce")
    for i = 1:length(obj)
        println(io, obj.time[i], " => ", obj.obs[i])
    end
end
function Base.show(io::IO, obj::ElimRange)
    print(io, "Elimination range: $(obj.kelstart) - $(obj.kelend) ")
    if length(obj.kelexcl) > 0
        print(io, "Exclusions: $(obj.kelexcl[i])")
        if length(obj.kelexcl) > 1 for i = 1:length(obj.kelexcl) print(io, ", $(obj.kelexcl[i])") end end
        print(io, ".")
    else
        print(io, "No exclusion.")
    end
end
function Base.show(io::IO, obj::KelData)
    m = copy(obj.s)
    m = hcat(m, obj.e)
    m = hcat(m, obj.a)
    m = hcat(m, obj.b)
    m = hcat(m, obj.r)
    m = hcat(m, obj.ar)
    println(io, "Elimination table:")
    print(io, m)
end
#-------------------------------------------------------------------------------
function ncarule!(data::DataFrame, conc::Symbol, time::Symbol, rule::LimitRule)
        sort!(data, [time])

        cmax, tmax, tmaxn = ctmax(data, conc, time)
        #NaN Rule
        if rule.nan !== NaN
            for i = 1:nrow(data)
                if data[i, conc] === NaN data[i, conc] = rule.nan end
            end
        end

        #LLOQ rule
        for i = 1:nrow(data)
            if data[i, conc] < rule.lloq
                if i <= tmaxn
                    data[i, conc] = rule.btmax
                else
                    data[i, conc] = rule.atmax
                end
            end
        end

        #NaN Remove rule
        if rule.rm
            filterv = Array{Int, 1}(undef, 0)
            for i = 1:nrow(data)
                if data[i, conc] === NaN push!(filterv, i) end
            end
            if length(filterv) > 0
                deleterows!(data, filterv)
            end
        end
end
function ncarule!(data::PKSubject, rule::LimitRule)
    #!!!!
end

    """
        cmax
        tmax
        tmaxn
    """
    function ctmax(time::Vector, obs::Vector)
        if length(time) != length(obs) throw(ArgumentError("Unequal vector length!")) end
        cmax  = maximum(obs[obs .!== NaN]);
        tmax  = NaN
        tmaxn = 0
        for i = 1:length(time)
            if obs[i] == cmax
                tmax  = time[i]
                tmaxn = i
                break
            end
        end
        return cmax, tmax, tmaxn
    end
    function ctmax(data::DataFrame, conc::Symbol, time::Symbol)
        return ctmax(data[!, time], data[!, conc])
    end
    function ctmax(data::PKSubject)
        return ctmax(data.time, data.obs)
    end
    function ctmax(data::PKSubject, dosetime)
        s     = 0
        for i = 1:length(data) - 1
            if dosetime >= data.time[i] && dosetime < data.time[i+1] s = i; break end
        end
        if dosetime == data.time[1] return ctmax(data.time, data.obs) end
        mask  = trues(length(data))
        cpredict   = linpredict(data.time[s], data.time[s+1], data.dosetime.time, data.obs[s], data.obs[s+1])
        mask[1:s] .= false
        cmax, tmax, tmaxn = ctmax(data.time[mask], data.obs[mask])
        if cmax > cpredict return cmax, tmax, tmaxn + s else return cpredict, data.dosetime.time, s end
    end
    function ctmax(data::PKSubject, dosetime, tau)
        if tau === NaN || tau <= 0 return ctmax(data, dosetime) end
        tautime = dosetime + tau
        s, e = ncarange(data, dosetime, tau)
        if dosetime == data.time[1] && tautime == data.time[end] return ctmax(data.time, data.obs) end
        mask           = trues(length(data))
        scpredict      = linpredict(data.time[s], data.time[s+1], data.dosetime.time, data.obs[s], data.obs[s+1])
        mask[1:s]     .= false

        if e < length(data)
            ecpredict      = linpredict(data.time[e], data.time[e+1], tautime, data.obs[e], data.obs[e+1])
            mask[e+1:end] .= false
        end
        cmax, tmax, tmaxn = ctmax(data.time[mask], data.obs[mask])
        if e < length(data)
            if cmax > max(scpredict, ecpredict)
                return cmax, tmax, tmaxn + s
            elseif scpredict >= max(cmax, ecpredict)
                return scpredict, data.dosetime.time, s
            else
                return ecpredict, tautime, e
            end
        elseif cmax > scpredict return max, tmax, tmaxn + s else return scpredict, data.dosetime.time, s end
    end
    """
    Range for AUC
        return start point number, end point number
    """
    function ncarange(data::PKSubject, dosetime, tau)
        tautime = dosetime + tau
        s     = 0
        e     = 0
        for i = 1:length(data) - 1
            if dosetime >= data.time[i] && dosetime < data.time[i+1] s = i; break end
        end
        if tautime >= data.time[end]
            e = length(data.time)
        else
            for i = s:length(data) - 1
                if tautime >= data.time[i] && tautime < data.time[i+1] e = i; break end
            end
        end
        return s, e
    end

    #---------------------------------------------------------------------------
    function slope(x::Vector, y::Vector)::Tuple{Number, Number, Number, Number}
        if length(x) != length(y) throw(ArgumentError("Unequal vector length!")) end
        n   = length(x)
        if n < 2 throw(ArgumentError("n < 2!")) end
        Σxy = sum(x .* y)
        Σx  = sum(x)
        Σy  = sum(y)
        Σx2 = sum(x .* x)
        Σy2 = sum(y .* y)
        a   = (n * Σxy - Σx * Σy)/(n * Σx2 - Σx^2)
        b   = (Σy * Σx2 - Σx * Σxy)/(n * Σx2 - Σx^2)
        r2  = (n * Σxy - Σx * Σy)^2/((n * Σx2 - Σx^2)*(n * Σy2 - Σy^2))
        if n > 2
            ar  = 1 - (1 - r2)*(n - 1)/(n - 2)
        else
            ar = NaN
        end
        return a, b, r2, ar
    end #end slope
        #Linear trapezoidal auc
    function linauc(t₁, t₂, c₁, c₂)::Number
        return (t₂-t₁)*(c₁+c₂)/2
    end
        #Linear trapezoidal aumc
    function linaumc(t₁, t₂, c₁, c₂)::Number
        return (t₂-t₁)*(t₁*c₁+t₂*c₂)/2
    end
        #Log trapezoidal auc
    function logauc(t₁, t₂, c₁, c₂)::Number
        return  (t₂-t₁)*(c₂-c₁)/log(c₂/c₁)
    end
        #Log trapezoidal aumc
    function logaumc(t₁, t₂, c₁, c₂)::Number
        return (t₂-t₁) * (t₂*c₂-t₁*c₁) / log(c₂/c₁) - (t₂-t₁)^2 * (c₂-c₁) / log(c₂/c₁)^2
    end
    #Intrapolation
        #linear prediction bx from ax, a1 < ax < a2
    function linpredict(a₁, a₂, ax, b₁, b₂)::Number
        return abs((ax - a₁) / (a₂ - a₁))*(b₂ - b₁) + b₁
    end

    function logtpredict(c₁, c₂, cx, t₁, t₂)::Number
        return log(cx/c₁)/log(c₂/c₁)*(t₂-t₁)+t₁
    end

    function logcpredict(t₁, t₂, tx, c₁, c₂)::Number
        return exp(log(c₁) + abs((tx-t₁)/(t₂-t₁))*(log(c₂) - log(c₁)))
    end
    function cpredict(t₁, t₂, tx, c₁, c₂, calcm)
        if calcm == :lint || c₂ >= c₁
            return linpredict(t₁, t₂, tx, c₁, c₂)
        else
            return logcpredict(t₁, t₂, tx, c₁, c₂)
        end
    end
    #---------------------------------------------------------------------------

    function aucpart(t₁, t₂, c₁, c₂, calcm, aftertmax)
        if calcm == :lint
            auc   =  linauc(t₁, t₂, c₁, c₂)    # (data[i - 1, conc] + data[i, conc]) * (data[i, time] - data[i - 1, time])/2
            aumc  = linaumc(t₁, t₂, c₁, c₂)    # (data[i - 1, conc]*data[i - 1, time] + data[i, conc]*data[i, time]) * (data[i, time] - data[i - 1, time])/2
        elseif calcm == :logt
            if aftertmax
                if c₁ > 0 && c₂ > 0
                    auc   =  logauc(t₁, t₂, c₁, c₂)
                    aumc  = logaumc(t₁, t₂, c₁, c₂)
                else
                    auc   =  linauc(t₁, t₂, c₁, c₂)
                    aumc  = linaumc(t₁, t₂, c₁, c₂)
                end
            else
                auc   =  linauc(t₁, t₂, c₁, c₂)
                aumc  = linaumc(t₁, t₂, c₁, c₂)
            end
        elseif calcm == :luld
            if c₁ > c₂ > 0
                auc   =  logauc(t₁, t₂, c₁, c₂)
                aumc  = logaumc(t₁, t₂, c₁, c₂)
            else
                auc   =  linauc(t₁, t₂, c₁, c₂)
                aumc  = linaumc(t₁, t₂, c₁, c₂)
            end
        elseif calcm == :luldt
            if aftertmax
                if c₁ > c₂ > 0
                    auc   =  logauc(t₁, t₂, c₁, c₂)
                    aumc  = logaumc(t₁, t₂, c₁, c₂)
                else
                    auc   =  linauc(t₁, t₂, c₁, c₂)
                    aumc  = linaumc(t₁, t₂, c₁, c₂)
                end
            else
                auc   =  linauc(t₁, t₂, c₁, c₂)
                aumc  = linaumc(t₁, t₂, c₁, c₂)
            end
        end
        return auc, aumc
    end

#-------------------------------------------------------------------------------
"""
    nca!(data::PKSubject; calcm = :lint, verbose = false, io::IO = stdout)

Pharmacokinetics non-compartment analysis for one PK subject.

calcm - calculation method;

- :lint  - Linear trapezoidal everywhere;
- :logt  - Log-trapezoidat rule after Tmax if c₁ > 0 and c₂ > 0, else Linear trapezoidal used;
- :luld  - Linear Up - Log Down everywhere if c₁ > c₂ > 0, else Linear trapezoidal used;
- :luldt - Linear Up - Log Down  after Tmax if c₁ > c₂ > 0, else Linear trapezoidal used;

intp - interpolation rule;
- :lint - linear interpolation;
- :logt - log interpolation;


verbose - print to out stream if "true";

- true;
- false.

io - out stream (stdout by default)
"""
function nca!(data::PKSubject; calcm = :lint, intp = :lint, verbose = false, io::IO = stdout)
        result   = Dict()
        #=
        result    = Dict(:Obsnum   => 0,   :Tmax   => 0,   :Tmaxn    => 0,   :Tlast   => 0,
        :Cmax    => 0,   :Cdose    => NaN, :Clast  => 0,   :Ctau     => NaN, :Ctaumin => NaN, :Cavg => NaN,
        :Swing   => NaN, :Swingtau => NaN, :Fluc   => NaN, :Fluctau  => NaN,
        :AUClast => 0,   :AUCall   => 0,   :AUCtau => 0,   :AUMClast => 0,   :AUMCall => 0,   :AUMCtau => 0,
        :MRTlast => 0,   :Kel      => NaN, :HL     => NaN,
        :Rsq     => NaN, :ARsq     => NaN, :Rsqn   => NaN,
        :AUCinf  => NaN, :AUCpct   => NaN)
        =#

        result[:Obsnum]  = length(data)
        result[:Cmax]    = maximum(data.obs)

        for i = result[:Obsnum]:-1:1
            if data.obs[i] > 0
                result[:Tlast]   = data.time[i]
                result[:Clast]   = data.obs[i]
                break
            end
        end
        result[:Cmax], result[:Tmax], result[:Tmaxn] = ctmax(data) #????!!!!
        #-----------------------------------------------------------------------
        if verbose
            println(io, "Non-compartmental Pharmacokinetic Analysis")
            printsortval(io, data.sort)
            println(io, "Settings:")
            println(io, "Method: $(calcm)")
        end
        #-----------------------------------------------------------------------
        #Areas
        aucpartl  = Array{Number, 1}(undef, result[:Obsnum] - 1)
        aumcpartl = Array{Number, 1}(undef, result[:Obsnum] - 1)
        #Calculate all UAC part based on data
        for i = 1:(result[:Obsnum] - 1)
            aucpartl[i], aumcpartl[i] = aucpart(data.time[i], data.time[i + 1], data.obs[i], data.obs[i + 1], calcm, i + 1 > result[:Tmaxn])
        end
        pmask   = Array{Bool, 1}(undef, result[:Obsnum] - 1)
        pmask  .= true

        if data.dosetime.time < data.time[1]
            #Dosetime < first time
        end
        #Find AUC part from dose time; exclude all before dosetime
        for i = 1:result[:Obsnum] - 1
            if  data.time[i] <= data.dosetime.time < data.time[i+1]
                if data.time[i] == data.dosetime.time
                    result[:Cdose] = data.obs[i]
                    break
                end
                result[:Cdose] = linpredict(data.time[i] , data.time[i+1], data.dosetime.time, data.obs[i], data.obs[i+1])
                aucpartl[i], aumcpartl[i] = aucpart(data.dosetime.time, data.time[i+1], result[:Cdose], data.obs[i+1], :lint, false) #? only :lint? always aftertmax = false?
                break
            else
                aucpartl[i]  = 0
                aumcpartl[i] = 0
                pmask[i]     = false
            end
        end
        #sum all AUC parts
        result[:AUCall]  = sum(aucpartl[pmask])
        result[:AUMCall] = sum(aumcpartl[pmask])
        #-----------------------------------------------------------------------
        if verbose
            aucpartlsum  = similar(aucpartl)
            aumcpartlsum = similar(aumcpartl)
            for i = 1:length(aucpartl)
                aucpartlsum[i]  = sum(aucpartl[1:i])
                aumcpartlsum[i] = sum(aumcpartl[1:i])
            end
            mx = hcat(data.time, data.obs, round.(vcat([0], aucpartl), digits = 3),  round.(vcat([0], aucpartlsum), digits = 3), round.(vcat([0], aumcpartl), digits = 3),  round.(vcat([0], aumcpartlsum), digits = 3))
            mx = vcat(permutedims(["Time", "Concentrtion", "AUC", "AUC (cumulate)", "AUMC", "AUMC (cumulate)"]), mx)
            printmatrix(io, mx)
            println(io, "")
        end
        #-----------------------------------------------------------------------
        #Exclude all AUC parts from observed concentation before 0 or less
        #Need elaborate!!!!
        for i = result[:Tmaxn]:result[:Obsnum] - 1
            if data.obs[i+1] <= 0 * data.obs[i+1] pmask[i:end] .= false; break end
        end
        result[:AUClast]  = sum(aucpartl[pmask])
        result[:AUMClast] = sum(aumcpartl[pmask])

        result[:MRTlast]       = result[:AUMClast] / result[:AUClast]
        #-----------------------------------------------------------------------
        #Elimination
        keldata                = KelData()
        if data.kelauto
            if result[:Obsnum] - result[:Tmaxn] >= 3
                logconc = log.(data.obs)
                cmask   = Array{Bool, 1}(undef, result[:Obsnum])
                for i  = result[:Tmaxn] + 1:result[:Obsnum] - 2
                    cmask                    .= false
                    cmask[i:result[:Obsnum]] .= true
                    sl = slope(data.time[cmask], logconc[cmask])
                    if sl[1] < 0 * sl[1]
                        push!(keldata, i, result[:Obsnum], sl[1], sl[2], sl[3], sl[4])
                    end
                end
            end
        else
            logconc = log.(data.obs)
            cmask   = Array{Bool, 1}(undef, result[:Obsnum])
            cmask  .= false
            cmask[data.kelrange.kelstart:data.kelrange.kelend] .= true
            cmask[data.kelrange.kelexcl] .= false
            sl = slope(data.time[cmask], logconc[cmask])
            push!(keldata, data.kelrange.kelstart, data.kelrange.kelend, sl[1], sl[2], sl[3], sl[4])
        end
        if  length(keldata) > 0
            result[:Rsq], result[:Rsqn] = findmax(keldata.ar)
            data.kelrange.kelstart   = keldata.s[result[:Rsqn]]
            result[:Kelstart]        = data.kelrange.kelstart
            data.kelrange.kelend     = keldata.e[result[:Rsqn]]
            result[:Kelend]          = data.kelrange.kelend
            data.keldata             = keldata
            result[:ARsq]            = keldata.ar[result[:Rsqn]]
            result[:Kel]             = abs(keldata.a[result[:Rsqn]])
            result[:LZ]              = keldata.a[result[:Rsqn]]
            result[:LZint]           = keldata.b[result[:Rsqn]]
            result[:Clast_pred]      = exp(result[:LZint] + result[:LZ]*result[:Tlast])
            result[:HL]              = LOG2 / result[:Kel]
            result[:AUCinf]          = result[:AUClast] + result[:Clast] / result[:Kel]
            result[:AUCpct]          = (result[:AUCinf] - result[:AUClast]) / result[:AUCinf] * 100.0
        end
        #-----------------------------------------------------------------------
        #Steady-state
        if data.dosetime.tau > 0
            ncas, ncae = ncarange(data, data.dosetime.time, data.dosetime.tau)
            tautime    = data.dosetime.time + data.dosetime.tau
            eaucpartl  = eaumcpartl = 0.0
            if tautime < data.time[end]
                #result[:Ctau] = linpredict(data.time[ncae] , data.time[ncae+1], tautime, data.obs[ncae], data.obs[ncae+1])
                if tautime > result[:Tmax] aftertmax = true else aftertmax = false end
                result[:Ctau] = cpredict(data.time[ncae], data.time[ncae+1], tautime, data.obs[ncae], data.obs[ncae+1], intp)
                aucpartl[ncae], aumcpartl[ncae] = aucpart(data.time[ncae], tautime, data.obs[ncae], result[:Ctau], calcm, aftertmax)
                if verbose
                    println(io, "Tau + dosetime is less then end time. Inrapolation used.")
                    println(io, "Interpolation between: $(data.time[ncae]) - $( data.time[ncae+1]), method: $(intp)")
                    println(io, "Ctau: $(result[:Ctau])")
                    println(io, "AUC final part: $(aucpartl[ncae])")
                end
                #remoove part after tau
                if ncae < result[:Obsnum] - 1 pmask[ncae+1:end] .= false end
            elseif tautime > data.time[end]
                #extrapolation
                result[:Ctau] = exp(result[:LZint] + result[:LZ] * tautime)
                eaucpartl, eaumcpartl = aucpart(data.time[ncae], tautime, data.obs[ncae], result[:Ctau], calcm, true)

                #!!!
            else
                result[:Ctau] = data.obs[end]
            end
            result[:Ctaumin]  = min(result[:Ctau], result[:Cdose], minimum(data.obs[ncas+1:ncae]))
            result[:AUCtau]   = sum(aucpartl[pmask])  + eaucpartl
            result[:AUMCtau]  = sum(aumcpartl[pmask]) + eaumcpartl
            result[:Cavg]     = result[:AUCtau]/data.dosetime.tau
            if result[:Ctaumin] != 0
                result[:Swing]    = (result[:Cmax] - result[:Ctaumin])/result[:Ctaumin]
            else
                result[:Swing]    = NaN
            end
            if result[:Ctau] != 0
                result[:Swingtau] = (result[:Cmax] - result[:Ctau])/result[:Ctau]
            else
                result[:Swingtau] = NaN
            end
            result[:Fluc]     = (result[:Cmax] - result[:Ctaumin])/result[:Cavg]
            result[:Fluctau]  = (result[:Cmax] - result[:Ctau])/result[:Cavg]
            #result[:AccInd]
        end
        #-----------------------------------------------------------------------
        return PKPDProfile(data, result; method = calcm)
    end
"""
    nca!(data::DataSet{PKSubject}; calcm = :lint, intp = :lint,
        verbose = false, io::IO = stdout)


Pharmacokinetics non-compartment analysis for PK subjects set.

calcm - calculation method;

- :lint  - Linear trapezoidal everywhere;
- :logt  - Log-trapezoidat rule after Tmax if c₁ > 0 and c₂ > 0, else Linear trapezoidal used;
- :luld  - Linear Up - Log Down everywhere if c₁ > c₂ > 0, else Linear trapezoidal used;
- :luldt - Linear Up - Log Down  after Tmax if c₁ > c₂ > 0, else Linear trapezoidal used;

intp - interpolation rule;
- :lint - linear interpolation;
- :logt - log interpolation;

verbose - print to out stream if "true";

- true;
- false.

io - out stream (stdout by default)
"""
function nca!(data::DataSet{PKSubject}; calcm = :lint, intp = :lint, verbose = false, io::IO = stdout)
        results = Array{PKPDProfile, 1}(undef, 0)
        for i = 1:length(data.data)
            push!(results, nca!(data.data[i]; calcm = calcm, intp = intp, verbose = verbose, io = io))

        end
        return DataSet(results)
end
    #---------------------------------------------------------------------------
"""
    nca!(data::PDSubject; verbose = false, io::IO = stdout)::PKPDProfile{PDSubject}

Pharmacodynamics non-compartment analysis for one PD subject.
"""
function nca!(data::PDSubject; verbose = false, io::IO = stdout)::PKPDProfile{PDSubject}
        result = Dict(:TH => NaN, :AUCABL => 0, :AUCBBL => 0, :AUCBLNET => 0,  :TABL => 0, :TBBL => 0)
        result[:Obsnum] = length(data)
        result[:TH] = data.th
        result[:BL] = data.bl
        result[:RMAX] = maximum(data.obs)
        result[:RMIN] = maximum(data.obs)
        for i = 2:obsnum(data) #BASELINE
            if data.obs[i - 1] <= result[:BL] && data.obs[i] <= result[:BL]
                result[:TBBL]   += data.time[i,] - data.time[i - 1]
                result[:AUCBBL] += linauc(data.time[i - 1], data.time[i], result[:BL] - data.obs[i - 1], result[:BL] - data.obs[i])
            elseif data.obs[i - 1] <= result[:BL] &&  data.obs[i] > result[:BL]
                tx      = linpredict(data.obs[i - 1], data.obs[i], result[:BL], data.time[i - 1], data.time[i])
                result[:TBBL]   += tx - data.time[i - 1]
                result[:TABL]   += data.time[i] - tx
                result[:AUCBBL] += (tx - data.time[i - 1]) * (result[:BL] - data.obs[i - 1]) / 2
                result[:AUCABL] += (data.time[i] - tx) * (data.obs[i] - result[:BL]) / 2 #Ok
            elseif data.obs[i - 1] > result[:BL] &&  data.obs[i] <= result[:BL]
                tx      = linpredict(data.obs[i - 1], data.obs[i], result[:BL], data.time[i - 1], data.time[i])
                result[:TBBL]   += data.time[i] - tx
                result[:TABL]   += tx - data.time[i - 1]
                result[:AUCBBL] += (data.time[i] - tx) * (result[:BL] - data.obs[i]) / 2
                result[:AUCABL] += (tx - data.time[i - 1]) * (data.obs[i - 1] - result[:BL]) / 2
            elseif data.obs[i - 1] > result[:BL] &&  data.obs[i] > result[:BL]
                result[:TABL]   += data.time[i] - data.time[i - 1]
                result[:AUCABL]     += linauc(data.time[i - 1], data.time[i], data.obs[i - 1] - result[:BL], data.obs[i] - result[:BL])
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
            for i = 2:length(data)
                if data.obs[i - 1] <= result[:TH] && data.obs[i] <= result[:TH]
                    result[:TBTH]   += data.time[i] - data.time[i - 1]
                    result[:AUCBTH] += linauc(data.time[i - 1], data.time[i], result[:TH] - data.obs[i - 1], result[:TH] - data.obs[i])
                elseif data.obs[i - 1] <= result[:TH] &&  data.obs[i] > result[:TH]
                    tx      = linpredict(data.obs[i - 1], data.obs[i], result[:TH], data.time[i - 1], data.time[i])
                    result[:TBTH]   += tx - data.time[i - 1]
                    result[:TATH]   += data.time[i] - tx
                    result[:AUCBTH] += (tx - data.time[i - 1]) * (result[:TH] - data.obs[i - 1]) / 2
                    result[:AUCATH] += (data.time[i] - tx) * (data.obs[i] - result[:TH]) / 2 #Ok
                elseif data.obs[i - 1] > result[:TH] &&  data.obs[i] <= result[:TH]
                    tx      = linpredict(data.obs[i - 1], data.obs[i], result[:TH], data.time[i - 1], data.time[i])
                    result[:TBTH]   += data.time[i] - tx
                    result[:TATH]   += tx - data.time[i - 1]
                    result[:AUCBTH] += (data.time[i] - tx) * (result[:TH] - data.obs[i]) / 2
                    result[:AUCATH] += (tx - data.time[i - 1]) * (data.obs[i - 1] - result[:TH]) / 2
                elseif data.obs[i - 1] > result[:TH] &&  data.obs[i] > result[:TH]
                    result[:TATH]   += data.time[i] - data.time[i - 1]
                    result[:AUCATH] += linauc(data.time[i - 1], data.time[i], data.obs[i - 1] - result[:TH], data.obs[i] - result[:TH])
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
"""
    nca!(data::PDSubject; verbose = false, io::IO = stdout)::PKPDProfile{PDSubject}

Pharmacodynamics non-compartment analysis for PD subjects set.
"""
function nca!(data::DataSet{PDSubject}; verbose = false, io::IO = stdout)
        results = Array{PKPDProfile, 1}(undef, 0)
        for i = 1:length(data.data)
            push!(results, nca!(data.data[i]))
        end
        return DataSet(results)
    end
    #---------------------------------------------------------------------------
"""
    pkimport(data::DataFrame, sort::Array, rule::LimitRule; conc::Symbol, time::Symbol)

Pharmacokinetics data import from DataFrame.

- data - sourece DataFrame;
- sort - sorting columns;
- rule - applied LimitRule.

- conc - concentration column;
- time - time column.
"""
function pkimport(data::DataFrame, sort::Array, rule::LimitRule; conc::Symbol, time::Symbol)
        sortlist = unique(data[:, sort])
        results  = Array{PKSubject, 1}(undef, 0)
        for i = 1:size(sortlist, 1) #For each line in sortlist
            datai = DataFrame(Time = eltype(data[!,time])[], Conc = eltype(data[!,conc])[])
            for c = 1:size(data, 1) #For each line in data
                if data[c, sort] == sortlist[i,:]
                    push!(datai, [data[c, time], data[c, conc]])
                end
            end
            if !(eltype(datai.Time) <: Real)
                datai.Time[.!isa.(datai.Time, Real)] .= NaN
                datai.Time = float.(datai.Time)
            end
            if !(eltype(datai.Conc) <: Real)
                datai.Conc[.!isa.(datai.Conc, Real)] .= NaN
                datai.Conc = float.(datai.Conc)
            end
            if rule.lloq !== NaN
                ncarule!(datai, :Conc, :Time, rule)
            else
                sort!(datai, :Time)
            end
            push!(results, PKSubject(datai[!, :Time], datai[!, :Conc], sort = Dict(sort .=> collect(sortlist[i,:]))))
        end
        return DataSet(results)
    end
"""
    pkimport(data::DataFrame, sort::Array; time::Symbol, conc::Symbol)

Pharmacokinetics data import from DataFrame.

- data - sourece DataFrame;
- sort - sorting columns.

- conc - concentration column;
- time - time column.
"""
function pkimport(data::DataFrame, sort::Array; time::Symbol, conc::Symbol)
        rule = LimitRule()
        return pkimport(data, sort, rule; conc = conc, time = time)
end
"""
    pkimport(data::DataFrame, sort::Array; time::Symbol, conc::Symbol)

Pharmacokinetics data import from DataFrame for one subject.

- data - sourece DataFrame;
- rule - applied LimitRule.

- conc - concentration column;
- time - time column.
"""
function pkimport(data::DataFrame, rule::LimitRule; time::Symbol, conc::Symbol)
        #rule = LimitRule()
        datai = ncarule!(copy(data[!,[time, conc]]), conc, time, rule)
        return DataSet([PKSubject(datai[!, time], datai[!, conc])])
end
"""
    pkimport(data::DataFrame, sort::Array; time::Symbol, conc::Symbol)

Pharmacokinetics data import from DataFrame for one subject.

- data - sourece DataFrame;

- conc - concentration column;
- time - time column.
"""
function pkimport(data::DataFrame; time::Symbol, conc::Symbol)
        datai = sort(data[!,[time, conc]], time)
        return DataSet([PKSubject(datai[!, time], datai[!, conc])])
end
    #---------------------------------------------------------------------------
"""
    pdimport(data::DataFrame, sort::Array; resp::Symbol, time::Symbol,
        bl::Real = 0, th::Real = NaN)

Pharmacodynamics data import from DataFrame.

- data - sourece DataFrame;
- sort - sorting columns;

- resp - responce column;
- time - time column
- bl - baseline;
- th - treashold.
"""
function pdimport(data::DataFrame, sort::Array; resp::Symbol, time::Symbol, bl::Real = 0, th::Real = NaN)
        sortlist = unique(data[:, sort])
        results  = Array{PDSubject, 1}(undef, 0)
        for i = 1:size(sortlist, 1) #For each line in sortlist
            datai = DataFrame(Time = Real[], Resp = Real[])
            for c = 1:size(data, 1) #For each line in data
                if data[c, sort] == sortlist[i,:]
                    push!(datai, [data[c, time], data[c, resp]])
                end
            end
            sort!(datai, :Time)
            push!(results, PDSubject(datai[!, :Time], datai[!, :Resp], bl, th, sort = Dict(sort .=> collect(sortlist[i,:]))))
        end
        return DataSet(results)
end
    function pdimport(data::DataFrame; resp::Symbol, time::Symbol, bl = 0, th = NaN)
        datai = sort(data[!,[time, resp]], time)
        return DataSet([PDSubject(datai[!, time], datai[!, resp], bl, th)])
    end
    #---------------------------------------------------------------------------
    function applyncarule(data::PKSubject, rule::LimitRule)
    end
    function applyncarule(data::DataSet{PKSubject}, rule::LimitRule)
    end
    function applyncarule(data::PKPDProfile{PKSubject}, rule::LimitRule)
    end
    function applyncarule(data::DataSet{PKPDProfile{PKSubject}}, rule::LimitRule)
    end
    function applyncarule!(data::PKSubject, rule::LimitRule)
    end
    function applyncarule!(data::DataSet{PKSubject}, rule::LimitRule)
    end
    function applyncarule!(data::PKPDProfile{PKSubject}, rule::LimitRule)
    end
    function applyncarule!(data::DataSet{PKPDProfile{PKSubject}}, rule::LimitRule)
    end
    #---------------------------------------------------------------------------
"""
    setelimrange!(data::PKSubject, range::ElimRange)

Set range for elimination parameters calculation for subject.

data - PK subject;
range - ElimRange object.
"""
function setelimrange!(data::PKSubject, range::ElimRange)
        if range.kelend > length(data) throw(ArgumentError("Kel endpoint out of range")) end
        data.kelrange = range
end

function setelimrange!(data::DataSet{PKSubject}, range::ElimRange)
        for i = 1:length(data)
            setelimrange!(data[i], range)
        end
        data
end
function setelimrange!(data::DataSet{PKSubject}, range::ElimRange, subj::Array{Int,1})
        for i in subj
            setelimrange!(data[i], range)
        end
        data
end
function setelimrange!(data::DataSet{PKSubject}, range::ElimRange, subj::Int)
        setelimrange!(data[subj], range)
        data
end
#-------------------------------------------------------------------------------
"""
    setkelauto!(data::PKSubject, kelauto::Bool)

Set range for elimination parameters calculation for subject.

data     - PK subject;
kelauto  - value.
"""
function setkelauto!(data::PKSubject, kelauto::Bool)
        data.kelauto = kelauto
end

function setkelauto!(data::DataSet{PKSubject}, kelauto::Bool)
        for i = 1:length(data)
            setkelauto!(data[i], kelauto)
        end
        data
end
function setkelauto!(data::DataSet{PKSubject}, kelauto::Bool, subj::Array{Int,1})
        for i in subj
            setkelauto!(data[i], kelauto)
        end
        data
end
function setkelauto!(data::DataSet{PKSubject}, kelauto::Bool, subj::Int)
        setkelauto!(data[subj], kelauto)
        data
end
#-------------------------------------------------------------------------------
"""
    applyelimrange!(data::PKPDProfile{PKSubject}, range::ElimRange)

Set range for elimination parameters calculation for profile's subject and recalculate NCA parameters.

data - PK subject;
range - ElimRange object.
"""
function applyelimrange!(data::PKPDProfile{PKSubject}, range::ElimRange)
        setelimrange!(data.subject, range)
        data.result = nca!(data.subject, calcm = data.method).result
        data
end
    function applyelimrange!(data::DataSet{PKPDProfile{PKSubject}}, range::ElimRange)
        for i = 1:length(data)
            applyelimrange!(data[i], range)
        end
        data
    end
    function applyelimrange!(data::DataSet{PKPDProfile{PKSubject}}, range::ElimRange, subj::Array{Int,1})
        for i = 1:length(data)
            if i ∈ subj applyelimrange!(data[i], range) end
        end
        data
    end
    function applyelimrange!(data::DataSet{PKPDProfile{PKSubject}}, range::ElimRange, subj::Int)
        applyelimrange!(data[subj], range)
        data
    end
    function applyelimrange!(data::DataSet{PKPDProfile{PKSubject}}, range::ElimRange, sort::Dict)
        for i = 1:length(data)
            if sort ∈ data[i].subject.sort
                applyelimrange!(data[i], range)
            end
        end
        data
    end
    #---------------------------------------------------------------------------
"""
    setdosetime!(data::PKSubject, dosetime::DoseTime)

Set dosing time and tau for subject.

data - PK subject;
dosetime - DoseTime object.
"""
function setdosetime!(data::PKSubject, dosetime::DoseTime)
    data.dosetime = dosetime
    data
end
function setdosetime!(data::DataSet{PKSubject}, dosetime::DoseTime)
    for i = 1:length(data)
        setdosetime!(data[i], dosetime)
    end
    data
end
function setdosetime!(data::DataSet{PKSubject}, dosetime::DoseTime, sort::Dict)
    for i = 1:length(data)
        if sort ∈ data[i].sort setdosetime!(data[i], dosetime) end
    end
    data
end
"""
    setth!(data::PDSubject, th::Real)

Set treashold for subject.

data - PK subject;
th - treashold.
"""
function setth!(data::PDSubject, th::Real)
        data.th = th
        data
end
function setth!(data::DataSet{PDSubject}, th::Real)
    for i = 1:length(data)
        setth!(data[i], th)
    end
    data
end
function setth!(data::DataSet{PDSubject}, th::Real, sort::Dict)
    for i = 1:length(data)
        if sort ∈ data[i].sort setth!(data[i], th) end
    end
    data
end

"""
    setbl!(data::PDSubject, bl::Real)

Set baseline for subject.

data - PK subject;
bl - baseline.
"""
function setbl!(data::PDSubject, bl::Real)
    data.bl = bl
    data
end
function setbl!(data::DataSet{PDSubject}, bl::Real)
    for i = 1:length(data)
        setbl!(data[i], bl)
    end
    data
end
function setbl!(data::DataSet{PDSubject}, bl::Real, sort::Dict)
    for i = 1:length(data)
        if sort ∈ data[i].sort setbl!(data[i], bl) end
    end
    data
end
#-------------------------------------------------------------------------------
"""
    keldata(data::PKPDProfile{PKSubject})

Return elimination data for PK subject in profile.
"""
function keldata(data::PKPDProfile{PKSubject})
    return data.subject.keldata
end
function keldata(data::DataSet{PKPDProfile{PKSubject}}, i::Int)
        return data[i].subject.keldata
end

"""
    DataFrames.DataFrame(data::DataSet{PKPDProfile}; unst = false, us = false)

Make datafrafe from PK/PD DataSet.

unst | us - unstack data;
"""
function DataFrames.DataFrame(data::DataSet{T}; unst = false, us = false) where T
        d = DataFrame(id = Int[], sortvar = Symbol[], sortval = Any[])
        for i = 1:length(data)
            if length(data[i].sort) > 0
                for s in data[i].sort
                    push!(d, [i, s[1], s[2]])
                end
            end
        end
        d   = unstack(d, :sortvar, :sortval)[!,2:end]
        df  = DataFrame()
        dfn = names(d)
        for i = 1:size(d,2)
            df[!,dfn[i]] = Array{eltype(d[!, dfn[i]]), 1}(undef, 0)
        end
        df = hcat(df, DataFrame(param = Any[], value = Real[]))
        for i = 1:size(d,1)
            a = Array{Any,1}(undef, size(d,2))
            copyto!(a, collect(d[i,:]))
            for p in data[i].result
                r = append!(copy(a), collect(p))
            push!(df, r)
            end
        end
        if unst || us
            return unstack(df, names(df)[end-1], names(df)[end])
        else
            return df
        end
end
