#Pharmacokinetics
#Makoid C, Vuchetich J, Banakar V. 1996-1999. Basic Pharmacokinetics.
#!!!

isnanormissing(x::Number) = isnan(x)
isnanormissing(x::Missing) = true


struct KelData
    s::Array{Int, 1}
    e::Array{Int, 1}
    a::Array{Float64, 1}
    b::Array{Float64, 1}
    r::Array{Float64, 1}
    ar::Array{Float64, 1}
    function KelData(s, e, a, b, r, ar)::KelData
        new(s, e, a, b, r, ar)::KelData
    end
    function KelData()::KelData
        new(Int[], Int[], Float64[], Float64[], Float64[], Float64[])::KelData
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


#Plasma PK subject
mutable struct PKSubject <: AbstractSubject
    time::Vector{Float64}
    obs::Vector{Float64}
    kelauto::Bool
    kelrange::ElimRange
    dosetime::DoseTime
    keldata::KelData
    sort::Dict
    function PKSubject(time::Vector, conc::Vector, kelauto::Bool, kelrange::ElimRange, dosetime::DoseTime, keldata::KelData; sort = Dict())
        new(time, conc, kelauto, kelrange, dosetime, keldata, sort)::PKSubject
    end
    function PKSubject(time::Vector, conc::Vector, kelauto::Bool, kelrange, dosetime; sort = Dict())
        PKSubject(time, conc, kelauto, kelrange, dosetime, KelData(); sort = sort)
    end
    function PKSubject(time::Vector, conc::Vector; sort = Dict())
        PKSubject(time, conc, true, ElimRange(), DoseTime(NaN, 0 * time[1]), KelData(); sort = sort)
    end
    function PKSubject(data; time::Symbol, conc::Symbol, timesort::Bool = false, sort = Dict())
        if timesort sort!(data, time) end
        #Check double time
        PKSubject(copy(data[!,time]), copy(data[!,conc]), true, ElimRange(), DoseTime(), KelData(); sort = sort)
    end
end

#Urine PK subject
mutable struct UPKSubject <: AbstractSubject
    stime::Vector
    etime::Vector
    conc::Vector
    vol::Vector
    time::Vector
    obs::Vector
    kelauto::Bool
    kelrange::ElimRange
    dosetime::DoseTime
    keldata::KelData
    sort::Dict
    function UPKSubject(stime::Vector, etime::Vector, conc::Vector, vol::Vector; sort = Dict())
        new(stime, etime, conc, vol, (stime .+ etime) / 2, (conc .* vol) ./  (etime .- stime), true, ElimRange(), DoseTime(NaN, 0 * stime[1]), KelData(), sort)::UPKSubject
    end
end

#PD subject
mutable struct PDSubject <: AbstractSubject
    time::Vector
    obs::Vector
    bl::Float64
    th::Float64
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
"""
    LimitRule(;lloq::Real = NaN, btmax::Real = NaN, atmax::Real = NaN, nan::Real = NaN, rm::Bool = false)


    LimitRule(lloq::Real, btmax::Real, atmax::Real, nan::Real, rm::Bool)

Rule for PK subject.

STEP 1 (NaN step): replace all NaN values with nan reyword value (if nan !== NaN);
STEP 2 (LLOQ step): replace values below lloq with btmax value if this value befor Tmax or with atmax if this value after Tmax (if lloq !== NaN);
STEP 3 (remove NaN): rm == true, then remove all NaN values;
"""
struct LimitRule{T <: Real}
    lloq::T
    btmax::Float64
    atmax::Float64
    nan::Float64
    rm::Bool
    function LimitRule(lloq::T, btmax, atmax, nan, rm::Bool) where T <: Real
        new{T}(lloq, btmax, atmax, nan, rm)::LimitRule
    end
    function LimitRule(lloq, btmax, atmax, nan)
        LimitRule(lloq, btmax, atmax, nan, false)
    end
    function LimitRule(lloq, btmax, atmax)
        LimitRule(lloq, btmax, atmax, NaN, true)
    end
    function LimitRule(;lloq = NaN, btmax = NaN, atmax = NaN, nan = NaN, rm::Bool = false)
        LimitRule(lloq, btmax, atmax, nan, rm)
    end
end


function Base.getindex(a::PKPDProfile{T}, s::Symbol) where T <: AbstractSubject
    return a.result[s]
end

#=
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
#=
function obsnum(data::T) where T <:  AbstractSubject
    return length(data.time)
end
function obsnum(keldata::KelData)
    return length(keldata.a)
end
=#
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
    println(io, "  Pharmacokinetic subject")
    println(io, "Observations: $(length(obj)); ", obj.dosetime)
    println(io,  obj.kelrange)
    println(io, "Time   Concentration")
    for i = 1:length(obj)
        println(io, obj.time[i], " => ", obj.obs[i])
    end
end
function Base.show(io::IO, obj::PKPDProfile)
    println(io, "PK/PD Profile")
    println(io, "Type: $(obj.subject)")
    println(io, "Method: $(obj.method)")
    print(io, "Results:")
    for k in keys(obj.result)
        print(io, " $(k)")
    end
    println(io, ";")
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
        print(io, "Exclusions: $(obj.kelexcl[1])")
        if length(obj.kelexcl) > 1 for i = 2:length(obj.kelexcl) print(io, ", $(obj.kelexcl[i])") end end
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
    printmatrix(io, m)
end
#-------------------------------------------------------------------------------
function notnanormissing(x)
    if x === NaN || x === missing return false else true end
end
    """
        cmax
        tmax
        tmaxn
    """
    function ctmax(time::Vector, obs::Vector)
        if length(time) != length(obs) throw(ArgumentError("Unequal vector length!")) end
        if length(obs[notnanormissing.(obs)]) > 0 cmax  = maximum(obs[notnanormissing.(obs)]) else cmax = NaN end
        tmax  = NaN
        tmaxn = 0
        for i = 1:length(time)
            if obs[i] !== missing
                if obs[i] == cmax
                    tmax  = time[i]
                    tmaxn = i
                    break
                end
            end
        end
        return cmax, tmax, tmaxn
    end
    #=
    function ctmax(data, conc::Symbol, time::Symbol)
        return ctmax(data[!, time], data[!, conc])
    end
    =#
    function ctmax(data::PKSubject)
        return ctmax(data.time, data.obs)
    end
#=
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
=#
function dropbeforedosetime(time::Vector, obs::Vector, dosetime::DoseTime)
    s = 0
    for i = 1:length(time)
        @inbounds if time[i] < dosetime.time || isnan(obs[i]) s = i else break end
    end
    return time[s+1:end], obs[s+1:end]
end
    """
    Range for AUC
        return start point number, end point number
    """
    function ncarange(time::Vector, dosetime, tau)
        tautime = dosetime + tau
        s     = 1
        e     = length(time)
        for i = 1:length(time) - 1
            if dosetime > time[i] && dosetime <= time[i+1] s = i+1;
                break
            end
        end
        if tautime >= time[end]
            e = length(time)
        else
            for i = s:length(time) - 1
                if tautime >= time[i] && tautime < time[i+1] e = i; break end
            end
        end
        return s, e
    end

    #---------------------------------------------------------------------------
    function slope(x::AbstractVector{T1}, y::AbstractVector{T2}) where T1 where T2
        if length(x) != length(y) throw(ArgumentError("Unequal vector length!")) end
        n   = length(x)
        if n < 2 throw(ArgumentError("n < 2!")) end
        T = promote_type(T1, T2)
        Σxy = Σx = Σy = Σx2 = Σy2 = zero(T)
        @inbounds for i = 1:n
            Σxy += x[i] * y[i]
            Σx  += x[i]
            Σy  += y[i]
            Σx2 += x[i]^2
            Σy2 += y[i]^2
        end
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
    function linauc(t₁, t₂, c₁, c₂)
        return (t₂-t₁)*(c₁+c₂)/2
    end
        #Linear trapezoidal aumc
    function linaumc(t₁, t₂, c₁, c₂)
        return (t₂-t₁)*(t₁*c₁+t₂*c₂)/2
    end
        #Log trapezoidal auc
    function logauc(t₁, t₂, c₁, c₂)
        return  (t₂-t₁)*(c₂-c₁)/log(c₂/c₁)
    end
        #Log trapezoidal aumc
    function logaumc(t₁, t₂, c₁, c₂)
        return (t₂-t₁) * (t₂*c₂-t₁*c₁) / log(c₂/c₁) - (t₂-t₁)^2 * (c₂-c₁) / log(c₂/c₁)^2
    end
    #Intrapolation
        #linear prediction bx from ax, a1 < ax < a2
    function linpredict(a₁, a₂, ax, b₁, b₂)
        return (ax - a₁) / (a₂ - a₁)*(b₂ - b₁) + b₁
    end

    function logtpredict(c₁, c₂, cx, t₁, t₂)
        return log(cx/c₁)/log(c₂/c₁)*(t₂-t₁)+t₁
    end

    function logcpredict(t₁, t₂, tx, c₁, c₂)
        return exp(log(c₁) + (tx-t₁)/(t₂-t₁)*(log(c₂) - log(c₁)))
    end
    function cpredict(t₁, t₂, tx, c₁, c₂, calcm)
        if calcm == :lint || c₂ >= c₁
            return linpredict(t₁, t₂, tx, c₁, c₂)
        else
            return logcpredict(t₁, t₂, tx, c₁, c₂)
        end
    end
#---------------------------------------------------------------------------
#=
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
        elseif calcm == :luld
            if c₁ > c₂ > 0
                auc   =  logauc(t₁, t₂, c₁, c₂)
                aumc  = logaumc(t₁, t₂, c₁, c₂)
            else
                auc   =  linauc(t₁, t₂, c₁, c₂)
                aumc  = linaumc(t₁, t₂, c₁, c₂)
            end
        end
        return auc, aumc
    end
=#
function aucpart(t₁, t₂, c₁::T1, c₂::T2, calcm, aftertmax, dosetime) where T1 where T2

    if calcm == :lint
        auc   =  linauc(t₁, t₂, c₁, c₂)
        aumc  = linaumc(t₁ - dosetime, t₂ - dosetime, c₁, c₂)
    elseif calcm == :logt && aftertmax && c₁ > zero(T1) && c₂ > zero(T2)
        auc   =  logauc(t₁, t₂, c₁, c₂)
        aumc  = logaumc(t₁ - dosetime, t₂ - dosetime, c₁, c₂)
    elseif calcm == :luldt && aftertmax && c₁ > c₂ > zero(T2)
        auc   =  logauc(t₁, t₂, c₁, c₂)
        aumc  = logaumc(t₁ - dosetime, t₂ - dosetime, c₁, c₂)
    elseif calcm == :luld &&  c₁ > c₂ > zero(T2)
        auc   =  logauc(t₁, t₂, c₁, c₂)
        aumc  = logaumc(t₁ - dosetime, t₂ - dosetime, c₁, c₂)
    else
        auc   =  linauc(t₁, t₂, c₁, c₂)
        aumc  = linaumc(t₁ - dosetime, t₂ - dosetime, c₁, c₂)
    end
    return auc, aumc
end
#-------------------------------------------------------------------------------
"""
    nca!(data::PKSubject; adm = :ev, calcm = :lint, intp = :lint,
        verbose = false, warn = true, io::IO = stdout)

Pharmacokinetics non-compartment analysis for one PK subject.

`adm` - administration:

* `:ev` - extravascular;
* `:iv` - intra vascular bolus;

`calcm` - calculation method;

- :lint  - Linear trapezoidal everywhere;
- :logt  - Log-trapezoidat rule after Tmax if c₁ > 0 and c₂ > 0, else Linear trapezoidal used;
- :luld  - Linear Up - Log Down everywhere if c₁ > c₂ > 0, else Linear trapezoidal used;
- :luldt - Linear Up - Log Down after Tmax if c₁ > c₂ > 0, else Linear trapezoidal used;

`intp` - interpolation rule;
- :lint - linear interpolation;
- :logt - log interpolation;


`verbose` - print to out stream if "true";

- true;
- false.

`warn` - warnings enabled;

`io` - out stream (stdout by default).
"""
function nca!(data::PKSubject; adm = :ev, calcm = :lint, intp = :lint, verbose = false, warn = true, io::IO = stdout)
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
        if warn
            if !allunique(data.time)
                @warn "Not all time values is unique!"
            end
        end
        #if  calcm == :logt / :luld / :luldt calculation method used - log interpolation using
        if calcm == :logt || calcm == :luld  || calcm == :luldt intp = :logt end

        #time             = data.time
        #obs              = data.obs
        time, obs        = dropbeforedosetime(data.time, data.obs, data.dosetime)
        #dropnan
        #nanlist  = findall(isnan, obs)
        #if length(nanlist) > 0
        #    deleteat!(time, nanlist)
        #    deleteat!(obs, nanlist)
        #end
        result[:Obsnum]  = length(obs)
        rmn = length(data.obs) - length(obs)
        if result[:Obsnum] < 2
            return PKPDProfile(data, result; method = calcm)
        end

        result[:Cmax]    = maximum(obs)
        doseaucpart      = 0.0
        doseaumcpart     = 0.0
        daucpartl        = 0.0
        daumcpartl       = 0.0

        for i = result[:Obsnum]:-1:1
            if obs[i] > 0
                result[:Tlast]   = time[i]
                result[:Clast]   = obs[i]
                break
            end
        end
        result[:Cmax], result[:Tmax], result[:Tmaxn] = ctmax(time, obs)
        #result[:Tmaxn] += rmn
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
            aucpartl[i], aumcpartl[i] = aucpart(time[i], time[i + 1], obs[i], obs[i + 1], calcm, i + 1 > result[:Tmaxn], data.dosetime.time)
        end
        #pmask   = Array{Bool, 1}(undef, result[:Obsnum] - 1)
        #pmask  .= true
        pmask   = trues(result[:Obsnum] - 1)

        ncas    = nothing
        ncae    = nothing

        # If TAU set, calculates start and edt timepoints for AUCtau
        if  data.dosetime.tau > 0
            ncas, ncae       = ncarange(time, data.dosetime.time, data.dosetime.tau)
            result[:Ctaumin] = minimum(view(obs, ncas:ncae))
        end

        #Elimination
        keldata                = KelData()
        if data.kelauto
            if result[:Obsnum] - result[:Tmaxn] >= 3
                logconc    = log.(obs)
                cmask      = trues(result[:Obsnum])
                kelexclrmn = data.kelrange.kelexcl .- rmn
                filter!(x-> 0 < x <= result[:Obsnum], kelexclrmn)
                view(cmask, kelexclrmn) .= false
                kellast = findlast(x -> x, cmask)
                if kellast - result[:Tmaxn] >= 3
                    if adm == :iv tmnc = 0 else tmnc = 1 end
                    for i  = result[:Tmaxn] + tmnc:kellast - 2
                        cmask                    .= false
                        #cmask[i:result[:Obsnum]] .= true
                        view(cmask, i:kellast) .= true
                        if length(data.kelrange.kelexcl) > 0
                            view(cmask, kelexclrmn) .= false
                        end
                        if sum(cmask) >= 3
                            sl = slope(time[cmask], logconc[cmask])
                            if sl[1] < 0
                                push!(keldata, i + rmn, kellast + rmn, sl[1], sl[2], sl[3], sl[4])
                            end
                        end
                    end
                end
            end
        else
            logconc = log.(obs)
            #cmask   = Array{Bool, 1}(undef, result[:Obsnum])
            #cmask  .= false
            cmask   = falses(result[:Obsnum])
            cmask[data.kelrange.kelstart-rmn:data.kelrange.kelend-rmn] .= true
            cmask[data.kelrange.kelexcl .- rmn] .= false
            sl = slope(view(time, cmask), view(logconc, cmask))
            push!(keldata, data.kelrange.kelstart, data.kelrange.kelend, sl[1], sl[2], sl[3], sl[4])
        end

        #Calcalation partial areas (doseaucpart, doseaumcpart)
        #Dosetime < first time
        if  data.dosetime.time == time[1]
            result[:Cdose] = obs[1]
        elseif data.dosetime.time < time[1]
            if adm == :iv
                if obs[1] > obs[2] > 0
                    result[:Cdose] = logcpredict(time[1], time[2], data.dosetime.time, obs[1], obs[2])
                else
                    result[:Cdose] = obs[findfirst(x->x>0, obs)]
                end
            else
                if  data.dosetime.tau > 0
                    result[:Cdose] = result[:Ctaumin]
                else
                    result[:Cdose] = 0
                end
            end
            doseaucpart, doseaumcpart  = aucpart(data.dosetime.time, time[1], result[:Cdose], obs[1], :lint, false, data.dosetime.time) #always :lint?
        else
            error("Some concentration before dose time after filtering!!!")
        end

        #sum all full AUC parts and part before dose
        result[:AUCall]  = sum(view(aucpartl, pmask))  + doseaucpart
        result[:AUMCall] = sum(view(aumcpartl, pmask)) + doseaumcpart
        #-----------------------------------------------------------------------
        #-----------------------------------------------------------------------
        #Find last mesurable concentration (>0) from start, all after excluded
        for i = result[:Tmaxn]:result[:Obsnum] - 1
            if obs[i+1] <= 0 pmask[i:end] .= false; break end
        end
        #=
        #Find last mesurable concentration (>0) from end, all after excluded
        for i = result[:Obsnum]:-1:1
            if obs[i] <= 0  pmask[i] = false else break end
        end
        =#
        result[:AUClast]       = sum(view(aucpartl, pmask)) + doseaucpart
        result[:AUMClast]      = sum(view(aumcpartl, pmask)) + doseaumcpart
        result[:MRTlast]       = result[:AUMClast] / result[:AUClast]


        result[:Cllast]           = data.dosetime.dose / result[:AUClast]
        #-----------------------------------------------------------------------
        #-----------------------------------------------------------------------

        if  length(keldata) > 0
            result[:ARsq], result[:Rsqn] = findmax(keldata.ar)
            data.kelrange.kelstart   = keldata.s[result[:Rsqn]]
            result[:Kelstart]        = data.kelrange.kelstart
            data.kelrange.kelend     = keldata.e[result[:Rsqn]]
            result[:Kelend]          = data.kelrange.kelend
            data.keldata             = keldata
            result[:Rsq]             = keldata.r[result[:Rsqn]]
            result[:Kel]             = abs(keldata.a[result[:Rsqn]])
            result[:LZ]              = keldata.a[result[:Rsqn]]
            result[:LZint]           = keldata.b[result[:Rsqn]]
            result[:Clast_pred]      = exp(result[:LZint] + result[:LZ]*result[:Tlast])
            result[:HL]              = LOG2 / result[:Kel]
            result[:AUCinf]          = result[:AUClast] + result[:Clast] / result[:Kel]
            result[:AUCinf_pred]     = result[:AUClast] + result[:Clast_pred] / result[:Kel]
            result[:AUCpct]          = (result[:AUCinf] - result[:AUClast]) / result[:AUCinf] * 100.0
            result[:AUMCinf]         = result[:AUMClast] + result[:Tlast] * result[:Clast] / result[:Kel] + result[:Clast] / result[:Kel] ^ 2
            result[:MRTinf]          = result[:AUMCinf] / result[:AUCinf]
            result[:Vzlast]           = data.dosetime.dose / result[:AUClast] / result[:Kel]
            result[:Vzinf]           = data.dosetime.dose / result[:AUCinf] / result[:Kel]
            result[:Clinf]           = data.dosetime.dose / result[:AUCinf]
            result[:Vssinf]           = result[:Clinf] * result[:MRTinf]

        else
            result[:Kel]             = NaN
            result[:HL]              = NaN
            result[:AUCinf]          = NaN
            result[:LZint]           = NaN
        end
        #-----------------------------------------------------------------------
        #Steady-state
        tautime    = data.dosetime.time + data.dosetime.tau
        if data.dosetime.tau > 0
            eaucpartl  = eaumcpartl = 0.0
            if tautime < time[end]
                if tautime > result[:Tmax] aftertmax = true else aftertmax = false end
                result[:Ctau] = cpredict(time[ncae], time[ncae + 1], tautime, obs[ncae], obs[ncae + 1], intp)
                eaucpartl, eaumcpartl = aucpart(time[ncae], tautime, obs[ncae], result[:Ctau], calcm, aftertmax, data.dosetime.time)
                #remoove part after tau
                if ncae < result[:Obsnum] - 1 pmask[ncae:end] .= false end
            elseif tautime > time[end] && result[:LZint] !== NaN
                #extrapolation
                result[:Ctau] = exp(result[:LZint] + result[:LZ] * tautime)
                eaucpartl, eaumcpartl = aucpart(time[ncae], tautime, obs[ncae], result[:Ctau], calcm, true, data.dosetime.time)
            else
                result[:Ctau] = obs[end]
            end
            result[:AUCtau]   = sum(view(aucpartl, pmask))  + eaucpartl  + doseaucpart
            result[:AUMCtau]  = sum(view(aumcpartl, pmask)) + eaumcpartl + doseaumcpart
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
            #If Kel calculated
            if result[:Kel] !== NaN
                result[:Accind]   = 1 / (1 - (exp(-result[:Kel] * data.dosetime.tau)))
            else
                result[:Accind]   = NaN
            end
            result[:MRTtauinf]       = (result[:AUMCtau] + data.dosetime.tau * (result[:AUCinf] - result[:AUCtau])) / result[:AUCtau]
            result[:Vztau]           = data.dosetime.dose / result[:AUCtau] / result[:Kel]
            result[:Cltau]           = data.dosetime.dose / result[:AUCtau]
            result[:Vsstau]          = result[:Cltau] * result[:MRTtauinf]
        end
        if verbose
            aucpartlsum  = similar(aucpartl)
            aumcpartlsum = similar(aumcpartl)
            for i = 1:length(aucpartl)
                aucpartlsum[i]  = sum(view(aucpartl, 1:i))
                aumcpartlsum[i] = sum(view(aumcpartl, 1:i))
            end
            astx    = Vector{String}(undef, length(time))
            astx[1] = ""
            for i = 1:length(pmask)
                if pmask[i] astx[i+1] = "Yes" else astx[i+1] = "No" end
            end

            aucpartlsum .+= doseaucpart

            aumcpartlsum .+= doseaumcpart
            if  doseaucpart > 0 || doseaumcpart > 0  astx[1] = "D" end
            mx = hcat(time, obs, round.(vcat([doseaucpart], aucpartl), digits = 3),  round.(vcat([doseaucpart], aucpartlsum), digits = 3), round.(vcat([doseaumcpart], aumcpartl), digits = 3),  round.(vcat([doseaumcpart], aumcpartlsum), digits = 3), astx)
            mx = vcat(permutedims(["Time", "Concentrtion", "AUC", "AUC (cumulate)", "AUMC", "AUMC (cumulate)", "Include"]), mx)
            printmatrix(io, mx)
            println(io, "")
            println(io, "Cdose: $(result[:Cdose]), Dose time: $(data.dosetime.time)")
            if data.dosetime.time > time[1]
                println("Dose AUC  part: $(doseaucpart)")
                println("Dose AUMC part: $(doseaumcpart)")
            end
            println(io, "")
            if tautime < time[end] && tautime > 0
                println(io, "Tau + dosetime is less then end time. Interpolation used.")
                println(io, "Interpolation between: $(time[ncae]) - $( time[ncae + 1]), method: $(intp)")
                println(io, "Ctau: $(result[:Ctau])")
                println(io, "AUC  final part: $(eaucpartl)")
                println(io, "AUMC final part: $(eaumcpartl)")
            end
        end
        result[:Tmaxn] += rmn #!
        #-----------------------------------------------------------------------
        return PKPDProfile(data, result; method = calcm)
    end
"""
    nca!(data::DataSet{PKSubject}; adm = :ev, calcm = :lint, intp = :lint,
        verbose = false, warn = true, io::IO = stdout, sort = nothing)

* `sort` - only for this type of subjects.

Pharmacokinetics non-compartment analysis for PK subjects DataSet.
"""
function nca!(data::DataSet{PKSubject}; adm = :ev, calcm = :lint, intp = :lint, verbose = false, warn = true, io::IO = stdout, sort = nothing)
        results = Array{PKPDProfile, 1}(undef, 0)
        if isnothing(sort)
            for i = 1:length(data.data)
                push!(results, nca!(data.data[i]; adm = adm, calcm = calcm, intp = intp, verbose = verbose, warn = warn, io = io))
            end
        else
            for i = 1:length(data.data)
                if sort ⊆  data.data[i].sort
                    push!(results, nca!(data.data[i]; adm = adm, calcm = calcm, intp = intp, verbose = verbose, warn = warn, io = io))
                end
            end
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
        result[:RMIN] = minimum(data.obs)
        for i = 2:length(data) #BASELINE
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
        return PKPDProfile(data, result; method = :lint)
    end
"""
    nca!(data::PDSubject; verbose = false, io::IO = stdout)::PKPDProfile{PDSubject}

Pharmacodynamics non-compartment analysis for PD subjects set.
"""
function nca!(data::DataSet{PDSubject}; verbose = false, io::IO = stdout)
        results = Array{PKPDProfile, 1}(undef, 0)
        for i = 1:length(data.data)
            push!(results, nca!(data.data[i]; verbose = verbose, io = io))
        end
        return DataSet(results)
end


"""
    nca!(data::UPKSubject; verbose = false, io::IO = stdout)::PKPDProfile{UPKSubject}

Pharmacodynamics non-compartment analysis for urine data.
"""
function nca!(data::UPKSubject; verbose = false, io::IO = stdout)::PKPDProfile{UPKSubject}

    result   = Dict()
    #time  = (data.stime .+ data.etime) / 2
    #exrate   = (data.conc .* data.vol) ./  (data.etime .- data.stime)
    result[:maxrate], result[:mTmax], result[:mTmaxn] = ctmax(data.time, data.obs)
    result[:ar]      = sum(data.conc .* data.vol)
    result[:volume]  = sum(data.vol)

    #result[:AUCrate]
    #result[:tmaxrate]
    #result[:arp]
    #result[:Kel]
    #result[:HL]
    return PKPDProfile(data, result; method = :upk)
end
"""
    nca!(data::DataSet{UPKSubject}; verbose = false, io::IO = stdout)

Pharmacodynamics non-compartment analysis for PD subjects set.
"""
function nca!(data::DataSet{UPKSubject}; verbose = false, io::IO = stdout)
    results = Vector{PKPDProfile}(undef, 0)
    for i = 1:length(data)
        push!(results, nca!(data[i]; verbose = verbose, io = io))
    end
    return DataSet(results)
end
#=
stime::Vector
etime::Vector
obs::Vector
vol::Vector
=#
#-------------------------------------------------------------------------------


function getdatai(data, sort, cols, func; sortby = nothing)
    sortlist = unique(data[!, sort])
    for i = 1:size(sortlist, 1)
        inds = Vector{Int}(undef, 0)
        for c = 1:size(data, 1) #For each line in data
            if data[c, sort] == sortlist[i,:]
                push!(inds, c)
            end
        end
        datai = data[inds, cols]
        if sortby !== nothing
            sort!(datai, sortby)
        end
        func(datai, sortlist[i,:])
    end
end

#-------------------------------------------------------------------------------

function tryfloatparse!(x)
    for i = 1:length(x)
        if typeof(x[i]) <: AbstractString
            try
                x[i] = parse(Float64, x[i])
            catch
                x[i] = NaN
            end
        else
            try
                if x[i] === missing
                    x[i] = NaN
                else
                    x[i] = float(x[i])
                end
            catch
                x[i] = NaN
            end
        end
    end
    return convert(Vector{Float64}, x)
end

"""
    pkimport(data, sort::Array, rule::LimitRule; conc::Symbol, time::Symbol)

Pharmacokinetics data import.

- data - sourece data;
- sort - sorting columns;
- rule - applied LimitRule.

- conc - concentration column;
- time - time column.
"""
function pkimport(data, sort::Array, rule::LimitRule; conc::Symbol, time::Symbol)::DataSet
    results  = Array{PKSubject, 1}(undef, 0)
    if !(eltype(data[!, conc]) <: Real) throw(ArgumentError("Type of concentration data not <: Real!"))end
    getdatai(data, sort, [conc, time], (x, y) -> push!(results, PKSubject(x[!, time], x[!, conc], sort = Dict(sort .=> collect(y)))); sortby = time)
    results = DataSet(results)
    applyncarule!(results, rule)
    return results
end
"""
    pkimport(data, sort::Array; time::Symbol, conc::Symbol)

Pharmacokinetics data import.

- data - sourece data;
- sort - sorting columns.

- conc - concentration column;
- time - time column.
"""
function pkimport(data, sort::Array; time::Symbol, conc::Symbol)::DataSet
    rule = LimitRule()
    return pkimport(data, sort, rule; conc = conc, time = time)
end
"""
    pkimport(data, sort::Array; time::Symbol, conc::Symbol)

Pharmacokinetics data import  for one subject.

- data - sourece data;
- rule - applied LimitRule.

- conc - concentration column;
- time - time column.
"""
function pkimport(data, rule::LimitRule; time::Symbol, conc::Symbol)::DataSet
        #rule = LimitRule()
    datai = ncarule!(copy(data[!,[time, conc]]), conc, time, rule)
    return DataSet([PKSubject(datai[!, time], datai[!, conc])])
end
"""
    pkimport(data, sort::Array; time::Symbol, conc::Symbol)

Pharmacokinetics data import  for one subject.

- data - sourece data;

- conc - concentration column;
- time - time column.
"""
function pkimport(data; time::Symbol, conc::Symbol)::DataSet
    datai = sort(data[!,[time, conc]], time)
    return DataSet([PKSubject(datai[!, time], datai[!, conc])])
end

"""
    pkimport(time::Vector, conc::Vector; sort = Dict())::PKSubject

Pharmacokinetics data import.

- time - time vector;
- conc - concentration vector.

"""
function pkimport(time::Vector, conc::Vector; sort = Dict())::PKSubject
     PKSubject(time, conc; sort = sort)
end
    #---------------------------------------------------------------------------
"""
    pdimport(data, sort::Array; resp::Symbol, time::Symbol,
        bl::Real = 0, th::Real = NaN)

Pharmacodynamics data import.

- data - sourece data;
- sort - sorting columns;

- resp - responce column;
- time - time column
- bl - baseline;
- th - treashold.
"""
function pdimport(data, sort::Array; resp::Symbol, time::Symbol, bl::Real = 0, th::Real = NaN)::DataSet
    sortlist = unique(data[:, sort])
    results  = Array{PDSubject, 1}(undef, 0)
    getdatai(data, sort, [resp, time], (x, y) -> push!(results, PDSubject(x[!, time], x[!, resp], bl, th, sort = Dict(sort .=> collect(y)))); sortby = time)
    return DataSet(results)
end

function pdimport(data; resp::Symbol, time::Symbol, bl = 0, th = NaN)
    datai = sort(data[!,[time, resp]], time)
    return DataSet([PDSubject(datai[!, time], datai[!, resp], bl, th)])
end
#-------------------------------------------------------------------------------
function upkimport(data, sort::Array; stime::Symbol, etime::Symbol, conc::Symbol, vol::Symbol)::DataSet
    results  = Array{UPKSubject, 1}(undef, 0)
    getdatai(data, sort, [stime, etime, conc, vol], (x, y) -> push!(results, UPKSubject(x[!, stime], x[!, etime], x[!, conc], x[!, vol]; sort =Dict(sort .=> collect(y)))); sortby = stime)
    results = DataSet(results)
    return results
end

#-------------------------------------------------------------------------------
"""
    applyncarule!(data::PKSubject, rule::LimitRule)

Apply rule to PK subject .

STEP 1 (NaN step): replace all NaN values with nan reyword value (if nan !== NaN);
STEP 2 (LLOQ step): replace values below lloq with btmax value if this value befor Tmax or with atmax if this value after Tmax (if lloq !== NaN);
STEP 3 (remove NaN): rm == true, then remove all NaN values;
"""
function applyncarule!(data::PKSubject, rule::LimitRule)
    cmax, tmax, tmaxn = ctmax(data)
    #NaN Rule
    obsn = length(data)
    if rule.nan !== NaN
        for i = 1:length(data)
            if data.obs[i] === NaN || data.obs[i] === missing
                data.obs[i] = rule.nan
            end
        end
    end
    #LLOQ rule
    if rule.lloq !== NaN
        for i = 1:obsn
            if data.obs[i] <= rule.lloq
                if i <= tmaxn
                    data.obs[i] = rule.btmax
                else
                    data.obs[i] = rule.atmax
                end
            end
        end
    end
    #NaN Remove rule
    if rule.rm
        inds = findall(isnanormissing, data.obs)
        deleteat!(data.time, inds)
        deleteat!(data.obs, inds)
        #filterv   = data.obs .!== NaN
        #data.time = data.time[filterv]
        #data.obs  = data.obs[filterv]
    end
end

function applyncarule!(data::DataSet{PKSubject}, rule::LimitRule; sort = nothing)
    if isnothing(sort)
        for i in data
            applyncarule!(i, rule)
        end
    else
        for i in data
            if sort ⊆ i.sort
                applyncarule!(i, rule)
            end
        end
    end
    data
end

#---------------------------------------------------------------------------
"""
    setelimrange!(data::PKSubject, range::ElimRange; kelauto = false)

Set range for elimination parameters calculation for subject.

data - PK subject;
range - ElimRange object;
kelauto - set kelauto value.

Set kelauto `false` if  kelend > kelstart > 0.
"""
function setelimrange!(data::PKSubject, range::ElimRange; kelauto = false)
    if range.kelend > length(data) throw(ArgumentError("Kel endpoint out of range")) end
    if range.kelend > range.kelstart > 0 setkelauto!(data, kelauto) end
    data.kelrange = range
end
"""
    setelimrange!(data::DataSet{PKSubject}, range::ElimRange)

Set range for elimination parameters calculation for DataSet.

data - DataSet of PK subject;
range - ElimRange object.
"""
function setelimrange!(data::DataSet{PKSubject}, range::ElimRange)
    for i = 1:length(data)
        setelimrange!(data[i], range)
    end
    data
end
"""
    setelimrange!(data::DataSet{PKSubject}, range::ElimRange, subj::Array{Int,1})

Set range for elimination parameters calculation for DataSet.

data - DataSet of PK subject;
range - ElimRange object;
subj - subjects.
"""
function setelimrange!(data::DataSet{PKSubject}, range::ElimRange, subj::Array{Int,1})
    for i in subj
        setelimrange!(data[i], range)
    end
    data
end
"""
    setelimrange!(data::DataSet{PKSubject}, range::ElimRange, subj::Int)

Set range for elimination parameters calculation for DataSet.

data - DataSet of PK subject;
range - ElimRange object;
subj - subject.
"""
function setelimrange!(data::DataSet{PKSubject}, range::ElimRange, subj::Int)
    setelimrange!(data[subj], range)
    data
end
"""
    setelimrange!(data::DataSet{PKSubject}, range::ElimRange, sort::Dict; kelauto = false)

Set range for elimination parameters calculation for DataSet.

data - DataSet of PK subject;
range - ElimRange object;
sort - subject sort.
"""
function setelimrange!(data::DataSet{PKSubject}, range::ElimRange, sort::Dict; kelauto = false)
    for i = 1:length(data)
        if sort ⊆ data[i].sort setelimrange!(data[i], range; kelauto = kelauto) end
    end
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
function getkelauto(data::PKSubject)
    return data.kelauto
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
function applyelimrange!(data::DataSet{PKPDProfile}, range::ElimRange)
    for i = 1:length(data)
        applyelimrange!(data[i], range)
    end
    data
end
function applyelimrange!(data::DataSet{PKPDProfile}, range::ElimRange, subj::Array{Int,1})
    for i = 1:length(data)
        if i ∈ subj applyelimrange!(data[i], range) end
    end
    data
end
function applyelimrange!(data::DataSet{PKPDProfile}, range::ElimRange, subj::Int)
    applyelimrange!(data[subj], range)
    data
end
function applyelimrange!(data::DataSet{PKPDProfile}, range::ElimRange, sort::Dict)
    for i = 1:length(data)
        if sort ⊆ data[i].subject.sort
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
"""
    setdosetime!(data::DataSet{PKSubject}, dosetime::DoseTime)

Set dosing time and tau for all subjects.

data - dataset of PK subjects;
dosetime - DoseTime object.
"""
function setdosetime!(data::DataSet{PKSubject}, dosetime::DoseTime)
    for i = 1:length(data)
        setdosetime!(data[i], dosetime)
    end
    data
end
"""
    setdosetime!(data::DataSet{PKSubject}, dosetime::DoseTime, sort::Dict)

Set dosing time and tau for  subjects in sort dict.

data - dataset of PK subjects;
dosetime - DoseTime object.
"""
function setdosetime!(data::DataSet{PKSubject}, dosetime::DoseTime, sort::Dict)
    for i = 1:length(data)
        if sort ⊆ data[i].sort setdosetime!(data[i], dosetime) end
    end
    data
end
"""
    findfirst(sort::Dict, data::DataSet{PKSubject})

The first item in the dataset that satisfies the condition - sort dictionary.
"""
function findfirst(sort::Dict, data::DataSet{PKSubject})
    for i = 1:length(data)
        if sort ⊆ data[i].sort
            return data[i]
        end
    end
end

"""
    dosetime(data::PKSubject)

Return dosing time and tau for subject.

data - PK subject.
"""
function dosetime(data::PKSubject)
    data.dosetime
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
        if sort ⊆ data[i].sort setth!(data[i], th) end
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
        if sort ⊆ data[i].sort setbl!(data[i], bl) end
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


#-------------------------------------------------------------------------------
# Param functions

function cmax(data::PKSubject)
    return maximum(data.obs)
end

function auc()
end
#----
function datatable(data; stack = true)
    if stack return datatable_st(data) else return datatable_unst(data) end
end
function datatable_st(data) where T

        sortcols = Vector{Any}(undef, 0)
        sorttype = Dict()
        for i = 1:length(data)
            if length(data[i].sort) > 0
                for k in keys(data[i].sort)
                    if k ∉ sortcols push!(sortcols, k) end
                    if k ∉ keys(sorttype) sorttype[k] = typeof(data[i].sort[k]) else sorttype[k] = promote_type(sorttype[k], typeof(data[i].sort[k])) end
                end
            end
        end
        pname = :Parameter
        vname = :Value

        for i = 1:length(sortcols)
            if sortcols[i] == pname pname = Symbol(string(pname)*"_") end
            if sortcols[i] == vname pname = Symbol(string(vname)*"_") end
        end
        t = Vector{Vector}(undef, length(sortcols))
        for i = 1:length(t)
            t[i] = Vector{sorttype[sortcols[i]]}(undef, 0)
        end
        nt = NamedTuple{Tuple(append!(copy(sortcols), [pname, vname]))}(Tuple(append!(copy(t), [Vector{Symbol}(undef, 0), Vector{Float64}(undef, 0)])))
        for i = 1:length(data)
            sort = Dict()
            for c = 1:length(sortcols)
                if sortcols[c] ∈ keys(data[i].sort) sort[sortcols[c]] = data[i].sort[sortcols[c]] else sort[sortcols[c]] = missing end
            end
            for k in keys(data[i].result)
                for s in keys(sort)
                    push!(nt[s], sort[s])
                end
                push!(nt[pname], k)
                push!(nt[vname], data[i].result[k])
            end
        end
        return nt
end
function datatable_unst(data)
end

"""
    auc_sparse(data::PKSubject)

AUC for sparse data.

```math
w_1 = (t_2 - t_1) / 2
w_j = (t_{j+1} - t_{j-1}) / 2  (2 \\leq j \\leq J - 1)
w_J = (t_J - t_{J-1}) / 2

AUC = \\sum_{j=1}^J \\mu_j w_j
```

where `math \\mu_j` is the mean drug concentration at time t.
"""
function auc_sparse(data::PKSubject)
    wts = Vector{Float64}(undef, length(data))
    wts[1]   = (data.time[2] - data.time[1])/2.0
    wts[end] = (data.time[end] - data.time[end-1])/2.0
    if length(data) > 2
        for i = 2:length(data)-1
            wts[i] = (data.time[i+1] - data.time[i-1])/2.0
        end
    end
    return data.obs'*wts
end
