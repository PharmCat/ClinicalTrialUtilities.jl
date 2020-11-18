abstract type AbstractTest end
abstract type AbstractConTab <: AbstractData end

struct ConTab{Int32, Int32} <: AbstractConTab
    tab::Matrix{Int}
    row::Vector
    col::Vector
    sort::Dict

    function ConTab(m::Matrix{T}) where T <: Int
        row = Vector{Symbol}(undef, size(m, 1)) .= Symbol("")
        col = Vector{Symbol}(undef, size(m, 2)) .= Symbol("")
        new{size(m, 1), size(m, 2)}(m, row, col, Dict())
    end
    function ConTab(m::Matrix{T}, row, col) where T <: Int
        new{size(m, 1), size(m, 2)}(m, row, col, Dict())
    end
    function ConTab(m::Matrix{T}, row, col, sort) where T <: Int
        new{size(m, 1), size(m, 2)}(m, row, col, sort)
    end
end

struct McnmConTab <: AbstractConTab
    tab::Matrix{Int}
    sort::Dict
    function McnmConTab(m::Matrix{Int})
        new(m, Dict())
    end
    function McnmConTab(m::Matrix{Int}, sort)
        new(m, sort)
    end
end

function Base.getproperty(obj::McnmConTab, attr::Symbol)
    if attr == :tab
        return getfield(obj, :tab)
    elseif attr == :a
        return getfield(obj, :tab)[1,1]
    elseif attr == :b
        return getfield(obj, :tab)[1,2]
    elseif attr == :c
        return getfield(obj, :tab)[2,1]
    elseif attr == :d
        return getfield(obj, :tab)[2,2]
    else
        error("No field!")
    end
end

function Base.getproperty(obj::ConTab{2, 2}, attr::Symbol)
    if attr == :tab
        return getfield(obj, :tab)
    elseif attr == :a
        return getfield(obj, :tab)[1,1]
    elseif attr == :b
        return getfield(obj, :tab)[1,2]
    elseif attr == :c
        return getfield(obj, :tab)[2,1]
    elseif attr == :d
        return getfield(obj, :tab)[2,2]
    else
        return getfield(obj, attr)
    end
end

struct Freque
    val
    n
end

"""
    freque(data::DataFrame; vars::Symbol, alpha = 0.05)::DataFrame

Frequencies.
"""
function freque(data::DataFrame; vars::Symbol, alpha = 0.05)::DataFrame
    result = DataFrame(value = Any[], n = Int[], p = Float64[], cil = Float64[], ciu = Float64[])
    list   = unique(data[:, vars])
    n      = length(data[:, vars])
    for i in list
        ne = count(x -> (x == i), data[:, vars])
        pe = ne/n
        ci = ClinicalTrialUtilities.propci(ne, n, alpha=alpha, method=:wald)
        push!(result, [i, ne, pe, ci.lower, ci.upper])
    end
    return result
end
"""
    contab(data::DataFrame; row::Symbol, col::Symbol, sort = Dict())::ConTab

Make contingency table.
"""
function contab(data::DataFrame; row::Symbol, col::Symbol, sort = Dict())::ConTab
    clist = unique(data[:, col])
    rlist = unique(data[:, row])
    cn    = length(clist)
    rn    = length(rlist)
    dfs   = Matrix{Int}(undef, rn, cn)
    for ri = 1:rn
        rowl  = data[data[:, row] .== rlist[ri], col]
        for ci = 1:cn
            cnt = count(x -> x == clist[ci], rowl)
            dfs[ri, ci] = cnt
        end
    end
    return ConTab(dfs, rlist, clist, sort)
end
"""
    contab(data::DataFrame, sort; row::Symbol, col::Symbol)

Make contingency tables set.
"""
function contab(data::DataFrame, sort; row::Symbol, col::Symbol)
    slist  = unique(data[:, sort])
    clist  = unique(data[:, col])
    rlist  = unique(data[:, row])
    rcv    = [row, col]
    result = Vector{ConTab}(undef, size(slist, 1))
    tempdf = DataFrame(row = Vector{eltype(data[!, row])}(undef, 0), col = Vector{eltype(data[!, col])}(undef, 0))
    rename!(tempdf, rcv)
    for si = 1:size(slist, 1)
        if size(tempdf, 1) > 0 deleterows!(tempdf, 1:size(tempdf, 1)) end
        for i = 1:size(data, 1)
            if data[i, sort] == slist[si, :]
                push!(tempdf, data[i, rcv])
            end
        end
        result[si] = contab(tempdf, row = row, col = col, sort = Dict(sort .=> collect(slist[si, :])))
    end
    return DataSet(result)
end

"""
    mcnmcontab
"""
function mcnmcontab(data::DataFrame; row::Symbol, col::Symbol, sort = Dict())::McnmConTab
    clist = unique(data[:, col])
    rlist = unique(data[:, row])
    cn    = length(clist)
    rn    = length(rlist)
    dfs   = Matrix{Int}(undef, rn, cn)
    for ri = 1:rn
        rowl  = data[data[:, row] .== rlist[ri], col]
        for ci = 1:cn
            cnt = count(x -> x == clist[ci], rowl)
            dfs[ri, ci] = cnt
        end
    end
    return McnmConTab(dfs, sort)
end
"""
    mcnmcontab
"""
function mcnmcontab(data::DataFrame, sort; row::Symbol, col::Symbol)
    slist  = unique(data[:, sort])
    clist  = unique(data[:, col])
    rlist  = unique(data[:, row])
    rcv    = [row, col]
    result = Vector{McnmConTab}(undef, size(slist, 1))
    tempdf = DataFrame(row = Vector{eltype(data[!, row])}(undef, 0), col = Vector{eltype(data[!, col])}(undef, 0))
    rename!(tempdf, rcv)
    for si = 1:size(slist, 1)
        if size(tempdf, 1) > 0 deleterows!(tempdf, 1:size(tempdf, 1)) end
        for i = 1:size(data, 1)
            if data[i, sort] == slist[si, :]
                push!(tempdf, data[i, rcv])
            end
        end
        result[si] = mcnmcontab(tempdf, row = row, col = col, sort = Dict(sort .=> collect(slist[si, :])))
    end
    return DataSet(result)
end


function pirson(a::Matrix{Int})
    n   = length(a[:,1])
    m   = length(a[1,:])
    tm  = sum(a, dims=1)[1,:]
    tn  = sum(a, dims=2)[:,1]
    num = sum(tm)
    ae  = Array{Float64, 2}(undef, n, m)
    for im = 1:m
        for in = 1:n
            ae[in, im] = tn[in]*tm[im]/num
        end
    end
    chsq  = sum(((a .- ae) .^2 ) ./ ae)
    chsqy = sum(((abs.((a .- ae)) .- 0.5) .^2 ) ./ ae)
    ml    = 2 * sum( a .* log.( a ./ ae ))
    df    = (n - 1)*(m - 1)
    ϕ     = sqrt(chsq / (num*(n - 1)*(m - 1)))
    C     = sqrt(chsq/(chsq+num))
    K     = sqrt(chsq/num/sqrt((n - 1)*(m - 1)))
    return chsq, chsqy, ml, 1-cdf(Chisq(df), chsq)
end
function pirson(a::ConTab)
    return (pirson(a.tab))
end

function fisher(a::Matrix{Int})
    dist  = Hypergeometric(sum(a[1, :]), sum(a[2, :]), sum(a[:, 1]))
    l = cdf(dist, a[1, 1])
    r = ccdf(dist, a[1, 1]-1)
    value = min(2 * min(r, l), 1.0)
    value, l, r
end
function fisher(t::ConTab{2, 2})
    fisher(t.tab)
end


function mcnmtest(a::Matrix{Int}; cc = false)
    if cc cc = 1 else cc = 0 end
    (abs(a[1,2] - a[2,1]) - cc) ^ 2 / (a[1,2] + a[2,1])
end


function mcnmtest(a::McnmConTab; cc = false)
    return mcnmtest(a.tab; cc = cc)
end
#=
function StatsBase.confint(t::ConTab{2, 2})
end
=#

struct McnmTest
    chisq
    function McnmTest(ct::McnmConTab)
        new(mcnmtest(ct.tab))
    end
end

struct MetaProp <: AbstractTest
    model::Symbol
    type::Symbol
    esti::Vector
    vari::Vector
    wts::Vector
    est::Real
    var::Real
    k::Int
    q::Real
    tausq::Real
    chisq::Real
    function MetaProp(model::Symbol, t::Symbol, esti::Vector, vari::Vector, wts::Vector, est::Real, var::Real, k::Int, q::Real, tausq, chisq)
        new(model, t, esti, vari, wts, est, var, k, q, tausq, chisq)
    end
end

function tabdiff(tab::McnmConTab, adj)
    a = tab.a + adj
    b = tab.b + adj
    c = tab.c + adj
    d = tab.d + adj
    n        = a + b + c + d
    return (a + b)/n - (a + c)/n
end
function tabdiffvar(tab::McnmConTab, adj)
    a = tab.a + adj
    b = tab.b + adj
    c = tab.c + adj
    d = tab.d + adj
    n        = a + b + c + d
    return ((a + d) * (b + c) + 4 * b * c) / n^3
end

function tabdiff(tab::ConTab{2,2}, adj)
    a = tab.a + adj
    b = tab.b + adj
    c = tab.c + adj
    d = tab.d + adj
    pt = a / (a + b)
    pc = c / (c + d)
    return  pt - pc
end
function tabdiffvar(tab::ConTab{2,2}, adj)
    a = tab.a + adj
    b = tab.b + adj
    c = tab.c + adj
    d = tab.d + adj
    pt = a / (a + b)
    pc = c / (c + d)
    return pt * (1 - pt) / (a + b) + pc * (1 - pc) / (c + d)
end

function tabor(tab::McnmConTab, adj)
    b = tab.b + adj
    c = tab.c + adj
    return b / c
end
function taborvar(tab::McnmConTab, adj)
    b = tab.b + adj
    c = tab.c + adj
    return  1 / b + 1 / c
end

function tabor(tab::ConTab{2,2}, adj)
    a = tab.a + adj
    b = tab.b + adj
    c = tab.c + adj
    d = tab.d + adj
    pt = a / (a + b)
    pc = c / (c + d)
    return (pt / (1 - pt)) / (pc / (1 - pc))
end
function taborvar(tab::ConTab{2,2}, adj)
    a = tab.a + adj
    b = tab.b + adj
    c = tab.c + adj
    d = tab.d + adj
    return 1 / a + 1 / b + 1 / c + 1 / d
end

function tabrr(tab::McnmConTab, adj)
    a = tab.a + adj
    b = tab.b + adj
    c = tab.c + adj
    return (a + b)/(a + c)
end
function tabrrvar(tab::McnmConTab, adj)
    a = tab.a + adj
    b = tab.b + adj
    c = tab.c + adj
    return (b + c) / (a + c) / (a + b)
end

function tabrr(tab::ConTab{2,2}, adj)
    a = tab.a + adj
    b = tab.b + adj
    c = tab.c + adj
    d = tab.d + adj
    pt = a / (a + b)
    pc = c / (c + d)
    return pt / pc
end
function tabrrvar(tab::ConTab{2,2}, adj)
    a = tab.a + adj
    b = tab.b + adj
    c = tab.c + adj
    d = tab.d + adj
    return 1 / a - 1 / (a + b) + 1 / c - 1 /(c + d)
end
#=
function cmhweight(tab::ConTab{2,2})
    return (tab.a + tab.b) * (tab.c + tab.d) / (tab.a + tab.b + tab.c + tab.d)

end
=#

"""
    metaprop(tab::Vector{T}; type::Symbol, model::Symbol = :fixed, zeroadj::Real = 0, tau::Symbol =:dl)::MetaProp where T <: AbstractConTab

Meta-analysis for 2x2 tables.

tab: vectro of ConTab{2,2} or McnmConTab;

Inverce Variance method used to get variance estimate for fixed effect.

type - type of measure:
- :rr
- :or/
- :diff

model:
- :fixed
- :random

zeroadj - zero adjustment value for all cells;

tau - τ² calculation method:
- :dl
- :ho
- :hm

"""
function metaprop(tab::Vector{T}; type::Symbol, model::Symbol = :fixed, zeroadj::Real = 0, tau::Symbol =:dl)::MetaProp where T <: AbstractConTab
    k = length(tab)

    if type == :diff
        θᵢ = tabdiff.(tab, zeroadj)
        vᵢ = tabdiffvar.(tab, zeroadj)
    elseif type == :or
        θᵢ = log.(tabor.(tab, zeroadj))
        vᵢ = taborvar.(tab, zeroadj)
    elseif type == :rr
        θᵢ = log.(tabrr.(tab, zeroadj))
        vᵢ = tabrrvar.(tab, zeroadj)
    else
        error("Type keyword unknown!")
    end

    wᵢ    = 1 ./ vᵢ
    v     = 1 / sum(wᵢ)

    θ     = sum(wᵢ .* θᵢ) / sum(wᵢ)
    chisq = sum(wᵢ .* (θᵢ .^ 2))
    q     = sum(wᵢ .* ((θᵢ .- θ) .^ 2))

        #Tau square for random effect

    if tau == :dl
        τ² = max(0, (q - (k - 1))/ (sum(wᵢ) - sum(wᵢ .^ 2) / sum(wᵢ)))
    elseif tau == :ho
        τ² = max(0, 1 / (k - 1) * sum((θᵢ .- mean(θᵢ)) .^ 2) - 1 / k * sum(vᵢ))
    elseif tau == :hm
        τ² = q ^ 2 / (2 * (k - 1) + q) / (sum(wᵢ) - sum(wᵢ .^ 2) / sum(wᵢ))
    elseif tau == :ml
            #in dev
        lnL(m, t, k, vi, yi) = - 1 / 2 * log(2pi) - 1 / 2 * sum(log.(vi .+ t)) - 1 / 2 * sum(((yi .- m) .^ 2)/(vi .+ t))
    else
        error("Tau keyword unknown!")
            #=
            w   = 1 / k * sum(wᵢ)
            s²  = 1 / (k - 1) * sum( wᵢ .^ 2 .- (k * w ^ 2))
            u   = (k - 1) * (w - s² / k / w)
            τ²  = max(0, (q - k + 1) / u)
            =#
    end

    if model == :fixed
        return MetaProp(model, type, θᵢ, vᵢ, wᵢ, θ, v, k, q, τ², chisq)
    end
    rwᵢ = 1 ./ (vᵢ .+ τ²)
    θ     = sum(rwᵢ .* θᵢ) / sum(rwᵢ)
    v     = 1 / sum(rwᵢ)
    return MetaProp(model, type, θᵢ, vᵢ, rwᵢ, θ, v, k, q, τ², chisq)
end
"""
    metaprop(data::DataSet{T}; type::Symbol, model::Symbol = :fixed, zeroadj::Real = 0, tau::Symbol =:dl)::MetaProp where T <: AbstractConTab

Meta-analysis for 2x2 tables.

data - DataSet of ConTab{2,2} or McnmConTab;
"""
function metaprop(data::DataSet{T}; type::Symbol, model::Symbol = :fixed, zeroadj::Real = 0, tau::Symbol =:dl)::MetaProp where T <: AbstractConTab
    metaprop(data.data; type = type, model = model, zeroadj = zeroadj, tau = tau)
end


function Base.show(io::IO, obj::MetaProp)
    println(io, "  Meta-analysis")
    println(io, "  -----------------------")
    println(io, "  Trial number (k): $(obj.k)")
    println(io, "  Model: $(obj.model)")
    println(io, "  Type: $(obj.type)")
    if obj.type == :diff
        println(io, "  Estimate: $(obj.est)")
        println(io, "  Variance: $(obj.var)")
    else
        println(io, "  Estimate: $(exp(obj.est))")
        println(io, "  Variance: $(exp(obj.var))")
    end

    println(io, "  Q: $(obj.q)")
    print(io,   "  τ²: $(obj.tausq)")
end
