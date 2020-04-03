abstract type AbstractTest end
abstract type AbstractConTab end

struct ConTab{Int32, Int32} <: AbstractConTab
    tab::Matrix{Real}
    row::Vector
    col::Vector

    function ConTab(m::Matrix{T}) where T <: Real
        row = Vector{Symbol}(undef, size(m, 1)) .= Symbol("")
        col = Vector{Symbol}(undef, size(m, 2)) .= Symbol("")
        new{size(m, 1), size(m, 2)}(m, row, col)
    end
    function ConTab(m::Matrix{T}, row, col) where T <: Real
        new{size(m, 1), size(m, 2)}(m, row, col)
    end
end

struct McnmConTab <: AbstractConTab
    tab::Matrix{Int}
    function McnmConTab(m::Matrix{Int})
        new(m)
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
        error("No field!")
    end
end

struct Freque
    val
    n
end


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

function contab(data::DataFrame; row::Symbol, col::Symbol)::ConTab
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
    return ConTab(dfs, rlist, clist)
end


function pirson(a::Matrix{Int})
    n   = length(a[:,1])
    m   = length(a[1,:])
    tm  = sum(a, dims=1)[1,:]
    tn  = sum(a, dims=2)[:,1]
    num = sum(tm)
    ae  = Array{Real, 2}(undef, n, m)
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

function fisher(a::Matrix{Int})
    dist  = Hypergeometric(sum(a[1, :]), sum(a[2, :]), sum(a[:, 1]))
    value = min(2 * min(cdf(dist, a[1, 1]), ccdf(dist, a[1, 1])), 1.0)
end


function fisher(t::ConTab{2, 2})
    fisher(t.tab)
end


function mcnmtest(a::Matrix{Int}; cc = false)
    if cc cc = 1 else cc = 0 end
    (abs(a[1,2] - a[2,1]) - cc) ^ 2 / (a[1,2] + a[2,1])
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

struct CMHTest <: AbstractTest
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
    function CMHTest(model::Symbol, t::Symbol, esti::Vector, vari::Vector, wts::Vector, est::Real, var::Real, k::Int, q::Real, tausq, chisq)
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


function cmh(tab::Vector{T}; type::Symbol, model::Symbol = :fixed, zeroadj::Real = 0, tau::Symbol =:dl)::CMHTest where T <: AbstractConTab
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
        return CMHTest(model, type, θᵢ, vᵢ, wᵢ, θ, v, k, q, τ², chisq)
    end

    rwᵢ = 1 ./ (vᵢ .+ τ²)
    θ     = sum(rwᵢ .* θᵢ) / sum(rwᵢ)
    v     = 1 / sum(rwᵢ)

    return CMHTest(model, type, θᵢ, vᵢ, rwᵢ, θ, v, k, q, τ², chisq)
end
