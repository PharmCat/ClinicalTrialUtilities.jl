struct Proportion{Bool} <: AbstractSimpleProportion
    x::Real
    n::Real
    val::Float64
    function Proportion(x::Int, n::Int)
        if x > n throw(ArgumentError("Error: X > N!")) end
        new{true}(x, n, x/n)::Proportion
    end
    function Proportion(val::Float64)
        new{false}(NaN, NaN, val)::Proportion
    end
end
#=
struct Probability{T} <: AbstractSimpleProportion where T <: AbstractFloat
    p::T
    function Probability(p::T) where T <: AbstractFloat
        if p ≥ 1.0 || p ≤ 0.0 throw(ArgumentError("Probability can't be ≥ 1.0 or ≤ 0.0!")) end
        new{T}(p)::Probability
    end
    function Probability(p::Proportion)
        prob = p.x/p.n
        new{typeof(prob)}(prob)::Probability
    end
end
=#
function getval(p::Proportion)
    return p.val
end
function Base.:+(a::Proportion{true}, b::Proportion{true})
    return Proportion(a.x + b.x, a.n + b.n)
end
#function getval(p::Probability)
#    return p.p
#end

struct DiffProportion{T1 <: AbstractSimpleProportion, T2 <: AbstractSimpleProportion} <: AbstractCompositeProportion
    a::T1
    b::T2
    function DiffProportion(a::T, b::T) where T <: AbstractSimpleProportion
        new{typeof(a), typeof(b)}(a, b)::DiffProportion
    end
    function DiffProportion(a::T1, b::T2) where T1 <: AbstractSimpleProportion where T2 <: Real
        b = Proportion(b)
        new{typeof(a), typeof(b)}(a, b)::DiffProportion
    end
    function DiffProportion(a::T1, b::T2) where T1 <: Real where T2 <: AbstractSimpleProportion
        a = Proportion(a)
        new{typeof(a), typeof(b)}(a, b)::DiffProportion
    end
    function DiffProportion(a::T1, b::T2) where T1 <: Real where T2 <: Real
        a = Proportion(a)
        b = Proportion(b)
        new{typeof(a), typeof(b)}(a, b)::DiffProportion
    end
end

function getval(p::DiffProportion)
    return getval(p.a) - getval(p.b)
end

struct OddRatio{T1 <: AbstractSimpleProportion, T2 <: AbstractSimpleProportion} <: AbstractCompositeProportion
    a::T1
    b::T2
    val::Float64
    function OddRatio(a::T, b::T) where T <: AbstractSimpleProportion
        new{typeof(a), typeof(b)}(a, b, getval(a)/getval(b))
    end
    function OddRatio(a::T1, b::T2) where T1 <: AbstractSimpleProportion where T2 <: Real
        b = Proportion(b)
        new{typeof(a), typeof(b)}(a, b, getval(a)/getval(b))::OddRatio
    end
    function OddRatio(a::T1, b::T2) where T1 <: Real where T2 <: AbstractSimpleProportion
        a = Proportion(a)
        new{typeof(a), typeof(b)}(a, b, getval(a)/getval(b))::OddRatio
    end
    function OddRatio(a::T1, b::T2) where T1 <: Real where T2 <: Real
        a = Proportion(a)
        b = Proportion(b)
        new{typeof(a), typeof(b)}(a, b, getval(a)/getval(b))::OddRatio
    end
    function OddRatio(val)
        a = Proportion(NaN)
        b = Proportion(NaN)
        new{typeof(a), typeof(b)}(NaN, NaN, val)
    end
end
function getval(p::OddRatio)
    return p.val
end
struct RiskRatio{T <: AbstractSimpleProportion} <: AbstractCompositeProportion
    a::T
    b::T
end

struct CoxHazardRatio <: AbstractParameter
    a::Float64
    b::Float64
    p::Float64
end

struct DiffCoxHazardRatio <: AbstractCompositeProportion
    a::CoxHazardRatio
    b::CoxHazardRatio
end
