struct Proportion <: AbstractSimpleProportion
    x::Int
    n::Int
    function Proportion(x::Int, n::Int)
        if x > n throw(ArgumentError("Error: X > N!")) end
        new(x, n)::Proportion
    end
end
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
function getval(p::Proportion)::Float64
    return p.x/p.n
end
function getval(p::Probability)::Float64
    return p.p
end

struct DiffProportion{S <: AbstractSimpleProportion, T <: Union{Proportion, Probability, Real}}  <: AbstractCompositeProportion{T}
    a::S
    b::T
end
struct OddRatio{T <: AbstractSimpleProportion} <: AbstractCompositeProportion{T}
    a::T
    b::T
end
struct RiskRatio{T <: AbstractSimpleProportion} <: AbstractCompositeProportion{T}
    a::T
    b::T
end
