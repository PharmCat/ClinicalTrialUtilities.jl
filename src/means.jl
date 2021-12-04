struct Mean{Bool} <: AbstractMean
    m::Float64
    sd::Float64
    n::Int
    function Mean(m, sd, n::Int)
        if sd ≤ 0.0 throw(ArgumentError("SD can't be ≤ 0.0!")) end
        new{true}(m, sd, n)::Mean
    end
    function Mean(m, sd)
        if sd ≤ 0.0 throw(ArgumentError("SD can't be ≤ 0.0!")) end
        new{false}(m, sd, 0)
    end
    function Mean(m::Real)
        new{false}(m, NaN, 0)
    end
    function Mean(v::Vector{T}) where T
        new{true}(mean(v), sqrt(var(v)), length(v))
    end
end

function getval(m::Mean)
    return m.m
end
"""
    Mean comparison hypothesis testing
    a - Test group value
    b - Reference group value {Mean} type, or reference value for one group test {Real}
"""
struct DiffMean{Bool} <: AbstractCompositeMean
    a::Mean
    b::Mean
    function DiffMean(a::Mean, b::Mean)
        if typeof(a) <: Mean{false} || typeof(b) <: Mean{false}  t = false else t = true end
        new{t}(a, b)
    end
end
