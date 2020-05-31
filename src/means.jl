struct Mean{Bool} <: AbstractMean
    m::Real
    sd::Real
    n::Int
    function Mean(m::Real, sd::Real, n::Int)
        if sd ≤ 0.0 throw(ArgumentError("SD can't be ≤ 0.0!")) end
        new{true}(m, sd, n)::Mean
    end
    function Mean(m::Real, sd::Real)
        if sd ≤ 0.0 throw(ArgumentError("SD can't be ≤ 0.0!")) end
        new{false}(m, sd, 0)
    end
    function Mean(m::Real)
        new{false}(m, NaN, 0)
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
