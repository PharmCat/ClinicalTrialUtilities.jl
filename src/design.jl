# Clinical Trial Utilities
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)

struct Parallel <: AbstractDesign
    df::Function
    bkni::Real
    sq::Int
    function Parallel()
        new(x -> x - 2, 1.0, 2)::Parallel
    end
end
struct OneGroup <: AbstractDesign
    df::Function
    bkni::Real
    sq::Int
    function OneGroup()
        new(x -> x - 1, 1, 1)::OneGroup
    end
end
struct Crossover{Symbol} <: AbstractDesign
    df::Function
    bkni::Real
    sq::Int
    type
    function Crossover(type::Symbol)
        return Design(type)
    end
    function Crossover(df::Function, bkni::Real, sq::Int, t::Symbol)
        new{t}(df, bkni, sq, t)
    end
end

function Design(type::Symbol)::AbstractDesign
    if type == :parallel
        return Parallel()
    elseif type == :d2x2
        return Crossover(x -> x - 2, 0.5, 2, type)
    elseif type == :d2x2x3
        return Crossover(x -> 2 * x - 3, 0.375, 2, type)
    elseif type == :d2x2x4
        return Crossover(x -> 3 * x - 4, 0.25, 2, type)
    elseif type == :d2x4x4
        return  Crossover(x -> 3 * x - 4, 0.0625, 4, type)
    elseif type == :d2x3x3
        return Crossover(x -> 2 * x - 3, 1/6, 3, type)
    elseif type == :d2x4x2
        return Crossover(x -> x - 2, 0.5, 4, type)
    elseif type == :d3x3
        return Crossover(x -> 2 * x - 4, 2/9, 3, type)
    elseif type == :d4x4
        return Crossover(x -> 3 * x - 6, 1/8, 4, type)
    elseif type == :d3x6x3
        return Crossover(x -> 2 * x - 4, 1/18, 6, type)
    else throw(ArgumentError("Design type not known!")) end
end

function designtype(d::Crossover)
    return d.type
end
function designtype(d::Parallel)
    return :parallel
end

@inline function sediv(d::AbstractDesign, n::Int)
    sqa   = Array{Float64, 1}(undef, d.sq)
    sqa  .= n÷d.sq
    for i = 1:n%d.sq
        sqa[i] += 1
    end
    return sqrt(sum(1 ./ sqa)*d.bkni)
end
