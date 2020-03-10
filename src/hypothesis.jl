# Clinical Trial Utilities
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)

struct Equivalence <:AbstractEquivalenceHypothesis
    llim::Real          #Lower lmit for Test group
    ulim::Real          #Upper lmit for Test group
    #diff::Real          #Margin difference
    function Equivalence(llim, ulim; bio::Bool = false)
        if llim == ulim throw(ArgumentError("llim == ulim!")) end
        if bio return Bioequivalence(llim, ulim) end
        new(llim, ulim)::Equivalence
    end
end
function mdiff(h::Equivalence)::AbstractFloat
    (h.ulim - h.llim)/2
end

struct Bioequivalence <: AbstractEquivalenceHypothesis
    llim::Real          #Lower lmit for Test group
    ulim::Real          #Upper lmit for Test group
end

struct Equality <: AbstractHypothesis
    val::Real
    function Equality()
        new(0)::Equality
    end
end
struct Superiority <: AbstractHypothesis
    llim::Real          #Lower lmit for Test group
    diff::Real          #Margin difference
end
struct McNemars <: AbstractHypothesis end
