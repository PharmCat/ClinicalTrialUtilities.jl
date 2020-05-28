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
function mdiff(h::Equivalence)
    (h.ulim - h.llim)/2
end
function refval(h::Equivalence)
    (h.ulim + h.llim)/2
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
    function Equality(val)
        new(val)::Equality
    end
end
function refval(h::Equality)
    h.val
end

struct Superiority <: AbstractHypothesis
    llim::Real          #Lower lmit for Test group
    diff::Real          #Margin difference
end
function refval(h::Superiority)
    h.llim - h.diff
end


function checkhyp(h::Superiority, ci::ConfInt)
    ci.lower > h.diff
end
function checkhyp(h::Equivalence, ci::ConfInt)
    ci.lower > h.llim && ci.upper < h.ulim
end
function checkhyp(h::Equality, ci::ConfInt)
    ci.lower > h.val || ci.upper < h.val
end


struct McNemars <: AbstractHypothesis end
