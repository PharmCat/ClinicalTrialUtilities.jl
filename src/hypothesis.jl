# Clinical Trial Utilities
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)

struct Equivalence <:AbstractEquivalenceHypothesis
    llim::Real          #Lower lmit for Test group
    ulim::Real          #Upper lmit for Test group
    alpha::Real
    #diff::Real          #Margin difference
    function Equivalence(llim, ulim, alpha; bio::Bool = false)
        if llim == ulim throw(ArgumentError("llim == ulim!")) end
        if bio return Bioequivalence(llim, ulim, alpha) end
        new(llim, ulim, alpha)::Equivalence
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
    alpha::Real
end

struct Equality <: AbstractHypothesis
    val::Real
    alpha::Real
    function Equality(alpha)
        new(0, alpha)::Equality
    end
    function Equality(val, alpha)
        new(val, alpha)::Equality
    end
end
function refval(h::Equality)
    h.val
end

struct Superiority <: AbstractHypothesis
    llim::Real          #Lower lmit for Test group
    diff::Real          #Margin difference
    alpha::Real
end
function refval(h::Superiority)
    h.llim - h.diff
end


function checkhyp(h::Superiority, param; method::Symbol = :default)
    ci = confint(param; level = 1.0 - h.alpha * 2, method = method)
    ci.lower > h.diff
end
function checkhyp(h::Equivalence, param; method::Symbol = :default)
    ci = confint(param; level = 1.0 - h.alpha * 2, method = method)
    ci.lower > h.llim && ci.upper < h.ulim
end
function checkhyp(h::Equality, param; method::Symbol = :default)
    ci = confint(param; level = 1.0 - h.alpha, method = method)
    ci.lower > h.val || ci.upper < h.val
end


struct McNemars <: AbstractHypothesis
    alpha::Real
end
