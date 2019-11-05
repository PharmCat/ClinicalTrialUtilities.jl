# Clinical Trial Utilities
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)

#ctsamplen, ctpower, besamplen, bepower

abstract type AbstractTask end
abstract type AbstractParameter end

abstract type AbstractProportion  <:  AbstractParameter end
abstract type AbstractMean  <:  AbstractParameter end

abstract type AbstractSimpleProportion <:  AbstractProportion end

abstract type AbstractCompositeProportion{T}  <:  AbstractProportion end
abstract type AbstractCompositeMean{T}  <:  AbstractMean end

abstract type AbstractObjective end
abstract type AbstractHypothesis end
abstract type AbstractEquivalenceHypothesis <: AbstractHypothesis end

abstract type AbstractDesign end

struct Parallel <: AbstractDesign
    df::Function
    bkni::Real
    sq::Int
    function Parallel()
        new(x -> x - 2, 1.0, 2)::Parallel
    end
end
struct Onegroup <: AbstractDesign
    df::Function
    bkni::Real
    sq::Int
    function Onegroup()
        new(x -> x - 1, 1, 1)::Onegroup
    end
end
struct Crossover <: AbstractDesign
    df::Function
    bkni::Real
    sq::Int
    function Crossover(type::Symbol)
        return Design(type)
    end
    function Crossover(df::Function, bkni::Real, sq::Int)
        new(df, bkni, sq)::Crossover
    end
end

function Design(type::Symbol)::AbstractDesign
    if type == :parallel
        return Parallel()
    elseif type == :d2x2
        return Crossover(x -> x - 2, 0.5, 2)
    elseif type == :d2x2x3
        return Crossover(x -> 2 * x - 3, 0.375, 2)
    elseif type == :d2x2x4
        return Crossover(x -> 3 * x - 4, 0.25, 2)
    elseif type == :d2x4x4
        return  Crossover(x -> 3 * x - 4, 0.0625, 4)
    elseif type == :d2x3x3
        return Crossover(x -> 2 * x - 3, 1/6, 3)
    elseif type == :d2x4x2
        return Crossover(x -> x - 2, 0.5, 4)
    elseif type == :d3x3
        return Crossover(x -> 2 * x - 4, 2/9, 3)
    elseif type == :d4x4
        return Crossover(x -> 3 * x - 6, 1/8, 4)
    elseif type == :d3x6x3
        return Crossover(x -> 2 * x - 4, 1/18, 6)
    else throw(ArgumentError("Design type not known!")) end
end

@inline function sediv(d::AbstractDesign, n::Int)
    sqa   = Array{Float64, 1}(undef, d.sq)
    sqa  .= n÷d.sq
    for i = 1:n%d.sq
        sqa[i] += 1
    end
    return sqrt(sum(1 ./ sqa)*d.bkni)
end

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
function mdiff(h::Equivalence)::Float64
    (h.ulim - h.llim)/2
end

struct Bioequivalence <: AbstractEquivalenceHypothesis
    llim::Real          #Lower lmit for Test group
    ulim::Real          #Upper lmit for Test group
end

struct Equality <: AbstractHypothesis
end
struct Superiority <: AbstractHypothesis
    llim::Real          #Lower lmit for Test group
    diff::Real          #Margin difference
end
struct McNemars <: AbstractHypothesis end

struct Power <: AbstractObjective
    val::Int
end
struct SampleSize <: AbstractObjective
    val::Float64
    function SampleSize(val)
        if val ≥ 1.0 || val ≤ 0.0 throw(ArgumentError("Beta ≥ 1.0 or ≤ 0.0!")) end
        new(val)::SampleSize
    end
end

struct Proportion <: AbstractSimpleProportion
    x::Int
    n::Int
end
struct Probability <: AbstractSimpleProportion
    p::Float64
    function Probability(p::Float64)
        if p ≥ 1.0 || p ≤ 0.0 throw(ArgumentError("Probability can't be ≥ 1.0 or ≤ 0.0!")) end
        new(p::Float64)::Probability
    end
end
function getval(p::Proportion)::Float64
    return p.x/p.n
end
function getval(p::Probability)::Float64
    return p.p
end
#=
struct Probability{T<: Float64, N <: Union{Int, Nothing}} <: AbstractParameter
    p::T
    n::N
    function Probability(p::T, n::N) where T <: Float64 where N <: Int
        new{T, N}(p, n)::Probability
    end
    function Probability(p::T, n::N = nothing) where T <: Float64 where N <: Nothing
        new{T, N}(p, nothing)::Probability
    end
end
=#
struct Mean{T <: Union{Int, Nothing}} <: AbstractMean
    m::Real
    sd::Real
    n::T
    function Mean(m::Real, sd::Real, n::T) where T <: Union{Int, Nothing}
        if sd ≤ 0.0 throw(ArgumentError("SD can't be ≤ 0.0!")) end
        new{T}(m, sd, n)::Mean
    end
    function Mean(m::Real, sd::Real)
        Mean(m, sd, nothing)
    end
end
"""
    Mean comparison hypothesis testing
    a - Test group value
    b - Reference group value {Mean} type, or reference value for one group test {Real}
"""
struct DiffMean{T <: Union{Mean, Real}} <: AbstractCompositeMean{T}
    a::Mean
    b::T
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

struct TaskResult{T <: AbstractTask}
    task::T
    method::Symbol
    result
end

struct TaskEstimate
    est
    pow
end

mutable struct CTask{T <: AbstractParameter, D <: AbstractDesign, H <: AbstractHypothesis, O <: AbstractObjective} <: AbstractTask
    param::T            #Testing parameter
    design::D           #Trial design
    hyp::H              #Hypothesis
    objective::O        #Objective (result)
    alpha::Real         #Alpha level
    k::Real             #Group coefficient
    function CTask(param::T, design::D, hyp::H, objective::O, alpha::Real, k::Real) where T <: AbstractParameter where D <: AbstractDesign where H <: AbstractHypothesis where O <: AbstractObjective
        #=
        if isa(hyp, Equality) && llim != ulim
            @warn "For Equality hypothesis llim and ulim can't be different! ulim set as llim!"
            ulim = llim
        end
        if k ≤ 0.0 throw(ArgumentError("Constant k can't be ≤ 0.0!")) end
        if alpha ≥ 1.0 || alpha ≤ 0.0 throw(ArgumentError("Alpha ≥ 1.0 or ≤ 0.0!")) end
        if isa(hyp, Equivalence) && diff ≤ 0 throw(ArgumentError("Diiference can't be ≤ 0.0 with Equivalence hypothesis!")) end
        =#
        new{T,D,H,O}(param, design, hyp, objective, alpha, k)::CTask{T,D,H,O}
    end
    function CTask(param::T, design::D, hyp::H, objective::O, alpha::Real) where T <: AbstractParameter where D <: AbstractDesign where H <: AbstractHypothesis where O <: AbstractObjective
        return CTask(param, design, hyp, objective, alpha, 1)
    end
end
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#main sample size function
function ctsamplen(;param=:notdef, type=:notdef, group=:notdef, alpha=0.05, beta=0.2, diff=0, sd=0, a=0, b=0, k=1, logscale=true)::TaskResult
    if alpha >= 1 || alpha <= 0 || beta >= 1 || beta <= 0  throw(ArgumentError("sampleSize: alpha and beta sould be > 0 and < 1.")) end
    if param == :prop && !(group == :one || group == :two || type == :mcnm)  throw(ArgumentError("sampleSize: group should be defined or mcnm type.")) end
    if sd ≤ 0 && param == :mean throw(ArgumentError("SD can't be ≤ 0.0!")) end
    if k <= 0 throw(ArgumentError("Constant k can't be ≤ 0.0!")) end
    if param == :mean
        if group == :one
            if type == :ea
                n = one_mean_equality(a, b, sd, alpha, beta)
                task = CTask(DiffMean(Mean(a, sd), b), Onegroup(), Equality(), SampleSize(beta), alpha)
            elseif type == :ei
                n = one_mean_equivalence(a, b, sd, diff, alpha, beta)
                task = CTask(DiffMean(Mean(a, sd), b),  Onegroup(), Equivalence(b-diff, b+diff), SampleSize(beta), alpha)
            elseif type == :ns
                n = one_mean_superiority(a, b, sd, diff, alpha, beta)
                task = CTask(DiffMean(Mean(a, sd), b), Onegroup(), Superiority(b + diff, diff), SampleSize(beta), alpha)
            else throw(ArgumentError("Keyword type unknown!")) end
        elseif group == :two
            if type == :ea
                n = two_mean_equality(a, b, sd, alpha, beta, k)
                task = CTask(DiffMean(Mean(a, sd), Mean(b, sd)), Parallel(), Equality(), SampleSize(beta), alpha, k)
            elseif type == :ei
                n = two_mean_equivalence(a, b, sd, diff, alpha, beta, k)
                task = CTask(DiffMean(Mean(a, sd), Mean(b, sd)), Parallel(), Equivalence(-diff, diff), SampleSize(beta), alpha, k)
            elseif type == :ns
                n = two_mean_superiority(a, b, sd, diff, alpha, beta, k)
                task = CTask(DiffMean(Mean(a, sd), Mean(b, sd)),  Parallel(), Superiority(diff, diff), SampleSize(beta), alpha, k)
            else throw(ArgumentError("Keyword type unknown!")) end
        else throw(ArgumentError("Keyword group unknown!")) end
    elseif param == :prop
        if 1 < a || a < 0 || 1 < b || b < 0 throw(ArgumentError("Keyword a or b out of the range!")) end
        if type == :mcnm
            n = mcnm(a, b, alpha, beta)
            task = CTask(DiffProportion(Probability(a), Probability(b)), Parallel(),  McNemars(), SampleSize(beta), alpha, k)
        else
            if group == :one
                if type == :ea
                    n = one_proportion_equality(a, b, alpha, beta)
                    task = CTask(DiffProportion(Probability(a), b), Onegroup(), Equality(), SampleSize(beta), alpha)
                elseif type == :ei
                    n = one_proportion_equivalence(a, b, diff, alpha, beta)
                    task = CTask(DiffProportion(Probability(a), b), Onegroup(), Equivalence(b-diff, b+diff), SampleSize(beta), alpha)
                elseif type == :ns
                    n = one_proportion_superiority(a, b, diff, alpha, beta)
                    task = CTask(DiffProportion(Probability(a), b), Onegroup(), Superiority(b+diff, diff), SampleSize(beta), alpha)
                else throw(ArgumentError("Keyword type unknown!")) end
            elseif group == :two
                if type == :ea
                    n = two_proportion_equality(a, b, alpha, beta, k)
                    task = CTask(DiffProportion(Probability(a), Probability(b)), Parallel(), Equality(), SampleSize(beta), alpha, k)
                elseif type == :ei
                    n = two_proportion_equivalence(a, b, diff, alpha, beta, k)
                    task = CTask(DiffProportion(Probability(a), Probability(b)), Parallel(), Equivalence(-diff, diff), SampleSize(beta), alpha, k)
                elseif type == :ns
                    n = two_proportion_superiority(a, b, diff, alpha, beta, k)
                    task = CTask(DiffProportion(Probability(a), Probability(b)), Parallel(), Superiority(diff, diff), SampleSize(beta), alpha, k)
                else throw(ArgumentError("Keyword type unknown!")) end
            else throw(ArgumentError("Keyword group unknown!")) end
        end
    elseif param == :or
        if type == :ea
            n = or_equality(a, b, alpha, beta, k)
            task = CTask(OddRatio(Probability(a), Probability(b)), Parallel(),  Equality(), SampleSize(beta), alpha, k)
        elseif type == :ei
            if !logscale diff = log(diff) end
            n = or_equivalence(a, b, diff, alpha, beta, k)
            task = CTask(OddRatio(Probability(a), Probability(b)), Parallel(), Equivalence(-diff, diff), SampleSize(beta), alpha, k)
        elseif type == :ns
            if !logscale diff = log(diff) end
            n = or_superiority(a, b, diff, alpha, beta, k)
            task = CTask(OddRatio(Probability(a), Probability(b)), Parallel(), Superiority(diff,  diff), SampleSize(beta), alpha, k)
        else throw(ArgumentError("Keyword type unknown!")) end
    else throw(ArgumentError("Keyword param unknown!")) end
    #TaskEstimate(task.design, n, k)
    return TaskResult(task, :chow, n)
end #sampleSize

function ctsamplen(t::CTask{T, D, H, O}) where T <: DiffMean{R} where R <: Real where D <: AbstractDesign where H <: Equality where O <: SampleSize
    return TaskResult(t, :chow, one_mean_equality(t.param.a.m, t.param.b, t.param.a.sd, t.alpha, t.objective.val))
end
function ctsamplen(t::CTask{T, D, H, O}) where T <: DiffMean{R} where R <: Real where D <: AbstractDesign where H <: Equivalence where O <: SampleSize
    return TaskResult(t, :chow, one_mean_equivalence(t.param.a.m, t.param.b, t.param.a.sd, mdiff(t.hyp), t.alpha, t.objective.val))
end
function ctsamplen(t::CTask{T, D, H, O}) where T <: DiffMean{R} where R <: Real where D <: AbstractDesign where H <: Superiority where O <: SampleSize
    return TaskResult(t, :chow, one_mean_superiority(t.param.a.m, t.param.b, t.param.a.sd, t.hyp.diff, t.alpha, t.objective.val))
end
#---
function ctsamplen(t::CTask{T, D, H, O}) where T <: DiffMean{M} where M <: AbstractMean where D <: AbstractDesign where H <: Equality where O <: SampleSize
    return TaskResult(t, :chow, two_mean_equality(t.param.a.m, t.param.b.m, t.param.a.sd, t.alpha, t.objective.val, t.k))
end
function ctsamplen(t::CTask{T, D, H, O}) where T <: DiffMean{M} where M <: AbstractMean  where D <: AbstractDesign where H <: Equivalence where O <: SampleSize
    return TaskResult(t, :chow, two_mean_equivalence(t.param.a.m, t.param.b.m, t.param.a.sd,  mdiff(t.hyp), t.alpha, t.objective.val, t.k))
end
function ctsamplen(t::CTask{T, D, H, O}) where T <: DiffMean{M} where M <: AbstractMean  where D <: AbstractDesign where H <: Superiority where O <: SampleSize
    return TaskResult(t, :chow, two_mean_superiority(t.param.a.m, t.param.b.m, t.param.a.sd, t.hyp.diff, t.alpha, t.objective.val, t.k))
end
#------
function ctsamplen(t::CTask{T, D, H, O}) where T <: DiffProportion{P, R} where P <: AbstractSimpleProportion where R <: Real  where D  <: AbstractDesign where H <: Equality where O <: SampleSize
    return TaskResult(t, :chow, one_proportion_equality(getval(t.param.a), t.param.b, t.alpha, t.objective.val))
end
function ctsamplen(t::CTask{T, D, H, O}) where T <: DiffProportion{P, R} where P <: AbstractSimpleProportion where R <: Real  where D <: AbstractDesign where H <: Equivalence where O <: SampleSize
    return TaskResult(t, :chow, one_proportion_equivalence(getval(t.param.a), t.param.b,  mdiff(t.hyp), t.alpha, t.objective.val))
end
function ctsamplen(t::CTask{T, D, H, O}) where T <: DiffProportion{P, R} where P <: AbstractSimpleProportion where R <: Real  where D <: AbstractDesign where H <: Superiority where O <: SampleSize
    return TaskResult(t, :chow, one_proportion_superiority(getval(t.param.a), t.param.b, t.hyp.diff, t.alpha, t.objective.val))
end
#---
function ctsamplen(t::CTask{T, D, H, O}) where T <: DiffProportion{P, P} where P <: AbstractSimpleProportion  where D <: AbstractDesign where H <: Equality where O <: SampleSize
    return TaskResult(t, :chow, two_proportion_equality(getval(t.param.a), getval(t.param.b), t.alpha, t.objective.val, t.k))
end
function ctsamplen(t::CTask{T, D, H, O}) where T <: DiffProportion{P, P} where P <: AbstractSimpleProportion  where D <: AbstractDesign where H <: Equivalence where O <: SampleSize
    return TaskResult(t, :chow, two_proportion_equivalence(getval(t.param.a), getval(t.param.b),  mdiff(t.hyp), t.alpha, t.objective.val, t.k))
end
function ctsamplen(t::CTask{T, D, H, O}) where T <: DiffProportion{P, P} where P <: AbstractSimpleProportion  where D <: AbstractDesign where H <: Superiority where O <: SampleSize
    return TaskResult(t, :chow, two_proportion_superiority(getval(t.param.a), getval(t.param.b), t.hyp.diff, t.alpha, t.objective.val, t.k))
end
#------
function ctsamplen(t::CTask{T, D, H, O}) where T <: OddRatio{P} where P <:  AbstractSimpleProportion  where D <: AbstractDesign where H <: Equality where O <: SampleSize
    return TaskResult(t, :chow, or_equality(getval(t.param.a), getval(t.param.b), t.alpha, t.objective.val, t.k))
end
function ctsamplen(t::CTask{T, D, H, O}) where T <: OddRatio{P} where P <:  AbstractSimpleProportion  where D <: AbstractDesign where H <: Equivalence where O <: SampleSize
    return TaskResult(t, :chow, or_equivalence(getval(t.param.a), getval(t.param.b),  mdiff(t.hyp), t.alpha, t.objective.val, t.k))
end
function ctsamplen(t::CTask{T, D, H, O}) where T <: OddRatio{P} where P <:  AbstractSimpleProportion  where D <: AbstractDesign where H <: Superiority where O <: SampleSize
    return TaskResult(t, :chow, or_superiority(getval(t.param.a), getval(t.param.b), t.hyp.diff, t.alpha, t.objective.val, t.k))
end
#------------
function ctsamplen(t::CTask{T, D, H, O}) where T <: DiffProportion{P, P} where P <: AbstractProportion  where D <: AbstractDesign where H <: McNemars where O <: SampleSize
    return TaskResult(t, :chow, mcnm(getval(t.param.a), getval(t.param.b), t.alpha, t.objective.val))
end

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#clinical trial power main function
function ctpower(;param=:notdef, type=:notdef, group=:notdef, alpha=0.05, logscale=true, diff=0, sd=0, a=0, b=0, n=0, k=1)::TaskResult
    if alpha >= 1 || alpha <= 0  throw(ArgumentError("Keyword alpha out of the range!")) end
    if param == :prop && !(group == :one || group == :two || type == :mcnm)  throw(ArgumentError("Keyword group or type defined incorrectly!")) end
    if sd == 0 && param == :mean throw(ArgumentError("Keyword sd = 0!")) end
    if k == 0 throw(ArgumentError("Keyword k = 0!")) end
    if n == 0 throw(ArgumentError("Keyword n = 0!")) end
    if param == :mean
        if group == :one
            if type == :ea
                pow =  one_mean_equality_pow(a, b, sd, alpha, n)
                task = CTask(DiffMean(Mean(a, sd), b), Parallel(), Equality(), Power(n), alpha)
            elseif type == :ei
                pow =  one_mean_equivalence_pow(a, b, sd, diff, alpha, n)
                task = CTask(DiffMean(Mean(a, sd), b), Parallel(),  Equivalence(b-diff, b+diff), Power(n), alpha)
            elseif type == :ns
                pow =  one_mean_superiority_pow(a, b, sd, diff, alpha, n)
                task = CTask(DiffMean(Mean(a, sd), b), Parallel(), Superiority(b+diff, diff), Power(n), alpha)
            else throw(ArgumentError("Keyword type unknown!")) end
        elseif group == :two
            if type == :ea
                pow =  two_mean_equality_pow(a, b, sd, alpha, n, k)
                task = CTask(DiffMean(Mean(a, sd), Mean(b, sd)), Parallel(), Equality(), Power(n), alpha, k)
            elseif type == :ei
                pow =  two_mean_equivalence_pow(a, b, sd, diff, alpha, n, k)
                task = CTask(DiffMean(Mean(a, sd), Mean(b, sd)), Parallel(), Equivalence(-diff, diff), Power(n), alpha, k)
            elseif type == :ns
                pow =  two_mean_superiority_pow(a, b, sd, diff, alpha, n, k)
                task = CTask(DiffMean(Mean(a, sd), Mean(b, sd)), Parallel(), Superiority(diff, diff), Power(n), alpha, k)
            else throw(ArgumentError("Keyword type unknown!")) end
        else throw(ArgumentError("Keyword group unknown!")) end
    elseif param == :prop
        if 1 < a || a < 0 || 1 < b || b < 0 throw(ArgumentError("Keyword a or b out of the range!")) end
        if type == :mcnm
            pow =  mcnm_pow(a, b, alpha, n)
            task = CTask(DiffProportion(Probability(a), Probability(b)), Parallel(), McNemars(), Power(n), alpha, k)
        else
            if group == :one
                if type == :ea
                    pow =  one_proportion_equality_pow(a, b, alpha, n)
                    task = CTask(DiffProportion(Probability(a), b), Parallel(), Equality(), Power(n), alpha)
                elseif type == :ei
                    pow =  one_proportion_equivalence_pow(a, b, diff, alpha, n)
                    task = CTask(DiffProportion(Probability(a), b), Parallel(), Equivalence(b-diff, b+diff), Power(n), alpha)
                elseif type == :ns
                    pow =  one_proportion_superiority_pow(a, b, diff, alpha, n)
                    task = CTask(DiffProportion(Probability(a), b), Parallel(), Superiority(b-diff, diff), Power(n), alpha)
                else throw(ArgumentError("Keyword type unknown!")) end
            elseif group == :two
                if type == :ea
                    pow =  two_proportion_equality_pow(a, b, alpha, n, k)
                    task = CTask(DiffProportion(Probability(a), Probability(b)), Parallel(), Equality(), Power(n), alpha, k)
                elseif type == :ei
                    pow =  two_proportion_equivalence_pow(a, b, diff, alpha, n, k)
                    task = CTask(DiffProportion(Probability(a), Probability(b)), Parallel(), Equivalence(-diff, +diff), Power(n), alpha, k)
                elseif type == :ns
                    pow =  two_proportion_superiority_pow(a, b, diff, alpha, n, k)
                    task = CTask(DiffProportion(Probability(a), Probability(b)), Parallel(), Superiority(diff, diff), Power(n), alpha, k)
                else throw(ArgumentError("Keyword type unknown!")) end
            else throw(ArgumentError("Keyword group unknown!")) end
        end
    elseif param == :or
        if type == :ea
            pow =  or_equality_pow(a, b, alpha, n, k)
            task = CTask(OddRatio(Probability(a), Probability(b)), Parallel(), Equality(), Power(n), alpha, k)
        elseif type == :ei
            if !logscale diff = log(diff) end
            pow =  or_equivalence_pow(a, b, diff, alpha, n, k)
            task = CTask(OddRatio(Probability(a), Probability(b)), Parallel(), Equivalence(-diff, +diff), Power(n), alpha, k)
        elseif type == :ns
            if !logscale diff = log(diff) end
            pow = or_superiority_pow(a, b, diff, alpha, n, k)
            task = CTask(OddRatio(Probability(a), Probability(b)), Parallel(), Superiority(diff, diff), Power(n), alpha, k)
        else throw(ArgumentError("Keyword type unknown!")) end
    else throw(ArgumentError("Keyword param unknown!")) end
    return TaskResult(task, :chow, pow)
end #ctpower
#---
function ctpower(t::CTask{T, D, H, O}) where T <: DiffMean{R} where R <: Real where D <: AbstractDesign where H <: Equality where O <: Power
    return TaskResult(t, :chow,  one_mean_equality_pow(t.param.a.m, t.param.b, t.param.a.sd, t.alpha, t.objective.val))
end
function ctpower(t::CTask{T, D, H, O}) where T <: DiffMean{R} where R <: Real where D <: AbstractDesign where H <: Equivalence where O <: Power
    return TaskResult(t, :chow,  one_mean_equivalence_pow(t.param.a.m, t.param.b, t.param.a.sd, mdiff(t.hyp), t.alpha, t.objective.val))
end
function ctpower(t::CTask{T, D, H, O}) where T <: DiffMean{R} where R <: Real where D <: AbstractDesign where H <: Superiority where O <: Power
    return TaskResult(t, :chow,  one_mean_superiority_pow(t.param.a.m, t.param.b, t.param.a.sd, t.hyp.diff, t.alpha, t.objective.val))
end
#------
function ctpower(t::CTask{T, D, H, O}) where T <: DiffMean{M} where M <: AbstractMean where D <: AbstractDesign where H <: Equality where O <: Power
    return TaskResult(t, :chow, two_mean_equality_pow(t.param.a.m, t.param.b.m, t.param.a.sd, t.alpha, t.objective.val, t.k))
end
function ctpower(t::CTask{T, D, H, O}) where T <: DiffMean{M} where M <: AbstractMean where D <: AbstractDesign where H <: Equivalence where O <: Power
    return TaskResult(t, :chow, two_mean_equivalence_pow(t.param.a.m, t.param.b.m, t.param.a.sd, mdiff(t.hyp), t.alpha, t.objective.val, t.k))
end
function ctpower(t::CTask{T, D, H, O}) where T <: DiffMean{M} where M <: AbstractMean where D <: AbstractDesign where H <: Superiority where O <: Power
    return TaskResult(t, :chow, two_mean_superiority_pow(t.param.a.m, t.param.b.m, t.param.a.sd, t.hyp.diff, t.alpha, t.objective.val, t.k))
end
#---
function ctpower(t::CTask{T, D, H, O}) where T <: DiffProportion{P, R} where P <: AbstractSimpleProportion where R <: Real where D <: AbstractDesign where H <: Equality where O <: Power
    return TaskResult(t, :chow, one_proportion_equality_pow(getval(t.param.a), t.param.b, t.alpha, t.objective.val))
end
function ctpower(t::CTask{T, D, H, O}) where T <: DiffProportion{P, R} where P <: AbstractSimpleProportion where R <: Real where D <: AbstractDesign where H <: Equivalence where O <: Power
    return TaskResult(t, :chow, one_proportion_equivalence_pow(getval(t.param.a), t.param.b, mdiff(t.hyp), t.alpha, t.objective.val))
end
function ctpower(t::CTask{T, D, H, O}) where T <: DiffProportion{P, R} where P <: AbstractSimpleProportion where R <: Real where D <: AbstractDesign where H <: Superiority where O <: Power
    return TaskResult(t, :chow, one_proportion_superiority_pow(getval(t.param.a), t.param.b, t.hyp.diff, t.alpha, t.objective.val))
end
#---
function ctpower(t::CTask{T, D, H, O}) where T <: DiffProportion{P, P} where P <: AbstractSimpleProportion where D <: AbstractDesign where H <: Equality where O <: Power
    return TaskResult(t, :chow, two_proportion_equality_pow(getval(t.param.a), getval(t.param.b), t.alpha, t.objective.val, t.k))
end
function ctpower(t::CTask{T, D, H, O}) where T <: DiffProportion{P, P} where P <: AbstractSimpleProportion where D <: AbstractDesign where H <: Equivalence where O <: Power
    return TaskResult(t, :chow, two_proportion_equivalence_pow(getval(t.param.a), getval(t.param.b), mdiff(t.hyp), t.alpha, t.objective.val, t.k))
end
function ctpower(t::CTask{T, D, H, O}) where T <: DiffProportion{P, P} where P <: AbstractSimpleProportion where D <: AbstractDesign where H <: Superiority where O <: Power
    return TaskResult(t, :chow, two_proportion_superiority_pow(getval(t.param.a), getval(t.param.b), t.hyp.diff, t.alpha, t.objective.val, t.k))
end
#------
function ctpower(t::CTask{T, D, H, O}) where T <: OddRatio{P} where P <:  AbstractSimpleProportion where D <: AbstractDesign where H <: Equality where O <: Power
    return TaskResult(t, :chow, or_equality_pow(getval(t.param.a), getval(t.param.b), t.alpha, t.objective.val, t.k))
end
function ctpower(t::CTask{T, D, H, O}) where T <: OddRatio{P} where P <:  AbstractSimpleProportion where D <: AbstractDesign where H <: Equivalence where O <: Power
    return TaskResult(t, :chow, or_equivalence_pow(getval(t.param.a), getval(t.param.b), mdiff(t.hyp), t.alpha, t.objective.val, t.k))
end
function ctpower(t::CTask{T, D, H, O}) where T <: OddRatio{P} where P <:  AbstractSimpleProportion where D <: AbstractDesign where H <: Superiority where O <: Power
    return TaskResult(t, :chow, or_superiority_pow(getval(t.param.a), getval(t.param.b), t.hyp.diff, t.alpha, t.objective.val, t.k))
end
#------------
function ctpower(t::CTask{T, D, H, O}) where T <: DiffProportion{P, P} where P <: AbstractProportion where D <: AbstractDesign  where H <: McNemars where O <: Power
    return TaskResult(t, :chow, mcnm_pow(getval(t.param.a), getval(t.param.b), t.alpha, t.objective.val))
end
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
function besamplen(;alpha::Real=0.05, beta::Real=0.2, theta0::Real=0.95, theta1::Real=0.8, theta2::Real=1.25, cv::Real=0.0, sd::Real=0.0, logscale::Bool=true, design::Symbol=:d2x2, method::Symbol=:owenq)::TaskResult
    if !(0 < beta < 1) throw(ArgumentError("Beta ≥ 1.0 or ≤ 0.0!")) end
    if !(0 < alpha < 1) throw(ArgumentError("Alfa ≥ 1.0 or ≤ 0.0!")) end
    if theta1 ≥ theta2  throw(ArgumentError("theta1 ≥ theta2!")) end
    if !logscale && cv > 0.0 && sd ≤ 0.0 throw(ArgumentError("Use sd instead cv for non-logscale parameters!")) end
    if !logscale && sd ≤ 0.0 throw(ArgumentError("sd ≤ 0.0!")) end
    if !logscale && cv > 0.0 @warn "sd and cv provided for non-logscale parameters, only sd used!" end
    if logscale && cv ≤ 0.0 && sd > 0.0 throw(ArgumentError("Use cv instead sd for non-logscale parameters!")) end
    if logscale && cv ≤ 0.0 throw(ArgumentError("cv ≤ 0"))  end
    if logscale && sd > 0.0 @warn "sd and cv provided for logscale parameters, only cv used!" end

    if logscale
        if theta1 ≤ 0 || theta2 ≤ 0 || theta0 ≤ 0 throw(ArgumentError("theta0 or theta1 or theta2 ≤ 0!")) end
        theta1 = log(theta1)
        theta2 = log(theta2)
        theta0 = log(theta0)
        sd     = cv2sd(cv)
    end

    task = CTask(DiffMean(Mean(0, sd), Mean(theta0, sd)), Crossover(design), Equivalence(theta1, theta2, bio=true), SampleSize(beta), alpha)

    estn, estpow = samplentostint(alpha, theta1, theta2, theta0, sd, beta, design, method)
    return TaskResult(task, :chow, estn)
end

function besamplen(t::CTask)

end

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
function bepower(;alpha::Real=0.05, theta1::Real=0.8, theta2::Real=1.25, theta0::Real=0.95, cv::Real=0.0, sd::Real=0.0, n::Int=0, logscale::Bool=true, design::Symbol=:d2x2, method::Symbol=:owenq)::TaskResult
    if n < 2 throw(ArgumentError("n < 2")) end

    if !(0 < alpha < 1) throw(ArgumentError("Alfa ≥ 1.0 or ≤ 0.0")) end
    if !logscale && cv > 0.0 && sd ≤ 0.0 throw(ArgumentError("Use sd instead cv for non-logscale parameters!")) end
    if !logscale && sd ≤ 0.0 throw(ArgumentError("sd ≤ 0.0!")) end
    if !logscale && cv > 0.0 @warn "sd and cv provided for non-logscale parameters, only sd used!" end
    if logscale && cv ≤ 0.0 && sd > 0.0 throw(ArgumentError("Use cv instead sd for non-logscale parameters!")) end
    if logscale && cv ≤ 0.0 throw(ArgumentError("cv ≤ 0"))  end
    if logscale && sd > 0.0 @warn "sd and cv provided for logscale parameters, only cv used!" end

    if logscale
        theta1 = log(theta1)
        theta2 = log(theta2)
        theta0 = log(theta0)
        sd     = cv2sd(cv)    # sqrt(ms)
    end

    task = CTask(DiffMean(Mean(theta0, sd), Mean(0, sd)), Crossover(design), Equivalence(theta1, theta2, bio=true), Power(n), alpha)

    pow =  powertostint(alpha,  theta1, theta2, theta0, sd, n, design, method)

    return TaskResult(task, method, pow)
end #bepower

function bepower(t::CTask; method::Symbol = :owenq)
    df    = t.design.df(t.objective.val)
    σ̵ₓ    = t.param.a.sd*sediv(t.design, t.objective.val)
    α     = t.alpha
    θ₁    = t.hyp.llim
    θ₂    = t.hyp.ulim
    δ     = t.param.a.m - t.param.b.m
    if method     == :owenq
        pow =    powertost_owenq(α, θ₁, θ₂, δ, σ̵ₓ, df)
    elseif method == :nct
        pow =      powertost_nct(α, θ₁, θ₂, δ, σ̵ₓ, df)
    elseif method == :mvt
        pow =      powertost_mvt(α, θ₁, θ₂, δ, σ̵ₓ, df) #not implemented
    elseif method == :shifted
        pow =  powertost_shifted(α, θ₁, θ₂, δ, σ̵ₓ, df)
    else
         throw(ArgumentError("method not known!"))
    end
    return TaskResult(t, method, pow)
end
