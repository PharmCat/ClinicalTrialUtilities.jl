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

struct Equivalence <: AbstractHypothesis
    bio::Bool
    function Equivalence(;bio::Bool = false)
        new(bio)::Equivalence
    end
end
struct Equality <: AbstractHypothesis end
struct Superiority <: AbstractHypothesis end
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

struct TOSTSampleSizeTask <: AbstractTask
    alpha::Real
    beta::Real
    gmr::Real
    llim::Real
    ulim::Real
    cv::Real
    logscale::Bool
    design::Symbol
    method::Symbol
end
struct TOSTPowerTask <: AbstractTask
    alpha::Real
    n::Real
    gmr::Real
    llim::Real
    ulim::Real
    cv::Real
    logscale::Bool
    design::Symbol
    method::Symbol
end
struct TaskResult{T <: AbstractTask}
    task::T
    method::Symbol
    result::Real
end

struct CTask{T <: AbstractParameter, H <: AbstractHypothesis, O <: AbstractObjective} <: AbstractTask
    param::T            #Testing parameter
    llim::Real          #Lower lmit for Test group
    ulim::Real          #Upper lmit for Test group
    diff::Real          #Margin difference
    alpha::Real         #Alpha level
    hyp::H              #Hypothesis
    k::Real             #Group coefficient
    objective::O        #Objective (result)
    function CTask(param::T, llim::Real, ulim::Real, diff::Real, alpha::Real, hyp::H, k::Real, objective::O) where T <: AbstractParameter where H <: AbstractHypothesis where O <: AbstractObjective
        if isa(hyp, Equality) && llim != ulim
            @warn "For Equality hypothesis llim and ulim can't be different! ulim set as llim!"
            ulim = llim
        end
        if k ≤ 0.0 throw(ArgumentError("Constant k can't be ≤ 0.0!")) end
        if alpha ≥ 1.0 || alpha ≤ 0.0 throw(ArgumentError("Alpha ≥ 1.0 or ≤ 0.0!")) end
        if isa(hyp, Equivalence) && diff ≤ 0 throw(ArgumentError("Diiference can't be ≤ 0.0 with Equivalence hypothesis!")) end
        new{T,H,O}(param, llim, ulim, diff, alpha, hyp, k, objective)::CTask{T,H,O}
    end
    function CTask(param::T, llim::Real, ulim::Real, alpha::Real, hyp::H, k::Real, objective::O) where T <: AbstractParameter where H <: AbstractHypothesis where O <: AbstractObjective
        return CTask(param, llim, ulim, NaN, alpha, hyp, k, objective)
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
                task = CTask(DiffMean(Mean(a, sd), b), b, b, alpha, Equality(), k, SampleSize(beta))
            elseif type == :ei
                n = one_mean_equivalence(a, b, sd, diff, alpha, beta)
                task = CTask(DiffMean(Mean(a, sd), b), b-diff, b+diff, diff, alpha, Equivalence(), k, SampleSize(beta))
            elseif type == :ns
                n = one_mean_superiority(a, b, sd, diff, alpha, beta)
                task = CTask(DiffMean(Mean(a, sd), b), b + diff, Inf, diff, alpha, Superiority(), k, SampleSize(beta))
            else throw(ArgumentError("Keyword type unknown!")) end
        elseif group == :two
            if type == :ea
                n = two_mean_equality(a, b, sd, alpha, beta, k)
                task = CTask(DiffMean(Mean(a, sd), Mean(b, sd)), 0, 0, alpha, Equality(), k, SampleSize(beta))
            elseif type == :ei
                n = two_mean_equivalence(a, b, sd, diff, alpha, beta, k)
                task = CTask(DiffMean(Mean(a, sd), Mean(b, sd)), -diff, diff, diff, alpha, Equivalence(), k, SampleSize(beta))
            elseif type == :ns
                n = two_mean_superiority(a, b, sd, diff, alpha, beta, k)
                task = CTask(DiffMean(Mean(a, sd), Mean(b, sd)), diff, Inf, diff, alpha, Superiority(), k, SampleSize(beta))
            else throw(ArgumentError("Keyword type unknown!")) end
        else throw(ArgumentError("Keyword group unknown!")) end
    elseif param == :prop
        if 1 < a || a < 0 || 1 < b || b < 0 throw(ArgumentError("Keyword a or b out of the range!")) end
        if type == :mcnm
            n = mcnm(a, b, alpha, beta)
            task = CTask(DiffProportion(Probability(a), Probability(b)), 0, 0, alpha, McNemars(), k, SampleSize(beta))
        else
            if group == :one
                if type == :ea
                    n = one_proportion_equality(a, b, alpha, beta)
                    task = CTask(DiffProportion(Probability(a), b), b, b, alpha, Equality(), 1, SampleSize(beta))
                elseif type == :ei
                    n = one_proportion_equivalence(a, b, diff, alpha, beta)
                    task = CTask(DiffProportion(Probability(a), b), b-diff, b+diff, diff, alpha, Equivalence(), 1, SampleSize(beta))
                elseif type == :ns
                    n = one_proportion_superiority(a, b, diff, alpha, beta)
                    task = CTask(DiffProportion(Probability(a), b), b+diff, Inf, diff, alpha, Superiority(), 1, SampleSize(beta))
                else throw(ArgumentError("Keyword type unknown!")) end
            elseif group == :two
                if type == :ea
                    n = two_proportion_equality(a, b, alpha, beta, k)
                    task = CTask(DiffProportion(Probability(a), Probability(b)), 0, 0, alpha, Equality(), k, SampleSize(beta))
                elseif type == :ei
                    n = two_proportion_equivalence(a, b, diff, alpha, beta, k)
                    task = CTask(DiffProportion(Probability(a), Probability(b)), -diff, diff, diff, alpha, Equivalence(), k, SampleSize(beta))
                elseif type == :ns
                    n = two_proportion_superiority(a, b, diff, alpha, beta, k)
                    task = CTask(DiffProportion(Probability(a), Probability(b)), diff, Inf, diff, alpha, Superiority(), k, SampleSize(beta))
                else throw(ArgumentError("Keyword type unknown!")) end
            else throw(ArgumentError("Keyword group unknown!")) end
        end
    elseif param == :or
        if type == :ea
            n = or_equality(a, b, alpha, beta, k)
            task = CTask(OddRatio(Probability(a), Probability(b)), 0, 0, alpha, Equality(), k, SampleSize(beta))
        elseif type == :ei
            if !logscale diff = log(diff) end
            n = or_equivalence(a, b, diff, alpha, beta, k)
            task = CTask(OddRatio(Probability(a), Probability(b)), -diff, diff, diff, alpha, Equivalence(), k, SampleSize(beta))
        elseif type == :ns
            if !logscale diff = log(diff) end
            n = or_superiority(a, b, diff, alpha, beta, k)
            task = CTask(OddRatio(Probability(a), Probability(b)), diff, Inf, diff, alpha, Superiority(), k, SampleSize(beta))
        else throw(ArgumentError("Keyword type unknown!")) end
    else throw(ArgumentError("Keyword param unknown!")) end
    return TaskResult(task, :chow, n)
end #sampleSize

function ctsamplen(t::CTask{T, H, O}) where T <: DiffMean{R} where R <: Real where H <: Equality where O <: SampleSize
    return TaskResult(t, :chow, one_mean_equality(t.param.a.m, t.param.b, t.param.a.sd, t.alpha, t.objective.val))
end
function ctsamplen(t::CTask{T, H, O}) where T <: DiffMean{R} where R <: Real where H <: Equivalence where O <: SampleSize
    return TaskResult(t, :chow, one_mean_equivalence(t.param.a.m, t.param.b, t.param.a.sd, t.diff, t.alpha, t.objective.val))
end
function ctsamplen(t::CTask{T, H, O}) where T <: DiffMean{R} where R <: Real where H <: Superiority where O <: SampleSize
    return TaskResult(t, :chow, one_mean_superiority(t.param.a.m, t.param.b, t.param.a.sd, t.diff, t.alpha, t.objective.val))
end
#---
function ctsamplen(t::CTask{T, H, O}) where T <: DiffMean{M} where M <: AbstractMean where H <: Equality where O <: SampleSize
    return TaskResult(t, :chow, two_mean_equality(t.param.a.m, t.param.b.m, t.param.a.sd, t.alpha, t.objective.val, t.k))
end
function ctsamplen(t::CTask{T, H, O}) where T <: DiffMean{M} where M <: AbstractMean where H <: Equivalence where O <: SampleSize
    return TaskResult(t, :chow, two_mean_equivalence(t.param.a.m, t.param.b.m, t.param.a.sd, t.diff, t.alpha, t.objective.val, t.k))
end
function ctsamplen(t::CTask{T, H, O}) where T <: DiffMean{M} where M <: AbstractMean where H <: Superiority where O <: SampleSize
    return TaskResult(t, :chow, two_mean_superiority(t.param.a.m, t.param.b.m, t.param.a.sd, t.diff, t.alpha, t.objective.val, t.k))
end
#------
function ctsamplen(t::CTask{T, H, O}) where T <: DiffProportion{P, R} where P <: AbstractSimpleProportion where R <: Real where H <: Equality where O <: SampleSize
    return TaskResult(t, :chow, one_proportion_equality(getval(t.param.a), t.param.b, t.alpha, t.objective.val))
end
function ctsamplen(t::CTask{T, H, O}) where T <: DiffProportion{P, R} where P <: AbstractSimpleProportion where R <: Real where H <: Equivalence where O <: SampleSize
    return TaskResult(t, :chow, one_proportion_equivalence(getval(t.param.a), t.param.b, t.diff, t.alpha, t.objective.val))
end
function ctsamplen(t::CTask{T, H, O}) where T <: DiffProportion{P, R} where P <: AbstractSimpleProportion where R <: Real where H <: Superiority where O <: SampleSize
    return TaskResult(t, :chow, one_proportion_superiority(getval(t.param.a), t.param.b, t.diff, t.alpha, t.objective.val))
end
#---
function ctsamplen(t::CTask{T, H, O}) where T <: DiffProportion{P, P} where P <: AbstractSimpleProportion where H <: Equality where O <: SampleSize
    return TaskResult(t, :chow, two_proportion_equality(getval(t.param.a), getval(t.param.b), t.alpha, t.objective.val, t.k))
end
function ctsamplen(t::CTask{T, H, O}) where T <: DiffProportion{P, P} where P <: AbstractSimpleProportion where H <: Equivalence where O <: SampleSize
    return TaskResult(t, :chow, two_proportion_equivalence(getval(t.param.a), getval(t.param.b), t.diff, t.alpha, t.objective.val, t.k))
end
function ctsamplen(t::CTask{T, H, O}) where T <: DiffProportion{P, P} where P <: AbstractSimpleProportion where H <: Superiority where O <: SampleSize
    return TaskResult(t, :chow, two_proportion_superiority(getval(t.param.a), getval(t.param.b), t.diff, t.alpha, t.objective.val, t.k))
end
#------
function ctsamplen(t::CTask{T, H, O}) where T <: OddRatio{P} where P <:  AbstractSimpleProportion where H <: Equality where O <: SampleSize
    return TaskResult(t, :chow, or_equality(getval(t.param.a), getval(t.param.b), t.alpha, t.objective.val, t.k))
end
function ctsamplen(t::CTask{T, H, O}) where T <: OddRatio{P} where P <:  AbstractSimpleProportion where H <: Equivalence where O <: SampleSize
    return TaskResult(t, :chow, or_equivalence(getval(t.param.a), getval(t.param.b), t.diff, t.alpha, t.objective.val, t.k))
end
function ctsamplen(t::CTask{T, H, O}) where T <: OddRatio{P} where P <:  AbstractSimpleProportion where H <: Superiority where O <: SampleSize
    return TaskResult(t, :chow, or_superiority(getval(t.param.a), getval(t.param.b), t.diff, t.alpha, t.objective.val, t.k))
end
#------------
function ctsamplen(t::CTask{T, H, O}) where T <: DiffProportion{P, P} where P <: AbstractProportion where H <: McNemars where O <: SampleSize
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
                task = CTask(DiffMean(Mean(a, sd), b), b, b, alpha, Equality(), 1, Power(n))
            elseif type == :ei
                pow =  one_mean_equivalence_pow(a, b, sd, diff, alpha, n)
                task = CTask(DiffMean(Mean(a, sd), b), b-diff, b+diff, diff, alpha, Equivalence(), 1, Power(n))
            elseif type == :ns
                pow =  one_mean_superiority_pow(a, b, sd, diff, alpha, n)
                task = CTask(DiffMean(Mean(a, sd), b), b+diff, Inf, diff, alpha, Superiority(), 1, Power(n))
            else throw(ArgumentError("Keyword type unknown!")) end
        elseif group == :two
            if type == :ea
                pow =  two_mean_equality_pow(a, b, sd, alpha, n, k)
                task = CTask(DiffMean(Mean(a, sd), Mean(b, sd)), 0, 0, alpha, Equality(), k, Power(n))
            elseif type == :ei
                pow =  two_mean_equivalence_pow(a, b, sd, diff, alpha, n, k)
                task = CTask(DiffMean(Mean(a, sd), Mean(b, sd)), -diff, diff, diff, alpha, Equivalence(), k, Power(n))
            elseif type == :ns
                pow =  two_mean_superiority_pow(a, b, sd, diff, alpha, n, k)
                task = CTask(DiffMean(Mean(a, sd), Mean(b, sd)), diff, Inf, diff, alpha, Superiority(), k, Power(n))
            else throw(ArgumentError("Keyword type unknown!")) end
        else throw(ArgumentError("Keyword group unknown!")) end
    elseif param == :prop
        if 1 < a || a < 0 || 1 < b || b < 0 throw(ArgumentError("Keyword a or b out of the range!")) end
        if type == :mcnm
            pow =  mcnm_pow(a, b, alpha, n)
            task = CTask(DiffProportion(Probability(a), Probability(b)), 0, 0, alpha, McNemars(), k, Power(n))
        else
            if group == :one
                if type == :ea
                    pow =  one_proportion_equality_pow(a, b, alpha, n)
                    task = CTask(DiffProportion(Probability(a), b), b, b, alpha, Equality(), 1, Power(n))
                elseif type == :ei
                    pow =  one_proportion_equivalence_pow(a, b, diff, alpha, n)
                    task = CTask(DiffProportion(Probability(a), b), b-diff, b+diff, diff, alpha, Equivalence(), 1, Power(n))
                elseif type == :ns
                    pow =  one_proportion_superiority_pow(a, b, diff, alpha, n)
                    task = CTask(DiffProportion(Probability(a), b), b-diff, Inf, diff, alpha, Superiority(), 1, Power(n))
                else throw(ArgumentError("Keyword type unknown!")) end
            elseif group == :two
                if type == :ea
                    pow =  two_proportion_equality_pow(a, b, alpha, n, k)
                    task = CTask(DiffProportion(Probability(a), Probability(b)), 0, 0, alpha, Equality(), k, Power(n))
                elseif type == :ei
                    pow =  two_proportion_equivalence_pow(a, b, diff, alpha, n, k)
                    task = CTask(DiffProportion(Probability(a), Probability(b)), -diff, +diff, diff, alpha, Equivalence(), k, Power(n))
                elseif type == :ns
                    pow =  two_proportion_superiority_pow(a, b, diff, alpha, n, k)
                    task = CTask(DiffProportion(Probability(a), Probability(b)), diff, Inf, diff, alpha, Superiority(), k, Power(n))
                else throw(ArgumentError("Keyword type unknown!")) end
            else throw(ArgumentError("Keyword group unknown!")) end
        end
    elseif param == :or
        if type == :ea
            pow =  or_equality_pow(a, b, alpha, n, k)
            task = CTask(OddRatio(Probability(a), Probability(b)), 0, 0, alpha, Equality(), k, Power(n))
        elseif type == :ei
            if !logscale diff = log(diff) end
            pow =  or_equivalence_pow(a, b, diff, alpha, n, k)
            task = CTask(OddRatio(Probability(a), Probability(b)), -diff, +diff, diff, alpha, Equivalence(), k, Power(n))
        elseif type == :ns
            if !logscale diff = log(diff) end
            pow = or_superiority_pow(a, b, diff, alpha, n, k)
            task = CTask(OddRatio(Probability(a), Probability(b)), diff, Inf, diff, alpha, Superiority(), k, Power(n))
        else throw(ArgumentError("Keyword type unknown!")) end
    else throw(ArgumentError("Keyword param unknown!")) end
    return TaskResult(task, :chow, pow)
end #ctpower
#---
function ctpower(t::CTask{T, H, O}) where T <: DiffMean{R} where R <: Real where H <: Equality where O <: Power
    return TaskResult(t, :chow,  one_mean_equality_pow(t.param.a.m, t.param.b, t.param.a.sd, t.alpha, t.objective.val))
end
function ctpower(t::CTask{T, H, O}) where T <: DiffMean{R} where R <: Real where H <: Equivalence where O <: Power
    return TaskResult(t, :chow,  one_mean_equivalence_pow(t.param.a.m, t.param.b, t.param.a.sd, t.diff, t.alpha, t.objective.val))
end
function ctpower(t::CTask{T, H, O}) where T <: DiffMean{R} where R <: Real where H <: Superiority where O <: Power
    return TaskResult(t, :chow,  one_mean_superiority_pow(t.param.a.m, t.param.b, t.param.a.sd, t.diff, t.alpha, t.objective.val))
end
#------
function ctpower(t::CTask{T, H, O}) where T <: DiffMean{M} where M <: AbstractMean where H <: Equality where O <: Power
    return TaskResult(t, :chow, two_mean_equality_pow(t.param.a.m, t.param.b.m, t.param.a.sd, t.alpha, t.objective.val, t.k))
end
function ctpower(t::CTask{T, H, O}) where T <: DiffMean{M} where M <: AbstractMean where H <: Equivalence where O <: Power
    return TaskResult(t, :chow, two_mean_equivalence_pow(t.param.a.m, t.param.b.m, t.param.a.sd, t.diff, t.alpha, t.objective.val, t.k))
end
function ctpower(t::CTask{T, H, O}) where T <: DiffMean{M} where M <: AbstractMean where H <: Superiority where O <: Power
    return TaskResult(t, :chow, two_mean_superiority_pow(t.param.a.m, t.param.b.m, t.param.a.sd, t.diff, t.alpha, t.objective.val, t.k))
end
#---
function ctpower(t::CTask{T, H, O}) where T <: DiffProportion{P, R} where P <: AbstractSimpleProportion where R <: Real where H <: Equality where O <: Power
    return TaskResult(t, :chow, one_proportion_equality_pow(getval(t.param.a), t.param.b, t.alpha, t.objective.val))
end
function ctpower(t::CTask{T, H, O}) where T <: DiffProportion{P, R} where P <: AbstractSimpleProportion where R <: Real where H <: Equivalence where O <: Power
    return TaskResult(t, :chow, one_proportion_equivalence_pow(getval(t.param.a), t.param.b, t.diff, t.alpha, t.objective.val))
end
function ctpower(t::CTask{T, H, O}) where T <: DiffProportion{P, R} where P <: AbstractSimpleProportion where R <: Real where H <: Superiority where O <: Power
    return TaskResult(t, :chow, one_proportion_superiority_pow(getval(t.param.a), t.param.b, t.diff, t.alpha, t.objective.val))
end
#---
function ctpower(t::CTask{T, H, O}) where T <: DiffProportion{P, P} where P <: AbstractSimpleProportion where H <: Equality where O <: Power
    return TaskResult(t, :chow, two_proportion_equality_pow(getval(t.param.a), getval(t.param.b), t.alpha, t.objective.val, t.k))
end
function ctpower(t::CTask{T, H, O}) where T <: DiffProportion{P, P} where P <: AbstractSimpleProportion where H <: Equivalence where O <: Power
    return TaskResult(t, :chow, two_proportion_equivalence_pow(getval(t.param.a), getval(t.param.b), t.diff, t.alpha, t.objective.val, t.k))
end
function ctpower(t::CTask{T, H, O}) where T <: DiffProportion{P, P} where P <: AbstractSimpleProportion where H <: Superiority where O <: Power
    return TaskResult(t, :chow, two_proportion_superiority_pow(getval(t.param.a), getval(t.param.b), t.diff, t.alpha, t.objective.val, t.k))
end
#------
function ctpower(t::CTask{T, H, O}) where T <: OddRatio{P} where P <:  AbstractSimpleProportion where H <: Equality where O <: Power
    return TaskResult(t, :chow, or_equality_pow(getval(t.param.a), getval(t.param.b), t.alpha, t.objective.val, t.k))
end
function ctpower(t::CTask{T, H, O}) where T <: OddRatio{P} where P <:  AbstractSimpleProportion where H <: Equivalence where O <: Power
    return TaskResult(t, :chow, or_equivalence_pow(getval(t.param.a), getval(t.param.b), t.diff, t.alpha, t.objective.val, t.k))
end
function ctpower(t::CTask{T, H, O}) where T <: OddRatio{P} where P <:  AbstractSimpleProportion where H <: Superiority where O <: Power
    return TaskResult(t, :chow, or_superiority_pow(getval(t.param.a), getval(t.param.b), t.diff, t.alpha, t.objective.val, t.k))
end
#------------
function ctpower(t::CTask{T, H, O}) where T <: DiffProportion{P, P} where P <: AbstractProportion where H <: McNemars where O <: Power
    return TaskResult(t, :chow, mcnm_pow(getval(t.param.a), getval(t.param.b), t.alpha, t.objective.val))
end
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
function besamplen(;alpha=0.05, beta=0.2, theta0=0.95, theta1=0.8, theta2=1.25, cv=0.0, logscale=true, design=:d2x2, method=:owenq)::TaskResult{TOSTSampleSizeTask}

    theta0 = convert(Float64, theta0); theta1 = convert(Float64, theta1); theta2 = convert(Float64, theta2); cv = convert(Float64, cv); alpha  = convert(Float64, alpha); beta = convert(Float64, beta)

    if cv <= 0 throw(ArgumentError("besamplen: cv can not be <= 0")) end
    if theta1 >= theta2  throw(ArgumentError("besamplen: theta1 should be < theta2")) end
    if alpha >= 1 || alpha <= 0 || beta >= 1 || beta <= 0 throw(ArgumentError("besamplen: alpha and beta shold be > 0 and < 1")) end

    task = TOSTSampleSizeTask(alpha, beta, theta0, theta1, theta2, cv, logscale, design, method)

    if logscale
        if theta1 < 0 || theta2 < 0 || theta0 < 0 throw(ArgumentError("besamplen: theta0, theta1, theta2 shold be > 0 and ")) end
        ltheta1 = log(theta1)
        ltheta2 = log(theta2)
        diffm   = log(theta0)
        sd      = cv2sd(cv)
    else
        ltheta1 = theta1
        ltheta2 = theta2
        diffm   = theta0
        sd      = cv
    end
    #values for approximate n
    td = (ltheta1 + ltheta2)/2
    rd = abs((ltheta1 - ltheta2)/2)

    #if rd <= 0 return false end
    d0 = diffm - td
    #approximate n
    n0::Int = convert(Int, ceil(two_mean_equivalence(0, d0, sd, rd, alpha, beta, 1)/2)*2)
    tp = 1 - beta  #target power
    if n0 < 4 n0 = 4 end
    if n0 > 5000 n0 = 5000 end
    pow = powertostint(alpha,  ltheta1, ltheta2, diffm, sd, n0, design, method)
    np::Int = 2
    powp::Float64 = pow
    if pow > tp
        while (pow > tp)
            np = n0
            powp = pow
            n0 = n0 - 2
            #pow = powerTOST(;alpha=alpha, logscale=false, theta1=ltheta1, theta2=ltheta2, theta0=diffm, cv=se, n=n0, design=design, method=method)
            if n0 < 4 break end #n0, pow end
            pow = powertostint(alpha,  ltheta1, ltheta2, diffm, sd, n0, design, method)
        end
        estpower = powp
        estn     = np
    elseif pow < tp
        while (pow < tp)
            np = n0
            powp = pow
            n0 = n0 + 2
            #pow = powerTOST(;alpha=alpha, logscale=false, theta1=ltheta1, theta2=ltheta2, theta0=diffm, cv=se, n=n0, design=design, method=method)
            pow = powertostint(alpha,  ltheta1, ltheta2, diffm, sd, n0, design, method)
            if n0 > 10000  break end # n0, pow end
        end
        estpower = pow
        estn     = n0
    else
        estpower = pow
        estn     = n0
    end
    return TaskResult(task, :chow, estn)
end
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
function bepower(;alpha=0.05, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.0, n=0, logscale=true, design=:d2x2, method=:owenq,  out=:num)::Float64
    if n < 2 throw(ArgumentError("powerTOST: n can't be < 2")) end
    if cv == 0.0 throw(ArgumentError("powerTOST: cv can't be equal to 0"))  end
    if !(0 < alpha < 1) throw(ArgumentError("powerTOST: alfa can't be > 1 or < 0")) end
    theta0   = convert(Float64, theta0)
    theta1   = convert(Float64, theta1)
    theta2   = convert(Float64, theta2)
    logscale = convert(Bool, logscale)
    cv       = convert(Float64, cv)
    n        = convert(Int, n)
    alpha    = convert(Float64, alpha)

    if logscale
        ltheta1 = log(theta1)
        ltheta2 = log(theta2)
        diffm   = log(theta0)
        sd      = cv2sd(cv)    # sqrt(ms)
    else
        ltheta1 = theta1;
        ltheta2 = theta2;
        diffm   = theta0;
        sd      = cv;
    end

    return powertostint(alpha,  ltheta1, ltheta2, diffm, sd, n, design, method)
end #bepower
