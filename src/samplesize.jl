# Clinical Trial Utilities
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)

#ctsamplen, ctpower, besamplen, bepower
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
    function Proportion(x::Int, n::Int)
        if x > n throw(ArgumentError("Error: X > N!")) end
        new(x, n)::Proportion
    end
end
struct Probability <: AbstractSimpleProportion
    p::Float64
    function Probability(p::Float64)
        if p ≥ 1.0 || p ≤ 0.0 throw(ArgumentError("Probability can't be ≥ 1.0 or ≤ 0.0!")) end
        new(p::Float64)::Probability
    end
    function Probability(p::Proportion)
        new(p.x/p.n)::Probability
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
"""
    ctsamplen(;param::Symbol, type::Symbol, group::Symbol = :notdef,
        alpha::Real = 0.05, beta::Real = 0.2, diff::Real = 0, sd::Real = 0,
        a::Real = 0, b::Real = 0,
        k::Real = 1, logscale::Bool = true)::TaskResult

Clinical trial sample size estimation.

**param (Parameter type):**
- :mean - Means;
- :prop - Proportions;
- :or - Odd Ratio;

**type (Hypothesis type):**
- :ea - Equality;
- :ei - Equivalencens;
- :ns - Non-Inferiority / Superiority (!one-sided hypothesis!);
- :mcnm - McNemar's Equality test;

**group (group num):**
- :one - One sample;
- :two - Two sample, result is for one group, second group size = n * k;

**alpha** - Alpha (o < α < 1)  (default=0.05);

**beta** - Beta (o < β < 1) (default=0.2); power = 1 - β;

**diff** - difference/equivalence margin/non-inferiority/superiority margin;

**sd** - Standard deviation (σ, for Means only);

**a** -  Group A  (μ₀/p₀) - Test group;

**b** -  Group B (μ₁/p₁) - Reference group or reference value;

**k** - Na/Nb (after sample size estimation second group size: Na=k*Nb, only for two sample design) (default=1);

**logscale** - diff is log transformed for OR:
- true  - diff is already in log-scale, no transformation required (default);
- false - diff is not in log-scale, will be transformed;
"""
function ctsamplen(;param::Symbol, type::Symbol, group::Symbol = :notdef, alpha::Real = 0.05, beta::Real = 0.2, diff::Real = 0, sd::Real = 0, a::Real = 0, b::Real = 0, k::Real = 1, logscale::Bool = true)::TaskResult
    if alpha >= 1 || alpha <= 0 || beta >= 1 || beta <= 0  throw(ArgumentError("sampleSize: alpha and beta sould be > 0 and < 1.")) end
    if param == :prop && !(group == :one || group == :two || type == :mcnm)  throw(ArgumentError("sampleSize: group should be defined or mcnm type.")) end
    if sd ≤ 0 && param == :mean throw(ArgumentError("SD can't be ≤ 0.0!")) end
    if k <= 0 throw(ArgumentError("Constant k can't be ≤ 0.0!")) end
    if param == :mean
        if group == :one
            if type == :ea
                n = one_mean_equality(a, b, sd, alpha, beta)
                task = CTask(DiffMean(Mean(a, sd), b), OneGroup(), Equality(), SampleSize(beta), alpha)
            elseif type == :ei
                n = one_mean_equivalence(a, b, sd, diff, alpha, beta)
                task = CTask(DiffMean(Mean(a, sd), b),  OneGroup(), Equivalence(b-diff, b+diff), SampleSize(beta), alpha)
            elseif type == :ns
                n = one_mean_superiority(a, b, sd, diff, alpha, beta)
                task = CTask(DiffMean(Mean(a, sd), b), OneGroup(), Superiority(b + diff, diff), SampleSize(beta), alpha)
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
                    task = CTask(DiffProportion(Probability(a), b), OneGroup(), Equality(), SampleSize(beta), alpha)
                elseif type == :ei
                    n = one_proportion_equivalence(a, b, diff, alpha, beta)
                    task = CTask(DiffProportion(Probability(a), b), OneGroup(), Equivalence(b-diff, b+diff), SampleSize(beta), alpha)
                elseif type == :ns
                    n = one_proportion_superiority(a, b, diff, alpha, beta)
                    task = CTask(DiffProportion(Probability(a), b), OneGroup(), Superiority(b+diff, diff), SampleSize(beta), alpha)
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
function pdiffnsn(a, b, diff; alpha = 0.05, beta = 0.2, k = 1.0)
    t = CTask(DiffProportion(Probability(a), Probability(b)), Parallel(), Superiority(diff, diff), SampleSize(beta), alpha, k)
    ctsamplen(t)
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
#------------
function ctsamplen(t::CTask{T, D, Bioequivalence, SampleSize}; method::Symbol = :owenq)  where T where D
    return besamplen(t; method = method)
end

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#clinical trial power main function
"""
    function ctpower(;param::Symbol, type::Symbol, group::Symbol = :notdef,
        alpha::Real = 0.05, n::Real = 0, diff::Real = 0, sd::Real = 0,
        a::Real = 0, b::Real = 0,
        k::Real = 1, logscale::Bool = true)::TaskResult

Clinical trial power estimation.

**param (Parameter type):**
- :mean - Means;
- :prop - Proportions;
- :or   - Odd Ratio;

**type (Hypothesis type):**
- :ea   - Equality;
- :ei   - Equivalence;
- :ns   - Non-Inferiority / Superiority;
- :mcnm - McNemar's Equality test;

**group (group num):**
- :one - one sample;
- :two - Two sample;

**alpha** - Alpha (0 < α < 1)  (default=0.05);

**n** - Subjects number;

**diff** - difference/equivalence margin/non-inferiority/superiority margin;

**sd** - Standard deviation (σ, for Means only);

**a** -  Group A  (μ₀/p₀) - Test group;

**b** -  Group B (μ₁/p₁) - Reference group or reference value;

**k** - Na/Nb (after sample size estimation second group size: Na=k*Nb, only for two sample design) (default=1);

**logscale** - diff is log transformed for OR:
- true  - diff is already in log-scale, no transformation required (default);
- false - diff is not in log-scale, will be transformed;
"""
function ctpower(;param::Symbol, type::Symbol, group::Symbol = :notdef, alpha::Real = 0.05, n::Real = 0, diff::Real = 0, sd::Real = 0, a::Real = 0, b::Real = 0,  k::Real = 1, logscale::Bool = true)::TaskResult
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
#------------

function ctpower(t::CTask{T, D, Bioequivalence, Power}; method::Symbol = :owenq)  where T where D
    return bepower(t; method = method)
end
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
"""
    besamplen(;alpha::Real=0.05, beta::Real=0.2,
        theta0::Real=0.95, theta1::Real=0.8, theta2::Real=1.25,
        cv::Real=0.0, sd::Real=0.0, design::Symbol=:d2x2,
        method::Symbol=:owenq, logscale::Bool=true)::TaskResult

Bioequivalence sample size estimation.

**alpha** - Alpha (o < α < 1)  (default=0.05);

**beta** - Beta (o < β < 1) (default=0.2); power = 1 - β;

**theta0** - T/R Ratio (default=0.95);

**theta1** - Lower Level (default=0.8);

**theta2** - Upper level (default=1.25);

**cv** - coefficient of variation;

**logscale** - theta1, theta2, theta0: if true - make log transformation (default true);

**design** - trial design;
- :parallel
- :d2x2 (default)
- :d2x2x4
- :d2x4x4
- :d2x3x3
- :d2x4x2
- :d3x3
- :d3x6x3

**method** - calculating method: Owen's Q Function | NonCentral T | Shifted;
- :owenq (default)
- :nct
- :shifted

"""
function besamplen(;alpha::Real=0.05, beta::Real=0.2, theta0::Real=0.95, theta1::Real=0.8, theta2::Real=1.25, cv::Real=0.0, sd::Real=0.0, design::Symbol=:d2x2, method::Symbol=:owenq, logscale::Bool=true)::TaskResult
    if !(0 < beta < 1) throw(ArgumentError("Beta ≥ 1.0 or ≤ 0.0!")) end
    if !(0 < alpha < 1) throw(ArgumentError("Alfa ≥ 1.0 or ≤ 0.0!")) end
    if !(theta2 > theta0 > theta1)  throw(ArgumentError("!(theta2 > thetao > theta1), check settings!")) end
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
        sd     = sdfromcv(cv)
    end

    task = CTask(DiffMean(Mean(theta0, sd), Mean(0, sd)), Crossover(design), Bioequivalence(theta1, theta2), SampleSize(beta), alpha)

    estn, estpow = samplentostint(alpha, theta1, theta2, theta0, sd, beta, design, method)
    return TaskResult(task, method, estn)
end

function besamplen(t::CTask; method::Symbol = :owenq)
    design = designtype(t.design)
    σ      = t.param.a.sd
    α      = t.alpha
    β      = t.objective.val
    θ₁     = t.hyp.llim
    θ₂     = t.hyp.ulim
    δ      = t.param.a.m - t.param.b.m
    estn, estpow = samplentostint(α, θ₁, θ₂, δ, σ, β, design, method)
    return TaskResult(t, method, estn)
end

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
"""
    bepower(;alpha::Real=0.05, theta1::Real=0.8, theta2::Real=1.25, theta0::Real=0.95,
        cv::Real=0.0, sd::Real=0.0, n::Int=0,
        design::Symbol=:d2x2, method::Symbol=:owenq, logscale::Bool=true)::TaskResult

Bioequivalence power estimation.

**alpha** - Alpha (0 < α < 1)  (default=0.05);

**theta1** - Lower Level (default=0.8);

**theta2** - Upper level (default=1.25);

**theta0** - T/R Ratio (default=0.95);

**cv** - coefficient of variation;

**n** - subject number;

**logscale** - theta1, theta2, theta0: if true - make log transformation (default true);

**design** - trial design;
- :parallel
- :d2x2 (default)
- :d2x2x4
- :d2x4x4
- :d2x3x3
- :d2x4x2
- :d3x3
- :d3x6x3

**method** - calculating method: Owen's Q Function | NonCentral T, Shifted;
- :owenq (default)
- :nct
- :shifted
"""
function bepower(;alpha::Real=0.05, theta1::Real=0.8, theta2::Real=1.25, theta0::Real=0.95, cv::Real=0.0, sd::Real=0.0, n::Int=0, design::Symbol=:d2x2, method::Symbol=:owenq, logscale::Bool=true)::TaskResult
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
        sd     = sdfromcv(cv)    # sqrt(ms)
    end

    task = CTask(DiffMean(Mean(theta0, sd), Mean(0, sd)), Crossover(design), Bioequivalence(theta1, theta2), Power(n), alpha)

    df         = task.design.df(task.objective.val)
    σ̵ₓ         = task.param.a.sd*sediv(task.design, task.objective.val)
    powertostf = powertostintf(method)
    pow        = powertostf(task.alpha, task.hyp.llim, task.hyp.ulim, task.param.a.m, σ̵ₓ, df)
    #pow =  powertostint(alpha,  theta1, theta2, theta0, sd, n, design, method)

    return TaskResult(task, method, pow)
end #bepower

function bepower(t::CTask; method::Symbol = :owenq)
    df         = t.design.df(t.objective.val)
    σ̵ₓ         = t.param.a.sd*sediv(t.design, t.objective.val)
    α          = t.alpha
    θ₁         = t.hyp.llim
    θ₂         = t.hyp.ulim
    δ          = t.param.a.m - t.param.b.m
    powertostf = powertostintf(method)
    pow        = powertostf(α, θ₁, θ₂, δ, σ̵ₓ, df)
    return TaskResult(t, method, pow)
end
