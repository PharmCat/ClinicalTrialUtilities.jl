# Clinical Trial Utilities
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)

#ctsamplen, ctpower, besamplen, bepower


struct Power <: AbstractPower
    val::Int
end
struct SampleSize{T} <: AbstractSampleSize where T <: AbstractFloat
    val::T
    function SampleSize(val::T) where T <: AbstractFloat
        if val ≥ 1.0 || val ≤ 0.0 throw(ArgumentError("Beta ≥ 1.0 or ≤ 0.0!")) end
        new{T}(val)::SampleSize
    end
end

function getval(o::T) where  T <:  AbstractObjective
    o.val
end


struct TaskResult{T <: AbstractTask}
    task::T
    method::Symbol
    res::Dict
    function TaskResult(t::T, method::Symbol, result::Dict) where T <: AbstractTask
        new{typeof(t)}(t, method, result)
    end
    function TaskResult(t::T, method::Symbol, result::R) where T <: AbstractTask where R <: Real
        new{typeof(t)}(t, method, Dict(:result => result))
    end
end

@inline function Base.getproperty(t::TaskResult, f::Symbol)
    if f == :result return t.res[:result] else return getfield(t, f) end
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
    k::Real             #Group coefficient
    function CTask(param::T, design::D, hyp::H, objective::O, k::Real) where T <: AbstractParameter where D <: AbstractDesign where H <: AbstractHypothesis where O <: AbstractObjective
        #=
        if isa(hyp, Equality) && llim != ulim
            @warn "For Equality hypothesis llim and ulim can't be different! ulim set as llim!"
            ulim = llim
        end
        if k ≤ 0.0 throw(ArgumentError("Constant k can't be ≤ 0.0!")) end
        if alpha ≥ 1.0 || alpha ≤ 0.0 throw(ArgumentError("Alpha ≥ 1.0 or ≤ 0.0!")) end
        if isa(hyp, Equivalence) && diff ≤ 0 throw(ArgumentError("Diiference can't be ≤ 0.0 with Equivalence hypothesis!")) end
        =#
        new{T,D,H,O}(param, design, hyp, objective, k)::CTask{T,D,H,O}
    end
    function CTask(param::T, design::D, hyp::H, objective::O) where T <: AbstractParameter where D <: AbstractDesign where H <: AbstractHypothesis where O <: AbstractObjective
        return CTask(param, design, hyp, objective, 1)
    end
end

function setobjval!(task::CTask{T,D,H,O}, val) where T where D where H where O <: AbstractPower
    task.objective = Power(val)
end
function setobjval!(task::CTask{T,D,H,O}, val) where T where D where H where O <: AbstractSampleSize
    task.objective = SampleSize(val)
end

"""
    ctask(;param::Symbol, hyp::Symbol, group::Symbol = :notdef, alpha::Real = 0.05, beta::Real = 0.2, k::Real = 1, kw...)::CTask
"""
function ctask(;param::Symbol, hyp::Symbol, group::Symbol = :notdef, alpha::Real = 0.05, k::Real = 1, kw...)::CTask
    kwk      = keys(kw)
    a        = kw[:a]
    b        = kw[:b]

    if :power in kwk && :n in kwk
         throw(ArgumentError("only beta or only n keyword should be used!"))
    elseif :beta in kwk
        objf = SampleSize
        objv = kw[:beta]
    elseif :n in kwk
        objf = Power
        objv = kw[:n]
    else
        throw(ArgumentError("beta or n keyword should be used!"))
    end

    if (hyp == :ei || hyp == :ns) && :diff ∉ kwk
        throw(ArgumentError("diff keyword should be specified when ei or ns hypothesis used!"))
    else
        diff     = kw[:diff]
    end

    if param == :mean

        sd = kw[:sd]
        if group == :one
            if hyp == :ea
                task = CTask(Mean(kw[:a], sd), OneGroup(), Equality(kw[:b], alpha), objf(objv))
            elseif hyp == :ei
                task = CTask(Mean(kw[:a], sd), OneGroup(), Equivalence(kw[:b] - diff, kw[:b] + diff, alpha), objf(objv))
            elseif hyp == :ns
                task = CTask(Mean(kw[:a], sd), OneGroup(), Superiority(kw[:b] + diff, diff, alpha), objf(objv))
            else throw(ArgumentError("Keyword hyp unknown!")) end
        elseif group == :two
            if hyp == :ea
                task = CTask(DiffMean(Mean(kw[:a], sd), Mean(kw[:b], sd)), Parallel(), Equality(alpha), objf(objv), k)
            elseif hyp == :ei
                task = CTask(DiffMean(Mean(kw[:a], sd), Mean(kw[:b], sd)), Parallel(), Equivalence(-diff, diff, alpha), objf(objv), k)
            elseif hyp == :ns
                task = CTask(DiffMean(Mean(kw[:a], sd), Mean(kw[:b], sd)),  Parallel(), Superiority(diff, diff, alpha), objf(objv), k)
            else throw(ArgumentError("Keyword hyp unknown!")) end
        elseif group == :co

        else throw(ArgumentError("Keyword group unknown!")) end
    elseif param == :prop
        if group == :one
            if hyp == :ea
                task = CTask(Proportion(kw[:a]), OneGroup(), Equality(kw[:b], alpha), objf(objv))
            elseif hyp == :ei
                task = CTask(Proportion(kw[:a]), OneGroup(), Equivalence(kw[:b] - diff,  kw[:b] + diff, alpha), objf(objv))
            elseif hyp == :ns
                task = CTask(Proportion(kw[:a]), OneGroup(), Superiority(kw[:b] + diff, diff, alpha), objf(objv))
            else throw(ArgumentError("Keyword hyp unknown!")) end
        elseif group == :two
            if hyp == :ea
                task = CTask(DiffProportion(Proportion(kw[:a]), Proportion(kw[:b])), Parallel(), Equality(alpha), objf(objv), k)
            elseif hyp == :ei
                task = CTask(DiffProportion(Proportion(kw[:a]), Proportion(kw[:b])), Parallel(), Equivalence(-diff, diff, alpha), objf(objv), k)
            elseif hyp == :ns
                task = CTask(DiffProportion(Proportion(kw[:a]), Proportion(kw[:b])), Parallel(), Superiority(diff, diff, alpha), objf(objv), k)
            else throw(ArgumentError("Keyword hyp unknown!")) end
        elseif group == :co
            if hyp == :mcnm
                task = CTask(DiffProportion(Proportion(kw[:a]), Proportion(kw[:b])), Parallel(),  McNemars(alpha), objf(objv), k)
            else throw(ArgumentError("Keyword hyp unknown!")) end
        else throw(ArgumentError("Keyword group unknown!")) end
    elseif param == :or

        logscale = kw[:logscale]
        if group == :two || group == :notdef
            if hyp == :ea
                task = CTask(OddRatio(Proportion(kw[:a]), Proportion(kw[:b])), Parallel(),  Equality(alpha), objf(objv), k)
            elseif hyp == :ei
                if !logscale diff = log(diff) end
                task = CTask(OddRatio(Proportion(kw[:a]), Proportion(kw[:b])), Parallel(), Equivalence(-diff, diff, alpha), objf(objv), k)
            elseif hyp == :ns
                if !logscale diff = log(diff) end
                task = CTask(OddRatio(Proportion(kw[:a]), Proportion(kw[:b])), Parallel(), Superiority(diff,  diff, alpha), objf(objv), k)
            else throw(ArgumentError("Keyword hyp unknown!")) end
        elseif group == :co

        else throw(ArgumentError("Keyword group unknown!")) end
    elseif param == :cox
        if hyp == :ea
            task = CTask(CoxHazardRatio(kw[:a], kw[:p]), Parallel(),  Equality(alpha), objf(objv), k)
        elseif hyp == :ei
            task = CTask(CoxHazardRatio(kw[:a], kw[:p]), Parallel(),  Equivalence(-diff, diff, alpha), objf(objv), k)
        elseif hyp == :ns
            task = CTask(CoxHazardRatio(kw[:a], kw[:p]), Parallel(),  Superiority(diff,  diff, alpha), objf(objv), k)
        else throw(ArgumentError("Keyword hyp unknown!")) end
    else throw(ArgumentError("Keyword param unknown!")) end
    return task
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
    if alpha >= 1 || alpha <= 0 || beta >= 1 || beta <= 0  throw(ArgumentError("ctsamplen: alpha and beta sould be > 0 and < 1.")) end
    if param == :prop && !(group == :one || group == :two || type == :mcnm)  throw(ArgumentError("ctsamplen: group should be defined or mcnm type.")) end
    if sd ≤ 0 && param == :mean throw(ArgumentError("SD can't be ≤ 0.0!")) end
    if k <= 0 throw(ArgumentError("Constant k can't be ≤ 0.0!")) end

    if type == :mcnm
        group = :co
    end

    kwd  = Dict(:beta => beta, :diff => diff, :sd => sd, :a => a, :b => b, :logscale => logscale)

    task = ctask(;param = param, hyp = type, group = group, alpha = alpha, k = k, kwd...)

    return ctsamplen(task)
end #sampleSize

function ctsamplen(t::CTask{T, D, H, O}) where T <: Mean where D <: AbstractDesign where H <: Equality where O <: AbstractSampleSize
    return TaskResult(t, :chow, one_mean_equality(t.param.m, refval(t.hyp), t.param.sd, t.hyp.alpha, t.objective.val))
end
function ctsamplen(t::CTask{T, D, H, O}) where T <: Mean where D <: AbstractDesign where H <: Equivalence where O <: AbstractSampleSize
    return TaskResult(t, :chow, one_mean_equivalence(t.param.m, refval(t.hyp), t.param.sd, mdiff(t.hyp), t.hyp.alpha, t.objective.val))
end
function ctsamplen(t::CTask{T, D, H, O}) where T <: Mean where D <: AbstractDesign where H <: Superiority where O <: AbstractSampleSize
    return TaskResult(t, :chow, one_mean_superiority(t.param.m, refval(t.hyp), t.param.sd, t.hyp.diff, t.hyp.alpha, t.objective.val))
end
#---
function ctsamplen(t::CTask{T, D, H, O}) where T <: DiffMean where D <: AbstractDesign where H <: Equality where O <: AbstractSampleSize
    return TaskResult(t, :chow, two_mean_equality(t.param.a.m, t.param.b.m, t.param.a.sd, t.hyp.alpha, t.objective.val, t.k))
end
function ctsamplen(t::CTask{T, D, H, O}) where T <: DiffMean  where D <: AbstractDesign where H <: Equivalence where O <: AbstractSampleSize
    return TaskResult(t, :chow, two_mean_equivalence(t.param.a.m, t.param.b.m, t.param.a.sd,  mdiff(t.hyp), t.hyp.alpha, t.objective.val, t.k))
end
function ctsamplen(t::CTask{T, D, H, O}) where T <: DiffMean  where D <: AbstractDesign where H <: Superiority where O <: AbstractSampleSize
    return TaskResult(t, :chow, two_mean_superiority(t.param.a.m, t.param.b.m, t.param.a.sd, t.hyp.diff, t.hyp.alpha, t.objective.val, t.k))
end
#------
function ctsamplen(t::CTask{T, D, H, O}) where T <: Proportion where D  <: OneGroup where H <: Equality where O <: AbstractSampleSize
    return TaskResult(t, :chow, one_proportion_equality(getval(t.param), refval(t.hyp), t.hyp.alpha, t.objective.val))
end
function ctsamplen(t::CTask{T, D, H, O}) where T <: Proportion   where D <: OneGroup where H <: Equivalence where O <: AbstractSampleSize
    return TaskResult(t, :chow, one_proportion_equivalence(getval(t.param), refval(t.hyp),  mdiff(t.hyp), t.hyp.alpha, t.objective.val))
end
function ctsamplen(t::CTask{T, D, H, O}) where T <: Proportion  where D <: OneGroup where H <: Superiority where O <: AbstractSampleSize
    return TaskResult(t, :chow, one_proportion_superiority(getval(t.param), refval(t.hyp), t.hyp.diff, t.hyp.alpha, t.objective.val))
end
#---
function ctsamplen(t::CTask{T, D, H, O}) where T <: DiffProportion{P, P} where P <: AbstractSimpleProportion  where D <: AbstractDesign where H <: Equality where O <: AbstractSampleSize
    return TaskResult(t, :chow, two_proportion_equality(getval(t.param.a), getval(t.param.b), t.hyp.alpha, t.objective.val, t.k))
end
function ctsamplen(t::CTask{T, D, H, O}) where T <: DiffProportion{P, P} where P <: AbstractSimpleProportion  where D <: AbstractDesign where H <: Equivalence where O <: AbstractSampleSize
    return TaskResult(t, :chow, two_proportion_equivalence(getval(t.param.a), getval(t.param.b),  mdiff(t.hyp), t.hyp.alpha, t.objective.val, t.k))
end
function pdiffnsn(a, b, diff; alpha = 0.05, beta = 0.2, k = 1.0)
    t = CTask(DiffProportion(Proportion(a), Proportion(b)), Parallel(), Superiority(diff, diff, alpha), SampleSize(beta), k)
    ctsamplen(t)
end
function ctsamplen(t::CTask{T, D, H, O}) where T <: DiffProportion{P, P} where P <: AbstractSimpleProportion  where D <: AbstractDesign where H <: Superiority where O <: AbstractSampleSize
    return TaskResult(t, :chow, two_proportion_superiority(getval(t.param.a), getval(t.param.b), t.hyp.diff, t.hyp.alpha, t.objective.val, t.k))
end
#------
function ctsamplen(t::CTask{T, D, H, O}) where T <: OddRatio  where D <: AbstractDesign where H <: Equality where O <: AbstractSampleSize
    return TaskResult(t, :chow, or_equality(getval(t.param.a), getval(t.param.b), t.hyp.alpha, t.objective.val, t.k))
end
function ctsamplen(t::CTask{T, D, H, O}) where T <: OddRatio  where D <: AbstractDesign where H <: Equivalence where O <: AbstractSampleSize
    return TaskResult(t, :chow, or_equivalence(getval(t.param.a), getval(t.param.b),  mdiff(t.hyp), t.hyp.alpha, t.objective.val, t.k))
end
function ctsamplen(t::CTask{T, D, H, O}) where T <: OddRatio  where D <: AbstractDesign where H <: Superiority where O <: AbstractSampleSize
    return TaskResult(t, :chow, or_superiority(getval(t.param.a), getval(t.param.b), t.hyp.diff, t.hyp.alpha, t.objective.val, t.k))
end
#------------
function ctsamplen(t::CTask{T, D, H, O}) where T <: DiffProportion{P, P} where P <: AbstractProportion  where D <: AbstractDesign where H <: McNemars where O <: AbstractSampleSize
    return TaskResult(t, :chow, mcnm(getval(t.param.a), getval(t.param.b), t.hyp.alpha, t.objective.val))
end
#------------
function ctsamplen(t::CTask{T, D, Bioequivalence, O}; method::Symbol = :owenq)  where T where D where O <: AbstractSampleSize
    return besamplen(t; method = method)
end
#------------

#------------
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

    if type == :mcnm
        group = :co
    end

    kwd  = Dict(:n => n, :diff => diff, :sd => sd, :a => a, :b => b, :logscale => logscale)

    task = ctask(;param = param, hyp = type, group = group, alpha = alpha, k = k, kwd...)

    return ctpower(task)
end #ctpower
#---
function ctpower(t::CTask{T, D, H, O}) where T <: Mean where D <: AbstractDesign where H <: Equality where O <: Power
    return TaskResult(t, :chow,  one_mean_equality_pow(t.param.m, refval(t.hyp), t.param.sd, t.hyp.alpha, t.objective.val))
end
function ctpower(t::CTask{T, D, H, O}) where T <: Mean where D <: AbstractDesign where H <: Equivalence where O <: Power
    return TaskResult(t, :chow,  one_mean_equivalence_pow(t.param.m, refval(t.hyp), t.param.sd, mdiff(t.hyp), t.hyp.alpha, t.objective.val))
end
function ctpower(t::CTask{T, D, H, O}) where T <: Mean where D <: AbstractDesign where H <: Superiority where O <: Power
    return TaskResult(t, :chow,  one_mean_superiority_pow(t.param.m, refval(t.hyp), t.param.sd, t.hyp.diff, t.hyp.alpha, t.objective.val))
end
#------
function ctpower(t::CTask{T, D, H, O}) where T <: DiffMean  where D <: AbstractDesign where H <: Equality where O <: Power
    return TaskResult(t, :chow, two_mean_equality_pow(t.param.a.m, t.param.b.m, t.param.a.sd, t.hyp.alpha, t.objective.val, t.k))
end
function ctpower(t::CTask{T, D, H, O}) where T <: DiffMean where D <: AbstractDesign where H <: Equivalence where O <: Power
    return TaskResult(t, :chow, two_mean_equivalence_pow(t.param.a.m, t.param.b.m, t.param.a.sd, mdiff(t.hyp), t.hyp.alpha, t.objective.val, t.k))
end
function ctpower(t::CTask{T, D, H, O}) where T <: DiffMean where D <: AbstractDesign where H <: Superiority where O <: Power
    return TaskResult(t, :chow, two_mean_superiority_pow(t.param.a.m, t.param.b.m, t.param.a.sd, t.hyp.diff, t.hyp.alpha, t.objective.val, t.k))
end
#---
function ctpower(t::CTask{T, D, H, O})  where T <: Proportion where D <: OneGroup where H <: Equality where O <: Power
    return TaskResult(t, :chow, one_proportion_equality_pow(getval(t.param), refval(t.hyp), t.hyp.alpha, t.objective.val))
end
function ctpower(t::CTask{T, D, H, O})  where T <: Proportion where D <: OneGroup where H <: Equivalence where O <: Power
    return TaskResult(t, :chow, one_proportion_equivalence_pow(getval(t.param), refval(t.hyp), mdiff(t.hyp), t.hyp.alpha, t.objective.val))
end
function ctpower(t::CTask{T, D, H, O})  where T <: Proportion where D <: OneGroup where H <: Superiority where O <: Power
    return TaskResult(t, :chow, one_proportion_superiority_pow(getval(t.param), refval(t.hyp), t.hyp.diff, t.hyp.alpha, t.objective.val))
end
#---
function ctpower(t::CTask{T, D, H, O}) where T <: DiffProportion{P, P} where P <: AbstractSimpleProportion where D <: Parallel where H <: Equality where O <: Power
    return TaskResult(t, :chow, two_proportion_equality_pow(getval(t.param.a), getval(t.param.b), t.hyp.alpha, t.objective.val, t.k))
end
function ctpower(t::CTask{T, D, H, O}) where T <: DiffProportion{P, P} where P <: AbstractSimpleProportion where D <: Parallel where H <: Equivalence where O <: Power
    return TaskResult(t, :chow, two_proportion_equivalence_pow(getval(t.param.a), getval(t.param.b), mdiff(t.hyp), t.hyp.alpha, t.objective.val, t.k))
end
function ctpower(t::CTask{T, D, H, O}) where T <: DiffProportion{P, P} where P <: AbstractSimpleProportion where D <: Parallel where H <: Superiority where O <: Power
    return TaskResult(t, :chow, two_proportion_superiority_pow(getval(t.param.a), getval(t.param.b), t.hyp.diff, t.hyp.alpha, t.objective.val, t.k))
end
#------
function ctpower(t::CTask{T, D, H, O}) where T <: OddRatio{P} where P <:  AbstractSimpleProportion where D <: AbstractDesign where H <: Equality where O <: Power
    return TaskResult(t, :chow, or_equality_pow(getval(t.param.a), getval(t.param.b), t.hyp.alpha, t.objective.val, t.k))
end
function ctpower(t::CTask{T, D, H, O}) where T <: OddRatio{P} where P <:  AbstractSimpleProportion where D <: AbstractDesign where H <: Equivalence where O <: Power
    return TaskResult(t, :chow, or_equivalence_pow(getval(t.param.a), getval(t.param.b), mdiff(t.hyp), t.hyp.alpha, t.objective.val, t.k))
end
function ctpower(t::CTask{T, D, H, O}) where T <: OddRatio{P} where P <:  AbstractSimpleProportion where D <: AbstractDesign where H <: Superiority where O <: Power
    return TaskResult(t, :chow, or_superiority_pow(getval(t.param.a), getval(t.param.b), t.hyp.diff, t.hyp.alpha, t.objective.val, t.k))
end
#------------
function ctpower(t::CTask{T, D, H, O}) where T <: DiffProportion{P, P} where P <: AbstractProportion where D <: AbstractDesign  where H <: McNemars where O <: Power
    return TaskResult(t, :chow, mcnm_pow(getval(t.param.a), getval(t.param.b), t.hyp.alpha, t.objective.val))
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

    task = CTask(DiffMean(Mean(theta0, sd), Mean(0, sd)), Crossover(design), Bioequivalence(theta1, theta2, alpha), SampleSize(beta))

    estn, estpow = samplentostint(alpha, theta1, theta2, theta0, sd, beta, design, method)
    return TaskResult(task, method, estn)
end

function besamplen(t::CTask; method::Symbol = :owenq)
    design = designtype(t.design)
    σ      = t.param.a.sd
    α      = t.hyp.alpha
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

    task = CTask(DiffMean(Mean(theta0, sd), Mean(0, sd)), Crossover(design), Bioequivalence(theta1, theta2, alpha), Power(n))

    df         = task.design.df(task.objective.val)
    σ̵ₓ         = task.param.a.sd*sediv(task.design, task.objective.val)
    powertostf = powertostintf(method)
    pow        = powertostf(task.hyp.alpha, task.hyp.llim, task.hyp.ulim, task.param.a.m, σ̵ₓ, df)
    #pow =  powertostint(alpha,  theta1, theta2, theta0, sd, n, design, method)

    return TaskResult(task, method, pow)
end #bepower

function bepower(t::CTask; method::Symbol = :owenq)
    df         = t.design.df(t.objective.val)
    σ̵ₓ         = t.param.a.sd*sediv(t.design, t.objective.val)
    α          = t.hyp.alpha
    θ₁         = t.hyp.llim
    θ₂         = t.hyp.ulim
    δ          = t.param.a.m - t.param.b.m
    powertostf = powertostintf(method)
    pow        = powertostf(α, θ₁, θ₂, δ, σ̵ₓ, df)
    return TaskResult(t, method, pow)
end
