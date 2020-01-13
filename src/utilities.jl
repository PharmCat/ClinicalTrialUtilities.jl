#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                             UTILITIES
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#CV2se
"""
    sdfromcv(cv::Real)::AbstractFloat

LnSD from CV.
"""
@inline function sdfromcv(cv::T)::AbstractFloat where T <: Real
    return sqrt(varfromcv(cv))
end
"""
    varfromcv(cv::Real)::AbstractFloat

LnVariance from CV.
"""
function varfromcv(cv::T)::AbstractFloat where T <: Real
     return log(1+cv^2)
end
"""
    cvfromvar(σ²::Real)::AbstractFloat

CV from variance.
"""
function cvfromvar(σ²::T)::AbstractFloat where T <: Real
    return sqrt(exp(σ²)-1)
end
#CV2se
"""
    cvfromsd(σ::Real)::AbstractFloat

CV from variance.
"""
function cvfromsd(σ::T)::AbstractFloat where T <: Real
    return sqrt(exp(σ^2)-1)
end

"""
    cvfromci(;alpha = 0.05, theta1 = 0.8, theta2 = 1.25, n, design=:d2x2, mso=false, cvms=false)::AbstractFloat

CV from bioequivalence confidence inerval.

**alpha** - Alpha (o < alpha < 1)  (default=0.05);

**beta** - Beta (o < beta < 1) (default=0.2); power = 1 - beta;

**theta1** - Lower Level (default=0.8);

**theta2** - Upper level (default=1.25);

**n** - subject n;

**design** - trial design;
- :parallel
- :d2x2 (default)
- :d2x2x4
- :d2x4x4
- :d2x3x3
- :d2x4x2
- :d3x3
- :d3x6x3

**mso**

Calculate MS only
- false(default)
- true


Deprecated:

**cvms**
Calculate CV and MS
- false(default)
- true
"""
function cvfromci(;alpha = 0.05, theta1 = 0.8, theta2 = 1.25, n, design=:d2x2, mso=false, cvms=false)::AbstractFloat
    d     = Design(design)
    df    = d.df(n)
    if df < 1 throw(ArgumentError("df < 1, check n!")) end
    sef = sediv(d, n)
    ms = ((log(theta2/theta1)/2/quantile(TDist(df), 1-alpha))/sef)^2
    if mso return ms end
    return cvfromvar(ms)
end
function cvfromci(theta1, theta2, n; alpha = 0.05, design=:d2x2, mso=false, cvms=false)
    return cvfromci(;alpha = alpha, theta1 = theta1, theta2 = theta2, n = n, design = design, mso = mso, cvms = cvms)
end

"""
    pooledcv(data::DataFrame; cv=:cv, df=:df, alpha=0.05, returncv=true)::ConfInt

Pooled CV from multiple sources.

**data**::DataFrame - Dataframe with CV data

**cv**::Symbol - CV column in dataframe

**df**::Symbol - DF column in dataframe

**alpha** - Alpha for var/cv confidence interval.

**returncv** - Return CV or var:

- true  - return cv
- false - return var

"""
function pooledcv(data::DataFrame; cv = :cv, df = :df, alpha::Real = 0.05, returncv::Bool = true)::ConfInt
    if isa(cv, String)  cv = Symbol(cv) end
    if isa(df, String)  df = Symbol(df) end
    return pooledcv(data[!, cv], data[!, df]; alpha = alpha, returncv = returncv)
end
"""
    pooledcv(cv::Vector, df::Vector; alpha = 0.05, returncv = true)::ConfInt

Pooled CV from multiple sources.

**cv** - CV Vector

**df** - DF Vector

**alpha** - Alpha for var/cv confidence interval.

**returncv** - Return CV or var:

- true  - return cv
- false - return var

"""
function pooledcv(cv::Vector, df::Vector; alpha = 0.05, returncv = true)::ConfInt
    tdf    = sum(df)
    result = sum(varfromcv.(cv) .* df)/tdf
    CHSQ   = Chisq(tdf)
    if returncv return ConfInt(cvfromvar(result*tdf/quantile(CHSQ, 1-alpha/2)), cvfromvar(result*tdf/quantile(CHSQ, alpha/2)), cvfromvar(result))
    else ConfInt(result*tdf/quantile(CHSQ, 1-alpha/2), result*tdf/quantile(CHSQ, alpha/2), result)
    end
end
"""
    pooledcv(cv::Vector, n::Vector, design::Vector; alpha = 0.05, returncv = true)::ConfInt

Pooled CV from multiple sources.

**cv** - CV Vector

**n** - n Vector

**design** - design Vector

**alpha** - Alpha for var/cv confidence interval.

**returncv** - Return CV or var:

- true  - return cv
- false - return var

"""
function pooledcv(cv::Vector, n::Vector, design::Vector; alpha = 0.05, returncv = true)::ConfInt
    d  = Design.(design)
    df = Array{Int, 1}(undef, length(d))
    for i = 1:length(d)
        df[i] = d[i].df(n[i])
    end
    return pooledcv(cv, df; alpha = alpha, returncv = returncv)
end
