# Clinical Trial Utilities
# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)

# ms = ss / df
# sd = ms^2
# se_diff = se = sd * sqrt((1/N1 + ... + 1/Nn)*bkni) = sqrt(ms*(1/N1 + ... + 1/Nn)*bkni)
# CI bounds = Diff +- t(df, alpha)*se


#powertostint
#powerTOSTOwenQ
#approxPowerTOST
#power1TOST
#approx2PowerTOST
#cv2sd
#cv2ms
#ms2cv
#sd2cv
#designProp
#ci2cv

function samplentostint(alpha, ltheta1, ltheta2, diffm, sd, beta, design, method)
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
    return estn, estpower
end

#  δ  - difference/theta0
#  σ  - SD
#  σ̵ₓ - SE/SEM
#  α  - alpha
#  θ₁ - theta1
#  θ₂ - theta2
function powertostint(α::Real,  θ₁::Real, θ₂::Real, δ::Real, σ::Real, n::Int, design::Symbol, method::Symbol)::Float64
    d     = Design(design) #dffunc if generic funtion with 1 arg return df
    df    = d.df(n)

    if df < 1 throw(ArgumentError("powertostint: df < 1")) end

    σ̵ₓ::Float64 = σ*sediv(d, n)

    if method     == :owenq
        return   powertost_owenq(α, θ₁, θ₂, δ, σ̵ₓ, df)
    elseif method == :nct
        return     powertost_nct(α, θ₁, θ₂, δ, σ̵ₓ, df)
    elseif method == :mvt
        return     powertost_mvt(α, θ₁, θ₂, δ, σ̵ₓ, df) #not implemented
    elseif method == :shifted
        return powertost_shifted(α, θ₁, θ₂, δ, σ̵ₓ, df)
    else
         throw(ArgumentError("method not known!"))
    end
end #powerTOST

#.power.TOST
function powertost_owenq(alpha::Real, ltheta1::Real, ltheta2::Real, diffm::Real, se::Real, df::Real)::Float64
    tval::Float64   = quantile(TDist(df), 1 - alpha)
    delta1::Float64 = (diffm-ltheta1)/se
    delta2::Float64 = (diffm-ltheta2)/se
    R::Float64      = (delta1 - delta2) * sqrt(df) / (tval + tval)
    if isnan(R) R   = 0 end
    if R <= 0 R     = Inf end
    # to avoid numerical errors in OwensQ implementation
    # 'shifted' normal approximation Jan 2015
    # former Julious formula (57)/(58) doesn't work
    if df > 10000
        #tval = qnorm(1-alpha)
        tval  = quantile(ZDIST, 1 - alpha)
        #p1   = pnorm(tval-delta1)
        p1    = cdf(ZDIST,  tval - delta1)
        #p2   = pnorm(-tval-delta2)
        p2    = cdf(ZDIST, -tval - delta2)
        pwr   = p2 - p1
        if pwr > 0 return pwr else return 0 end
    elseif df >= 5000
        # approximation via non-central t-distribution
        return powertost_nct(alpha, ltheta1, ltheta2, diffm, se, df)
    end
    p1  = owensq(df, tval, delta1, 0.0, R)
    p2  = owensq(df,-tval, delta2, 0.0, R)
    pwr = p2 - p1
    if pwr > 0 return pwr else return 0 end
end #powerTOSTOwenQ

#------------------------------------------------------------------------------
# approximation based on non-central t
# .approx.power.TOST - PowerTOST
function powertost_nct(alpha::Real, ltheta1::Real, ltheta2::Real, diffm::Real, se::Real, df::Real)::Float64
    tdist           = TDist(df)
    tval::Float64   = quantile(tdist, 1-alpha)
    delta1::Float64 = (diffm-ltheta1)/se
    delta2::Float64 = (diffm-ltheta2)/se
    pow             = cdf(NoncentralT(df,delta2), -tval) - cdf(NoncentralT(df,delta1), tval)
    if pow > 0 return pow else return 0 end
end #approxPowerTOST

#.power.1TOST
function powertost_mvt(alpha::Real, ltheta1::Real, ltheta2::Real, diffm::Real, se::Real, df::Real)::Float64
    throw(ArgumentError("Method not implemented!"))
    #Method ON MULTIVARIATE t AND GAUSS PROBABILITIES IN R not implemented
    # Distributions.MvNormal - in plan
    # see  Distributions.jl/src/multivariate/mvtdist.jl
    # Multivariate t-distribution
    # Generic multivariate t-distribution class
    # mvt = MvTDist()
end

#.approx2.power.TOST
function powertost_shifted(alpha::Real, ltheta1::Real, ltheta2::Real, diffm::Real, se::Real, df::Real)::Float64
    tdist           = TDist(df)
    tval::Float64   = quantile(tdist, 1-alpha)
    delta1::Float64 = (diffm-ltheta1)/se
    delta2::Float64 = (diffm-ltheta2)/se
    if isnan(delta1) delta1 = 0 end
    if isnan(delta2) delta2 = 0 end
    pow = cdf(tdist,-tval-delta2) - cdf(tdist, tval-delta1)
    if pow > 0 return pow else return 0 end
end #approx2PowerTOST

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                             UTILITIES
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#CV2se
@inline function cv2sd(cv::Real)::Float64
    return sqrt(cv2ms(cv))
end

@inline function cv2ms(cv::Real)::Float64
     return log(1+cv^2)
end
@inline function ms2cv(ms::Real)::Float64
    return sqrt(exp(ms)-1)
end
#CV2se
@inline function sd2cv(sd::Real)::Float64
    return sqrt(exp(sd^2)-1)
end

function ci2cv(;alpha = 0.05, theta1 = 0.8, theta2 = 1.25, n, design=:d2x2, mso=false, cvms=false)::Union{Float64, Tuple{Float64, Float64}}
    d     = Design(design)
    df    = d.df(n)
    if df < 1 throw(ArgumentError("df < 1, check n!")) end

    sef = sediv(d, n)

    ms = ((log(theta2/theta1)/2/quantile(TDist(df), 1-alpha))/sef)^2
    if cvms return ms2cv(ms), ms end
    if mso return ms end
    return ms2cv(ms)
end

function pooledcv(data::DataFrame; cv=:cv, df=:df, alpha=0.05, returncv=true)::ConfInt
    if isa(cv, String)  cv = Symbol(cv) end
    if isa(df, String)  df = Symbol(df) end
    tdf = sum(data[:, df])
    result = sum(cv2ms.(data[:, cv]) .* data[:, df])/tdf
    CHSQ = Chisq(tdf)
    if returncv return ConfInt(ms2cv(result*tdf/quantile(CHSQ, 1-alpha/2)), ms2cv(result*tdf/quantile(CHSQ, alpha/2)), ms2cv(result))
    else ConfInt(result*tdf/quantile(CHSQ, 1-alpha/2), result*tdf/quantile(CHSQ, alpha/2), result)
    end
end
