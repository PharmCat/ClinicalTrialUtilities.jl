# Clinical Trial Utilities
# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)

# ms = ss / df
# sd = ms^2
# se_diff = se = sd * sqrt((1/N1 + ... + 1/Nn)*bkni) = sqrt(ms*(1/N1 + ... + 1/Nn)*bkni)
# CI bounds = Diff +- t(df, alpha)*se


function samplentostint(α::T, θ₁::T, θ₂::T, δ::T, σ::T, β::T, design::Symbol, method::Symbol)::Tuple{Float64, Float64} where T <: AbstractFloat
     #values for approximate n
    td = (θ₁ + θ₂) / 2
    rd = abs(θ₁ - θ₂) / 2

    d0 = δ - td
    # approximate n
    n₀::Int = convert(Int, ceil(two_mean_equivalence(zero(T), d0, sqrt(σ^2 * 2), rd, α, β, one(T))) * 2)
    tp = one(typeof(β)) - β   #target power
    if n₀ < 4 n₀ = 4 end
    if n₀ > 5000 n₀ = 5000 end

    d     = Design(design) #dffunc if generic funtion with 1 arg return df
    df    = d.df(n₀)
    σ̵ₓ    = σ * sediv(d, n₀)
    if df < 1 throw(ArgumentError("powertostint: df < 1")) end

    powertostf = powertostintf(method) #PowerTOST function

    pow  = powertostf(α, θ₁, θ₂, δ, σ̵ₓ, df)
    np   = 2
    powp = pow
    if pow > tp
        while (pow > tp)
            np   = n₀
            powp = pow
            n₀   = n₀ - 2
            if n₀ < 4 break end #n0, pow end
            df    = d.df(n₀)
            σ̵ₓ    = σ * sediv(d, n₀)
            pow  = powertostf(α, θ₁, θ₂, δ, σ̵ₓ, df)
        end
        estpower = powp
        estn     = np
    elseif pow < tp
        while (pow < tp)
            np   = n₀
            powp = pow
            n₀   = n₀ + 2
            df    = d.df(n₀)
            σ̵ₓ    = σ * sediv(d, n₀)
            pow = powertostf(α, θ₁, θ₂, δ, σ̵ₓ, df)
            if n₀ > 10000  break end # n0, pow end
        end
        estpower = pow
        estn     = n₀
    else
        estpower = pow
        estn     = n₀
    end
    return estn, estpower
end

#  δ  - difference/theta0
#  σ  - SD
#  σ̵ₓ - SE/SEM
#  α  - alpha
#  θ₁ - theta1
#  θ₂ - theta2

function powertostintf(method::Symbol)::Function
    if method     == :owenq
        return   powertost_owenq
    elseif method == :nct
        return     powertost_nct
    elseif method == :shifted
        return powertost_shifted
    else
         throw(ArgumentError("method not known!"))
    end
end 

function powertost_owenq(α::T, θ₁::T, θ₂::T, δ::T, σ̵ₓ::T, df::Real)::T where T <: AbstractFloat
    tval            = quantile(TDist(df), one(T) - α)
    delta1          = (δ - θ₁) / σ̵ₓ
    delta2          = (δ - θ₂) / σ̵ₓ
    R               = (delta1 - delta2) * sqrt(df) / (tval + tval)
    if isnan(R) R   = zero(eltype(R)) end
    if R <= zero(T) R   = Inf end
    # to avoid numerical errors in OwensQ implementation
    # 'shifted' normal approximation Jan 2015
    # former Julious formula (57)/(58) doesn't work
    if df > 10000
        #tval = qnorm(1-alpha)
        tval  = quantile(ZDIST, one(T) - α)
        return max(zero(T), (cdf(ZDIST, -tval - delta2) - cdf(ZDIST,  tval - delta1)))
        #if pwr > 0 return pwr else return 0 end
    elseif df >= 5000
        # approximation via non-central t-distribution
        return powertost_nct(α, θ₁, θ₂, δ, σ̵ₓ, df)
    else
        return max(zero(T), (owensq(df, -tval, delta2, zero(T), R) - owensq(df, tval, delta1, zero(T), R)))
    end
end 

#------------------------------------------------------------------------------
# approximation based on non-central t
function powertost_nct(α::T, θ₁::T, θ₂::T, δ::T, σ̵ₓ::T, df::Real)::T where T <: AbstractFloat
    tdist           = TDist(df)
    tval            = quantile(tdist, one(T) - α)
    delta1          = (δ - θ₁)/σ̵ₓ
    delta2          = (δ - θ₂)/σ̵ₓ
    pow             = cdf(NoncentralT(df, delta2), -tval) - cdf(NoncentralT(df, delta1), tval)
    return max(zero(T), pow)
end 

function powertost_shifted(α::T, θ₁::T, θ₂::T, δ::T, σ̵ₓ::T, df::Real)::T where T <: AbstractFloat
    tdist           = TDist(df)
    tval            = quantile(tdist, one(T) - α)
    delta1          = (δ - θ₁) / σ̵ₓ
    delta2          = (δ - θ₂) / σ̵ₓ
    if isnan(delta1) delta1 = zero(T) end
    if isnan(delta2) delta2 = zero(T) end
    pow = cdf(tdist,-tval - delta2) - cdf(tdist, tval - delta1)
    return max(zero(T), pow)
end 
