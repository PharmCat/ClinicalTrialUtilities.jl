# Clinical Trial Utilities
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)
# Calculation based on Chow S, Shao J, Wang H. 2008. Sample Size Calculations in Clinical Research. 2nd Ed. Chapman & Hall/CRC Biostatistics Series.
# OwensQ function rewrited from https://github.com/Detlew/PowerTOST by Detlew Labes, Helmut Schuetz, Benjamin Lang
# Licence: GNU Affero General Public License v3.0

#Ref: Chow S, Shao J, Wang H. 2008. Sample Size Calculations in Clinical Research. 2nd Ed. Chapman & Hall/CRC Biostatistics Series.
#Sample Size
#Compare Means
#One Sample
# μ₀ = m₀ = Test
# m₁ = Reference value
# σ  = SD - Standart deviation
# δ  = difference μ
function one_mean_equality(μ₀::Real, μ₁::Real, σ::Real, α::Float64, β::Float64)::Float64
    return ((quantile(ZDIST, 1 - α / 2) + quantile(ZDIST, 1 - β)) * σ / (μ₀ - μ₁))^2
end
function one_mean_equivalence(μ₀::Real, μ₁::Real, σ::Real, δ::Real, α::Float64, β::Float64)::Float64
    return (σ*(quantile(ZDIST, 1 - α) + quantile(ZDIST, 1 - β/2))/(δ - abs(μ₀ - μ₁)))^2
end
function one_mean_superiority(μ₀::Real, μ₁::Real, σ::Real, δ::Real, α::Float64, β::Float64)::Float64 #Non-inferiority / Superiority
    return (σ*(quantile(ZDIST, 1 - α) + quantile(ZDIST, 1 - β))/(μ₀ - μ₁ - δ))^2
end
#Two Sample
# μ₀ = μA - Group A; μ₁ = μB - Group B
function two_mean_equality(μ₀::Real, μ₁::Real, σ::Real, α::Float64, β::Float64, k::Real)
    return (1 + 1 / k)*(σ*(quantile(ZDIST, 1 - α / 2) + quantile(ZDIST, 1 - β))/(μ₀ - μ₁))^2
end
function two_mean_equivalence(μ₀::Real, μ₁::Real, σ::Real, δ::Real, α::Float64, β::Float64, k::Real)
    return (1 + 1 / k)*(σ*(quantile(ZDIST, 1 - α) + quantile(ZDIST, 1 - β / 2))/(abs(μ₀ - μ₁) - δ))^2
end
function twoSampleMeanNS(μ₀::Real, μ₁::Real, σ::Real, δ::Real; alpha=0.05, beta=0.2, k=1) #Non-inferiority / Superiority
    return (1 + 1 / k)*(σ*(quantile(ZDIST, 1 - alpha) + quantile(ZDIST, 1 - beta))/(μ₀ - μ₁ - δ))^2
end
#Compare Proportion
#One Sample
function one_proportion_equality(p₀::Float64, p₁::Float64, α::Float64, β::Float64)::Float64
    return p₀*(1 - p₀)*((quantile(ZDIST, 1 - α / 2)+quantile(ZDIST, 1 - β))/(p₀ - p₁))^2
end
function oneProportionEquivalence(p₀::Float64, p₁::Float64, δ::Real; alpha=0.05, beta=0.2)
    return p₀*(1 - p₀)*((quantile(ZDIST, 1 - alpha)+quantile(ZDIST, 1 - beta / 2))/(abs(p₀ - p₁) - δ))^2
end
function oneProportionNS(p₀::Float64, p₁::Float64, δ::Real; alpha=0.05, beta=0.2)
    return p₀*(1 - p₀)*((quantile(ZDIST, 1 - alpha)+quantile(ZDIST, 1 - beta))/(p₀ - p₁ - δ))^2
end
#Two Sample
function twoProportionEquality(p₀::Float64, p₁::Float64; alpha=0.05, beta=0.2, k=1)
    return (p₀*(1-p₀)/k+p₁*(1-p₁))*((quantile(ZDIST, 1-alpha/2)+quantile(ZDIST, 1 - beta))/(p₀-p₁))^2
end
function twoProportionEquivalence(p₀::Float64, p₁::Float64, δ::Real; alpha=0.05, beta=0.2, k=1)
    return (p₀*(1-p₀)/k+p₁*(1-p₁))*((quantile(ZDIST, 1-alpha)+quantile(ZDIST, 1 - beta/2))/(δ - abs(p₀-p₁)))^2
end
function twoProportionNS(p₀::Float64, p₁::Float64, δ::Real; alpha=0.05, beta=0.2, k=1)
    return (p₀*(1-p₀)/k+p₁*(1-p₁))*((quantile(ZDIST, 1-alpha)+quantile(ZDIST, 1 - beta))/(p₀-p₁-δ))^2
end
#Odd ratio Chow S, Shao J, Wang H. 2008. Sample Size Calculations in Clinical Research. 2nd Ed. Chapman & Hall/CRC Biostatistics Series.
function orEquality(p₀::Float64, p₁::Float64; alpha=0.05, beta=0.2, k=1)
    OR=p₀*(1-p₁)/p₁/(1-p₀)
    return (1/k/p₀/(1-p₀)+1/p₁/(1-p₁))*((quantile(ZDIST, 1-alpha/2)+quantile(ZDIST, 1 - beta))/log(OR))^2
end
function orEquivalence(p₀::Float64, p₁::Float64, δ::Real; alpha=0.05, beta=0.2, k=1)
    OR=p₀*(1-p₁)/p₁/(1-p₀)
    return (1/k/p₀/(1-p₀)+1/p₁/(1-p₁))*((quantile(ZDIST, 1-alpha)+quantile(ZDIST, 1 - beta/2))/(log(OR)-δ))^2
end
function orNS(p₀::Float64, p₁::Float64, δ::Real; alpha=0.05, beta=0.2, k=1)
    OR=p₀*(1-p₁)/p₁/(1-p₀)
    return (1/k/p₀/(1-p₀)+1/p₁/(1-p₁))*((quantile(ZDIST, 1-alpha)+quantile(ZDIST, 1 - beta))/(log(OR)-δ))^2
end

#ClinicalTrialUtilities.mcnm(0.45, 0.05, alpha=0.1, beta=0.1)
#22.80907052752994
#Connor R. J. 1987. Sample size for testing differences in proportions for the paired-sample design. Biometrics 43(1):207-211. page 209.
function mcnm(p10, p01; alpha=0.05, beta=0.2)
    pdisc=p10+p01
    pdiff=p10-p01
    return ((quantile(ZDIST, 1-alpha/2)*sqrt(pdisc)+quantile(ZDIST, 1 - beta)*sqrt(pdisc-pdiff^2))/pdiff)^2
end

function mcnmP(p10, p01, n; alpha=0.05)
    pdisc=p10+p01
    pdiff=p10-p01
    x1=(pdiff*sqrt(n)-quantile(ZDIST, 1-alpha/2)*sqrt(pdisc))/sqrt(pdisc-pdiff^2);
    x2=(-pdiff*sqrt(n)-quantile(ZDIST, 1-alpha/2)*sqrt(pdisc))/sqrt(pdisc-pdiff^2);
    return cdf(ZDIST, x1)+cdf(ZDIST, x2)
end
#-------------------------------------------------------------------------------
# Power Section
# Mean
# One
function oneSampleMeanEqualityP(μ₀::Real, μ₁::Real, σ::Real, n; alpha=0.05)
    z = (μ₁-μ₀)/σ*sqrt(n)
    return cdf(ZDIST, z - quantile(ZDIST, 1-alpha/2)) + cdf(ZDIST, - z - quantile(ZDIST, 1-alpha/2))
end
function oneSampleMeanEquivalenceP(μ₀::Real, μ₁::Real, σ::Real, δ::Real, n; alpha=0.05)
    z=(abs(μ₁-μ₀)-δ)/σ*sqrt(n)
    return 2*(cdf(ZDIST,z-quantile(ZDIST,1-alpha))+cdf(ZDIST,-z-quantile(ZDIST,1-alpha)))-1
end
function oneSampleMeanNSP(μ₀::Real, μ₁::Real, σ::Real, δ::Real, n; alpha=0.05) #Non-inferiority / Superiority
    z=(μ₁-μ₀-δ)/σ*sqrt(n)
    return cdf(ZDIST, z-quantile(ZDIST, 1-alpha))+cdf(ZDIST, -z-quantile(ZDIST, 1 - alpha))
end
#Two
function twoSampleMeanEqualityP(μ₀::Real, μ₁::Real, σ::Real, n; alpha=0.05, k=1)
    z=(μ₀-μ₁)/(σ*sqrt((1+1/k)/n))
    return cdf(ZDIST, z - quantile(ZDIST, 1-alpha/2)) + cdf(ZDIST, - z - quantile(ZDIST, 1-alpha/2))
end
function twoSampleMeanEquivalenceP(μ₀::Real, μ₁::Real, σ::Real, δ::Real, n; alpha=0.05, k=1)
    z=(abs(μ₀-μ₁)-δ)/(σ*sqrt((1+1/k)/n))
    return 2*(cdf(ZDIST,z-quantile(ZDIST,1-alpha))+cdf(ZDIST,-z-quantile(ZDIST,1-alpha)))-1
end
function twoSampleMeanNSP(μ₀::Real, μ₁::Real, σ::Real, δ::Real, n; alpha=0.05, k=1) #Non-inferiority / Superiority
    z=(μ₀-μ₁-δ)/(σ*sqrt((1+1/k)/n))
    return cdf(ZDIST, z-quantile(ZDIST, 1-alpha))+cdf(ZDIST, -z-quantile(ZDIST, 1 - alpha))
end
#Compare Proportion
#One Sample
function oneProportionEqualityP(p₀::Float64, p₁::Float64, n; alpha=0.05)
    z=(p₁-p₀)/sqrt(p₁*(1-p₁)/n)
    return cdf(ZDIST, z - quantile(ZDIST, 1-alpha/2)) + cdf(ZDIST, - z - quantile(ZDIST, 1-alpha/2))
end
function oneProportionEquivalenceP(p₀::Float64, p₁::Float64, δ::Real, n; alpha=0.05)
    z=(abs(p₁-p₀)-δ)/sqrt(p₁*(1-p₁)/n)
    return 2*(cdf(ZDIST,z-quantile(ZDIST,1-alpha))+cdf(ZDIST,-z-quantile(ZDIST,1-alpha)))-1
end
function oneProportionNSP(p₀::Float64, p₁::Float64, δ::Real, n; alpha=0.05)
    z=(p₁-p₀-δ)/sqrt(p₁*(1-p₁)/n)
    return cdf(ZDIST, z-quantile(ZDIST, 1-alpha)) + cdf(ZDIST, -z-quantile(ZDIST, 1 - alpha))
end

#Two Sample
function twoProportionEqualityP(p₀::Float64, p₁::Float64, n; alpha=0.05, k=1)
    z=(p₀-p₁)/sqrt(p₀*(1-p₀)/n/k+p₁*(1-p₁)/n)
    return cdf(ZDIST, z - quantile(ZDIST, 1-alpha/2)) + cdf(ZDIST, - z - quantile(ZDIST, 1-alpha/2))
end
function twoProportionEquivalenceP(p₀::Float64, p₁::Float64, δ::Real, n; alpha=0.05, k=1)
    z=(abs(p₀-p₁)-δ)/sqrt(p₀*(1-p₀)/n/k+p₁*(1-p₁)/n)
    return 2*(cdf(ZDIST,z-quantile(ZDIST,1-alpha)) + cdf(ZDIST,-z-quantile(ZDIST,1-alpha)))-1
end
function twoProportionNSP(p₀::Float64, p₁::Float64, δ::Real, n; alpha=0.05, k=1)
    z=(p₀-p₁-δ)/sqrt(p₀*(1-p₀)/n/k+p₁*(1-p₁)/n)
    return cdf(ZDIST, z-quantile(ZDIST, 1-alpha)) + cdf(ZDIST, -z-quantile(ZDIST, 1 - alpha))
end
#OR
function orEqualityP(p₀::Float64, p₁::Float64, n; alpha=0.05, k=1)
    OR=p₀*(1-p₁)/p₁/(1-p₀)
    z=log(OR)*sqrt(n)/sqrt(1/(k*p₀*(1-p₀))+1/(p₁*(1-p₁)))
    return cdf(ZDIST, z - quantile(ZDIST, 1-alpha/2)) + cdf(ZDIST, - z - quantile(ZDIST, 1-alpha/2))
end
function orEquivalenceP(p₀::Float64, p₁::Float64, δ::Real, n; alpha=0.05, k=1)
    OR=p₀*(1-p₁)/p₁/(1-p₀)
    z=(abs(log(OR))-δ)*sqrt(n)/sqrt(1/(k*p₀*(1-p₀))+1/(p₁*(1-p₁)))
    return  2*(cdf(ZDIST,z-quantile(ZDIST,1-alpha)) + cdf(ZDIST,-z-quantile(ZDIST,1-alpha)))-1
end
function orNSP(p₀::Float64, p₁::Float64, δ::Real, n; alpha=0.05, k=1)
    OR=p₀*(1-p₁)/p₁/(1-p₀)
    z=(log(OR)-δ)*sqrt(n)/sqrt(1/(k*p₀*(1-p₀))+1/(p₁*(1-p₁)))
    return cdf(ZDIST, z-quantile(ZDIST, 1-alpha)) + cdf(ZDIST, -z-quantile(ZDIST, 1 - alpha))
end
