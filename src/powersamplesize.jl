# Clinical Trial Utilities
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)
# Calculation based on Chow S, Shao J, Wang H. 2008. Sample Size Calculations in Clinical Research. 2nd Ed. Chapman & Hall/CRC Biostatistics Series.
# OwensQ function rewrited from https://github.com/Detlew/PowerTOST by Detlew Labes, Helmut Schuetz, Benjamin Lang
# Licence: GNU Affero General Public License v3.0

# μ₀ = Test
# μ₁ = Reference value
# p₀ = Test Proportion
# p₁ = Reference Proportion
# σ  = SD - Standart deviation
# δ  = difference μ
# α  = Alpha
# β  = Beta
# n  = Subject number
# k  = group coefficient

#Ref: Chow S, Shao J, Wang H. 2008. Sample Size Calculations in Clinical Research. 2nd Ed. Chapman & Hall/CRC Biostatistics Series.
#Sample Size
#Compare Means
#One Sample
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
function two_mean_equality(μ₀::Real, μ₁::Real, σ::Real, α::Float64, β::Float64, k::Real)::Float64
    return (1 + 1 / k)*(σ*(quantile(ZDIST, 1 - α / 2) + quantile(ZDIST, 1 - β))/(μ₀ - μ₁))^2
end
function two_mean_equivalence(μ₀::Real, μ₁::Real, σ::Real, δ::Real, α::Float64, β::Float64, k::Real)::Float64
    return (1 + 1 / k)*(σ*(quantile(ZDIST, 1 - α) + quantile(ZDIST, 1 - β / 2))/(abs(μ₀ - μ₁) - δ))^2
end
function two_mean_superiority(μ₀::Real, μ₁::Real, σ::Real, δ::Real, α::Float64, β::Float64, k::Real)::Float64 #Non-inferiority / Superiority
    return (1 + 1 / k)*(σ*(quantile(ZDIST, 1 - α) + quantile(ZDIST, 1 - β))/(μ₀ - μ₁ - δ))^2
end
#Compare Proportion
#One Sample
function one_proportion_equality(p₀::Float64, p₁::Float64, α::Float64, β::Float64)::Float64
    return p₀*(1 - p₀)*((quantile(ZDIST, 1 - α / 2)+quantile(ZDIST, 1 - β))/(p₀ - p₁))^2
end
function one_proportion_equivalence(p₀::Float64, p₁::Float64, δ::Real, α::Float64, β::Float64)::Float64
    return p₀*(1 - p₀)*((quantile(ZDIST, 1 - α)+quantile(ZDIST, 1 - β / 2))/(abs(p₀ - p₁) - δ))^2
end
function one_proportion_superiority(p₀::Float64, p₁::Float64, δ::Real, α::Float64, β::Float64)::Float64
    return p₀*(1 - p₀)*((quantile(ZDIST, 1 - α)+quantile(ZDIST, 1 - β))/(p₀ - p₁ - δ))^2
end
#Two Sample
function two_proportion_equality(p₀::Float64, p₁::Float64, α::Float64, β::Float64, k::Real)::Float64
    return (p₀*(1-p₀)/k+p₁*(1-p₁))*((quantile(ZDIST, 1-α/2)+quantile(ZDIST, 1 - β))/(p₀-p₁))^2
end
function two_proportion_equivalence(p₀::Float64, p₁::Float64, δ::Real, α::Float64, β::Float64, k::Real)::Float64
    return (p₀*(1-p₀)/k+p₁*(1-p₁))*((quantile(ZDIST, 1-α)+quantile(ZDIST, 1 - β/2))/(δ - abs(p₀-p₁)))^2
end
function two_proportion_superiority(p₀::Float64, p₁::Float64, δ::Real, α::Float64, β::Float64, k::Real)::Float64
    return (p₀*(1-p₀)/k+p₁*(1-p₁))*((quantile(ZDIST, 1-α)+quantile(ZDIST, 1 - β))/(p₀-p₁-δ))^2
end
#Odd ratio Chow S, Shao J, Wang H. 2008. Sample Size Calculations in Clinical Research. 2nd Ed. Chapman & Hall/CRC Biostatistics Series.
function or_equality(p₀::Float64, p₁::Float64, α::Float64, β::Float64, k::Real)::Float64
    OR=p₀*(1-p₁)/p₁/(1-p₀)
    return (1/k/p₀/(1-p₀)+1/p₁/(1-p₁))*((quantile(ZDIST, 1-α/2)+quantile(ZDIST, 1 - β))/log(OR))^2
end
function or_equivalence(p₀::Float64, p₁::Float64, δ::Real, α::Float64, β::Float64, k::Real)::Float64
    OR=p₀*(1-p₁)/p₁/(1-p₀)
    return (1/k/p₀/(1-p₀)+1/p₁/(1-p₁))*((quantile(ZDIST, 1-α)+quantile(ZDIST, 1 - β/2))/(log(OR)-δ))^2
end
function or_superiority(p₀::Float64, p₁::Float64, δ::Real, α::Float64, β::Float64, k::Real)::Float64
    OR=p₀*(1-p₁)/p₁/(1-p₀)
    return (1/k/p₀/(1-p₀)+1/p₁/(1-p₁))*((quantile(ZDIST, 1-α)+quantile(ZDIST, 1 - β))/(log(OR)-δ))^2
end

#Connor R. J. 1987. Sample size for testing differences in proportions for the paired-sample design. Biometrics 43(1):207-211. page 209.
function mcnm(p10, p01, α::Float64, β::Float64)::Float64
    pdisc=p10+p01
    pdiff=p10-p01
    return ((quantile(ZDIST, 1-α/2)*sqrt(pdisc)+quantile(ZDIST, 1 - β)*sqrt(pdisc-pdiff^2))/pdiff)^2
end

function mcnm_pow(p10, p01, n, α::Float64)::Float64
    pdisc=p10+p01
    pdiff=p10-p01
    x1=(pdiff*sqrt(n)-quantile(ZDIST, 1-α/2)*sqrt(pdisc))/sqrt(pdisc-pdiff^2);
    x2=(-pdiff*sqrt(n)-quantile(ZDIST, 1-α/2)*sqrt(pdisc))/sqrt(pdisc-pdiff^2);
    return cdf(ZDIST, x1)+cdf(ZDIST, x2)
end
#-------------------------------------------------------------------------------
# Power Section
# Mean
# One
function one_mean_equality_pow(μ₀::Real, μ₁::Real, σ::Real, n::Int, α::Float64)::Float64
    z = (μ₀-μ₁)/σ*sqrt(n)
    return cdf(ZDIST, z - quantile(ZDIST, 1 - α / 2)) + cdf(ZDIST, - z - quantile(ZDIST, 1 - α / 2))
end
function one_mean_equivalence_pow(μ₀::Real, μ₁::Real, σ::Real, δ::Real, n::Int, α::Float64)::Float64
    z=(abs(μ₀-μ₁)-δ)/σ*sqrt(n)
    return 2*(cdf(ZDIST,z-quantile(ZDIST,1-α))+cdf(ZDIST,-z-quantile(ZDIST,1-α)))-1
end
function one_mean_superiority_pow(μ₀::Real, μ₁::Real, σ::Real, δ::Real, n::Int, α::Float64)::Float64 #Non-inferiority / Superiority
    z=(μ₀ - μ₁ - δ) / σ*sqrt(n)
    return cdf(ZDIST, z-quantile(ZDIST, 1 - α))+cdf(ZDIST, -z-quantile(ZDIST, 1 - α))
end
#Two
function two_mean_equality_pow(μ₀::Real, μ₁::Real, σ::Real, n::Int, α::Float64, k::Real)::Float64
    z=(μ₀ - μ₁)/(σ*sqrt((1 + 1 / k)/n))
    return cdf(ZDIST, z - quantile(ZDIST, 1 - α / 2)) + cdf(ZDIST, - z - quantile(ZDIST, 1 - α / 2))
end
function two_mean_equivalence_pow(μ₀::Real, μ₁::Real, σ::Real, δ::Real, n::Int, α::Float64, k::Real)::Float64
    z=(abs(μ₀-μ₁)-δ)/(σ*sqrt((1+1/k)/n))
    return 2*(cdf(ZDIST,z-quantile(ZDIST,1-α))+cdf(ZDIST,-z-quantile(ZDIST,1-α)))-1
end
function two_mean_superiority_pow(μ₀::Real, μ₁::Real, σ::Real, δ::Real, n::Int, α::Float64, k::Real)::Float64 #Non-inferiority / Superiority
    z=(μ₀-μ₁-δ)/(σ*sqrt((1+1/k)/n))
    return cdf(ZDIST, z-quantile(ZDIST, 1-α))+cdf(ZDIST, -z-quantile(ZDIST, 1 - α))
end
#Compare Proportion
#One Sample
function one_proportion_equality_pow(p₀::Float64, p₁::Float64, n::Int, α::Float64)::Float64
    z=(p₀-p₁)/sqrt(p₀*(1-p₀)/n)
    return cdf(ZDIST, z - quantile(ZDIST, 1-α/2)) + cdf(ZDIST, - z - quantile(ZDIST, 1-α/2))
end
function oneProportionEquivalenceP(p₀::Float64, p₁::Float64, δ::Real, n::Int, α::Float64)::Float64
    z=(abs(p₀-p₁)-δ)/sqrt(p₀*(1-p₀)/n)
    return 2*(cdf(ZDIST,z-quantile(ZDIST,1-α))+cdf(ZDIST,-z-quantile(ZDIST,1-α)))-1
end
function one_proportion_superiority_pow(p₀::Float64, p₁::Float64, δ::Real, n, α::Float64)::Float64
    z=(p₀-p₁-δ)/sqrt(p₀*(1-p₀)/n)
    return cdf(ZDIST, z-quantile(ZDIST, 1-α)) + cdf(ZDIST, -z-quantile(ZDIST, 1 - α))
end

#Two Sample
function two_proportion_equality_pow(p₀::Float64, p₁::Float64, n::Int, α::Float64, k::Real)::Float64
    z=(p₀-p₁)/sqrt(p₀*(1-p₀)/n/k+p₁*(1-p₁)/n)
    return cdf(ZDIST, z - quantile(ZDIST, 1-α/2)) + cdf(ZDIST, - z - quantile(ZDIST, 1-α/2))
end
function two_proportion_equivalence_pow(p₀::Float64, p₁::Float64, δ::Real, n::Int, α::Float64, k::Real)::Float64
    z=(abs(p₀-p₁)-δ)/sqrt(p₀*(1-p₀)/n/k+p₁*(1-p₁)/n)
    return 2*(cdf(ZDIST,z-quantile(ZDIST,1-α)) + cdf(ZDIST,-z-quantile(ZDIST,1-α)))-1
end
function two_proportion_superiority_pow(p₀::Float64, p₁::Float64, δ::Real, n::Int, α::Float64, k::Real)::Float64
    z=(p₀-p₁-δ)/sqrt(p₀*(1-p₀)/n/k+p₁*(1-p₁)/n)
    return cdf(ZDIST, z-quantile(ZDIST, 1-α)) + cdf(ZDIST, -z-quantile(ZDIST, 1 - α))
end
#OR
function or_equality_pow(p₀::Float64, p₁::Float64, n::Int, α::Float64, k::Real)::Float64
    OR=p₀*(1-p₁)/p₁/(1-p₀)
    z=log(OR)*sqrt(n)/sqrt(1/(k*p₀*(1-p₀))+1/(p₁*(1-p₁)))
    return cdf(ZDIST, z - quantile(ZDIST, 1-α/2)) + cdf(ZDIST, - z - quantile(ZDIST, 1-α/2))
end
function or_equivalence_pow(p₀::Float64, p₁::Float64, δ::Real, n::Int, α::Float64, k::Real)::Float64
    OR=p₀*(1-p₁)/p₁/(1-p₀)
    z=(abs(log(OR))-δ)*sqrt(n)/sqrt(1/(k*p₀*(1-p₀))+1/(p₁*(1-p₁)))
    return  2*(cdf(ZDIST,z-quantile(ZDIST,1-α)) + cdf(ZDIST,-z-quantile(ZDIST,1-α)))-1
end
function or_superiority_pow(p₀::Float64, p₁::Float64, δ::Real, n::Int, α::Float64, k::Real)::Float64
    OR=p₀*(1-p₁)/p₁/(1-p₀)
    z=(log(OR)-δ)*sqrt(n)/sqrt(1/(k*p₀*(1-p₀))+1/(p₁*(1-p₁)))
    return cdf(ZDIST, z-quantile(ZDIST, 1-α)) + cdf(ZDIST, -z-quantile(ZDIST, 1 - α))
end
