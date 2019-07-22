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
# m0 = μ0; m1 = μ

function oneSampleMeanEquality(m0, m1, sd; alpha=0.05, beta=0.2)::Float64
    return ((quantile(ZDIST, 1-alpha/2) + quantile(ZDIST, 1-beta))*sd/(m1-m0))^2
end
function oneSampleMeanEquivalence(m0, m1, sd, diff; alpha=0.05, beta=0.2)
    return (sd*(quantile(ZDIST, 1-alpha) + quantile(ZDIST, 1 - beta/2))/(diff-abs(m1-m0)))^2
end
function oneSampleMeanNS(m0, m1, sd, diff; alpha=0.05, beta=0.2) #Non-inferiority / Superiority
    return (sd*(quantile(ZDIST, 1-alpha) + quantile(ZDIST, 1 - beta))/(m1 - m0 - diff))^2
end
#Two Sample
# m0 = μA - Group A; m1 = μB - Group B
function twoSampleMeanEquality(m0, m1, sd; alpha=0.05, beta=0.2, k=1)
    return (1+1/k)*((quantile(ZDIST, 1-alpha/2) + quantile(ZDIST, 1-beta))*sd/(m0-m1))^2
end
function twoSampleMeanEquivalence(m0, m1, sd, diff; alpha=0.05, beta=0.2, k=1)
    return (1+1/k)*(sd*(quantile(ZDIST, 1-alpha) + quantile(ZDIST, 1 - beta/2))/(abs(m0-m1)-diff))^2
end
function twoSampleMeanNS(m0, m1, sd, diff; alpha=0.05, beta=0.2, k=1) #Non-inferiority / Superiority
    return (1+1/k)*(sd*(quantile(ZDIST, 1-alpha) + quantile(ZDIST, 1 - beta))/(m0 - m1 - diff))^2
end
#Compare Proportion
#One Sample
function oneProportionEquality(p0, p1; alpha=0.05, beta=0.2)
    return p1*(1-p1)*((quantile(ZDIST, 1-alpha/2)+quantile(ZDIST, 1 - beta))/(p1-p0))^2
end
function oneProportionEquivalence(p0, p1, diff; alpha=0.05, beta=0.2)
    return p1*(1-p1)*((quantile(ZDIST, 1-alpha)+quantile(ZDIST, 1 - beta/2))/(abs(p1-p0)-diff))^2
end
function oneProportionNS(p0, p1, diff; alpha=0.05, beta=0.2)
    return p1*(1-p1)*((quantile(ZDIST, 1-alpha)+quantile(ZDIST, 1 - beta))/(p1-p0-diff))^2
end
#Two Sample
function twoProportionEquality(p0, p1; alpha=0.05, beta=0.2, k=1)
    return (p0*(1-p0)/k+p1*(1-p1))*((quantile(ZDIST, 1-alpha/2)+quantile(ZDIST, 1 - beta))/(p0-p1))^2
end
function twoProportionEquivalence(p0, p1, diff; alpha=0.05, beta=0.2, k=1)
    return (p0*(1-p0)/k+p1*(1-p1))*((quantile(ZDIST, 1-alpha)+quantile(ZDIST, 1 - beta/2))/(abs(p0-p1)-diff))^2
end
function twoProportionNS(p0, p1, diff; alpha=0.05, beta=0.2, k=1)
    return (p0*(1-p0)/k+p1*(1-p1))*((quantile(ZDIST, 1-alpha)+quantile(ZDIST, 1 - beta))/(p0-p1-diff))^2
end
#Odd ratio Chow S, Shao J, Wang H. 2008. Sample Size Calculations in Clinical Research. 2nd Ed. Chapman & Hall/CRC Biostatistics Series.
function orEquality(p0, p1; alpha=0.05, beta=0.2, k=1)
    OR=p0*(1-p1)/p1/(1-p0)
    return (1/k/p0/(1-p0)+1/p1/(1-p1))*((quantile(ZDIST, 1-alpha/2)+quantile(ZDIST, 1 - beta))/log(OR))^2
end
function orEquivalence(p0, p1, diff; alpha=0.05, beta=0.2, k=1, logdiff=false)
    if !logdiff diff=log(diff) end
    OR=p0*(1-p1)/p1/(1-p0)
    return (1/k/p0/(1-p0)+1/p1/(1-p1))*((quantile(ZDIST, 1-alpha)+quantile(ZDIST, 1 - beta/2))/(log(OR)-diff))^2
end
function orNS(p0, p1, diff; alpha=0.05, beta=0.2, k=1, logdiff=false)
    if !logdiff diff=log(diff) end
    OR=p0*(1-p1)/p1/(1-p0)
    return (1/k/p0/(1-p0)+1/p1/(1-p1))*((quantile(ZDIST, 1-alpha)+quantile(ZDIST, 1 - beta))/(log(OR)-diff))^2
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
function oneSampleMeanEqualityP(m0,m1,sd, n; alpha=0.05)
    z = (m1-m0)/sd*sqrt(n)
    return cdf(ZDIST, z - quantile(ZDIST, 1-alpha/2)) + cdf(ZDIST, - z - quantile(ZDIST, 1-alpha/2))
end
function oneSampleMeanEquivalenceP(m0, m1, sd, diff, n; alpha=0.05)
    z=(abs(m1-m0)-diff)/sd*sqrt(n)
    return 2*(cdf(ZDIST,z-quantile(ZDIST,1-alpha))+cdf(ZDIST,-z-quantile(ZDIST,1-alpha)))-1
end
function oneSampleMeanNSP(m0, m1, sd, diff, n; alpha=0.05) #Non-inferiority / Superiority
    z=(m1-m0-diff)/sd*sqrt(n)
    return cdf(ZDIST, z-quantile(ZDIST, 1-alpha))+cdf(ZDIST, -z-quantile(ZDIST, 1 - alpha))
end
#Two
function twoSampleMeanEqualityP(m0, m1, sd, n; alpha=0.05, k=1)
    z=(m0-m1)/(sd*sqrt((1+1/k)/n))
    return cdf(ZDIST, z - quantile(ZDIST, 1-alpha/2)) + cdf(ZDIST, - z - quantile(ZDIST, 1-alpha/2))
end
function twoSampleMeanEquivalenceP(m0, m1, sd, diff, n; alpha=0.05, k=1)
    z=(abs(m0-m1)-diff)/(sd*sqrt((1+1/k)/n))
    return 2*(cdf(ZDIST,z-quantile(ZDIST,1-alpha))+cdf(ZDIST,-z-quantile(ZDIST,1-alpha)))-1
end
function twoSampleMeanNSP(m0, m1, sd, diff, n; alpha=0.05, k=1) #Non-inferiority / Superiority
    z=(m0-m1-diff)/(sd*sqrt((1+1/k)/n))
    return cdf(ZDIST, z-quantile(ZDIST, 1-alpha))+cdf(ZDIST, -z-quantile(ZDIST, 1 - alpha))
end
#Compare Proportion
#One Sample
function oneProportionEqualityP(p0, p1, n; alpha=0.05)
    z=(p1-p0)/sqrt(p1*(1-p1)/n)
    return cdf(ZDIST, z - quantile(ZDIST, 1-alpha/2)) + cdf(ZDIST, - z - quantile(ZDIST, 1-alpha/2))
end
function oneProportionEquivalenceP(p0, p1, diff, n; alpha=0.05)
    z=(abs(p1-p0)-diff)/sqrt(p1*(1-p1)/n)
    return 2*(cdf(ZDIST,z-quantile(ZDIST,1-alpha))+cdf(ZDIST,-z-quantile(ZDIST,1-alpha)))-1
end
function oneProportionNSP(p0, p1, diff, n; alpha=0.05)
    z=(p1-p0-diff)/sqrt(p1*(1-p1)/n)
    return cdf(ZDIST, z-quantile(ZDIST, 1-alpha)) + cdf(ZDIST, -z-quantile(ZDIST, 1 - alpha))
end

#Two Sample
function twoProportionEqualityP(p0, p1, n; alpha=0.05, k=1)
    z=(p0-p1)/sqrt(p0*(1-p0)/n/k+p1*(1-p1)/n)
    return cdf(ZDIST, z - quantile(ZDIST, 1-alpha/2)) + cdf(ZDIST, - z - quantile(ZDIST, 1-alpha/2))
end
function twoProportionEquivalenceP(p0, p1, diff, n; alpha=0.05, k=1)
    z=(abs(p0-p1)-diff)/sqrt(p0*(1-p0)/n/k+p1*(1-p1)/n)
    return 2*(cdf(ZDIST,z-quantile(ZDIST,1-alpha)) + cdf(ZDIST,-z-quantile(ZDIST,1-alpha)))-1
end
function twoProportionNSP(p0, p1, diff, n; alpha=0.05, k=1)
    z=(p0-p1-diff)/sqrt(p0*(1-p0)/n/k+p1*(1-p1)/n)
    return cdf(ZDIST, z-quantile(ZDIST, 1-alpha)) + cdf(ZDIST, -z-quantile(ZDIST, 1 - alpha))
end
#OR
function orEqualityP(p0, p1, n; alpha=0.05, k=1)
    OR=p0*(1-p1)/p1/(1-p0)
    z=log(OR)*sqrt(n)/sqrt(1/(k*p0*(1-p0))+1/(p1*(1-p1)))
    return cdf(ZDIST, z - quantile(ZDIST, 1-alpha/2)) + cdf(ZDIST, - z - quantile(ZDIST, 1-alpha/2))
end
function orEquivalenceP(p0, p1, diff, n; alpha=0.05, k=1, logdiff=true)
    if !logdiff diff=log(diff) end
    OR=p0*(1-p1)/p1/(1-p0)
    z=(abs(log(OR))-diff)*sqrt(n)/sqrt(1/(k*p0*(1-p0))+1/(p1*(1-p1)))
    return  2*(cdf(ZDIST,z-quantile(ZDIST,1-alpha)) + cdf(ZDIST,-z-quantile(ZDIST,1-alpha)))-1
end
function orNSP(p0, p1, diff, n; alpha=0.05, k=1, logdiff=true)
    if !logdiff diff=log(diff) end
    OR=p0*(1-p1)/p1/(1-p0)
    z=(log(OR)-diff)*sqrt(n)/sqrt(1/(k*p0*(1-p0))+1/(p1*(1-p1)))
    return cdf(ZDIST, z-quantile(ZDIST, 1-alpha)) + cdf(ZDIST, -z-quantile(ZDIST, 1 - alpha))
end
