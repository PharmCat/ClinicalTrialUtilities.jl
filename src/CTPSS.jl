# Clinical Trial Power and Sample Size calculation
# Version: 0.1.1
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat
# Calculation based on Chow S, Shao J, Wang H. 2008. Sample Size Calculations in Clinical Research. 2nd Ed. Chapman & Hall/CRC Biostatistics Series.
# OwensQ function rewrited from https://github.com/Detlew/PowerTOST by Detlew Labes, Helmut Schuetz, Benjamin Lang
# Licence: GNU Affero General Public License v3.0
__precompile__(true)
module CTPSS
using Distributions
using Rmath #should be rewrited
using QuadGK
#using SpecialFunctions
import SpecialFunctions.lgamma
include("OwensQ.jl")
include("PowerTOST.jl")

export sampleSize
export owensQ
export ParamSet
export sampleSizeParam

    const ZDIST = Normal()

    mutable struct ParamSet
        param::String
        type::String
        group::String
        alpha::Float32
        beta::Float32
        sd::Float32
        a::Float32
        b::Float32
        k::Float32
    end
    function sampleSizeParam(x::ParamSet)
        return sampleSize(param=x.param, type=x.type, group=x.group, alpha=x.alpha, beta=x.beta, sd=x.sd, a=x.a, b=x.b, k=x.k)
    end

    function sampleSize(;param="mean", type="ea", group="one", alpha=0.05, beta=0.2, diff=0, sd=0, a=0, b=0, k=1)
        if alpha >= 1 || alpha <= 0 || beta >= 1 || beta <= 0 return false end
        if (type == "ei" || type == "ns") && diff == 0 return false end
        if sd == 0 && param == "mean" return false end
        if k == 0 return false end
        if param == "mean"
            if group == "one"
                if type == "ea"
                    return oneSampleMeanEquality(a, b, sd, alpha=alpha, beta=beta)
                elseif type == "ei"
                    return oneSampleMeanEquivalence(a, b, sd, diff, alpha=alpha, beta=beta)
                elseif type == "ns"
                    return oneSampleMeanNS(a, b, sd, diff, alpha=alpha, beta=beta)
                else return false end
            elseif group == "two"
                if type == "ea"
                    return twoSampleMeanEquality(m0, m1, sd, alpha=alpha, beta=beta, k=k)
                elseif type == "ei"
                    return twoSampleMeanEquivalence(m0, m1, sd, diff, alpha=alpha, beta=beta, k=k)
                elseif type == "ns"
                    return twoSampleMeanNS(m0, m1, sd, diff, alpha=alpha, beta=beta, k=k)
                else return false end
            else return false end
        elseif param == "prop"
            if 1 > a < 0 || 1 > b < 0 return false end
            if group == "one"
                if type == "ea"
                    return oneProportionEquality(a, b; alpha=alpha, beta=beta)
                elseif type == "ei"
                    return oneProportionEquivalence(a, b, diff; alpha=alpha, beta=beta)
                elseif type == "ns"
                    return oneProportionNS(a, b, diff; alpha=alpha, beta=beta)
                else return false end
            elseif group == "two"
                if type == "ea"
                    return twoProportionEquality(a, b; alpha=alpha, beta=beta, k=k)
                elseif type == "ei"
                    return twoProportionEquivalence(a, b, diff; alpha=alpha, beta=beta, k=k)
                elseif type == "ns"
                    return twoProportionNS(a, b, diff; alpha=alpha, beta=beta, k=k)
                else return false end
            else return false end
        elseif param == "or"
            if type == "ea"
                return orEquality(a, b; alpha=alpha, beta=beta, k=k)
            elseif type == "ei"
                return orEquivalence(a, b, diff; alpha=alpha, beta=beta, k=k, logdiff=true)
            elseif type == "ns"
                return orNS(a, b, diff; alpha=alpha, beta=beta, k=k, logdiff=true)
            else return false end
        else return false end
    end

    #Ref: Chow S, Shao J, Wang H. 2008. Sample Size Calculations in Clinical Research. 2nd Ed. Chapman & Hall/CRC Biostatistics Series.

    #Sample Size
    #Compare Means
    #One Sample
    # m0 = μ0; m1 = μ
    function oneSampleMeanEquality(m0, m1, sd; alpha=0.05, beta=0.2)
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
        return (p1*(1-p1)/k+p0*(1-p0))*((quantile(ZDIST, 1-alpha/2)+quantile(ZDIST, 1 - beta))/(p1-p0))^2
    end

    function twoProportionEquivalence(p0, p1, diff; alpha=0.05, beta=0.2, k=1)
        return (p1*(1-p1)/k+p0*(1-p0))*((quantile(ZDIST, 1-alpha)+quantile(ZDIST, 1 - beta/2))/(abs(p1-p0)-diff))^2
    end

    function twoProportionNS(p0, p1, diff; alpha=0.05, beta=0.2, k=1)
        return (p1*(1-p1)/k+p0*(1-p0))*((quantile(ZDIST, 1-alpha)+quantile(ZDIST, 1 - beta))/(p1-p0-diff))^2
    end

    function orEquality(p0, p1; alpha=0.05, beta=0.2, k=1)
        OR=p1*(1-p0)/p0/(1-p1)
        return (1/k/p1/(1-p1)+1/p0/(1-p0))*((quantile(ZDIST, 1-alpha/2)+quantile(ZDIST, 1 - beta))/log(OR))^2
    end

    function orEquivalence(p0, p1, diff; alpha=0.05, beta=0.2, k=1, logdiff=true)
        if !logdiff diff=log(diff) end
        OR=p1*(1-p0)/p0/(1-p1)
        return (1/k/p1/(1-p1)+1/p0/(1-p0))*((quantile(ZDIST, 1-alpha)+quantile(ZDIST, 1 - beta/2))/(log(OR)-diff))^2
    end

    function orNS(p0, p1, diff; alpha=0.05, beta=0.2, k=1, logdiff=true)
        if !logdiff diff=log(diff) end
        OR=p1*(1-p0)/p0/(1-p1)
        return (1/k/p1/(1-p1)+1/p0/(1-p0))*((quantile(ZDIST, 1-alpha)+quantile(ZDIST, 1 - beta))/(log(OR)-diff))^2
    end
#p01=0.45
#p10=0.05
#alpha=0.05*2 (one side)
#beta=0.10
#pdisc=p10+p01
#pdiff=p10-p01
#(n=((qnorm(1-alpha/2)*sqrt(pdisc)+qnorm(1-beta)*sqrt(pdisc-pdiff^2))/pdiff)^2)
#22.80907
#CTPSS.mcnm(0.45, 0.05, alpha=0.1, beta=0.1)
#22.80907052752994
    function mcnm(p10, p01; alpha=0.05, beta=0.2)
        pdisc=p10+p01
        pdiff=p10-p01
        return ((quantile(ZDIST, 1-alpha/2)*sqrt(pdisc)+quantile(ZDIST, 1 - beta)*sqrt(pdisc-pdiff^2))/pdiff)^2
    end
#-------------------------------------------------------------------------------

end # module
