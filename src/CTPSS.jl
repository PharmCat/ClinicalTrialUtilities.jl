# Clinical Trial Power and Sample Size calculation
# Version: 0.1.4
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat
# OwensQ/PowerTOST functions rewrited from https://github.com/Detlew/PowerTOST by Detlew Labes, Helmut Schuetz, Benjamin Lang
# Licence: GNU Affero General Public License v3.0
# Reference:
# Calculation based on Chow S, Shao J, Wang H. 2008. Sample Size Calculations in Clinical Research. 2nd Ed. Chapman & Hall/CRC Biostatistics Series.
# Connor R. J. 1987. Sample size for testing differences in proportions for the paired-sample design. Biometrics 43(1):207-211. page 209.
# Owen, D B (1965) "A Special Case of a Bivariate Non-central t-Distribution" Biometrika Vol. 52, pp.437-446.
# FORTRAN code in the References and matlab code given on https://people.sc.fsu.edu/~jburkardt/m_src/asa076/asa076.html by J. Burkhardt, license GNU LGPL
# D.B. Owen "Tables for computing bivariate normal Probabilities" The Annals of Mathematical Statistics, Vol. 27 (4) Dec. 1956, pp. 1075-1090
# matlab code given on https://people.sc.fsu.edu/~jburkardt/m_src/owens/owens.html by J. Burkhardt license GNU LGPL
# If you want to check and get R code you can find some here: http://powerandsamplesize.com/Calculators/
__precompile__(true)
module CTPSS
using Distributions
using Rmath #should be rewrited
using QuadGK
#using SpecialFunctions
import SpecialFunctions.lgamma
include("OwensQ.jl")
include("PowerTOST.jl")
include("PowerSampleSize.jl")

export sampleSize
export ctPower
export owensQ
export owensT
export powerTOST
export beSampleN
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

    #main sample size function
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
                    return twoSampleMeanEquality(a, b, sd, alpha=alpha, beta=beta, k=k)
                elseif type == "ei"
                    return twoSampleMeanEquivalence(a, b, sd, diff, alpha=alpha, beta=beta, k=k)
                elseif type == "ns"
                    return twoSampleMeanNS(a, b, sd, diff, alpha=alpha, beta=beta, k=k)
                else return false end
            else return false end
        elseif param == "prop"
            if 1 < a || a < 0 || 1 < b || b < 0 return false end
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
    end #sampleSize

    #clinical trial power main function
    function ctPower(;param="mean", type="ea", group="one", alpha=0.05, diff=0, sd=0, a=0, b=0, n=0, k=1)
        if alpha >= 1 || alpha <= 0  return false end
        if (type == "ei" || type == "ns") && diff == 0 return false end
        if sd == 0 && param == "mean" return false end
        if k == 0 return false end
        if n == 0 return false end

        if param == "mean"
            if group == "one"
                if type == "ea"
                    return oneSampleMeanEqualityP(a, b, sd, n, alpha=alpha)
                elseif type == "ei"
                    return oneSampleMeanEquivalenceP(a, b, sd, diff, n, alpha=alpha)
                elseif type == "ns"
                    return oneSampleMeanNSP(a, b, sd, diff, n, alpha=alpha)
                else return false end
            elseif group == "two"
                if type == "ea"
                    return twoSampleMeanEqualityP(a, b, sd, n, alpha=alpha, k=k)
                elseif type == "ei"
                    return twoSampleMeanEquivalenceP(a, b, sd, diff, n, alpha=alpha, k=k)
                elseif type == "ns"
                    return twoSampleMeanNSP(a, b, sd, diff, n, alpha=alpha, k=k)
                else return false end
            else return false end
        elseif param == "prop"
            if 1 < a || a < 0 || 1 < b || b < 0 return false end
            if group == "one"
                if type == "ea"
                    return oneProportionEqualityP(a, b, n; alpha=alpha)
                elseif type == "ei"
                    return oneProportionEquivalenceP(a, b, diff, n; alpha=alpha)
                elseif type == "ns"
                    return oneProportionNSP(a, b, diff, n; alpha=alpha)
                else return false end
            elseif group == "two"
                if type == "ea"
                    return twoProportionEqualityP(a, b, n; alpha=alpha, k=k)
                elseif type == "ei"
                    return twoProportionEquivalenceP(a, b, diff, n; alpha=alpha, k=k)
                elseif type == "ns"
                    return twoProportionNSP(a, b, diff, n; alpha=alpha, k=k)
                else return false end
            else return false end
        elseif param == "or"
            if type == "ea"
                return orEqualityP(a, b, n; alpha=alpha, k=k)
            elseif type == "ei"
                return orEquivalenceP(a, b, diff, n; alpha=alpha, k=k, logdiff=true)
            elseif type == "ns"
                return orNSP(a, b, diff, n; alpha=alpha, k=k, logdiff=true)
            else return false end
        else return false end
    end #ctPower
#-------------------------------------------------------------------------------

    function beSampleN(;theta0=0.95, theta1=0.8, theta2=1.25, cv=0, alpha=0.05, beta=0.2, logscale=true, design="2x2", method="owenq", txt=true)
        if cv <= 0 return false end
        if alpha >= 1 || alpha <= 0 || beta >= 1 || beta <= 0 return false end

        if logscale
            if theta1 < 0 || theta2 < 0 || theta0 < 0 return false end
            ltheta1 = log(theta1)
            ltheta2 = log(theta2)
            diffm   = log(theta0)
            se      = cv2se(cv)
        else
            ltheta1 = theta1
            ltheta2 = theta2
            diffm   = theta0
            se      = cv
        end

        #values for approximate n
        td = (ltheta1 + ltheta2)/2
        rd = abs((ltheta1 - ltheta2)/2)

        if rd <= 0 return false end

        d0 = diffm - td

        n0 = twoSampleMeanEquivalence(0, d0, se, rd, alpha=alpha, beta=beta) #approximate n
        n0 = ceil(n0/2)*2
        tp = 1 - beta  #target power

        if n0 > 5000 n0 = 5000 end
        pow = powerTOST(;alpha=alpha, logscale=false, theta1=ltheta1, theta2=ltheta2, theta0=diffm, cv=se, n=n0, design=design, method=method)
        np = 0
        powp = 0
        if pow > tp
            while (pow > tp)
                np = n0
                powp = pow
                n0 = n0 - 2
                pow = powerTOST(;alpha=alpha, logscale=false, theta1=ltheta1, theta2=ltheta2, theta0=diffm, cv=se, n=n0, design=design, method=method)
                if n0 < 4 return n0, pow end
            end
            estpower = powp
            estn     = np
        elseif pow < tp
            while (pow < tp)
                np = n0
                powp = pow
                n0 = n0 + 2
                pow = powerTOST(;alpha=alpha, logscale=false, theta1=ltheta1, theta2=ltheta2, theta0=diffm, cv=se, n=n0, design=design, method=method)
                if n0 > 10000 return n0, pow end
            end
            estpower = pow
            estn     = n0
        else return n0, pow end

        return estn, estpower

    end
end # module
