# Clinical Trial Utilities
# Version: 0.1.6
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)
# OwensQ/PowerTOST functions rewrited from https://github.com/Detlew/PowerTOST by Detlew Labes, Helmut Schuetz, Benjamin Lang
# Licence: GNU Affero General Public License v3.0
# Reference:
# Calculation based on Chow S, Shao J, Wang H. 2008. Sample Size Calculations in Clinical Research. 2nd Ed. Chapman & Hall/CRC Biostatistics Series.
# Connor R. J. 1987. Sample size for testing differences in proportions for the paired-sample design. Biometrics 43(1):207-211. page 209.
# Owen, D B (1965) "A Special Case of a Bivariate Non-central t-Distribution" Biometrika Vol. 52, pp.437-446.
# FORTRAN code by J. Burkhardt, license GNU LGPL
# D.B. Owen "Tables for computing bivariate normal Probabilities" The Annals of Mathematical Statistics, Vol. 27 (4) Dec. 1956, pp. 1075-1090
# matlab code  by J. Burkhardt license GNU LGPL
# If you want to check and get R code you can find some here: http://powerandsamplesize.com/Calculators/
__precompile__(true)
module ClinicalTrialUtilities
using Distributions
#using Rmath #should be rewrited
using QuadGK
#using SpecialFunctions
import SpecialFunctions.lgamma
#Owen function calc: owensQ, owensQo, ifun1, owensTint2, owensT, tfn
include("OwensQ.jl")
#powerTOST calc: powerTOST, powerTOSTint, powerTOSTOwenQ, approxPowerTOST, power1TOST, approx2PowerTOST, cv2se, designProp
include("PowerTOST.jl")
#Sample sise and Power atomic functions
include("PowerSampleSize.jl")
#Main sample size and power functions: sampleSize, ctPower, beSampleN
include("SampleSize.jl")
#Confidence interval calculation
include("CI.jl")
#info function
include("Info.jl")

#Sample size
export sampleSize
export beSampleN
#Power
export ctPower
export powerTOST
#Other
export owensQ
export owensT
#Structs - should be rewrited or deleted
export ParamSet
export sampleSizeParam
#Exceptions
export CTUException
#Confidence Intervals for Proportions and Means
export CI

#Exceptions
struct CTUException <: Exception
    n::Int
    var::String
end
Base.showerror(io::IO, e::CTUException) = print("CTU Exception code: ", e.n, " Message: ", e.var)

    const ZDIST = Normal()
#-------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------
end # module end
