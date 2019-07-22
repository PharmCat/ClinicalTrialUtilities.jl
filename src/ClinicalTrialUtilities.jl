# Clinical Trial Utilities
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
using Distributions, StatsBase, Statistics, Random
using QuadGK
using DataFrames
import SpecialFunctions.lgamma

#Exceptions
struct CTUException <: Exception
    n::Int
    var::String
end

Base.showerror(io::IO, e::CTUException) = print("CTU Exception code: ", e.n, " Message: ", e.var);
const ZDIST  = Normal()
const LOG2   = log(2)
const PI2    = π * 2.0
const PI2INV = 1.0 / (π * 2.0)
const VERSION = "0.1.14"
#Exceptions

struct ConfInt
    lower::Float64
    upper::Float64
    estimate::Float64
end

function getindex(a::ConfInt, b::Int64)
    if b == 1
        return a.lower
    elseif b == 2
        return a.upper
    elseif b == 3
        return a.estimate
    else
        throw(ArgumentError("Index should be in 1:3"))
    end
end

struct NCA
    result::DataFrame
    elimination::DataFrame
    pd::DataFrame
    settings::DataFrame
    textout::String
    errorlog::String
    errors::Array
end

export CTUException, ConfInt, NCA

#Owen function calc: owensQ, owensQo, ifun1, owensTint2, owensT, tfn
include("owensq.jl")
#powerTOST calc: powerTOST, powertostint, powerTOSTOwenQ, approxPowerTOST, power1TOST, approx2PowerTOST, cv2se, designProp
include("powertost.jl")
#Sample sise and Power atomic functions
include("powersamplesize.jl")
#Main sample size and power functions: sampleSize, ctPower, beSampleN
include("samplesize.jl")
#Confidence interval calculation
include("CI.jl")
#Simulations
include("SIM.jl")
#PK
include("PK.jl")
#info function
include("info.jl")
#Descriptive statistics
include("descriptives.jl")
#Frequency
include("freque.jl")
#Export
include("Export.jl")
#Randomization
include("randomization.jl")

#Sample size
export ctSampleN, beSampleN
#Power
export ctPower, bePower
#Utils
export ci2cv, pooledCV
#Other
export descriptives, owensQ, owensT
#Mudules
export SIM, CI, PK



#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------



end # module end
