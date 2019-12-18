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
using Distributions, StatsBase, Random, Roots, QuadGK, DataFrames, RecipesBase, Plots
import SpecialFunctions
import Base.show
import Base.showerror
import Base.getindex
import Base.length
import Base.in
import Base.iterate
import StatsBase.confint
import DataFrames.DataFrame


function lgamma(x)
    return SpecialFunctions.logabsgamma(x)[1]
end

const ZDIST  = Normal()
const LOG2   = log(2)
const PI2    = π * 2.0
const PI2INV = 0.5 / π
const VERSION = "0.2.0"
#Exceptions

include("abstracttypes.jl")
#Owen function calc: owensQ, owensQo, ifun1, owensTint2, owensT, tfn
include("owensq.jl")
#powerTOST calc: powerTOST, powertostint, powerTOSTOwenQ, approxPowerTOST, power1TOST, approx2PowerTOST, cv2se, designProp
include("powertost.jl")
#Sample sise and Power atomic functions
include("powersamplesize.jl")
#Main sample size and power functions: sampleSize, ctpower, besamplen
include("samplesize.jl")
#Confidence interval calculation
include("ci.jl")
#Simulations
include("sim.jl")
#info function
include("info.jl")
#DataSET
include("dataset.jl")
#Descriptive statistics
include("descriptives.jl")
#PK
include("pk.jl")
#Frequency
include("freque.jl")
#Export
include("export.jl")
#Randomization
include("randomization.jl")
#Utilities
include("utilities.jl")
#Show
include("show.jl")

#Deprecated
include("deprecated.jl")

const CTU = ClinicalTrialUtilities
#Types
export CTU, ConfInt,
#Task
CTask,
#Sample size
ctsamplen,
besamplen,
#Power
ctpower,
bepower,
#Utils
cvfromci,
cvfromvar,
pooledcv,
#Other
descriptive,
freque,
contab,
htmlexport,
#CI
confint,
propci,
diffpropci,
rrpropci,
orpropci,
meanci,
diffmeanci,
diffcmhci,
orcmhci,
rrcmhci,
#Randomization
randomtable,
randomseq,
#Pharmacokinetics
nca!,
pkimport,
pdimport,
DataFrame,
ElimRange,
DoseTime,
LimitRule,
keldata,
# Simulation
ctsim



#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------



end # module end
