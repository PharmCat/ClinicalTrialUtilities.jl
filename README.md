# ClinicalTrialUtilities

 Clinical trial related calculation: descriptive statistics, power and sample size calculation, power simulations, confidence interval, pharmacokinetics/pharmacodynamics parameters calculation.

[![Build Status](https://travis-ci.com/PharmCat/ClinicalTrialUtilities.jl.svg?branch=master)](https://travis-ci.com/PharmCat/ClinicalTrialUtilities.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/35f8b5vq259sbssg?svg=true)](https://ci.appveyor.com/project/PharmCat/clinicaltrialutilities-jl)
[![codecov](https://codecov.io/gh/PharmCat/ClinicalTrialUtilities.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/PharmCat/ClinicalTrialUtilities.jl)
[![Coverage Status](https://coveralls.io/repos/github/PharmCat/ClinicalTrialUtilities.jl/badge.svg?branch=master)](https://coveralls.io/github/PharmCat/ClinicalTrialUtilities.jl?branch=master)

## Description

The package is designed to perform calculations related to the planning and analysis of the results of clinical trials. The package includes the basic functions described below, as well as a few modules to perform specific calculations.

## Content

- [Installation](#Installation)
- [Basic functions](#Basic)
- [Usage](#Usage)
  - [descriptives](#descriptives)
  - [ctSampleN](#ctSampleN)
  - [ctPower](#ctPower)
  - [besamplen](#besamplen)
  - [bePower](#bePower)
  - [ci2cv](#ci2cv)
  - [pooledCV](#pooledCV)
- [Submodules](#Submodules)
  - [CI](https://github.com/PharmCat/ClinicalTrialUtilities.jl/blob/master/doc/CI.md)
    - [oneProp](https://github.com/PharmCat/ClinicalTrialUtilities.jl/blob/master/doc/CI.md#oneProp)
    - [oneMean](https://github.com/PharmCat/ClinicalTrialUtilities.jl/blob/master/doc/CI.md#oneMean)
    - [twoProp](https://github.com/PharmCat/ClinicalTrialUtilities.jl/blob/master/doc/CI.md#twoProp)
    - [twoMeans](https://github.com/PharmCat/ClinicalTrialUtilities.jl/blob/master/doc/CI.md#twoMeans)
    - [cmh](https://github.com/PharmCat/ClinicalTrialUtilities.jl/blob/master/doc/CI.md#cmh)
  - [PK](https://github.com/PharmCat/ClinicalTrialUtilities.jl/blob/master/doc/PK.md)
    - [nca](https://github.com/PharmCat/ClinicalTrialUtilities.jl/blob/master/doc/PK.md#nca)
  - [SIM](https://github.com/PharmCat/ClinicalTrialUtilities.jl/blob/master/doc/SIM.md)
    - [bePower](https://github.com/PharmCat/ClinicalTrialUtilities.jl/blob/master/doc/SIM.md#bePower)
    - [ctPropPower](https://github.com/PharmCat/ClinicalTrialUtilities.jl/blob/master/doc/SIM.md#ctPropPower)
    - [ctPropSampleN](https://github.com/PharmCat/ClinicalTrialUtilities.jl/blob/master/doc/SIM.md#ctPropSampleN)
    - [ctMeansPower](https://github.com/PharmCat/ClinicalTrialUtilities.jl/blob/master/doc/SIM.md#ctMeansPower)
    - [ctMeansPowerFS](https://github.com/PharmCat/ClinicalTrialUtilities.jl/blob/master/doc/SIM.md#ctMeansPowerFS)
- [Types](#Types)
- [Examples](#Examples)
- [Other](#Other)
- [Dependencies](#Dependencies)
- [Copyrights&References](#Copyrights&References)

### <a name="Installation">Installation</a>
```
using Pkg; Pkg.add("ClinicalTrialUtilities");
```
or
```
using Pkg; Pkg.clone("https://github.com/PharmCat/ClinicalTrialUtilities.jl.git");
```

And then to perform tests:

```
Pkg.test("ClinicalTrialUtilities");
```

### <a name="Basic">Basic functions</a>

- [Descriptive statistics](#descriptives)

- [Clinical trial sample size estimation](#ctSampleN)

- [Clinical trial power estimation](#ctPower)

- [Iterative sample size estimation for bioequivalence trials](#besamplen)

- [Power estimation for bioequivalence trials](#bePower)

- [CV from CI for bioequivalence trials](#ci2cv)

- [Pooled CV from multiple sources](#pooledCV)

### <a name="Usage">Usage</a>

**NB! Hypothesis types:**

- :ea - Equality: two-sided;
- :ei - Equivalencens: two one-sided hypothesis;
- :ns - Non-Inferiority / Superiority: one-sided hypothesis, for some cases you should use two-sided hypothesis for  Non-Inferiority/Superiority, you can use alpha/2 for this;

### <a name="descriptives">descriptives</a>

Descriptive statistics.

```
descriptives(data::DataFrame; sort = NaN, vars = NaN, stats = [:n, :mean, :sd, :sem, :uq, :median, :lq])::DataFrame
```

### <a name="ctSampleN">ctSampleN</a>

Sample size estimation for clinical trial.

```
ctsamplen(;param=:notdef, type=:notdef, group=:notdef, alpha=0.05, beta=0.2, diff=0, sd=0, a=0, b=0, k=1, logdiff=false, out=:num)
```

**param (Parameter type):**
- :mean - Means;
- :prop - Proportions;
- :or - Odd Ratio;

**type (Hypothesis type):**
- :ea - Equality;
- :ei - Equivalencens;
- :ns - Non-Inferiority / Superiority (!one-sided hypothesis!);
- :mcnm - McNemar's Equality test;

**group (group num):**
- :one - One sample;
- :two - Two sample, result is for one group, second group size = n * k;

**alpha** - Alpha (o < alpha < 1)  (default=0.05);

**beta** - Beta (o < beta < 1) (default=0.2); power = 1 - beta;

**diff** - difference/equivalence margin/non-inferiority/superiority margin;

**sd** - Standard deviation (σ, for Means only);

**a** - Null Hypothesis mean (μ0), Group A;

**b** - True mean (μ) for one sample / Group B for two sample design;

**k** - Na/Nb (after sample size estimation second group size: Na=k*Nb, only for two sample design) (default=1);

**logdiff** - diff is log transformed for OR:
- false (default, diff would be transformed)
- true

**out** - output type:
- :num   - numeric (default);
- :str   - String variable with text output;
- :vstr  - numeric and String variable;
- :print - print to console;

### <a name="ctPower">ctPower</a>

Power estimation for clinical trials.

```
ctPower(;param=:notdef, type=:notdef, group=:notdef, alpha=0.05, logdiff=false, diff=0, sd=0, a=0, b=0, n=0, k=1,  out=:num)
```

**param (Parameter type):**
- :mean - Means;
- :prop - Proportions;
- :or   - Odd Ratio;

**type (Hypothesis type):**
- :ea   - Equality;
- :ei   - Equivalence;
- :ns   - Non-Inferiority / Superiority;
- :mcnm - McNemar's Equality test;

**group (group num):**
- :one - one sample;
- :two - Two sample;

**alpha** - Alpha (0<alpha<1)  (default=0.05);

**n** - Subjects number;

**diff** - difference/equivalence margin/non-inferiority/superiority margin;

**sd** - Standard deviation (σ, for Means only);

**a** - Null Hypothesis mean (μ0), Group A;

**b** - True mean (μ) for one sample / Group B for two sample design;

**k** - Na/Nb (after sample size estimation second group size: Na=k*Nb, only for two sample design) (default=1);

**logdiff** - diff is log transformed for OR:
- false (default, diff would be transformed)
- true

**out** - output type:
- :num   - numeric (default);
- :str   - String variable with text output;
- :vstr  - numeric and String variable;
- :print - print to console;

### <a name="besamplen">besamplen</a>

Sample size estimation for bioequivalence study (iterative procedure).

```
besamplen(;alpha=0.05, beta=0.2, theta0=0.95, theta1=0.8, theta2=1.25, cv=0.0, logscale=true, design=:d2x2, method=:owenq,  out=:num)
```

**alpha** - Alpha (o < alpha < 1)  (default=0.05);

**beta** - Beta (o < beta < 1) (default=0.2); power = 1 - beta;

**theta0** - T/R Ratio (default=0.95);

**theta1** - Lower Level (default=0.8);

**theta2** - Upper level (default=1.25);

**cv** - coefficient of variation;

**logscale** - theta1, theta2, theta0: if true - make log transformation (default true);

**design** - trial design;
- :parallel
- :d2x2 (default)
- :d2x2x4
- :d2x4x4
- :d2x3x3
- :d2x4x2
- :d3x3
- :d3x6x3

**method** - calculating method: Owen's Q Function | NonCentral T | Shifted;
- :owenq (default)
- :nct
- :shifted

**out** - output type:
- :num   - numeric (default);
- :str   - String variable with text output;
- :vstr  - numeric and String variable;
- :print - print to console;

### <a name="bePower">bePower</a>

Power estimation for bioequivalence trials.

```
bePower(;alpha=0.05, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.0, n=0, logscale=true, design=:d2x2, method=:owenq,  out=:num)
```

**alpha** - Alpha (0 < alpha < 1)  (default=0.05);

**theta1** - Lower Level (default=0.8);

**theta2** - Upper level (default=1.25);

**theta0** - T/R Ratio (default=0.95);

**cv** - coefficient of variation;

**n** - subject number;

**logscale** - theta1, theta2, theta0: if true - make log transformation (default true);

**design** - trial design;
- :parallel
- :d2x2 (default)
- :d2x2x4
- :d2x4x4
- :d2x3x3
- :d2x4x2
- :d3x3
- :d3x6x3

**method** - calculating method: Owen's Q Function | NonCentral T, Shifted;
- :owenq (default)
- :nct
- :shifted

**out** - output type:
- :num   - numeric (default);
- :str   - String variable with text output;
- :vstr  - numeric and String variable;
- :print - print to console;

### <a name="ci2cv">ci2cv</a>

Take CV from known CI and subject number.

```
ci2cv(;alpha = 0.05, theta1 = 0.8, theta2 = 1.25, n, design=:d2x2, mso=false, cvms=false)
```

**alpha** - Alpha (o < alpha < 1)  (default=0.05);

**beta** - Beta (o < beta < 1) (default=0.2); power = 1 - beta;

**theta1** - Lower Level (default=0.8);

**theta2** - Upper level (default=1.25);

**n** - subject n;

**design** - trial design;
- :parallel
- :d2x2 (default)
- :d2x2x4
- :d2x4x4
- :d2x3x3
- :d2x4x2
- :d3x3
- :d3x6x3

**mso**

Calculate MS only
- false(default)
- true

**cvms**

Calculate CV and MS
- false(default)
- true

### <a name="pooledCV">pooledCV</a>

Get pooled CV from multiple sources.

```
pooledCV(data::DataFrame; cv=:cv, df=:df, alpha=0.05, returncv=true)::ConfInt
```

**data**::DataFrame - Dataframe with CV data

**cv**::Symbol - CV column in dataframe

**df**::Symbol - DF column in dataframe

**alpha** - Alpha for var/cv confidence interval.

**returncv** - Return CV or var:

- true  - return cv
- false - return var

### <a name="Submodules">Submodules</a>

* Confidence interval calculation - [Doc](https://github.com/PharmCat/ClinicalTrialUtilities.jl/blob/master/doc/CI.md)
  * oneProp
  * oneMeans
  * twoProp
  * twoMeans
  * cmh

* Pharmacokinetics calculation - [Doc](https://github.com/PharmCat/ClinicalTrialUtilities.jl/blob/mater/doc/PK.md)
  * nca

* Simulations - [Doc](https://github.com/PharmCat/ClinicalTrialUtilities.jl/blob/master/doc/SIM.md)
  * bePower
  * bePowerSIM
  * ctPropPower
  * ctPropSampleN
  * ctMeansPower
  * ctMeansPowerFS

### <a name="Types">Types</a>

Confidence intervals

```
struct ConfInt
    lower::Float64
    upper::Float64
    estimate::Float64
end
```

Pharmacokinetics noncompartment analysis output

```
struct NCA
    result::DataFrame
    elimination::DataFrame
    settings::DataFrame
    textout::String
    errorlog::String
    errors::Array
end
```

### <a name="Examples">Examples</a>

```
#Sample size for one proportion equality
ctsamplen(param=:prop, type=:ea, group=:one, a=0.3, b=0.5)
#Equivalence for two means
ctsamplen(param=:mean, type=:ei, group=:two, diff=0.3, sd=1, a=0.3, b=0.5)
#Odd ratio non-inferiority
ctsamplen(param=:or, type=:ns, diff=-0.1, a=0.3, b=0.5, k=2)
#Odd ratio equality
ctsamplen(param=:or, type=:ea, a=0.3, b=0.5, k=2)

#Power
ctPower(param=:mean, type=:ea, group=:one, a=1.5, b=2, sd=1,n=32, alpha=0.05)

#Bioequivalence sample size
besamplen(alpha=0.05,  theta1=0.8, theta2=1.25, theta0=0.95, cv=0.15, method=:owenq)
besamplen(cv=0.20, method=:nct)
besamplen(cv=0.347, design=:parallel,  out=:print)
besamplen(cv=0.40)
n, p, s = besamplen(cv=0.347, design=:d2x2x4, method=:nct, out=:vstr)

#Bioequivalence power for 2x2 design, default method - OwensQ
bePower(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.2, n=20, design=:d2x2, method=:owenq)
#Same
bePower(alpha=0.05, cv=0.2, n=20, design=:d2x2)
#Bioequivalence power for cv 14%, 21 subjects, default OwensQ method, logscale false
bePower(alpha=0.1, logscale=false, theta1=-0.1, theta2=0.1, theta0=0, cv=0.14, n=21)
#Bioequivalence power for cv 14%, 21 subjects, shifted method, logscale false
bePower(alpha=0.1, logscale=false, theta1=-0.1, theta2=0.1, theta0=0, cv=0.14, n=21, method=:shifted)
#Simple notations
bePower(cv=0.4, n=35, design=:d2x4x4)
bePower(cv=0.14, n=21)

#CV from CI
ci2cv(;alpha = 0.05, theta1 = 0.9, theta2 = 1.25, n=30, design=:d2x2x4)

#Polled CV
data = DataFrame(cv = Float64[], df = Int[])
push!(data, (0.12, 12))
push!(data, (0.2, 20))
push!(data, (0.25, 30))
pooledCV(data; cv=:cv, df=:df, alpha=0.05, returncv=true)

```

### <a name="Other">Other functions</a>

- Owen's T function:

```
owensT(h::Float64,a::Float64)::Float64
```

- Owen's Q function (a,b always should be >= 0):

```
owensQ(nu, t::Float64, delta::Float64, a::Float64, b::Float64)::Float64
```

### <a name="Dependencies">Dependencies</a>

 - Distributions
 - StatsBase
 - Statistics
 - QuadGK
 - SpecialFunctions
 - Random
 - Roots
 - DataFrames

 ### <a name="Copyrights&References">Copyrights&References</a>

 > Clinical Trial Utilities

 > Author: Vladimir Arnautov aka PharmCat

 > Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)

 > OwensQ/PowerTOST functions rewritten from https://github.com/Detlew/PowerTOST by Detlew Labes, Helmut Schuetz, Benjamin Lang

 > Calculation based on Chow S, Shao J, Wang H. 2008. Sample Size Calculations in Clinical Research. 2nd Ed. Chapman & Hall/CRC Biostatistics Series.

 > Connor R. J. 1987. Sample size for testing differences in proportions for the paired-sample design. Biometrics 43(1):207-211. page 209.

 >Phillips KF.Power of the Two One-Sided Tests Procedure in BioequivalenceJ Pharmacokin Biopharm. 1990;18(2):137–44. doi: 10.1007/BF01063556

 >Diletti D, Hauschke D, Steinijans VW.Sample Size Determination for Bioequivalence Assessment by Means of Confidence IntervalsInt J Clin Pharmacol Ther Toxicol. 1991;29(1):1–8.

 > Owen, D B (1965) "A Special Case of a Bivariate Non-central t-Distribution" Biometrika Vol. 52, pp.437-446.

 > FORTRAN code by J. Burkhardt, license GNU LGPL

 > D.B. Owen "Tables for computing bivariate normal Probabilities" The Annals of Mathematical Statistics, Vol. 27 (4) Dec. 1956, pp. 1075-1090

 > matlab code  by J. Burkhardt license GNU LGPL

 > If you want to check and get R code - you can find some here: http://powerandsamplesize.com/Calculators/

 > Some ideas was taken from R project packages:

 > PropCIs by Ralph Scherer https://cran.r-project.org/web/packages/PropCIs/index.html

 > pairwiseCI by Frank Schaarschmidt, Daniel Gerhard  https://CRAN.R-project.org/package=pairwiseCI

 > binGroup by Boan Zhang, Christopher Bilder, Brad Biggerstaff, Frank Schaarschmidt Brianna Hitt https://CRAN.R-project.org/package=binGroup

 > proportion by M.Subbiah, V.Rajeswaran https://CRAN.R-project.org/package=proportion

 > binom by Sundar Dorai-Raj https://CRAN.R-project.org/package=binom

 > DescTools https://CRAN.R-project.org/package=DescTools

 > ORCI by Libo Sun https://CRAN.R-project.org/package=ORCI

 > metafor by Wolfgang Viechtbauer https://cran.r-project.org/package=metafor
