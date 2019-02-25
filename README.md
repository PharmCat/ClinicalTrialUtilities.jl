# ClinicalTrialUtilities
 Clinical Trial Power and Sample Size calculation

Version:0.1.7

Author: Vladimir Arnautov

2019 &copy; Vladimir Arnautov

[![Build Status](https://travis-ci.com/PharmCat/ClinicalTrialUtilities.jl.svg?branch=master)](https://travis-ci.com/PharmCat/ClinicalTrialUtilities.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/35f8b5vq259sbssg?svg=true)](https://ci.appveyor.com/project/PharmCat/clinicaltrialutilities-jl)
[![codecov](https://codecov.io/gh/PharmCat/ClinicalTrialUtilities.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/PharmCat/ClinicalTrialUtilities.jl)
[![Coverage Status](https://coveralls.io/repos/github/PharmCat/ClinicalTrialUtilities.jl/badge.svg?branch=master)](https://coveralls.io/github/PharmCat/ClinicalTrialUtilities.jl?branch=master)
### Dependencies:

 - Distributions
 - QuadGK
 - SpecialFunctions

### Install:
```
using Pkg
Pkg.add("ClinicalTrialUtilities");
```
or
```
using Pkg
Pkg.clone("https://github.com/PharmCat/ClinicalTrialUtilities.jl.git");
```
### Functions:

Sample size calculation:
sampleSize (;param, type, group, alpha, beta, diff, sd, a, b, k, out="num|str|vstr|print")

Clinical trial power estimation:
ctPower(;param, type, group, alpha, n, diff, sd, a, b, k, out)

Power calculation for TOST (for Bioequivalence trials):
powerTOST(;theta0, theta1, theta2, cv, n, alpha, logscale, method,  design)

Iterative sample size calculation for Bioequivalence trials:
beSampleN(;theta0, theta1, theta2, cv, alpha, beta, logscale, method, design, out)

Owen's T function:
owensT(h::Float64,a::Float64)::Float64

Owen's Q function (a,b always should be >= 0):
owensQ(nu, t::Float64, delta::Float64, a::Float64, b::Float64)::Float64

### Usage

#### sampleSize
```
using ClinicalTrialUtilities
sampleSize(param="mean|prop|or", type="ea|ei|ns|mcnm", group="one|two", alpha=0.05, beta=0.2, diff=1, sd=2, a=1, b=2, k=1)

```
**param (Parameter type):**
- mean - Means;
- prop - Proportions;
- or - Odd Ratio;

**type (Hypothesis type):**
- ea - Equality;
- ei - Equivalencens;
- ns - Non-Inferiority / Superiority;
- mcnm - McNemar's Equality test;

**group (group num):**
- one - One sample;
- two - Two sample, result is for one group, second group size = n * k;

**alpha** - Alpha (o < alpha < 1)  (default=0.05);

**beta** - Beta (o < beta < 1) (default=0.2); power = 1 - beta;

**diff** - difference/equivalence margin/non-inferiority/superiority margin;

**sd** - Standard deviation (σ, for Means only);

**a** - Null Hypothesis mean (μ0), Group A;

**b** - True mean (μ) for one sample / Group B for two sample design;

**k** - Na/Nb (after sample size estimation second group size: Na=k*Nb, only for two sample design) (default=1);

**out** - output type:
- num   - numeric (default);
- str   - String variable with text output;
- vstr  - numeric and String variable;
- print - print to console;

#### ctPower
```
using ClinicalTrialUtilities
ctPower(param="mean|prop|or", type="ea|ei|ns|mcnm", group="one|two", alpha=0.05, n=30, diff=1, sd=2, a=1, b=2, k=1)
```

**param (Parameter type):**
- mean - Means;
- prop - Proportions;
- or - Odd Ratio;

**type (Hypothesis type):**
- ea - Equality;
- ei - Equivalence;
- ns - Non-Inferiority / Superiority;
- mcnm - McNemar's Equality test;

**group (group num):**
- one - one sample;
- two - Two sample;

**alpha** - Alpha (0<alpha<1)  (default=0.05);

**n** - Subjects number;

**diff** - difference/equivalence margin/non-inferiority/superiority margin;

**sd** - Standard deviation (σ, for Means only);

**a** - Null Hypothesis mean (μ0), Group A;

**b** - True mean (μ) for one sample / Group B for two sample design;

**k** - Na/Nb (after sample size estimation second group size: Na=k*Nb, only for two sample design) (default=1);

**out** - output type:
- num   - numeric (default);
- str   - String variable with text output;
- vstr  - numeric and String variable;
- print - print to console;

#### powerTOST

```
using ClinicalTrialUtilities
powerTOST(alpha=0.05, logscale=true|false, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.15, n=36, method="owenq|nct|shifted", design="parallel|2x2|2x2x3|2x2x4|2x4x4|2x3x3")
```
**logscale** - theta1, theta2, theta0: if true - make log transformation (default true);

**alpha** - Alpha (0 < alpha < 1)  (default=0.05);

**theta1** - Lower Level (default=0.8);

**theta2** - Upper level (default=1.25);

**theta0** - T/R Ratio (default=0.95);

**cv** - coefficient of variation;

**n** - subject number;

**method** - calculating method: Owen's Q Function | NonCentral T, Shifted;
- owenq (default)
- nct
- shifted

**design** - trial design;
- parralel
- 2x2 (default)
- 2x2x4
- 2x4x4
- 2x3x3

#### beSampleN

Using for bioequivalence study.

```
using ClinicalTrialUtilities
beSampleN(alpha=0.05, logscale=true|false, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.15, method="owenq|nct|shifted", design="parallel|2x2|2x2x3|2x2x4|2x4x4|2x3x3")
```
**logscale** - theta1, theta2, theta0: if true - make log transformation (default true);

**alpha** - Alpha (o < alpha < 1)  (default=0.05);

**beta** - Beta (o < beta < 1) (default=0.2); power = 1 - beta;

**theta1** - Lower Level (default=0.8);

**theta2** - Upper level (default=1.25);

**theta0** - T/R Ratio (default=0.95);

**cv** - coefficient of variation;

**method** - calculating method: Owen's Q Function | NonCentral T | Shifted;
- owenq (default)
- nct
- shifted

**design** - trial design;
- parralel
- 2x2 (default)
- 2x2x4
- 2x4x4
- 2x3x3

**out** - output type:
- num   - numeric (default);
- str   - String variable with text output;
- vstr  - numeric and String variable;
- print - print to console;

### Examples:

```
sampleSize(param="prop", type="ea", group="one", a=0.3, b=0.5)
sampleSize(param="mean", type="ei", group="two", diff=0.3, sd=1, a=0.3, b=0.5)
sampleSize(param="or", type="ns", diff=-0.1, a=0.3, b=0.5, k=2)
sampleSize(param="or", type="ea", a=0.3, b=0.5, k=2)

powerTOST(alpha=0.1, logscale=false, theta1=-0.1, theta2=0.1, theta0=0, cv=0.14, n=21, method="shifted")
powerTOST(alpha=0.1, logscale=false, theta1=-0.1, theta2=0.1, theta0=0, cv=0.14, n=21)
powerTOST(cv=0.4, n=35, design="2x4x4")
powerTOST(cv=0.14, n=21)

beSampleN(alpha=0.05,  theta1=0.8, theta2=1.25, theta0=0.95, cv=0.15, method="owenq")
beSampleN(cv=0.20, method="nct")
beSampleN(cv=0.347, design="parallel",  out="print")
beSampleN(cv=0.40)

n, p, s = beSampleN(cv=0.347, design="2x2x4", method="nct", out="vstr")
```

### Confidence Interval Submodule



### ToDo:

 - atomic function interface with struct
 - cvfromci ()
 - Simulations
