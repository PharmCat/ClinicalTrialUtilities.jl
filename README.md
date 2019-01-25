# CTPSS
 Clinical Trial Power and Sample Size calculation

Version:

Author: Vladimir Arnautov

[![Build Status](https://travis-ci.com/PharmCat/CTPSS.svg?branch=master)](https://travis-ci.com/PharmCat/CTPSS)
[![Build status](https://ci.appveyor.com/api/projects/status/c7x06t501eqjvd6s?svg=true)](https://ci.appveyor.com/project/PharmCat/ctpss)
[![codecov](https://codecov.io/gh/PharmCat/CTPSS/branch/master/graph/badge.svg)](https://codecov.io/gh/PharmCat/CTPSS)
[![Coverage Status](https://coveralls.io/repos/github/PharmCat/CTPSS/badge.svg?branch=master)](https://coveralls.io/github/PharmCat/CTPSS?branch=master)
### Dependencies:

 - Distributions
 - QuadGK
 - SpecialFunctions
 - Rmath

### Install:
```
using Pkg
Pkg.add("git://github.com/PharmCat/CTPSS.git")
```

### Function

sampleSize (;param, type, group, alpha, beta, diff, sd, a, b, k)

Sample size calculation.

### Usage:
```
using CTPSS
sampleSize(param="mean|prop|or", type="ea|ei|ns", group="one|two", alpha=0.05, beta=0.2, diff=1, sd=2, a=1, b=2, k=1)

```

**param (Parameter type):**
- mean - Means (default);
- prop - Proportions;
- or - Odd Ratio;

**type (Hypothesis type):**
- ea - Equality  (default);
- ei - Equivalencens;
- ns - Non-Inferiority / Superiority;

**group (Group num):**
- one - one sample  (default);
- two - Two sample;

**alpha** - Alpha (o < alpha < 1)  (default=0.05);
**beta** - Beta (o < beta < 1) (default=0.2); power = 1 - beta;

**diff** - difference /equivalencens margin/non-inferiority/superiority margin;
**sd** - Standart deviation (σ, for Means only);
**a** - Null Hypothesis mean (μ0), Group A;
**b** - True mean (μ) for one sample / Group B for two sample design;
**k** - Na/Nb (after sample size estimation second grop size: Na = k * Nb, only for two sample design) (default=1);

### Function

powerTOST(;alpha, logscale, theta1, theta2, theta0, cv, n, method)

Power calculation for TOST (for Bioequivalence trials).

### Usage

```
using CTPSS
powerTOST(alpha=0.05, logscale=true|false, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.15, n=36, method="owenq|nct|shifted")
```
**logscale** - theta1, theta2, theta0: if true - make log transformation;

**alpha** - Alpha (o < alpha < 1)  (default=0.05);

**beta** - Beta (o < beta < 1) (default=0.2); power = 1 - beta;

**theta1** - Lower Level;

**theta2** - Upper level;

**theta0** - T/R Ratio;

**cv** - coefficient of variation;

**n** - subject number;

**method** - calculating method: Oqen'sQ Function | NonCentral T, Shifted;

###Example:

```
sampleSize(param="prop", type="ea", group="one", a=0.3, b=0.5)
sampleSize(param="mean", type="ei", group="two", diff=0.3, sd=1, a=0.3, b=0.5)
sampleSize(param="or", type="ns", diff=-0.1, a=0.3, b=0.5, k=2)
sampleSize(param="or", type="ea", a=0.3, b=0.5, k=2)

powerTOST(alpha=0.1, logscale=false, theta1=-0.1, theta2=0.1, theta0=0, cv=0.14, n=21, method="shifted")
powerTOST(alpha=0.1, logscale=false, theta1=-0.1, theta2=0.1, theta0=0, cv=0.14, n=21)
powerTOST(cv=0.14, n=21)
```

### ToDo:

 - powerEst ()
 - knownDesign () from PowerTost
 - sampleSizeBE ()
 - cvfromci ()
 - Simulations.
