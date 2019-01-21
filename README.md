# CTPSS
 Clinical Trial Power and Sample Size calculation

Alpha version

Author: Vladimir Arnautov aka PharmCat

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

### Functions:

sampleSize (;param, type, group, alpha, beta, diff, sd, a, b, k)

### Usage:
```
using CTPSS
sampleSize(param="mean|prop|or", type="ea|ei|ns", group="one|two", alpha=0.05, beta=0.2, diff=1, sd=2, a=1, b=2, k=1)

```

**param (Parameter type):**
- mean - Means;
- prop - Proportions;
- or - Odd Ratio;

**type (Hypothesis type):**
- ea - Equality;
- ei - Equivalencens;
- ns - Non-Inferiority / Superiority;

**group (Group num):**
- one - one sample;
- two - Two sample;

**alpha** - Alpha (o < alpha < 1);
**beta** - Beta (o < beta < 1); power = 1 - beta;

**diff** - difference /equivalencens margin/non-inferiority/superiority margin;
**sd** - Standart deviation (σ, for Means only);
**a** - Null Hypothesis mean (μ0), Group A;
**b** - True mean (μ) for one sample / Group B for two sample design;
**k** - Na/Nb (after sample size estimation second grop size: Na = k * Nb, only for two sample design);

### ToDo:

 - powerEst ()
 - owensQ(nu, t, delta, a, b)
 - knownDesign () from PowerTost
 -  sampleSizeBE ()
 - powerBE ()
 - cvfromci ()
 - Simulations.
