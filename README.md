# ClinicalTrialUtilities
 Clinical Trial Power and Sample Size calculation

Version:0.1.5

Author: Vladimir Arnautov

2019 &copy; Vladimir Arnautov

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

sampleSize (;param, type, group, alpha, beta, diff, sd, a, b, k, out="num|str|vstr|print") - Sample size calculation.

ctPower(;param, type, group, alpha, n, diff, sd, a, b, k, out="num|str|vstr|print") - Clinical trial power estimation.

powerTOST(;theta0=0.95, theta1=0.8, theta2=1.25, cv, n, alpha=0.5, logscale=true, method="owenq",  design="parallel|2x2|2x2x3|2x2x4|2x4x4|2x3x3") - Power calculation for TOST (for Bioequivalence trials).

beSampleN(;theta0=0.95, theta1=0.8, theta2=1.25, cv, alpha=0.05, beta=0.2, logscale=true, method="owenq", design="parallel|2x2|2x2x3|2x2x4|2x4x4|2x3x3", out="num|str|vstr|print") - iterative sample size calculation for Bioequivalence trials.

owensT(h::Float64,a::Float64)::Float64 - Owen's T function.

owensQ(nu, t, delta, a, b) - Owen's Q function; a,b always should be >= 0.

### Usage

#### sampleSize
```
using ClinicalTrialUtilities
sampleSize(param="mean|prop|or", type="ea|ei|ns|mcnm", group="one|two", alpha=0.05, beta=0.2, diff=1, sd=2, a=1, b=2, k=1)

```
**param (Parameter type):**
- mean - Means (default);
- prop - Proportions;
- or - Odd Ratio;

**type (Hypothesis type):**
- ea - Equality  (default);
- ei - Equivalencens;
- ns - Non-Inferiority / Superiority;
- mcnm - McNemar's Equality test;

**group (Group num):**
- one - One sample  (default);
- two - Two sample, result is for one group, second group size = n * k;

**alpha** - Alpha (o < alpha < 1)  (default=0.05);

**beta** - Beta (o < beta < 1) (default=0.2); power = 1 - beta;

**diff** - difference / equivalence margin/non-inferiority/superiority margin;

**sd** - Standard deviation (σ, for Means only);

**a** - Null Hypothesis mean (μ0), Group A;

**b** - True mean (μ) for one sample / Group B for two sample design;

**k** - Na/Nb (after sample size estimation second group size: Na = k * Nb, only for two sample design) (default=1);

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
- mean - Means (default);
- prop - Proportions;
- or - Odd Ratio;

**type (Hypothesis type):**
- ea - Equality  (default);
- ei - Equivalence;
- ns - Non-Inferiority / Superiority;
- mcnm - McNemar's Equality test;

**group (Group num):**
- one - one sample  (default);
- two - Two sample;

**alpha** - Alpha (o < alpha < 1)  (default=0.05);

**n** - Subjects number;

**diff** - difference /equivalence margin/non-inferiority/superiority margin;

**sd** - Standard deviation (σ, for Means only);

**a** - Null Hypothesis mean (μ0), Group A;

**b** - True mean (μ) for one sample / Group B for two sample design;

**k** - Na/Nb (after sample size estimation second group size: Na = k * Nb, only for two sample design) (default=1);

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
**logscale** - theta1, theta2, theta0: if true - make log transformation;

**alpha** - Alpha (o < alpha < 1)  (default=0.05);

**theta1** - Lower Level (default=0.8);

**theta2** - Upper level (default=1.25);

**theta0** - T/R Ratio (default=0.95);

**cv** - coefficient of variation;

**n** - subject number;

**method** - calculating method: Owen's Q Function | NonCentral T, Shifted;
- owenq
- nct
- shifted

#### beSampleN

Using for bioequivalence study.

```
using ClinicalTrialUtilities
beSampleN(alpha=0.05, logscale=true|false, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.15, method="owenq|nct|shifted", design="parallel|2x2|2x2x3|2x2x4|2x4x4|2x3x3")
```
**logscale** - theta1, theta2, theta0: if true - make log transformation;

**alpha** - Alpha (o < alpha < 1)  (default=0.05);

**beta** - Beta (o < beta < 1) (default=0.2); power = 1 - beta;

**theta1** - Lower Level (default=0.8);

**theta2** - Upper level (default=1.25);

**theta0** - T/R Ratio (default=0.95);

**cv** - coefficient of variation;

**method** - calculating method: Owen's Q Function | NonCentral T | Shifted;
- owenq
- nct
- shifted

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

### ToDo:

 - atomic function interface with struct
 - cvfromci ()
 - Simulations
