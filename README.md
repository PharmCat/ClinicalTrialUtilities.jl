# ClinicalTrialUtilities

 Clinical trial related calculation: descriptive statistics, power and sample size calculation, power simulations, confidence interval, pharmacokinetics/pharmacodynamics parameters calculation. This program comes with absolutely no warranty. No liability is accepted for any loss and risk to public health resulting from use of this software.

[![Build Status](https://travis-ci.com/PharmCat/ClinicalTrialUtilities.jl.svg?branch=master)](https://travis-ci.com/PharmCat/ClinicalTrialUtilities.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/35f8b5vq259sbssg?svg=true)](https://ci.appveyor.com/project/PharmCat/clinicaltrialutilities-jl)
[![codecov](https://codecov.io/gh/PharmCat/ClinicalTrialUtilities.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/PharmCat/ClinicalTrialUtilities.jl)
[![Coverage Status](https://coveralls.io/repos/github/PharmCat/ClinicalTrialUtilities.jl/badge.svg?branch=master)](https://coveralls.io/github/PharmCat/ClinicalTrialUtilities.jl?branch=master)
[![Latest docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://pharmcat.github.io/ClinicalTrialUtilities.jl/dev/)

## Description

The package is designed to perform calculations related to the planning and analysis of the results of clinical trials. The package includes the basic functions described below, as well as a few modules to perform specific calculations.

### <a name="Installation">Installation</a>
```
using Pkg; Pkg.add("ClinicalTrialUtilities");
```

### <a name="Features">Main features</a>

- SampleSize calculation
- Power calculation
- Confidence Interval calculation
- NCA Pharmacokinetics parameters calculation
- Randomization


### <a name="Examples">Examples</a>

#### SampleSize

**NB! Hypothesis types:**

- :ea - Equality: two-sided;
- :ei - Equivalencens: two one-sided hypothesis;
- :ns - Non-Inferiority / Superiority: one-sided hypothesis, for some cases you should use two-sided hypothesis for  Non-Inferiority/Superiority, you can use alpha/2 for this;

```
#Sample size for one proportion equality
ctsamplen(param=:prop, type=:ea, group=:one, a=0.3, b=0.5)

#Equivalence for two means
ctsamplen(param=:mean, type=:ei, group=:two, diff=0.3, sd=1, a=0.3, b=0.5)

#Odd ratio non-inferiority
ctsamplen(param=:or, type=:ns, diff=-0.1, a=0.3, b=0.5, k=2)

#Odd ratio equality
ctsamplen(param=:or, type=:ea, a=0.3, b=0.5, k=2)
```

#### Bioequivalence sample size
```
besamplen(alpha=0.05,  theta1=0.8, theta2=1.25, theta0=0.95, cv=0.15, method=:owenq)
besamplen(cv=0.20, method=:nct)
besamplen(cv=0.347, design=:parallel)
besamplen(cv=0.40)
besamplen(cv=0.347, design=:d2x2x4, method=:nct)
```

#### Power
```
ctpower(param=:mean, type=:ea, group=:one, a=1.5, b=2, sd=1,n=32, alpha=0.05)
```

#### Bioequivalence power
```
#2x2 design, default method - OwensQ
bepower(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.2, n=20, design=:d2x2, method=:owenq)

#Same
bepower(alpha=0.05, cv=0.2, n=20, design=:d2x2)

#Bioequivalence power for cv 14%, 21 subjects, default OwensQ method, logscale false
bepower(alpha=0.1, logscale=false, theta1=-0.1, theta2=0.1, theta0=0, cv=0.14, n=21)

#Bioequivalence power for cv 14%, 21 subjects, shifted method, logscale false
bepower(alpha=0.1, logscale=false, theta1=-0.1, theta2=0.1, theta0=0, cv=0.14, n=21, method=:shifted)

#Simple notations
bepower(cv=0.4, n=35, design=:d2x4x4)
bepower(cv=0.14, n=21)
```

#### Bioequivalence CV from CI
```
cvfromci(;alpha = 0.05, theta1 = 0.9, theta2 = 1.25, n=30, design=:d2x2x4)
```

#### Polled CV
```
data = DataFrame(cv = Float64[], df = Int[])
push!(data, (0.12, 12))
push!(data, (0.2, 20))
push!(data, (0.25, 30))
pooledcv(data; cv=:cv, df=:df, alpha=0.05, returncv=true)


pooledcv([0.12, 0.2, 0.25], [14, 22, 32], [:d2x2, :d2x2, :d2x2])

```

#### NCA
```
using CSV, DataFrames, ClinicalTrialUtilities
pkdata2 = CSV.File("pkdata.csv") |> DataFrame
pkds    = pkimport(pkdata2, [:Subject, :Formulation]; time = :Time, conc = :Concentration)
pk      = nca!(pkds)
ncadf   = DataFrame(pk; unst = true)
ds      = ClinicalTrialUtilities.descriptive(ncadf, stats = [:n, :mean, :sd], sort = [:Formulation])
dsdf    = ClinicalTrialUtilities.DataFrame(ds; unst = true)

```


### <a name="Copyrights">Copyrights</a>

Clinical Trial Utilities

Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)

If you want to check and get R code for power/sample size estimation, you can find examples here: http://powerandsamplesize.com/Calculators/
