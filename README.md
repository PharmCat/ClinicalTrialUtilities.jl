# ClinicalTrialUtilities

 Clinical trial related calculation:  power and sample size calculation, randomization. This program comes with absolutely no warranty. No liability is accepted for any loss and risk to public health resulting from use of this software.

![Tier 1](https://github.com/PharmCat/ClinicalTrialUtilities.jl/workflows/Tier%201/badge.svg)

[![codecov](https://codecov.io/gh/PharmCat/ClinicalTrialUtilities.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/PharmCat/ClinicalTrialUtilities.jl)

[![Latest docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://pharmcat.github.io/ClinicalTrialUtilities.jl/dev/)


### <a name="Installation">Installation</a>
```
using Pkg; Pkg.add("ClinicalTrialUtilities");
```

### <a name="Features">Main features</a>

- Clinical trial sample size calculation
- Power calculation
- Randomization


### <a name="Examples">Examples</a>

#### SampleSize

**NB! Hypothesis types:**

- :ea - Equality: two-sided;
- :ei - Equivalencens: two one-sided hypothesis (TOST);
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

#### Randomization
```
using DataFrames, ClinicalTrialUtilities
rt = ClinicalTrialUtilities.randomtable(;blocksize = 4, subject = 32, group = 2, ratio = [1,1], grseq = ["TR", "RT"], seed = 36434654652452)
```

#### Confidence Intervals

Proportion CI moved to [MetidaFreq.jl](https://github.com/PharmCat/MetidaFreq.jl)


#### NCA

NCA moved to [MetidaNCA.jl](https://github.com/PharmCat/MetidaNCA.jl)

### <a name="Copyrights">Copyrights</a>

Clinical Trial Utilities

Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)

If you want to check and get R code for power/sample size estimation, you can find examples here: http://powerandsamplesize.com/Calculators/
