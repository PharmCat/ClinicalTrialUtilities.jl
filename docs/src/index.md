# ClinicalTrialUtilities

 Clinical trial related calculation: power and sample size calculation, randomization.

[![Build Status](https://travis-ci.com/PharmCat/ClinicalTrialUtilities.jl.svg?branch=master)](https://travis-ci.com/PharmCat/ClinicalTrialUtilities.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/35f8b5vq259sbssg?svg=true)](https://ci.appveyor.com/project/PharmCat/clinicaltrialutilities-jl)
[![codecov](https://codecov.io/gh/PharmCat/ClinicalTrialUtilities.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/PharmCat/ClinicalTrialUtilities.jl)
[![Coverage Status](https://coveralls.io/repos/github/PharmCat/ClinicalTrialUtilities.jl/badge.svg?branch=master)](https://coveralls.io/github/PharmCat/ClinicalTrialUtilities.jl?branch=master)
[![Latest docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://pharmcat.github.io/ClinicalTrialUtilities.jl/dev/)

## Description

The package is designed to perform calculations related to the planning of clinical trials.
## Installation
```
using Pkg; Pkg.add("ClinicalTrialUtilities");
```

### Pharmacodynamics

Further development of NCA PK/PD will be based on [MetidaNCA.jl](https://github.com/PharmCat/MetidaNCA.jl) package.

All NCA PK/PD functions moved to [MetidaNCA](https://github.com/PharmCat/MetidaNCA.jl).

Descriptive statistics moved to  [MetidaStats](https://github.com/PharmCat/MetidaStats.jl).

Confidence intervals moved to  [MetidaFreq](https://github.com/PharmCat/MetidaFreq.jl).

Simulations removed.


## Note

**NB! Hypothesis types:**

When ctpower/ctsamplen used:

- :ea - Equality: two-sided hypothesis;
- :ei - Equivalencens: two one-sided hypothesis;
- :ns - Non-Inferiority / Superiority: one-sided hypothesis, for some cases you should use two-sided hypothesis for  Non-Inferiority/Superiority, you can use alpha/2 for this;

## Invalidated version v0.3.0 (removed from release page)

  Wrong Equivalence Hypothesis alpha level! Use 0.2.7, 0.3.1 version or higher.

## Contents

```@contents
Pages = [

        "samplesize.md",
        "power.md",
        "random.md",
        "utils.md",
        "examples.md",
        "export.md",
        "ref.md"]
Depth = 4
```

## Copyrights


Clinical Trial Utilities

Copyright © 2019 Vladimir Arnautov aka PharmCat <mail@pharmcat.net>

If you want to check and get R code for power/sample size estimation, you can find examples here: http://powerandsamplesize.com/Calculators/
