# ClinicalTrialUtilities

**VERSION 0.2.0 INCOMATIBLE WITH 0.1.x**

 Clinical trial related calculation: descriptive statistics, power and sample size calculation, power simulations, confidence interval, pharmacokinetics/pharmacodynamics parameters calculation.

[![Build Status](https://travis-ci.com/PharmCat/ClinicalTrialUtilities.jl.svg?branch=master)](https://travis-ci.com/PharmCat/ClinicalTrialUtilities.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/35f8b5vq259sbssg?svg=true)](https://ci.appveyor.com/project/PharmCat/clinicaltrialutilities-jl)
[![codecov](https://codecov.io/gh/PharmCat/ClinicalTrialUtilities.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/PharmCat/ClinicalTrialUtilities.jl)
[![Coverage Status](https://coveralls.io/repos/github/PharmCat/ClinicalTrialUtilities.jl/badge.svg?branch=master)](https://coveralls.io/github/PharmCat/ClinicalTrialUtilities.jl?branch=master)
[![Latest docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://pharmcat.github.io/ClinicalTrialUtilities.jl/dev/)

## Description

The package is designed to perform calculations related to the planning and analysis of the results of clinical trials. The package includes the basic functions described below, as well as a few modules to perform specific calculations.

## Installation
```
using Pkg; Pkg.add("ClinicalTrialUtilities");
```

## Note

**NB! Hypothesis types:**

When ctpower/ctsamplen used:

- :ea - Equality: two-sided hypothesis;
- :ei - Equivalencens: two one-sided hypothesis;
- :ns - Non-Inferiority / Superiority: one-sided hypothesis, for some cases you should use two-sided hypothesis for  Non-Inferiority/Superiority, you can use alpha/2 for this;

## Contents

```@contents
Pages = [
        "index.md",
        "samplesize.md",
        "power.md",
        "ci.md",
        "ds.md",
        "nca.md",
        "pk.md",
        "pd.md",
        "random.md",
        "sim.md",
        "utils.md",
        "examples.md",
        "export.md",
        "validation.md",
        "ref.md"]
Depth = 3
```

## Copyrights

Clinical Trial Utilities

Copyright © 2019 Vladimir Arnautov aka PharmCat <mail@pharmcat.net>

If you want to check and get R code for power/sample size estimation, you can find examples here: http://powerandsamplesize.com/Calculators/
