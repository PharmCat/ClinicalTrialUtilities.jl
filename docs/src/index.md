# ClinicalTrialUtilities

**VERSION 0.2.0 INCOMATIBLE WITH 0.1.x**

 Clinical trial related calculation: descriptive statistics, power and sample size calculation, power simulations, confidence interval, pharmacokinetics/pharmacodynamics parameters calculation.

[![Build Status](https://travis-ci.com/PharmCat/ClinicalTrialUtilities.jl.svg?branch=master)](https://travis-ci.com/PharmCat/ClinicalTrialUtilities.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/35f8b5vq259sbssg?svg=true)](https://ci.appveyor.com/project/PharmCat/clinicaltrialutilities-jl)
[![codecov](https://codecov.io/gh/PharmCat/ClinicalTrialUtilities.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/PharmCat/ClinicalTrialUtilities.jl)
[![Coverage Status](https://coveralls.io/repos/github/PharmCat/ClinicalTrialUtilities.jl/badge.svg?branch=master)](https://coveralls.io/github/PharmCat/ClinicalTrialUtilities.jl?branch=master)


```@meta
CurrentModule = ClinicalTrialUtilities
DocTestSetup = quote
    using StatsBase
end
```

## Description

The package is designed to perform calculations related to the planning and analysis of the results of clinical trials. The package includes the basic functions described below, as well as a few modules to perform specific calculations.

```@contents
Pages = ["samplesize.md", "power.md", "ci.md", "ds.md", "nca.md", "pk.md", "pd.md", "random.md", "sim.md", "utils.md", "examples.md", "export.md", "validation.md", "ref.jl"]
Depth = 2
```
