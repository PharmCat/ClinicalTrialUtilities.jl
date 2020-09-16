# SampleSize

Sample size estimation.

## Introduction

### Hypothesis

#### Equality

Two-sided hypothesis:

H₀: A - B = 0
Hₐ: A - B ≠ 0

100\*(1 - α)% two-sided confidence interval used to demonstrate it.

#### Equivalence | Bioequivalence

Two-sided hypothesis or two one-sided hypothesis (TOST):

H₀: |A − B| ≥ δ
Hₐ: |A − B| < δ

or

H₀: A - B ≤ δˡ || A - B ≥ δᵘ
Hₐ: δˡ < A - B < δᵘ

100\*(1 - 2\*α)% two-sided confidence interval used to demonstrate it.

See Walker&Nowacki, 2011, Understanding Equivalence and Noninferiority Testing.

*Chow at al. Sample Size Calculations in Clinical Research 3-rd ed. for sample size calculation used alpha level for testing corresponding 100\*(1 - α)% two-sided confidence interval*

For bioequivalence clinical trial 90% *(100\*(1 - 2\*α)%)* two-sided confidence interval is used (α = 0.05), but for therapeutic equivalence 95% two-sided confidence interval should be used (in this case use α = 0.025).

#### Non-inferiority | Superiority

One-sided hypothesis:

H₀: A − B ≤ δ
Hₐ: A − B > δ

100\*(1 - α)% one-sided confidence interval used to demonstrate it.

In clinical research α = 0.025 or less should be used, or 100*(1 - α)% two-sided confidence interval should be used.

### ctsamplen
```@docs
ClinicalTrialUtilities.ctsamplen
```

### besamplen
```@docs
ClinicalTrialUtilities.besamplen
```
