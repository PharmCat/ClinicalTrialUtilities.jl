# Confidence Intervals

## Proportions

### One proportion

P = n / x

where n - number of outcomes; x - number of observations.

### Absolute risk difference

Diff(𝛿) = P₁ - P₂ = n₁ / x₁ - n₂ / x₂

### Risk Ratio

RR = P₁ / P₂ = (n₁ / x₁) - (n₂ / x₂)

### Odd Ratio

OR = (n₁ / (x₁ - n₁)) - (n₂ / (x₂ - n₂))

### propci
```@docs
ClinicalTrialUtilities.propci
```

### diffpropci
```@docs
ClinicalTrialUtilities.diffpropci
```

### rrpropci
```@docs
ClinicalTrialUtilities.rrpropci
```

### orpropci
```@docs
ClinicalTrialUtilities.orpropci
```

## Means

### meanci
```@docs
ClinicalTrialUtilities.meanci
```

### diffmeanci
```@docs
ClinicalTrialUtilities.diffmeanci
```

## Cochran–Mantel–Haenszel confidence intervals

Table cell map:


| group   | outcome 1 | outcome 2 |
|---------|-----------|-----------|
| group 1 |     a     |     b     |
| group 2 |     c     |     d     |


### diffcmhci
```@docs
ClinicalTrialUtilities.diffcmhci
```

### orcmhci
```@docs
ClinicalTrialUtilities.orcmhci
```

### rrcmhci
```@docs
ClinicalTrialUtilities.rrcmhci
```
