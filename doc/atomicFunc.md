

μ₀ - Test Mean
μ₁ - Reference Mean
p₀ - Test Proportion
p₁ - Reference Proportion
σ  - SD
δ  - Margin Difference
α  - Alpha (Type I Error)
β  - Beta  (Type II Error)

# Sample Size

## Means

### One Sample Means

```
one_mean_equality(μ₀::Real, μ₁::Real, σ::Real, α::Float64, β::Float64)::Float64

one_mean_equivalence(μ₀::Real, μ₁::Real, σ::Real, δ::Real, α::Float64, β::Float64)::Float64

one_mean_superiority(μ₀::Real, μ₁::Real, σ::Real, δ::Real, α::Float64, β::Float64)::Float64
```

### Two Sample Means

```
two_mean_equality(μ₀::Real, μ₁::Real, σ::Real, α::Float64, β::Float64, k::Real)::Float64

two_mean_equivalence(μ₀::Real, μ₁::Real, σ::Real, δ::Real, α::Float64, β::Float64, k::Real)::Float64

two_mean_superiority(μ₀::Real, μ₁::Real, σ::Real, δ::Real, α::Float64, β::Float64, k::Real)::Float64
```

## Compare Proportion

### One Sample Proportion

```
one_proportion_equality(p₀::Float64, p₁::Float64, α::Float64, β::Float64)::Float64

one_proportion_equivalence(p₀::Float64, p₁::Float64, δ::Real, α::Float64, β::Float64)::Float64

one_proportion_superiority(p₀::Float64, p₁::Float64, δ::Real, α::Float64, β::Float64)::Float64
```

### Two Sample Proportion

```
two_proportion_equality(p₀::Float64, p₁::Float64, α::Float64, β::Float64, k::Real)::Float64

two_proportion_equivalence(p₀::Float64, p₁::Float64, δ::Real, α::Float64, β::Float64, k::Real)::Float64

two_proportion_superiority(p₀::Float64, p₁::Float64, δ::Real, α::Float64, β::Float64, k::Real)::Float64
```

## Odd ratio

```
or_equality(p₀::Float64, p₁::Float64, α::Float64, β::Float64, k::Real)::Float64

or_equivalence(p₀::Float64, p₁::Float64, δ::Real, α::Float64, β::Float64, k::Real)::Float64

or_superiority(p₀::Float64, p₁::Float64, δ::Real, α::Float64, β::Float64, k::Real)::Float64
```

## McNemar's Equality test

```
mcnm(p10::Float64, p01::Float64, α::Float64, β::Float64)::Float64
```

# Power

## Means

### One Sample Means

```
one_mean_equality_pow(μ₀::Real, μ₁::Real, σ::Real, n::Int, α::Float64)::Float64

one_mean_equivalence_pow(μ₀::Real, μ₁::Real, σ::Real, δ::Real, n::Int, α::Float64)::Float64

one_mean_superiority_pow(μ₀::Real, μ₁::Real, σ::Real, δ::Real, n::Int, α::Float64)::Float64
```

### Two Sample Means

```
two_mean_equality_pow(μ₀::Real, μ₁::Real, σ::Real, n::Int, α::Float64, k::Real)::Float64

two_mean_equivalence_pow(μ₀::Real, μ₁::Real, σ::Real, δ::Real, n::Int, α::Float64, k::Real)::Float64

two_mean_superiority_pow(μ₀::Real, μ₁::Real, σ::Real, δ::Real, n::Int, α::Float64, k::Real)::Float64
```

## Compare Proportion

### One Sample proportion

```
one_proportion_equality_pow(p₀::Float64, p₁::Float64, n::Int, α::Float64)::Float64

one_proportion_equivalence_pow(p₀::Float64, p₁::Float64, δ::Real, n::Int, α::Float64)::Float64

one_proportion_superiority_pow(p₀::Float64, p₁::Float64, δ::Real, n, α::Float64)::Float64
```

### Two Sample Proportion

```
two_proportion_equality_pow(p₀::Float64, p₁::Float64, n::Int, α::Float64, k::Real)::Float64

two_proportion_equivalence_pow(p₀::Float64, p₁::Float64, δ::Real, n::Int, α::Float64, k::Real)::Float64

two_proportion_superiority_pow(p₀::Float64, p₁::Float64, δ::Real, n::Int, α::Float64, k::Real)::Float64
```

## Odd Ratio

```
or_equality_pow(p₀::Float64, p₁::Float64, n::Int, α::Float64, k::Real)::Float64

or_equivalence_pow(p₀::Float64, p₁::Float64, δ::Real, n::Int, α::Float64, k::Real)::Float64

or_superiority_pow(p₀::Float64, p₁::Float64, δ::Real, n::Int, α::Float64, k::Real)::Float64
```

## McNemar's Equality test

```
mcnm_pow(p10::Float64, p01::Float64, n::Int, α::Float64)::Float64
```
