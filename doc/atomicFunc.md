# Sample Size

## Means

### One Sample Means

m0 = μ0 - null H; m1 = μT - Test Group

oneSampleMeanEquality(m0, m1, sd; alpha=0.05, beta=0.2)::Float64

oneSampleMeanEquivalence(m0, m1, sd, diff; alpha=0.05, beta=0.2)

oneSampleMeanNS(m0, m1, sd, diff; alpha=0.05, beta=0.2) #Non-inferiority / Superiority

### Two Sample Means

m0 = μA - Group A; m1 = μB - Group B

twoSampleMeanEquality(m0, m1, sd; alpha=0.05, beta=0.2, k=1)

twoSampleMeanEquivalence(m0, m1, sd, diff; alpha=0.05, beta=0.2, k=1)

twoSampleMeanNS(m0, m1, sd, diff; alpha=0.05, beta=0.2, k=1) #Non-inferiority / Superiority

## Compare Proportion

### One Sample Proportion

oneProportionEquality(p0, p1; alpha=0.05, beta=0.2)

oneProportionEquivalence(p0, p1, diff; alpha=0.05, beta=0.2)

oneProportionNS(p0, p1, diff; alpha=0.05, beta=0.2)

### Two Sample Proportion

twoProportionEquality(p0, p1; alpha=0.05, beta=0.2, k=1)

twoProportionEquivalence(p0, p1, diff; alpha=0.05, beta=0.2, k=1)

twoProportionNS(p0, p1, diff; alpha=0.05, beta=0.2, k=1)

## Odd ratio

orEquality(p0, p1; alpha=0.05, beta=0.2, k=1)

orEquivalence(p0, p1, diff; alpha=0.05, beta=0.2, k=1, logdiff=true)

orNS(p0, p1, diff; alpha=0.05, beta=0.2, k=1, logdiff=true)


## McNemar's Equality test

mcnm(p10, p01; alpha=0.05, beta=0.2)


# Power

## Means

### One Sample Means

oneSampleMeanEqualityP(m0,m1,sd, n; alpha=0.05)

oneSampleMeanEquivalenceP(m0, m1, sd, diff, n; alpha=0.05)

oneSampleMeanNSP(m0, m1, sd, diff, n; alpha=0.05) #Non-inferiority / Superiority

### Two Sample Means

twoSampleMeanEqualityP(m0, m1, sd, n; alpha=0.05, k=1)

twoSampleMeanEquivalenceP(m0, m1, sd, diff, n; alpha=0.05, k=1)

twoSampleMeanNSP(m0, m1, sd, diff, n; alpha=0.05, k=1) #Non-inferiority / Superiority

## Compare Proportion

### One Sample proportion

oneProportionEqualityP(p0, p1, n; alpha=0.05)

oneProportionEquivalenceP(p0, p1, diff, n; alpha=0.05)

oneProportionNSP(p0, p1, diff, n; alpha=0.05)

### Two Sample Proportion

twoProportionEqualityP(p0, p1, n; alpha=0.05, k=1)

twoProportionEquivalenceP(p0, p1, diff, n; alpha=0.05, k=1)

twoProportionNSP(p0, p1, diff, n; alpha=0.05, k=1)

## Odd Ratio

orEqualityP(p0, p1, n; alpha=0.05, k=1)

orEquivalenceP(p0, p1, diff, n; alpha=0.05, k=1, logdiff=true)

orNSP(p0, p1, diff, n; alpha=0.05, k=1, logdiff=true)

## McNemar's Equality test

mcnmP(p10, p01, n; alpha=0.05)
