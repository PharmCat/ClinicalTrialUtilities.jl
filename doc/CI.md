# Confidence Intervals

CI submodule needs for clinical trial simulations with binomial outcomes, but can easily used for confidence interval evaluating in common practice.

## Functions

One proportion

oneProp(x, n; alpha, method)

One mean

oneMean(m, s, n, alpha; method)

Two proportions

twoProp(x1, n1, x2, n2; alpha, type, method)

Two means

twoMeans(m1, s1, n1, m2, s2, n2; alpha, type, method)

Cochran–Mantel–Haenszel confidence intervals.

cmh(data::DataFrame; a = :a, b = :b, c = :c, d = :d, type = :diff, alpha = 0.05, method)

Return result as ConfInt(lower::Float64, upper::Float64, estimate::Float64)

### oneProp

Confidence interval calculation for one proportion.

```
oneProp(x::Int, n::Int; alpha=0.05, method=[:wilson|:wilsoncc|:cp|:soc|:blaker|:arcsine|:wald])
```

**x** - number of positive outcomes;

**n** - number of subjects;

**alpha** - alpha value (all CI computed as two-sided);

**method** - computation method;

- wilson - Wilson's confidence interval (CI) for a single proportion (wilson score);
- wilsoncc - Wilson's CI with continuity correction (CC);
- cp - Clopper-Pearson exact CI;
- soc - SOC: Second-Order corrected CI;
- blaker - Blaker exact CI for discrete distributions;
- arcsine - Arcsine CI;
- wald - Wald CI without CC;

### twoProp

Confidence interval calculation for two proportion.

```
 twoProp(x1::Int, n1::Int, x2::Int, n2::Int; alpha=0.05, type=[:diff|:rr|:or], method::Symbol)
 ```

 **x1** - number of positive outcomes in group 1;

 **n1** - number of subjects in group 1;

 **x2** - number of positive outcomes in group 2;

 **n2** - number of subjects in group 2;

 **alpha** - alpha value (all CI computed as two-sided);

 **type** - comparation type: difference, relative risk, odd ratio;

- diff - Difference;
- rr - Risk ratio;
- or - Odd ratio;

 **method** - computation method;

Difference:

- nhs - Newcombes Hybrid (wilson) Score interval for the difference of proportions
- nhscc - Newcombes Hybrid Score CC
- ac - Agresti-Caffo interval for the difference of proportions
- mn - Method of Mee 1984 with Miettinen and Nurminen modification
- mee - Mee maximum likelihood method
- wald - Wald CI
- waldcc - Wald CC CI

Risk ratio

- cli - Crude log interval
- mover - Method of variance estimates recovery

Odd ratio

- mn - Miettinen-Nurminen CI
- woolf - Woolf logit CI
- awoolf - Adjusted Woolf interval (Gart adjusted logit)

### oneMean

Confidence interval for one mean.

...


### twoMeans

Confidence interval for mean difference.

...


### cmh

Cochran–Mantel–Haenszel confidence intervals.


> cmh(data::DataFrame; a = :a, b = :b, c = :c, d = :d, alpha = 0.05, type = [**:diff**|:or|:rr], method = :default, logscale = false)::ConfInt


 **data**- dataframe with 4 columns, each line represent 2X2 table

**a** **b** **c** **d** - dataframe table names (number of subjects in 2X2 table):

|         | outcome 1 | outcome 2 |
|---------|-----------|-----------|
| group 1 |     a     |     b     |
| group 2 |     c     |     d     |

**alpha** - alpha value (all CI computed as two-sided);

**type** - estimation type
- diff - risk difference (default)
- or - odd ratio
- rr - risk ratio

**logscale** - return CI in log scale.

- false (default)
- true

For example, let we have 2 sources:

1:

|         | outcome 1 | outcome 2 |
|---------|-----------|-----------|
| group 1 |     8     |     98    |
| group 2 |     5     |    115    |

2:

|         | outcome 1 | outcome 2 |
|---------|-----------|-----------|
| group 1 |     22    |     76    |
| group 2 |     16    |     69    |


Dataframe construction:

```
data = DataFrame(a = Int[], b = Int[], c = Int[], d = Int[])
  push!(data, (8, 98, 5, 115))
  push!(data, (22, 76, 16, 69))
```

For risk difference use:

```
ci = ClinicalTrialUtilities.CI.cmh(data, alpha = 0.1)
```

## Examples:

```
ClinicalTrialUtilities.CI.oneProp(81, 263, alpha=0.05, method=:wilson)

ClinicalTrialUtilities.CI.twoProp(7, 34, 1, 34; alpha=0.05, type=:diff, method=:nhs)

ClinicalTrialUtilities.CI.twoProp(2, 14, 1, 11; alpha=0.05, type=:or, method=:woolf)

ClinicalTrialUtilities.CI.twoMeans(30, 10, 30, 40, 12, 35, alpha=0.05, method=:ev)

```

## References:

Wilson, E.B. (1927) Probable inference, the law of succession, and statistical inference. J. Amer. Stat. Assoc22, 209–212.

Clopper, C. and Pearson, E.S. (1934) The use of confidence or fiducial limits illustrated in the case of the binomial. Biometrika 26, 404–413.

Blaker, H. (2000). Confidence curves and improved exact confidence intervals for discrete distributions. Canadian Journal of Statistics28 (4), 783–798.

T. Tony Cai One-sided confidence intervals in discrete distributions. doi:10.1016/j.jspi.2004.01.00.

Newcombe, R. G. (1998) Two-sided confidence intervals for the single proportion: comparison of seven methods. Statistics in Medicine. 17 (8): 857–872. doi:10.1002/(SICI)1097-0258(19980430)17:8<857::AID-SIM777>3.0.CO;2-E. PMID 9595616.

Miettinen O. S., Nurminen M. (1985) Comparative analysis of two rates. Statistics in Medicine4,213–226

Lawson, R (2005):Smallsample confidence intervals for the odds ratio.  Communication in Statistics Simulation andComputation, 33, 1095-1113.

Woolf, B. (1955). On estimating the relation between blood group and disease. Annals of human genetics, 19(4):251-253.

Gart, JJand Nam, J (1988): Approximate interval estimation of the ratio of binomial parameters: Areview and corrections for skewness. Biometrics 44, 323-338.

Donner, A. and Zou, G. (2012). Closed-form confidence intervals for functions of the normal mean and standard deviation. Statistical Methods in Medical Research, 21(4):347-359.

Newcombe RG (1998): Interval Estimation for the Difference Between Independent Proportions: Comparison of Eleven Methods. Statistics in Medicine 17, 873-890.

Agresti A, Caffo B., “Simple and effective confidence intervals for proportions and differences of proportions result from adding two successes and two failures”, American Statistician 54: 280–288 (2000)

Mee RW (1984) Confidence bounds for the difference between two probabilities,Biometrics40:1175-1176

Rothman, K. J., Greenland, S., & Lash, T. L. (2008). Modern epidemiology (3rd ed.). Philadelphia: Lippincott Williams & Wilkins.
