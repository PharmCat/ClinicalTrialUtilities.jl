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

Return result as ConfInt(lower::Float64, upper::Float64, estimate::Float64)

### oneProp

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
- soc - SOC: Second-Order corrected CI
- blaker - Blaker exact CI for discrete distributions;
- arcsine - Arcsine CI;
- wald - Wald CI without CC;

### twoProp

```
 twoProp(x1::Int, n1::Int, x2::Int, n2::Int; alpha=0.05, type=[:diff|:rr|:or], method::Symbol)
 ```

 **x1** - number of positive outcomes in group 1;

 **n1** - number of subjects in group 1;

 **x2** - number of positive outcomes in group 2;

 **n2** - number of subjects in group 2;

 **alpha** - alpha value (all CI computed as two-sided);

 **type** - comparation type: difference, relative risk, odd ratio;

 **method** - computation method;


## References:

Wilson, E.B. (1927) Probable inference, the law of succession, and statistical inference. J. Amer. Stat. Assoc22, 209–212.

Clopper, C. and Pearson, E.S. (1934) The use of confidence or fiducial limits illustrated in the case of the binomial. Biometrika 26, 404–413.

Blaker, H. (2000). Confidence curves and improved exact confidence intervals for discrete distributions. Canadian Journal of Statistics28 (4), 783–798.

T. Tony Cai One-sided confidence intervals in discrete distributions. doi:10.1016/j.jspi.2004.01.00.

Newcombe, R. G. (1998) Two-sided confidence intervals for the single proportion: comparison of seven methods. Statistics in Medicine. 17 (8): 857–872. doi:10.1002/(SICI)1097-0258(19980430)17:8<857::AID-SIM777>3.0.CO;2-E. PMID 9595616.
