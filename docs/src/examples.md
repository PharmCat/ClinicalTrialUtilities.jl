## Confidence intervals

```
propci(81, 263, alpha=0.05, method=:wilson)

diffpropci(7, 34, 1, 34; alpha=0.05, method=:nhs)

orpropci(2, 14, 1, 11; alpha=0.05, method=:woolf)

rrpropci(2, 14, 1, 11; alpha=0.05)

diffmeanci(30, 10, 30, 40, 12, 35, alpha=0.05, method=:ev)

```

## Sample size

### Equivalence for two means

```
ctsamplen(param=:mean, type=:ei, group=:two, diff=0.3, sd=1, a=0.3, b=0.5)
```

### One proportion equality

```
ctsamplen(param=:prop, type=:ea, group=:one, a=0.3, b=0.5)
```

### Odd ratio non-inferiority

```
ctsamplen(param=:or, type=:ns, diff=-0.1, a=0.3, b=0.5, k=2)
```

### Odd ratio equality

```
ctsamplen(param=:or, type=:ea, a=0.3, b=0.5, k=2)
```

### Bioequivalence

```
besamplen(alpha=0.05,  theta1=0.8, theta2=1.25, theta0=0.95, cv=0.15, method=:owenq)

besamplen(cv=0.20, method=:nct)

besamplen(cv=0.347, design=:parallel)

besamplen(cv=0.40)

besamplen(cv=0.347, design=:d2x2x4, method=:nct)
```

#Power

### Equality for one mean

```
ctpower(param=:mean, type=:ea, group=:one, a=1.5, b=2, sd=1,n=32, alpha=0.05)
```

### Bioequivalence

```
#Bioequivalence power for 2x2 design, default method - OwensQ
bepower(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.2, n=20, design=:d2x2, method=:owenq)

#The same
bepower(alpha=0.05, cv=0.2, n=20, design=:d2x2)

#Bioequivalence power for cv 14%, 21 subjects, default OwensQ method, logscale false
bepower(alpha=0.1, logscale=false, theta1=-0.1, theta2=0.1, theta0=0, cv=0.14, n=21)

#Bioequivalence power for cv 14%, 21 subjects, shifted method, logscale false
bepower(alpha=0.1, logscale=false, theta1=-0.1, theta2=0.1, theta0=0, cv=0.14, n=21, method=:shifted)

#Simple notations
bepower(cv=0.4, n=35, design=:d2x4x4)
bepower(cv=0.14, n=21)
```

## Utilities

### CV from CI

```
cvfromci(;alpha = 0.05, theta1 = 0.9, theta2 = 1.25, n=30, design=:d2x2x4)
```

### Polled CV

```
data = DataFrame(cv = Float64[], df = Int[])
push!(data, (0.12, 12))
push!(data, (0.2, 20))
push!(data, (0.25, 30))
pooledcv(data; cv=:cv, df=:df, alpha=0.05, returncv=true)
```
