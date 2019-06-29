# Pharmacokinetics

Pharmacokinetics module.

## Functions

```
nca(data; conc=:Concentration, time=:Time, sort = NaN, calcm = :lint)
```

**conc**::Symbol - concentration column of dataframe data;

**time**::Symbol - time column of dataframe data;

**sort**::Symbol[] - sorting columns in dataframe (NaN by default in only one profile);

**calcm**::Symbol - calculation method:

- :lint - linear trapezoidal;
- :logt - linear trapezoidal before tmax and linear up log down after;
- :lulg - linear up log down: if c2 >= c1 - linear, else log;


### Simple example

```
data = DataFrame(Concentration = Float64[], Time = Float64[], Subject = String[], Formulation = String[])
  push!(data, (0.0, 0, "1", "T"))
  push!(data, (0.2, 1, "1", "T"))
  push!(data, (0.3, 2, "1", "T"))
  push!(data, (0.4, 3, "1", "T"))
  push!(data, (0.3, 4, "1", "T"))
  push!(data, (0.2, 5, "1", "T"))
  push!(data, (0.1, 6, "1", "T"))

pk = ClinicalTrialUtilities.PK.nca(data)
```

### Example with sorting

```
data = DataFrame(Concentration = Float64[], Time = Float64[], Subject = String[], Formulation = String[])
  push!(data, (0.0, 0, "1", "T"))
  push!(data, (0.2, 1, "1", "T"))
  push!(data, (0.3, 2, "1", "T"))
  push!(data, (0.4, 3, "1", "T"))
  push!(data, (0.3, 4, "1", "T"))
  push!(data, (0.2, 5, "1", "T"))
  push!(data, (0.1, 6, "1", "T"))
  push!(data, (0.0, 0, "2", "T"))
  push!(data, (0.3, 1, "2", "T"))
  push!(data, (0.4, 2, "2", "T"))
  push!(data, (0.5, 3, "2", "T"))
  push!(data, (0.3, 4, "2", "T"))
  push!(data, (0.1, 5, "2", "T"))
  push!(data, (0.05, 6, "2", "T"))
  push!(data, (0.1, 0, "1", "R"))
  push!(data, (0.2, 1, "1", "R"))
  push!(data, (0.9, 2, "1", "R"))
  push!(data, (0.8, 3, "1", "R"))
  push!(data, (0.3, 4, "1", "R"))
  push!(data, (0.2, 5, "1", "R"))
  push!(data, (0.1, 6, "1", "R"))
pk = ClinicalTrialUtilities.PK.nca(data; sort=[:Formulation, :Subject], calcm = :logt)

```
