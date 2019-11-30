# Non-compartment analysis

NCA analysis based on following steps:

1. Loadind data to DataFrame;
2. Constructing subjects list with pkimport/pdimport;
3. Run NCA;
4. Exporting to DataFrame;
5. Descriptive statistics / HTML export.


### nca!
```@docs
ClinicalTrialUtilities.nca!
```

### Scenario

1 Loading DataFrame

```julia
using CSV, DataFrames, ClinicalTrialUtilities
pkdata2 = CSV.File("pkdata.csv") |> DataFrame
```
2 Subject list

```julia
  pkds = pkimport(pkdata2, [:Subject, :Formulation]; time = :Time, conc = :Concentration)
```

3 NCA  analysis with default settings

```julia
  pk   = nca!(pkds)
```

4 Exporting

```julia
  ncadf   = DataFrame(pk; unst = true)
```

5 Descriptive statistics

```julia
  ds   = ClinicalTrialUtilities.descriptive(ncadf, stats = [:n, :mean, :sd], sort = [:Formulation])
```

6 Exporting  

```julia
  dsdf   = ClinicalTrialUtilities.DataFrame(ds; unst = true)
```
