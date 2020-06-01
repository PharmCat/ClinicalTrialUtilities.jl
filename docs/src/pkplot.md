### PK plots

#### pkplot
```@docs
ClinicalTrialUtilities.pkplot
```

#### Examples


```@example pkplots
using ClinicalTrialUtilities, DataFrames, CSV, Plots
pkdatapath = joinpath(dirname(pathof(ClinicalTrialUtilities)))*"/../test/csv/pkdata2.csv"
pkdata     = CSV.File(pkdatapath) |> DataFrame
pkds       = ClinicalTrialUtilities.pkimport(pkdata, [:Subject, :Formulation]; conc = :Concentration, time = :Time)
plot1      = pkplot(pkds[1], legend = false);
#savefig("pkplot1.svg"); # hide
plot2      = pkplot(pkds;  pagesort = [:Formulation], typesort = [:Subject])[1];
#savefig("pkplot2.svg"); nothing # hide
```

Plot for subject:

```@example pkplots
plot1
```

First plot for DataSet:

```@example pkplots
plot2
```
