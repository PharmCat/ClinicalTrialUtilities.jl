# Non-compartment analysis

NCA analysis based on following steps:

1. Loadind data to DataFrame;
2. Constructing subjects list with pkimport/pdimport;
3. Run NCA;
4. Exporting to DataFrame;
5. Descriptive statistics / HTML export.

### Description

```julia
include("ncatable.jl")
```
#### AUC

```math
AUC = \sum_{n=1}^N AUC_{n}
```

Where AUCn - partial AUC.

Linear trapezoidal rule

```math
AUC\mid_{t_1}^{t_2} = \delta t \times \frac{C_1 + C_2}{2}

AUMC\mid_{t_1}^{t_2} = \delta t \times \frac{t_1 \times C_1 + t_2 \times C_2}{2}
```

Logarithmic trapezoidal rule

```math
AUC\mid_{t_1}^{t_2} =   \delta t \times \frac{ C_2 - C_1}{ln(C_2/C_1)}

AUMC\mid_{t_1}^{t_2} = \delta t \times \frac{t_2 \times C_2 - t_1 \times C_1}{ln(C_2/C_1)} -  \delta t^2 \times \frac{ C_2 - C_1}{ln(C_2/C_1)^2}
```

Linear interpolation rule

```math
C_x = C_1 + \frac{(t_x-t_1)\times(C_2 - C_1)}{t_2 - t_1}
```

Logarithmic interpolation rule

```math
C_x = exp\left(ln(C_1) + \frac{(t_x-t_1)\times(ln(C_2) - ln(C_1))}{t_2 - t_1}\right)
```

#### \lambda_z - elimination constant

#### HL

```math
HL = ln(2) / \lambda_z
```

#### AUCinf

```math
AUCinf = AUClast + \frac{Clast}{\lambda_z}
```

#### AUMCinf

```math
AUMCinf =  AUMClast + \frac{tlast\times Clast}{\lambda_z} + \frac{Clast}{\lambda_z^2}
```

#### Accumulation index

```math
Accind = \frac{1}{1 - exp(-\lambda_z \tau)}
```

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
