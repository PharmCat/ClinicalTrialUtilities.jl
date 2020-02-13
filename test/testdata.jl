#Simple PK DataFrame
pkdata         = CSV.File(path*"/csv/pkdata.csv") |> DataFrame
#Simple PK DataFrame
nanpkdata         = CSV.File(path*"/csv/nanpk.csv") |> DataFrame
#Simple PD DataFrame
pddata         = CSV.File(path*"/csv/pddata.csv") |> DataFrame
# Multiple subjects PK DataFrame
pkdata2        = CSV.File(path*"/csv/pkdata2.csv") |> DataFrame
#Glucose2
#Pinheiro, J. C. and Bates, D. M. (2000), Mixed-Effects Models in S and S-PLUS, Springer, New York. (Appendix A.10)
#Hand, D. and Crowder, M. (1996), Practical Longitudinal Data Analysis, Chapman and Hall, London.
glucose2       = CSV.File(path*"/csv/glucose2.csv") |> DataFrame
#Simple urine PK DataFrame
upkdata         = CSV.File(path*"/csv/upkdata.csv") |> DataFrame
#Simple frequency dataset
freqdat        = CSV.File(path*"/csv/freqdat.csv") |> DataFrame
#Small dataset with negative and zero value
negdat         = CSV.File(path*"/csv/negdat.csv") |> DataFrame
#Dataset foe descriptive statistics check
descriptivedat = CSV.File(path*"/csv/descriptivedat.csv") |> DataFrame
