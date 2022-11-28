
#Glucose2
#Pinheiro, J. C. and Bates, D. M. (2000), Mixed-Effects Models in S and S-PLUS, Springer, New York. (Appendix A.10)
#Hand, D. and Crowder, M. (1996), Practical Longitudinal Data Analysis, Chapman and Hall, London.
glucose2       = CSV.File(path*"/csv/glucose2.csv") |> DataFrame
#theo
#https://cran.r-project.org/web/packages/PKNCA/vignettes/Example-theophylline.html
theodata       = CSV.File(path*"/csv/theo.csv") |> DataFrame

#Simple frequency dataset
freqdat        = CSV.File(path*"/csv/freqdat.csv") |> DataFrame
#Small dataset with negative and zero value
negdat         = CSV.File(path*"/csv/negdat.csv") |> DataFrame
#Dataset foe descriptive statistics check
descriptivedat = CSV.File(path*"/csv/descriptivedat.csv") |> DataFrame
#Meta-analysis
metadf = CSV.File(path*"/csv/meta.csv") |> DataFrame
