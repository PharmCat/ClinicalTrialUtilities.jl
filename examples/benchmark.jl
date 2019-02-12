using BenchmarkTools
#using ClinicalTrialUtilities

b = @benchmark ClinicalTrialUtilities.sampleSize(param="mean", type="ea", group="one", alpha=0.05, beta=0.2, sd=1, a=1.5, b=2, k=1)

display(b[:])

b = @benchmark ClinicalTrialUtilities.owensQ(4,100.0,30.0,0.0,0.8)

display(b[:])

b = @benchmark ClinicalTrialUtilities.powerTOST(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.2, n=20, design="2x2", method="owenq")

display(b[:])

b = @benchmark ClinicalTrialUtilities.beSampleN(;theta0=1.0, theta1=0.8, theta2=1.25, cv=0.3, alpha=0.05, beta=0.1, logscale=true, method="owenq")

display(b[:])
