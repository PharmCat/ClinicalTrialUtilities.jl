#=
println(" ---------------------------------- ")
@testset "#14 Simulations         " begin

t      = ClinicalTrialUtilities.bepower(cv=0.2, n=20).task
result = ClinicalTrialUtilities.besim(t; nsim = 100, seed = 1234, rng = StableRNG(0))
@test result ≈ 0.8 atol=0.1
Base.show(io, t)

t      = ClinicalTrialUtilities.ctask(param=:mean, hyp=:ns, group=:two, alpha=0.05, n = 108, sd=10, diff=5, a=7, b=0, k=1)
result = ClinicalTrialUtilities.ctsim(t; nsim = 1000, seed = 1234, rng = StableRNG(0))
@test result == 0.439
Base.show(io, t)
#=
julia> ClinicalTrialUtilities.ctpower(t)
            Power Estimation
-----------------------------------------
  Parameter type: Mean Difference
  Design: Parallel
  Hypothesis: Superiority/Non-Inferiority
  H₀: A − B ≤ 5.0
  Hₐ: A − B > 5.0
  Alpha: 0.05
  N: 108
  A Mean(SD): 7 ± 10
  B Mean(SD): 0 ± 10
-----------------------------------------
Power: 0.431398
=#

 t      = ClinicalTrialUtilities.ctask(param=:mean, hyp=:ea, group=:two, alpha=0.05, n = 108, sd=5, diff=0, a=5, b=6, k=1)
 result = ClinicalTrialUtilities.ctsim(t; nsim = 1000, seed = 1234, rng = StableRNG(0))
 @test result == 0.323
 Base.show(io, t)

#=
julia> ClinicalTrialUtilities.ctpower(t)
            Power Estimation
-----------------------------------------
  Parameter type: Mean Difference
  Design: Parallel
  Hypothesis: Equality
  H₀: A - B = 0
  Hₐ: A - B ≠ 0
  Alpha: 0.05
  N: 108
  A Mean(SD): 5 ± 5
  B Mean(SD): 6 ± 5
-----------------------------------------
Power: 0.312274
=#

t      = ClinicalTrialUtilities.ctask(param=:mean, hyp=:ei, group=:two, alpha=0.05, n = 108, sd=10, diff=5, a=5, b=4, k=1)
result = ClinicalTrialUtilities.ctsim(t; nsim = 1000, seed = 1234, rng = StableRNG(0))
@test result == 0.884
Base.show(io, t)
#=
julia> ClinicalTrialUtilities.ctpower(t)
            Power Estimation
-----------------------------------------
  Parameter type: Mean Difference
  Design: Parallel
  Hypothesis: Equivalence
  H₀: |A − B| ≥ 5.0
  Hₐ: |A − B| < 5.0
  Alpha: 0.05
  N: 108
  A Mean(SD): 5 ± 10
  B Mean(SD): 4 ± 10
-----------------------------------------
Power: 0.902674
=#

#---

t      = ClinicalTrialUtilities.CTask(
ClinicalTrialUtilities.DiffProportion(ClinicalTrialUtilities.Proportion(30, 100), ClinicalTrialUtilities.Proportion(40, 100)),
ClinicalTrialUtilities.Parallel(),
ClinicalTrialUtilities.Superiority(-0.15, -0.15, 0.05),
ClinicalTrialUtilities.Power(100), 1.0)
result = ClinicalTrialUtilities.ctsim(t; nsim = 1000, method = :nhs, seed = 1234, rng = StableRNG(0))
@test result == 0.214
Base.show(io, t)

#=
julia> ClinicalTrialUtilities.ctpower(t)
            Power Estimation
-----------------------------------------
  Parameter type: Proportion Difference
  Design: Parallel
  Hypothesis: Superiority/Non-Inferiority
  H₀: A − B ≤ -0.15
  Hₐ: A − B > -0.15
  Alpha: 0.05
  N: 100
  A: 30/100
  B: 40/100
-----------------------------------------
Power: 0.192613
=#

t      = ClinicalTrialUtilities.CTask(
ClinicalTrialUtilities.DiffProportion(ClinicalTrialUtilities.Proportion(75, 100), ClinicalTrialUtilities.Proportion(80, 100)),
ClinicalTrialUtilities.Parallel(),
ClinicalTrialUtilities.Equivalence(-0.2, 0.2, 0.05),
ClinicalTrialUtilities.Power(132), 1.0)
result = ClinicalTrialUtilities.ctsim(t; nsim = 1000, seed = 1234, rng = StableRNG(0))
@test result == 0.897
Base.show(io, t)

#=
julia> ClinicalTrialUtilities.ctpower(t)
            Power Estimation
-----------------------------------------
  Parameter type: Proportion Difference
  Design: Parallel
  Hypothesis: Equivalence
  H₀: |A − B| ≥ 0.2
  Hₐ: |A − B| < 0.2
  Alpha: 0.05
  N: 132
  A: 75/100
  B: 80/100
-----------------------------------------
Power: 0.899422
=#

t      = ClinicalTrialUtilities.CTask(
ClinicalTrialUtilities.DiffProportion(ClinicalTrialUtilities.Proportion(30, 100), ClinicalTrialUtilities.Proportion(40, 100)),
ClinicalTrialUtilities.Parallel(),
ClinicalTrialUtilities.Equality(0, 0.05),
ClinicalTrialUtilities.Power(100), 1.0)
result = ClinicalTrialUtilities.ctsim(t; nsim = 1000, method = :mn, seed = 1234, rng = StableRNG(0))
@test result == 0.315
Base.show(io, t)

#=
julia> ClinicalTrialUtilities.ctpower(t)
            Power Estimation
-----------------------------------------
  Parameter type: Proportion Difference
  Design: Parallel
  Hypothesis: Equality
  H₀: A - B = 0
  Hₐ: A - B ≠ 0
  Alpha: 0.05
  N: 100
  A: 30/100
  B: 40/100
-----------------------------------------
Power: 0.319724
=#

#---

t      = ClinicalTrialUtilities.ctask(param=:prop, hyp=:ei, group=:two, alpha=0.05, n = 135, diff = 0.15,  a=0.45, b=0.4, k=1)
end

=#