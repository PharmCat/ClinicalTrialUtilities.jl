println(" ---------------------------------- ")
@testset "  Simulations           " begin

    t      = ClinicalTrialUtilities.bepower(cv=0.2, n=20).task
    result = ClinicalTrialUtilities.besim(t; nsim = 100, seed = 1234)
    @test result == 0.83
    Base.show(io, t)

    t      = ClinicalTrialUtilities.CTask(
    ClinicalTrialUtilities.DiffProportion(ClinicalTrialUtilities.Proportion(30, 100), ClinicalTrialUtilities.Proportion(40, 100)),
    ClinicalTrialUtilities.Parallel(),
    ClinicalTrialUtilities.Superiority(-0.15, -0.15, 0.05),
    ClinicalTrialUtilities.Power(100), 1.0)
    result = ClinicalTrialUtilities.ctsim(t; nsim = 1000, method = :nhs, seed = 1234)
    @test result == 0.169
    Base.show(io, t)

    t      = ClinicalTrialUtilities.ctask(param=:mean, hyp=:ei, group=:two, alpha=0.05, n = 108, sd=10, diff=5, a=5, b=4, k=1)
    result = ClinicalTrialUtilities.ctsim(t; nsim = 1000, seed = 1234)
    @test result == 0.809
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
    Power: 0.804524
    =#



end
