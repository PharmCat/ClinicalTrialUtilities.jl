println(" ---------------------------------- ")
@testset "  Simulations           " begin

    t      = ClinicalTrialUtilities.bepower(cv=0.2, n=20).task
    result = ClinicalTrialUtilities.ctsim(t; nsim = 100, seed = 1234)
    @test result == 0.83
    Base.show(io, t)

    t      = ClinicalTrialUtilities.CTask(
    ClinicalTrialUtilities.DiffProportion(ClinicalTrialUtilities.Proportion(30, 100), ClinicalTrialUtilities.Proportion(40, 100)),
    ClinicalTrialUtilities.Parallel(),
    ClinicalTrialUtilities.Superiority(-0.15, -0.15),
    ClinicalTrialUtilities.Power(100), 0.05, 1.0)
    result = ClinicalTrialUtilities.ctsim(t, nsim = 1000, method = :nhs, seed = 1234)
    @test result == 0.169
    Base.show(io, t)

    t      = ClinicalTrialUtilities.CTask(
    ClinicalTrialUtilities.DiffMean(ClinicalTrialUtilities.Mean(10, 2), ClinicalTrialUtilities.Mean(12, 2)),
    ClinicalTrialUtilities.Parallel(),
    ClinicalTrialUtilities.Equivalence(-1.15, 1.15),
    ClinicalTrialUtilities.Power(100), 0.05, 1.0)
    result = ClinicalTrialUtilities.ctsim(t, nsim = 1000, seed = 1234)
    @test result == 0.691
    Base.show(io, t)
end
