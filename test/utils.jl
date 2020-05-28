println(" ---------------------------------- ")
@testset "#5  Utilities           " begin

@testset "  ci2cv Test            " begin
    #cvms = ClinicalTrialUtilities.cvfromci(;alpha = 0.05, theta1 = 0.9, theta2 = 1.25, n=30, design=:d2x2x4, cvms=true)
    #@test cvms[1] ≈ 0.583175066140736 && cvms[2] ≈ 0.29273913226894244
    @test ClinicalTrialUtilities.cvfromci(;alpha = 0.05, theta1 = 0.9, theta2 = 1.25, n=30, design=:d2x2x4, mso=true) ≈ 0.29273913226894244
    @test ClinicalTrialUtilities.cvfromci(;alpha = 0.05, theta1 = 0.9, theta2 = 1.25, n=30) ≈ 0.387417014838382
    @test ClinicalTrialUtilities.cvfromci(;alpha = 0.05, theta1 = 0.8, theta2 = 1.25, n=30) ≈ 0.5426467811605913
    @test ClinicalTrialUtilities.cvfromci(;alpha = 0.05, theta1 = 1.01, theta2 = 1.21, n=31, design=:d2x2) ≈ 0.21151405971696524

    @test ClinicalTrialUtilities.cvfromci(0.9, 1.25, 30) ≈ 0.387417014838382
end

@testset "  pooledcv            " begin
    data = DataFrame(cv = Float64[], df = Int[])
    push!(data, (0.12, 12))
    push!(data, (0.2, 20))
    push!(data, (0.25, 30))
    ci = ClinicalTrialUtilities.pooledcv(data; cv="cv", df="df")
    @test ci.lower    ≈ 0.18145259424967664 atol=1E-15
    @test ci.upper    ≈ 0.2609307413637307 atol=1E-15
    @test ci.estimate ≈ 0.21393949168210136 atol=1E-15

    ci = ClinicalTrialUtilities.pooledcv(data; cv="cv", df="df", returncv=false)
    @test ci.lower    ≈ 0.032394625994562894 atol=1E-15
    @test ci.upper    ≈ 0.06586718662236608 atol=1E-15
    @test ci.estimate ≈ 0.04475355764465427 atol=1E-15

    ci = ClinicalTrialUtilities.pooledcv([0.12, 0.2, 0.25], [12, 20, 30])
    @test ci.lower    ≈ 0.18145259424967664 atol=1E-15
    @test ci.upper    ≈ 0.2609307413637307 atol=1E-15
    @test ci.estimate ≈ 0.21393949168210136 atol=1E-15

    ci = ClinicalTrialUtilities.pooledcv([0.12, 0.2, 0.25], [14, 22, 32], [:d2x2, :d2x2, :d2x2])
    @test ci.lower    ≈ 0.18145259424967664 atol=1E-15
    @test ci.upper    ≈ 0.2609307413637307 atol=1E-15
    @test ci.estimate ≈ 0.21393949168210136 atol=1E-15
end

@test ClinicalTrialUtilities.cvfromsd(ClinicalTrialUtilities.sdfromcv(0.2)) ≈ 0.2

@test ClinicalTrialUtilities.cvfromvar(ClinicalTrialUtilities.varfromcv(0.2)) ≈ 0.2
end
