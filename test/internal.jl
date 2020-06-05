
@testset "Internal functions test " begin
@testset "  tfn function        " begin
    @test ClinicalTrialUtilities.tfn(1.0,2.0) ≈  0.07846821 atol=1E-8
    @test ClinicalTrialUtilities.tfn(0.1,10.0) > 0 #Not validated with PowerTOST
    @test ClinicalTrialUtilities.tfn(0.1,10E20) > 0 #Not validated with PowerTOST
end
@testset "  owensTint2 function " begin
    @test round(ClinicalTrialUtilities.owenstint2(1.0, 3.0, 20.0, 3.0), digits=7) ≈ 0.4839414
end
@testset "  owensqo function    " begin
    @test ClinicalTrialUtilities.owensqo(1 ,2.0,1.0,1.0;a=0.0) ≈ 0.321429    atol=1E-6
    @test ClinicalTrialUtilities.owensqo(2 ,1.0,0.5,0.2;a=0.0) ≈ 0.006781741 atol=1E-9
    @test ClinicalTrialUtilities.owensqo(4 ,2.0,1.0,1.0;a=0.0) ≈ 0.03739024  atol=1E-8
    @test ClinicalTrialUtilities.owensqo(7 ,2.0,1.0,1.0;a=0.0) ≈ 0.001888241 atol=1E-9
    @test ClinicalTrialUtilities.owensqo(3 ,2.0,1.0,Inf;a=0.0) ≈ 0.7436299   atol=1E-7
end
@testset "  owensq  function    " begin
    @test ClinicalTrialUtilities.owensq(4 ,100.0,40.0,0.0,Inf) ≈ 0.9584071  atol=1E-7
    @test ClinicalTrialUtilities.owensq(1 ,1.0,1.0,0.0,Inf)    ≈ 0.42202    atol=1E-5
    @test ClinicalTrialUtilities.owensq(4 ,100.0,30.0,0.0,0.8) ≈ 0.02702275 atol=1E-8
    @test ClinicalTrialUtilities.owensq(1,100.0,40.0,0.0,1.0)  ≈ 0.3718607  atol=1E-7
    @test ClinicalTrialUtilities.owensq(4 ,100.0,40.0,0.0,Inf) ≈ 0.9584071  atol=1E-7
    @test ClinicalTrialUtilities.owensq(1 ,1.0,1.0,0.0,Inf)    ≈ 0.42202    atol=1E-5
    @test ClinicalTrialUtilities.owensq(4 ,100.0,30.0,0.0,0.8) ≈ 0.02702275 atol=1E-8
    @test ClinicalTrialUtilities.owensq(1,100.0,40.0,0.0,1.0)  ≈ 0.3718607  atol=1E-7
end
@testset "  powerTOSTOwenQ      " begin
    @test ClinicalTrialUtilities.powertost_owenq(0.05,0.1,0.4,0.05,0.11,23.0) ≈ 0.00147511 atol=1E-8
end
@testset "  approxPowerTOST     " begin
    @test ClinicalTrialUtilities.powertost_nct(0.05,0.4,0.9,0.05,0.11,23.0) ≈ 1.076964e-06 atol=1E-12
    @test ClinicalTrialUtilities.powertost_nct(0.05,1.0,1.0,0.5,0.2,100.0) == 0
end
@testset "  approx2PowerTOST    " begin
    @test ClinicalTrialUtilities.powertost_shifted(0.05,0.1,1.0,0.5,0.2,1000.0) ≈ 0.4413917 atol=1E-7
end
@testset "  owenst              " begin
    @test ClinicalTrialUtilities.owenst(1.0,Inf)   ≈ 0.07932763  atol=1E-8
    @test ClinicalTrialUtilities.owenst(-1.0,Inf)  ≈ 0.07932763  atol=1E-8
    @test ClinicalTrialUtilities.owenst(1.0,-Inf)  ≈ -0.07932763 atol=1E-8
    @test ClinicalTrialUtilities.owenst(-1.0,-Inf) ≈ -0.07932763 atol=1E-8
    @test ClinicalTrialUtilities.owenst(1.0, 0.5)  ≈ 0.0430647   atol=1E-8
    @test ClinicalTrialUtilities.owenst(1.0,2.0)   ≈ 0.07846819  atol=1E-8
    @test ClinicalTrialUtilities.owenst(Inf, 1.0)   == 0
end

@testset "  designProp          " begin
    d = ClinicalTrialUtilities.Design(:parallel)
    @test d.df(30) ≈ 28 && d.bkni ≈ 1.0 && d.sq ≈ 2
    d = ClinicalTrialUtilities.Design(:d2x2)
    @test d.df(31) ≈ 29 && d.bkni ≈ 0.5 && d.sq ≈ 2
    d = ClinicalTrialUtilities.Design(:d2x2x4)
    @test d.df(31) ≈ 89 && d.bkni ≈ 0.25 && d.sq ≈ 2
    d = ClinicalTrialUtilities.Design(:d2x4x4)
    @test d.df(31) ≈ 89 && d.bkni ≈ 0.0625 && d.sq ≈ 4
    d = ClinicalTrialUtilities.Design(:d2x2x3)
    @test d.df(31) ≈ 59 && d.bkni ≈ 0.375 && d.sq ≈ 2
    d = ClinicalTrialUtilities.Design(:d2x3x3)
    @test d.df(31) ≈ 59 && d.bkni ≈ 1/6 && d.sq ≈ 3

end

@testset "  Proportion         " begin
    p1 = ClinicalTrialUtilities.DiffProportion(ClinicalTrialUtilities.Proportion(0.2), 0.1)
    p2 = ClinicalTrialUtilities.DiffProportion(0.2, ClinicalTrialUtilities.Proportion(0.1))
    p3 = ClinicalTrialUtilities.DiffProportion(0.2, 0.1)
    @test ClinicalTrialUtilities.getval(p1) == ClinicalTrialUtilities.getval(p2) == ClinicalTrialUtilities.getval(p3)
    #@test ClinicalTrialUtilities.getval(p1 + p2) == (p1.x + p2.x)/(p1.n + p2.n)
end

@testset "  Means        " begin
    m1 = ClinicalTrialUtilities.Mean(10, 2, 100)
    m2 = ClinicalTrialUtilities.Mean(10, 2)
    m3 = ClinicalTrialUtilities.Mean(10)
    @test ClinicalTrialUtilities.getval(m1) == ClinicalTrialUtilities.getval(m2) == ClinicalTrialUtilities.getval(m3)
end

@testset "  ConfInt         " begin
    ci = ClinicalTrialUtilities.ConfInt(1, 3, 2, 0.05)
    @test ci[1] == 1
    @test ci[2] == 3
    @test ci[3] == 2
    @test ci[4] == 0.05
end

@testset "  DataSet         " begin
    pkds = ClinicalTrialUtilities.pkimport(pkdata, [:Subject, :Formulation]; conc = :Concentration, time = :Time)
    @test eltype(pkds) == eltype(pkds.data)
    a = pkds[:Subject => 1]
    b = pkds[Dict(:Subject => 1)]
    @test a === b

    b = findall(pkds, Dict(:Formulation => "T"))
    @test length(b) == 2

    deleteat!(pkds, Dict(:Formulation => "T"))
    @test length(pkds) == 1

    pkds = ClinicalTrialUtilities.pkimport(pkdata, [:Subject, :Formulation]; conc = :Concentration, time = :Time)
    deleteat!(pkds, 1)
    @test length(pkds) == 2

    pkds = ClinicalTrialUtilities.pkimport(pkdata, [:Subject, :Formulation]; conc = :Concentration, time = :Time)
    deleteat!(pkds, [1,2])
    @test length(pkds) == 1
end
end
