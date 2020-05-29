
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
end
