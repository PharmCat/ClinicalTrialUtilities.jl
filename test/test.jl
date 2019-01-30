# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)


println(" ---------------------------------- ")
println(" ---------   START TEST   --------- ")
println(" ---------------------------------- ")
@testset "sampleSize Test       " begin
    @test ceil(ClinicalTrialUtilities.sampleSize(param="mean", type="ea", group="one", alpha=0.05, beta=0.2, sd=1, a=1.5, b=2, k=1)) == 32
    @test ceil(ClinicalTrialUtilities.sampleSize(param="mean", type="ei", group="one", alpha=0.05, beta=0.2, sd=0.1, diff=0.05, a=2, b=2, k=1)) == 35
    @test ceil(ClinicalTrialUtilities.sampleSize(param="mean", type="ns", group="one", alpha=0.05, beta=0.2, sd=1, diff=-0.5, a=1.5, b=2, k=1)) == 7
    @test ceil(ClinicalTrialUtilities.sampleSize(param="mean", type="ea", group="two", alpha=0.05, beta=0.2, sd=10, a=5, b=10, k=1)) == 63
    @test ceil(ClinicalTrialUtilities.sampleSize(param="mean", type="ei", group="two", alpha=0.05, beta=0.2, sd=10, diff=5, a=5, b=4, k=1)) == 108
    @test ceil(ClinicalTrialUtilities.sampleSize(param="mean", type="ns", group="two", alpha=0.05, beta=0.2, sd=10, diff=5, a=5, b=5, k=1)) == 50
    @test ceil(ClinicalTrialUtilities.sampleSize(param="prop", type="ea", group="one", alpha=0.05, beta=0.2, a=0.3, b=0.5)) == 50
    @test ceil(ClinicalTrialUtilities.sampleSize(param="prop", type="ei", group="one", alpha=0.05, beta=0.2, diff=0.2, a=0.6, b=0.6)) == 52
    @test ceil(ClinicalTrialUtilities.sampleSize(param="prop", type="ns", group="one", alpha=0.05, beta=0.2, diff=-0.1, a=0.3, b=0.5)) == 18
    @test ceil(ClinicalTrialUtilities.sampleSize(param="prop", type="ea", group="two", alpha=0.05, beta=0.2, a=0.65, b=0.85)) == 70
    @test ceil(ClinicalTrialUtilities.sampleSize(param="prop", type="ei", group="two", alpha=0.05, beta=0.2, diff=0.05, a=0.65, b=0.85)) == 136
    @test ceil(ClinicalTrialUtilities.sampleSize(param="prop", type="ns", group="two", alpha=0.05, beta=0.2, diff=-0.1, a=0.85, b=0.65)) == 25
    @test ceil(ClinicalTrialUtilities.sampleSize(param="or", type="ea",  alpha=0.05, beta=0.2, a=0.4, b=0.25)) == 156
    @test ceil(ClinicalTrialUtilities.sampleSize(param="or", type="ei",  alpha=0.05, beta=0.2, diff=0.5, a=0.25, b=0.25)) == 366
    @test ceil(ClinicalTrialUtilities.sampleSize(param="or", type="ns",  alpha=0.05, beta=0.2, diff=0.2, a=0.4, b=0.25)) == 242

    @test ceil(ClinicalTrialUtilities.mcnm(0.45,0.05)) == 23

    @test round(sampleSize(param="mean", type="ns", group="two", alpha=0.05, beta=0.2, diff=1, sd=20, a=1, b=2), digits=12) ≈ 1236.511446403953
end
println(" ---------------------------------- ")
@testset "sampleSizeParam Test  " begin
    vals = ClinicalTrialUtilities.ParamSet("mean","ea","one",0.05,0.2,1.0,1.5,2,1)
    @test ceil(ClinicalTrialUtilities.sampleSizeParam(vals)) == 32
end
println(" ---------------------------------- ")
@testset "sampleSizeParam Errors" begin

    @test !ClinicalTrialUtilities.sampleSize(param="mean", type="ea", group="one", alpha=0.05, beta=0.2, diff=1, sd=1, a=0, b=0, k=0)
    @test !ClinicalTrialUtilities.sampleSize(param="mean", type="ea", group="one", alpha=2, beta=0.2, diff=1, sd=1, a=1, b=1, k=1)
    @test !ClinicalTrialUtilities.sampleSize(param="mean", type="ea", group="one", alpha=-1, beta=0.2, diff=1, sd=1, a=1, b=1, k=1)
    @test !ClinicalTrialUtilities.sampleSize(param="mean", type="ea", group="one", alpha=0.5, beta=1, diff=1, sd=1, a=1, b=1, k=1)
    @test !ClinicalTrialUtilities.sampleSize(param="mean", type="ea", group="one", alpha=0.5, beta=0, diff=1, sd=1, a=1, b=1, k=1)
    @test !ClinicalTrialUtilities.sampleSize(param="mean", type="ei", group="one", alpha=0.5, beta=0.2, diff=0, sd=1, a=1, b=1, k=1)
    @test !ClinicalTrialUtilities.sampleSize(param="mean", type="ns", group="one", alpha=0.5, beta=0.2, diff=1, sd=0, a=1, b=1, k=1)

    @test !ClinicalTrialUtilities.sampleSize(param="mean", type="ea", group="", alpha=0.5, beta=0.2, diff=1, sd=1, a=1, b=1, k=1)
    @test !ClinicalTrialUtilities.sampleSize(param="mean", type="", group="one", alpha=0.5, beta=0.2, diff=1, sd=1, a=1, b=1, k=1)
    @test !ClinicalTrialUtilities.sampleSize(param="mean", type="", group="two", alpha=0.5, beta=0.2, diff=1, sd=1, a=1, b=1, k=1)
    @test !ClinicalTrialUtilities.sampleSize(param="", type="", group="one", alpha=0.5, beta=0.2, diff=1, sd=1, a=1, b=1, k=1)
    @test !ClinicalTrialUtilities.sampleSize(param="prop", type="", group="one", alpha=0.05, beta=0.2, diff=1, a=0.5, b=0.5, k=1)
    @test !ClinicalTrialUtilities.sampleSize(param="prop", type="", group="two", alpha=0.05, beta=0.2, diff=1, a=0.5, b=0.5, k=1)
    @test !ClinicalTrialUtilities.sampleSize(param="prop", type="ea", group="", alpha=0.05, beta=0.2, diff=1, a=0.5, b=0.5, k=1)

    @test !ClinicalTrialUtilities.sampleSize(param="prop", type="ea", group="one", alpha=0.05, beta=0.2, diff=1, sd=1, a=-1, b=0, k=1)
    @test !ClinicalTrialUtilities.sampleSize(param="prop", type="ei", group="one", alpha=0.05, beta=0.2, diff=1, sd=1, a=0.4, b=2, k=1)

    @test !ClinicalTrialUtilities.sampleSize(param="or", type="", group="",  diff=1, a=0.5, b=0.5, k=1)
end

println(" ---------------------------------- ")
@testset "sampleSizePower       " begin

    @test ClinicalTrialUtilities.oneSampleMeanEqualityP(1.5,2,1,32;alpha=0.05) ≈ 0.8074304194325561 ≈ ClinicalTrialUtilities.ctPower(param="mean", type="ea", group="one", a=1.5, b=2, sd=1,n=32, alpha=0.05)
    @test ClinicalTrialUtilities.oneSampleMeanEquivalenceP(2, 2, 0.1, 0.05, 35; alpha=0.05) ≈ 0.8108839754376387  ≈ ClinicalTrialUtilities.ctPower(param="mean", type="ei", group="one", a=2, b=2, sd=0.1, diff=0.05, n=35, alpha=0.05)
    @test ClinicalTrialUtilities.oneSampleMeanNSP(1.5, 2, 1, -.5, 7; alpha=0.05) ≈ 0.8415707712023641 ≈ ClinicalTrialUtilities.ctPower(param="mean", type="ns", group="one", a=1.5, b=2, sd=1, diff=-0.5, n=7, alpha=0.05)
    @test ClinicalTrialUtilities.twoSampleMeanEqualityP(5, 10, 10, 63; alpha=0.05, k=1) ≈ 0.8013023941055788 ≈ ClinicalTrialUtilities.ctPower(param="mean", type="ea", group="two", a=5, b=10, sd=10, n=63, alpha=0.05)
    @test ClinicalTrialUtilities.twoSampleMeanEquivalenceP(5, 4, 10, 5,  108; alpha=0.05, k=1) ≈ 0.80452354556060  ≈ ClinicalTrialUtilities.ctPower(param="mean", type="ei", group="two", a=5, b=4, sd=10, diff=5, n=108, alpha=0.05)
    @test ClinicalTrialUtilities.twoSampleMeanNSP(5, 5, 10, 5, 50; alpha=0.05, k=1) ≈ 0.8037819415575257 ≈ ClinicalTrialUtilities.ctPower(param="mean", type="ns", group="two", a=5, b=5, sd=10, diff=5, n=50, alpha=0.05)
    @test ClinicalTrialUtilities.oneProportionEqualityP(0.3, 0.5, 50; alpha=0.05) ≈ 0.8074304194325561 ≈ ClinicalTrialUtilities.ctPower(param="prop", type="ea", group="one", a=0.3, b=0.5, n=50, alpha=0.05)
    @test ClinicalTrialUtilities.oneProportionEquivalenceP(0.6, 0.6, 0.2, 52; alpha=0.05) ≈ 0.8060834056011101 ≈ ClinicalTrialUtilities.ctPower(param="prop", type="ei", group="one", a=0.6, b=0.6, diff=0.2, n=52, alpha=0.05)
    @test ClinicalTrialUtilities.oneProportionNSP(0.3, 0.5, -0.1, 18; alpha=0.05) ≈ 0.8161481827204281 ≈ ClinicalTrialUtilities.ctPower(param="prop", type="ns", group="one", a=0.3, b=0.5, diff=-0.1, n=18, alpha=0.05)
    @test ClinicalTrialUtilities.twoProportionEqualityP(0.65, 0.85, 70; alpha=0.05, k=1) ≈ 0.8019139068528 ≈ ClinicalTrialUtilities.ctPower(param="prop", type="ea", group="two", a=0.65, b=0.85, n=70, alpha=0.05)
    @test ClinicalTrialUtilities.twoProportionEquivalenceP(0.65, 0.85, 0.05, 136; alpha=0.05, k=1) ≈ 0.8033294052407269  ≈ ClinicalTrialUtilities.ctPower(param="prop", type="ei", group="two", a=0.65, b=0.85, diff=0.05, n=136, alpha=0.05)
    @test ClinicalTrialUtilities.twoProportionNSP(0.85, 0.65, -0.1, 25; alpha=0.05, k=1) ≈ 0.808599833380679 ≈ ClinicalTrialUtilities.ctPower(param="prop", type="ns", group="two", a=0.85, b=0.65, diff=-0.1, n=25, alpha=0.05)
    @test ClinicalTrialUtilities.orEqualityP(0.4, 0.25, 156; alpha=0.05, k=1) ≈ 0.8020239054864792 ≈ ClinicalTrialUtilities.ctPower(param="or", type="ea", a=0.4, b=0.25, n=156, alpha=0.05)
    @test ClinicalTrialUtilities.orEquivalenceP(0.25, 0.25, 0.5, 366; alpha=0.05, k=1, logdiff=true) ≈ 0.8008593380478983  ≈ ClinicalTrialUtilities.ctPower(param="or", type="ei", a=0.25, b=0.25, diff=0.5, n=366, alpha=0.05)
    @test ClinicalTrialUtilities.orNSP(0.4, 0.25, 0.2, 242; alpha=0.05, k=1, logdiff=true) ≈ 0.8007200876001626  ≈ ClinicalTrialUtilities.ctPower(param="or", type="ns", a=0.4, b=0.25, diff=0.2, n=242, alpha=0.05)

end
println(" ---------------------------------- ")

@testset "PowerTOST Test        " begin

    @test round(ClinicalTrialUtilities.tfn(1,2), digits=8) ≈  0.07846821
    @test round(ClinicalTrialUtilities.tfn(0.1,10), digits=8) > 0 #Not validated with PowerTOST
    @test round(ClinicalTrialUtilities.tfn(0.1,10E20), digits=8) > 0 #Not validated with PowerTOST

    @test round(ClinicalTrialUtilities.owensTint2(1, 3, 20, 3), digits=7) ≈ 0.4839414

    @test round(ClinicalTrialUtilities.owensQo(1,2,1,1;a=0), digits=6) ≈ 0.321429
    @test round(ClinicalTrialUtilities.owensQo(2,1,0.5,0.2;a=0), digits=9) ≈ 0.006781741
    @test round(ClinicalTrialUtilities.owensQo(4,2,1,1;a=0), digits=8) ≈ 0.03739024
    @test round(ClinicalTrialUtilities.owensQo(7,2,1,1;a=0), digits=9) ≈ 0.001888241
    @test round(ClinicalTrialUtilities.owensQo(3,2,1,Inf;a=0), digits=7) ≈ 0.7436299

    @test round(ClinicalTrialUtilities.owensQ(4,100,40,0,Inf), digits=7) ≈ 0.9584071
    @test round(ClinicalTrialUtilities.owensQ(1,1,1,0,Inf), digits=5) ≈  0.42202
    @test round(ClinicalTrialUtilities.owensQ(4,100,30,0,0.8), digits=8) ≈ 0.02702275
    @test round(ClinicalTrialUtilities.owensQ(1,100,40,0,1), digits=7) ≈ 0.3718607

    @test round(ClinicalTrialUtilities.powerTOSTOwenQ(0.05,0.1,0.4,0.05,0.11,23), digits=8) ≈ 0.00147511

    @test round(ClinicalTrialUtilities.approxPowerTOST(0.05,0.4,0.9,0.05,0.11,23), digits=12) ≈ 1.076964e-06
    @test       ClinicalTrialUtilities.approxPowerTOST(0.05,1.0,1.0,0.5,0.2,100) == 0

    @test round(ClinicalTrialUtilities.approx2PowerTOST(0.05,0.1,1.0,0.5,0.2,1000), digits=7) ≈ 0.4413917

    @test round(ClinicalTrialUtilities.powerTOST(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.2, n=20, design="2x2", method="owenq"), digits=7) ≈ 0.8346802
    @test round(ClinicalTrialUtilities.powerTOST(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.2, n=10, design="2x2", method="nct"), digits=7) ≈ 0.4316618
    @test round(ClinicalTrialUtilities.powerTOST(alpha=0.1, logscale=false, theta1=-0.1, theta2=0.1, theta0=0, cv=0.14, n=21, design="2x2", method="shifted"), digits=7) ≈ 0.6626132
    @test round(ClinicalTrialUtilities.powerTOST(alpha=0.05, logscale=false, theta1=-0.1, theta2=0.1, theta0=0, cv=0.14, n=30, design="2x2", method="nct"), digits=7) ≈ 0.7079951
    @test round(ClinicalTrialUtilities.powerTOST(alpha=0.0000001, logscale=false, theta1=-0.1, theta2=0.1, theta0=0, cv=1, n=10000, design="2x2", method="owenq"), digits=7) ≈ 0.9380914
    @test round(ClinicalTrialUtilities.powerTOST(alpha=0.0001, logscale=false, theta1=-0.1, theta2=0.1, theta0=0, cv=1, n=3500, design="2x2", method="owenq"), digits=7) ≈ 0.3545904
    @test round(powerTOST(alpha=0.00000005, logscale=false, theta1=-0.1, theta2=0.1, theta0=0, cv=1.5, n=20000, design="2x2", method="owenq"), digits=7) ≈ 0.8197361

    @test round(ClinicalTrialUtilities.owensT(1,Inf), digits=8)   ≈  0.07932763
    @test round(ClinicalTrialUtilities.owensT(-1,Inf), digits=8)  ≈ 0.07932763
    @test round(ClinicalTrialUtilities.owensT(1,-Inf), digits=8)  ≈ -0.07932763
    @test round(ClinicalTrialUtilities.owensT(-1,-Inf), digits=8) ≈ -0.07932763
    @test round(ClinicalTrialUtilities.owensT(1, 0.5), digits=8) ≈0.0430647
    @test round(ClinicalTrialUtilities.owensT(1,2), digits=8) ≈ 0.07846819
    @test       ClinicalTrialUtilities.owensT(Inf, 1) == 0
end

@testset "PowerTOST Errors      " begin
    @test !ClinicalTrialUtilities.powerTOST(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.2, n=20,  method="")
    @test !ClinicalTrialUtilities.powerTOST(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.2, n=20, design="2x2", method="mvt")

end
println(" ---------------------------------- ")
@testset "beSampleN Test        " begin

    n, p = ClinicalTrialUtilities.beSampleN(;theta0=0.95, theta1=0.8, theta2=1.25, cv=0.2, alpha=0.05, beta=0.2, logscale=true, method="owenq")
    @test n == 20 && round(p, digits=7) == 0.8346802

    n, p = ClinicalTrialUtilities.beSampleN(;theta0=0.95, theta1=0.8, theta2=1.25, cv=0.3, alpha=0.05, beta=0.2, logscale=true, method="owenq")
    @test n == 40 && round(p, digits=7) == 0.8158453

    n, p = ClinicalTrialUtilities.beSampleN(;theta0=1.0, theta1=0.8, theta2=1.25, cv=0.3, alpha=0.05, beta=0.1, logscale=true, method="owenq")
    @test n == 40 && round(p, digits=7) == 0.9095603

    n, p = ClinicalTrialUtilities.beSampleN(;theta0=1.05, theta1=0.8, theta2=1.25, cv=0.4, alpha=0.05, beta=0.15, logscale=true, method="nct")
    @test n == 74 && round(p, digits=7) == 0.8558178

    n, p = ClinicalTrialUtilities.beSampleN(;theta0=1.05, theta1=0.9, theta2=1.25, cv=0.4, alpha=0.05, beta=0.15, logscale=true, method="nct")
    @test n == 108 && round(p, digits=7) == 0.8506248

    n, p = ClinicalTrialUtilities.beSampleN(;theta0=1.05, theta1=0.8, theta2=1.2, cv=0.5, alpha=0.05, beta=0.2, logscale=true, method="nct")
    @test n == 158 && round(p, digits=7) == 0.8039191

    n, p = ClinicalTrialUtilities.beSampleN(;theta0=1.05, theta1=0.8, theta2=1.25, cv=0.8, alpha=0.05, beta=0.2, logscale=true, method="shifted")
    @test n == 210 && round(p, digits=7) == 0.8012471

    n, p = ClinicalTrialUtilities.beSampleN(;theta0=0, theta1=-0.2, theta2=0.2, cv=0.5, alpha=0.05, beta=0.2, logscale=false, method="owenq")
    @test n == 110 && round(p, digits=7) == 0.8074124

    n, p = ClinicalTrialUtilities.beSampleN(;theta0=0, theta1=-0.2, theta2=0.2, cv=2, alpha=0.05, beta=0.2, logscale=false, method="owenq")
    @test n == 1716 && round(p, digits=7) == 0.8005618

    n, p = ClinicalTrialUtilities.beSampleN(;theta0=0, theta1=-0.2, theta2=0.2, cv=2, alpha=0.001, beta=0.2, logscale=false, method="owenq")
    @test n == 3828 && round(p, digits=7) == 0.8001454

    n, p = ClinicalTrialUtilities.beSampleN(;theta0=0, theta1=-0.2, theta2=0.2, cv=2, alpha=0.01, beta=0.01, logscale=false, method="owenq")
    @test n == 4810 && round(p, digits=7) == 0.9900151

    n, p = beSampleN( cv=0.347)
    @test n == 52 && round(p, digits=7) == 0.8136415

end
