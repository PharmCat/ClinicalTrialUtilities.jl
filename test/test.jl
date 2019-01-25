


println(" ---------------------------------- ")
println(" ---------   START TEST   --------- ")
println(" ---------------------------------- ")
@testset "sampleSize Test       " begin
    @test ceil(CTPSS.sampleSize(param="mean", type="ea", group="one", alpha=0.05, beta=0.2, sd=1, a=1.5, b=2, k=1)) == 32
    @test ceil(CTPSS.sampleSize(param="mean", type="ei", group="one", alpha=0.05, beta=0.2, sd=0.1, diff=0.05, a=2, b=2, k=1)) == 35
    @test ceil(CTPSS.sampleSize(param="mean", type="ns", group="one", alpha=0.05, beta=0.2, sd=1, diff=-0.5, a=1.5, b=2, k=1)) == 7
    @test ceil(CTPSS.sampleSize(param="mean", type="ea", group="two", alpha=0.05, beta=0.2, sd=10, a=5, b=10, k=1)) == 63
    @test ceil(CTPSS.sampleSize(param="mean", type="ei", group="two", alpha=0.05, beta=0.2, sd=10, diff=5, a=5, b=4, k=1)) == 108
    @test ceil(CTPSS.sampleSize(param="mean", type="ns", group="two", alpha=0.05, beta=0.2, sd=10, diff=5, a=5, b=5, k=1)) == 50
    @test ceil(CTPSS.sampleSize(param="prop", type="ea", group="one", alpha=0.05, beta=0.2, a=0.3, b=0.5)) == 50
    @test ceil(CTPSS.sampleSize(param="prop", type="ei", group="one", alpha=0.05, beta=0.2, diff=0.2, a=0.6, b=0.6)) == 52
    @test ceil(CTPSS.sampleSize(param="prop", type="ns", group="one", alpha=0.05, beta=0.2, diff=-0.1, a=0.3, b=0.5)) == 18
    @test ceil(CTPSS.sampleSize(param="prop", type="ea", group="two", alpha=0.05, beta=0.2, a=0.65, b=0.85)) == 70
    @test ceil(CTPSS.sampleSize(param="prop", type="ei", group="two", alpha=0.05, beta=0.2, diff=0.05, a=0.65, b=0.85)) == 136
    @test ceil(CTPSS.sampleSize(param="prop", type="ns", group="two", alpha=0.05, beta=0.2, diff=-0.1, a=0.85, b=0.65)) == 25
    @test ceil(CTPSS.sampleSize(param="or", type="ea",  alpha=0.05, beta=0.2, a=0.4, b=0.25)) == 156
    @test ceil(CTPSS.sampleSize(param="or", type="ei",  alpha=0.05, beta=0.2, diff=0.5, a=0.25, b=0.25)) == 366
    @test ceil(CTPSS.sampleSize(param="or", type="ns",  alpha=0.05, beta=0.2, diff=0.2, a=0.4, b=0.25)) == 242

    @test ceil(CTPSS.mcnm(0.45,0.05)) == 23

    @test round(sampleSize(param="mean", type="ns", group="two", alpha=0.05, beta=0.2, diff=1, sd=20, a=1, b=2), digits=12) ≈ 1236.511446403953
end
println(" ---------------------------------- ")
@testset "sampleSizeParam Test  " begin
    vals = CTPSS.ParamSet("mean","ea","one",0.05,0.2,1.0,1.5,2,1)
    @test ceil(CTPSS.sampleSizeParam(vals)) == 32
end
println(" ---------------------------------- ")
@testset "sampleSizeParam Errors" begin

    @test !CTPSS.sampleSize(param="mean", type="ea", group="one", alpha=0.05, beta=0.2, diff=1, sd=1, a=0, b=0, k=0)
    @test !CTPSS.sampleSize(param="mean", type="ea", group="one", alpha=2, beta=0.2, diff=1, sd=1, a=1, b=1, k=1)
    @test !CTPSS.sampleSize(param="mean", type="ea", group="one", alpha=-1, beta=0.2, diff=1, sd=1, a=1, b=1, k=1)
    @test !CTPSS.sampleSize(param="mean", type="ea", group="one", alpha=0.5, beta=1, diff=1, sd=1, a=1, b=1, k=1)
    @test !CTPSS.sampleSize(param="mean", type="ea", group="one", alpha=0.5, beta=0, diff=1, sd=1, a=1, b=1, k=1)
    @test !CTPSS.sampleSize(param="mean", type="ei", group="one", alpha=0.5, beta=0.2, diff=0, sd=1, a=1, b=1, k=1)
    @test !CTPSS.sampleSize(param="mean", type="ns", group="one", alpha=0.5, beta=0.2, diff=1, sd=0, a=1, b=1, k=1)
    @test !CTPSS.sampleSize(param="mean", type="ea", group="", alpha=0.5, beta=0.2, diff=1, sd=1, a=1, b=1, k=1)
    @test !CTPSS.sampleSize(param="mean", type="", group="one", alpha=0.5, beta=0.2, diff=1, sd=1, a=1, b=1, k=1)
    @test !CTPSS.sampleSize(param="", type="", group="one", alpha=0.5, beta=0.2, diff=1, sd=1, a=1, b=1, k=1)
    @test !CTPSS.sampleSize(param="prop", type="ea", group="one", alpha=0.05, beta=0.2, diff=1, sd=1, a=-1, b=0, k=1)
    @test !CTPSS.sampleSize(param="prop", type="ei", group="one", alpha=0.05, beta=0.2, diff=1, sd=1, a=0.4, b=2, k=1)
end
println(" ---------------------------------- ")

@testset "PowerTOST Test        " begin

    @test round(CTPSS.tfn(1,2), digits=8) ≈  0.07846821
    @test round(CTPSS.owensT(1,2), digits=8) ≈ 0.07846819
    @test round(CTPSS.owensTint2(1, 3, 20, 3), digits=7) ≈ 0.4839414
    @test round(CTPSS.owensQo(2,1,0.5,0.2;a=0), digits=9) ≈ 0.006781741
    @test round(CTPSS.owensQo(1,2,1,1;a=0), digits=6) ≈ 0.321429
    @test round(CTPSS.owensQo(3,2,1,Inf;a=0), digits=7) ≈ 0.7436299
    @test round(CTPSS.owensQ(4,100,40,0,Inf), digits=7) ≈ 0.9584071
    @test round(CTPSS.owensQ(4,100,30,0,0.8), digits=8) ≈ 0.02702275
    @test round(CTPSS.powerTOSTOwenQ(0.05,0.1,0.4,0.05,0.11,23), digits=8) ≈ 0.00147511
    @test round(CTPSS.approxPowerTOST(0.05,0.4,0.9,0.05,0.11,23), digits=12) ≈ 1.076964e-06
    @test CTPSS.approxPowerTOST(0.05,1.0,1.0,0.5,0.2,100) == 0
    @test round(CTPSS.approx2PowerTOST(0.05,0.1,1.0,0.5,0.2,1000), digits=7) ≈ 0.4413917
    @test round(CTPSS.powerTOST(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.2, n=20, design="2x2", method="owenq"), digits=7) ≈ 0.8346802
    @test round(CTPSS.powerTOST(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.2, n=10, design="2x2", method="nct"), digits=7) ≈ 0.4316618
    @test round(CTPSS.powerTOST(alpha=0.1, logscale=false, theta1=-0.1, theta2=0.1, theta0=0, cv=0.14, n=21, design="2x2", method="shifted"), digits=7) ≈ 0.6626132
    @test round(CTPSS.powerTOST(alpha=0.05, logscale=false, theta1=-0.1, theta2=0.1, theta0=0, cv=0.14, n=30, design="2x2", method="nct"), digits=7) ≈ 0.7079951
    @test round(CTPSS.powerTOST(alpha=0.0000001, logscale=false, theta1=-0.1, theta2=0.1, theta0=0, cv=1, n=10000, design="2x2", method="owenq"), digits=7) ≈ 0.9380914
    @test round(CTPSS.powerTOST(alpha=0.0001, logscale=false, theta1=-0.1, theta2=0.1, theta0=0, cv=1, n=3500, design="2x2", method="owenq"), digits=7) ≈ 0.3545904

    @test round(CTPSS.owensT(1,Inf), digits=8)   ≈  0.07932763
    @test round(CTPSS.owensT(-1,Inf), digits=8)  ≈ 0.07932763
    @test round(CTPSS.owensT(1,-Inf), digits=8)  ≈ -0.07932763
    @test round(CTPSS.owensT(-1,-Inf), digits=8) ≈ -0.07932763
    @test CTPSS.owensT(Inf, 1) == 0
end

@testset "PowerTOST Errors      " begin
    @test !CTPSS.powerTOST(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.2, n=20,  method="")


end
