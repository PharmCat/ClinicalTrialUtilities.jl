using CTPSS
using Test

    @testset "Sample Size Test" begin

    n = ceil(sampleSize(param="mean", type="ea", group="one", alpha=0.05, beta=0.2, sd=1, a=1.5, b=2, k=1))
    if n == 32 prinln "OK!"

    n = ceil(sampleSize(param="mean", type="ei", group="one", alpha=0.05, beta=0.2, sd=0.1, diff=0.05, a=2, b=2, k=1))
    if n == 35 prinln "OK!"

    end
