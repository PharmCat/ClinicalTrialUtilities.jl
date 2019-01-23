
using Test
import CTPSS.sampleSize

    @testset "Test" begin
        vals = CTPSS.ParamSet("mean","ea","one",0.05,0.2,1.0,1.5,2,1)
        @test ceil(CTPSS.sampleSizeParam(vals)) == 32
        #@test ceil(sampleSize(param="mean", type="ea", group="one", alpha=0.05, beta=0.2, sd=1, a=1.5, b=2, k=1)) == 32
        #@test ceil(sampleSize(param="mean", type="ei", group="one", alpha=0.05, beta=0.2, sd=0.1, diff=0.05, a=2, b=2, k=1)) == 35
    end
