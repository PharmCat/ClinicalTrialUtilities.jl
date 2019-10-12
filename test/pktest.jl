println(" ---------------------------------- ")
@testset "  PK                  " begin


    #pk = ClinicalTrialUtilities.PK.nca(data; conc = :Concentration, sort=[:Formulation, :Subject])
    #@test pk.result.AUCinf[1] ≈ 1.63205 atol=1E-5
    #@test pk.result.Cmax[1] ≈ 0.4 atol=1E-5
    #@test pk.result.MRTlast[1] ≈ 3.10345 atol=1E-5
    #@test pk.result.Tmax[1] ≈ 3.0 atol=1E-5

    pkds = ClinicalTrialUtilities.PK.pkimport(pkdata, [:Subject, :Formulation]; conc = :Concentration, time = :Time)
    pk   = ClinicalTrialUtilities.PK.nca!(pkds)
    @test pk[1, :AUCinf]  ≈ 1.63205 atol=1E-5
    @test pk[1, :Cmax]    ≈ 0.4 atol=1E-5
    @test pk[1, :MRTlast] ≈ 3.10345 atol=1E-5
    @test pk[1, :Tmax]    ≈ 3.0 atol=1E-5

    pk   = ClinicalTrialUtilities.PK.nca!(pkds; calcm = :logt)
    @test pk[1, :AUClast]  ≈ 1.43851 atol=1E-5
    @test pk[1, :AUMClast] ≈ 4.49504 atol=1E-5
    pk   = ClinicalTrialUtilities.PK.nca!(pkds; calcm = :luld)
    @test pk[1, :AUClast]  ≈ 1.43851 atol=1E-5
    @test pk[1, :AUMClast] ≈ 4.49504 atol=1E-5
    #pk = ClinicalTrialUtilities.PK.nca(data; conc = :Concentration, sort=[:Formulation, :Subject], bl = 2.0)

    #pk = ClinicalTrialUtilities.PK.nca(data[1:7,:]; conc = :Concentration, sort=[:Formulation, :Subject])
    #@test pk.result.AUCinf[1] ≈ 1.63205 atol=1E-5
    #@test pk.result.Cmax[1] ≈ 0.4 atol=1E-5
    #@test pk.result.MRTlast[1] ≈ 3.10345 atol=1E-5
    #@test pk.result.Tmax[1] ≈ 3.0 atol=1E-5

    pkds = ClinicalTrialUtilities.PK.pkimport(pkdata[1:7,:], [:Subject, :Formulation]; conc = :Concentration, time = :Time)
    pk   = ClinicalTrialUtilities.PK.nca!(pkds)
    @test pk[1, :AUCinf]  ≈ 1.63205 atol=1E-5
    @test pk[1, :Cmax]    ≈ 0.4 atol=1E-5
    @test pk[1, :MRTlast] ≈ 3.10345 atol=1E-5
    @test pk[1, :Tmax]    ≈ 3.0 atol=1E-5

    #pk = ClinicalTrialUtilities.PK.nca(data[1:7,:]; conc = :Concentration)
    #@test pk.result.AUCinf[1] ≈ 1.63205 atol=1E-5
    #@test pk.result.Cmax[1] ≈ 0.4 atol=1E-5
    #@test pk.result.MRTlast[1] ≈ 3.10345 atol=1E-5
    #@test pk.result.Tmax[1] ≈ 3.0 atol=1E-5

    #pkds = ClinicalTrialUtilities.PK.pkimport(pkdata, [:Subject, :Formulation]; conc = :Concentration, time = :Time)


    #=
    pk = ClinicalTrialUtilities.PK.nca(data; conc = :Concentration, sort=[:Formulation, :Subject], calcm = :logt)
    @test pk.result.AUClast[1] ≈ 1.43851 atol=1E-5
    @test pk.result.AUMClast[1] ≈ 4.49504 atol=1E-5

    pk = ClinicalTrialUtilities.PK.nca(data; conc = :Concentration, sort=[:Formulation, :Subject], calcm = :luld)
    @test pk.result.AUClast[1] ≈ 1.43851 atol=1E-5
    @test pk.result.AUMClast[1] ≈ 4.49504 atol=1E-5
    =#

    pdds = ClinicalTrialUtilities.PK.pdimport(pddata; time=:time, resp=:effect, bl = 3.0)
    pd   = ClinicalTrialUtilities.PK.nca!(pdds)
    @test pd[1,:AUCABL]   ≈ 7.38571428571429 atol=1E-5
    @test pd[1,:AUCBBL]   ≈ 8.73571428571429 atol=1E-5
    ClinicalTrialUtilities.PK.setth!(pdds, 1.5)
    pd   = ClinicalTrialUtilities.PK.nca!(pdds)
    @test pd[1,:AUCATH]   ≈ 13.9595238095238 atol=1E-5
    @test pd[1,:AUCBTH]   ≈ 1.80952380952381 atol=1E-5
    @test pd[1,:TABL]     ≈ 3.48095238095238 atol=1E-5
    @test pd[1,:TBBL]     ≈ 5.51904761904762 atol=1E-5
    @test pd[1,:TATH]     ≈ 5.76190476190476 atol=1E-5
    @test pd[1,:TBTH]     ≈ 3.23809523809524 atol=1E-5
    @test pd[1,:AUCBLNET] ≈ -1.35 atol=1E-5
    @test pd[1,:AUCTHNET] ≈ 12.15 atol=1E-5

    pdds = ClinicalTrialUtilities.PK.pdimport(data, [:Formulation, :Subject]; time=:Time, resp=:Concentration, bl=0.2, th=0.3)
    pd   = ClinicalTrialUtilities.PK.nca!(pdds)
    @test pd[2, :AUCDBLTH] ≈ 0.3416666666666665 atol=1E-5
    ClinicalTrialUtilities.PK.setbl!(pdds, 0.3)
    ClinicalTrialUtilities.PK.setth!(pdds, 0.2)
    pd   = ClinicalTrialUtilities.PK.nca!(pdds)
    @test pd[3, :AUCDBLTH] ≈ 0.3428571428571429 atol=1E-5

    #=
    df = CSV.read(IOBuffer(pddata)) |> DataFrame
    pk = ClinicalTrialUtilities.PK.nca(df; effect=:effect, time=:time, bl = 3.0)
    @test pk.result.AUCABL[1] ≈ 7.38571428571429 atol=1E-5
    @test pk.result.AUCBBL[1] ≈ 8.73571428571429 atol=1E-5

    pk = ClinicalTrialUtilities.PK.nca(df; effect=:effect, time=:time, bl = 3.0, th = 1.5)
    @test pk.result.AUCATH[1] ≈ 13.9595238095238 atol=1E-5
    @test pk.result.AUCBTH[1] ≈ 1.80952380952381 atol=1E-5
    @test pk.result.TABL[1] ≈ 3.48095238095238 atol=1E-5
    @test pk.result.TBBL[1] ≈ 5.51904761904762 atol=1E-5
    @test pk.result.TATH[1] ≈ 5.76190476190476 atol=1E-5
    @test pk.result.TBTH[1] ≈ 3.23809523809524 atol=1E-5
    @test pk.result.AUCBLNET[1] ≈ -1.35 atol=1E-5
    @test pk.result.AUCTHNET[1] ≈ 12.15 atol=1E-5

    pd = ClinicalTrialUtilities.PK.nca(data; effect = :Concentration, sort=[:Formulation, :Subject], bl=0.2, th=0.3).result
    @test pd[2, :AUCDBLTH] ≈ 0.3416666666666665 atol=1E-5

    pd = ClinicalTrialUtilities.PK.nca(data; effect = :Concentration, sort=[:Formulation, :Subject], bl=0.3, th=0.2).result
    @test pd[3, :AUCDBLTH] ≈ 0.3428571428571429 atol=1E-5
    =#
end
