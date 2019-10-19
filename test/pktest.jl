println(" ---------------------------------- ")
@testset "  PK                  " begin


    #pk = ClinicalTrialUtilities.PK.nca(data; conc = :Concentration, sort=[:Formulation, :Subject])
    #@test pk.result.AUCinf[1] ≈ 1.63205 atol=1E-5
    #@test pk.result.Cmax[1] ≈ 0.4 atol=1E-5
    #@test pk.result.MRTlast[1] ≈ 3.10345 atol=1E-5
    #@test pk.result.Tmax[1] ≈ 3.0 atol=1E-5

    pkds = ClinicalTrialUtilities.pkimport(pkdata, [:Subject, :Formulation]; conc = :Concentration, time = :Time)
    pk   = ClinicalTrialUtilities.nca!(pkds)
    @test pk[1, :AUCinf]  ≈ 1.63205 atol=1E-5
    @test pk[1, :Cmax]    ≈ 0.4 atol=1E-5
    @test pk[1, :MRTlast] ≈ 3.10345 atol=1E-5
    @test pk[1, :Tmax]    ≈ 3.0 atol=1E-5

    pk   = ClinicalTrialUtilities.nca!(pkds; calcm = :logt)
    @test pk[1, :AUClast]  ≈ 1.43851 atol=1E-5
    @test pk[1, :AUMClast] ≈ 4.49504 atol=1E-5
    pk   = ClinicalTrialUtilities.nca!(pkds; calcm = :luld)
    @test pk[1, :AUClast]  ≈ 1.43851 atol=1E-5
    @test pk[1, :AUMClast] ≈ 4.49504 atol=1E-5

    pkds = ClinicalTrialUtilities.pkimport(pkdata[1:7,:], [:Subject, :Formulation]; conc = :Concentration, time = :Time)
    pk   = ClinicalTrialUtilities.nca!(pkds)
    @test pk[1, :AUCinf]  ≈ 1.63205 atol=1E-5
    @test pk[1, :Cmax]    ≈ 0.4 atol=1E-5
    @test pk[1, :MRTlast] ≈ 3.10345 atol=1E-5
    @test pk[1, :Tmax]    ≈ 3.0 atol=1E-5

    pkds = ClinicalTrialUtilities.pkimport(pkdata2, [:Subject, :Formulation]; conc = :Concentration, time = :Time)
    pk   = ClinicalTrialUtilities.nca!(pkds)

    pkds = ClinicalTrialUtilities.pkimport(glucose2, [:Subject, :Date]; conc = :glucose, time = :Time)
    pk   = ClinicalTrialUtilities.nca!(pkds)

    pdds = ClinicalTrialUtilities.pdimport(pddata; time=:time, resp=:effect, bl = 3.0)
    pd   = ClinicalTrialUtilities.nca!(pdds)
    @test pd[1,:AUCABL]   ≈ 7.38571428571429 atol=1E-5
    @test pd[1,:AUCBBL]   ≈ 8.73571428571429 atol=1E-5
    ClinicalTrialUtilities.setth!(pdds, 1.5)
    pd   = ClinicalTrialUtilities.nca!(pdds)
    @test pd[1,:AUCATH]   ≈ 13.9595238095238 atol=1E-5
    @test pd[1,:AUCBTH]   ≈ 1.80952380952381 atol=1E-5
    @test pd[1,:TABL]     ≈ 3.48095238095238 atol=1E-5
    @test pd[1,:TBBL]     ≈ 5.51904761904762 atol=1E-5
    @test pd[1,:TATH]     ≈ 5.76190476190476 atol=1E-5
    @test pd[1,:TBTH]     ≈ 3.23809523809524 atol=1E-5
    @test pd[1,:AUCBLNET] ≈ -1.35 atol=1E-5
    @test pd[1,:AUCTHNET] ≈ 12.15 atol=1E-5

    pdds = ClinicalTrialUtilities.pdimport(pkdata, [:Formulation, :Subject]; time=:Time, resp=:Concentration, bl=0.2, th=0.3)
    pd   = ClinicalTrialUtilities.nca!(pdds)
    @test pd[2, :AUCDBLTH] ≈ 0.3416666666666665 atol=1E-5
    ClinicalTrialUtilities.setbl!(pdds, 0.3)
    ClinicalTrialUtilities.setth!(pdds, 0.2)
    pd   = ClinicalTrialUtilities.nca!(pdds)
    @test pd[3, :AUCDBLTH] ≈ 0.3428571428571429 atol=1E-5

end
