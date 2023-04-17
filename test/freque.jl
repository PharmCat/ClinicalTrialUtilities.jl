#=
println(" ---------------------------------- ")
@testset "#8  Frequency           " begin
    df   = freqdat
    ctab =  ClinicalTrialUtilities.contab(df, row = :row, col = :col)
    @test ctab.tab == [9 8; 5 21]

    frtab =  DataFrame(ClinicalTrialUtilities.freque(df; vars=:row, alpha = 0.05))
    @test frtab[1,2] == 17

    @test ClinicalTrialUtilities.fisher(ClinicalTrialUtilities.ConTab([12 23; 22 33]))[2] ≈ 0.3752579856319276

    @test ClinicalTrialUtilities.pirson(ClinicalTrialUtilities.ConTab([12 23; 22 33]))[4] ≈ 0.5856943077831229


    ctds = ClinicalTrialUtilities.contab(metadf, [:trial]; row = :group, col = :result)
    cmht = ClinicalTrialUtilities.metaprop(ctds, type = :or, model = :fixed, zeroadj = 0.0, tau = :ho)
    @test cmht.est   ≈ 0.4164682333774169 atol=1E-5
    @test cmht.var   ≈ 0.13107613010773247 atol=1E-5
    @test cmht.tausq ≈ 0.5554302118884364 atol=1E-5
    println(io, cmht)

    cmht = ClinicalTrialUtilities.metaprop(ctds, type = :or, model = :random, zeroadj = 0.0, tau = :dl)
    @test cmht.est   ≈ 0.46774668790076773 atol=1E-5
    @test cmht.var   ≈ 0.2500964678287197 atol=1E-5
    @test cmht.tausq ≈ 0.3284301973957149 atol=1E-5
    println(io, cmht)

    #Not validated
    cmht = ClinicalTrialUtilities.metaprop(ctds, type = :diff, model = :random, zeroadj = 0.0, tau = :ho)
    cmht = ClinicalTrialUtilities.metaprop(ctds, type = :diff, model = :random, zeroadj = 0.0, tau = :dl)
    cmht = ClinicalTrialUtilities.metaprop(ctds, type = :diff, model = :random, zeroadj = 0.0, tau = :hm)

    cmht = ClinicalTrialUtilities.metaprop(ctds, type = :rr, model = :random, zeroadj = 0.0, tau = :ho)
    cmht = ClinicalTrialUtilities.metaprop(ctds, type = :rr, model = :random, zeroadj = 0.0, tau = :dl)
    cmht = ClinicalTrialUtilities.metaprop(ctds, type = :rr, model = :random, zeroadj = 0.0, tau = :hm)

    # MCNM Test
    ctds = ClinicalTrialUtilities.mcnmcontab(metadf, [:trial]; row = :group, col = :result)

    cmht = ClinicalTrialUtilities.metaprop(ctds, type = :or, model = :fixed, zeroadj = 0.0, tau = :ho)
    println(io, cmht)

    cmht = ClinicalTrialUtilities.metaprop(ctds, type = :or, model = :random, zeroadj = 0.0, tau = :ho)
    cmht = ClinicalTrialUtilities.metaprop(ctds, type = :or, model = :random, zeroadj = 0.0, tau = :dl)
    cmht = ClinicalTrialUtilities.metaprop(ctds, type = :diff, model = :random, zeroadj = 0.0, tau = :hm)

    println(io, cmht)




end
=#