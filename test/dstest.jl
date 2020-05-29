println(" ---------------------------------- ")
@testset "#7  Descriptives        " begin
    df = descriptivedat
    ds = ClinicalTrialUtilities.descriptive(df, stats = :all, sort = [:C1, :C2], vars=[:P1, :P2])
    @test ds[1,:n]        ≈ 20 atol=1E-5
    @test ds[1,:min]      ≈ 11.30162773 atol=1E-5
    @test ds[1,:max]      ≈ 13561.95328 atol=1E-5
    @test ds[1,:range]    ≈ 13550.65165227 atol=1E-5
    @test ds[1,:mean]     ≈ 1620.9521026095 atol=1E-5
    @test ds[1,:var]      ≈ 14034427.3245356 atol=1E-5
    @test ds[1,:sd]       ≈ 3746.25510670797 atol=1E-5
    @test ds[1,:sem]      ≈ 837.688107965475 atol=1E-5
    @test ds[1,:cv]       ≈ 231.114485164431 atol=1E-5
    @test ds[1,:harmmean] ≈ 43.7567236075722 atol=1E-5
    @test ds[1,:geomean]  ≈ 205.844666821181 atol=1E-5
    @test ds[1,:geovar]   ≈ 4.91093447695305 atol=1E-5
    @test ds[1,:geosd]    ≈ 2.21606283235676 atol=1E-5
    @test ds[1,:geocv]    ≈ 1160.88856294209 atol=1E-5
    @test ds[1,:skew]     ≈ 2.63411684077405 atol=1E-5
    @test ds[1,:kurt]     ≈ 5.20836944102797 atol=1E-5
    @test ds[1,:uq]       ≈ 1130.162773 atol=1E-5
    @test ds[1,:median]   ≈ 135.6195328 atol=1E-5
    @test ds[1,:lq]       ≈ 31.022968125 atol=1E-5   #SPSS 30.853444  Phoenix WNL 30.683919295
    @test ds[1,:iqr]      ≈ 1099.139804875 atol=1E-5 #Same as above
    @test ds[1,:mode]     ≈ 645.8072989 atol=1E-5

    ds = ClinicalTrialUtilities.descriptive(df, stats = [:geosd], sort = [:C1, :C2], vars=[:P1, :P2])
    @test ds[1,:geosd]    ≈ 2.21606283235676 atol=1E-5

    ds = ClinicalTrialUtilities.descriptive(df, stats = [:mean, :geomean, :geocv], vars=[:P1, :P2, :C3])
    @test ds[1,:mean]     ≈ 2181.5170114224 atol=1E-5
    @test ds[1,:geomean]  ≈ 4.86763332369581 atol=1E-5
    @test ds[1,:geocv]    ≈ 124574201599.317 atol=1E-2

    ds = ClinicalTrialUtilities.descriptive(df, stats = [:n, :ses, :sek],  vars=[:P1, :P2, :C3])
    @test ds[3,:ses]  ≈ 0.171925 atol=1E-6
    @test ds[3,:sek]  ≈ 0.342202 atol=1E-6

    ds = ClinicalTrialUtilities.descriptive(df[:,:C3], stats = [:n, :ses, :sek])
    @test ds[:ses]  ≈ 0.171925 atol=1E-6
    @test ds[:sek]  ≈ 0.342202 atol=1E-6

    ds = ClinicalTrialUtilities.descriptive(df, sort=[:C3], vars=[:P1, :P2])
    @test ds[3,:mean]     ≈ 51.35 atol=1E-5
    @test ds[3,:sem]      ≈ 48.65 atol=1E-5
    @test ds[3,:median]   ≈ 51.35 atol=1E-3

    df = negdat
    ds = ClinicalTrialUtilities.descriptive(df, stats = :all, vars=[:P1])
    @test ds[1,:harmmean] === NaN
    @test ds[1,:geomean]  === NaN
    @test ds[1,:geovar]   === NaN
    @test ds[1,:geosd]    === NaN
    @test ds[1,:geocv]    === NaN
end
