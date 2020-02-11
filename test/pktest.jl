println(" ---------------------------------- ")
@testset "  PK                    " begin


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

    pkds = ClinicalTrialUtilities.pkimport(pkdata[1:7,2], pkdata[1:7,1])
    pk   = ClinicalTrialUtilities.nca!(pkds)
    @test pk[:AUClast]    ≈ 1.4499999999999997

    pkds = ClinicalTrialUtilities.pkimport(pkdata2, [:Subject, :Formulation]; conc = :Concentration, time = :Time)
    pk   = ClinicalTrialUtilities.nca!(pkds)
    df   = DataFrame(pk; unst = true)
    sort!(df, :Subject)

    # Linear-trapezoidal rule (linear interpolation)
    #AUC last
    @test round.(df[!, :AUClast], sigdigits = 6) == round.([9585.4218
    10112.176
    5396.5498
    9317.8358
    9561.26
    6966.598
    7029.5735
    7110.6745
    8315.0803
    5620.8945], sigdigits = 6)

    #AUMC last
    @test round.(df[!, :AUMClast], sigdigits = 6) == round.([333582.48
    298701.39
    186032.06
    313955.9
    315181.56
    226977.06
    219797.71
    240526.05
    277613.98
    154893.06], sigdigits = 6)

    #Cmax
    @test df[!, :Cmax] == [190.869
    261.177
    105.345
    208.542
    169.334
    154.648
    153.254
    138.327
    167.347
    125.482]

    #Clast
    @test df[!, :Clast] == [112.846
    85.241
    67.901
    97.625
    110.778
    69.501
    58.051
    74.437
    93.44
    42.191]

    #Adjusted R sq
    @test round.(df[!, :ARsq], digits = 6) == round.([0.71476928
    0.99035145
    0.77630678
    0.83771737
    0.82891994
    0.92517856
    0.96041642
    0.92195356
    0.92130684
    0.86391165], digits = 6)

    #Kel
    @test round.(df[!, :Kel], sigdigits = 6) == round.([0.0033847439
    0.014106315
    0.0032914304
    0.0076953442
    0.0068133279
    0.0076922807
    0.012458956
    0.0089300798
    0.0056458649
    0.017189737], sigdigits = 6)

    @test round.(df[!, :AUCinf], sigdigits = 6) == round.([42925.019
    16154.93
    26026.183
    22004.078
    25820.275
    16001.76
    11688.953
    15446.21
    24865.246
    8075.3242], sigdigits = 6)


    @test df[!, :Tmax] == [1
    1
    1.5
    1
    4
    2.5
    2.5
    4
    3
    2]

    @test round.(df[!, :MRTlast], digits = 6) == round.([34.801023
    29.538786
    34.472406
    33.69408
    32.964438
    32.58076
    31.267574
    33.826053
    33.386807
    27.556657], digits = 6)

    @test round.(df[!, :Clast_pred], sigdigits = 6) == round.([117.30578
    82.53669
    66.931057
    100.76793
    105.29832
    71.939942
    61.172702
    75.604277
    93.761762
    38.810857], sigdigits = 6)


    # Linear-trapezoidal Log Interpolation

    ClinicalTrialUtilities.setdosetime!(pkds, ClinicalTrialUtilities.DoseTime(dose = 120, time = 0, tau = 12))
    pk   = ClinicalTrialUtilities.nca!(pkds, intp = :logt)
    df   = DataFrame(pk; unst = true)
    sort!(df, :Subject)


    @test round.(df[!, :AUCtau], sigdigits = 6) == round.([1670.1018
    2380.2695
    980.34575
    1711.0358
    1738.46
    1409.998
    1436.5595
    1105.0705
    1638.1903
    1293.7065], sigdigits = 6)

    # Linear-trapezoidal Linear Interpolation

    ClinicalTrialUtilities.setdosetime!(pkds, ClinicalTrialUtilities.DoseTime(dose = 120, time = 2, tau = 10))
    pk   = ClinicalTrialUtilities.nca!(pkds)
    df   = DataFrame(pk; unst = true)
    sort!(df, :Subject)

    @test round.(df[!, :AUCtau], sigdigits = 6) == round.([
    1367.7388
    2043.0158
    838.8295
    1444.7985
    1618.6205
    1224.6608
    1265.7013
    880.311
    1417.942
    1167.911], sigdigits = 6)

    # --------------------------------------------------------------------------
    # Linear Up Log Down - Log Interpolation
    # TAU 0 - 36
    ClinicalTrialUtilities.setdosetime!(pkds, ClinicalTrialUtilities.DoseTime(dose = 100, time = 0, tau = 36))
    pk   = ClinicalTrialUtilities.nca!(pkds, calcm = :luldt, intp = :logt)
    df   = DataFrame(pk; unst = true)
    sort!(df, :Subject)

    @test round.(df[!, :Ctau], sigdigits = 6) == round.([135.2155
    132.9739623
    75.47731257
    135.0568587
    133.7790224
    96.70617845
    99.89068375
    103.1329035
    113.3749669
    72.57153458], sigdigits = 6)

    #AUCtau
    #=
    4947.052768
    6265.637822
    2863.392745
    5008.644865
    5415.246398
    3882.203934
    4096.937951 (!)
    3802.30601
    4531.779637
    3744.768624
    =#

    # TAU 0.25 - 9
    ClinicalTrialUtilities.setdosetime!(pkds, ClinicalTrialUtilities.DoseTime(dose = 100, time = 0.25, tau = 9))
    pk   = ClinicalTrialUtilities.nca!(pkds, calcm = :luldt, intp = :logt)
    df   = DataFrame(pk; unst = true)
    sort!(df, :Subject)

    @test round.(df[!, :Ctau], sigdigits = 6) == round.([137.721875
    204.9625055
    84.915625
    141.9969549
    160.0570781
    106.9862605
    123.572375
    84.9595
    142.566125
    122.7865], sigdigits = 6)

    #AUC Tau (!!!)
    #=
    1268.275631
    1831.820528
    754.6493604
    1336.480932
    1310.904516
    1114.240351
    1079.366854
    766.620245
    1219.631864
    970.3062692
    =#

    # TAU 0.0 - 100
    ClinicalTrialUtilities.setdosetime!(pkds, ClinicalTrialUtilities.DoseTime(dose = 100, time = 0.0, tau = 100))
    pk   = ClinicalTrialUtilities.nca!(pkds, calcm = :luldt, intp = :logt)
    df   = DataFrame(pk; unst = true)
    sort!(df, :Subject)

    @test round.(df[!, :Ctau], sigdigits = 6) == round.([106.6989373
    55.60460942
    61.0383917
    81.23535402
    87.01010737
    58.00027625
    43.15724396
    58.8781831
    80.05171762
    23.98401112], sigdigits = 6)

#=
    12646.63632
    11996.6718
    7195.902904
    11794.11692
    12274.83395
    8729.151856
    8395.400098 (!)
    8930.999936 (!)
    10727.4135 (!)
    6389.420453
=#

    # Log-trapezoidal Linear Interpolation

    pkds = ClinicalTrialUtilities.pkimport(pkdata2, [:Subject, :Formulation]; conc = :Concentration, time = :Time)
    pk   = ClinicalTrialUtilities.nca!(pkds, calcm = :logt, intp = :lint)
    df   = DataFrame(pk; unst = true)
    sort!(df, :Subject)

    @test round.(df[!, :AUClast], sigdigits = 6) == round.([9572.8582
    10054.037
    5391.5322
    9296.2179
    9518.6531
    6948.5757
    6987.0645
    7064.7816
    8298.9634
    5485.6538], sigdigits = 6)

    @test round.(df[!, :AUCinf], sigdigits = 6) == round.([42912.456
    16096.791
    26021.165
    21982.46
    25777.668
    15983.737
    11646.444
    15400.317
    24849.129
    7940.0834], sigdigits = 6)

    pk   = ClinicalTrialUtilities.nca!(pkds, calcm = :luldt, io = io, verbose = true)
    df   = DataFrame(pk; unst = true)
    sort!(df, :Subject)

    @test round.(df[!, :AUClast], sigdigits = 6) == round.([9573.810558691312
    10054.286478059563
    5392.457219413793
    9297.096334450325
    9519.181808199797
    6948.985621117448
    6988.960344867885
    7073.306755718137
    8303.373085532965
    5486.838889441992], sigdigits = 6)
    #---------------------------------------------------------------------------
    # Elimination range check
    ClinicalTrialUtilities.setkelauto!(pkds, false)
    ClinicalTrialUtilities.setelimrange!(pkds, ClinicalTrialUtilities.ElimRange(11, 15))
    pk   = ClinicalTrialUtilities.nca!(pkds)
    df   = DataFrame(pk; unst = true)
    sort!(df, :Subject)

    @test round.(df[!, :AUCinf], sigdigits = 6) == round.([655061.1994219155
    15395.482644599213
    22227.779049490928
    53919.75612131339
    21428.55291240885
    21778.845190342177
    15100.132319764292
    25511.20156654014
    27923.624963124363
    7386.834076155846], sigdigits = 6)

    # glucose2

    pkds = ClinicalTrialUtilities.pkimport(glucose2, [:Subject, :Date]; conc = :glucose, time = :Time)
    pk   = ClinicalTrialUtilities.nca!(pkds)
    df   = DataFrame(pk; unst = true)
    p    = ClinicalTrialUtilities.plot(pkds; pagesort = [:Date], typesort = [:Subject])
    @test length(p) == 2
    print(io, pk[1].subject.keldata)
end


println(" ---------------------------------- ")
@testset "  PD                    " begin

    #PD
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

    print(io, pdds)
    print(io, pdds[1])
    print(io, pd)
    print(io, pd[1])

end
