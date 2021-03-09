println(" ---------------------------------- ")
@testset "#10  Pharmacokinetics   " begin

    # Basic tests
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

    #setkelauto! getkelauto
    ClinicalTrialUtilities.setkelauto!(pkds, false, [1,2])
    ClinicalTrialUtilities.setkelauto!(pkds, true, 1)
    @test ClinicalTrialUtilities.getkelauto(pkds[1]) == true

    ClinicalTrialUtilities.applyelimrange!(pk, ClinicalTrialUtilities.ElimRange(4,7))
    ClinicalTrialUtilities.applyelimrange!(pk, ClinicalTrialUtilities.ElimRange(5,7), [1,2,3])
    ClinicalTrialUtilities.applyelimrange!(pk, ClinicalTrialUtilities.ElimRange(4,7), 1)
    ClinicalTrialUtilities.applyelimrange!(pk, ClinicalTrialUtilities.ElimRange(5,7), Dict(:Subject => 2))
    @test pk[1].subject.kelrange.kelstart == 4
    @test pk[2].subject.kelrange.kelstart == 5

    pkds = ClinicalTrialUtilities.pkimport(pkdata[1:7,:], [:Subject, :Formulation]; conc = :Concentration, time = :Time)
    pk   = ClinicalTrialUtilities.nca!(pkds)
    @test pk[1, :AUCinf]  ≈ 1.63205 atol=1E-5
    @test pk[1, :Cmax]    ≈ 0.4 atol=1E-5
    @test pk[1, :MRTlast] ≈ 3.10345 atol=1E-5
    @test pk[1, :Tmax]    ≈ 3.0 atol=1E-5

    pkds = ClinicalTrialUtilities.pkimport(pkdata[1:7,2], pkdata[1:7,1])
    pk   = ClinicalTrialUtilities.nca!(pkds)
    @test pk[:AUClast]    ≈ 1.4499999999999997

#Program:
    #AUClast
    #AUCinf
    #AUCpct
    #Tmax
    #Cmax
    #Clast
    #Adjusted R sq
    #Kel
    #HL
    #AUMClast
    #MRTlast
    #Clast_pred

    # --------------------------------------------------------------------------
    # Linear Trapezoidal Linear Interpolation
    # calcm = :lint, intp = :lint

    pkds = ClinicalTrialUtilities.pkimport(pkdata2, [:Subject, :Formulation]; conc = :Concentration, time = :Time)
    pk   = ClinicalTrialUtilities.nca!(pkds)
    df   = DataFrame(pk; unst = true)
    sort!(df, :Subject)

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
    #AUCinf
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
    #AUCpct
    @test round.(df[!, :AUCpct], sigdigits = 5) == round.([77.669383
    37.405019
    79.26492
    57.65405
    62.969953
    56.463551
    39.861391
    53.964925
    66.559429
    30.394194], sigdigits = 5)
    #Tmax
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
    #HL
    @test round.(df[!, :HL], sigdigits = 5) == round.([204.78571
    49.137367
    210.59148
    90.073577
    101.73401
    90.10945
    55.634451
    77.619371
    122.77077
    40.323315], sigdigits = 5)
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
    #MRTlast
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
    #Clast_pred
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


    # --------------------------------------------------------------------------
    # Linear Log Trapezoidal
    # Log-trapezoidat rule after Tmax if c₁ > 0 and c₂ > 0, else Linear trapezoidal used;


    pk   = ClinicalTrialUtilities.nca!(pkds, calcm = :logt,  intp = :logt, io = io, verbose = true)
    df   = DataFrame(pk; unst = true)
    sort!(df, :Subject)
    #AUClast
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
    #AUCinf
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
    #AUCpct
    @test round.(df[!, :AUCpct], sigdigits = 5) == round.([77.692122
    37.540119
    79.280204
    57.710748
    63.074033
    56.527216
    40.006884
    54.12574
    66.602599
    30.911888], sigdigits = 5)
    #Tmax
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
    @test round.(df[!, :Kel], sigdigits = 6) == round.([0.003384744
    0.014106315
    0.00329143
    0.007695344
    0.006813328
    0.007692281
    0.012458956
    0.00893008
    0.005645865
    0.017189737], sigdigits = 6)
    #HL
    @test round.(df[!, :HL], sigdigits = 5) == round.([204.78571
    49.137367
    210.59148
    90.073577
    101.73401
    90.10945
    55.634451
    77.619371
    122.77077
    40.323315], sigdigits = 5)

    #Clast_pred
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


    # --------------------------------------------------------------------------
    # Linear Up Log Down
    # Linear Up Log Down everywhere if c₁ > c₂ > 0, else Linear trapezoidal used;

    pk   = ClinicalTrialUtilities.nca!(pkds, calcm = :luld,  intp = :logt, io = io, verbose = true)
    df   = DataFrame(pk; unst = true)
    sort!(df, :Subject)

    #AUClast
    @test round.(df[!, :AUClast], sigdigits = 6) == round.([9573.8106
    10054.286
    5392.4572
    9297.0963
    9519.1809
    6948.9856
    6988.7726
    7073.0922
    8303.3586
    5486.8389], sigdigits = 6)
    #AUCinf
    @test round.(df[!, :AUCinf], sigdigits = 6) == round.([42913.408
    16097.041
    26022.09
    21983.338
    25778.196
    15984.147
    11648.152
    15408.628
    24853.524
    7941.2686], sigdigits = 6)
    #AUCpct
    #Tmax
    @test round.(df[!, :Tmax], sigdigits = 6) == round.([1
    1
    1.5
    1
    4
    2.5
    2.5
    4
    3
    2], sigdigits = 6)
    #Cmax
    @test round.(df[!, :Cmax], sigdigits = 6) == round.([190.869
    261.177
    105.345
    208.542
    169.334
    154.648
    153.254
    138.327
    167.347
    125.482], sigdigits = 6)
    #Clast
    @test round.(df[!, :Clast], sigdigits = 6) == round.([112.846
    85.241
    67.901
    97.625
    110.778
    69.501
    58.051
    74.437
    93.44
    42.191], sigdigits = 6)
    #Adjusted R sq
    @test round.(df[!, :ARsq], sigdigits = 6) == round.([0.71476928
    0.99035145
    0.77630678
    0.83771737
    0.82891994
    0.92517856
    0.96041642
    0.92195356
    0.92130684
    0.86391165], sigdigits = 6)
    #Kel
    @test round.(df[!, :Kel], sigdigits = 6) == round.([0.003384744
    0.014106315
    0.00329143
    0.007695344
    0.006813328
    0.007692281
    0.012458956
    0.00893008
    0.005645865
    0.017189737], sigdigits = 6)
    #HL
    @test round.(df[!, :HL], sigdigits = 5) == round.([204.78571
    49.137367
    210.59148
    90.073577
    101.73401
    90.10945
    55.634451
    77.619371
    122.77077
    40.323315], sigdigits = 5)
    #AUMClast
    #MRTlast
    #Clast_pred
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

    # --------------------------------------------------------------------------
    # Linear Trapezoidal Linear/Log Interpolation.
    # calcm = :lint, intp = :logt
    # Linear trapezoidal everywhere + log interpolation

    pkds = ClinicalTrialUtilities.pkimport(pkdata2, [:Subject, :Formulation]; conc = :Concentration, time = :Time)
    pk   = ClinicalTrialUtilities.nca!(pkds; intp = :logt)
    df   = DataFrame(pk; unst = true)
    sort!(df, :Subject)

    #AUClast
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


    # --------------------------------------------------------------------------
    # Not present in Phoenix
    # Linear Up Log Down after Tmax if c₁ > c₂ > 0, else Linear trapezoidal used;

    pk   = ClinicalTrialUtilities.nca!(pkds, calcm = :luldt,  intp = :logt, io = io, verbose = true)
    df   = DataFrame(pk; unst = true)
    sort!(df, :Subject)

    #AUClast
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
    #AUCinf
    @test round.(df[!, :AUCinf], sigdigits = 6) == round.([42913.407881300096
    16097.04111262767
    26022.090028134975
    21983.33845321829
    25778.19670343543
    15984.14736468632
    11648.33951823218
    15408.84256127436
    24853.53879891218
    7941.268553852982], sigdigits = 6)


    # --------------------------------------------------------------------------
    # NOT USED
    # Log-trapezoidal Linear Interpolation
    # if  calcm == :logt / :luld / :luldt calculation method used - log interpolation used when possible
    #=
    ClinicalTrialUtilities.setdosetime!(pkds, ClinicalTrialUtilities.DoseTime(dose = 100, time = 0, tau = 0))
    pk   = ClinicalTrialUtilities.nca!(pkds, calcm = :logt, intp = :lint)
    df   = DataFrame(pk; unst = true)
    sort!(df, :Subject)
    #AUClast
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
    #AUCinf
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
    #AUCpct
    @test round.(df[!, :AUCpct], sigdigits = 6) == round.([77.692122
    37.540119
    79.280204
    57.710748
    63.074033
    56.527216
    40.006884
    54.12574
    66.602599
    30.911888], sigdigits = 6)
    #Tmax
    @test round.(df[!, :Tmax], sigdigits = 6) == round.([1
    1
    1.5
    1
    4
    2.5
    2.5
    4
    3
    2], sigdigits = 6)
    #Cmax
    @test round.(df[!, :Cmax], sigdigits = 6) == round.([190.869
    261.177
    105.345
    208.542
    169.334
    154.648
    153.254
    138.327
    167.347
    125.482], sigdigits = 6)
    #Clast
    @test round.(df[!, :Clast], sigdigits = 6) == round.([112.846
    85.241
    67.901
    97.625
    110.778
    69.501
    58.051
    74.437
    93.44
    42.191], sigdigits = 6)
    #Adjusted R sq :ARsq
    @test round.(df[!, :ARsq], sigdigits = 6) == round.([0.71476928
    0.99035145
    0.77630678
    0.83771737
    0.82891994
    0.92517856
    0.96041642
    0.92195356
    0.92130684
    0.86391165], sigdigits = 6)
    #Kel
    @test round.(df[!, :Kel], sigdigits = 6) == round.([0.003384744
    0.014106315
    0.00329143
    0.007695344
    0.006813328
    0.007692281
    0.012458956
    0.00893008
    0.005645865
    0.017189737], sigdigits = 6)
    #HL
    @test round.(df[!, :HL], sigdigits = 5) == round.([204.78571
    49.137367
    210.59148
    90.073577
    101.73401
    90.10945
    55.634451
    77.619371
    122.77077
    40.323315], sigdigits = 5)
    =#


    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # Steady state or Tau defined
    # --------------------------------------------------------------------------
    # Linear-trapezoidal, Log Interpolation

    ClinicalTrialUtilities.setdosetime!(pkds, ClinicalTrialUtilities.DoseTime(dose = 120, time = 0, tau = 12))
    pk   = ClinicalTrialUtilities.nca!(pkds, intp = :logt)
    df   = DataFrame(pk; unst = true)
    sort!(df, :Subject)

    #AUCtau
    #Ctau
    #Cavg
    #Swingtau
    #Fluctau
    #Accind

    #AUCtau
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

    # --------------------------------------------------------------------------
    # Linear-trapezoidal Linear Interpolation
    # TAU 2 - 10

    ClinicalTrialUtilities.setdosetime!(pkds, ClinicalTrialUtilities.DoseTime(dose = 120, time = 2, tau = 10))
    pk   = ClinicalTrialUtilities.nca!(pkds)
    df   = DataFrame(pk; unst = true)
    sort!(df, :Subject)

    #AUCtau
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
    # Linear Up Log Down, Log Interpolation
    # TAU 0 - 36

    ClinicalTrialUtilities.setdosetime!(pkds, ClinicalTrialUtilities.DoseTime(dose = 100, time = 0, tau = 36))
    pk   = ClinicalTrialUtilities.nca!(pkds, calcm = :luld, intp = :logt)
    df   = DataFrame(pk; unst = true)
    sort!(df, :Subject)

    #AUClast
    @test round.(df[!, :AUClast], sigdigits = 6) == round.([9573.810559
    10054.28648
    5392.457219
    9297.096334
    9519.180874
    6948.985621
    6988.772632
    7073.092214
    8303.358586
    5486.838889], sigdigits = 6)
    #AUCtau
    @test round.(df[!, :AUCtau], sigdigits = 6) == round.([4947.052768
    6265.637822
    2863.392745
    5008.644865
    5415.246398
    3882.203934
    4096.937951
    3802.30601
    4531.779637
    3744.768624], sigdigits = 6)
    #Ctau
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
    #Cavg
    @test round.(df[!, :Cavg], sigdigits = 6) == round.([137.4181324
    174.0454951
    79.53868737
    139.129024
    150.423511
    107.8389982
    113.803832
    105.6196114
    125.8827677
    104.0213507], sigdigits = 6)
    #Swingtau
    @test round.(df[!, :Swingtau], sigdigits = 6) == round.([0.411591127
    0.964121363
    0.39571742
    0.544105216
    0.265773938
    0.599153255
    0.534217149
    0.341249934
    0.476048942
    0.729080151], sigdigits = 6)
    #Fluctau
    @test round.(df[!, :Fluctau] * 100, sigdigits = 6) == round.([40.49938608
    73.660647
    37.5511445
    52.81798086
    23.6365827
    53.72993308
    46.8906146
    33.32155462
    42.87483828
    50.8650052], sigdigits = 6)
    #Accind
    @test round.(df[!, :Accind], sigdigits = 6) == round.([8.716910716
    2.511311393
    8.949296376
    4.132742748
    4.597396017
    4.1341712
    2.766795121
    3.637329819
    5.43694762
    2.167194353], sigdigits = 6)

    # Steady state, partial areas
    # --------------------------------------------------------------------------
    # Linear Up Log Downl, Log Interpolation
    # TAU 0.25 - 9

    ClinicalTrialUtilities.setdosetime!(pkds, ClinicalTrialUtilities.DoseTime(dose = 100, time = 0.25, tau = 9))
    pk   = ClinicalTrialUtilities.nca!(pkds, calcm = :luld, intp = :logt, io = io, verbose = true)
    df   = DataFrame(pk; unst = true)
    sort!(df, :Subject)

    #AUClast
    @test round.(df[!, :AUClast], sigdigits = 6) == round.([9566.596809
    10054.28648
    5392.457219
    9297.096334
    9519.180874
    6948.985621
    6988.772632
    7058.818964
    8302.368086
    5486.838889], sigdigits = 6)
    #AUCtau
    @test round.(df[!, :AUCtau], sigdigits = 6) == round.([1268.275631
    1831.820528
    754.6493604
    1336.480932
    1310.904516
    1114.240351
    1079.366854
    766.620245
    1219.631864
    970.3062692], sigdigits = 6)
    #Ctau
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
    #Cavg
    @test round.(df[!, :Cavg], sigdigits = 6) == round.([140.9195145
    203.5356142
    83.84992893
    148.4978814
    145.6560573
    123.8044835
    119.9296505
    85.18002722
    135.5146515
    107.8118077], sigdigits = 6)
    #Swing tau
    @test round.(df[!, :Swingtau], sigdigits = 6) == round.([0.38590184
    0.27426721
    0.240584404
    0.468637128
    0.057960085
    0.445494022
    0.240196282
    0.628152237
    0.173820219
    0.021952739], sigdigits = 6)
    #Fluctuation%_Tau
    @test round.(df[!, :Fluctau] * 100, sigdigits = 6) == round.([37.71452462
    27.61899665
    24.36421266
    44.8121175
    6.369060133
    38.4975876
    24.74919661
    62.65259796
    18.28649133
    2.500189968], sigdigits = 6)
    #Accind
    @test round.(df[!, :Accind], sigdigits = 6) == round.([33.3295745
    8.38726982
    34.26016611
    14.94451581
    16.81301567
    14.95026394
    9.427514161
    12.94903929
    20.18432092
    6.976692409], sigdigits = 6)

    #=
    #Kel
    @test round.(df[!, :Kel], sigdigits = 6) == round.([0.003384744
    0.014106315
    0.00329143
    0.007695344
    0.006813328
    0.007692281
    0.012458956
    0.00893008
    0.005645865
    0.017189737], sigdigits = 6)
    =#

    # TAU 0.0 - 100

    ClinicalTrialUtilities.setdosetime!(pkds, ClinicalTrialUtilities.DoseTime(dose = 100, time = 0.0, tau = 100))
    pk   = ClinicalTrialUtilities.nca!(pkds, calcm = :luld, intp = :logt)
    df   = DataFrame(pk; unst = true)
    sort!(df, :Subject)

    #AUClast
    @test round.(df[!, :AUClast], sigdigits = 6) == round.([9573.810559
    10054.28648
    5392.457219
    9297.096334
    9519.180874
    6948.985621
    6988.772632
    7073.092214
    8303.358586
    5486.838889], sigdigits = 6)
    #AUCtau
    @test round.(df[!, :AUCtau], sigdigits = 6) == round.([12646.63632
    11996.6718
    7195.902904
    11794.11692
    12274.83395
    8729.151856
    8395.400098
    8930.999936
    10727.4135
    6389.420453], sigdigits = 6)
    #Ctau
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
    @test round.(df[!, :Cavg], sigdigits = 6) == round.([126.4663632
    119.966718
    71.95902904
    117.9411692
    122.7483395
    87.29151856
    83.95400098
    89.30999936
    107.274135
    63.89420453], sigdigits = 6)
    #Swing tau
    @test round.(df[!, :Swingtau], sigdigits = 6) == round.([0.78885568
    3.697038658
    0.725880992
    1.567133516
    0.94614172
    1.666332128
    2.551060864
    1.349376165
    1.090486063
    4.231902178], sigdigits = 6)


    # If TAU = Tlast then AUCtau = AUClast
    ClinicalTrialUtilities.setdosetime!(pkds, ClinicalTrialUtilities.DoseTime(dose = 100, time = 0, tau = 72))
    pk   = ClinicalTrialUtilities.nca!(pkds, calcm = :luld, intp = :logt)
    df   = DataFrame(pk; unst = true)
    sort!(df, :Subject)
    @test df[!, :AUCtau] == df[!, :AUClast]



    # Utilities
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
    p    = ClinicalTrialUtilities.pkplot(pkds; pagesort = [:Date], typesort = [:Subject])
    @test length(p) == 2
    print(io, pk[1].subject.keldata)
    p    = ClinicalTrialUtilities.pkplot(pkds; pagesort = nothing, typesort = [:Subject])
    @test length(p) == 1
    p    = ClinicalTrialUtilities.pkplot(pkds; pagesort = nothing, typesort = nothing, legend = false, xlabel = "L1", ylabel = "L2", xlims = (0,12))
    @test length(p) == 1
    p    = ClinicalTrialUtilities.pkplot(pkds; pagesort = [:Date], typesort = nothing)
    @test length(p) == 2
    p    = ClinicalTrialUtilities.pkplot(pkds[1])

    # The Theoph dataset: 132 observations from 12 subjects
    # see also  https://github.com/asancpt/NonCompart-tests/blob/master/docs/validation_0.4.4.pdf
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6989226/

    pkds = ClinicalTrialUtilities.pkimport(theodata, [:Subject]; conc = :conc, time = :Time)
    pk   = ClinicalTrialUtilities.nca!(pkds; calcm = :luld)
    df   = DataFrame(pk; unst = true)
    @test round.(df[!, :AUClast], sigdigits = 6) == round.([ 147.2347485
  88.7312755
  95.8781978
 102.6336232
 118.1793538
  71.6970150
  87.9692274
  86.8065635
  83.9374360
 135.5760701
  77.8934723
 115.2202082], sigdigits = 6)
    #=
    @test round.(df[!, :AUCinf], sigdigits = 6) == round.([214.9236316
  97.3779346
 106.1276685
 114.2162046
 136.3047316
  81.74333453119011 #82.1758833 #6 Rsq&
 100.9876292
 102.1533003
  97.52000174742838 #97.5200039
 167.8600307
  86.9026173
 125.8315397], sigdigits = 6)
 =#
    # NaN PK LimitRule test
    nanpkdata.Concentration = ClinicalTrialUtilities.tryfloatparse!(nanpkdata.Concentration)
    pkds = ClinicalTrialUtilities.pkimport(nanpkdata, [:Subject, :Formulation]; conc = :Concentration, time = :Time)
    rule = ClinicalTrialUtilities.LimitRule(0, 0, NaN, 0, true)
    ClinicalTrialUtilities.applyncarule!(pkds, rule)
    @test length(pkds[1]) == 5
    @test pkds[3].obs[1]  == 0.0
    @test pkds[2].obs[5]  == 0.1
    pkds = ClinicalTrialUtilities.pkimport(nanpkdata, [:Subject, :Formulation]; conc = :Concentration, time = :Time)
    rule = ClinicalTrialUtilities.LimitRule(0, 0, NaN, -1, false)
    ClinicalTrialUtilities.applyncarule!(pkds, rule)
    @test pkds[3].obs[6]  === NaN

    ClinicalTrialUtilities.setdosetime!(pkds, ClinicalTrialUtilities.DoseTime(100,0), Dict(:Subject => 1, :Formulation => "R"))
    @test pkds[3].dosetime.dose == 100

    pks = ClinicalTrialUtilities.findfirst(Dict(:Subject => 1, :Formulation => "R"), pkds)
end

println(" ---------------------------------- ")
@testset "#11 Urine PK            " begin
    upk = ClinicalTrialUtilities.upkimport(upkdata, [:subj]; stime = :st, etime = :et, conc = :conc, vol = :vol)
    unca = ClinicalTrialUtilities.nca!(upk[1])
    unca = ClinicalTrialUtilities.nca!(upk)
    @test unca[1][:maxrate] ≈ 4.0
    @test unca[1][:mTmax]   ≈ 1.5
    @test unca[1][:ar]      ≈ 16.0
    @test unca[1][:volume]  ≈ 11.0
     p    = ClinicalTrialUtilities.pkplot(upk, ylabel = "Excretion")
end
println(" ---------------------------------- ")
@testset "#12  Pharmacodynamics   " begin


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
println(" ---------------------------------- ")
@testset "#13  Sparse PK          " begin
    pks = ClinicalTrialUtilities.PKSubject(sparse_pk; conc = :Concentration, time = :Time)
    auc = ClinicalTrialUtilities.auc_sparse(pks)
    @test auc ≈ 1.35 atol=1E-5
end
