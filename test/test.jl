# Clinical Trial Utilities
# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)
using Distributions, Random, DataFrames, CSV, Test, Plots, StatsBase, StableRNGs

path    = dirname(@__FILE__)
io      = IOBuffer();

include("testdata.jl")

@testset "  Info:                 " begin
    ClinicalTrialUtilities.info()
    ClinicalTrialUtilities.citation(io = io)
    ClinicalTrialUtilities.licence(io = io)
end

println(" ---------------------------------- ")
println(" ---------   START TEST   --------- ")
println(" ---------------------------------- ")
println(" ---------------------------------- ")
println("")
include("internal.jl")

include("ctsamplen.jl")

println(" ---------------------------------- ")
@testset "#3  besamplen Test      " begin
    t = ClinicalTrialUtilities.besamplen(;theta0=0.95, theta1=0.8, theta2=1.25, cv=0.2, alpha=0.05, beta=0.2, logscale=true, method=:owenq)
    @test ClinicalTrialUtilities.besamplen(t.task).result == 20
    @test ClinicalTrialUtilities.ctsamplen(t.task).result == 20
    Base.show(io, t)
    Base.show(io, t.task)
    #@test n == 20 && round(p, digits=7) == 0.8346802
    t = ClinicalTrialUtilities.besamplen(;theta0=0.95, theta1=0.8, theta2=1.25, cv=0.3, alpha=0.05, beta=0.2, logscale=true, method=:owenq)
    @test ClinicalTrialUtilities.besamplen(t.task).result == 40
    Base.show(io, t)
    Base.show(io, t.task)
    #@test n == 40 && round(p, digits=7) == 0.8158453
    t = ClinicalTrialUtilities.besamplen(;theta0=1.0, theta1=0.8, theta2=1.25, cv=0.3, alpha=0.05, beta=0.1, logscale=true, method=:owenq)
    @test ClinicalTrialUtilities.besamplen(t.task).result  == 40
    Base.show(io, t)
    Base.show(io, t.task)
    #@test n == 40 && round(p, digits=7) == 0.9095603
    t = ClinicalTrialUtilities.besamplen(;theta0=1.05, theta1=0.8, theta2=1.25, cv=0.4, alpha=0.05, beta=0.15, logscale=true, method=:nct)
    @test ClinicalTrialUtilities.besamplen(t.task).result == 74
    Base.show(io, t)
    Base.show(io, t.task)
    #@test n == 74 && round(p, digits=7) == 0.8558178
    t = ClinicalTrialUtilities.besamplen(;theta0=1.05, theta1=0.9, theta2=1.25, cv=0.4, alpha=0.05, beta=0.15, logscale=true, method=:nct)
    @test ClinicalTrialUtilities.besamplen(t.task).result == 108
    Base.show(io, t)
    Base.show(io, t.task)
    #@test n == 108 && round(p, digits=7) == 0.8506248
    t = ClinicalTrialUtilities.besamplen(;theta0=1.05, theta1=0.8, theta2=1.2, cv=0.5, alpha=0.05, beta=0.2, logscale=true, method=:nct)
    @test ClinicalTrialUtilities.besamplen(t.task).result == 158
    Base.show(io, t)
    Base.show(io, t.task)
    #@test n == 158 && round(p, digits=7) == 0.8039191
    t = ClinicalTrialUtilities.besamplen(;theta0=1.05, theta1=0.8, theta2=1.25, cv=0.8, alpha=0.05, beta=0.2, logscale=true, method=:shifted)
    @test ClinicalTrialUtilities.besamplen(t.task).result == 210
    Base.show(io, t)
    Base.show(io, t.task)
    #@test n == 210 && round(p, digits=7) == 0.8012471
    t = ClinicalTrialUtilities.besamplen(;theta0=0.0, theta1=-0.2, theta2=0.2, sd=0.5, alpha=0.05, beta=0.2, logscale=false, method=:owenq)
    @test ClinicalTrialUtilities.besamplen(t.task).result == 110
    Base.show(io, t)
    Base.show(io, t.task)
    #@test n == 110 && round(p, digits=7) == 0.8074124
    t = ClinicalTrialUtilities.besamplen(;theta0=0.0, theta1=-0.2, theta2=0.2, sd=2.0, alpha=0.05, beta=0.2, logscale=false, method=:owenq)
    @test ClinicalTrialUtilities.besamplen(t.task).result == 1716
    Base.show(io, t)
    Base.show(io, t.task)
    #@test n == 1716 && round(p, digits=7) == 0.8005618
    t = ClinicalTrialUtilities.besamplen(;theta0=0.0, theta1=-0.2, theta2=0.2, sd=2.0, alpha=0.001, beta=0.2, logscale=false, method=:owenq)
    @test ClinicalTrialUtilities.besamplen(t.task).result == 3828
    Base.show(io, t)
    Base.show(io, t.task)
    #@test n == 3828 && round(p, digits=7) == 0.8001454
    t = ClinicalTrialUtilities.besamplen(;theta0=0.0, theta1=-0.2, theta2=0.2, sd=2.0, alpha=0.01, beta=0.01, logscale=false, method=:owenq)
    @test ClinicalTrialUtilities.besamplen(t.task).result == 4810
    Base.show(io, t)
    Base.show(io, t.task)
    #@test n == 4810 && round(p, digits=7) == 0.9900151
    t = ClinicalTrialUtilities.besamplen(cv=0.347)
    @test ClinicalTrialUtilities.besamplen(t.task).result == 52
    Base.show(io, t)
    Base.show(io, t.task)
    #@test n == 52 && round(p, digits=7) == 0.8136415
    t = ClinicalTrialUtilities.besamplen(;theta0=1.05, theta1=0.9, theta2=1.25, cv=0.0001, alpha=0.05, beta=0.15, logscale=true, method=:nct, design=:parallel)
    @test ClinicalTrialUtilities.besamplen(t.task).result == 4
    Base.show(io, t)
    Base.show(io, t.task)
    #@test n == 4 && p ≈ 1.0
    t = ClinicalTrialUtilities.besamplen(;theta0=1.0, theta1=0.95, theta2=1.5, cv=0.8, alpha=0.0000001, beta=0.0001, logscale=true, method=:shifted, design=:d2x2x4)
    @test ClinicalTrialUtilities.besamplen(t.task).result == 10002
    Base.show(io, t)
    Base.show(io, t.task)
    #@test n == 10002 && p ≈ 0.9818179411719451
    t = ClinicalTrialUtilities.besamplen(;theta0=1.05, theta1=0.8, theta2=1.25, cv=0.8, alpha=0.05, beta=0.2, logscale=true, method=:shifted, design=:d2x2x4)
    @test ClinicalTrialUtilities.besamplen(t.task).result == 106
    Base.show(io, t)
    Base.show(io, t.task)
    #@test n == 106 && p ≈ 0.8060551186037984
    t = ClinicalTrialUtilities.besamplen(;theta0=0.95, theta1=0.8, theta2=1.25, cv=0.35, alpha=0.1, beta=0.1, logscale=true, method=:shifted, design=:parallel)
    @test ClinicalTrialUtilities.besamplen(t.task).result == 106
    Base.show(io, t)
    Base.show(io, t.task)
    #@test n == 106 && p ≈ 0.9013894463164578

    #CHECK with R
    t = ClinicalTrialUtilities.besamplen(;cv=0.35, design=:d4x4)
    @test ClinicalTrialUtilities.besamplen(t.task).result == 52
    Base.show(io, t)
    Base.show(io, t.task)
end

println(" ---------------------------------- ")
@testset "#4  bepower Test        " begin
    #parallel
    #1
    t = ClinicalTrialUtilities.bepower(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.3, n=31, design=:parallel, method=:owenq)
    @test ClinicalTrialUtilities.bepower(t.task).result ≈ 0.2949476 atol=1E-7
    @test ClinicalTrialUtilities.ctpower(t.task).result ≈ 0.2949476 atol=1E-7
    Base.show(io, t)
    Base.show(io, t.task)
    #2
    t = ClinicalTrialUtilities.bepower(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.3, n=32, design=:parallel, method=:owenq)
    @test ClinicalTrialUtilities.bepower(t.task).result ≈ 0.3166927 atol=1E-7
    Base.show(io, t)
    Base.show(io, t.task)
    #2x2
    #3
    t = ClinicalTrialUtilities.bepower(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.2, n=20, design=:d2x2, method=:owenq)
    @test ClinicalTrialUtilities.bepower(t.task).result         ≈ 0.8346802 atol=1E-7
    Base.show(io, t)
    Base.show(io, t.task)
    #4
    t = ClinicalTrialUtilities.bepower(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.2, n=10, design=:d2x2, method=:nct)
    @test ClinicalTrialUtilities.bepower(t.task, method=:nct).result        ≈ 0.4316618 atol=1E-7
    Base.show(io, t)
    Base.show(io, t.task)
    #5
    t = ClinicalTrialUtilities.bepower(alpha=0.1, logscale=false, theta1=-0.1, theta2=0.1, theta0=0.0, sd=0.14, n=21, design=:d2x2, method=:shifted)
    @test ClinicalTrialUtilities.bepower(t.task, method=:shifted).result       ≈ 0.6626132 atol=1E-7
    Base.show(io, t)
    Base.show(io, t.task)
    #6
    t = ClinicalTrialUtilities.bepower(alpha=0.05, logscale=false, theta1=-0.1, theta2=0.1, theta0=0.0, sd=0.14, n=30, design=:d2x2, method=:nct)
    @test ClinicalTrialUtilities.bepower(t.task, method=:nct).result          ≈ 0.7079951 atol=1E-7
    Base.show(io, t)
    Base.show(io, t.task)
    #7
    t = ClinicalTrialUtilities.bepower(alpha=0.0000001, logscale=false, theta1=-0.1, theta2=0.1, theta0=0.0, sd=1.0, n=10000, design=:d2x2, method=:owenq)
    @test ClinicalTrialUtilities.bepower(t.task).result ≈ 0.9380914 atol=1E-7
    Base.show(io, t)
    Base.show(io, t.task)
    #8
    t = ClinicalTrialUtilities.bepower(alpha=0.0001, logscale=false, theta1=-0.1, theta2=0.1, theta0=0.0, sd=1.0, n=3500, design=:d2x2, method=:owenq)
    @test ClinicalTrialUtilities.bepower(t.task).result     ≈ 0.3545904 atol=1E-7
    Base.show(io, t)
    Base.show(io, t.task)
    #9
    t = ClinicalTrialUtilities.bepower(alpha=0.00000005, logscale=false, theta1=-0.1, theta2=0.1, theta0=0.0, sd=1.5, n=20000, design=:d2x2, method=:owenq)
    @test ClinicalTrialUtilities.bepower(t.task).result ≈ 0.8197361 atol=1E-7
    Base.show(io, t)
    Base.show(io, t.task)
    #10
    t = ClinicalTrialUtilities.bepower(alpha=0.1, logscale=false, theta1=-0.1, theta2=0.1, theta0=0.0, sd=0.14, n=4, design=:d2x2, method=:shifted)
    @test ClinicalTrialUtilities.bepower(t.task, method=:shifted).result ≈ 0.0
    Base.show(io, t)
    Base.show(io, t.task)
    #11
    t = ClinicalTrialUtilities.bepower(alpha=0.1, logscale=false, theta1=-0.1, theta2=0.1, theta0=0.0, sd=0.02, n=3, design=:d2x2, method=:shifted)
    @test ClinicalTrialUtilities.bepower(t.task, method=:shifted).result ≈ 0.7738659 atol=1E-7
    Base.show(io, t)
    Base.show(io, t.task)
    #
    t = ClinicalTrialUtilities.bepower(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.2, n=27, design=:d2x2, method=:owenq)
    @test ClinicalTrialUtilities.bepower(t.task).result ≈ 0.9264365737448076
    t = ClinicalTrialUtilities.bepower(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.2, n=29, design=:d2x2, method=:owenq)
    @test ClinicalTrialUtilities.bepower(t.task).result ≈ 0.941900827163551
    t = ClinicalTrialUtilities.bepower(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.2, n=31, design=:d2x2, method=:owenq)
    @test ClinicalTrialUtilities.bepower(t.task).result ≈ 0.9542152686694777

    #2x2x4
    t = ClinicalTrialUtilities.bepower(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.4, n=35, design=:d2x2x4)
    @test ClinicalTrialUtilities.bepower(t.task).result ≈ 0.829747  atol=1E-6
    t = ClinicalTrialUtilities.bepower(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=1, n=35, design=:d2x2x4)
    @test ClinicalTrialUtilities.bepower(t.task).result ≈ 0.014249535210231756  atol=1E-6
    #2x4x4
    t = ClinicalTrialUtilities.bepower(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.4, n=35, design=:d2x4x4)
    @test ClinicalTrialUtilities.bepower(t.task).result ≈ 0.8291076  atol=1E-7
    t = ClinicalTrialUtilities.bepower(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.4, n=34, design=:d2x4x4)
    @test ClinicalTrialUtilities.bepower(t.task).result ≈ 0.8180596  atol=1E-7
    t = ClinicalTrialUtilities.bepower(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.4, n=33, design=:d2x4x4)
    @test ClinicalTrialUtilities.bepower(t.task).result ≈ 0.8069565  atol=1E-7
    #2x3x3
    t = ClinicalTrialUtilities.bepower(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.4, n=32, design=:d2x3x3)
    @test ClinicalTrialUtilities.bepower(t.task).result ≈ 0.5976873  atol=1E-7
    t = ClinicalTrialUtilities.bepower(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.4, n=31, design=:d2x3x3)
    @test ClinicalTrialUtilities.bepower(t.task).result ≈ 0.579468  atol=1E-6
    t = ClinicalTrialUtilities.bepower(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.4, n=30, design=:d2x3x3)
    @test ClinicalTrialUtilities.bepower(t.task).result ≈ 0.5614358  atol=1E-7
    #3x3
    t = ClinicalTrialUtilities.bepower(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.4, n=30, design=:d3x3)
    @test ClinicalTrialUtilities.bepower(t.task).result ≈ 0.3847067  atol=1E-7
    #3x6x3
    t = ClinicalTrialUtilities.bepower(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.4, n=30, design=:d3x6x3)
    @test ClinicalTrialUtilities.bepower(t.task).result ≈ 0.3847067  atol=1E-7
    #2x4x2
    t = ClinicalTrialUtilities.bepower(alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.4, n=30, design=:d2x4x2)
    @test ClinicalTrialUtilities.bepower(t.task).result ≈ 0.0001785665  atol=1E-10
    #@test ClinicalTrialUtilities.ctpower(t.task).result ≈ 0.0001785665  atol=1E-10
end

include("utils.jl")

include("citest.jl")

include("dstest.jl")

include("freque.jl")


println(" ---------------------------------- ")
@testset "#9  Random              " begin
    rng = StableRNG(0)
    @test ClinicalTrialUtilities.randomseq(seed = 1234, rng = rng) == [1, 2, 1, 2, 2, 1, 1, 2, 2, 1]
    rdf   = ClinicalTrialUtilities.randomtable(seed = 1234, rng = rng)
    @test rdf[:Group] == [1, 2, 1, 2, 2, 1, 1, 2, 2, 1]
end


#PK
include("pktest.jl")


include("sim.jl")


println(" ---------------------------------- ")
@testset "  Scenario              " begin
    #Pharmacokinetics statistics
    pkds = ClinicalTrialUtilities.pkimport(pkdata2, [:Subject, :Formulation]; time = :Time, conc = :Concentration)
    pk   = ClinicalTrialUtilities.nca!(pkds)
    df   = ClinicalTrialUtilities.DataFrame(pk; unst = true)
    ds   = ClinicalTrialUtilities.descriptive(df, stats = [:n, :mean, :sd], sort = [:Formulation])
    df   = ClinicalTrialUtilities.DataFrame(ds; unst = true)
    Base.show(io, pkds)
    Base.show(io, pkds[1])
    Base.show(io, pk)
    Base.show(io, pk[1])
    Base.show(io, ds)
    Base.show(io, ds[1])


    @test ds[Dict(:Variable=>:AUClast,:Formulation=>"R")][:mean] ≈ 7431.283916666667

    pdds = ClinicalTrialUtilities.pdimport(pkdata2, [:Subject, :Formulation]; time = :Time, resp = :Concentration)
    pd   = ClinicalTrialUtilities.nca!(pdds)
    df   = ClinicalTrialUtilities.DataFrame(pd; unst = true)
    ds   = ClinicalTrialUtilities.descriptive(df, stats = :default, sort = [:Formulation])

    @test ds[Dict(:Variable=>:AUCBLNET,:Formulation=>"T")][:mean] ≈ 8607.09
    @test ds[(:Variable=>:AUCBLNET,:Formulation=>"R")][:sd] ≈ 1919.7954123986283

end

println(" ---------------------------------- ")
@testset "  HTML export           " begin

    pkds = ClinicalTrialUtilities.pkimport(pkdata2, [:Subject, :Formulation]; time = :Time, conc = :Concentration)
    pk   = ClinicalTrialUtilities.nca!(pkds)
    df   = ClinicalTrialUtilities.DataFrame(pk; unst = true)
    ClinicalTrialUtilities.htmlexport(df; io = io, sort = :Subject, rspan=:all, title = "PK Parameters", dict = :pk)
    @test true
end

println(" ---------------------------------- ")
@testset "  Errors                " begin

    #ERROR: ArgumentError: sampleSize: alpha and beta sould be > 0 and < 1.
    @test_throws ArgumentError ClinicalTrialUtilities.ctsamplen(param=:mean, type=:ea, group=:one, alpha=2, beta=0.2, diff=1, sd=1, a=1, b=1, k=1)
    #ERROR: ArgumentError: Diiference can't be ≤ 0.0 with Equivalence hypothesis!
    @test_throws ArgumentError ClinicalTrialUtilities.ctsamplen(param=:mean, type=:ei, group=:one, alpha=0.5, beta=0.2, diff=0, sd=1, a=1, b=1, k=1)
    #ERROR: ArgumentError: SD can't be ≤ 0.0!
    @test_throws ArgumentError ClinicalTrialUtilities.ctsamplen(param=:mean, type=:ns, group=:one, alpha=0.5, beta=0.2, diff=1, sd=0, a=1, b=1, k=1)
    #ERROR: ArgumentError: Constant k can't be ≤ 0.0!
    @test_throws ArgumentError ClinicalTrialUtilities.ctsamplen(param=:mean, type=:ea, group=:one, alpha=0.05, beta=0.2, diff=1, sd=1, a=0, b=0, k=0)
    #ERROR: ArgumentError: Design not known!
    @test_throws ArgumentError ClinicalTrialUtilities.Design(:ddd)

    #ERROR: ArgumentError: !(theta2 > thetao > theta1), check settings!
    @test_throws ArgumentError ClinicalTrialUtilities.besamplen(;cv=0.35, theta1 = 1.0)
    #ERROR: ArgumentError: Beta ≥ 1.0 or ≤ 0.0!
    @test_throws ArgumentError ClinicalTrialUtilities.besamplen(;cv=0.35, beta = 1.0)

    @test_throws ArgumentError ClinicalTrialUtilities.besamplen(;cv=0.35, alpha = 1.0)

    @test_throws ArgumentError ClinicalTrialUtilities.besamplen(;cv=-0.35)
    #CI
    ci = ClinicalTrialUtilities.ConfInt(1, 3, 2, 0.05)
    @test_throws ArgumentError ci[5]
    @test_throws ArgumentError ci = ClinicalTrialUtilities.propci(38, 100, alpha=0.05, method=:waldd)
    @test_throws ArgumentError ci = ClinicalTrialUtilities.diffpropci(30, 100, 40, 90; alpha=-0.05, method=:mee)
    @test_throws ArgumentError ci = ClinicalTrialUtilities.diffpropci(30, 100, 40, 90; alpha=0.05, method=:meee)
    @test_throws ArgumentError ci = ClinicalTrialUtilities.rrpropci(30, 100, 40, 90; alpha=0.05, method=:mnnn)
    @test_throws ArgumentError ci = ClinicalTrialUtilities.orpropci(30, 100, 40, 90; alpha=0.05,  method=:awoolfff)

    #DS
    pkds = ClinicalTrialUtilities.pkimport(pkdata, [:Subject, :Formulation]; conc = :Concentration, time = :Time)
    @test_throws BoundsError pkds[20]
    @test_throws ArgumentError pkds[:Subject => 20]
    @test_throws ArgumentError pkds[Dict(:Subject => 20)]


    #Type check
    ClinicalTrialUtilities.besamplen(;theta0=0, theta1=-1, theta2=1, sd=2, logscale = false)
end

println(" ---------------------------------- ")
@testset "  Experimental          " begin
    @test ClinicalTrialUtilities.bartlettstest(10, 100, 12, 50)[1] ≈ 0.4578814034969132

    @test  ClinicalTrialUtilities.meanratiot(10, 2, 100, 12, 3, 120, 0.2, 0.05)[1] ≈ 0.7623819718221434

    @test  ClinicalTrialUtilities.diffmcnmwaldci(10, 20, 11, 15; alpha = 0.05).lower ≈ -0.029553406052837572

    @test  ClinicalTrialUtilities.diffmcnmwaldccci(10, 20, 11, 15; alpha = 0.05).lower ≈ -0.04741054890998043

    @test  ClinicalTrialUtilities.mcnmtest(ClinicalTrialUtilities.McnmConTab([12 23; 22 33]); cc = false) ≈ 0.022222222222222223

    #= /*SAS Tesing code*/
data effect;
length Group $12 Response $3;
input Group Response N;
datalines;
T          Yes  75
T          No   25
R         Yes  80
R         No   20
;

proc freq data=effect order=data;
   weight N;
   tables Group*Response / riskdiff (EQUIV MARGIN=0.08) alpha=0.05;
run;

proc freq data=effect order=data;
   weight N;
   tables Group*Response / riskdiff (EQUAL(NULL=0)) alpha=0.05;
run;

proc freq data=effect order=data;
   weight N;
   tables Group*Response / riskdiff (NONINF MARGIN=0.15) alpha=0.05;
run;
=#
    param = ClinicalTrialUtilities.DiffProportion(ClinicalTrialUtilities.Proportion(75, 100), ClinicalTrialUtilities.Proportion(80, 100))
    hyp   = ClinicalTrialUtilities.Equivalence(-0.08, 0.08, 0.05)
    hyp2 = ClinicalTrialUtilities.Equality(0.0, 0.05)
    hyp3 = ClinicalTrialUtilities.Superiority(-0.15, 0.05)
    @test  ClinicalTrialUtilities.pval(param, hyp; method = :wald, atol = 1E-6) ≈ 0.3054046630859374
    @test  ClinicalTrialUtilities.pval(param, hyp2; method = :wald, atol = 1E-6) ≈ 0.39634132385253895
    @test  ClinicalTrialUtilities.pval(param, hyp3; method = :wald, atol = 1E-6) ≈ 0.04490566253662108


end
#println(" ---------------------------------- ")
#@testset "  Tpl                 " begin
#end
