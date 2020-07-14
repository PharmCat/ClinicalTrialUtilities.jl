println(" ---------------------------------- ")
@testset "#1  ctsamplen Test      " begin
    #Sample Size Calculations in Clinical Research, Second Edition, Shein-Chung Chow, Ph.D., 2008, International Standard Book Number‑13: 978‑1‑58488‑982‑3
    #1
    t = ClinicalTrialUtilities.ctsamplen(param=:mean, type=:ea, group=:one, alpha=0.05, beta=0.2, sd=1, a=2, b=1.5)
    @test ceil(t.result) == ceil(ClinicalTrialUtilities.ctsamplen(t.task).result) == 32
    Base.show(io, t)
    Base.show(io, t.task)
    #2
    t = ClinicalTrialUtilities.ctsamplen(param=:mean, type=:ei, group=:one, alpha = 0.05 / 2, beta=0.2, sd=0.1, diff=0.05, a=2, b=2)
    @test ceil(t.result) == ceil(ClinicalTrialUtilities.ctsamplen(t.task).result) == 35
    Base.show(io, t)
    #3
    t = ClinicalTrialUtilities.ctsamplen(param=:mean, type=:ns, group=:one, alpha=0.05, beta=0.2, sd=1, diff=-0.5, a=2, b=1.5)
    @test ceil(t.result) == ceil(ClinicalTrialUtilities.ctsamplen(t.task).result) == 7
    Base.show(io, t)
    #4
    t = ClinicalTrialUtilities.ctsamplen(param=:mean, type=:ea, group=:two, alpha=0.05, beta=0.2, sd=10, a=5, b=10, k=1)
    @test ceil(t.result) == ceil(ClinicalTrialUtilities.ctsamplen(t.task).result) == 63
    Base.show(io, t)
    #5
    t = ClinicalTrialUtilities.ctsamplen(param=:mean, type=:ei, group=:two, alpha=0.05 / 2, beta=0.2, sd=10, diff=5, a=5, b=4, k=1)
    @test ceil(t.result) == ceil(ClinicalTrialUtilities.ctsamplen(t.task).result) == 108
    Base.show(io, t)
    #6
    t = ClinicalTrialUtilities.ctsamplen(param=:mean, type=:ns, group=:two, alpha=0.05, beta=0.2, sd=10, diff=5, a=5, b=5, k=1)
    @test ceil(t.result) == ceil(ClinicalTrialUtilities.ctsamplen(t.task).result) == 50
    Base.show(io, t)
    #7
    t = ClinicalTrialUtilities.ctsamplen(param=:prop, type=:ea, group=:one, alpha=0.05, beta=0.2, a=0.5, b=0.3)
    @test ceil(t.result) == ceil(ClinicalTrialUtilities.ctsamplen(t.task).result) == 50
    Base.show(io, t)
    #8
    t = ClinicalTrialUtilities.ctsamplen(param=:prop, type=:ei, group=:one, alpha=0.05 / 2, beta=0.2, diff=0.2, a=0.6, b=0.6)
    @test ceil(t.result) == ceil(ClinicalTrialUtilities.ctsamplen(t.task).result) == 52
    Base.show(io, t)
    #9
    t = ClinicalTrialUtilities.ctsamplen(param=:prop, type=:ns, group=:one, alpha=0.05, beta=0.2, diff=-0.1, a=0.5, b=0.3)
    @test ceil(t.result) == ceil(ClinicalTrialUtilities.ctsamplen(t.task).result) == 18
        Base.show(io, t)
    #10    p.92
    t = ClinicalTrialUtilities.ctsamplen(param=:prop, type=:ea, group=:two, alpha=0.05, beta=0.2, a=0.65, b=0.85)
    @test ceil(t.result) == ceil(ClinicalTrialUtilities.ctsamplen(t.task).result) == 70
    Base.show(io, t)
    #11    p.93
    t = ClinicalTrialUtilities.ctsamplen(param=:prop, type=:ei, group=:two, alpha=0.05 / 2, beta=0.2, diff=0.2, a=0.75, b=0.80)
    @test ceil(t.result) == ceil(ClinicalTrialUtilities.ctsamplen(t.task).result) == 133
    Base.show(io, t)
    #12    p.92
    t = ClinicalTrialUtilities.ctsamplen(param=:prop, type=:ns, group=:two, alpha=0.05, beta=0.2, diff=-0.1, a=0.85, b=0.65)
    @test ceil(t.result) == ceil(ClinicalTrialUtilities.ctsamplen(t.task).result) == 25
    @test ceil(ClinicalTrialUtilities.pdiffnsn(0.85, 0.65, -0.1).result) == 25
    Base.show(io, t)
    #13    p.108
    t = ClinicalTrialUtilities.ctsamplen(param=:or, type=:ea,  alpha=0.05, beta=0.2, a=0.4, b=0.25)
    @test ceil(t.result) == ceil(ClinicalTrialUtilities.ctsamplen(t.task).result) == 156
    Base.show(io, t)
    #14.1    p.109
    t = ClinicalTrialUtilities.ctsamplen(param=:or, type=:ei,  alpha=0.05 / 2, beta=0.2, diff=0.5, a=0.25, b=0.25)
    @test ceil(t.result) == ceil(ClinicalTrialUtilities.ctsamplen(t.task).result) == 366
    Base.show(io, t)
    #15.1    p.108
    t = ClinicalTrialUtilities.ctsamplen(param=:or, type=:ns,  alpha=0.05, beta=0.2, diff=0.2, a=0.4, b=0.25)
    @test ceil(t.result) == ceil(ClinicalTrialUtilities.ctsamplen(t.task).result) == 242
    Base.show(io, t)
    #14.2
    t = ClinicalTrialUtilities.ctsamplen(param=:or, type=:ei,  alpha=0.05 / 2, beta=0.2, diff=exp(0.5), a=0.25, b=0.25, logscale = false)
    @test ceil(t.result) == 366
    Base.show(io, t)
    #15.2
    t = ClinicalTrialUtilities.ctsamplen(param=:or, type=:ns,  alpha=0.05, beta=0.2, diff=exp(0.2), a=0.4, b=0.25, logscale = false)
    @test ceil(t.result) == 242
    Base.show(io, t)
    #16
    t = ClinicalTrialUtilities.ctsamplen(param=:prop, type=:mcnm, a=0.45, b=0.05)
    @test ceil(t.result) == ceil(ClinicalTrialUtilities.ctsamplen(t.task).result) == 23
    Base.show(io, t)

    #COX
    #ClinicalTrialUtilities.cox_superiority(2, 1, 0.8, 0.05/2, 0.2, 1) ≈  81.68207407358884
    #ClinicalTrialUtilities.cox_equivalence(1, exp(0.5), 0.8, 0.05, 0.2, 1) ≈ 171.27694701335946
    #ClinicalTrialUtilities.cox_equality(2, 1, 0.8, 0.05, 0.2, 1) ≈ 81.68207407358884

    #Additional
    #Different type input
    @test ClinicalTrialUtilities.ctsamplen(param=:mean, type=:ns, group=:two, alpha=0.05, beta=0.2, diff=1, sd=20.0, a=1, b=2).result ≈ 1236.511446403953 atol=1E-12

    @test ceil(ClinicalTrialUtilities.ctsamplen(param=:prop, type=:ns, group=:two, alpha=0.05, beta=0.2, diff=0.05, a=0.85, b=0.65).result) == 98
end

println(" ---------------------------------- ")
@testset "#2  ctpower Test        " begin

    #1
    t = ClinicalTrialUtilities.ctpower(param=:mean, type=:ea, group=:one, a=2, b=1.5, sd=1,n=32, alpha=0.05)
    @test ClinicalTrialUtilities.ctpower(t.task).result ≈ 0.8074304194325561  atol=1E-7
    Base.show(io, t)
    Base.show(io, t.task)
    #2
    t = ClinicalTrialUtilities.ctpower(param=:mean, type=:ei, group=:one, a=2, b=2, sd=0.1, diff=0.05, n=35, alpha=0.05 / 2)
    @test ClinicalTrialUtilities.ctpower(t.task).result ≈ 0.8108839754376387  atol=1E-7
    Base.show(io, t)
    #3
    t = ClinicalTrialUtilities.ctpower(param=:mean, type=:ns, group=:one, a=2, b=1.5, sd=1, diff=-0.5, n=7, alpha=0.05)
    @test ClinicalTrialUtilities.ctpower(t.task).result ≈ 0.8415707712023641  atol=1E-7
    Base.show(io, t)
    #4
    t = ClinicalTrialUtilities.ctpower(param=:mean, type=:ea, group=:two, a=5, b=10, sd=10, n=63, alpha=0.05)
    @test ClinicalTrialUtilities.ctpower(t.task).result ≈ 0.8013023941055788  atol=1E-7
    Base.show(io, t)
    #5
    t = ClinicalTrialUtilities.ctpower(param=:mean, type=:ei, group=:two, a=5, b=4, sd=10, diff=5, n=108, alpha=0.05 / 2)
    @test ClinicalTrialUtilities.ctpower(t.task).result ≈ 0.80452354556060    atol=1E-7
    Base.show(io, t)
    #6
    t = ClinicalTrialUtilities.ctpower(param=:mean, type=:ns, group=:two, a=5, b=5, sd=10, diff=5, n=50, alpha=0.05)
    @test ClinicalTrialUtilities.ctpower(t.task).result ≈ 0.8037819415575257  atol=1E-7
    Base.show(io, t)
    #7
    t = ClinicalTrialUtilities.ctpower(param=:prop, type=:ea, group=:one, a=0.5, b=0.3, n=50, alpha=0.05)
    @test ClinicalTrialUtilities.ctpower(t.task).result ≈ 0.8074304194325561  atol=1E-7
    Base.show(io, t)
    #8
    t = ClinicalTrialUtilities.ctpower(param=:prop, type=:ei, group=:one, a=0.6, b=0.6, diff=0.2, n=52, alpha=0.05 / 2)
    @test ClinicalTrialUtilities.ctpower(t.task).result ≈ 0.8060834056011101  atol=1E-7
    Base.show(io, t)
    #9
    t = ClinicalTrialUtilities.ctpower(param=:prop, type=:ns, group=:one, a=0.5, b=0.3, diff=-0.1, n=18, alpha=0.05)
    @test ClinicalTrialUtilities.ctpower(t.task).result ≈ 0.8161481827204281  atol=1E-7
    Base.show(io, t)
    #10
    t = ClinicalTrialUtilities.ctpower(param=:prop, type=:ea, group=:two, a=0.65, b=0.85, n=70, alpha=0.05)
    @test ClinicalTrialUtilities.ctpower(t.task).result ≈ 0.8019139068528     atol=1E-7
    Base.show(io, t)
    #11
    t = ClinicalTrialUtilities.ctpower(param=:prop, type=:ei, group=:two, a=0.65, b=0.85, diff=0.05, n=136, alpha=0.05 / 2)
    @test ClinicalTrialUtilities.ctpower(t.task).result ≈ 0.8033294052407269  atol=1E-7
    Base.show(io, t)
    #12
    t = ClinicalTrialUtilities.ctpower(param=:prop, type=:ns, group=:two, a=0.85, b=0.65, diff=-0.1, n=25, alpha=0.05)
    @test ClinicalTrialUtilities.ctpower(t.task).result ≈ 0.808599833380679   atol=1E-7
    Base.show(io, t)
    #13
    t = ClinicalTrialUtilities.ctpower(param=:or, type=:ea, a=0.4, b=0.25, n=156, alpha=0.05)
    @test ClinicalTrialUtilities.ctpower(t.task).result ≈ 0.8020239054864792  atol=1E-7
    Base.show(io, t)
    #14
    t = ClinicalTrialUtilities.ctpower(param=:or, type=:ei, a=0.25, b=0.25, diff=0.5, n=366, alpha=0.05 / 2)
    @test ClinicalTrialUtilities.ctpower(t.task).result ≈ 0.8008593380478983  atol=1E-7
    Base.show(io, t)
    #15
    t = ClinicalTrialUtilities.ctpower(param=:or, type=:ns, a=0.4, b=0.25, diff=0.2, n=242, alpha=0.05)
    @test ClinicalTrialUtilities.ctpower(t.task).result ≈ 0.8007200876001626  atol=1E-7
    Base.show(io, t)
    #16
    t = ClinicalTrialUtilities.ctpower(param=:prop, type=:mcnm, a=0.45, b=0.05, n=23, alpha=0.1)
    @test ClinicalTrialUtilities.ctpower(t.task).result ≈ 0.9023805225074147  atol=1E-7
    Base.show(io, t)

    t = ClinicalTrialUtilities.ctpower(param=:or, type=:ns, a=0.4, b=0.25, diff=exp(0.2), n=242, alpha=0.05, logscale = false)
    @test ClinicalTrialUtilities.ctpower(t.task).result ≈ 0.8007200876001626  atol=1E-7
    Base.show(io, t)
end
