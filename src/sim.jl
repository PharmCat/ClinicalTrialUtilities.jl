
"""
    ctsim(t::CTask{T, D, Bioequivalence, Power}; nsim = 100, seed=0)  where T where D

Bioequivalence power simulation.

"""
function ctsim(t::CTask{T, D, Bioequivalence, Power}; nsim = 100, seed=0)  where T where D
    if seed != 0  Random.seed!(seed) end
    df      = t.design.df(t.objective.val)
    sef     = sediv(t.design, t.objective.val)
    CHSQ    = Chisq(df)
    tval    = quantile(TDist(df), 1 - t.alpha)
    σ̵ₓ      = sef * t.param.a.sd
    pow     = 0
    for i = 1:nsim
        smean   = rand(ZDIST) * σ̵ₓ + (t.param.a.m - t.param.b.m)
        σ²      = rand(CHSQ)  * t.param.a.sd  ^ 2 / df
        hw      = tval * sqrt(σ²) * sef
        θ₁      = smean - hw
        θ₂      = smean + hw
        if θ₁ > t.hyp.llim && θ₂ < t.hyp.ulim pow += 1 end
    end
    return pow/nsim
end
