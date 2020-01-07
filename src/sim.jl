
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

# Proportion, Superiority, Power
"""
    ctsim(t::CTask{DiffProportion{P, P}, D, H, Power}; nsim = 100, seed=0, method = :default, dropout = 0.0)  where P <: Proportion where D where H <: Superiority

Proportion difference power simulation.

"""
function ctsim(t::CTask{DiffProportion{P, P}, D, H, Power}; nsim = 100, seed=0, method = :default, dropout = 0.0)  where P <: Proportion where D where H <: Superiority
    if seed != 0  Random.seed!(seed) end

    pow     = 0
    D₁      = Binomial(t.param.a.n, getval(t.param.a))
    D₂      = Binomial(t.param.b.n, getval(t.param.b))
    for i = 1:nsim
        n1      = rand(D₁)
        n2      = rand(D₂)
        ci      = confint(DiffProportion(Proportion(n1, t.param.a.n), Proportion(n2, t.param.b.n)); level = 1.0 - t.alpha * 2, method = method)
        if ci.lower > t.hyp.diff pow += 1 end
    end
    return pow/nsim
end
