
"""
    besim(t::CTask{T, D, Bioequivalence, Power}; nsim = 100, seed=0)  where T where D

Bioequivalence power simulation.

"""
function besim(t::CTask{P, D, Bioequivalence, Power}; nsim::Int = 100, rsabe::Bool = false, rsabeconst = 0.760, seed = 0, rng = MersenneTwister())  where P where D
    if seed != 0  Random.seed!(rng, seed) end
    df      = t.design.df(t.objective.val)
    sef     = sediv(t.design, t.objective.val)
    CHSQ    = Chisq(df)
    tval    = quantile(TDist(df), 1 - t.hyp.alpha)
    σ̵ₓ      = sef * t.param.a.sd
    pow     = 0
    for i = 1:nsim
        smean   = rand(rng, ZDIST) * σ̵ₓ + (t.param.a.m - t.param.b.m)
        σ²      = rand(rng, CHSQ)  * t.param.a.sd  ^ 2 / df
        σ       = sqrt(σ²)
        hw      = tval * σ * sef
        θ₁      = smean - hw
        θ₂      = smean + hw
        if rsabe
            cv = cvfromsd(σ)
            if cv > 0.30
                if θ₁ > (- rsabeconst * σ) && θ₂ < (rsabeconst * σ) pow += 1 end
            else
                if θ₁ > t.hyp.llim && θ₂ < t.hyp.ulim pow += 1 end
            end
        else
            if θ₁ > t.hyp.llim && θ₂ < t.hyp.ulim pow += 1 end
        end
    end
    return pow/nsim
end

# Proportion, Superiority, Power
"""
    ctsim(t::CTask{DiffProportion{P, P}, D, H, Power}; nsim = 100, seed=0, method = :default, dropout = 0.0)  where P <: Proportion where D where H <: Superiority

Proportion difference power simulation.

"""
function ctsim(t::CTask{DiffProportion{P1, P2}, D, H, Power}; nsim = 100, method = :default, dropout = 0.0, seed = 0, rng = MersenneTwister())  where P1 <: Proportion where P2 <: Proportion where D where H <: AbstractHypothesis
    if seed != 0  Random.seed!(rng, seed) end
    n₁      = getval(t.objective)
    n₂      = Int(ceil(getval(t.objective) / t.k))
    pow     = 0
    D₁      = Binomial(n₁, getval(t.param.a))
    D₂      = Binomial(n₂, getval(t.param.b))
    for i = 1:nsim
        xn₁      = rand(rng, D₁)
        xn₂      = rand(rng, D₂)
        if  checkhyp(t.hyp, DiffProportion(Proportion(xn₁, n₁), Proportion(xn₂, n₂)); method = method) pow += 1 end
    end
    return pow/nsim
end


#=
function ctsimadapt(t::CTask{DiffProportion{P, P}, D, H, Power}; maxsize = 200, nsim = 100, seed=0, method = :default, dropout = 0.0)  where P <: Proportion where D where H
    if seed != 0  Random.seed!(seed) end
    n₁      = getval(t.objective)
    n₂      = Int(ceil(getval(t.objective) / t.k))
    pow     = 0
    D₁      = Binomial(n₁, getval(t.param.a))
    D₂      = Binomial(n₂, getval(t.param.b))
    for i = 1:nsim
        n1      = rand(D₁)
        n2      = rand(D₂)
        dp      = DiffProportion(Proportion(n1, getval(t.objective)), Proportion(n2, getval(t.objective)))

        if  checkhyp(t.hyp, dp; method = method)
            pow += 1
        else
            #Recalc
            task     = CTask(dp, t.design, t.hyp, SampleSize(0.2), t.k)
            #Recalc for total n1
            recalc   = ctsamplen(task)
            #println(recalc)
            if recalc.result > maxsize continue end
            n1t      = Int(ceil(recalc.result))
            #Calc for n2
            n2t      = Int(ceil(n1t / t.k))

            #Add
            n1a      = n1t - n₁
            n2a      = n2t - n₂

            if n1a <= 0 || n2a <= 0
                continue
            end
            println("New N: $n1t")

            n1p      = rand(Binomial(n1a, getval(dp.a)))
            n2p      = rand(Binomial(n2a, getval(dp.b)))

            dp      = DiffProportion(Proportion(n1p, n1t), Proportion(n2p, n2t))

            if  checkhyp(t.hyp, dp; method = method)
                pow += 1
            end
        end
    end
    return pow/nsim
end
=#
#=
#Methods for adaptive design simulation
function checksatage(task, diffp; method)
end
function futilitie(task, diffp; method)
end
=#



#Means
"""
    ctsim(t::CTask{DiffMean{M, M}, D, H, Power}; nsim = 100, seed=0, method = :default, dropout = 0.0)  where M <: Mean where D where H

Means power simulation.

"""
function ctsim(t::CTask{P, D, H, Power}; nsim = 100, method = :default, dropout = 0.0, seed = 0, rng = MersenneTwister())  where P <: DiffMean where D where H
    if seed != 0  Random.seed!(rng, seed) end
    n₁      = getval(t.objective)
    n₂      = Int(ceil(getval(t.objective) / t.k))
    pow     = 0

    for i = 1:nsim
        n1      = rand(rng, ZDIST, n₁) .* t.param.a.sd .+ getval(t.param.a)
        n2      = rand(rng, ZDIST, n₂) .* t.param.b.sd .+ getval(t.param.b)
        if  checkhyp(t.hyp, DiffMean(Mean(n1), Mean(n2)); method = method) pow += 1 end
    end
    return pow/nsim
end
