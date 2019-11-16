# Clinical Trial Utilities
# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)

module SIM
    using Distributions
    using Random

    import ..ctsamplen
    import ..varfromcv
    import ..sediv
    import ..Design
    import ..ZDIST
    import ..ConfInt
    #import ..CTUException
    import ..twoprop
    import ..diffmeanci
    import ..AbstractTask, ..CTask, ..AbstractParameter, ..AbstractObjective, ..AbstractHypothesis, ..SampleSize, ..Power, ..Proportion, ..Probability, ..DiffProportion
    import ..TaskResult

    function simdist(t::CTask{DiffProportion{Probability}, H, O}) where H <: AbstractHypothesis where O <: AbstractObjective
        return (Binomial(t.objective.val, t.param.a.p), Binomial(t.objective.val, t.param.b.p))
    end

    function randvar(dist::Binomial)
        return Proportion(rand(dist), dist.n)
    end
    function randvar(dist::Tuple{Vararg{Binomial}})
        a = Array{Proportion, 1}(undef, length(dist))
        for i = 1:length(dist)
            a[i] = Proportion(rand(dist[i]), dist[i].n)
        end
        return Tuple(a)
    end

    function confint(p1::Proportion, p2::Proportion; alpha = 0.05, method=:wald)
        return twoprop(p1.x, p1.n, p2.x, p2.n; alpha=alpha, type=:diff, method = method)
    end
    function confint(p::Tuple{Proportion, Proportion}; alpha = 0.05, method=:wald)
        return twoprop(p[1].x, p[1].n, p[2].x, p[2].n; alpha=alpha, type=:diff, method = method)
    end

    function ctpowersim(task::T; simnum=5, seed=0, f = (cl, tl, cu, tu) -> cl > tl && cu < tu) where T <: Task
        rng     = MersenneTwister(1234)
        if seed == 0  Random.seed!(rng) else Random.seed!(seed) end
        pow     = 0
        nsim    = 10^simnum
        d       = simdist(task)
        for i=1:nsim
            x = randvar(d)
            ci = confint(x; alpha = task.alpha)
            if f(ci.lower, task.llim, ci.upper, task.ulim) pow += 1 end
        end
        return TaskResult(task, :sim, pow/nsim)
    end





    function bepowersim(;alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.0, n=0, simnum=5, seed=0)
        if alpha <= 0.0 || alpha >= 1.0  throw(CTUException(1111,"SIM.bepower: alpha should be > 0 and < 1")) end
        if cv <= 0.0   throw(CTUException(1112,"SIM.bepower: cv should be > 0")) end
        if n <= 0   throw(CTUException(1113,"SIM.bepower: n should be > 0")) end
        if simnum <= 0   throw(CTUException(1114,"SIM.bepower: simnum should be > 0")) end

        d     = Design(:d2x2)                                   #dffunc if generic funtion with 1 arg return df
        df    = d.df(n)

        sef   = sediv(d, n)                                 #for se
        #ms    = log(1+cv^2)
        if logscale
            ltheta1 = log(theta1); ltheta2 = log(theta2); diffm = log(theta0); ms = log(1+cv^2);
        else
            ltheta1 = theta1; ltheta2 = theta2; diffm = theta0; ms = cv*cv;
        end
        return bepowerSIM(ltheta1, ltheta2, ms, diffm, df, sef, 10^simnum, alpha, seed=seed)
    end

    function bepowerSIM(theta1, theta2, ms, mean, df, sef, nsim, alpha; seed=0)
        rng  = MersenneTwister(1234)
        if seed == 0  Random.seed!(rng) else Random.seed!(seed) end
        CHSQ    = Chisq(df)
        tval    = quantile(TDist(df), 1-alpha)
        se      = sef*sqrt(ms)                                                  #se_diff
        pow     = 0
        for i = 1:nsim
            smean   = rand(ZDIST)*se+mean
            sms     = rand(CHSQ)*ms/df
            hw      = tval*sqrt(sms)*sef
            lwr     = smean - hw
            upr     = smean + hw
            if lwr > theta1 && upr < theta2 pow += 1 end
        end
        return pow/nsim
    end


    function ctpowersim(p1, n1, p2, n2, ref; alpha=0.05, type=:notdef, citype=:notdef, method=:notdef, simnum=5, seed=0)
        if type == :notdef || method == :notdef throw(CTUException(1115,"ctPropPower: type or method not defined.")) end
        rng = MersenneTwister(1234)
        if seed == 0  Random.seed!(rng) else Random.seed!(seed) end
        BIN1 = Binomial(n1, p1)
        BIN2 = Binomial(n2, p2)
        pow     = 0
        nsim = 10^simnum
        for i=1:nsim
            x1 = rand(BIN1)
            x2 = rand(BIN2)
            ci = twoprop(x1, n1, x2, n2; alpha=alpha, type=citype, method=method)
            if type == :ns
                if ci.lower > ref pow += 1 end
            elseif type == :ei
                if ci.lower > ref[1] &&  ci.upper < ref[2] pow += 1 end
            end
        end
        return pow/nsim
    end

    function ctsamplensim(p1, p2, ref; alpha=0.05, beta=0.2, type=:notdef, citype =:notdef, method=:notdef, simnum=5, seed=0)
        if type == :notdef || citype == :notdef || method == :notdef throw(CTUException(1116,"ctPropSampleN: type or method not defined.")) end
            st::Int = sn::Int = 10
        if citype == :diff
            st = sn = ceil(ctsamplen(param=:prop, type=type, group=:two, alpha=alpha, beta=beta, diff=ref, a=p1, b=p2).result)
        elseif citype == :or
            st = sn = ceil(ctsamplen(param=:or, type=type, group=:two, alpha=alpha/2, beta=beta, diff=ref, a=p1, b=p2).result)
        end

        pow = ctPropPower(p1, sn, p2, sn, ref; alpha=alpha, type=type, citype=citype, method=method, simnum=simnum, seed=seed)
        conu(x, y) = x < 1 - y
        cond(x, y) = x > 1 - y
        if pow < 1-beta
            inc = 1
            con = conu
        else
            con = cond
            inc = -1
        end
        pown = pow
        snn  = sn
        while con(pown, beta)
            pow = pown
            sn  = snn
            snn = snn + inc
            pown = ctPropPower(p1, snn, p2, snn, ref; alpha=alpha, type=type, citype=citype,method=method, simnum=simnum, seed=seed)
        end
        return snn, pown, sn, pow
    end

    function ctMeansPower(m1, s1, n1, m2, s2, n2, ref; alpha=0.05, method=:notdef, simnum=5, seed=0)
        #if type == :notdef || method == :notdef throw(CTUException(1115,"ctPropPower: type or method not defined.")) end

        rng = MersenneTwister(1234)
        if seed == 0  Random.seed!(rng) else Random.seed!(seed) end
        se1    = sqrt(s1/(n1-1))
        se2    = sqrt(s2/(n2-1))
        CHSQ1  = Chisq(n1-1)
        CHSQ2  = Chisq(n2-1)
        pow    = 0
        nsim   = 10^simnum
        for i=1:nsim
            sm1    = rand(ZDIST)*se1+m1
            sm2    = rand(ZDIST)*se2+m2
            ss1    = rand(CHSQ1)*s1/(n1-1)
            ss2    = rand(CHSQ2)*s2/(n2-1)
            ci = diffmeanci(sm1, ss1, n1, sm2, ss2, n2; alpha=alpha,  method=:ev)
            if ci.lower > ref pow += 1 end
        end
        return pow/nsim
    end

    function ctMeansPowerFS(m1, s1, n1, m2, s2, n2, ref;alpha=0.05, method=:notdef, simnum::Real=5, seed=0)
        #if type == :notdef || method == :notdef throw(CTUException(1115,"ctPropPower: type or method not defined.")) end
        if method == :notdef throw(CTUException(1117,"ctMeansPowerFS: method not defined.")) end
        rng = MersenneTwister(1234)
        if seed == 0  Random.seed!(rng) else Random.seed!(seed) end

        sd1    = sqrt(s1)
        sd2    = sqrt(s2)
        pow    = 0
        nsim   = convert(Int, 10^simnum)
        for i=1:nsim
            set1   = (rand(ZDIST, n1)*sd1).+m1
            set2   = (rand(ZDIST, n2)*sd2).+m2
            ci     = diffmeanci(mean(set1), var(set1), n1, mean(set2), var(set2), n2; alpha=alpha,  method=method)
            if ci.lower > ref pow += 1 end
        end
        return pow/nsim
    end



end #end module
