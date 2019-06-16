# Clinical Trial Utilities
# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)

module SIM
    using Distributions
    using Random
    import ..ctSampleN
    import ..designProp
    import ..cv2ms
    import ..ZDIST
    import ..CTUException
    import ..CI.twoProp
    import ..CI.twoMeans

    export bePower, ctPropPower, ctPropSampleN, ctMeansPower, ctMeansPowerFS

    function bePower(;alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.0, n=0, simnum=5, seed=0)
        if alpha <= 0.0 || alpha >= 1.0  throw(CTUException(1111,"SIM.bePower: alpha should be > 0 and < 1")) end
        if cv <= 0.0   throw(CTUException(1112,"SIM.bePower: cv should be > 0")) end
        if n <= 0   throw(CTUException(1113,"SIM.bePower: n should be > 0")) end
        if simnum <= 0   throw(CTUException(1114,"SIM.bePower: simnum should be > 0")) end
        dffunc, bkni, seq = designProp(:d2x2)                                   #dffunc if generic funtion with 1 arg return df
        df    = dffunc(n)
        sqa   = Array{Float64, 1}(undef, seq)
        sqa  .= n÷seq
        for i = 1:n%seq
            sqa[i] += 1
        end
        sef   = sqrt(sum(1 ./ sqa)*bkni)                                        #for se
        #ms    = log(1+cv^2)
        if logscale
            ltheta1 = log(theta1); ltheta2 = log(theta2); diffm = log(theta0); ms = log(1+cv^2);
        else
            ltheta1 = theta1; ltheta2 = theta2; diffm = theta0; ms = cv*cv;
        end
        return bePowerSIM(ltheta1, ltheta2, ms, diffm, df, sef, 10^simnum, alpha, seed=seed)
    end

    function bePowerSIM(theta1, theta2, ms, mean, df, sef, nsim, alpha; seed=0)
        rng = MersenneTwister(1234)
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

    function beSampleSetGen()
    end

    function bePowerFSS()
    end

    function twoStageBEPower()
    end

    function beReplPowerFSS()
    end

    function ctPropPower(p1, n1, p2, n2, ref; alpha=0.05, type=:notdef, citype=:notdef, method=:notdef, simnum=5, seed=0)
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
            ci = twoProp(x1, n1, x2, n2; alpha=alpha, type=citype, method=method)
            if type == :ns
                if ci.lower > ref pow += 1 end
            elseif type == :ei
                if ci.lower > ref[1] &&  ci.upper < ref[2] pow += 1 end
            end
        end
        return pow/nsim
    end

    function ctPropSampleN(p1, p2, ref; alpha=0.05, beta=0.2, type=:notdef, citype =:notdef, method=:notdef, simnum=5, seed=0)
        if type == :notdef || citype == :notdef || method == :notdef throw(CTUException(1116,"ctPropSampleN: type or method not defined.")) end
            st::Int = sn::Int = 10
        if citype == :diff
            st = sn = ceil(ctSampleN(param=:prop, type=type, group=:two, alpha=alpha, beta=beta, diff=ref, a=p1, b=p2))
        elseif citype == :or
            st = sn = ceil(ctSampleN(param=:or, type=type, group=:two, alpha=alpha, beta=beta, diff=ref, a=p1, b=p2, logdiff = false))
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
            ci = twoMeans(sm1, ss1, n1, sm2, ss2, n2; alpha=alpha,  method=:ev)
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
            ci     = twoMeans(mean(set1), var(set1), n1, mean(set2), var(set2), n2; alpha=alpha,  method=method)
            if ci.lower > ref pow += 1 end
        end
        return pow/nsim
    end



end #end module
