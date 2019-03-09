# Clinical Trial Utilities
# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)

module SIM
    using Distributions, Random
    import ..designProp
    import ..cv2ms
    import ..ZDIST
    import ..CTUException
    import ..CI.twoProp
    import ..CI.twoMeans

    export bePower, ctBinPower


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
        for i=1:nsim
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

    function ctPropPower(p1, n1, p2, n2, diff; alpha=0.05, type=:or, method=:mn, simnum=5, seed=0)
        rng = MersenneTwister(1234)
        if seed == 0  Random.seed!(rng) else Random.seed!(seed) end

        BIN1 = Binomial(n1, p1)
        BIN2 = Binomial(n2, p2)
        pow     = 0
        nsim = 10^simnum

        for i=1:nsim
            x1 = rand(BIN1)
            x2 = rand(BIN2)
            ci = twoProp(x1, n1, x2, n2; alpha=alpha, type=type, method=method)
            if ci.lower > diff pow += 1 end
        end
        return pow/nsim
    end

    function ctMeansPower()

    end



end #end module
