# Clinical Trial Utilities
# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)

module SIM
    using Distributions
    import ClinicalTrialUtilities.designProp
    import ClinicalTrialUtilities.ZDIST
    import ClinicalTrialUtilities.CTUException


    function bePower(;alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.0, n=0)
        dffunc, bkni, seq = designProp("2x2")                                   #dffunc if generic funtion with 1 arg return df
        df    = dffunc(n)
        sqa   = Array{Float64, 1}(undef, seq)
        sqa  .= n÷seq
        for i = 1:n%seq
            sqa[i] += 1
        end
        sef   = sqrt(sum(1 ./ sqa)*bkni)                                        #for se
        ms    = log(1+cv^2)                                                     #var

        return bePowerSIM(log(0.8), log(1.25), ms, log(1.0), df, sef, 10^6, 0.05)
    end

    function bePowerSIM(theta1, theta2, ms, mean, df, sef, nsim, alpha)
        CHSQ    = Chisq(df)
        tval    = quantile(TDist(df), 1-alpha)
        se      = sef*sqrt(ms)                                                #se_diff
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

end #end module
