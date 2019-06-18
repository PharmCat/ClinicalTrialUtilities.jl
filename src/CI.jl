# Clinical Trial Utilities
# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)
# Some ideas was taken from R project packages:
# PropCIs by Ralph Scherer https://cran.r-project.org/web/packages/PropCIs/index.html
# pairwiseCI by Frank Schaarschmidt, Daniel Gerhard  https://CRAN.R-project.org/package=pairwiseCI
# binGroup by Boan Zhang, Christopher Bilder, Brad Biggerstaff, Frank Schaarschmidt Brianna Hitt https://CRAN.R-project.org/package=binGroup
# proportion by M.Subbiah, V.Rajeswaran https://CRAN.R-project.org/package=proportion
# binom by Sundar Dorai-Raj https://CRAN.R-project.org/package=binom
# DescTools https://CRAN.R-project.org/package=DescTools

module CI
    using Distributions, Roots, DataFrames
    import ..ZDIST
    import ..CTUException
    import ..ConfInt
    #const ZDIST = Normal()
    export oneProp, oneMeans, twoProp, twoMeans, cmh

    function oneProp(x::Int, n::Int; alpha=0.05, method=:wilson)
        if method==:wilson
            return propWilsonCI(x, n, alpha)
        elseif method==:wilsoncc
            return propWilsonCCCI(x, n, alpha)
        elseif method==:cp
            return propCPCI(x, n, alpha)
        elseif method==:soc
            return propSOCCI(x, n, alpha)
        elseif method==:blaker
            return propBlakerCI(x, n, alpha)
        elseif method==:arcsine
            return propARCCI(x, n, alpha)
        elseif method==:wald
            return propWaldCI(x, n, alpha)
        else
            throw(CTUException(1301,"oneProp: no such method."))
        end
    end

    function oneMean(m,s,n,alpha; method=:norm)
        if method==:norm
            meanNormCI(m,s,n,alpha)
        elseif method==:tdist
            meanTdistCI(m,s,n,alpha)
        end
    end

    function twoProp(x1::Int, n1::Int, x2::Int, n2::Int; alpha=0.05, type::Symbol, method::Symbol)::ConfInt

        if type==:diff
            if method ==:nhs
                return propDiffNHSCI(x1, n1, x2, n2, alpha)
            elseif method ==:nhscc
                return propDiffNHSCCCI(x1, n1, x2, n2, alpha)
            elseif method ==:ac
                return propDiffACCI(x1, n1, x2, n2, alpha)
            elseif method ==:mn
                return propDiffMNCI(x1, n1, x2, n2, alpha)
            elseif method ==:mee
                return propDiffMeeCI(x1, n1, x2, n2, alpha)
            elseif method ==:mee2
                return propDiffFMCI(x1, n1, x2, n2, alpha)
            elseif method ==:wald
                return propDiffWaldCI(x1, n1, x2, n2, alpha)
            elseif method ==:waldcc
                return propDiffWaldCCCI(x1, n1, x2, n2, alpha)
            end
        elseif type==:rr
            if method ==:cli
                return propRRCLICI(x1, n1, x2, n2, alpha)
            elseif method ==:mover
                return  propRRMOVERCI(x1, n1, x2, n2, alpha)
            end
        elseif type==:or
            if method==:mn
                return propORCI(x1, n1, x2, n2, alpha)
            elseif method==:awoolf
                return propORaWoolfCI(x1, n1, x2, n2, alpha)
            elseif method==:woolf
                return propORWoolfCI(x1, n1, x2, n2, alpha)
            end
        end
    end #twoProp

    function twoMeans(m1::Real, s1::Real, n1::Real, m2::Real, s2::Real, n2::Real; alpha::Real=0.05, type=:diff, method=:notdef)::ConfInt
        if type==:diff
            if method == :ev
                return meanDiffEV(m1::Real, s1::Real, n1::Real, m2::Real, s2::Real, n2::Real, alpha::Real)
            elseif method == :uv
                return meanDiffUV(m1::Real, s1::Real, n1::Real, m2::Real, s2::Real, n2::Real, alpha::Real)
            end
        end
    end #twoMeans

    #-----------------------------PROPORTIONS-----------------------------------

    #Wilson’s confidence interval for a single proportion, wilson score
    #Wilson, E.B. (1927) Probable inference, the law of succession, and statistical inferenceJ. Amer.Stat. Assoc22, 209–212
    function propWilsonCI(x::Int, n::Int, alpha::Float64)::ConfInt
        z = abs(quantile(ZDIST, 1-alpha/2))
        p = x/n
        b = z*sqrt((p*(1-p)+(z^2)/(4*n))/n)/(1+(z^2)/n)
        m = (p+(z^2)/(2*n))/(1+(z^2)/n)
        return ConfInt(m - b,m + b,m)
    end
    #Wilson CC
    #Newcombe, R. G. (1998). "Two-sided confidence intervals for the single proportion: comparison of seven methods". Statistics in Medicine. 17 (8): 857–872. doi:10.1002/(SICI)1097-0258(19980430)17:8<857::AID-SIM777>3.0.CO;2-E. PMID 9595616
    function propWilsonCCCI(x::Int, n::Int, alpha::Float64)::ConfInt
        z = abs(quantile(ZDIST, 1-alpha/2))
        p = x/n
        l = (2*n*p+z*z-1-z*sqrt(z*z-2-1/n+4*p*(n*(1-p)+1)))/2/(n+z*z)
        u = (2*n*p+z*z+1+z*sqrt(z*z+2-1/n+4*p*(n*(1-p)-1)))/2/(n+z*z)
        return ConfInt(min(p, l), max(p, u), p)
    end

    #Clopper-Pearson exatct CI
    #Clopper, C. and Pearson, E.S. (1934) The use of confidence or fiducial limits illustrated in the caseof the binomial.Biometrika26, 404–413.
    function propCPCI(x::Int, n::Int, alpha::Float64)::ConfInt
        if x==0
            ll = 0.0
            ul = 1.0-(alpha/2)^(1/n)
        elseif x==n
            ul = 1.0
            ll = (alpha/2)^(1/n)
        else
            ll = 1/(1+(n-x+1)/(x*quantile(FDist(2*x, 2*(n-x+1)), alpha/2)))
            ul = 1/(1+(n-x) / ((x+1)*quantile(FDist(2*(x+1), 2*(n-x)), 1-alpha/2)))
        end
        return ConfInt(ll, ul, x/n)
    end
    #Blaker CI
    #Blaker, H. (2000). Confidence curves and improved exact confidence intervals for discrete distributions,Canadian Journal of Statistics28 (4), 783–798
    function propBlakerCI(x::Int, n::Int, alpha::Float64)::ConfInt
        tol = 1E-5; lower = 0; upper = 1;
        if n != 0
            lower = quantile(Beta(x, n-x+1), alpha/2)
            while acceptbin(x, n, lower+tol) < alpha
                lower +=tol
            end
        end
        if x != n
            upper = quantile(Beta(x+1, n-x), 1-alpha/2)
            while acceptbin(x, n, upper-tol) < alpha
                upper -=tol
            end
        end
        return ConfInt(lower,upper, x/n)
    end
    @inline function acceptbin(x::Int, n::Int, p::Float64)
        BIN = Binomial(n,p)
        p1 = 1-cdf(BIN,x-1)
        p2 =   cdf(BIN,x)
        a1 = p1 + cdf(BIN, quantile(BIN,p1)-1)
        a2 = p2+1-cdf(BIN, quantile(BIN,1-p2))
        return min(a1,a2)
    end
    #Wald CI
    function propWaldCI(x::Int, n::Int, alpha::Float64)::ConfInt
        p=x/n
        b = quantile(ZDIST, 1-alpha/2)*sqrt(p*(1-p)/n)
        return ConfInt(p-b,p+b,p)
    end
    #SOC  Second-Order corrected
    #T. Tony Cai One-sided con&dence intervals in discrete distributions doi:10.1016/j.jspi.2004.01.00
    function propSOCCI(x::Int, n::Int, alpha::Float64)::ConfInt
        p  = x/n
        k  = quantile(ZDIST, 1-alpha/2)
        k2 = k^2
        η  = k2/3+1/6
        γ1 = -(k2*13/18+17/18)
        γ2 = k2/18+7/36
        m  = (x+η)/(n+2*η)
        b  = k*sqrt(p*(1-p)+(γ1*p*(1-p)+γ2)/n)/sqrt(n)
        return ConfInt(m-b, m+b, p)
    end
    #Arcsine
    function propARCCI(x::Int, n::Int, alpha::Float64)::ConfInt
        q = quantile(ZDIST, 1-alpha/2)
        p = x/n
        z = q/(2*sqrt(n))
        return ConfInt(sin(asin(sqrt(p))-z)^2, sin(asin(sqrt(p))+z)^2, p)
    end
    #--------------------------------OR-----------------------------------------
    #Cornfield, J. (1956) A statistical problem arising from retrospective studies.  In Neyman J. (ed.),Proceedings of the third Berkeley Symposium on Mathematical Statistics and Probability4,  pp.135–148.
    #Miettinen O. S., Nurminen M. (1985) Comparative analysis of two rates.Statistics in Medicine4,213–226
    #Agresti, A. 2002. Categorical Data Analysis. Wiley, 2nd Edition.
    #MN Score
    function propORCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
        p1 = x1/n1
        p2 = x2/n2
        if (x1==0 && x2==0) || (x1==n1 && x2==n2)
            ul    = Inf
            ll    = 0.0
        elseif x1==0 || x2==n2
            ll    = 0.0
            theta = 0.01/n2
            ul    = limit(x1,n1,x2,n2,alpha, theta, 1)
        elseif x1==n1 || x2 == 0
            ul    = Inf
            theta = 100*n1
            ll    = limit(x1,n1,x2,n2,alpha, theta, 0)
        else
            theta = p1/(1-p1)/(p2/(1-p2))/1.1
            ll    = limit(x1,n1,x2,n2,alpha,theta,0)
            theta = p1/(1-p1)/(p2/(1-p2))*1.1
            ul    = limit(x1,n1,x2,n2,alpha,theta,1)
        end
        return ConfInt(ll, ul, (p1/(1-p1))/(p2/(1-p2)))
    end
    @inline function limit(x1, n1, x2, n2, alpha, lim, t)
        z  = quantile(Chisq(1), 1-alpha)
        ci::Float64 = 0.0
        px = x1/n1
        score = 0.0
        while score < z
            a = n2*(lim-1)
            b = n1*lim+n2-(x1+x2)*(lim-1)
            c = -(x1+x2)
            p2d = (-b+sqrt(b^2-4*a*c))/(2*a)
            p1d = p2d*lim/(1+p2d*(lim-1))
            score = ((n1*(px-p1d))^2)*(1/(n1*p1d*(1-p1d))+1/(n2*p2d*(1-p2d)))*(n1+n2-1)/(n1+n2)
            ci = lim
            if t==0 lim = ci/1.001 else lim = ci*1.001 end
        end
        return ci
    end #limit

    #Adjusted Woolf interval (Gart adjusted logit) Lawson, R (2005):Smallsample confidence intervals for the odds ratio.  Communication in Statistics Simulation andComputation, 33, 1095-1113.
    function propORaWoolfCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
            xa = x1 + 0.5
            xb = n1 - x1 + 0.5
            xc = x2 + 0.5
            xd = n2 - x2 + 0.5
            est = xa*xd/xc/xb
            estI = log(est)
            stde = sqrt(1/xa + 1/xb + 1/xc + 1/xd)
            z   = quantile(ZDIST, 1-alpha/2)
            return ConfInt(exp(estI - z*stde), exp(estI + z*stde), est)
    end
    #Woolf logit
    #Woolf, B. (1955). On estimating the relation between blood group and disease. Annals of human genetics, 19(4):251-253.
    function propORWoolfCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
            xa = x1
            xb = n1 - x1
            xc = x2
            xd = n2 - x2
            estimate = xa*xd/xc/xb
            estI = log(estimate)
            stde = sqrt(1/xa + 1/xb + 1/xc + 1/xd)
            z   = quantile(ZDIST, 1-alpha/2)
            return ConfInt(exp(estI - z*stde), exp(estI + z*stde), estimate)
    end

    #Method of variance estimates recovery
    #Donner, A. and Zou, G. (2012). Closed-form confidence intervals for functions of the normal mean and standard deviation. Statistical Methods in Medical Research, 21(4):347-359.
    #function propRRMOVERCI()
    #end

    #?Agresti independence-smoothed logit

    #?Cornfield exact confidence interval
    #https://rdrr.io/cran/ORCI/man/Cornfieldexact.CI.html



    #------------------------------DIFF-----------------------------------------
    #Wald
    function propDiffWaldCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
        p1  = x1/n1
        p2  = x2/n2
        est = p1-p2
        q   = quantile(ZDIST, 1 - alpha/2)
        stderr = sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)
        return ConfInt(est - q*stderr, est + q*stderr, est)
    end
    #Wald CC
    function propDiffWaldCCCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
        p1  = x1/n1
        p2  = x2/n2
        est = p1-p2
        cc  = 0.5*(1/n1+1/n2)
        q   = quantile(ZDIST, 1 - alpha/2)
        stderr = sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)
        return ConfInt(est - q*stderr - cc, est + q*stderr + cc, est)
    end

    #Newcombes Hybrid (wilson) Score interval for the difference of proportions
    #Newcombe 1998
    #10
    #Newcombe RG (1998): Interval Estimation for the Difference Between Independent Proportions: Comparison of Eleven Methods. Statistics in Medicine 17, 873-890.
    function propDiffNHSCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
        p1  = x1/n1
        p2  = x2/n2
        est = p1-p2
        q   = quantile(ZDIST, 1 - alpha/2)
        ci1 = propWilsonCI(x1, n1, alpha)
        ci2 = propWilsonCI(x2, n2, alpha)
        return ConfInt(est-q*sqrt(ci1.lower*(1-ci1.lower)/n1+ci2.upper*(1-ci2.upper)/n2),
                       est+q*sqrt(ci1.upper*(1-ci1.upper)/n1+ci2.lower*(1-ci2.lower)/n2), est)
    end

    function propDiffNHSCCCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
        p1  = x1/n1
        p2  = x2/n2
        est = p1-p2
        q   = quantile(ZDIST, 1 - alpha/2)
        ci1 = propWilsonCCCI(x1, n1, alpha)
        ci2 = propWilsonCCCI(x2, n2, alpha)
        return ConfInt(est-sqrt((p1-ci1.lower)^2 + (ci2.upper-p2)^2),
                       est+sqrt((ci1.upper - p1)^2 + (p2 - ci2.lower)^2), est)
    end

    #Agresti-Caffo interval for the difference of proportions
    #13
    #Agresti A, Caffo B., “Simple and effective confidence intervals for proportions and differences of proportions result from adding two successes and two failures”, American Statistician 54: 280–288 (2000)
    function propDiffACCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
        p1       = x1/n1
        p2       = x2/n2
        est      = p1-p2
        q        = quantile(ZDIST, 1 - alpha/2)
        p1I      = (x1+1)/(n1+2)
        p2I      = (x2+1)/(n2+2)
        n1I      = n1+2
        n2I      = n2+2
        estI     = p1I-p2I
        stderr   = sqrt(p1I*(1-p1I)/n1I+p2I*(1-p2I)/n2I)
        return ConfInt(estI-q*stderr, estI+q*stderr, est)
    end
    #Method of Mee 1984 with Miettinen and Nurminen modification n / (n - 1) Newcombe 1998
    #Score intervals for the difference of two binomial proportions
    #6
    function propDiffMNCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt  #With correction factor
        p1    = x1/n1
        p2    = x2/n2
        z     = quantile(Chisq(1), 1-alpha)
        score = 0
        proot = p1 - p2
        dp    = 1 - proot
        up2   = 0
        while dp > 0.0000001 && abs(z-score) > 0.000001
            dp    = dp/2
            up2   = proot + dp
            score = mnzstat(p1, n1, p2, n2, up2)
            if score < z proot = up2 end
        end
        score = 0
        proot = p1-p2
        dp    = 1 + proot
        low2  = 0
        while dp > 0.0000001 && abs(z-score) > 0.000001
            dp    = dp/2
            low2  = proot - dp
            score = mnzstat(p1, n1, p2, n2, low2)
            if score < z proot = low2 end
        end
        return ConfInt(low2, up2, p1-p2)
    end

    @inline function mnzstat(p1::Float64, n1::Int, p2::Float64, n2::Int, delta::Float64)::Float64
        diff = p1-p2-delta
        n    = n1+n2
        return diff*diff/(n/(n-1)*fmvar(p1, n1, p2, n2, delta))
    end
    #Mee
    #Mee RW (1984) Confidence bounds for the difference between two probabilities,Biometrics40:1175-1176
    function propDiffMeeCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)
        p1   = x1/n1
        p2   = x2/n2
        est  = p1 - p2
        f(x) = fmpval(p1, n1, p2, n2, est, x) - alpha
        return ConfInt(find_zero(f, (-1.0+1e-6, est-1e-6)), find_zero(f, ( est+1e-6, 1.0-1e-6)), est)

    end
    @inline function fmvar(p1::Float64, n1::Int, p2::Float64, n2::Int, delta::Float64)::Float64
        if p1-p2-delta == 0 return 0.0 end
        theta = n2/n1
        a = 1 + theta
        b = -(1 + theta + p1 + theta * p2 + delta*(theta + 2))
        c = delta^2 + delta*(2 * p1 + theta + 1) + p1 + theta * p2
        d = -p1*delta*(1 + delta)
        v = (b/a/3)^3 - b*c/(6*a*a) + d/2/a
        u = sign(v) * sqrt((b/3/a)^2 - c/3/a)
        w = (pi + acos(v/u^3))/3
        p1n = 2*u*cos(w) - b/3/a
        p2n = p1n - delta
        return p1n * (1-p1n)/n1 + p2n * (1 - p2n)/n2
    end
    @inline function fmpval(p1, n1, p2, n2, est, delta)::Float64
        z = (est-delta)/sqrt(fmvar(p1, n1, p2, n2, delta))
        p = cdf(ZDIST, z)
        return 2*min(1-p, p)
    end
    #FM2 - Faster than propDiffMeeCI
    #Farrington, C. P. and Manning, G. (1990), “Test Statistics and Sample Size Formulae for Comparative Binomial Trials with Null Hypothesis of Non-zero Risk Difference or Non-unity Relative Risk,” Statistics in Medicine, 9, 1447–1454
    function propDiffFMCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
        p1  = x1/n1
        p2  = x2/n2
        est = p1 - p2
        z   = quantile(ZDIST, 1 - alpha/2)
        f(x) = fmpval2(p1, n1, p2, n2, est, x) - z
        return ConfInt(find_zero(f, (-1.0+1e-6, est-1e-6), atol=1E-5), find_zero(f, (est+1e-6, 1.0-1e-6), atol=1E-5), est)
    end
    @inline function fmpval2(p1, n1, p2, n2, est, delta)
        return abs((est-delta)/sqrt(fmvar(p1, n1, p2, n2, delta)))
    end
    #--------------------------------RR-----------------------------------------
    #Miettinen-Nurminen Score interval
    #Not implemented
    #function propRRMNCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64) end

    #Crude log interval
    #Gart, JJand Nam, J (1988): Approximate interval estimation of the ratio of binomial parameters: Areview and corrections for skewness. Biometrics 44, 323-338.
    function propRRCLICI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
        x1I  = x1+0.5
        x2I  = x2+0.5
        n1I  = n1+0.5
        n2I  = n2+0.5
        estI = log((x1I/n1I)/(x2I/n2I))
        stderrlog = sqrt(1/x2I+1/x1I-1/n2I-1/n1I)
        estimate  = (x1/n1)/(x2/n2)
        Z         =  quantile(ZDIST,1-alpha/2)
        return ConfInt(exp(estI-Z*stderrlog), exp(estI+Z*stderrlog), estimate)
    end
    #Method of variance estimates recovery (Donner, Zou, 2012)
    function propRRMOVERCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
        p1       = (x1/n1)
        p2       = (x2/n2)
        estimate = (x1/n1)/(x2/n2)
        Z        = quantile(ZDIST, 1-alpha/2)
        wilci1   = propWilsonCI(x1, n1, alpha)
        wilci2   = propWilsonCI(x2, n2, alpha)
        lower    = (p1*p2-sqrt((p1*p2)^2 - wilci1.lower*wilci2.upper*(2*p1-wilci1.lower)*(2*p2-wilci2.upper)))/(wilci2.upper*(2*p2 - wilci2.upper))
        upper    = (p1*p2+sqrt((p1*p2)^2 - wilci1.upper*wilci2.lower*(2*p1-wilci1.upper)*(2*p2-wilci2.lower)))/(wilci2.lower*(2*p2 - wilci2.lower))
        return ConfInt(lower, upper, estimate)
    end
    #-------------------------------MEANS---------------------------------------

    #Normal
    function meanNormCI(m,s,n,alpha)::ConfInt
        e = quantile(ZDIST, 1-alpha/2)*sqrt(s/n)
        return ConfInt(m-e, m+e, m)
    end
    #T Distribution
    function meanTdistCI(m,s,n,alpha)::ConfInt
        e = quantile(TDist(n-1), 1-alpha/2)*sqrt(s/n)
        return ConfInt(m-e, m+e, m)
    end
    #mean diff equal var
    function meanDiffEV(m1::Real, s1::Real, n1::Real, m2::Real, s2::Real, n2::Real, alpha::Real)::ConfInt
        diff   = m1 - m2
        stddev = sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
        stderr = stddev * sqrt(1/n1 + 1/n2)
        d      = stderr*quantile(TDist(n1+n2-2), 1-alpha/2)
        return ConfInt(diff-d, diff+d, diff)
    end
    function meanDiffEV(a1::AbstractVector{T}, a2::AbstractVector{S}, alpha::Real)::ConfInt where {T<:Real,S<:Real}
        return meanDiffEV(mean(a1), var(a1), length(a1), mean(a2), var(a2), length(a2), alpha)
    end
    #mean diff unequal var
    #Two sample t-test (unequal variance)
    #Welch-Satterthwaite df
    function meanDiffUV(m1::Real, s1::Real, n1::Real, m2::Real, s2::Real, n2::Real, alpha::Real)::ConfInt
        diff   = m1 - m2
        v      = (s1/n1+s2/n2)^2/(s1^2/n1^2/(n1-1)+s2^2/n2^2/(n2-1))
        stderr = sqrt(s1/n1 + s2/n2)
        d      = stderr*quantile(TDist(v), 1-alpha/2)
        return ConfInt(diff-d, diff+d, diff)
    end
    function meanDiffUV(a1::AbstractVector{T}, a2::AbstractVector{S}, alpha::Real)::ConfInt where {T<:Real,S<:Real}
        return meanDiffUV(mean(a1), var(a1), length(a1), mean(a2), var(a2), length(a2), alpha)
    end

    #Not validated
    function bartlettsTest(s1::Real, n1::Real, s2::Real, n2::Real)
        n = n1 + n2
        k = 2
        s = ((n1 - 1)*s1 + (n2 - 1)*s2)/(n - k)
        x = ((n-k)*log(s) - ((n1-1)*log(s1)+(n2-1)*log(s2)))/(1+(1/(n1-1) + 1/(n2-1) - 1/(n-k))/3/(k-1))
        return 1-cdf(Chisq(k-1),x), x
    end

    #Cochran–Mantel–Haenszel confidence intervals
    #p275, Rothman, K. J., Greenland, S., & Lash, T. L. (2008). Modern epidemiology (3rd ed.). Philadelphia: Lippincott Williams & Wilkins.
    #1 equation in: Sato, Greenland, & Robins (1989)
    #2 equation in: Greenland & Robins (1985)
    # metafor: Meta-Analysis Package for R - Wolfgang Viechtbauer
    function cmh(data::DataFrame; a = :a, b = :b, c = :c, d = :d, alpha = 0.05, type = :diff, method = :default)::ConfInt
        if type == :diff
            if method == :default method = :sato end #default method

            n1 = data[a] + data[b]
            n2 = data[c] + data[d]
            N  = data[a] + data[b] + data[c] + data[d]
            est = sum(data[a] .* (n2 ./ N) - data[c] .* (n1 ./ N)) / sum(n1 .* (n2 ./ N))
            #1
            se   = sqrt((est * (sum(data[c] .* (n1 ./ N) .^ 2 - data[a] .* (n2 ./ N) .^2 + (n1 ./ N) .* (n2 ./ N) .* (n2-n1) ./ 2)) + sum(data[a] .* (n2 - data[c]) ./ N + data[c] .* (n1 - data[a]) ./ N)/2) / sum(n1 .* (n2 ./ N)) .^ 2) # equation in: Sato, Greenland, & Robins (1989)
            #2
            #se   = sqrt(sum(((data[a]./N.^2).*data[b].*(n2.^2./n1) + (data[c]./N.^2).*data[d].*(n1.^2./n2))) / sum(n1.*(n2./N))^2)
            #zval = est/se
            #pval = 2*(1-cdf(Normal(), abs(zval)))
            z = quantile(ZDIST, 1 - alpha/2)
            return ConfInt(est - z*se, est + z*se, est)
        elseif type == :or
            #...
        elseif type == :rr
            #...
        end
    end

end #end module CI
