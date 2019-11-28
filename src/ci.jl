# Clinical Trial Utilities
# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)
# Some ideas was taken from R project packages:
# PropCIs by Ralph Scherer https://cran.r-project.org/web/packages/PropCIs/index.html
# pairwiseCI by Frank Schaarschmidt, Daniel Gerhard  https://CRAN.R-project.org/package=pairwiseCI
# binGroup by Boan Zhang, Christopher Bilder, Brad Biggerstaff, Frank Schaarschmidt Brianna Hitt https://CRAN.R-project.org/package=binGroup
# proportion by M.Subbiah, V.Rajeswaran https://CRAN.R-project.org/package=proportion
# binom by Sundar Dorai-Raj https://CRAN.R-project.org/package=binom
# DescTools https://CRAN.R-project.org/package=DescTools
# metafor by Wolfgang Viechtbauer https://cran.r-project.org/package=metafor

"""
    StatsBase.confint(param::Proportion; level = 0.95, method = :default)::ConfInt

# Computation methods:

- wilson - Wilson's confidence interval (CI) for a single proportion (wilson score) (default);
- wilsoncc - Wilson's CI with continuity correction (CC);
- cp - Clopper-Pearson exact CI;
- soc - SOC: Second-Order corrected CI;
- blaker - Blaker exact CI for discrete distributions;
- arcsine - Arcsine CI;
- wald - Wald CI without CC;
"""
function StatsBase.confint(param::Proportion; level = 0.95, method = :default)::ConfInt
    propci(param.x, param.n; alpha = 1 - level, method = method)
end
function StatsBase.confint(param::DiffProportion{P, P}; level = 0.95, method = :default)::ConfInt  where P <: Proportion
    diffpropci(param.a.x, param.a.n, param.b.x, param.b.n; alpha = 1 - level, method = method)
end
function StatsBase.confint(param::OddRatio{P}; level = 0.95, method = :default)::ConfInt where P <: Proportion
    orpropci(param.a.x, param.a.n, param.b.x, param.b.n; alpha = 1 - level, method = method)
end
function StatsBase.confint(param::RiskRatio{P}; level = 0.95, method = :default)::ConfInt where P <: Proportion
    rrpropci(param.a.x, param.a.n, param.b.x, param.b.n; alpha = 1 - level, method = method)
end
function StatsBase.confint(param::Mean; level = 0.95, method = :default)::ConfInt
    meanci(param.m, param.s, param.n; alpha = 1 - level, method = method)
end
function StatsBase.confint(param::DiffMean{T}; level = 0.95, method = :default)::ConfInt where T <: Mean
    diffmeanci(param.a.m, param.a.s, param.a.n, param.b.m, param.b.s, param.b.n; alpha = 1 - level, method = method)
end
"""
    propci(x::Int, n::Int; alpha=0.05, method = :default)::ConfInt

Confidence interval for proportion.

# Computation methods:

- :wilson | :default - Wilson's confidence interval (CI) for a single proportion (wilson score);
- :wilsoncc - Wilson's CI with continuity correction (CC);
- :cp - Clopper-Pearson exact CI;
- :soc - SOC: Second-Order corrected CI;
- :blaker - Blaker exact CI for discrete distributions;
- :arcsine - Arcsine CI;
- :wald - Wald CI without CC;
- :waldcc - Wald CI with CC;
"""
function propci(x::Int, n::Int; alpha::Real = 0.05, method = :default)::ConfInt
    if alpha >= 1.0 || alpha <= 0.0
        @warn "Alpha >= 1.0 or <= 0.0"
        return ConfInt(0,0,x/n)
    end
    if method == :wilson || method == :default
        return propwilsonci(x, n, alpha)
    elseif method==:wilsoncc
        return propwilsonccci(x, n, alpha)
    elseif method==:cp
        return propcpci(x, n, alpha)
    elseif method==:soc
        return propsocci(x, n, alpha)
    elseif method==:blaker
        return propblakerci(x, n, alpha)
    elseif method==:arcsine
        return proparcci(x, n, alpha)
    elseif method==:wald
        return propwaldci(x, n, alpha)
    elseif method==:waldcc
        return propwaldcicc(x, n, alpha)
    else
        throw(ArgumentError("unknown method!"))
    end
end
"""
    diffpropci(x1::Int, n1::Int, x2::Int, n2::Int; alpha::Real = 0.05, method::Symbol = :default)::ConfInt

Confidence interval for proportion difference.

# Computation methods:

- :nhs - Newcombes Hybrid (wilson) Score interval for the difference of proportions;
- :nhscc - Newcombes Hybrid Score CC;
- :ac - Agresti-Caffo interval for the difference of proportions;
- :mn | :default - Method of Mee 1984 with Miettinen and Nurminen modification;
- :mee | :fm - Mee maximum likelihood method;
- :wald - Wald CI without CC;
- :waldcc - Wald CI with CC;
"""
function diffpropci(x1::Int, n1::Int, x2::Int, n2::Int; alpha::Real = 0.05, method::Symbol = :default)::ConfInt
    if alpha >= 1.0 || alpha <= 0.0
        throw(ArgumentError("Alpha shold be > 0.0 and < 1.0"))
    end
    if method == :nhs
        return propDiffNHSCI(x1, n1, x2, n2, alpha)
    elseif method == :nhscc
        return propDiffNHSCCCI(x1, n1, x2, n2, alpha)
    elseif method == :ac
        return propDiffACCI(x1, n1, x2, n2, alpha)
    elseif method == :mn || method == :default
        return propDiffMNCI(x1, n1, x2, n2, alpha)
    elseif method == :mee2
        return propDiffMeeCI(x1, n1, x2, n2, alpha)
    elseif method == :mee || method == :fm
        return propDiffFMCI(x1, n1, x2, n2, alpha)
    elseif method == :wald
        return propDiffWaldCI(x1, n1, x2, n2, alpha)
    elseif method == :waldcc
        return propDiffWaldCCCI(x1, n1, x2, n2, alpha)
    else
        throw(ArgumentError("Method unknown!"))
    end
end
"""
    rrpropci(x1::Int, n1::Int, x2::Int, n2::Int; alpha::Real = 0.05, method::Symbol = :default)::ConfInt

Confidence interval for relative risk.

# Computation methods:

- :mn | :default - Miettinen-Nurminen Score interval;
- :cli | :walters - Crude log interval;
- :li | :katz - Log interval for the risk ratio;
- :mover - Method of variance estimates recovery;
"""
function rrpropci(x1::Int, n1::Int, x2::Int, n2::Int; alpha::Real = 0.05, method::Symbol = :default)::ConfInt
    if alpha >= 1.0 || alpha <= 0.0 throw(ArgumentError("Alpha shold be > 0.0 and < 1.0")) end
    if method==:mn || method == :default
        return propRRMNCI(x1, n1, x2, n2, alpha)
    elseif method == :cli || method == :walters
        return propRRCLICI(x1, n1, x2, n2, alpha)
    elseif method == :li || method == :katz
        return proprrkatzci(x1, n1, x2, n2, alpha)
    elseif method ==:mover
        return  propRRMOVERCI(x1, n1, x2, n2, alpha)
    else
        throw(ArgumentError("Method unknown!"))
    end
end
"""
    orpropci(x1::Int, n1::Int, x2::Int, n2::Int; alpha::Real = 0.05, method::Symbol = :default)::ConfInt

Confidence interval for odd ratio.

# Computation methods:

- :mn - Miettinen-Nurminen CI (deprecated);
- :mn2 | :default - Miettinen-Nurminen CI;
- :woolf - Woolf logit CI;
- :awoolf | :gart - Adjusted Woolf interval (Gart adjusted logit);
- :mover - Method of variance estimates recovery;
"""
function orpropci(x1::Int, n1::Int, x2::Int, n2::Int; alpha::Real = 0.05, method::Symbol = :default)::ConfInt
    if alpha >= 1.0 || alpha <= 0.0 throw(ArgumentError("Alpha shold be > 0.0 and < 1.0")) end
    if method==:mn
        return propORMNCI(x1, n1, x2, n2, alpha)
    elseif method==:awoolf || method==:gart
        return propORaWoolfCI(x1, n1, x2, n2, alpha)
    elseif method==:woolf
        return propORWoolfCI(x1, n1, x2, n2, alpha)
    elseif method==:mover
        return propORMOVERCI(x1, n1, x2, n2, alpha)
    elseif method==:mn2 || method == :default
        return propORCI(x1, n1, x2, n2, alpha)
    else
        throw(ArgumentError("Method unknown!"))
    end
end
"""
    meanci(m::Real, s::Real, n::Int; alpha::Real = 0.05, method=:default)::ConfInt

Confidence interval for mean.

# Computation methods:

- :norm - Normal distribution (default);
- :tdist - T Distribution.

"""
function meanci(m::Real, s::Real, n::Int; alpha::Real = 0.05, method=:default)::ConfInt
        if method==:norm || method == :default
            return meanNormCI(m, s, n, alpha)
        elseif method == :tdist
            return meanTdistCI(m, s, n, alpha)
        end
end
"""
    diffmeanci(m1::Real, s1::Real, n1::Real, m2::Real, s2::Real, n2::Real; alpha::Real = 0.05, method::Symbol = :default)::ConfInt

Confidence interval for mead difference.

- :ev - equal variance (default);
- :uv - unequal variance with Welch-Satterthwaite df correction.
"""
function diffmeanci(m1::Real, s1::Real, n1::Real, m2::Real, s2::Real, n2::Real; alpha::Real = 0.05, method::Symbol = :default)::ConfInt
        if method == :ev || method == :default
            return meanDiffEV(m1::Real, s1::Real, n1::Real, m2::Real, s2::Real, n2::Real, alpha::Real)
        elseif method == :uv
            return meanDiffUV(m1::Real, s1::Real, n1::Real, m2::Real, s2::Real, n2::Real, alpha::Real)
        end
end

    #-----------------------------PROPORTIONS-----------------------------------

    #Wilson’s confidence interval for a single proportion, wilson score
    #Wilson, E.B. (1927) Probable inference, the law of succession, and statistical inferenceJ. Amer.Stat. Assoc22, 209–212
    function propwilsonci(x::Int, n::Int, alpha::Float64)::ConfInt
        z = abs(quantile(ZDIST, 1-alpha/2))
        p = x/n
        b = z*sqrt((p*(1-p)+(z^2)/(4*n))/n)/(1+(z^2)/n)
        m = (p+(z^2)/(2*n))/(1+(z^2)/n)
        return ConfInt(m - b, m + b, m, alpha)
    end
    #Wilson CC
    #Newcombe, R. G. (1998). "Two-sided confidence intervals for the single proportion: comparison of seven methods". Statistics in Medicine. 17 (8): 857–872. doi:10.1002/(SICI)1097-0258(19980430)17:8<857::AID-SIM777>3.0.CO;2-E. PMID 9595616
    function propwilsonccci(x::Int, n::Int, alpha::Float64)::ConfInt
        z = abs(quantile(ZDIST, 1-alpha/2))
        p = x/n
        l = (2*n*p+z*z-1-z*sqrt(z*z-2-1/n+4*p*(n*(1-p)+1)))/2/(n+z*z)
        u = (2*n*p+z*z+1+z*sqrt(z*z+2-1/n+4*p*(n*(1-p)-1)))/2/(n+z*z)
        return ConfInt(min(p, l), max(p, u), p, alpha)
    end

    #Clopper-Pearson exatct CI
    #Clopper, C. and Pearson, E.S. (1934) The use of confidence or fiducial limits illustrated in the caseof the binomial.Biometrika26, 404–413.
    function propcpci(x::Int, n::Int, alpha::Float64)::ConfInt
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
        return ConfInt(ll, ul, x/n, alpha)
    end
    #Blaker CI
    #Blaker, H. (2000). Confidence curves and improved exact confidence intervals for discrete distributions,Canadian Journal of Statistics28 (4), 783–798
    function propblakerci(x::Int, n::Int, alpha::Float64)::ConfInt
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
        return ConfInt(lower,upper, x/n, alpha)
    end
    @inline function acceptbin(x::Int, n::Int, p::Float64)::Float64
        BIN = Binomial(n,p)
        p1 = 1-cdf(BIN,x-1)
        p2 =   cdf(BIN,x)
        a1 = p1 + cdf(BIN, quantile(BIN,p1)-1)
        a2 = p2+1-cdf(BIN, quantile(BIN,1-p2))
        return min(a1,a2)
    end
    #Wald CI
    function propwaldci(x::Int, n::Int, alpha::Float64)::ConfInt
        p=x/n
        b = quantile(ZDIST, 1-alpha/2)*sqrt(p*(1-p)/n)
        return ConfInt(p-b, p+b, p, alpha)
    end
    function propwaldcicc(x::Int, n::Int, alpha::Float64)::ConfInt
        p=x/n
        b = quantile(ZDIST, 1-alpha/2)*sqrt(p*(1-p)/n)
        cc = 0.5/n
        return ConfInt(p-b-cc, p+b+cc, p, alpha)
    end
    #SOC  Second-Order corrected
    #T. Tony Cai One-sided confdence intervals in discrete distributions doi:10.1016/j.jspi.2004.01.00
    function propsocci(x::Int, n::Int, alpha::Float64)::ConfInt
        p  = x/n
        k  = quantile(ZDIST, 1-alpha/2)
        k2 = k^2
        η  = k2/3+1/6
        γ1 = -(k2*13/18+17/18)
        γ2 = k2/18+7/36
        m  = (x+η)/(n+2*η)
        b  = k*sqrt(p*(1-p)+(γ1*p*(1-p)+γ2)/n)/sqrt(n)
        return ConfInt(m-b, m+b, p, alpha)
    end
    #Arcsine
    function proparcci(x::Int, n::Int, alpha::Float64)::ConfInt
        q = quantile(ZDIST, 1-alpha/2)
        p = x/n
        z = q/(2*sqrt(n))
        return ConfInt(sin(asin(sqrt(p))-z)^2, sin(asin(sqrt(p))+z)^2, p, alpha)
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
        return ConfInt(ll, ul, (p1/(1-p1))/(p2/(1-p2)), alpha)
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

    @inline function mlemnor(φ::Float64, x1::Int, n1::Int, x2::Int, n2::Int)::Tuple{Float64,Float64}
        a  = n2*(φ-1)
        b  = φ*n1+n2-(x1+x2)*(φ-1)
        c  = -(x1 + x2)
        p2 = (-b+sqrt(b*b-4*a*c))/a/2
        p1 = p2*φ/(1+p2*(φ-1))
        return p1, p2
    end
    @inline function mnorval(φ::Float64, x1::Int, n1::Int, x2::Int, n2::Int, z)::Float64
        p1 = x1/n1
        pmle1, pmle2 = mlemnor(φ, x1, n1, x2, n2)
        return (n1*(p1-pmle1))^2 * (1/(n1*pmle1*(1-pmle1)) + 1/(n2*pmle2*(1-pmle2))) / ((n1+n2)/(n1+n2-1))  - z
    end
    function propORMNCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
        z        = quantile(Chisq(1), 1-alpha)
        fmnor(x) = mnorval(x, x1, n1, x2, n2, z)
        if (x1==0 && x2==0) || (x1==n1 && x2==n2)
            return ConfInt(0.0, Inf, NaN)
        elseif x1==0 || x2==n2
            return ConfInt(0.0, find_zero(fmnor, 1e-6, atol=1E-6), 0.0)
        elseif x1==n1 || x2 == 0
            return ConfInt(find_zero(fmnor, 1e-6, atol=1E-6), Inf, Inf)
        else
            estimate = (x1/(n1-x1))/(x2/(n2-x2))
            return ConfInt(find_zero(fmnor, 1e-6, atol=1E-6), find_zero(fmnor, estimate+1e-6, atol=1E-6), estimate, alpha)
        end
    end

    #Adjusted Woolf interval (Gart adjusted logit) Lawson, R (2005):Smallsample confidence intervals for the odds ratio.  Communication in Statistics Simulation andComputation, 33, 1095-1113.
    function propORaWoolfCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
            xa        = x1 + 0.5
            xb        = n1 - x1 + 0.5
            xc        = x2 + 0.5
            xd        = n2 - x2 + 0.5
            estimate  = xa*xd/xc/xb
            estI      = log(estimate)
            stde      = sqrt(1/xa + 1/xb + 1/xc + 1/xd)
            z         = quantile(ZDIST, 1-alpha/2)
            return ConfInt(exp(estI - z*stde), exp(estI + z*stde), estimate, alpha)
    end
    #Woolf logit
    #Woolf, B. (1955). On estimating the relation between blood group and disease. Annals of human genetics, 19(4):251-253.
    function propORWoolfCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
            xa       = x1
            xb       = n1 - x1
            xc       = x2
            xd       = n2 - x2
            estimate = xa*xd/xc/xb
            estI     = log(estimate)
            stde     = sqrt(1/xa + 1/xb + 1/xc + 1/xd)
            z        = quantile(ZDIST, 1-alpha/2)
            return ConfInt(exp(estI - z*stde), exp(estI + z*stde), estimate, alpha)
    end

    #Method of variance estimates recovery
    #Donner, A. and Zou, G. (2012). Closed-form confidence intervals for functions of the normal mean and standard deviation. Statistical Methods in Medical Research, 21(4):347-359.
    function propORMOVERCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
        p1       = (x1/(n1-x1))
        p2       = (x2/(n2-x2))
        estimate = p1/p2
        Z        = quantile(ZDIST, 1-alpha/2)
        wilci1   = propwilsonci(x1, n1, alpha)
        wilci2   = propwilsonci(x2, n2, alpha)
        wilci1   = ConfInt(wilci1.lower/(1-wilci1.lower), wilci1.upper/(1-wilci1.upper), estimate)
        wilci2   = ConfInt(wilci2.lower/(1-wilci2.lower), wilci2.upper/(1-wilci2.upper), estimate)
        lower    = (p1*p2-sqrt((p1*p2)^2 - wilci1.lower*wilci2.upper*(2*p1-wilci1.lower)*(2*p2-wilci2.upper)))/(wilci2.upper*(2*p2 - wilci2.upper))
        upper    = (p1*p2+sqrt((p1*p2)^2 - wilci1.upper*wilci2.lower*(2*p1-wilci1.upper)*(2*p2-wilci2.lower)))/(wilci2.lower*(2*p2 - wilci2.lower))
        return ConfInt(lower, upper, estimate, alpha)
    end

    #?Agresti independence-smoothed logit

    #?Cornfield exact confidence interval
    #https://rdrr.io/cran/ORCI/man/Cornfieldexact.CI.html



    #------------------------------DIFF-----------------------------------------
    #Wald
    function propDiffWaldCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
        p1       = x1/n1
        p2       = x2/n2
        estimate = p1-p2
        z        = quantile(ZDIST, 1 - alpha/2)
        stderr   = sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)
        return ConfInt(estimate - z*stderr, estimate + z*stderr, estimate, alpha)
    end
    #Wald CC
    function propDiffWaldCCCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
        p1       = x1/n1
        p2       = x2/n2
        estimate = p1-p2
        cc       = 0.5*(1/n1+1/n2)
        z        = quantile(ZDIST, 1 - alpha/2)
        stderr   = sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)
        return ConfInt(estimate - z*stderr - cc, estimate + z*stderr + cc, estimate, alpha)
    end

    #Newcombes Hybrid (wilson) Score interval for the difference of proportions
    #Newcombe 1998
    #10
    #Newcombe RG (1998): Interval Estimation for the Difference Between Independent Proportions: Comparison of Eleven Methods. Statistics in Medicine 17, 873-890.
    function propDiffNHSCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
        p1       = x1/n1
        p2       = x2/n2
        estimate = p1-p2
        z        = quantile(ZDIST, 1 - alpha/2)
        ci1      = propwilsonci(x1, n1, alpha)
        ci2      = propwilsonci(x2, n2, alpha)
        return ConfInt(estimate-z*sqrt(ci1.lower*(1-ci1.lower)/n1+ci2.upper*(1-ci2.upper)/n2),
                       estimate+z*sqrt(ci1.upper*(1-ci1.upper)/n1+ci2.lower*(1-ci2.lower)/n2), estimate, alpha)
    end

    function propDiffNHSCCCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
        p1       = x1/n1
        p2       = x2/n2
        estimate = p1-p2
        z        = quantile(ZDIST, 1 - alpha/2)
        ci1      = propwilsonccci(x1, n1, alpha)
        ci2      = propwilsonccci(x2, n2, alpha)
        return ConfInt(estimate-sqrt((p1-ci1.lower)^2 + (ci2.upper-p2)^2),
                       estimate+sqrt((ci1.upper - p1)^2 + (p2 - ci2.lower)^2), estimate, alpha)
    end

    #Agresti-Caffo interval for the difference of proportions
    #13
    #Agresti A, Caffo B., “Simple and effective confidence intervals for proportions and differences of proportions result from adding two successes and two failures”, American Statistician 54: 280–288 (2000)
    function propDiffACCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
        p1       = x1/n1
        p2       = x2/n2
        estimate = p1-p2
        z        = quantile(ZDIST, 1 - alpha/2)
        p1I      = (x1+1)/(n1+2)
        p2I      = (x2+1)/(n2+2)
        n1I      = n1+2
        n2I      = n2+2
        estI     = p1I-p2I
        stderr   = sqrt(p1I*(1-p1I)/n1I+p2I*(1-p2I)/n2I)
        return ConfInt(estI-z*stderr, estI+z*stderr, estimate, alpha)
    end
    #Method of Mee 1984 with Miettinen and Nurminen modification n / (n - 1) Newcombe 1998
    #Score intervals for the difference of two binomial proportions
    #6
    function propDiffMNCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
        ci       = propDiffNHSCCCI(x1, n1, x2, n2, alpha) #approx zero bounds
        p1       = x1/n1
        p2       = x2/n2
        estimate = p1 - p2
        z        = quantile(Chisq(1), 1-alpha)
        fmnd(x)  = mndiffval(p1, n1, p2, n2, estimate, x) - z
        if fmnd(ci.lower)*fmnd(estimate-1e-8) < 0.0
            ll = ci.lower
            lu = estimate-1e-8
        else
            ll = -1.0+1e-8
            lu = ci.lower
        end
        if fmnd(ci.upper)*fmnd(estimate+1e-8) < 0.0
            ul = estimate+1e-8
            uu = ci.upper
        else
            ul = ci.upper
            uu = 1.0-1e-8
        end
        return ConfInt(find_zero(fmnd, (ll, lu), atol=1E-6), find_zero(fmnd, (ul, uu), atol=1E-6), estimate, alpha)
    end
    @inline function mndiffval(p1::Float64, n1::Int, p2::Float64, n2::Int, estimate::Float64, Δ::Float64)::Float64
        return (estimate-Δ)^2/((n1+n2)/(n1+n2-1)*mlemndiff(p1, n1, p2, n2, Δ))
    end
    @inline function mlemndiff(p1::Float64, n1::Int, p2::Float64, n2::Int, Δ::Float64)::Float64
        if p1-p2-Δ == 0 return 0.0 end
        theta = n2/n1
        a = 1 + theta
        b = -(1 + theta + p1 + theta * p2 + Δ*(theta + 2))
        c = Δ^2 + Δ*(2 * p1 + theta + 1) + p1 + theta * p2
        d = -p1*Δ*(1 + Δ)
        v = (b/a/3)^3 - b*c/(6*a*a) + d/2/a
        u = sign(v) * sqrt((b/3/a)^2 - c/3/a)
        w = (pi + acos(v/u^3))/3
        p1n = 2*u*cos(w) - b/3/a
        p2n = p1n - Δ
        return p1n * (1-p1n)/n1 + p2n * (1 - p2n)/n2
    end
    #Mee
    #Mee RW (1984) Confidence bounds for the difference between two probabilities,Biometrics40:1175-1176
    #Farrington, C. P. and Manning, G. (1990), “Test Statistics and Sample Size Formulae for Comparative Binomial Trials with Null Hypothesis of Non-zero Risk Difference or Non-unity Relative Risk,” Statistics in Medicine, 9, 1447–1454
    #FM - Faster than propDiffMeeCI
    function propDiffFMCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
        p1       = x1/n1
        p2       = x2/n2
        estimate = p1 - p2
        z        = quantile(Chisq(1), 1-alpha)
        f(x)     = fmpval(p1, n1, p2, n2, estimate, x) - z
        return ConfInt(find_zero(f, (-1.0+1e-8, estimate-1e-8), atol=1E-6),
                       find_zero(f, (estimate+1e-8, 1.0-1e-8), atol=1E-6), estimate, alpha)
    end
    @inline function fmpval(p1::Float64, n1::Int, p2::Float64, n2::Int, estimate::Float64, Δ::Float64)
        return abs((estimate-Δ)^2/mlemndiff(p1, n1, p2, n2, Δ))
    end
    function propDiffMeeCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
        p1        = x1/n1
        p2        = x2/n2
        estimate  = p1 - p2
        f(x)      = fmpval2(p1, n1, p2, n2, estimate, x) - alpha
        return ConfInt(find_zero(f, (-1.0+1e-8, estimate-1e-8), atol=1E-6), find_zero(f, ( estimate+1e-6, 1.0-1e-6), atol=1E-6), estimate, alpha)
    end
    @inline function fmpval2(p1, n1, p2, n2, estimate, delta)::Float64
        z = (estimate - delta)/sqrt(mlemndiff(p1, n1, p2, n2, delta))
        p = cdf(ZDIST, z)
        return 2*min(1 - p, p)
    end


    # Brown, Li's Jeffreys ?
    #=
    function propDiffJeffCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
        p1   = (x1 + 0.5) / (n1+1)
        p2   = (x2 + 0.5) / (n2+1)
        v    = p1*(1-p1)/n1 + p2*(1-p2)/n2
        z    = quantile(ZDIST, 1 - alpha/2)
        t    = z * sqrt(v)
        estI = p1-p2
        return ConfInt(max(-1, estI - t), min(1, estI + t), x1/n1-x2/n2)
    end
    =#

    #--------------------------------RR-----------------------------------------
    #Miettinen-Nurminen Score interval
    #Miettinen, O. and Nurminen, M. (1985), Comparative analysis of two rates. Statist. Med., 4: 213-226. doi:10.1002/sim.4780040211
    @inline function mlemnrr(φ::Float64, x1::Int, n1::Int, x2::Int, n2::Int)::Tuple{Float64,Float64}
        a = (n1+n2)*φ
        b = -(φ*(x1+n2)+x2+n1)
        c = x1 + x2
        p2 = (-b-sqrt(b*b-4*a*c))/2/a
        p1 = p2*φ
        return p1, p2
    end
    @inline function mnrrval(φ::Float64, x1::Int, n1::Int, x2::Int, n2::Int, z::Float64)::Float64
        p1 = x1/n1
        p2 = x2/n2
        pmle1, pmle2 = mlemnrr(φ, x1, n1, x2, n2)
        return (p1 - φ*p2)^2/((pmle1*(1-pmle1)/n1 + φ*φ*pmle2*(1-pmle2)/n2)*((n1+n2-1)/(n1+n2)))-z
    end
    function propRRMNCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
        z        = quantile(Chisq(1), 1 - alpha)
        fmnrr(x) = mnrrval(x, x1, n1, x2, n2, z)
        if (x1==0 && x2==0) || (x1==n1 && x2==n2)
            return ConfInt(0.0, Inf, NaN)
        elseif x1==0 || x2==n2
            return ConfInt(0.0, find_zero(fmnrr, 1e-8, atol=1E-6), 0.0, alpha)
        elseif x1==n1 || x2 == 0
            return ConfInt(find_zero(fmnrr, 1e-8, atol=1E-6), Inf, Inf, alpha)
        else
            estimate = (x1/n1)/(x2/n2)
            return ConfInt(find_zero(fmnrr, 1e-8, atol=1E-6), find_zero(fmnrr, estimate+1e-6, atol=1E-6), estimate, alpha)
        end
    end #propRRMNCI

    #Katz D, Baptista J, Azen SP and Pike MC. Obtaining confidence intervals for the risk ratio in cohort studies. Biometrics 1978; 34: 469–474
    function proprrkatzci(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
        estimate  = (x1/n1)/(x2/n2)
        estI      = log(estimate)
        stderrlog = sqrt(1/x2+1/x1-1/n2-1/n1)
        Z         = quantile(ZDIST,1-alpha/2)
        return ConfInt(exp(estI-Z*stderrlog), exp(estI+Z*stderrlog), estimate, alpha)
    end

    #Crude log interval
    #walters
    #Gart, JJand Nam, J (1988): Approximate interval estimation of the ratio of binomial parameters: Areview and corrections for skewness. Biometrics 44, 323-338.
    function propRRCLICI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
        x1I       = x1+0.5
        x2I       = x2+0.5
        n1I       = n1+0.5
        n2I       = n2+0.5
        estI      = log((x1I/n1I)/(x2I/n2I))
        stderrlog = sqrt(1/x2I+1/x1I-1/n2I-1/n1I)
        estimate  = (x1/n1)/(x2/n2)
        Z         =  quantile(ZDIST,1-alpha/2)
        return ConfInt(exp(estI-Z*stderrlog), exp(estI+Z*stderrlog), estimate, alpha)
    end
    #Method of variance estimates recovery (Donner, Zou, 2012)
    function propRRMOVERCI(x1::Int, n1::Int, x2::Int, n2::Int, alpha::Float64)::ConfInt
        p1       = (x1/n1)
        p2       = (x2/n2)
        estimate = (x1/n1)/(x2/n2)
        Z        = quantile(ZDIST, 1-alpha/2)
        wilci1   = propwilsonci(x1, n1, alpha)
        wilci2   = propwilsonci(x2, n2, alpha)
        lower    = (p1*p2-sqrt((p1*p2)^2 - wilci1.lower*wilci2.upper*(2*p1-wilci1.lower)*(2*p2-wilci2.upper)))/(wilci2.upper*(2*p2 - wilci2.upper))
        upper    = (p1*p2+sqrt((p1*p2)^2 - wilci1.upper*wilci2.lower*(2*p1-wilci1.upper)*(2*p2-wilci2.lower)))/(wilci2.lower*(2*p2 - wilci2.lower))
        return ConfInt(lower, upper, estimate, alpha)
    end
    #-------------------------------MEANS---------------------------------------

    #Normal
    function meanNormCI(m::Real, s::Real, n::Real, alpha::Real)::ConfInt
        e = quantile(ZDIST, 1-alpha/2)*sqrt(s/n)
        return ConfInt(m-e, m+e, m, alpha)
    end
    #T Distribution
    function meanTdistCI(m::Real, s::Real, n::Real, alpha::Real)::ConfInt
        e = quantile(TDist(n-1), 1-alpha/2)*sqrt(s/n)
        return ConfInt(m-e, m+e, m, alpha)
    end
    #mean diff equal var
    function meanDiffEV(m1::Real, s1::Real, n1::Real, m2::Real, s2::Real, n2::Real, alpha::Real)::ConfInt
        diff   = m1 - m2
        stddev = sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
        stderr = stddev * sqrt(1/n1 + 1/n2)
        d      = stderr*quantile(TDist(n1+n2-2), 1-alpha/2)
        return ConfInt(diff-d, diff+d, diff, alpha)
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
        return ConfInt(diff-d, diff+d, diff, alpha)
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

"""
    diffcmhci(data::DataFrame; a = :a, b = :b, c = :c, d = :d, alpha = 0.05, method = :default)::ConfInt

Cochran–Mantel–Haenszel confidence intervals for proportion difference.

**data**- dataframe with 4 columns, each line represent 2X2 table

**a** **b** **c** **d** - dataframe table names (number of subjects in 2X2 table):

"""
    function diffcmhci(data::DataFrame; a = :a, b = :b, c = :c, d = :d, alpha = 0.05, method = :default)::ConfInt
        n1 = data[:, a] + data[:, b]
        n2 = data[:, c] + data[:, d]
        N  = data[:, a] + data[:, b] + data[:, c] + data[:, d]
        z = quantile(ZDIST, 1 - alpha/2)

        estimate = sum(data[:, a] .* (n2 ./ N) - data[:, c] .* (n1 ./ N)) / sum(n1 .* (n2 ./ N))
        #1
        if method == :sato || method == :default
            se   = sqrt((estimate * (sum(data[:, c] .* (n1 ./ N) .^ 2 - data[:, a] .* (n2 ./ N) .^2 + (n1 ./ N) .* (n2 ./ N) .* (n2-n1) ./ 2)) + sum(data[:, a] .* (n2 - data[:, c]) ./ N + data[:, c] .* (n1 - data[:, a]) ./ N)/2) / sum(n1 .* (n2 ./ N)) .^ 2) # equation in: Sato, Greenland, & Robins (1989)
        #2
        elseif method == :gr
            se   = sqrt(sum(((data[:, a] ./N .^2) .* data[:, b] .* (n2 .^2 ./ n1) + (data[:, c] ./N .^2) .* data[:, d] .* (n1 .^2 ./ n2))) / sum(n1 .*(n2 ./ N))^2) # equation in: Greenland & Robins (1985)
        else
            throw(ArgumentError("Method unknown!"))
        end
        return ConfInt(estimate - z*se, estimate + z*se, estimate, alpha)
    end

    """
        orcmhci(data::DataFrame; a = :a, b = :b, c = :c, d = :d, alpha = 0.05, logscale = false)::ConfInt

    Cochran–Mantel–Haenszel confidence intervals for odd ratio.
    """
    function orcmhci(data::DataFrame; a = :a, b = :b, c = :c, d = :d, alpha = 0.05, logscale = false)::ConfInt
        #n1 = data[:, a] + data[:, b]
        #n2 = data[:, c] + data[:, d]
        N  = data[:, a] + data[:, b] + data[:, c] + data[:, d]
        z = quantile(ZDIST, 1 - alpha/2)

        Pi = (data[:, a] ./ N) + (data[:, d] ./ N)
        Qi = (data[:, b] ./ N) + (data[:, c] ./ N)
        Ri = (data[:, a] ./ N) .* data[:, d]
        Si = (data[:, b] ./ N) .* data[:, c]
        R  = sum(Ri)
        S  = sum(Si)
        if R == 0 || S == 0 return ConfInt(NaN, NaN, NaN) end
        estimate = log(R/S)
        se  = sqrt(1/2 * (sum(Pi .* Ri)/R^2 + sum(Pi .* Si + Qi .* Ri)/(R*S) + sum(Qi .* Si)/S^2)) # based on Robins et al. (1986)
        #zval= est / se
        #pval= 2*(1-cdf(Normal(), abs(zval)))
        if logscale return ConfInt(estimate - z*se, estimate + z*se, estimate, alpha) else return ConfInt(exp(estimate - z*se), exp(estimate + z*se), exp(estimate), alpha) end
    end

    """
        rrcmhci(data::DataFrame; a = :a, b = :b, c = :c, d = :d, alpha = 0.05, logscale = false)::ConfInt

    Cochran–Mantel–Haenszel confidence intervals for risk ratio.
    """
    function rrcmhci(data::DataFrame; a = :a, b = :b, c = :c, d = :d, alpha = 0.05, logscale = false)::ConfInt
        n1 = data[:, a] + data[:, b]
        n2 = data[:, c] + data[:, d]
        N  = data[:, a] + data[:, b] + data[:, c] + data[:, d]
        z = quantile(ZDIST, 1 - alpha/2)
            #...
        R = sum(data[:, a] .* (n2 ./ N))
        S = sum(data[:, c] .* (n1 ./ N))
        if sum(data[:, a]) == 0 || sum(data[:, c]) == 0 return ConfInt(NaN, NaN, NaN) end
        estimate = log(R/S)
        se  = sqrt(sum(((n1 ./ N) .* (n2 ./ N) .* (data[:, a] + data[:, c]) - (data[:, a] ./ N) .* data[:, c])) / (R*S))
        #zval= est / se
        #pval= 2*(1-cdf(Normal(), abs(zval)))
        if logscale return ConfInt(estimate  - z*se, estimate  + z*se, estimate , alpha) else return ConfInt(exp(estimate  - z*se), exp(estimate  + z*se), exp(estimate ), alpha) end
    end
