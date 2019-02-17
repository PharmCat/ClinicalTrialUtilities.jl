
module CI
    using Distributions
    import ClinicalTrialUtilities.ZDIST
    import ClinicalTrialUtilities.CTUException

    struct ConfInt
        lower::Float64
        upper::Float64
        estimate::Float64
    end

    function oneProp(n, na; alpha=0.05, method="wilson")
        if method=="wilson"
            return propWilsonCI(n, na, alpha)
        elseif method=="cp"
            return propCPCI(n, na, alpha)
        elseif method=="soc"
            return propSOCCI(n, na, alpha)
        elseif method=="blaker"
            return propBlakerCI(n, na, alpha)
        end
    end

    function oneMeans(m,s,n,alpha; method="norm")
        if method=="norm"
            meanNormCI(m,s,n,alpha)
        elseif method=="tdist"
            meanTdistCI(m,s,n,alpha)
        end
    end

    function twoProp(n1, n1a, n2, n2a; alpha=0.05, type="diff", method="wilson")

        if type=="diff"
            if method == ""
            end
        elseif type=="rr"
            return
        elseif type=="or"
            return
        end
    end

    function twoMeans()

    end

    #Wilson’s confidence interval for a single proportion
    #Wilson, E.B. (1927) Probable inference, the law of succession, and statistical inferenceJ. Amer.Stat. Assoc22, 209–212
    function propWilsonCI(n, na, alpha)
        z = abs(quantile(ZDIST, 1-alpha/2))
        p = n/na
        b = (z*((p*(1-p)+(z^2)/(4*na))/na)^(1/2))/(1+(z^2)/na)
        m = (p+(z^2)/(2*na))/(1+(z^2)/na)
        return ConfInt(m - b,m + b,m)
    end
    #Clopper-Pearson exatct CI
    #Clopper, C. and Pearson, E.S. (1934) The use of confidence or fiducial limits illustrated in the caseof the binomial.Biometrika26, 404–413.
    function propCPCI(n, na, alpha)
        if n==0
            ll=0.0
            ul=1.0-(alpha/2)^(1/na)
        elseif n==na
            ul=1.0
            ll=(alpha/2)^(1/na)
        else
            ll = 1/(1+(na-n+1)/(n*quantile(FDist(2*n, 2*(na-n+1)), alpha/2)))
            ul = 1/(1+(na-n) / ((n+1)*quantile(FDist(2*(n+1), 2*(na-n)), 1-alpha/2)))
        end
        return ll, ul
    end
    #Blaker, H. (2000). Confidence curves and improved exact confidence intervals for discrete distribu-tions,Canadian Journal of Statistics28 (4), 783–798
    function propBlakerCI(n, na, alpha)
        tol = 1E-5
        lower = 0
        upper = 1
        if n != 0
            lower = quantile(Beta(n, na-n+1), alpha/2)
            while acceptbin(n, na, lower+tol) < alpha
                lower +=tol
            end
        end
        if n != na
            upper = quantile(Beta(n+1, na-n), 1-alpha/2)
            while acceptbin(n, na, upper-tol) < alpha
                upper -=tol
            end
        end
        return ConfInt(lower,upper, n/na)
    end
    function acceptbin(n,na,p)
        BIN = Binomial(na,p)
        p1 = 1-cdf(BIN,n-1)
        p2 =   cdf(BIN,n)
        a1 = p1 + cdf(BIN, quantile(BIN,p1)-1)
        a2 = p2+1-cdf(BIN, quantile(BIN,1-p2))
        return min(a1,a2)
    end
    #Wald
    function propWaldCI(n, na, alpha)
        p=n/na
        b = quantile(ZDIST, 1-alpha/2)*sqrt(p*(1-p)/na)
        return ConfInt(p-b,p+b,p)
    end
    #SOC  Second-Order corrected
    #T. Tony Cai One-sided con&dence intervals in discretedistributions doi:10.1016/j.jspi.2004.01.00
    #not clear implementation
    function propSOCCI(n, na, alpha)
        p=n/na
        k=quantile(ZDIST, 1-alpha/2)
        k2=k^2
        η=k2/3+1/6
        γ1=-(k2*13/18+17/18)
        γ2= k2/18+7/36
        m=(n+η)/(na+2*η)
        b=k*sqrt(p*(1-p)+(γ1*p*(1-p)+γ2)/na)/sqrt(na)
        return ConfInt(m-b, m+b, p)
    end

    #Normal
    function meanNormCI(m,s,n,alpha)
        e = quanlile(ZDIST, 1-alpha/2)*s/sqrt(n)
        return m-e, m+e
    end
    #T Distribution
    function meanTdistCI(m,s,n,alpha)
        e = quanlile(TDist(n-1), 1-alpha/2)*s/sqrt(n)
        return m-e, m+e
    end
end
