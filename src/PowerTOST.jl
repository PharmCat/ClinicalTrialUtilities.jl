
function powerTOST(;alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0, n=0, design="2x2", method="owenq")

    #if isnan(cv) || isnan(n) return false end
    if n < 2 return false end
    if cv == 0 return false end
    if !(0 < alpha < 1) return false end

    df = n-2
    n1  = n ÷ 2
    n2  = n1 + n % 2
    sef = sqrt((1/n1+1/n2)*0.5)

    if df < 1 return false end

    if logscale
        ltheta1 = log(theta1)
        ltheta2 = log(theta2)
        diffm   = log(theta0)
        se      = cv2se(cv)
    else
        ltheta1 = theta1
        ltheta2 = theta2
        diffm   = theta0
        se      = cv
    end

    sem = se*sef

    if method=="owenq"
        return powerTOSTOwenQ(alpha,ltheta1,ltheta2,diffm,sem,df)
    elseif method=="nct"
        return approxPowerTOST(alpha,ltheta1,ltheta2,diffm,sem,df)
    elseif method=="mvt"
        return power1TOST(alpha,ltheta1,ltheta2,diffm,sem,df) #not implemented
    elseif method=="shifted"
        return approx2PowerTOST(alpha,ltheta1,ltheta2,diffm,sem,df)
    else
        return false
    end
end #powerTOST

#.power.TOST
function powerTOSTOwenQ(alpha,ltheta1,ltheta2,diffm,sem,df)
    tval = qt(1-alpha, df)

    delta1 = (diffm-ltheta1)/sem
    delta2 = (diffm-ltheta2)/sem

    R = (delta1-delta2)*sqrt(df)/(tval+tval) #not clear R implementation of vector form

    if isnan(R) R = 0 end

    if R <= 0 R = Inf end

    # to avoid numerical errors in OwensQ implementation
    # 'shifted' normal approximation Jan 2015
    # former Julious formula (57)/(58) doesn't work
    if df > 10000
        tval = qnorm(1-alpha)
        p1   = pnorm(tval-delta1)
        p2   = pnorm(-tval-delta2)

        pwr = p2-p1
        if pwr > 0 return pwr else return 0 end
    elseif df >= 5000
        # approximation via non-central t-distribution
        return approxPowerTOST(alpha,ltheta1,ltheta2,diffm,sem,df)
    end
    # get correct values for p1 and p2:
    # p1 is for left hypothesis (right-tailed)
    # -> critical value is always first component
    # p2 for right hypothesis (left-tailed)
    # -> critical value is first component if alpha is 1-dim,
    #    second component if alpha is 2-dim
    # NOT VECTORIZED

    p1 = owensQ(df, tval, delta1, 0, R)
    p2 = owensQ(df,-tval, delta2, 0, R)

    pwr = p2 - p1
    if pwr > 0 return pwr else return 0 end

end #powerTOSTOwenQ

#------------------------------------------------------------------------------
# 'raw' approximate power function without any error checks,
# approximation based on non-central t
# this vectorizes ok
# .approx.power.TOST - PowerTOST
function approxPowerTOST(alpha,ltheta1,ltheta2,diffm,sem,df)
    #pt(x,df,ncp)
    #NCT=NoncentralT(df,ncp)
    #cdf(NCT,x)
    tval = qt(1-alpha,df)
    delta1 = (diffm-ltheta1)/sem
    delta2 = (diffm-ltheta2)/sem
    pow = cdf(NoncentralT(df,delta2), -tval) - cdf(NoncentralT(df,delta1), tval)
    if pow > 0 return pow else return 0 end
end #approxPowerTOST

#.power.1TOST
function power1TOST(alpha,ltheta1,ltheta2,diffm,sem,df)
    return false #Method ON MULTIVARIATE t AND GAUSS PROBABILITIES IN R not implemented
    # Distributions.MvNormal - in plan
    # see  Distributions.jl/src/multivariate/mvtdist.jl
    # Multivariate t-distribution
    # Generic multivariate t-distribution class
    # mvt = MvTDist()

end

#.approx2.power.TOST
function approx2PowerTOST(alpha,ltheta1,ltheta2,diffm,sem,df)
    tval = qt(1-alpha, df)
    delta1 = (diffm-ltheta1)/sem
    delta2 = (diffm-ltheta2)/sem
    if isnan(delta1) delta1 = 0 end
    if isnan(delta2) delta2 = 0 end
    #pt(x, df)
    #TDIST=TDist(df)
    #cdf(TDIST, x)
    #pow <- pt(-tval[length(tval)]-delta2, df) - pt(tval[1]-delta1, df)
    tdist=TDist(df)
    pow = cdf(tdist,-tval-delta2) - cdf(tdist,tval-delta1)
    if pow > 0 return pow else return 0 end
end #approx2PowerTOST

#CV2se
function cv2se(cv)
    return sqrt(log(1+cv^2))
end
