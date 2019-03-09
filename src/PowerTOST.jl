# Clinical Trial Utilities
# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)

# ms = ss / df
# sd = ms^2
# se_diff = se = sd * sqrt((1/N1 + ... + 1/Nn)*bkni) = sqrt(ms*(1/N1 + ... + 1/Nn)*bkni)
# CI bounds = Diff +- t(df, alpha)*se

function powerTOST(;alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.0, n=0, design=:d2x2, method=:owenq,  out=:num)
    if n < 2 throw(CTUException(1021,"powerTOST: n can not be < 2")) end
    if cv == 0 throw(CTUException(1022,"powerTOST: cv can not be equal to 0"))  end
    if !(0 < alpha < 1) throw(CTUException(1023,"powerTOST: alfa can not be > 1 or < 0")) end
    theta0   = convert(Float64, theta0)
    theta1   = convert(Float64, theta1)
    theta2   = convert(Float64, theta2)
    logscale = convert(Bool, logscale)
    cv       = convert(Float64, cv)
    n        = convert(Int, n)
    alpha    = convert(Float64, alpha)

    if logscale
        ltheta1 = log(theta1)
        ltheta2 = log(theta2)
        diffm   = log(theta0)
        sd      = cv2sd(cv)    # sqrt(ms)
    else
        ltheta1 = theta1;
        ltheta2 = theta2;
        diffm   = theta0;
        sd      = cv;
    end

    return powerTOSTint(alpha,  ltheta1, ltheta2, diffm, sd, n, design, method)
end

function powerTOSTint(alpha::Float64,  ltheta1::Float64, ltheta2::Float64, diffm::Float64, sd::Float64, n::Int, design::Symbol, method::Symbol)::Float64

    dffunc, bkni, seq = designProp(design) #dffunc if generic funtion with 1 arg return df
    df    = dffunc(n)
    sqa   = Array{Float64, 1}(undef, seq)
    sqa  .= n÷seq
    for i = 1:n%seq
        sqa[i] += 1
    end
    sef = sqrt(sum(1 ./ sqa)*bkni)

    if df < 1 throw(CTUException(1024,"powerTOSTint: df < 1")) end

    se::Float64 = sd*sef

    if method     == :owenq
        return powerTOSTOwenQ(alpha,ltheta1,ltheta2,diffm,se,df)
    elseif method == :nct
        return approxPowerTOST(alpha,ltheta1,ltheta2,diffm,se,df)
    elseif method == :mvt
        return power1TOST(alpha,ltheta1,ltheta2,diffm,se,df) #not implemented
    elseif method == :shifted
        return approx2PowerTOST(alpha,ltheta1,ltheta2,diffm,se,df)
    else
         throw(CTUException(1025,"powerTOST: method not known!"))
    end
end #powerTOST

#.power.TOST
function powerTOSTOwenQ(alpha,ltheta1,ltheta2,diffm,se,df)
    tval::Float64   = quantile(TDist(df), 1-alpha) #qt(1-alpha, df)
    delta1::Float64 = (diffm-ltheta1)/se
    delta2::Float64 = (diffm-ltheta2)/se
    R::Float64      = (delta1-delta2)*sqrt(df)/(tval+tval) #not clear R implementation of vector form
    if isnan(R) R   = 0 end
    if R <= 0 R     = Inf end

    # to avoid numerical errors in OwensQ implementation
    # 'shifted' normal approximation Jan 2015
    # former Julious formula (57)/(58) doesn't work
    if df > 10000
        #tval = qnorm(1-alpha)
        tval  = quantile(ZDIST, 1-alpha)
        #p1   = pnorm(tval-delta1)
        p1    = cdf(ZDIST, tval-delta1)
        #p2   = pnorm(-tval-delta2)
        p2    = cdf(ZDIST, -tval-delta2)
        pwr   = p2-p1
        if pwr > 0 return pwr else return 0 end
    elseif df >= 5000
        # approximation via non-central t-distribution
        return approxPowerTOST(alpha,ltheta1,ltheta2,diffm,se,df)
    end
    p1  = owensQ(df, tval, delta1, 0.0, R)
    p2  = owensQ(df,-tval, delta2, 0.0, R)
    pwr = p2 - p1
    if pwr > 0 return pwr else return 0 end
end #powerTOSTOwenQ

#------------------------------------------------------------------------------
# 'raw' approximate power function without any error checks,
# approximation based on non-central t
# .approx.power.TOST - PowerTOST
function approxPowerTOST(alpha,ltheta1,ltheta2,diffm,se,df)
    tdist           = TDist(df)
    tval::Float64   = quantile(tdist, 1-alpha)
    delta1::Float64 = (diffm-ltheta1)/se
    delta2::Float64 = (diffm-ltheta2)/se
    pow             = cdf(NoncentralT(df,delta2), -tval) - cdf(NoncentralT(df,delta1), tval)
    if pow > 0 return pow else return 0 end
end #approxPowerTOST

#.power.1TOST
function power1TOST(alpha,ltheta1,ltheta2,diffm,se,df)
    throw(CTUException(1000,"Method not implemented!"))
    #Method ON MULTIVARIATE t AND GAUSS PROBABILITIES IN R not implemented
    # Distributions.MvNormal - in plan
    # see  Distributions.jl/src/multivariate/mvtdist.jl
    # Multivariate t-distribution
    # Generic multivariate t-distribution class
    # mvt = MvTDist()
end

#.approx2.power.TOST
function approx2PowerTOST(alpha,ltheta1,ltheta2,diffm,se,df)
    tdist           = TDist(df)
    tval::Float64   = quantile(tdist, 1-alpha)
    delta1::Float64 = (diffm-ltheta1)/se
    delta2::Float64 = (diffm-ltheta2)/se
    if isnan(delta1) delta1 = 0 end
    if isnan(delta2) delta2 = 0 end
    pow = cdf(tdist,-tval-delta2) - cdf(tdist,tval-delta1)
    if pow > 0 return pow else return 0 end
end #approx2PowerTOST

#CV2se
@inline function cv2sd(cv::Float64)::Float64
    return sqrt(log(1+cv^2))
end

@inline function cv2ms(cv::Float64)::Float64
    return log(1+cv^2)
end
@inline function ms2cv(ms::Float64)::Float64
    return sqrt(exp(ms)-1)
end
#CV2se
@inline function sd2cv(sd::Float64)::Float64
    return sqrt(exp(sd^2)-1)
end

function designProp(type::Symbol)
    if type == :parallel
        function f1(n) n - 2 end
        return f1, 1.0, 2
    elseif type == :d2x2
        function f2(n) n - 2 end
        return f2, 0.5, 2
    elseif type == :d2x2x3
        return function f3(n) 2*n - 3 end, 0.375, 2
    elseif type == :d2x2x4
        return function f4(n) 3*n - 4 end, 0.25, 2
    elseif type == :d2x4x4
        return function f5(n) 3*n - 4 end, 0.0625, 4
    elseif type == :d2x3x3
        return function f6(n) 2*n - 3 end, 1/6, 3
    else throw(CTUException(1031,"designProp: design not known!")) end
end
