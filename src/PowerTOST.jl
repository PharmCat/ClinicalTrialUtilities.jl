# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)

function powerTOST(;alpha=0.05, logscale=true, theta1=0.8, theta2=1.25, theta0=0.95, cv=0.0, n=0, design="2x2", method="owenq")
    if n < 2 return throw(CTUException(1021,"powerTOST: n can not be < 2")) end
    if cv == 0 throw(CTUException(1022,"powerTOST: cv can not be equal to 0"))  end
    if !(0 < alpha < 1) return throw(CTUException(1023,"powerTOST: alfa can not be > 1 or < 0")) end
    theta0 = convert(Float64, theta0)
    theta1 = convert(Float64, theta1)
    theta2 = convert(Float64, theta2)
    logscale = convert(Bool, logscale)
    cv = convert(Float64, cv)
    n = convert(Int, n)
    alpha = convert(Float64, alpha)

    powerTOSTint(alpha, logscale, theta1, theta2, theta0, cv, n, design, method)
end

function powerTOSTint(alpha::Float64, logscale::Bool, theta1::Float64, theta2::Float64, theta0::Float64, cv::Float64, n::Int, design::String, method::String)::Float64

    #if isnan(cv) || isnan(n) return false end

    #df::Int = n-2

    df, bkni, seq = designProp(n, design)
    sqa   = Array{Float64, 1}(undef, seq)
    sqa .= n÷seq
    for i = 1:n%seq
        sqa[i] += 1
    end

    sef = sqrt(sum(1 ./ sqa)*bkni)

    if df < 1 return  throw(CTUException(1024,"powerTOST: df < 1")) end

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

    sem::Float64 = se*sef

    if method=="owenq"
        return powerTOSTOwenQ(alpha,ltheta1,ltheta2,diffm,sem,df)
    elseif method=="nct"
        return approxPowerTOST(alpha,ltheta1,ltheta2,diffm,sem,df)
    elseif method=="mvt"
        return power1TOST(alpha,ltheta1,ltheta2,diffm,sem,df) #not implemented
    elseif method=="shifted"
        return approx2PowerTOST(alpha,ltheta1,ltheta2,diffm,sem,df)
    else
         throw(CTUException(1025,"powerTOST: method "*method*" not known!"))
    end
end #powerTOST

#.power.TOST
function powerTOSTOwenQ(alpha,ltheta1,ltheta2,diffm,sem,df)
    tval::Float64 = qt(1-alpha, df)

    delta1::Float64 = (diffm-ltheta1)/sem
    delta2::Float64 = (diffm-ltheta2)/sem

    R::Float64 = (delta1-delta2)*sqrt(df)/(tval+tval) #not clear R implementation of vector form

    if isnan(R) R = 0 end

    if R <= 0 R = Inf end

    # to avoid numerical errors in OwensQ implementation
    # 'shifted' normal approximation Jan 2015
    # former Julious formula (57)/(58) doesn't work
    if df > 10000
        #tval = qnorm(1-alpha)
        tval  = quantile(ZDIST, 1-alpha)
        #p1   = pnorm(tval-delta1)
        p1    = cdf(ZDIST, tval-delta1)
        #p2    = pnorm(-tval-delta2)
        p2    = cdf(ZDIST, -tval-delta2)

        pwr = p2-p1
        if pwr > 0 return pwr else return 0 end
    elseif df >= 5000
        # approximation via non-central t-distribution
        return approxPowerTOST(alpha,ltheta1,ltheta2,diffm,sem,df)
    end
    # NOT VECTORIZED

    p1 = owensQ(df, tval, delta1, 0.0, R)
    p2 = owensQ(df,-tval, delta2, 0.0, R)

    pwr = p2 - p1
    if pwr > 0 return pwr else return 0 end

end #powerTOSTOwenQ

#------------------------------------------------------------------------------
# 'raw' approximate power function without any error checks,
# approximation based on non-central t
# .approx.power.TOST - PowerTOST
function approxPowerTOST(alpha,ltheta1,ltheta2,diffm,sem,df)

    tval::Float64 = qt(1-alpha,df)
    delta1::Float64 = (diffm-ltheta1)/sem
    delta2::Float64 = (diffm-ltheta2)/sem
    pow = cdf(NoncentralT(df,delta2), -tval) - cdf(NoncentralT(df,delta1), tval)
    if pow > 0 return pow else return 0 end
end #approxPowerTOST

#.power.1TOST
function power1TOST(alpha,ltheta1,ltheta2,diffm,sem,df)

    throw(CTUException(1000,"Method not implemented!"))
    #Method ON MULTIVARIATE t AND GAUSS PROBABILITIES IN R not implemented
    # Distributions.MvNormal - in plan
    # see  Distributions.jl/src/multivariate/mvtdist.jl
    # Multivariate t-distribution
    # Generic multivariate t-distribution class
    # mvt = MvTDist()

end

#.approx2.power.TOST
function approx2PowerTOST(alpha,ltheta1,ltheta2,diffm,sem,df)
    tval::Float64 = qt(1-alpha, df)
    delta1::Float64 = (diffm-ltheta1)/sem
    delta2::Float64 = (diffm-ltheta2)/sem
    if isnan(delta1) delta1 = 0 end
    if isnan(delta2) delta2 = 0 end
    tdist=TDist(df)
    pow = cdf(tdist,-tval-delta2) - cdf(tdist,tval-delta1)
    if pow > 0 return pow else return 0 end
end #approx2PowerTOST

#CV2se
function cv2se(cv::Float64)::Float64
    return sqrt(log(1+cv^2))
end

function designProp(n::Int, type::String)
    if type == "parallel"
        return n - 2, 1.0, 2
    elseif type == "2x2"
        return n - 2, 0.5, 2
    elseif type == "2x2x3"
        return 2*n - 3, 0.375, 2
    elseif type == "2x2x4"
        return 3*n - 4, 0.25, 2
    elseif type == "2x4x4"
        return 3*n - 4, 0.0625, 4
    elseif type == "2x3x3"
        return 2*n - 3, 1/6, 3
    else throw(CTUException(1031,"designProp: design not known!")) end
end
