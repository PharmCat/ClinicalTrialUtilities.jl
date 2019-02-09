# Clinical Trial Power and Sample Size calculation
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)


#main sample size function
function sampleSize(;param="", type="", group="", alpha=0.05, beta=0.2, diff=0, sd=0, a=0, b=0, k=1, out="num")
    if alpha >= 1 || alpha <= 0 || beta >= 1 || beta <= 0 return false end
    if (type == "ei" || type == "ns") && diff == 0 return false end
    if param == "prop" && !(group == "one" || group == "two" || type == "mcnm")  return false end
    if sd == 0 && param == "mean" return false end
    if k == 0 return false end
    if param == "mean"
        if group == "one"
            if type == "ea"
                n = oneSampleMeanEquality(a, b, sd, alpha=alpha, beta=beta)
            elseif type == "ei"
                n = oneSampleMeanEquivalence(a, b, sd, diff, alpha=alpha, beta=beta)
            elseif type == "ns"
                n = oneSampleMeanNS(a, b, sd, diff, alpha=alpha, beta=beta)
            else return false end
        elseif group == "two"
            if type == "ea"
                n = twoSampleMeanEquality(a, b, sd, alpha=alpha, beta=beta, k=k)
            elseif type == "ei"
                n = twoSampleMeanEquivalence(a, b, sd, diff, alpha=alpha, beta=beta, k=k)
            elseif type == "ns"
                n = twoSampleMeanNS(a, b, sd, diff, alpha=alpha, beta=beta, k=k)
            else return false end
        else return false end
    elseif param == "prop"
        if 1 < a || a < 0 || 1 < b || b < 0 return false end
        if type == "mcnm"
            n = mcnm(a, b; alpha=alpha, beta=beta)
        else
            if group == "one"
                if type == "ea"
                    n = oneProportionEquality(a, b; alpha=alpha, beta=beta)
                elseif type == "ei"
                    n = oneProportionEquivalence(a, b, diff; alpha=alpha, beta=beta)
                elseif type == "ns"
                    n = oneProportionNS(a, b, diff; alpha=alpha, beta=beta)
                else return false end
            elseif group == "two"
                if type == "ea"
                    n = twoProportionEquality(a, b; alpha=alpha, beta=beta, k=k)
                elseif type == "ei"
                    n = twoProportionEquivalence(a, b, diff; alpha=alpha, beta=beta, k=k)
                elseif type == "ns"
                    n = twoProportionNS(a, b, diff; alpha=alpha, beta=beta, k=k)
                else return false end
            else return false end
        end
    elseif param == "or"
        if type == "ea"
            n = orEquality(a, b; alpha=alpha, beta=beta, k=k)
        elseif type == "ei"
            n = orEquivalence(a, b, diff; alpha=alpha, beta=beta, k=k, logdiff=true)
        elseif type == "ns"
            n = orNS(a, b, diff; alpha=alpha, beta=beta, k=k, logdiff=true)
        else return false end
    else return false end

    if out == "num" return n
    else
        params = "" #string parameter type
        if param == "mean" params = "Mean" elseif param== "prop" params = "Proportion" elseif param == "or" params = "Odd Ration" end
        if group == "one" groups = "One" elseif  group == "two" groups = "Two" else groups = "NA" end
        if type == "ea" hyps = "Equality" elseif type == "ei" hyps = "Equivalence" elseif type == "ns" hyps = "Non-Infriority/Superiority" elseif type == "mcnm" hyps = "McNemar's Equality test" end
        output  = "----------------------------------------\n"
        output *= "         Sample Size Estimation        \n"
        output *= "----------------------------------------\n"
        output *= "  Paramether type: "*params*"\n"
        output *= "  Groups: "*groups*"\n"
        output *= "  Hypothesis: "*hyps*"\n"
        output *= "----------------------------------------\n"
        output *= "  Alpha: "*string(alpha)*" Beta: "*string(beta)*"\n"
        if group == "two" output *= "  Constant k: "*string(k)*"\n" end
        output *= "----------------------------------------\n"
        if param == "mean"
            output *= "  SD: "*string(sd)*"\n"
        end
        if group == "one"
            output *= "  Null Value: "*string(a)*" "
            output *= "Test Value: "*string(b)*"\n"
        else
            output *= "  Group A Value: "*string(a)*" "
            output *= "Group B Value: "*string(b)*"\n"
        end
        if (type == "ei" || type == "ns")
            output *= "  Difference: "*string(diff)*"\n"
        end
        output *= "----------------------------------------\n"
        if group == "one"
            output *= "  Estimate: "*string(ceil(n))*"\n"
        else
            output *= "  Group A: "*string(ceil(n*k))*"  "
            output *= "Group B: "*string(ceil(n))*"\n"
            output *= "  Total: "*string(ceil(n)+ceil(n*k))*"\n"
        end
        #output *= "Estimate group A: "*string(ceil(n))*"\n"
        #output *= "Total: "*string(ceil(n))*"\n"
        output *= "----------------------------------------\n"
        if out == "str"
            return output
        elseif out == "print"
            print(output)
            return
        elseif out == "vstr"
            return n, output
        end
    end
    return n
end #sampleSize

#clinical trial power main function
function ctPower(;param="", type="", group="", alpha=0.05, diff=0, sd=0, a=0, b=0, n=0, k=1)
    if alpha >= 1 || alpha <= 0  return false end
    if (type == "ei" || type == "ns") && diff == 0 return false end
    if param == "prop" && !(group == "one" || group == "two" || type == "mcnm")  return false end
    if sd == 0 && param == "mean" return false end
    if k == 0 return false end
    if n == 0 return false end
    if param == "mean"
        if group == "one"
            if type == "ea"
                return oneSampleMeanEqualityP(a, b, sd, n, alpha=alpha)
            elseif type == "ei"
                return oneSampleMeanEquivalenceP(a, b, sd, diff, n, alpha=alpha)
            elseif type == "ns"
                return oneSampleMeanNSP(a, b, sd, diff, n, alpha=alpha)
            else return false end
        elseif group == "two"
            if type == "ea"
                return twoSampleMeanEqualityP(a, b, sd, n, alpha=alpha, k=k)
            elseif type == "ei"
                return twoSampleMeanEquivalenceP(a, b, sd, diff, n, alpha=alpha, k=k)
            elseif type == "ns"
                return twoSampleMeanNSP(a, b, sd, diff, n, alpha=alpha, k=k)
            else return false end
        else return false end
    elseif param == "prop"
        if 1 < a || a < 0 || 1 < b || b < 0 return false end
        if type == "mcnm"
            return mcnmP(a, b, n; alpha=alpha)
        else
            if group == "one"
                if type == "ea"
                    return oneProportionEqualityP(a, b, n; alpha=alpha)
                elseif type == "ei"
                    return oneProportionEquivalenceP(a, b, diff, n; alpha=alpha)
                elseif type == "ns"
                    return oneProportionNSP(a, b, diff, n; alpha=alpha)
                else return false end
            elseif group == "two"
                if type == "ea"
                    return twoProportionEqualityP(a, b, n; alpha=alpha, k=k)
                elseif type == "ei"
                    return twoProportionEquivalenceP(a, b, diff, n; alpha=alpha, k=k)
                elseif type == "ns"
                    return twoProportionNSP(a, b, diff, n; alpha=alpha, k=k)
                else return false end
            else return false end
        end
    elseif param == "or"
        if type == "ea"
            return orEqualityP(a, b, n; alpha=alpha, k=k)
        elseif type == "ei"
            return orEquivalenceP(a, b, diff, n; alpha=alpha, k=k, logdiff=true)
        elseif type == "ns"
            return orNSP(a, b, diff, n; alpha=alpha, k=k, logdiff=true)
        else return false end
    else return false end
end #ctPower
#-------------------------------------------------------------------------------
function beSampleN(;theta0=0.95, theta1=0.8, theta2=1.25, cv=0.0, alpha=0.05, beta=0.2, logscale::Bool=true, design::String="2x2", method::String="owenq", txt=true)
    theta0 = convert(Float64, theta0)
    theta1 = convert(Float64, theta1)
    theta2 = convert(Float64, theta2)
    cv = convert(Float64, cv)
    alpha = convert(Float64, alpha)
    beta = convert(Float64, beta)
    if cv <= 0 return false end
    if alpha >= 1 || alpha <= 0 || beta >= 1 || beta <= 0 return false end

    if logscale
        if theta1 < 0 || theta2 < 0 || theta0 < 0 return false end
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
    #values for approximate n
    td = (ltheta1 + ltheta2)/2
    rd = abs((ltheta1 - ltheta2)/2)

    if rd <= 0 return false end

    d0 = diffm - td

    #approximate n
    n0::Int = convert(Int, ceil(twoSampleMeanEquivalence(0, d0, se, rd, alpha=alpha, beta=beta)/2)*2)
    tp = 1 - beta  #target power

    if n0 > 5000 n0 = 5000 end
    pow = powerTOST(;alpha=alpha, logscale=false, theta1=ltheta1, theta2=ltheta2, theta0=diffm, cv=se, n=n0, design=design, method=method)
    np = 0
    powp::Float64 = 0.0
    if pow > tp
        while (pow > tp)
            np = n0
            powp = pow
            n0 = n0 - 2
            #pow = powerTOST(;alpha=alpha, logscale=false, theta1=ltheta1, theta2=ltheta2, theta0=diffm, cv=se, n=n0, design=design, method=method)
            pow = powerTOSTint(alpha, false, ltheta1, ltheta2, diffm, se, n0, design, method)
            if n0 < 4 return n0, pow end
        end
        estpower = powp
        estn     = np
    elseif pow < tp
        while (pow < tp)
            np = n0
            powp = pow
            n0 = n0 + 2
            #pow = powerTOST(;alpha=alpha, logscale=false, theta1=ltheta1, theta2=ltheta2, theta0=diffm, cv=se, n=n0, design=design, method=method)
            pow = powerTOSTint(alpha, false, ltheta1, ltheta2, diffm, se, n0, design, method)

            if n0 > 10000 return n0, pow end
        end
        estpower = pow
        estn     = n0
    else return n0, pow end
    return estn, estpower
end