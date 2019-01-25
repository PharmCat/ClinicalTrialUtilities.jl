# Clinical Trial Power and Sample Size calculation
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat
# Calculation based on Chow S, Shao J, Wang H. 2008. Sample Size Calculations in Clinical Research. 2nd Ed. Chapman & Hall/CRC Biostatistics Series.
# OwensQ function rewrited from https://github.com/Detlew/PowerTOST by Detlew Labes, Helmut Schuetz, Benjamin Lang
# Licence: GNU Affero General Public License v3.0

# general OwensQ function
# https://github.com/Detlew/PowerTOST/blob/master/R/OwensQ.R#L52

function owensQ(nu, t, delta, a, b)
    if a==b return(0) end
    if nu < 29 && abs(delta) > 37.62
        if isinf(b)
            return quadgk(x -> ifun1(x, nu, t, delta), 0, 1, rtol=1.0E-8)[1] #ifun1 described in OwensQ.jl
        else
            return owensQo(nu,t,delta,b) #described in OwensQ.jl
        end
    else
        if isinf(b)
            #45 of OwensQ
            return cdf(NoncentralT(nu,delta),t)
        else
            integral = quadgk(x -> ifun1(x,nu,t,delta, b=b), 0, 1)[1]
            #58 of OwensQ
            return cdf(NoncentralT(nu,delta),t)-integral
        end
    end
end #owensQ

# ----------------------------------------------------------------------------
# Owen's Q-function
# Calculates Owen's Q-function via repeated integration by parts
# formulas as given in
# Owen, D B (1965)
# "A Special Case of a Bivariate Non-central t-Distribution"
# Biometrika Vol. 52, pp.437-446.
# Port from PowerTOST - dlabes Mar 2012

#OwensQOwen
function owensQo(nu,t,delta,b;a=0)
    if nu < 1 return false end
    if a != 0 return false end
    if isinf(b) return cdf(NoncentralT(nu,delta),t) end
    if isinf(delta) delta = sign(delta)*1e20 end
    A = t/sqrt(nu)
    B = nu/(nu + t*t)
    upr = nu-2

    #some code abs L vector H vector M vector
    av = Array{Float64}(undef, nu)
    for i = 1:length(av)
        if i==1||i==2 av[i] = 1 else av[i] = 1/((i-2)*av[i-1]) end
    end

    if (upr-1) > 0 ll = upr-1 else ll = 0 end
    L  = Array{Float64}(undef, ll)
    if isfinite(b)
         for i=1:length(L)
             if i==1 L[1] = 0.5*A*B*b*dnorm(b)*dnorm(A*b-delta)
             else L[i] = av[i+3]*b*L[i-1] end
         end
    end

    if (upr+1)>0 ll = upr+1 else ll = 0 end
    H = Array{Float64}(undef, ll)
    if isfinite(b)
        for i = 1:length(H)
            if i==1 H[1] = -dnorm(b)*pnorm(A*b-delta)
            else H[i] = av[i+1]*b*H[i-1] end
        end
    end
    #pass1
    M = Array{Float64}(undef, ll)
    sB = sqrt(B)

    for i = 1:length(M)
        if i==1 M[1] = A*sB*dnorm(delta*sB)*(pnorm(delta*A*sB)-pnorm((delta*A*B-b)/sB)) end
        if i==2 M[2] = B*(delta*A*M[1]+A*dnorm(delta*sB)*(dnorm(delta*A*sB)-dnorm((delta*A*B-b)/sB))) end
        if i>2 M[i] = ((i-2)/(i-1))*B*(av[i-1]*delta*A*M[i-1]+M[i-2]) - L[i-2] end
    end
    #pass
    sumt = 0
    if nu%2 > 0
        if upr>=1
            for i = 1:2:upr
                sumt = sumt + M[i+1] + H[i+1]
            end
        end
        qv = pnorm(b) - 2*owensT(b,A-delta/b) -
                        2*owensT(delta*sB, (delta*A*B-b)/B/delta) +
                        2*owensT(delta*sB, A) - (delta>=0) + 2*sumt
    else
        if upr>=0
            for i = 0:2:upr
                sumt = sumt + M[i+1] + H[i+1]
            end
        end
        qv = pnorm(-delta)+sqrt(2*π)*sumt
    end
    #
end #OwensQo

#-------------------------------------------------------------------------------
# functions bellow used in OwensQ in CTPSS.jl
#PowerTost OwensQ #37 and impl b for #52-54
function ifun1(x, nu, t, delta; b=0)
    return owensTint2(b+x/(1-x), nu, t, delta)*1/(1-x)^2
end
#.Q.integrand
function owensTint2(x, nu, t, delta)
    if x == 0 return 0 end
    Qconst = -((nu/2.0)-1.0)*log(2.0)-lgamma(nu/2.0)
    return sign(x)^(nu-1)*pnorm(t*x/sqrt(nu)-delta)*exp((nu-1)*log(abs(x))-0.5*x^2+Qconst)
end
#-------------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Owen's T-function according to algorithm AS76 and remarks AS R65, AS R80
# R port of the FORTRAN code in the References and matlab code given on
# https://people.sc.fsu.edu/~jburkardt/m_src/asa076/asa076.html
# by J. Burkhardt, license GNU LGPL
#
# no trouble with integrate()
# arguments must be scalars!
# rewrite from PowerTOST
function owensT(h,a)
    #not implemented
    epsilon = eps()
    # special cases
    # D.B. Owen
    # "Tables for computing bivariate normal Probabilities"
    # The Annals of Mathematical Statistics, Vol. 27 (4) Dec. 1956, pp. 1075-1090
    if abs(a) < epsilon || !isfinite(h) || abs(1-abs(a)) < epsilon || abs(h) < epsilon || !isfinite(abs(a))
        if abs(a)< epsilon return 0 end
        if !isfinite(h) return 0 end
        if abs(1-abs(a)) < epsilon return sign(a)*0.5*pnorm(h)*(1-pnorm(h)) end
        if abs(h)< epsilon return atan(a)/2/pi end
        if !isfinite(abs(a))
            if h<0 tha=pnorm(h)/2 else tha = (1-pnorm(h))/2 end
            return sign(a)*tha
        end
    end
    aa = abs(a)
    if aa <= 1.0
        tha = tfn(h, a)
        return tha
    else
        ah  = aa * h
        gh  = pnorm(h)
        gah = pnorm(ah)
        tha = 0.5*(gh + gah) - gh*gah - tfn(ah, 1.0/aa)
    end
    if a < 0.0 return -tha else return tha end
end #OwensT(h,a)

# Julia implemetation from PowerTOST
# Auxillary function
# R port of the matlab code given on
# https://people.sc.fsu.edu/~jburkardt/m_src/owens/owens.html
# by J. Burkhardt license GNU LGPL
#
# is called as tfn(h, a) if a<=1
# otherwise as tfn(a*h, 1/a)
function tfn(x, fx)
    ng  = 5
    r   = [0.1477621, 0.1346334, 0.1095432, 0.0747257, 0.0333357]
    u   = [0.0744372, 0.2166977, 0.3397048, 0.4325317, 0.4869533]
    tp  = 1/2/pi
    tv1 = eps();
    tv2 = 15.0
    tv3 = 15.0
    tv4 = 1.0E-05
    if tv2 < abs(x) return 0 end
    xs  = -0.5*x*x
    x2  = fx
    fxs = fx * fx
    if tv3 <= (log(1.0 + fxs) - xs*fxs)
        x1  = 0.5 * fx
        fxs = 0.25 * fxs
        while true
            rt  = fxs + 1.0
            x2  = x1 + (xs*fxs+tv3-log(rt))/(2.0*x1*(1.0/rt-xs))
            fxs = x2 * x2
            if abs(x2 - x1) < tv4 break end
            x1 = x2
        end
    end
        # 10 point Gaussian quadrature.
        # original via loop
        rt = 0
        for i in 1:ng
            r1 = 1.0 + fxs*(0.5 + u[i])^2
            r2 = 1.0 + fxs*(0.5 - u[i])^2
            rt = rt + r[i]*(exp(xs*r1)/r1 + exp(xs*r2)/r2)
        end
    return rt*x2*tp
end #tfn(x, fx)

"""
#-------------------------------------------------------------------------------
#                        DEPRECATED
#OT_integrand
function owensTint(x, h)
    return exp(-0.5*h^2*(1+x^2))/(1+x^2)/2/π
end
#OwensT_old
function owensTo(h, a)
    return quadgk(x -> owensTint(x, h),0, a)[1]
end
#-------------------------------------------------------------------------------
"""
