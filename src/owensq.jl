# Clinical Trial Utilities
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)
# Calculation based on Chow S, Shao J, Wang H. 2008. Sample Size Calculations in Clinical Research. 2nd Ed. Chapman & Hall/CRC Biostatistics Series.
# OwensQ function rewrited from https://github.com/Detlew/PowerTOST by Detlew Labes, Helmut Schuetz, Benjamin Lang
# Licence: GNU Affero General Public License v3.0

# general OwensQ function
# https://github.com/Detlew/PowerTOST/blob/master/R/OwensQ.R#L52


#function owensqf(nu::Float64, t::Float64, delta::Float64, a::Float64, b::Float64)::Float64
#    sqrt(2pi)/gamma(nu/2.0)/2.0^((nu-2.0)/2.0)*quadgk(x -> cdf(Normal(), t * x / sqrt(nu) - delta) * x ^ (nu - 1.0) * pdf(Normal(), x), a, b, rtol=1.0E-8)[1]
#end

function owensq(nu::Real, t::Real, delta::Real, a::Real, b::Real)::AbstractFloat
    if a < 0 return  throw(CTUException(1011,"owensQ: a can not be < 0")) end
    if a==b return(0) end
    if a > b return throw(CTUException(1012,"owensQ: a can not be > b")) end
    if a > 0 return owensq(nu, t, delta, 0, b) - owensq(nu, t, delta, 0, a)  end #not effective - double integration
    if nu < 29 && abs(delta) > 37.62
        if isinf(b)
            return quadgk(x -> ifun1(x, nu, t, delta), 0, 1, rtol=1.0E-8)[1]
        else
            return owensqo(nu,t,delta,b)
        end
    else
        if isinf(b)
            #45 of OwensQ
            return cdf(NoncentralT(nu,delta),t)
        else
            integral = quadgk(x -> ifun1(x,nu,t,delta, b=b), 0, 1, rtol=1.0E-8)[1]
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
function owensqo(nu::Real, t::Real, delta::Real, b::Real; a::Real = 0.0)::AbstractFloat
    if nu < 1  throw(CTUException(1001,"owensQo: nu can not be < 1")) end
    if a != 0.0 throw(CTUException(1002,"owensQo: a can not be not 0")) end
    if isinf(b) return cdf(NoncentralT(nu,delta),t) end
    if isinf(delta) delta = sign(delta)*1e20 end
    A   = t/sqrt(nu)
    B   = nu/(nu + t*t)
    upr = nu-2
    #some code abs L vector H vector M vector
    av  = Array{Float64}(undef, nu)
    for i = 1:length(av)
        if i==1||i==2 av[i] = 1 else av[i] = 1 / ((i - 2) * av[i-1]) end
    end
    if (upr-1) > 0 ll = upr-1 else ll = 0 end
    L  = Array{Float64}(undef, ll)
    if isfinite(b)
         for i=1:length(L)
             if i==1 L[1] = 0.5 * A * B * b * pdf(ZDIST, b) * pdf(ZDIST, A * b - delta)
             else L[i] = av[i+3]*b*L[i-1] end
         end
    end
    if (upr+1)>0 ll = upr+1 else ll = 0 end
    H = Array{Float64}(undef, ll)
    if isfinite(b)
        for i = 1:length(H)
            if i==1 H[1] = -pdf(ZDIST, b)*cdf(ZDIST, A * b - delta)
            else H[i] = av[i+1] * b * H[i-1] end
        end
    end
    #pass1
    M   = Array{Float64}(undef, ll)
    sB  = sqrt(B)
    for i = 1:length(M)
        if i==1 M[1] = A * sB * pdf(ZDIST, delta * sB) * (cdf(ZDIST, delta * A * sB)-cdf(ZDIST, (delta * A * B - b) / sB)) end
        if i==2 M[2] = B * (delta * A * M[1] + A * pdf(ZDIST, delta * sB) * (pdf(ZDIST, delta * A * sB)-pdf(ZDIST,(delta * A * B - b) / sB))) end
        if i>2 M[i] = ((i - 2) / (i - 1)) * B * (av[i-1] * delta * A * M[i-1] + M[i-2]) - L[i-2] end
    end
    #pass
    sumt  = 0.0;
    if nu%2 > 0
        if upr>=1
            for i = 1:2:upr
                sumt = sumt + M[i + 1] + H[i + 1]
            end
        end
        qv = cdf(ZDIST, b) - 2 * owenst(b, A - delta / b) -
                        2 * owenst(delta * sB, (delta * A * B - b) / B / delta) +
                        2 * owenst(delta * sB, A) - (delta >= 0) + 2 * sumt
    else
        if upr>=0
            for i = 0:2:upr
                sumt = sumt + M[i+1] + H[i+1]
            end
        end
        qv = cdf(ZDIST,-delta)+sqrt(2*π)*sumt
    end
    #
end #OwensQo

#-------------------------------------------------------------------------------
# functions bellow used in owensQ
#PowerTost owensQ #37 and impl b for #52-54
@inline function ifun1(x::Real, nu::Real, t::Real, delta::Real; b=0.0::Real)::AbstractFloat
    return owenstint2(b + x / (1 - x), nu, t, delta) * 1 / (1 - x) ^ 2
end
#.Q.integrand
@inline function owenstint2(x::Real, nu::Real, t::Real, delta::Real)::AbstractFloat
    if x == 0 return 0 end
    return sign(x) ^ (nu - 1.0) * cdf(ZDIST, t * x / sqrt(nu) - delta) *
    exp((nu - 1.0) * log(abs(x)) - 0.5 * x ^ 2 -
    ((nu * 0.5) - 1.0) * LOG2 - lgamma(nu * 0.5))
end
#-------------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Owen's T-function according to algorithm AS76 and remarks AS R65, AS R80
# R port of the FORTRAN code in the References and matlab code given on
# https://people.sc.fsu.edu/~jburkardt/m_src/asa076/asa076.html
# by J. Burkhardt, license GNU LGPL
# rewrite from PowerTOST
function owenst(h::Real, a::Real)::AbstractFloat
    #not implemented
    epsilon = eps()
    # special cases
    # D.B. Owen
    # "Tables for computing bivariate normal Probabilities"
    # The Annals of Mathematical Statistics, Vol. 27 (4) Dec. 1956, pp. 1075-1090
    if abs(a) < epsilon || !isfinite(h) || abs(1-abs(a)) < epsilon || abs(h) < epsilon || !isfinite(abs(a))
        if abs(a)< epsilon return 0 end
        if !isfinite(h) return 0 end
        if abs(1-abs(a)) < epsilon return sign(a) * 0.5 * cdf(ZDIST, h) * (1 - cdf(ZDIST, h)) end
        if abs(h) < epsilon return atan(a) * 0.5 / pi end
        if !isfinite(abs(a))
            if h < 0 tha=cdf(ZDIST, h) * 0.5  else tha = (1 - cdf(ZDIST, h)) * 0.5 end
            return sign(a) * tha
        end
    end
    aa = abs(a)
    if aa <= 1.0
        return tfn(h, a)
    else
        ah   = aa * h
        gh   = cdf(ZDIST, h)
        gah  = cdf(ZDIST, ah)
        tha  = 0.5 * (gh + gah) - gh * gah - tfn(ah, 1.0 / aa)
    end
    if a < 0.0 return -tha else return tha end
end #OwensT(h,a)

# Julia implemetation from PowerTOST
# Auxillary function
# R port of the matlab code given on
# https://people.sc.fsu.edu/~jburkardt/m_src/owens/owens.html
# https://people.sc.fsu.edu/~jburkardt/c_src/asa076/asa076.c
# https://gist.github.com/mtrbean/464859cdf09b260d5ea6
# by J. Burkhardt license GNU LGPL
# is called as tfn(h, a) if a<=1
# otherwise as tfn(a*h, 1/a)
@inline function tfn(x, fx)
    ng  = 5
    r   = [0.14776211235737646, 0.13463335965499826, 0.10954318125799103, 0.07472567457529028, 0.033335672154344104]
    u   = [0.07443716949081561, 0.2166976970646236, 0.3397047841495122, 0.43253168334449227, 0.48695326425858587]
    tv1 = eps();
    tv2  = tv3 = 15.0
    tv4  = 1.0E-05
    if tv2 < abs(x) return 0 end
    xs   = -0.5 * x * x
    x2   = fx
    fxs  = fx * fx
    if tv3 <= (log(1.0 + fxs) - xs * fxs)
        x1  = 0.5 * fx
        fxs = 0.25 * fxs
        while true
            rt  = fxs + 1.0
            x2  = x1 + (xs * fxs + tv3 - log(rt)) / (2.0 * x1 * (1.0 / rt - xs))
            fxs = x2 * x2
            if abs(x2 - x1) < tv4 break end
            x1 = x2
        end
    end
        #Vectorized:
    r1 = @. 1.0 + fxs * (0.5 + u) ^ 2
    r2 = @. 1.0 + fxs * (0.5 - u) ^ 2
    rt = @. r * (exp(xs * r1) / r1 + exp(xs * r2) / r2)

    return sum(rt) * x2 * PI2INV
end #tfn(x, fx)
