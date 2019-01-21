
#OwensQOwen PowerTost
# ----------------------------------------------------------------------------
# Owen's Q-function
# Calculates Owen's Q-function via repeated integration by parts
# formulas as given in
# Owen, D B (1965)
# "A Special Case of a Bivariate Non-central t-Distribution"
# Biometrika Vol. 52, pp.437-446.
# Port from PowerTOST - dlabes Mar 2012
function OwensQo(nu,t,delta,b;a=0)
    if nu < 1 return false end
    if a != 0 return false end
    if isinf(b) return cdf(NoncentralT(nu,delta),t) end
    if isinf(delta) delta = sign(delta)*1e20 end

    A = t/sqrt(nu)
    B = nu/(nu + t*t)
    upr = nu-2

    #some code abs L vector H vector M vector
    av = Array{Float64}(undef, nu)

    println("hhh")
    for i = 1:length(av)
        #println(i)
        if i==1||i==2 av[i] = 1 else av[i] = 1/((i-2)*av[i-1]) end
    end

    if (upr-1)>0 ll = upr-1 else ll = 0 end

    L  = Array{Float64}(undef, ll)

    if isfinite(b)
         for i=1:length(L)
             if (i==1) L[1] = A*B*b*dnorm(b)*dnorm(A*b-delta)
             else L[i] = av[i+3]*b*L[k-1] end
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

    M = Array{Float64}(undef, ll)
    sB = sqrt(B)
    for i = 1:length(M)
        if i==1 M[1] = A*sB*dnorm(delta*sB)*(pnorm(delta*A*sB)-pnorm((delta*A*B-b)/sB)) end
        if i==2 M[2] = B*(delta*A*M[1]+A*dnorm(delta*sB)*(dnorm(delta*A*sB)-dnorm((delta*A*B-b)/sB))) end
        if i>2 M[i] = ((i-2)/(i-1))*(av[i-1]*delta*A*M[i-1]+M[i-2]) - L[i-2] end
    end

    sumt = 0
    if nu%2 > 0
        if upr>=1
            for i = 1:upr:2
                sumt = sumt + M[i+1] + M[i+1]
            end
        end

        qv = pnorm(b) - 2*OwensT(b,A-delta/b)-
                        2*OwensT(delta*sB, (delta*A*B-b)/B/delta)+
                        2*OwensT(delta*sB, A) - (delta>=0) + 2*sumt

    else
        if upr>=0
            for i = 0:upr:2
                sumt = sumt + M[i+1] + M[i+1]
            end
        end
        qv = pnorm(-delta)+sqrt(2*π)*sumt
    end




end #OwensQo

#-------------------------------------------------------------------------------
#                        DEPRECATED
#OT_integrand
function OwensTint(x, h)
    return exp(-0.5*h^2*(1+x^2))/(1+x^2)/2/π
end
#OwensT_old
function OwensTo(h, a)
    return quadgk(x -> OwensTint(x, h),0, a)[1]
end
#-------------------------------------------------------------------------------


#.Q.integrand
function OwensTint2(x, nu, t, delta)
    Qconst = -(nu/2-1)*log(2)-lgamma(nu/2)
    return sign(x)^(nu-1)*pnorm(t*x/sqrt(nu)-delta)*exp((nu-1)*log(abs(x))-0.5*x^2+Qconst)
end
#PowerTost OwensQ #37 and impl b for #52-54
function ifun1(x, nu, t, delta; b=0)
    return OwensTint2(b+x/(1-x), nu, t, delta)/(1-y)^2
end

# ----------------------------------------------------------------------------
# Owen's T-function according to algorithm AS76 and remarks AS R65, AS R80
# R port of the FORTRAN code in the References and matlab code given on
# https://people.sc.fsu.edu/~jburkardt/m_src/asa076/asa076.html
# by J. Burkhardt, license GNU LGPL
#
# no trouble with integrate()
# arguments must be scalars!
# rewrite from PowerTOST
function OwensT(h,a)
    #not implemented
end
