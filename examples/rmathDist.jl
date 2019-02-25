# Clinical Trial Power and Sample Size calculation
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)
# Calculation based on Chow S, Shao J, Wang H. 2008. Sample Size Calculations in Clinical Research. 2nd Ed. Chapman & Hall/CRC Biostatistics Series.
# OwensQ function rewrited from https://github.com/Detlew/PowerTOST by Detlew Labes, Helmut Schuetz, Benjamin Lang
# Licence: GNU Affero General Public License v3.0

#This is comparison of Rmath and Distributions library approach
"""
pt(x, df)
TDIST=TDist(df)
cdf(TDIST, x)

pt(x,df,ncp)
NCT=NoncentralT(df,ncp)
cdf(NCT,x)


qnorm(1-0.05)
[1] 1.644854
quantile(Normal(), 1-0.05)
1.6448536269514717

pnorm(1)
[1] 0.8413447
cdf(Normal(),1)
0.841344746068543


> dnorm(1)
[1] 0.2419707
dnorm(1)
0.24197072451914337
pdf(Normal(), 1)
0.24197072451914337

qf(0.2,44,66)
quantile(FDist(44,66),0.2)

> pbinom(30,100,0.3)
[1] 0.5491236
julia> cdf(Binomial(100,0.3),30)
0.5491236007687905

> qbinom(0.05,100,0.3)
[1] 23
julia> quantile(Binomial(100,0.3),0.05)
23

> qchisq(0.95,1)
[1] 3.841459
quantile(Chisq(1), 0.95)
3.841458820694124

"""
