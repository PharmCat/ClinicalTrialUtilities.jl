# Clinical Trial Power and Sample Size calculation
# Version: 0.0.1
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat
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

"""
