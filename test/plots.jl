# Clinical Trial Power and Sample Size calculation
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)

# This script contains some graphics examples
# Plots package require, not in dependencies


using Plots

# owensQ function
function graphOwensQ()
    m = Array{Float64}(undef, 120, 2)
    for a=1:120
        m[a,1] = a/20
        m[a,2] = ClinicalTrialUtilities.owensQ(7, 9, 0.5, 0, m[a,1])
    end
    plot(m[:,1],m[:,2],linewidth=2,title="owensQ", label=["Function"])
end
