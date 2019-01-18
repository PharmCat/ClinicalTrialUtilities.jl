module CTPSS
using Distributions

    function OneSampleMean(m0, m1, sd; alpha=0.05, beta=0.8)

        ZDIST = Normal()
        return (quantile(ZDIST, 1-alpha/2) + quantile(ZDIST, beta))^2*sd^2/(m0-m1)^2
    end

end # module
