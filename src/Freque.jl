

function freque(data::DataFrame; vars::Symbol, alpha = 0.05)::DataFrame
    result = DataFrame(value = Any[], n = Int[], p = Float64[], cil = Float64[], ciu = Float64[])
    list = unique(data[vars])

    n = length(data[vars])
    for i in list
        ne = count(x -> (x == i), data[vars])
        pe = ne/n
        ci = ClinicalTrialUtilities.CI.oneProp(ne, n, alpha=alpha, method=:wald)
        push!(result, [i, ne, pe, ci.lower, ci.upper])
    end
    return result
end
