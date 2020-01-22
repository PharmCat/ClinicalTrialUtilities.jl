# Clinical Trial Utilities
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)
#=
struct MyType
    time
    obs
end

@recipe function f(::Type{MyType}, mw::MyType)
    x = mw.time
    y = mw.obs
end
=#

@userplot PKPlot

@recipe function f(subj::PKPlot)
    x, y = subj.args

    seriestype  --> :line
    xlabel      --> "Time"
    ylabel      --> "Concentration"
    link        --> :both
    legend      --> false
    grid        --> false
    #ticks       := [nothing :auto nothing]
    xlims       --> (minimum(x), maximum(x)),
    ylims       --> (0, maximum(y)*1.1)
    seriescolor --> :blue
    markershape --> :circle
    markersize  --> 2
    markercolor --> :match
    msalpha     --> 0
    (x, y)
end

function plot(subj::T; title = nothing) where T <: AbstractSubject

    if title === nothing
        title = ""
        if length(subj.sort) > 0
            for (k, v) in subj.sort
                title *= "$(k) = $(v); "
            end
            title = title[1:end-2]
        end
        return pkplot(subj.time, subj.obs;
            title = title,
            )
        end
function plot!(subj::T) where T <: AbstractSubject
    return pkplot!(subj.time, subj.obs;)
end

function plot(data::DataSet{T}) where T <: AbstractSubject

end
