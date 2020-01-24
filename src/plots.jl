# Clinical Trial Utilities
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)

struct PKPlotStyle
    linestyle::Symbol
    linecolor::Symbol
    markershape::Symbol
    markercolor::Symbol
end

const PKPLOTSTYLE = [
PKPlotStyle(:solid, :blue, :circle, :blue),
PKPlotStyle(:solid, :red, :utriangle, :red),
PKPlotStyle(:solid, :green, :diamond, :green),
PKPlotStyle(:solid, :yellow, :rect, :yellow),
PKPlotStyle(:solid, :gray, :xcross, :gray),
PKPlotStyle(:solid, :cyan, :cross, :cyan),
PKPlotStyle(:solid, :gold, :star, :gold),
PKPlotStyle(:solid, :magenta, :pentagon, :magenta),
PKPlotStyle(:solid, :purple, :heptagon, :purple),
PKPlotStyle(:solid, :indigo, :octagon, :indigo),
PKPlotStyle(:dot, :blue, :utriangle, :blue),
PKPlotStyle(:dot, :red, :circle, :red),
PKPlotStyle(:dot, :green, :rect, :green),
PKPlotStyle(:dot, :yellow, :diamond, :yellow),
PKPlotStyle(:dot, :gray, :cross, :gray),
PKPlotStyle(:dot, :cyan, :xcross, :cyan),
PKPlotStyle(:dot, :gold, :pentagon, :gold),
PKPlotStyle(:dot, :magenta, :star, :magenta),
PKPlotStyle(:dot, :purple, :octagon, :purple),
PKPlotStyle(:dot, :indigo, :heptagon, :indigo),
]

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
end
function plot!(subj::T) where T <: AbstractSubject
    return pkplot!(subj.time, subj.obs;)
end

function usort(data::DataSet{T}, list) where T <: AbstractData
    dl = Vector{Dict}(undef, 0)
    for i in data
        subd = Dict(k => i.sort[k] for k in list)
        if subd ∉ dl push!(dl, subd) end
    end
    dl
end

function plot(data::DataSet{T}; pagesort = nothing, typesort = nothing) where T <: AbstractSubject
    typedict = usort(data, typesort)
    pagedict = usort(data, pagesort)
    plots    = Vector{Any}(undef, 0)
    p        = nothing
    for i in pagedict
        fst = true
        for d in data
            if i ∈ d.sort
                #println(i , " --- ", d.sort)
                if fst
                    #p = pkplot(d.time, d.obs)
                    p = plot(d)
                    fst = false
                else
                    #pkplot!(p, d.time, d.obs)
                    plot!(d)
                end
            end
        end
        push!(plots, p)
    end
    plots
end
