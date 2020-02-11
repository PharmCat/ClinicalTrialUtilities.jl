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
PKPlotStyle(:solid, :magenta, :pentagon, :magenta),
PKPlotStyle(:solid, :purple, :heptagon, :purple),
PKPlotStyle(:solid, :indigo, :octagon, :indigo),
PKPlotStyle(:solid, :gold, :star, :gold),
PKPlotStyle(:solid, :yellow, :rect, :yellow),
PKPlotStyle(:solid, :gray, :xcross, :gray),
PKPlotStyle(:solid, :cyan, :cross, :cyan),
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

function randomstyle()
    linestyle   = ([:solid, :dash, :dot, :dashdot, :dashdotdot])[sample(1:5)]
    linecolor   = ([:blue, :red, :green, :yellow, :gray, :cyan, :gold, :magenta, :purple, :indigo])[sample(1:10)]
    markershape = ([:circle, :rect, :star5, :diamond, :hexagon, :cross, :xcross, :utriangle, :dtriangle, :rtriangle, :ltriangle, :pentagon, :heptagon, :octagon, :star4, :star6, :star7, :star8, :vline, :hline])[sample(1:20)]
    markercolor = ([:blue, :red, :green, :yellow, :gray, :cyan, :gold, :magenta, :purple, :indigo])[sample(1:10)]
    return PKPlotStyle(linestyle, linecolor, markershape, markercolor)
end

@userplot PKPlot

@recipe function f(subj::PKPlot)
    x, y = subj.args

    seriestype  --> :line
    xlabel      --> "Time"
    ylabel      --> "Concentration"
    link        --> :both
    legend      --> true
    grid        --> false
    #ticks       := [nothing :auto nothing]
    xlims       --> (minimum(x), maximum(x)),
    ylims       --> (0, maximum(y)*1.1)
    seriescolor --> :blue
    markershape --> :circle
    markersize  --> 3
    markercolor --> :match
    msalpha     --> 0
    (x, y)
end

function plotlabel(d)
    title = ""
    if length(d) > 0
        for (k, v) in d
            title *= "$(k) = $(v); "
        end
        title = title[1:length(title)]
    end
    return title
end

function plot(subj::T; title = nothing, legend = true, label = "AUTO", xlims = nothing, plotstyle::PKPlotStyle = PKPLOTSTYLE[1]) where T <: AbstractSubject

    if title === nothing
        title = plotlabel(subj.sort)
    end
    if xlims === nothing xlims = (minimum(subj.time), maximum(subj.time)) end
    p = pkplot(subj.time, subj.obs;
        title       = title,
        linestyle   = plotstyle.linestyle,
        linecolor   = plotstyle.linecolor,
        markershape = plotstyle.markershape,
        markercolor = plotstyle.markercolor,
        legend      = legend,
        label       = label,
        xlims       = xlims,
        )
    return p

end
function plot!(subj::T; legend = true, label = "AUTO", xlims = nothing,  plotstyle::PKPlotStyle = PKPLOTSTYLE[1]) where T <: AbstractSubject
    if xlims === nothing xlims = (minimum(subj.time), maximum(subj.time)) end
    p = pkplot!(subj.time, subj.obs;
        linestyle   = plotstyle.linestyle,
        linecolor   = plotstyle.linecolor,
        markershape = plotstyle.markershape,
        markercolor = plotstyle.markercolor,
        legend      = legend,
        label       = label,
        xlims       = xlims,
        )
    return p
end

function usort(data::DataSet{T}, list) where T <: AbstractData
    dl = Vector{Dict}(undef, 0)
    for i in data
        subd = Dict(k => i.sort[k] for k in list)
        if subd ∉ dl push!(dl, subd) end
    end
    dl
end

function pageplot(pagedatatmp, styledict, typesort; title = "", legend = true, xlims = nothing)
    utypes    = keys(styledict)
    labels    = Vector{String}(undef, 0)
    fst       = true
    p         = nothing
    for u in utypes
        for d in pagedatatmp
            if u ∈ d.sort
                label = "AUTO"
                if typesort !== nothing
                    tempsd = Dict(k => d.sort[k] for k in typesort)
                    style  = styledict[tempsd]
                    label  = plotlabel(tempsd)
                    if label ∉ labels
                        push!(labels, label)
                    else
                        label = ""
                    end
                end
                if fst
                    p = plot(d; title = title, legend = legend, xlims = xlims, plotstyle = style, label = label)
                    fst = false
                else
                    plot!(                  d; legend = legend, xlims = xlims, plotstyle = style, label = label)
                end
            end
        end
    end
    p
end
function pageplot(pagedatatmp; title = "", legend = true, xlims = nothing)
    fst       = true
    p         = nothing
    for d in pagedatatmp
        label = "AUTO"
        if fst
            p = plot(d; title = title, legend = legend, xlims = xlims, plotstyle = PKPLOTSTYLE[1], label = label)
            fst = false
        else
            plot!(                  d; legend = legend, xlims = xlims, plotstyle = PKPLOTSTYLE[1], label = label)
        end
    end
    p
end

function plot(data::DataSet{T}; title = nothing, legend = true, xlims = nothing, pagesort = nothing, typesort = nothing) where T <: AbstractSubject
    # Style types

    styledict  = nothing
    #style      = PKPLOTSTYLE[1]
    if pagesort !== nothing
        if !isa(pagesort, Array) pagesort = [pagesort] end
    end
    if typesort !== nothing
        styledict  = Dict()
        if !isa(typesort, Array) typesort = [typesort] end
        utypes = usort(data, typesort)
        for i = 1:length(utypes)
            if i <= 20
                styledict[utypes[i]] = PKPLOTSTYLE[i]
            else
                styledict[utypes[i]] = randomstyle()
            end
        end
    end
    autotitle = false
    if title === nothing autotitle = true end
    #---------------------------------------------------------------------------
    # PAGES
    plots    = Vector{Any}(undef, 0)
    p        = nothing
    if pagesort !== nothing
        pagedict   = usort(data, pagesort)
        for i in pagedict
            if autotitle
                title  = plotlabel(i)
            end
            pagedatatmp = Vector{eltype(data)}(undef, 0)
            for d in data
                if i ∈ d.sort
                    push!(pagedatatmp, d)
                end
            end
            if typesort !== nothing
                p = pageplot(pagedatatmp, styledict, typesort; title = title, legend = legend, xlims = xlims)
            else
                p = pageplot(pagedatatmp; title = title, legend = legend, xlims = xlims)
            end
            push!(plots, p)
        end
    else
        if typesort !== nothing
            p = pageplot(data, styledict, typesort; title = title, legend = legend, xlims = xlims)
        else
            p = pageplot(data; title = title, legend = legend, xlims = xlims)
        end
        push!(plots, p)
    end
    plots
end
