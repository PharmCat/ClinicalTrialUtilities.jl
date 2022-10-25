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

function randomstyle(rng)
    linestyle   = ([:solid, :dash, :dot, :dashdot, :dashdotdot])[sample(rng, 1:5)]
    linecolor   = ([:blue, :red, :green, :yellow, :gray, :cyan, :gold, :magenta, :purple, :indigo])[sample(rng, 1:10)]
    markershape = ([:circle, :rect, :star5, :diamond, :hexagon, :cross, :xcross, :utriangle, :dtriangle, :rtriangle, :ltriangle, :pentagon, :heptagon, :octagon, :star4, :star6, :star7, :star8, :vline, :hline])[sample(rng, 1:20)]
    markercolor = ([:blue, :red, :green, :yellow, :gray, :cyan, :gold, :magenta, :purple, :indigo])[sample(rng, 1:10)]
    return PKPlotStyle(linestyle, linecolor, markershape, markercolor)
end

@userplot PKPlot
@userplot PKElimpPlot

function luceil(x)
    fl = Int(floor(log10(x)))
    if fl < 0 fl = 0 end
    ceil(x/10^fl)*10^fl
end

@recipe function f(subj::PKPlot; lcd = 5, tcd = 6)
    x, y = subj.args

    if isa(lcd, Real)
        lc = luceil(maximum(x->isnan(x) ? -Inf : x, y)/lcd)
        yt = 0:lc:lc*lcd
    elseif isa(lcd, StepRange)
        yt = lcd
    elseif lcd == :all
        yt = y
    end

    if isa(tcd, Real)
        tc = luceil(maximum(x)/tcd)
        xt = 0:tc:tc*tcd
    elseif isa(tcd, StepRange)
        xt = tcd
    elseif tcd == :all
        xt = x
    end

    widen             --> true
    seriestype        --> :line
    xguide            --> "Time"
    link              --> :both
    legend            --> true
    grid              --> true
    gridstyle         --> :auto
    #ticks       := [nothing :auto nothing]
    #xlims             --> (minimum(x), maximum(x)*1.1)
    #ylims             --> (0, maximum(y)*1.1)

    if !isa(lcd, Symbol) || lcd != :auto
        yticks            --> yt
    end
    if !isa(tcd, Symbol) || tcd != :auto
        xticks            --> xt
    end

    seriescolor       --> :blue
    markershape       --> :circle
    markersize        --> 3
    markercolor       --> :match
    markerstrokealpha --> 0
    (x, y)
end

@recipe function f(subj::PKElimpPlot)
    x, y = subj.args
    seriestype        --> :line
    legend            --> false
    markersize        --> 0
    markerstrokealpha --> 0
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

"""
    pkplot(subj::T; plotstyle::PKPlotStyle = PKPLOTSTYLE[1], kwargs...) where T <: AbstractSubject

Plot for subject

* subj - subject;
* plotstyle - styles for plots.

"""
function pkplot(subj::T; plotstyle::PKPlotStyle = PKPLOTSTYLE[1], ls = false, elim = false, xticksn = 6, yticksn = 5, kwargs...) where T <: AbstractSubject
    time = subj.time
    obs  = subj.obs
    kwargs = Dict{Symbol, Any}(kwargs)
    k = keys(kwargs)

    if !(:title in k)
        kwargs[:title] = plotlabel(subj.sort)
    end
    if !(:xlims in k)
        kwargs[:xlims] = (minimum(subj.time), maximum(subj.time))
    end
    if !(:legend in k)
        kwargs[:legend] = true
    end
    if !(:ylabel in k)
        kwargs[:ylabel] = "Concentration"
    end
    if :yscale in k
        if kwargs[:yscale] in [:ln, :log, :log2, :log10] ls = true end
    end
    if ls == true
        inds = findall(x->x>0, subj.obs)
        time = subj.time[inds]
        obs = log.(subj.obs[inds])
        if (:ylims in k)
            kwargs[:ylims] = (0,log(kwargs[:ylims][2]))
        end
    end

    p = pkplot(time, obs;  lcd = yticksn, tcd = xticksn,
        linestyle   = plotstyle.linestyle,
        linecolor   = plotstyle.linecolor,
        markershape = plotstyle.markershape,
        markercolor = plotstyle.markercolor,
        kwargs...
        )
    if elim
        if length(subj.keldata) > 0
            rsq, rsqn = findmax(subj.keldata.ar)
            lz        = subj.keldata.a[rsqn]
            lzint     = subj.keldata.b[rsqn]
            ts        = subj.time[subj.kelrange.kelstart]
            te        = subj.time[subj.kelrange.kelend]

            if ls true
                x = [ts,te]
                y = lzint .+ lz .* x
            else
                x = collect(ts:(te-ts)/100:te)
                y = exp.(lzint .+ lz .* x)
            end
            pkelimpplot!(p, x, y)
        end
    end
    return p

end

function pkplot!(subj::T; plotstyle::PKPlotStyle = PKPLOTSTYLE[1], ls = false, elim = false, xticksn = :auto, yticksn = :auto, kwargs...) where T <: AbstractSubject
    time = subj.time
    obs = subj.obs
    kwargs = Dict{Symbol, Any}(kwargs)
    k = keys(kwargs)

    if !(:xlims in k)
        kwargs[:xlims] = (minimum(subj.time), maximum(subj.time))
    end
    if !(:legend in k)
        kwargs[:legend] = true
    end
    if :yscale in k
        if kwargs[:yscale] in [:ln, :log, :log2, :log10] ls = true end
    end
    if ls == true
        inds = findall(x->x>0, subj.obs)
        time = subj.time[inds]
        obs = log.(subj.obs[inds])
        if (:ylims in k)
            kwargs[:ylims] = (0,log(kwargs[:ylims][2]))
        end
    end


    p = pkplot!(time, obs;  lcd = yticksn, tcd = xticksn,
        linestyle   = plotstyle.linestyle,
        linecolor   = plotstyle.linecolor,
        markershape = plotstyle.markershape,
        markercolor = plotstyle.markercolor,
        kwargs...
        )
    if :blline in k
        if kwargs[:blline] && isa(subj, PDSubject)
            hline!(p, [subj.bl], color = :red, label = "bl")
        end
    end
    if :thline in k
        if kwargs[:thline] && isa(subj, PDSubject)
            hline!(p, [subj.th], color = :blue, label = "th")
        end
    end
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

function pageplot(pagedatatmp, styledict, typesort, utypes; xticksn = 6, yticksn = 5, kwargs...)
    kwargs = Dict{Symbol, Any}(kwargs)
    k = keys(kwargs)

    if !(:title in k)
        kwargs[:title] = ""
    end
    if !(:xlims in k)
        kwargs[:xlims] = (minimum(subj.time), maximum(subj.time))
    end
    if !(:legend in k)
        kwargs[:legend] = true
    end
    if !(:ylabel in k)
        kwargs[:ylabel] = "Concentration"
    end

    #utypes    = keys(styledict)

    labels    = Vector{String}(undef, 0)
    fst       = true
    p         = nothing
    ylm = ylm2 = 0.0
    for u in utypes
        for d in pagedatatmp
            if u ⊆ d.sort
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

                kwargs[:label] = label
                if fst
                    ylm = maximum(x->isnan(x) ? -Inf : x, d.obs)
                    #MNIMUM
                    p = pkplot(d; plotstyle = style, xticksn = xticksn, yticksn = yticksn, ylims = (0, ylm*1.2), kwargs...)
                    fst = false
                else
                    ylm2 = maximum(x->isnan(x) ? -Inf : x, d.obs)
                    #MNIMUM
                    if ylm2 > ylm ylm = ylm2 end
                    pkplot!(   d; plotstyle = style, xticksn = xticksn, yticksn = yticksn, ylims = (0, ylm*1.2), kwargs...)
                end
            end
        end
    end
    if :vline in k
        vline!(p, kwargs[:vline])
    end
    if :hline in k
        hline!(p, kwargs[:hline])
    end
    p
end
function pageplot(pagedatatmp; elim = false,  xticksn = 6, yticksn = 5,  kwargs...)
    kwargs = Dict{Symbol, Any}(kwargs)
    k = keys(kwargs)

    if !(:title in k)
        kwargs[:title] = ""
    end
    if !(:xlims in k)
        kwargs[:xlims] = (minimum(subj.time), maximum(subj.time))
    end
    if !(:legend in k)
        kwargs[:legend] = true
    end
    if !(:ylabel in k)
        kwargs[:ylabel] = "Concentration"
    end

    fst       = true
    p         = nothing
    if length(pagedatatmp) > 1 elim = false end
    ylm = ylm2 = 0.0
    for d in pagedatatmp
        kwargs[:label] = "AUTO"
        if fst
            ylm = maximum(x->isnan(x) ? -Inf : x, d.obs)
            #MNIMUM
            p = pkplot(d; plotstyle = PKPLOTSTYLE[1], elim = elim, xticksn = xticksn, yticksn = yticksn, ylims = (0, ylm*1.2), kwargs...)
            fst = false
        else
            ylm2 = maximum(x->isnan(x) ? -Inf : x, d.obs)
            #MNIMUM
            if ylm2 > ylm ylm = ylm2 end
            pkplot!(   d; plotstyle = PKPLOTSTYLE[1], xticksn = xticksn, yticksn = yticksn, ylims = (0, ylm*1.2), kwargs...)
        end
    end
    if :vline in k
        vline!(p, kwargs[:vline])
    end
    if :hline in k
        hline!(p, kwargs[:hline])
    end
    p
end

"""
    pkplot(data::DataSet{T};  pagesort::Union{Nothing, Vector{Symbol}} = nothing, typesort::Union{Nothing, Vector{Symbol}} = nothing, kwargs...) where T <: AbstractSubject

Plot for subjects in dataset.

* data - subjects dataset;
* pagesort - subject page groupping;
* typesort - subject sorting within page;
"""
function pkplot(data::DataSet{T};  pagesort::Union{Nothing, Vector{Symbol}} = nothing, typesort::Union{Nothing, Vector{Symbol}} = nothing, elim = false, xticksn = :auto, yticksn = :auto, kwargs...) where T <: AbstractSubject
    # Style types
    kwargs = Dict{Symbol, Any}(kwargs)
    k = keys(kwargs)

    if !(:title in k)
        kwargs[:title] = :AUTO
    end
    if !(:xlims in k)
        kwargs[:xlims] = nothing
    end
    if !(:legend in k)
        kwargs[:legend] = true
    end
    if !(:ylabel in k)
        kwargs[:ylabel] = "Concentration"
    end

    styledict  = nothing
    local utypes
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
                styledict[utypes[i]] = randomstyle(MersenneTwister(34534))
            end
        end
    end
    autotitle = false
    if kwargs[:title] == :AUTO autotitle = true end
    #---------------------------------------------------------------------------
    # PAGES
    plots    = Vector{Any}(undef, 0)
    p        = nothing
    if pagesort !== nothing
        pagedict   = usort(data, pagesort)
        for i in pagedict
            if autotitle
                kwargs[:title]  = plotlabel(i)
            end
            pagedatatmp = Vector{eltype(data)}(undef, 0)
            for d in data
                if i ⊆ d.sort
                    push!(pagedatatmp, d)
                end
            end
            if typesort !== nothing
                p = pageplot(pagedatatmp, styledict, typesort, utypes; xticksn = xticksn, yticksn = yticksn,  kwargs...)
            else
                p = pageplot(pagedatatmp; elim = elim,  xticksn = xticksn, yticksn = yticksn, kwargs...)
            end
            push!(plots, p)
        end
    else
        if kwargs[:title] == :AUTO kwargs[:title] = "" end
        if typesort !== nothing
            p = pageplot(data, styledict, typesort, utypes; xticksn = xticksn, yticksn = yticksn,  kwargs...)
        else
            p = pageplot(data; elim = elim, xticksn = xticksn, yticksn = yticksn,  kwargs...)
        end
        return p
    end
    plots
end
