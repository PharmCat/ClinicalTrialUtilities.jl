#deprecated


#=
function ncarule!(data::DataFrame, conc::Symbol, time::Symbol, rule::LimitRule)
        sort!(data, [time])

        cmax, tmax, tmaxn = ctmax(data, conc, time)
        #NaN Rule
        if rule.nan !== NaN
            for i = 1:nrow(data)
                if data[i, conc] === NaN data[i, conc] = rule.nan end
            end
        end

        #LLOQ rule
        for i = 1:nrow(data)
            if data[i, conc] < rule.lloq
                if i <= tmaxn
                    data[i, conc] = rule.btmax
                else
                    data[i, conc] = rule.atmax
                end
            end
        end

        #NaN Remove rule
        if rule.rm
            filterv = Array{Int, 1}(undef, 0)
            for i = 1:nrow(data)
                if data[i, conc] === NaN push!(filterv, i) end
            end
            if length(filterv) > 0
                deleterows!(data, filterv)
            end
        end
end
function ncarule!(data::PKSubject, rule::LimitRule)
    #!!!!
end
=#

#=
function ctmax(data::PKSubject, dosetime, tau)
    if tau === NaN || tau <= 0 return ctmax(data, dosetime) end
    tautime = dosetime + tau
    s, e = ncarange(data, dosetime, tau)
    if dosetime == data.time[1] && tautime == data.time[end] return ctmax(data.time, data.obs) end
    mask           = trues(length(data))
    scpredict      = linpredict(data.time[s], data.time[s+1], data.dosetime.time, data.obs[s], data.obs[s+1])
    mask[1:s]     .= false

    if e < length(data)
        ecpredict      = linpredict(data.time[e], data.time[e+1], tautime, data.obs[e], data.obs[e+1])
        mask[e+1:end] .= false
    end
    cmax, tmax, tmaxn = ctmax(data.time[mask], data.obs[mask])
    if e < length(data)
        if cmax > max(scpredict, ecpredict)
            return cmax, tmax, tmaxn + s
        elseif scpredict >= max(cmax, ecpredict)
            return scpredict, data.dosetime.time, s
        else
            return ecpredict, tautime, e
        end
    elseif cmax > scpredict return max, tmax, tmaxn + s else return scpredict, data.dosetime.time, s end
end
=#

#=
function designinfo(d::OneGroup)::String
        return "One group design"
end
function designinfo(d::Parallel)::String
        return "Two parallel groups design"
end
function designinfo(d::Crossover)::String
        return "Crossover design"
end
=#

#=
function cmh(data::DataFrame; a = :a, b = :b, c = :c, d = :d, alpha = 0.05, type = :diff, method = :default, logscale = false)::ConfInt
    n1 = data[:, a] + data[:, b]
    n2 = data[:, c] + data[:, d]
    N  = data[:, a] + data[:, b] + data[:, c] + data[:, d]
    z = quantile(ZDIST, 1 - alpha/2)
    if type == :diff
        if method == :default method = :sato end #default method

        est = sum(data[:, a] .* (n2 ./ N) - data[:, c] .* (n1 ./ N)) / sum(n1 .* (n2 ./ N))
        #1
        if method == :sato
            se   = sqrt((est * (sum(data[:, c] .* (n1 ./ N) .^ 2 - data[:, a] .* (n2 ./ N) .^2 + (n1 ./ N) .* (n2 ./ N) .* (n2-n1) ./ 2)) + sum(data[:, a] .* (n2 - data[:, c]) ./ N + data[:, c] .* (n1 - data[:, a]) ./ N)/2) / sum(n1 .* (n2 ./ N)) .^ 2) # equation in: Sato, Greenland, & Robins (1989)
        #2
        elseif method == :gr
            se   = sqrt(sum(((data[:, a] ./N .^2) .* data[:, b] .* (n2 .^2 ./ n1) + (data[:, c] ./N .^2) .* data[:, d] .* (n1 .^2 ./ n2))) / sum(n1 .*(n2 ./ N))^2) # equation in: Greenland & Robins (1985)
        end
        #zval = est/se
        #pval = 2*(1-cdf(Normal(), abs(zval)))
        return ConfInt(est - z*se, est + z*se, est, alpha)
    elseif type == :or
        #...
        Pi = (data[:, a] ./ N) + (data[:, d] ./ N)
        Qi = (data[:, b] ./ N) + (data[:, c] ./ N)
        Ri = (data[:, a] ./ N) .* data[:, d]
        Si = (data[:, b] ./ N) .* data[:, c]
        R  = sum(Ri)
        S  = sum(Si)
        if R == 0 || S == 0 return ConfInt(NaN, NaN, NaN) end

        est = log(R/S)
        se  = sqrt(1/2 * (sum(Pi .* Ri)/R^2 + sum(Pi .* Si + Qi .* Ri)/(R*S) + sum(Qi .* Si)/S^2)) # based on Robins et al. (1986)
        #zval= est / se
        #pval= 2*(1-cdf(Normal(), abs(zval)))
        if logscale return ConfInt(est - z*se, est + z*se, est, alpha) else return ConfInt(exp(est - z*se), exp(est + z*se), exp(est), alpha) end
    elseif type == :rr
        #...
        R = sum(data[:, a] .* (n2 ./ N))
        S = sum(data[:, c] .* (n1 ./ N))

        if sum(data[:, a]) == 0 || sum(data[:, c]) == 0 return ConfInt(NaN, NaN, NaN) end

        est = log(R/S)
        se  = sqrt(sum(((n1 ./ N) .* (n2 ./ N) .* (data[:, a] + data[:, c]) - (data[:, a] ./ N) .* data[:, c])) / (R*S))
        #zval= est / se
        #pval= 2*(1-cdf(Normal(), abs(zval)))
        if logscale return ConfInt(est - z*se, est + z*se, est, alpha) else return ConfInt(exp(est - z*se), exp(est + z*se), exp(est), alpha) end
    end
end

function powertostint(α::Real,  θ₁::Real, θ₂::Real, δ::Real, σ::Real, n::Int, design::Symbol, method::Symbol)::Float64
    d     = Design(design) #dffunc if generic funtion with 1 arg return df
    df    = d.df(n)
    if df < 1 throw(ArgumentError("powertostint: df < 1")) end
    σ̵ₓ::Float64 = σ*sediv(d, n)
    if method     == :owenq
        return   powertost_owenq(α, θ₁, θ₂, δ, σ̵ₓ, df)
    elseif method == :nct
        return     powertost_nct(α, θ₁, θ₂, δ, σ̵ₓ, df)
    elseif method == :mvt
        return     powertost_mvt(α, θ₁, θ₂, δ, σ̵ₓ, df) #not implemented
    elseif method == :shifted
        return powertost_shifted(α, θ₁, θ₂, δ, σ̵ₓ, df)
    else
         throw(ArgumentError("method not known!"))
    end
end


function twoprop(x1::Int, n1::Int, x2::Int, n2::Int; alpha=0.05, type::Symbol, method::Symbol)::ConfInt
    if alpha >= 1.0 || alpha <= 0.0 throw(ArgumentError("Alpha shold be > 0.0 and < 1.0")) end
    if type==:diff
        if method ==:nhs
            return propdiffnhsci(x1, n1, x2, n2, alpha)
        elseif method ==:nhscc
            return propdiffnhsccci(x1, n1, x2, n2, alpha)
        elseif method ==:ac
            return propdiffacci(x1, n1, x2, n2, alpha)
        elseif method ==:mn
            return propdiffmnci(x1, n1, x2, n2, alpha)
        elseif method ==:mee2
            return propdiffmeeci(x1, n1, x2, n2, alpha)
        elseif method ==:mee || method == :fm
            return propdifffmci(x1, n1, x2, n2, alpha)
        elseif method ==:wald
            return propdiffwaldci(x1, n1, x2, n2, alpha)
        elseif method ==:waldcc
            return propdiffwaldccci(x1, n1, x2, n2, alpha)
        end
    elseif type==:rr
        if method==:mn
            return proprrmnci(x1, n1, x2, n2, alpha)
        elseif method == :cli || method == :walters
            return proprrclici(x1, n1, x2, n2, alpha)
        elseif method == :li || method == :katz
            return proprrkatzci(x1, n1, x2, n2, alpha)
        elseif method ==:mover
            return  proprrmoverci(x1, n1, x2, n2, alpha)
        end
    elseif type==:or
        if method==:mn
            return propormnci(x1, n1, x2, n2, alpha)
        elseif method==:awoolf || method==:gart
            return proporawoolfci(x1, n1, x2, n2, alpha)
        elseif method==:woolf
            return proporwoolfci(x1, n1, x2, n2, alpha)
        elseif method==:mover
            return propormoverci(x1, n1, x2, n2, alpha)
        elseif method==:mn2
            return proporci(x1, n1, x2, n2, alpha)
        end
    end
end #twoProp


function designProp(type::Symbol)::Tuple{Function, Float64, Int}
    if type == :parallel
        #function f1(n) n - 2 end
        #return f1, 1.0, 2
        return x -> x - 2.0, 1.0, 2
    elseif type == :d2x2
        #function f2(n) n - 2 end
        #return f2, 0.5, 2
        return x -> x - 2.0, 0.5, 2
    elseif type == :d2x2x3
        #return function f3(n) 2*n - 3 end, 0.375, 2
        return x -> 2.0 * x - 3.0, 0.375, 2
    elseif type == :d2x2x4
        #return function f4(n) 3*n - 4 end, 0.25, 2
        return x -> 3.0 * x - 4.0, 0.25, 2
    elseif type == :d2x4x4
        #return function f5(n) 3*n - 4 end, 0.0625, 4
        return x -> 3.0 * x - 4.0, 0.0625, 4
    elseif type == :d2x3x3
        #return function f6(n) 2*n - 3 end, 1/6, 3
        return x -> 2.0 * x - 3.0, 1/6, 3
    elseif type == :d2x4x2
        #return function f7(n) n - 2 end, 0.5, 4
        return x -> x - 2.0, 0.5, 4
    elseif type == :d3x3
        #return function f8(n) 2*n - 4 end, 2/9, 3
        return x -> 2.0 * x - 4.0, 2/9, 3
    elseif type == :d3x6x3
        #return function f9(n) 2*n - 4 end, 1/18, 6
        return x -> 2.0 * x - 4.0, 1/18, 6
    else throw(CTUException(1031,"designProp: design not known!")) end
end

#Exceptions
struct CTUException <: Exception
    n::Int
    var::String
end

Base.showerror(io::IO, e::CTUException) = print("CTU Exception code: ", e.n, " Message: ", e.var);

#Descriptives
function descriptive_deprecated(data::DataFrame;
    sort::Union{Symbol, Array{T, 1}, Nothing} = nothing,
    vars::Union{Symbol, Array{T, 1},  Nothing} = nothing,
    stats::Union{Symbol, Array{T, 1}} = [:n, :mean, :sd, :sem, :uq, :median, :lq])::DataFrame where T <: Union{Symbol, String}
    #Filtering
    dfnames = names(data) # Col names of dataframe
    #Filter sort
    if isa(sort, Array)
        sort = Symbol.(sort)
        filter!(x->x in dfnames, sort)
        if length(sort) == 0 sort = nothing end
    elseif isa(sort, Symbol)
        if any(x -> x == sort, dfnames)
            sort = [sort]
        else
            sort =  nothing
        end
    else
        sort =  nothing
    end
    #Filters vars
    if isa(vars, Symbol) vars = [vars] end
    if isa(vars, Array)
        vars = Symbol.(vars)
        filter!(x->x in dfnames, vars)
        filter!(x-> !(x in sort), vars)
        if length(vars) > 0
            for i = 1:length(vars)
                if !(eltype(data[:, vars[1]]) <: Real) deleteat!(vars, i) end
            end
        end
        if length(vars) == 0 error("No variables of type Real[] in dataset found! Check vars array!") end
    else
        vars = Array{Symbol,1}(undef, 0)
        for i in dfnames
            if eltype(data[:, i]) <: Real push!(vars, i) end
        end
        if sort !== nothing filter!(x-> !(x in sort), vars) end
        if length(vars) == 0 error("Real[] columns not found!") end
    end
    #Filter statistics array
    if stats == :all
        stats = [:n, :min, :max, :range, :mean, :var, :sd, :sem, :cv, :harmmean, :geomean, :geovar, :geosd, :geocv, :skew, :kurt, :uq, :median, :lq, :iqr, :mode]
    elseif isa(stats, Array)
        stats = Symbol.(stats)
        filter!(x->x in [:n, :min, :max, :range, :mean, :var, :sd, :sem, :cv, :harmmean, :geomean, :geovar, :geosd, :geocv, :skew, :ses, :kurt, :sek, :uq, :median, :lq, :iqr, :mode], stats)
        if length(stats) == 0
            stats = [:n, :mean, :sd, :sem, :uq, :median, :lq]
            @warn "Error in stats, default used."
        end
    else
        stats = [:n, :mean, :sd, :sem, :uq, :median, :lq]
        @warn "Unknown stats, default used."
    end
    #End filtering

    #construct datasets
    dfs   = DataFrame(vars = Symbol[]) #Init DataFrames
    sdata = DataFrame()                #Temp dataset for sorting

    if isa(sort, Array)
        for i in sort                  #if sort - make sort columns in dataset
            dfs[:, i] = eltype(data[:, i])[]
        end
    end
    for i in stats                     #make columns for statistics
        dfs[:, i] = Real[]
    end
    for i in vars                      #var columns for temp dataset
        sdata[:, i] = Real[]
    end
    if !isa(sort, Array)               #if no sort - for vars without sorting
        for v in vars
            sarray = Array{Any,1}(undef, 0)
            deleteat!(data[:, v], findall(x->x === NaN || x === nothing, data[:, v]))
            if length(data[:, v]) > 1
                push!(sarray, v)
                append!(sarray, descriptive_(data[:, v], stats))
                push!(dfs, sarray)
            end
        end
        return dfs
    end

    #If sort exist
    sortlist = unique(data[:, sort])
    sort!(sortlist, tuple(sort))
    for i = 1:size(sortlist, 1) #For each line in sortlist
            if size(sdata, 1) > 0 deleterows!(sdata, 1:size(sdata, 1)) end
            for c = 1:size(data, 1) #For each line in data
                if data[c, sort] == sortlist[i,:]
                    push!(sdata, data[c, vars])  #tmp ds constructing
                end
            end
        for v in vars #For each variable in list
            sarray = Array{Any,1}(undef, 0)
            push!(sarray, v)
            for c in sort
                push!(sarray, sortlist[i, c])
            end
            append!(sarray, descriptive_(sdata[:, v], stats))
            push!(dfs, sarray)
        end
    end
    return dfs
end
=#
#=
struct Probability{T<: Float64, N <: Union{Int, Nothing}} <: AbstractParameter
    p::T
    n::N
    function Probability(p::T, n::N) where T <: Float64 where N <: Int
        new{T, N}(p, n)::Probability
    end
    function Probability(p::T, n::N = nothing) where T <: Float64 where N <: Nothing
        new{T, N}(p, nothing)::Probability
    end
end
=#
#=
function nca(data::DataFrame; conc = NaN, effect = NaN, time=:Time, sort = NaN, calcm = :lint, bl::Float64 = NaN, th::Float64 = NaN)::NCA
    columns  = DataFrames.names(data)
    errorlog = ""
    errorn   = 0
    errors   = Int[]
    if typeof(findfirst(x ->  x  == conc, columns)) == Nothing && typeof(findfirst(x ->  x  == effect, columns)) == Nothing
        errorn+=1
        errorlog *= string(errorn)*": Concentration column not found\n"
        push!(errors, 1)
    end
    if typeof(findfirst(x ->  x  == time, columns)) == Nothing
        errorn+=1
        errorlog *= string(errorn)*": Time column not found\n"
        push!(errors, 2)
    end
    if  !(calcm == :lint || calcm == :logt || calcm == :luld)
        errorn+=1
        errorlog *= string(errorn)*": Calculation method not found\n"
        push!(errors, 3)
    end
    if sort !== NaN
        if isa(sort, Array)
            for i = 1:length(sort)
                if typeof(findfirst(x ->  x  == sort[i], columns)) == Nothing
                    errorn+=1
                    errorlog *= string(errorn)*": Sort column not found\n"
                    push!(errors, 4)
                end
            end
        else
            errorn+=1
            errorlog *= string(errorn)*": Sort is not array\n"
            push!(errors, 5)
        end
    end
    if errorn > 0 return NCA(DataFrame(), DataFrame(), DataFrame(), "", errorlog, errors) end
    if isa(conc, String)  conc = Symbol(conc) end
    if isa(time, String)  time = Symbol(time) end

    res  = DataFrame()
    elim = DataFrame()

    if sort === NaN
        if effect !== NaN
            res   = ncapd(data; conc = effect, time = time, calcm = calcm, bl = bl, th = th)
        elseif conc !== NaN
            res, rsqn, elim = nca_(data, conc, time, calcm)
        end
    else
        #PK NCA
        if conc !== NaN
            sortlist = unique(data[:, sort])
            res  = DataFrame(AUClast = Float64[], Cmax = Float64[], Tmax = Float64[], AUMClast = Float64[], MRTlast = Float64[], Kel = Float64[], HL = Float64[], Rsq = Float64[], AUCinf = Float64[], AUCpct = Float64[])
            elim = DataFrame(Start = Int[], End = Int[], b = Float64[], a = Float64[], Rsq = Float64[])
            for i = 1:size(sortlist, 1) #For each line in sortlist
                datai = DataFrame(Concentration = Float64[], Time = Float64[])
                for c = 1:size(data, 1) #For each line in data
                    if data[c, sort] == sortlist[i,:]
                        push!(datai, [data[c, conc], data[c, time]])
                    end
                end
                ncares = nca_(datai, :Concentration, :Time, calcm)
                append!(res, ncares[1])
                if ncares[2] !== NaN
                    append!(elim, DataFrame(ncares[3][ncares[2], :]))
                else
                    append!(elim, DataFrame(Start = [0], End = [0], b = [NaN], a = [NaN], Rsq = [NaN]))
                end
            end
            elim = hcat(sortlist, elim)
        #PD NCA
        elseif effect !== NaN && bl !== NaN
            sortlist = unique(data[:, sort])
            res      = DataFrame(RMAX = Float64[], TH = Float64[], BL = Float64[], AUCABL = Float64[], AUCBBL = Float64[], AUCATH = Float64[], AUCBTH = Float64[], AUCBLNET = Float64[], AUCTHNET = Float64[], AUCDBLTH = Float64[], TABL = Float64[], TBBL = Float64[], TATH = Float64[], TBTH = Float64[])
            for i = 1:size(sortlist, 1) #For each line in sortlist
                datai = DataFrame(Concentration = Float64[], Time = Float64[])
                for c = 1:size(data, 1) #For each line in data
                    if data[c, sort] == sortlist[i,:]
                        #println(conc, time)
                        push!(datai, [data[c, effect], data[c, time]])
                    end
                end
                pdres  = ncapd(datai; conc = :Concentration, time = :Time, calcm = calcm, bl = bl, th = th)
                append!(res, pdres)
            end
        end
        res  = hcat(sortlist, res)
    end
    return NCA(res, elim, DataFrame(), "", "", [0])
end
=#
#=
function nca_(data::DataFrame, conc::Symbol, time::Symbol, calcm = :lint)
    if length(unique(data[:,time])) != length(data[:,time])
        return  DataFrame(AUClast = [NaN], Cmax = [NaN], Tmax = [NaN], AUMClast = [NaN], MRTlast = [NaN], Kel = [NaN], HL = [NaN], Rsq = [NaN], AUCinf = [NaN], AUCpct = [NaN])
    end
    ncares = Dict(:Nums => 0, :Cmax => NaN, :Tmax => NaN, :Tmaxn => 0, :Tlast => NaN, :Clast => NaN, :AUClast => NaN, :AUMClast => NaN, :MRTlast => NaN, :Kel => NaN, :HL => NaN, :Rsq => NaN, :AUCinf => NaN, :AUCpct => NaN)

    pklog = "NCA analysis: \n"
    pklog *= "Concentration column: "*string(conc)*"\n"
    pklog *= "Time column: "*string(time)*"\n"
    #sort!(data, [time])
    ncarule!(data, conc, time, LimitRule(eps(), 0, NaN))

    n::Int     = nrow(data)
    cmax = maximum(data[:, conc])               #Cmax
    clast =  data[n, conc]
    pklog *= "Cmax = "*string(cmax)*"\n"
    auc   = 0                                  #AUClast
    aumc  = 0                                 #AUMClast
    mrt   = 0
    tmax  = NaN                               #Tmax
    tmaxn::Int = 0                             #Time point Tmax
    kel   = NaN
    hl    = NaN
    rsq   = NaN
    rsqn  = NaN
    aucinf = NaN
    aucinfpct = NaN
    mrtlast = NaN
# --- Tmax ---
    cmax, tmax, tmaxn = ctmax(data, conc, time)
    # --- AUC AUMC ---
    for i = 2:nrow(data)
        if calcm == :lint
            auc_  = linauc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])    # (data[i - 1, conc] + data[i, conc]) * (data[i, time] - data[i - 1, time])/2
            aumc_ = linaumc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])  # (data[i - 1, conc]*data[i - 1, time] + data[i, conc]*data[i, time]) * (data[i, time] - data[i - 1, time])/2
        elseif calcm == :logt
            if i > tmaxn
                if data[i, conc] < data[i - 1, conc] &&  data[i, conc] > 0 &&  data[i - 1, conc] > 0
                    auc_  = logauc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                    aumc_ = logaumc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                else
                    auc_  = linauc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                    aumc_ = linaumc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                end
            else
                auc_ = linauc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                aumc_ = linaumc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
            end
        elseif calcm == :luld
            if data[i, conc] < data[i - 1, conc] &&  data[i, conc] > 0 &&  data[i - 1, conc] > 0
                auc_  = logauc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                aumc_ = logaumc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
            else
                auc_  = linauc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                aumc_ = linaumc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
            end
        end
        auc   += auc_
        aumc  += aumc_
        pklog *= "AUC("*string(i-1)*") = "*string(auc_)*"; Sum = "*string(auc)*"\n"
    end
    # --- MRT ---
    mrt    = aumc / auc
    pklog *= "AUClast = "*string(auc)*"\n"

    # --- Kel, HL, ets
    if n - tmaxn >= 3
        keldata = DataFrame(Start = Int[], End = Int[], b = Float64[], a = Float64[], Rsq = Float64[])
        logconc = log.(data[:, conc])
        for i::Int = tmaxn+1:n-2
            sl = slope(data[i:n,time], logconc[i:n])
            if sl[2] < 0
                push!(keldata, [i, n, sl[1], sl[2], sl[3]])
            end
        end
        if nrow(keldata) > 0
            rsq, rsqn = findmax(keldata[:, :Rsq])
            kel = abs(keldata[rsqn,:a])
            hl  = LOG2/kel
            aucinf = auc + clast/kel
            aucinfpct = (aucinf - auc) / aucinf * 100.0
        end
    else
        keldata = nothing
    end
    # arsq = 1 - (1-rsq)*(rsqn-1)/(rsqn-2)
    return DataFrame(AUClast = [auc], Cmax = [cmax], Tmax = [tmax], AUMClast = [aumc], MRTlast = [mrt], Kel = [kel], HL = [hl], Rsq = [rsq], AUCinf = [aucinf], AUCpct = [aucinfpct]), rsqn, keldata
end
=#
#=
function ncapd(data::DataFrame; conc::Symbol, time::Symbol, calcm::Symbol = :lint, bl::Real, th::Real)

    aucabl = 0
    aucbbl = 0
    aucath = 0
    aucbth = 0
    aucdblth = 0
    tabl     = 0
    tbbl     = 0
    tath     = 0
    tbth     = 0
    for i = 2:nrow(data)
        #BASELINE
        if data[i - 1, conc] <= bl && data[i, conc] <= bl
            tbbl   += data[i, time] - data[i - 1, time]
            aucbbl += linauc(data[i - 1, time], data[i, time], bl - data[i - 1, conc], bl - data[i, conc])
        elseif data[i - 1, conc] <= bl &&  data[i, conc] > bl
            tx      = linpredict(data[i - 1, conc], data[i, conc], bl, data[i - 1, time], data[i, time])
            tbbl   += tx - data[i - 1, time]
            tabl   += data[i, time] - tx
            aucbbl += (tx - data[i - 1, time]) * (bl - data[i - 1, conc]) / 2
            aucabl += (data[i, time] - tx) * (data[i, conc] - bl) / 2 #Ok
        elseif data[i - 1, conc] > bl &&  data[i, conc] <= bl
            tx      = linpredict(data[i - 1, conc], data[i, conc], bl, data[i - 1, time], data[i, time])
            tbbl   += data[i, time] - tx
            tabl   += tx - data[i - 1, time]
            aucbbl += (data[i, time] - tx) * (bl - data[i, conc]) / 2
            aucabl += (tx - data[i - 1, time]) * (data[i - 1, conc] - bl) / 2
        elseif data[i - 1, conc] > bl &&  data[i, conc] > bl
            tabl   += data[i, time] - data[i - 1, time]
            aucabl += linauc(data[i - 1, time], data[i, time], data[i - 1, conc] - bl, data[i, conc] - bl)
        end
        #THRESHOLD
        if th !== NaN
            if data[i - 1, conc] <= th && data[i, conc] <= th
                tbth   += data[i, time] - data[i - 1, time]
                aucbth += linauc(data[i - 1, time], data[i, time], th - data[i - 1, conc], th - data[i, conc])
            elseif data[i - 1, conc] <= th &&  data[i, conc] > th
                tx      = linpredict(data[i - 1, conc], data[i, conc], th, data[i - 1, time], data[i, time])
                tbth   += tx - data[i - 1, time]
                tath   += data[i, time] - tx
                aucbth += (tx - data[i - 1, time]) * (th - data[i - 1, conc]) / 2
                aucath += (data[i, time] - tx) * (data[i, conc] - th) / 2 #Ok
            elseif data[i - 1, conc] > th &&  data[i, conc] <= th
                tx      = linpredict(data[i - 1, conc], data[i, conc], th, data[i - 1, time], data[i, time])
                tbth   += data[i, time] - tx
                tath   += tx - data[i - 1, time]
                aucbth += (data[i, time] - tx) * (th - data[i, conc]) / 2
                aucath += (tx - data[i - 1, time]) * (data[i - 1, conc] - th) / 2
            elseif data[i - 1, conc] > th &&  data[i, conc] > th
                tath   += data[i, time] - data[i - 1, time]
                aucath += linauc(data[i - 1, time], data[i, time], data[i - 1, conc] - th, data[i, conc] - th)
            end
        end
    end
    if bl > th
        aucdblth = aucath - aucabl
    elseif bl < th
        aucdblth = aucabl - aucath
    end
    return DataFrame(RMAX = [maximum(data[!, conc])], TH = [th], BL = [bl], AUCABL = [aucabl], AUCBBL = [aucbbl], AUCATH = [aucath], AUCBTH = [aucbth], AUCBLNET = [aucabl-aucbbl], AUCTHNET = [aucath-aucbth], AUCDBLTH = [aucdblth], TABL = [tabl], TBBL = [tbbl], TATH = [tath], TBTH = [tbth])
end
=#
#=
struct NCA
    result::DataFrame
    elimination::DataFrame
    settings::DataFrame
    textout::String
    errorlog::String
    errors::Array
end
=#
