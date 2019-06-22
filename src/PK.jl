module PK

import ..NCA

export nca

using DataFrames

    function nca(data; conc=:Concentration, time=:Time, sort = NaN, calcm = :lint, stack = false)
        columns  = DataFrames.names(data)
        errorlog = ""
        errorn   = 0
        errors   = Int[]

        if typeof(findfirst(x ->  x  == conc, columns)) == Nothing
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

        if sort === NaN
            nca, rsqn, elim = nca_(data, conc, time, calcm)
        else
            sortlist = unique(data[sort])
            nca  = DataFrame(AUClast = Float64[], Cmax = Float64[], Tmax = Float64[], AUMClast = Float64[], MRTlast = Float64[], Kel = Float64[], HL = Float64[], Rsq = Float64[], AUCinf = Float64[], AUCpct = Float64[])
            elim = DataFrame(Start = Int[], End = Int[], b = Float64[], a = Float64[], Rsq = Float64[])
            for i = 1:size(sortlist, 1) #For each line in sortlist
                datai = DataFrame(Concentration = Float64[], Time = Float64[])
                for c = 1:size(data, 1) #For each line in data
                    if data[c, sort] == sortlist[i,:]
                        push!(datai, [data[c, conc], data[c, time]])
                    end
                end
                ncares = nca_(datai, :Concentration, :Time, calcm)
                append!(nca, ncares[1])
                append!(elim, DataFrame(ncares[3][ncares[2],:]))
            end
            nca  = hcat(sortlist, nca)
            elim = hcat(sortlist, elim)
        end
        return NCA(nca, elim, DataFrame(), "", "", [0])
    end

    function nca_(data, conc, time, calcm = :lint)
        if length(unique(data[:,time])) != length(data[:,time])
            return  DataFrame(AUClast = [NaN], Cmax = [NaN], Tmax = [NaN], AUMClast = [NaN], MRTlast = [NaN], Kel = [NaN], HL = [NaN], Rsq = [NaN], AUCinf = [NaN], AUCpct = [NaN])
        end
        pklog = "NCA analysis: \n"
        pklog *= "Concentration column: "*string(conc)*"\n"
        pklog *= "Time column: "*string(time)*"\n"
        sort!(data, [time])
        pklog *= "Sorting by Time... \n"
        n::Int     = nrow(data)
        cmax = maximum(data[conc])               #Cmax
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
        aucinf = NaN
        aucinfpct = NaN
        mrtlast = NaN
    # --- Tmax ---
        for i = 1:nrow(data)
            if data[i, conc] == cmax
                tmax = data[i, time]
                tmaxn = i
                pklog *= "Tmax = "*string(tmax)*"\n"
                break
            end
        end


        pklog *= "AUC calculation:\n"
        # --- AUC AUMC ---
        for i = 2:nrow(data)
            if calcm == :lint
                auc_ = linauc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])    # (data[i - 1, conc] + data[i, conc]) * (data[i, time] - data[i - 1, time])/2
                aumc_ = linaumc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])  # (data[i - 1, conc]*data[i - 1, time] + data[i, conc]*data[i, time]) * (data[i, time] - data[i - 1, time])/2
            elseif calcm == :logt
                if i > tmaxn
                    if data[i, conc] < data[i - 1, conc] &&  data[i, conc] > 0 &&  data[i - 1, conc] > 0
                        auc_ = logauc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                        aumc_ = logaumc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                    else
                        auc_ = linauc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                        aumc_ = linaumc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                    end
                else
                    auc_ = linauc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                    aumc_ = linaumc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                end
            elseif calcm == :luld
                if data[i, conc] < data[i - 1, conc] &&  data[i, conc] > 0 &&  data[i - 1, conc] > 0
                    auc_ = logauc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                    aumc_ = logaumc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                else
                    auc_ = linauc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                    aumc_ = linaumc(data[i - 1, time], data[i, time], data[i - 1, conc], data[i, conc])
                end
            end
            auc += auc_
            aumc += aumc_
            pklog *= "AUC("*string(i-1)*") = "*string(auc_)*"; Sum = "*string(auc)*"\n"
        end
        # --- MRT ---
        mrt = aumc / auc
        pklog *= "AUClast = "*string(auc)*"\n"

        # --- Kel, HL, ets
        if n - tmaxn >= 3
            keldata = DataFrame(Start = Int[], End = Int[], b = Float64[], a = Float64[], Rsq = Float64[])
            logconc = log.(data[conc])
            for i::Int = tmax+1:n-2
                sl = slope(data[i:n,time], logconc[i:n])
                if sl[2] < 0
                    push!(keldata, [i, n, sl[1], sl[2], sl[3]])
                end
            end
            if nrow(keldata) > 0
                rsq, rsqn = findmax(keldata[:Rsq])
                kel = abs(keldata[rsqn,:a])
                hl  = 0.6931471805599453/kel
                aucinf = auc + clast/kel
                aucinfpct = (aucinf - auc) / aucinf * 100.0
            end
        end

        return DataFrame(AUClast = [auc], Cmax = [cmax], Tmax = [tmax], AUMClast = [aumc], MRTlast = [mrt], Kel = [kel], HL = [hl], Rsq = [rsq], AUCinf = [aucinf], AUCpct = [aucinfpct]), rsqn, keldata
    end

    function slope(x,y)
        n   = length(x)
        Σxy = sum(x .* y)
        Σx  = sum(x)
        Σy  = sum(y)
        Σx2 = sum( x .* x)
        Σy2 = sum( y .* y)

        a   = (Σy * Σx2 - Σx * Σxy)/(n * Σx2 - Σx^2)
        b   = (n * Σxy - Σx * Σy)/(n * Σx2 - Σx^2)
        r2  = (n * Σxy - Σx * Σy)^2/((n * Σx2 - Σx^2)*(n * Σy2 - Σy^2))
        return a, b, r2
    end #end slope
        #Linear trapezoidal auc
    function linauc(t1, t2, c1, c2)
        return (t2-t1)*(c1+c2)/2
    end
        #Linear trapezoidal aumc
    function linaumc(t1, t2, c1, c2)
        return (t2-t1)*(t1*c1+t2*c2)/2
    end
        #Log trapezoidal auc
    function logauc(t1, t2, c1, c2)
        return  (t2-t1)*(c2-c1)/log(c2/c1)
    end
        #Log trapezoidal aumc
    function logaumc(t1, t2, c1, c2)
        return (t2-t1) * (t2*c2-t1*c1) / log(c2/c1) - (t2-t1)^2 * (c2-c1) / log(c2/c1)^2
    end
end #end module
