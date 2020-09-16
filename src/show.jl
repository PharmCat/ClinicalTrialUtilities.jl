#-------------------------------------------------------------------------------
#Task
function Base.show(io::IO, obj::CTask)
        println(io, obj.param)
        println(io, "  Design: $(obj.design)")

        #println(io, "  Alpha: $(obj.hyp.alpha)")
        println(io, "  Hypothesis: $(obj.hyp)")
        println(io, "  K: $(obj.k)")
        print(io,   "  Objective: $(obj.objective)")
end

#-------------------------------------------------------------------------------
#Designs
function Base.show(io::IO, obj::Crossover{:d2x2})
        print(io, "2x2")
end
function Base.show(io::IO, obj::Crossover{:d2x2x3})
        print(io, "2x2x3")
end
function Base.show(io::IO, obj::Crossover{:d2x2x4})
        print(io, "2x2x4")
end
function Base.show(io::IO, obj::Crossover{:d2x4x4})
        print(io, "2x4x4")
end
function Base.show(io::IO, obj::Crossover{:d2x3x3})
        print(io, "2x3x3")
end
function Base.show(io::IO, obj::Crossover{:d3x3})
        print(io, "3x3")
end
function Base.show(io::IO, obj::Crossover{:d4x4})
        print(io, "4x4")
end
function Base.show(io::IO, obj::Crossover{:d3x6x3})
        print(io, "3x6x3")
end
function Base.show(io::IO, obj::Parallel)
        print(io, "Parallel")
end
function Base.show(io::IO, obj::OneGroup)
        print(io, "OneGroup")
end

#-------------------------------------------------------------------------------
# Objectives
function Base.show(io::IO, obj::SampleSize)
        print(io, "SampleSize (β: $(obj.val))")
end
function Base.show(io::IO, obj::Power)
        print(io, "Power (n: $(obj.val))")
end

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#Task result
function Base.show(io::IO, obj::TaskResult{CT}) where CT <:  CTask{T, D, H, O} where T where D where H where O <: Union{SampleSize, Power}
        println(io, objectivename(obj.task.objective))
        println(io,"-----------------------------------------")
        println(io,"  Parameter type: ",  paramname(obj.task.param))
        println(io,"  Design: $(obj.task.design)")
        println(io,"  Hypothesis: ", obj.task.hyp)
        #println(io,"  Lower limit: ", round(obj.task.llim, sigdigits = 4))
        #println(io,"  Upper limit: ", round(obj.task.ulim, sigdigits = 4))
        #println(io,"  Alpha: ", obj.task.alpha)
        showobjective(io, obj.task.objective)
        showparams(io, obj.task)
        #println(io, obj.task.param)
        println(io,"-----------------------------------------")
        showresult(io, obj)
end
function Base.show(io::IO, obj::TaskResult{CT}) where CT <:  CTask{T, D, Bioequivalence, O} where T where D where O <: Union{SampleSize, Power}
        println(io, objectivename(obj.task.objective))
        println(io,"-----------------------------------------")
        println(io,"  Parameter type: ",  paramname(obj.task.param))
        println(io,"  Design: $(obj.task.design)")
        println(io,"  Hypothesis: ", obj.task.hyp)
        #println(io,"  Lower limit: ", round(obj.task.llim, sigdigits = 4))
        #println(io,"  Upper limit: ", round(obj.task.ulim, sigdigits = 4))
        #println(io,"  Alpha: ", obj.task.alpha)
        showobjective(io, obj.task.objective)
        println(io, "  A/B = $(round(exp(obj.task.param.a.m - obj.task.param.b.m), sigdigits = 4))")
        println(io, "  σ  = $(round(obj.task.param.a.sd, sigdigits = 4))")
        println(io, "  CV  = $(round(cvfromsd(obj.task.param.a.sd), sigdigits = 4))")
        println(io,"-----------------------------------------")
        showresult(io, obj)
end

function objectivename(o::O)::String where O <: AbstractObjective
        if isa(o, SampleSize) return "         Sample Size Estimation         "
        elseif isa(o, Power)  return "            Power Estimation            "
        else return "NA" end

end
function showobjective(io, o::SampleSize)
        println(io, "  Beta: $(o.val) (Power: $((1-o.val)*100)%)")
end
function showobjective(io, o::Power)
        println(io, "  N: ", o.val)
end
function hypname(h::T)::String where T <: AbstractHypothesis
        if isa(h, Equivalence) return "Equivalence"
        elseif isa(h, Equality) return "Equality"
        elseif isa(h, Superiority) return "Superiority/Non-Inferiority"
        elseif isa(h, McNemars) return "McNemar's Equality test"
        else return "NA" end
end
function paramname(p::T)::String where T <: AbstractParameter
        if isa(p, DiffProportion) return "Proportion Difference"
        elseif isa(p, OddRatio) return "Odd Ratio"
        elseif isa(p, RiskRatio) return "Risk Ratio"
        elseif isa(p, DiffMean) return "Mean Difference"
        elseif isa(p, Mean) return "One Mean"
        elseif isa(p, Proportion) return "One Proportion"
        else return "NA" end
end
function groupnum(p::T)::String where T <: AbstractParameter
        if typeof(p) <: AbstractCompositeProportion{R} where R <: Real || typeof(p) <: AbstractCompositeMean{R} where R <: Real
                return "One"
        else
                return "Two"
        end
end

function showresult(io, obj::TaskResult{CT}) where CT <:  CTask{T, D, H, O} where T where D <: OneGroup where H where O <: SampleSize
        println(io, "Sample size: $(ceil(obj.result))")
end
function showresult(io, obj::TaskResult{CT}) where CT <:  CTask{T, D, H, O} where T where D <: Crossover where H where O <: SampleSize
        println(io, "Sample size: $(ceil(obj.result))")
end
function showresult(io, obj::TaskResult{CT}) where CT <:  CTask{T, D, H, O} where T where D <: Crossover where H <: Bioequivalence where O <: SampleSize
        println(io, "Sample size: $(ceil(obj.result)) (Power: $(round(obj.res[:pow], sigdigits = 4)))")
end
function showresult(io, obj::TaskResult{CT})  where CT <:  CTask{T, D, H, O} where T where D <: Parallel where H where O <: SampleSize
        println(io, "Sample size (k=$(obj.task.k)):")
        println(io, "  Group A: ", ceil(obj.result * obj.task.k), "  Group B: ", ceil(obj.result))
        print(io,   "  Total: ", (ceil(obj.result) + ceil(obj.result*obj.task.k)))
end
function showresult(io, obj::TaskResult{CT}) where CT <:  CTask{T, D, H, O} where T where D where H where O <: Power
        println(io, "Power: ", round(obj.result, sigdigits = 6))
end

function showparams(io, task::T) where T <: CTask{P, D, H, O} where P where D <: OneGroup where H where O
        println(io,"  A: ", task.param)
        println(io,"  B: ", refval(task.hyp))
end
function showparams(io, task::T) where T <: CTask{P, D, H, O} where P  where D where H where O
        println(io, task.param)
end
#-------------------------------------------------------------------------------
#Parameters
function Base.show(io::IO, p::Proportion{true})
        print(io, p.x, "/", p.n)
end
function Base.show(io::IO, p::Proportion{false})
        print(io, p.val)
end
#=
function Base.show(io::IO, p::Probability)
        print(io, p.p)
end
=#
#=
function Base.show(io::IO, dp::DiffProportion{Probability, R}) where R <: Real
        println(io, "  A   = ", dp.a.p)
        print(io,   "  Ref = ", dp.b)
end

function Base.show(io::IO, dp::DiffProportion{Probability, Probability})
        println(io, "  A = ", dp.a.p)
        print(io,   "  B = ", dp.b.p)
end
function Base.show(io, dp::DiffProportion{Proportion})::String
        println(io,"  A: ", dp.a.x,"/",dp.a.n)
        println(io,"  B: ", dp.b.x,"/",dp.b.n)
end
=#
function Base.show(io::IO, dp::T) where T <: AbstractCompositeProportion
        println(io,"  A: ", dp.a)
        print(io,"  B: ", dp.b)
end
#=
function Base.show(io::IO, dp::T) where T <: Union{DiffProportion{P, P}, OddRatio{P}, RiskRatio{P}} where P <: Probability
        println(io,"  A: ", dp.a.p)
        print(io,"  B: ", dp.b.p)
end
=#
function Base.show(io::IO, dm::DiffMean)
        println(io,"  A ", dm.a)
        print(io,"  B ", dm.b)
end
#=
function Base.show(io::IO, dm::DiffMean{false})
        println(io,"  A: ", dm.a.m, " ± ", round(dm.a.sd, sigdigits = 4))
        print(io,"  B: ", dm.b)
end
=#
function Base.show(io::IO, m::Mean{true})
        print(io,"Mean ± SD (N): $(m.m) ± $(m.sd) ($(m.n))")
end
function Base.show(io::IO, m::Mean{false})
        if m.sd === NaN
                print(io,"Mean: ", m.m)
        else
                print(io,"Mean(SD): ", m.m, " ± ", m.sd)
        end
end

#-------------------------------------------------------------------------------
# Hypothesis
function Base.show(io::IO, h::Equality)
        println(io,"Equality")
        println(io,"  H₀: A - B = 0")
        println(io,"  Hₐ: A - B ≠ 0")
        print(io,  "  Alpha: $(h.alpha)")
end
function Base.show(io::IO, h::Equivalence)
        println(io,"Equivalence")
        println(io,"  H₀: |A − B| ≥ $(round(mdiff(h), sigdigits = 4))")
        println(io,"  Hₐ: |A − B| < $(round(mdiff(h), sigdigits = 4))")
        print(io,  "  Alpha: $(h.alpha)")
end
function Base.show(io::IO, h::Bioequivalence)
        println(io,"BioEquivalence")
        println(io,"  H₀: A/B ≤ $(round(exp(h.llim), sigdigits = 4)) || A/B ≥ $(round(exp(h.ulim), sigdigits = 4))")
        println(io,"  Hₐ: $(round(exp(h.llim), sigdigits = 4)) < A/B < $(round(exp(h.ulim), sigdigits = 4))")
        print(io,  "  Alpha: $(h.alpha)")
end
function Base.show(io::IO, h::Superiority)
        println(io,"Superiority/Non-Inferiority")
        println(io,"  H₀: A − B ≤ $(round(h.diff, sigdigits = 4))")
        println(io,"  Hₐ: A − B > $(round(h.diff, sigdigits = 4))")
        print(io,  "  Alpha: $(h.alpha)")
end
function Base.show(io::IO, h::McNemars)
        println(io,"McNemar's Equality test")
        print(io,  "  Alpha: $(h.alpha)")
end
function Base.show(io::IO, e::TaskEstimate)
        for i = 1:length(e)
                prinln("Group $i: $(ceil(e.est))")
        end
end

#-------------------------------------------------------------------------------
# Confidence intervals
function Base.show(io::IO, obj::ConfInt)
        print(io, "Estimate: $(obj.estimate) ($(obj.lower) - $(obj.upper))")
end

#-------------------------------------------------------------------------------
# OUTPUT
function addspace(s::String, n::Int; first = false)::String
    if n > 0
        for i = 1:n
            if first s = Char(' ') * s else s = s * Char(' ') end
        end
    end
    return s
end
function printmatrix(io::IO, m::Matrix)
    sm = string.(m)
    lv = maximum(length.(sm), dims = 1)
    for r = 1:size(sm, 1)
        for c = 1:size(sm, 2)
            print(io, addspace(sm[r,c], lv[c] - length(sm[r,c]))*"   ")
        end
        println(io, "")
    end
end
function printsortval(io::IO, d::Dict)
   println(io, "Sort:")
   for v in d
      print(io, " $(v[1]) => $(v[2]),")
   end
   println(io, "\b ")
end
#-------------------------------------------------------------------------------
