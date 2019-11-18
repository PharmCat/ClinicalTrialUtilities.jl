
function Base.show(io::IO, obj::CTask)
        println(io, obj.param)
        println(io, obj.llim)
        println(io, obj.ulim)
        println(io, obj.alpha)
        println(io, obj.hyp)
        println(io, obj.k)
        println(io, obj.objective)
end

function Base.show(io::IO, obj::TaskResult{CT}) where CT <:  CTask{T, D, H, O} where T where D where H where O <: Union{SampleSize, Power}
        println(io, objectivename(obj.task.objective))
        println(io,"-----------------------------------------")
        println(io,"  Parameter type: ",  paramname(obj.task.param))
        println(io,"  Design: ",  designinfo(obj.task.design))
        println(io,"  Hypothesis: ", obj.task.hyp)
        #println(io,"  Lower limit: ", round(obj.task.llim, sigdigits = 4))
        #println(io,"  Upper limit: ", round(obj.task.ulim, sigdigits = 4))
        println(io,"  Alpha: ", obj.task.alpha)
        showobjective(io, obj.task.objective)
        println(io, obj.task.param)
        println(io,"-----------------------------------------")
        showresult(io, obj)
end
function Base.show(io::IO, obj::TaskResult{CT}) where CT <:  CTask{T, D, Bioequivalence, O} where T where D where O <: Union{SampleSize, Power}
        println(io, objectivename(obj.task.objective))
        println(io,"-----------------------------------------")
        println(io,"  Parameter type: ",  paramname(obj.task.param))
        println(io,"  Design: ",  designinfo(obj.task.design))
        println(io,"  Hypothesis: ", obj.task.hyp)
        #println(io,"  Lower limit: ", round(obj.task.llim, sigdigits = 4))
        #println(io,"  Upper limit: ", round(obj.task.ulim, sigdigits = 4))
        println(io,"  Alpha: ", obj.task.alpha)
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
        elseif isa(p, Probability) return "One Proportion"
        else return "NA" end
end
function groupnum(p::T)::String where T <: AbstractParameter
        if typeof(p) <: AbstractCompositeProportion{R} where R <: Real || typeof(p) <: AbstractCompositeMean{R} where R <: Real
                return "One"
        else
                return "Two"
        end
end
function designinfo(d::OneGroup)::String
        return "One group design"
end
function designinfo(d::Parallel)::String
        return "Two parallel groups design"
end
function designinfo(d::Crossover)::String
        return "Crossover design"
end

function showresult(io, obj::TaskResult{CT}) where CT <:  CTask{T, D, H, O} where T where D <: OneGroup where H where O <: SampleSize
        println(io, "Sample size: $(ceil(obj.result))")
end
function showresult(io, obj::TaskResult{CT}) where CT <:  CTask{T, D, H, O} where T where D <: Crossover where H where O <: SampleSize
        println(io, "Sample size: $(ceil(obj.result))")
end
function showresult(io, obj::TaskResult{CT})  where CT <:  CTask{T, D, H, O} where T where D <: Parallel where H where O <: SampleSize
        println(io, "Sample size (k=$(obj.task.k)):")
        println(io, "  Group A: ", ceil(obj.result * obj.task.k), "  Group B: ", ceil(obj.result))
        print(io,   "  Total: ", (ceil(obj.result) + ceil(obj.result*obj.task.k)))
end
function showresult(io, obj::TaskResult{CT}) where CT <:  CTask{T, D, H, O} where T where D where H where O <: Power
        println(io, "Power: ", round(obj.result, sigdigits = 6))
end

function Base.show(io::IO, p::Proportion)
        print(io,"  Proportion: ", p.x, "/", p.n)
end
function Base.show(io::IO, p::Probability)
        print(io,"  Probability: ", p.p)
end

function Base.show(io::IO, dp::DiffProportion{Probability, R}) where R <: Real
        println(io, "  A   = ", dp.a.p)
        print(io,   "  Ref = ", dp.b)
end

function Base.show(io::IO, dp::DiffProportion{Probability, Probability})
        println(io, "  A = ", dp.a.p)
        print(io,   "  B = ", dp.b.p)
end
#=
function Base.show(io, dp::DiffProportion{Proportion})::String
        println(io,"  A: ", dp.a.x,"/",dp.a.n)
        println(io,"  B: ", dp.b.x,"/",dp.b.n)
end
=#
function Base.show(io::IO, dp::T) where T <: Union{DiffProportion{P, P}, OddRatio{P}, RiskRatio{P}} where P <: Proportion
        println(io,"  A: ", dp.a.x,"/",dp.a.n)
        print(io,"  B: ", dp.b.x,"/",dp.b.n)
end
function Base.show(io::IO, dp::T) where T <: Union{DiffProportion{P, P}, OddRatio{P}, RiskRatio{P}} where P <: Probability
        println(io,"  A: ", dp.a.p)
        print(io,"  B: ", dp.b.p)
end
function Base.show(io::IO, dm::DiffMean{T}) where T <: AbstractMean
        println(io,"  A: ", dm.a.m, " ± ", round(dm.a.sd, sigdigits = 4))
        print(io,"  B: ", dm.b.m, " ± ", round(dm.b.sd, sigdigits = 4))
end
function Base.show(io::IO, dm::DiffMean{T}) where T <: Real
        println(io,"  A: ", dm.a.m, " ± ", round(dm.a.sd, sigdigits = 4))
        print(io,"  Ref: ", dm.b)
end
function Base.show(io::IO, m::Mean{Nothing})
        print(io,"  Mean(SD): ", m.m, " ± ", m.sd)
end

function Base.show(io::IO, h::Equality)
        println(io,"Equality")
        println(io,"  H₀: A = B")
        print(io,  "  Hₐ: A ≠ B")
end
function Base.show(io::IO, h::Equivalence)
        println(io,"Equivalence")
        println(io,"  H₀: |A − B| ≥ $(round(mdiff(h), sigdigits = 4))")
        print(io,  "  Hₐ: |A − B| < $(round(mdiff(h), sigdigits = 4))")
end
function Base.show(io::IO, h::Bioequivalence)
        println(io,"BioEquivalence")
        println(io,"  H₀: A/B ≤ $(round(exp(h.llim), sigdigits = 4)) || A/B ≥ $(round(exp(h.ulim), sigdigits = 4))")
        print(io,  "  Hₐ: $(round(exp(h.llim), sigdigits = 4)) < A/B < $(round(exp(h.ulim), sigdigits = 4))")
end
function Base.show(io::IO, h::Superiority)
        println(io,"Superiority/Non-Inferiority")
        println(io,"  H₀: A − B ≤ $(round(h.diff, sigdigits = 4))")
        print(io,  "  Hₐ: A − B > $(round(h.diff, sigdigits = 4))")
end
function Base.show(io::IO, h::McNemars)
        println(io,"McNemar's Equality test")
end
function Base.show(io::IO, e::TaskEstimate)
        for i = 1:length(e)
                prinln("Group $i: $(ceil(e.est))")
        end
end

function Base.show(io::IO, obj::ConfInt)
        print(io, "Estimate: $(obj.estimate) ($(obj.lower) - $(obj.upper))")
end
