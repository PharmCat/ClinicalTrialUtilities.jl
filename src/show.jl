
function Base.show(io::IO, obj::CTask)
        println(io, obj.param)
        println(io, obj.llim)
        println(io, obj.ulim)
        println(io, obj.alpha)
        println(io, obj.hyp)
        println(io, obj.k)
        println(io, obj.objective)
end

function Base.show(io::IO, obj::TaskResult{CT}) where CT <:  CTask{T, H, O} where T where H where O <: Union{SampleSize, Power}
        println(io, objectivename(obj.task.objective))
        println(io,"-----------------------------------------")
        println(io,"  Parameter type: ",  paramname(obj.task.param))
        println(io,"  Hypothesis: ", hypname(obj.task.hyp))
        println(io,"  Lower limit: ", round(obj.task.llim, sigdigits = 4))
        println(io,"  Upper limit: ", round(obj.task.ulim, sigdigits = 4))
        println(io,"  Alpha: ", obj.task.alpha)
        showobjective(io, obj.task.objective)
        println(io, obj.task.param)
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

function showresult(io, obj)
        println(io, "Estimate: ")
        if isa(obj.task.objective, SampleSize)
                if typeof(obj.task.param) <: AbstractTwoProportion || isa(obj.task.param, DiffMean)
                        println(io,"  K: ", obj.task.k)
                        println(io,"  Group A: ", ceil(obj.result * obj.task.k), "  Group B: ", ceil(obj.result))
                        print(io,"  Total: ", (ceil(obj.result) + ceil(obj.result*obj.task.k)))
                else
                        print(io,"  Total: ", ceil(obj.result))
                end
        elseif  isa(obj.task.objective, Power)
                println(io, "Estimate: ", round(obj.result, sigdigits = 6))
        end
end

function swowtost(io, obj)
        println(io,"           TOST sample size             ")
        println(io,"----------------------------------------")
        println(io,"  Design: ", obj.task.design)
        println(io,"  Method: ", obj.task.method)
        println(io,"  Logscale: ", obj.task.logscale)
        println(io,"----------------------------------------")
        println(io,"  Alpha: ", obj.task.alpha)
        println(io,"  Beta: ", obj.task.beta)
        println(io,"----------------------------------------")
        println(io,"  Lower limit: ", obj.task.llim)
        println(io,"  Upper limit: ", obj.task.ulim)
        println(io,"  GMR: ", obj.task.gmr)
        println(io,"  CV: ", obj.task.cv)
        println(io,"----------------------------------------")
        println(io,"  Sample Size: ", obj.result)
end

function Base.show(io::IO, p::Proportion)
        print(io,"  Proportion: ", p.x, "/", p.n)
end
function Base.show(io::IO, p::Probability)
        print(io,"  Probability: ", p.p)
end

function Base.show(io::IO, dp::DiffProportion{Probability})
        println(io, "  A: ", dp.a.p)
        print(io, "  B: ", dp.b.p)
end
#=
function Base.show(io, dp::DiffProportion{Proportion})::String
        println(io,"  A: ", dp.a.x,"/",dp.a.n)
        println(io,"  B: ", dp.b.x,"/",dp.b.n)
end
=#
function Base.show(io::IO, dp::T) where T <: Union{DiffProportion{P}, OddRatio{P}, RiskRatio{P}} where P <: Proportion
        println(io,"  A: ", dp.a.x,"/",dp.a.n)
        print(io,"  B: ", dp.b.x,"/",dp.b.n)
end
function Base.show(io::IO, dp::T) where T <: Union{DiffProportion{P}, OddRatio{P}, RiskRatio{P}} where P <: Probability
        println(io,"  A: ", dp.a.p)
        print(io,"  B: ", dp.b.p)
end
function Base.show(io::IO, dm::DiffMean)
        println(io,"  A: ", dm.a.m, " ± ", dm.a.sd)
        print(io,"  B: ", dm.b.m, " ± ", dm.b.sd)
end
function Base.show(io::IO, m::Mean{Nothing})
        print(io,"  Mean(SD): ", m.m, " ± ", m.sd)
end

function Base.show(io::IO, h::Equality)
        println(io,"  Equality (Two-Sided):")
        println(io,"  H₀: A = B")
        print(io,"  Hₐ: A ≠ B")
end
function Base.show(io::IO, h::Equivalence)
        println(io,"  Equivalence (One-Sided):")
        println(io,"  H₀: |A − B| ≥ δ")
        print(io,"  Hₐ: |A − B| < δ")
end
function Base.show(io::IO, h::Superiority)
        println(io,"  Superiority/Non-Inferiority (One-Sided):")
        println(io,"  H₀: A − B ≤ δ")
        print(io,"  Hₐ: A − B > δ")
end
