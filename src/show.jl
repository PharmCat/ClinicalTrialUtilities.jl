

function Base.show(io::IO, obj::TaskResult{T}) where T <: Union{SampleSizeTask, PowerTask, TOSTSampleSizeTask}
        if isa(obj.task, SampleSizeTask)
        println(io,"         Sample Size Estimation         ")
        elseif isa(obj.task, PowerTask)
        println(io,"            Power Estimation            ")
        elseif isa(obj.task, TOSTSampleSizeTask)
        swowtost(io, obj)
        return nothing
        end
        println(io,"----------------------------------------")
        showmod1(io, obj)
        println(io,"----------------------------------------")
        showmod2(io, obj)
        println(io,"----------------------------------------")
        showmod3(io, obj)
        println(io,"----------------------------------------")
        showmod4(io, obj)
end

function showmod1(io, obj)
        str::String = ""
        if obj.task.param == :mean str = "Mean" elseif obj.task.param == :prop str = "Proportion" elseif obj.task.param == :or str = "Odd Ratio" end
        println(io,"  Parameter type: ", str)
        if obj.task.group == :one str = "One" elseif  obj.task.group == :two str = "Two" else str = "NA" end
        println(io,"  Groups: ", str)
        if obj.task.type == :ea str = "Equality" elseif obj.task.type == :ei str = "Equivalence" elseif obj.task.type == :ns str = "Non-Inferiority/Superiority" elseif obj.task.type == :mcnm str = "McNemar's Equality test" elseif obj.task.type == :co str = "Crossover" end
        println(io,"  Hypothesis: ", str)
end
function showmod2(io, obj::TaskResult{SampleSizeTask})
        println(io,"  Alpha: ", obj.task.alpha, " Beta: ", obj.task.beta)
        if obj.task.group == :two
        println(io,"  Constant k: ", obj.task.k)
        end
end
function showmod2(io, obj::TaskResult{PowerTask})
        println(io,"  Alpha: ", obj.task.alpha, " N: ", obj.task.n)
        if obj.task.group == :two
        println(io,"  Constant k: ", obj.task.k)
        end
end
function showmod3(io, obj)
        if obj.task.param == :mean
        println(io,"  SD: ",obj.task.sd)
        end
        if obj.task.group == :one
        println(io,"  Null Value: ", obj.task.a, " Test Value: ", obj.task.b)
        else
        println(io,"  Group A Value: ", obj.task.a, " Group B Value: ", obj.task.b)
        end
        if (obj.task.type == :ei || obj.task.type == :ns)
        println(io,"  Difference: ", obj.task.diff)
        end
end
function showmod4(io, obj::TaskResult{SampleSizeTask})
        if obj.task.group == :one
        println(io,"  Number estimate: ", ceil(obj.result))
        else
        println(io,"  Group A: ", ceil(obj.result * obj.task.k), "  Group B: ", ceil(obj.result))
        println(io,"  Total: ", (ceil(obj.result) + ceil(obj.result*obj.task.k)))
        end
end
function showmod4(io, obj::TaskResult{PowerTask})
        println(io,"  Power estimate: ", round(obj.result, sigdigits = 6))
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
