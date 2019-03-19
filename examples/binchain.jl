using Distributions, Random, Distributed, MultivariateStats, Plots, DataFrames, GLM

function binChain()

    rng = MersenneTwister(1234)
    Random.seed!(rng)

    BER1 = Bernoulli(0.75)

    obs = Array{Float64, 1}(undef,0)
    cum = Array{Float64, 1}(undef,0)
    n   = Array{Float64, 1}(undef,0)

    for i=1:10000
        push!(obs, rand(BER1))
        push!(cum, mean(obs))
        push!(n, i)
    end

    df = DataFrames.DataFrame(x = n, y = cum);

    reg =  glm(@formula(y ~ x), df, Normal(), wts=n)

    coefs = coef(reg)
    a = llsq(n[:,:], cum[:,:])
    f(x) = coefs[2]*x+coefs[1]
    #println("Linear Val: "*string(a[2]))
    #println("Mid Linear Val: "*string(a[1]*100+a[2]))
    #println("WLinear Val: "*string(coefs[1]))
    println("Mid WLinear Val: "*string(coefs[2]*100+coefs[1]))
    println("Val: "*string(cum[length(cum)]))
    #println(a[2])
    #cum[20:length(cum)]

    plot(cum, seriestype=:scatter, title="P", marker = (:hexagon, 1, 0.6, :red, stroke(0)), legend=false)
    plot!(f, seriestype=:line, marker = (:hexagon, 1, 0.6, :blue, stroke(0)), legend=false)
    #print(cum)

    #if coefs[2] > 1.111 || abs(cum[length(cum)]-0.75) > abs(coefs[2]*100+coefs[1]-0.75)

    #    return true else return false end
end
#rat=0
#for i=1:200
#    global rat
#    if binChain() rat+=1 end
#end
#print(rat)
binChain()
